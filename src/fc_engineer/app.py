"""Flet local desktop GUI and state controller for Fc Engineering Studio."""

import atexit
import copy
import logging
import math
import os
import threading
from typing import Any, Dict, List, Optional, Tuple

import flet as ft
import flet.canvas as cv
import pyperclip

from fc_engineer import config
from fc_engineer.core import CONSTANT_DOMAIN_EU_RANGES
from fc_engineer.designer import (
    ALLOWED_TEMPLATE_KEYS_BY_TYPE,
    DEFAULT_TEMPLATE_KEYS,
    DOMAIN_STYLE_MAP,
    MAX_CHAINS,
    MAX_COORD,
    MAX_CUSTOM_SEQ_LENGTH,
    MAX_DOMAINS_PER_CHAIN,
    MAX_METADATA_LENGTH,
    MAX_PAIRINGS,
    VALID_DOMAIN_TYPES,
    VALID_PAIRING_TYPES,
    AntibodyProject,
    DomainBlock,
    InterChainPairing,
    PolypeptideChain,
    SchematicScene,
    add_file_custom_mutation,
    analyze_chain_developability,
    assemble_chain_sequence,
    build_format_template,
    build_preset_map,
    build_schematic_scene,
    diff_chain_sequences,
    reflow_chain_coordinates,
    remove_chain_from_project,
    remove_domain_from_chain,
    remove_preset_from_chain,
    resolve_domain_sequence,
    route_mutations_to_chain,
)
from fc_engineer.exporters import (
    format_project_fasta,
    load_ab_file,
    save_ab_file,
    save_fasta_file,
    save_schematic_png,
)

logger = logging.getLogger(__name__)

_FORMAT_LABELS: List[Tuple[str, str]] = [
    ("igg", "Monospecific IgG"),
    ("crossmab", "CrossMab"),
    ("igg_scfv", "IgG-scFv"),
    ("scfv_fc", "scFv-Fc"),
    ("dvd_ig", "DVD-Ig"),
    ("bite", "BiTE"),
    ("custom", "Custom Construct"),
]


class AntibodyStudioController:
    """Decoupled workspace controller managing antibody design state and security lifecycle."""

    def __init__(
        self,
        initial_format: str = "igg",
        isotype: str = "igg1",
        allotype: Optional[str] = None,
    ) -> None:
        self._clipboard_lock = threading.RLock()
        self.copied_fasta: str = ""
        self._clipboard_timer: Optional[threading.Timer] = None
        self._clipboard_gen: int = 0
        self.apply_to_all_heavy: bool = False
        self.selected_chain_id: Optional[str] = None
        self.selected_domain_id: Optional[str] = None
        self.project: AntibodyProject = build_format_template(
            initial_format, isotype=isotype, allotype=allotype
        )
        self._reset_selection()

    def _reset_selection(self) -> None:
        if self.project.chains:
            first_chain = self.project.chains[0]
            self.selected_chain_id = first_chain.id
            self.selected_domain_id = (
                first_chain.domains[0].id if first_chain.domains else None
            )
        else:
            self.selected_chain_id = None
            self.selected_domain_id = None

    def _get_chain(self, chain_id: str) -> PolypeptideChain:
        for chain in self.project.chains:
            if chain.id == chain_id:
                return chain
        raise ValueError(f"Chain '{chain_id}' not found.")

    def _get_domain(
        self, chain_id: str, domain_id: str
    ) -> Tuple[PolypeptideChain, DomainBlock, int]:
        chain = self._get_chain(chain_id)
        for idx, dom in enumerate(chain.domains):
            if dom.id == domain_id:
                return chain, dom, idx
        raise ValueError(f"Domain '{domain_id}' not found on chain '{chain_id}'.")

    def select_domain(
        self, chain_id: Optional[str], domain_id: Optional[str] = None
    ) -> None:
        """Select the active chain/domain inspected in the right panel."""
        if chain_id is None:
            self.selected_chain_id = None
            self.selected_domain_id = None
            return
        chain = self._get_chain(chain_id)
        self.selected_chain_id = chain.id
        if domain_id is not None:
            for dom in chain.domains:
                if dom.id == domain_id:
                    self.selected_domain_id = dom.id
                    return
            raise ValueError(f"Domain '{domain_id}' not found on chain '{chain_id}'.")
        self.selected_domain_id = chain.domains[0].id if chain.domains else None

    def select_format(
        self,
        format_type: str,
        isotype: str = "igg1",
        allotype: Optional[str] = None,
    ) -> AntibodyProject:
        """Switch to a new antibody format template and clear clipboard state."""
        new_proj = build_format_template(format_type, isotype=isotype, allotype=allotype)
        try:
            self.project = new_proj
            self._reset_selection()
        finally:
            self.clear_clipboard()
        return self.project

    def rename_project(self, name: str) -> None:
        """Update the project display name after control-character and length checks."""
        if not isinstance(name, str):
            raise ValueError("Project name must be a string.")
        if any(ord(ch) < 0x20 or ord(ch) == 0x7F for ch in name):
            raise ValueError("Project name contains disallowed control characters.")
        stripped = name.strip()
        if not stripped:
            raise ValueError("Project name cannot be empty.")
        if len(stripped) > MAX_METADATA_LENGTH:
            raise ValueError(
                f"Project name exceeds maximum length of {MAX_METADATA_LENGTH} characters."
            )
        self.project.name = stripped

    def add_chain(self, name: Optional[str] = None) -> PolypeptideChain:
        """Append a new empty PolypeptideChain using the lowest unused Chain_{n} ID."""
        if len(self.project.chains) >= MAX_CHAINS:
            raise ValueError(f"Maximum of {MAX_CHAINS} chains reached.")
        existing_ids = {c.id for c in self.project.chains}
        n = 1
        while f"Chain_{n}" in existing_ids:
            n += 1
        chain_id = f"Chain_{n}"

        if name is None:
            resolved_name = chain_id
        else:
            if not isinstance(name, str):
                raise ValueError("Chain name must be a string.")
            if any(ord(ch) < 0x20 or ord(ch) == 0x7F for ch in name):
                raise ValueError("Chain name contains disallowed control characters.")
            resolved_name = name.strip()
            if not resolved_name or len(resolved_name) > MAX_METADATA_LENGTH:
                raise ValueError("Invalid chain name length.")

        chain = PolypeptideChain(
            id=chain_id,
            name=resolved_name,
            domains=[],
            x=40.0,
            y=40.0 + len(self.project.chains) * 96.0,
        )
        self.project.chains.append(chain)
        self.selected_chain_id = chain.id
        self.selected_domain_id = None
        return chain

    def remove_chain(self, chain_id: str) -> None:
        """Remove a polypeptide chain and update current selection safely."""
        self._get_chain(chain_id)
        remove_chain_from_project(self.project, chain_id)
        if self.selected_chain_id == chain_id:
            self._reset_selection()

    def _next_domain_id(self, chain_id: str, domain_type: str) -> str:
        existing_ids = {d.id for c in self.project.chains for d in c.domains}
        prefix = f"{chain_id.lower()}_{domain_type.lower()}"
        n = 1
        while f"{prefix}_{n}" in existing_ids:
            n += 1
        return f"{prefix}_{n}"

    def add_domain_to_chain(
        self, chain_id: str, domain_type: str, index: Optional[int] = None
    ) -> DomainBlock:
        """Insert or append a DomainBlock onto a chain and reflow unpinned coords."""
        chain = self._get_chain(chain_id)
        if domain_type not in VALID_DOMAIN_TYPES:
            raise ValueError(f"Unsupported domain type '{domain_type}'.")
        if len(chain.domains) >= MAX_DOMAINS_PER_CHAIN:
            raise ValueError(
                f"Chain '{chain_id}' reached maximum of {MAX_DOMAINS_PER_CHAIN} domains."
            )

        is_const = domain_type in CONSTANT_DOMAIN_EU_RANGES
        iso_map = config.SEQUENCES.get("igg1", {})
        default_allotype = next(iter(iso_map)) if (is_const and iso_map) else ""
        dom = DomainBlock(
            id=self._next_domain_id(chain.id, domain_type),
            domain_type=domain_type,
            isotype="igg1",
            allotype=default_allotype,
            template_key=DEFAULT_TEMPLATE_KEYS.get(domain_type, "") if not is_const else "",
        )
        if index is None:
            chain.domains.append(dom)
        else:
            clamped = max(0, min(len(chain.domains), index))
            chain.domains.insert(clamped, dom)

        reflow_chain_coordinates(chain)
        self.selected_chain_id = chain.id
        self.selected_domain_id = dom.id
        return dom

    def reorder_domain(
        self, chain_id: str, from_index: int, to_index: int
    ) -> None:
        """Reorder a domain within a single chain and reflow unpinned coords."""
        chain = self._get_chain(chain_id)
        n = len(chain.domains)
        if not (0 <= from_index < n) or not (0 <= to_index < n):
            raise ValueError("Domain reorder index out of bounds.")
        dom = chain.domains.pop(from_index)
        chain.domains.insert(to_index, dom)
        reflow_chain_coordinates(chain)
        self.selected_chain_id = chain.id
        self.selected_domain_id = dom.id

    def move_domain_across_chains(
        self,
        source_chain_id: str,
        domain_id: str,
        target_chain_id: str,
        target_index: Optional[int] = None,
    ) -> None:
        """Move or reorder a DomainBlock within the same chain or onto another chain."""
        src_chain, dom, src_idx = self._get_domain(source_chain_id, domain_id)
        dst_chain = self._get_chain(target_chain_id)

        if source_chain_id == target_chain_id:
            src_chain.domains.pop(src_idx)
            if target_index is None:
                src_chain.domains.append(dom)
            else:
                adj = target_index - 1 if target_index > src_idx else target_index
                clamped = max(0, min(len(src_chain.domains), adj))
                src_chain.domains.insert(clamped, dom)
            reflow_chain_coordinates(src_chain)
            self.selected_chain_id = src_chain.id
            self.selected_domain_id = dom.id
            return

        if len(dst_chain.domains) >= MAX_DOMAINS_PER_CHAIN:
            raise ValueError(
                f"Target chain '{target_chain_id}' reached maximum of {MAX_DOMAINS_PER_CHAIN} domains."
            )

        src_chain.domains.pop(src_idx)
        # Prune any pairings referencing the old (source_chain_id, domain_id) endpoint
        self.project.pairings = [
            p
            for p in self.project.pairings
            if not (
                (p.chain_a_id == source_chain_id and p.domain_a_id == domain_id)
                or (p.chain_b_id == source_chain_id and p.domain_b_id == domain_id)
            )
        ]

        if target_index is None:
            dst_chain.domains.append(dom)
        else:
            clamped = max(0, min(len(dst_chain.domains), target_index))
            dst_chain.domains.insert(clamped, dom)

        reflow_chain_coordinates(src_chain)
        reflow_chain_coordinates(dst_chain)
        self.selected_chain_id = dst_chain.id
        self.selected_domain_id = dom.id

    def remove_domain(self, chain_id: str, domain_id: str) -> None:
        """Delete a domain from a chain and clean up any orphaned pairings."""
        chain, _, _ = self._get_domain(chain_id, domain_id)
        remove_domain_from_chain(self.project, chain_id, domain_id)
        if self.selected_chain_id == chain_id and self.selected_domain_id == domain_id:
            self.selected_domain_id = chain.domains[0].id if chain.domains else None

    def move_domain(
        self, chain_id: str, domain_id: str, x: float, y: float
    ) -> None:
        """Reposition a domain block manually on the 2D canvas and pin it."""
        _, dom, _ = self._get_domain(chain_id, domain_id)
        if (
            isinstance(x, bool)
            or isinstance(y, bool)
            or not isinstance(x, (int, float))
            or not isinstance(y, (int, float))
        ):
            raise ValueError("Domain coordinates must be numeric.")
        fx, fy = float(x), float(y)
        if not math.isfinite(fx) or not math.isfinite(fy):
            raise ValueError("Domain coordinates must be finite.")
        dom.x = max(-MAX_COORD, min(MAX_COORD, fx))
        dom.y = max(-MAX_COORD, min(MAX_COORD, fy))
        dom.pinned = True

    def update_domain_properties(
        self, chain_id: str, domain_id: str, **kwargs: Any
    ) -> List[str]:
        """Update allowed DomainBlock attributes with rollback on resolution failure."""
        allowed = {"isotype", "allotype", "template_key", "custom_sequence", "pinned"}
        unknown = set(kwargs.keys()) - allowed
        if unknown:
            raise ValueError(f"Unsupported domain properties: {', '.join(sorted(unknown))}")
        if "pinned" in kwargs and not isinstance(kwargs["pinned"], bool):
            raise ValueError("pinned must be a boolean.")

        _, dom, _ = self._get_domain(chain_id, domain_id)
        snapshot = (
            dom.isotype,
            dom.allotype,
            dom.template_key,
            dom.custom_sequence,
            dom.pinned,
        )

        if "isotype" in kwargs and kwargs["isotype"] != dom.isotype and "allotype" not in kwargs:
            dom.allotype = ""

        for key, val in kwargs.items():
            setattr(dom, key, val)

        _, errs = resolve_domain_sequence(dom, self.get_preset_map())
        if errs:
            (
                dom.isotype,
                dom.allotype,
                dom.template_key,
                dom.custom_sequence,
                dom.pinned,
            ) = snapshot
            return errs
        return []

    def _resolve_mutation_target_chains(
        self, chain_id: str, apply_to_all_heavy: Optional[bool]
    ) -> List[PolypeptideChain]:
        selected = self._get_chain(chain_id)
        use_all = (
            self.apply_to_all_heavy
            if apply_to_all_heavy is None
            else bool(apply_to_all_heavy)
        )
        if not use_all:
            return [selected]
        heavy = [
            c
            for c in self.project.chains
            if any(d.domain_type in ("CH2", "CH3") for d in c.domains)
        ]
        return heavy if heavy else [selected]

    def toggle_preset(
        self,
        chain_id: str,
        preset_label: str,
        enabled: bool,
        apply_to_all_heavy: Optional[bool] = None,
    ) -> List[str]:
        """Enable or disable an Fc mutation preset atomically across target chain(s)."""
        targets = self._resolve_mutation_target_chains(chain_id, apply_to_all_heavy)
        preset_map = self.get_preset_map()

        if not enabled:
            for c in targets:
                remove_preset_from_chain(c, preset_label)
            return []

        if preset_label not in preset_map:
            return [f"Unknown mutation preset '{preset_label}'."]
        mut_str = preset_map[preset_label]

        # Atomic dry-run across copies of all target chains
        for c in targets:
            errs = route_mutations_to_chain(
                copy.deepcopy(c), mut_str, preset_label=preset_label, preset_map=preset_map
            )
            if errs:
                return errs

        for c in targets:
            route_mutations_to_chain(
                c, mut_str, preset_label=preset_label, preset_map=preset_map
            )
        return []

    def apply_custom_mutation(
        self,
        chain_id: str,
        mutants_str: str,
        apply_to_all_heavy: Optional[bool] = None,
    ) -> List[str]:
        """Apply ad-hoc EU point mutations atomically across target chain(s)."""
        targets = self._resolve_mutation_target_chains(chain_id, apply_to_all_heavy)
        preset_map = self.get_preset_map()

        for c in targets:
            errs = route_mutations_to_chain(
                copy.deepcopy(c), mutants_str, preset_label="", preset_map=preset_map
            )
            if errs:
                return errs

        for c in targets:
            route_mutations_to_chain(
                c, mutants_str, preset_label="", preset_map=preset_map
            )
        return []

    def register_file_custom_mutation(self, label: str, mutants_str: str) -> List[str]:
        """Register a project-scoped custom mutation definition."""
        return add_file_custom_mutation(self.project, label, mutants_str)

    def add_pairing(
        self,
        chain_a_id: str,
        domain_a_id: str,
        chain_b_id: str,
        domain_b_id: str,
        pairing_type: str,
    ) -> InterChainPairing:
        """Create a validated inter-chain or intra-chain pairing between two domains."""
        if len(self.project.pairings) >= MAX_PAIRINGS:
            raise ValueError(f"Maximum of {MAX_PAIRINGS} pairings reached.")
        if pairing_type not in VALID_PAIRING_TYPES:
            raise ValueError(f"Invalid pairing type '{pairing_type}'.")
        self._get_domain(chain_a_id, domain_a_id)
        self._get_domain(chain_b_id, domain_b_id)
        if (chain_a_id, domain_a_id) == (chain_b_id, domain_b_id):
            raise ValueError("Cannot pair a domain with itself.")

        existing_ids = {p.id for p in self.project.pairings}
        n = 1
        while f"pair_{n}" in existing_ids:
            n += 1
        pairing = InterChainPairing(
            id=f"pair_{n}",
            chain_a_id=chain_a_id,
            domain_a_id=domain_a_id,
            chain_b_id=chain_b_id,
            domain_b_id=domain_b_id,
            pairing_type=pairing_type,
        )
        self.project.pairings.append(pairing)
        return pairing

    def remove_pairing(self, pairing_id: str) -> None:
        """Remove an existing pairing by ID."""
        before = len(self.project.pairings)
        self.project.pairings = [p for p in self.project.pairings if p.id != pairing_id]
        if len(self.project.pairings) == before:
            raise ValueError(f"Pairing '{pairing_id}' not found.")

    def get_preset_map(self) -> Dict[str, str]:
        return build_preset_map(self.project)

    def export_fasta(self) -> Tuple[str, List[str]]:
        return format_project_fasta(self.project)

    def save_fasta(self, path: str) -> None:
        save_fasta_file(path, self.project)

    def save_ab(self, path: str) -> None:
        save_ab_file(path, self.project)

    def load_ab(self, path: str) -> AntibodyProject:
        loaded = load_ab_file(path)
        try:
            self.project = loaded
            self._reset_selection()
        finally:
            self.clear_clipboard()
        return self.project

    def save_png(self, path: str, scale: float = 2.0) -> None:
        save_schematic_png(path, self.project, scale=scale)

    def _wipe_os_clipboard_if_matching(self, target_text: str) -> None:
        if not target_text:
            return
        try:
            if pyperclip.paste() == target_text:
                pyperclip.copy("")
                logger.info("Clipboard automatically cleared for security.")
        except Exception as e:
            config.log_sanitized_error(logger, "Error clearing clipboard", e)

    def clear_clipboard(
        self, content_to_clear: Optional[str] = None, gen: Optional[int] = None
    ) -> None:
        """Thread-safe conditional clipboard wipe honoring generation tokens and matching content."""
        with self._clipboard_lock:
            if gen is not None and gen != self._clipboard_gen:
                return
            if content_to_clear is not None and self.copied_fasta != content_to_clear:
                return

            if gen is None:
                self._clipboard_gen += 1
                if self._clipboard_timer is not None:
                    self._clipboard_timer.cancel()
                    self._clipboard_timer = None
            else:
                self._clipboard_timer = None

            target = content_to_clear if content_to_clear is not None else self.copied_fasta
            try:
                if target:
                    self._wipe_os_clipboard_if_matching(target)
            finally:
                self.copied_fasta = ""

    def copy_text_secure(self, text: str, label: str = "FASTA") -> bool:
        """Copy sequence text to OS clipboard and arm a 30-second auto-clear timer."""
        if not isinstance(text, str) or not text:
            return False

        with self._clipboard_lock:
            if self._clipboard_timer is not None:
                self._clipboard_timer.cancel()
                self._clipboard_timer = None
            self._clipboard_gen += 1
            prior_text = self.copied_fasta
            if prior_text:
                self._wipe_os_clipboard_if_matching(prior_text)
            self.copied_fasta = ""

            try:
                pyperclip.copy(text)
                self.copied_fasta = text
                gen = self._clipboard_gen
                timer = threading.Timer(
                    30.0, self.clear_clipboard, args=(text, gen)
                )
                timer.daemon = True
                timer.start()
                self._clipboard_timer = timer
                logger.info("Audit: Copied %s data to OS clipboard.", label)
                return True
            except Exception as e:
                config.log_sanitized_error(logger, "Error copying to clipboard", e)
                self.copied_fasta = ""
                return False

    def shutdown(self) -> None:
        """Cancel active timers and purge any copied sequence data on exit."""
        with self._clipboard_lock:
            self._clipboard_gen += 1
            if self._clipboard_timer is not None:
                self._clipboard_timer.cancel()
                self._clipboard_timer = None
            self.clear_clipboard()


def _build_canvas_shapes(scene: SchematicScene) -> List[Any]:
    shapes: List[Any] = []
    for pedge in scene.pairing_edges:
        if pedge.chain_a_id == pedge.chain_b_id and pedge.pairing_type == "VH_VL":
            mid_x = (pedge.x1 + pedge.x2) / 2.0
            ctrl_y = min(pedge.y1, pedge.y2) - 30.0
            shapes.append(
                cv.Path(
                    elements=[
                        cv.Path.MoveTo(pedge.x1, pedge.y1),
                        cv.Path.QuadraticTo(mid_x, ctrl_y, pedge.x2, pedge.y2),
                    ],
                    paint=ft.Paint(
                        color="#0E7490",
                        stroke_width=2,
                        style=ft.PaintingStyle.STROKE,
                        stroke_dash_pattern=[6, 4],
                    ),
                )
            )
        elif pedge.pairing_type == "HINGE_DISULFIDE":
            dx = pedge.x2 - pedge.x1
            dy = pedge.y2 - pedge.y1
            dist = math.hypot(dx, dy) or 1.0
            nx, ny = (-dy / dist) * 3.0, (dx / dist) * 3.0
            for sign in (-1.0, 1.0):
                shapes.append(
                    cv.Line(
                        pedge.x1 + nx * sign,
                        pedge.y1 + ny * sign,
                        pedge.x2 + nx * sign,
                        pedge.y2 + ny * sign,
                        paint=ft.Paint(color="#B45309", stroke_width=2),
                    )
                )
        else:
            shapes.append(
                cv.Line(
                    pedge.x1,
                    pedge.y1,
                    pedge.x2,
                    pedge.y2,
                    paint=ft.Paint(
                        color="#64748B",
                        stroke_width=2,
                        stroke_dash_pattern=[6, 4],
                    ),
                )
            )

    for edge in scene.edges:
        shapes.append(
            cv.Line(
                edge.x1,
                edge.y1,
                edge.x2,
                edge.y2,
                paint=ft.Paint(color="#475569", stroke_width=2.5),
            )
        )
    return shapes


def build_studio_ui(
    page: ft.Page, controller: Optional[AntibodyStudioController] = None
) -> AntibodyStudioController:
    """Build the 3-pane Flet Antibody & Bispecific Domain Engineering Studio."""
    ctrl = controller if controller is not None else AntibodyStudioController()
    page.title = "Fc Engineering Studio"
    page.padding = 10
    page.theme_mode = ft.ThemeMode.LIGHT
    page.window.prevent_close = True

    file_picker = ft.FilePicker()
    if hasattr(page, "services") and file_picker not in page.services:
        page.services.append(file_picker)

    def notify(msg: str, is_error: bool = False) -> None:
        snack = ft.SnackBar(
            content=ft.Text(msg),
            bgcolor="#DC2626" if is_error else "#0F766E",
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(snack)

    async def on_window_event(e: ft.WindowEvent) -> None:
        if e.type != ft.WindowEventType.CLOSE:
            return
        try:
            ctrl.shutdown()
        finally:
            await page.window.destroy()

    if hasattr(page, "window") and page.window is not None:
        page.window.on_event = on_window_event

    # Dynamic containers updated via render_workspace()
    preset_filter_state = {"query": ""}
    pairing_state: Dict[str, Optional[str]] = {
        "a": None,
        "b": None,
        "type": "VH_VL",
    }

    canvas_bg = cv.Canvas(shapes=[], expand=True)
    canvas_overlay = ft.Stack(controls=[], width=1100, height=560)
    chain_strip_col = ft.Column(spacing=6, scroll=ft.ScrollMode.AUTO)
    pairing_controls_col = ft.Column(spacing=6)
    inspector_col = ft.Column(spacing=8, scroll=ft.ScrollMode.AUTO, expand=True)

    project_name_tf = ft.TextField(
        label="Project Name",
        value=ctrl.project.name,
        width=220,
        dense=True,
        max_length=MAX_METADATA_LENGTH,
        input_filter=ft.InputFilter(allow=False, regex_string=r"[\r\n]"),
    )

    format_dd = ft.Dropdown(
        label="Format Template",
        value=ctrl.project.format_type,
        width=185,
        dense=True,
        options=[ft.dropdown.Option(k, text=lbl) for k, lbl in _FORMAT_LABELS],
    )

    def handle_drop_accept(
        target_chain_id: str, target_index: Optional[int], e: ft.DragTargetEvent
    ) -> None:
        if getattr(e, "src_id", None) is None:
            return
        src = page.get_control(e.src_id)
        if src is None:
            return
        payload = getattr(src, "data", None)
        if not isinstance(payload, str):
            return
        try:
            if payload.startswith("palette:"):
                dtype = payload.split(":", 1)[1]
                ctrl.add_domain_to_chain(target_chain_id, dtype, index=target_index)
            elif payload.startswith("chip:"):
                parts = payload.split(":", 2)
                if len(parts) == 3:
                    ctrl.move_domain_across_chains(
                        parts[1], parts[2], target_chain_id, target_index=target_index
                    )
            render_workspace()
        except Exception as err:
            notify(str(err), is_error=True)

    def make_slot_accept(chain_id: str, slot_index: Optional[int]):
        def _on_accept(e: ft.DragTargetEvent) -> None:
            handle_drop_accept(chain_id, slot_index, e)

        return _on_accept

    def rebuild_canvas_only() -> None:
        scene = build_schematic_scene(ctrl.project, ctrl.get_preset_map())
        canvas_bg.shapes = _build_canvas_shapes(scene)
        max_x, max_y = scene.bounds[2], scene.bounds[3]
        canvas_w = max(860.0, max_x + 48.0)
        canvas_h = max(420.0, max_y + 48.0)
        canvas_overlay.width = canvas_w
        canvas_overlay.height = canvas_h

        controls: List[ft.Control] = [
            ft.Container(content=canvas_bg, width=canvas_w, height=canvas_h)
        ]

        for node in scene.nodes:
            is_sel = (
                ctrl.selected_chain_id == node.chain_id
                and ctrl.selected_domain_id == node.domain_id
            )
            border_color = "#FACC15" if is_sel else ("#EF4444" if node.error_flag else "#0F172A")
            border_width = 2.5 if is_sel else 1.2

            dom_lookup: Optional[DomainBlock] = None
            for c in ctrl.project.chains:
                if c.id == node.chain_id:
                    for d in c.domains:
                        if d.id == node.domain_id:
                            dom_lookup = d
                            break

            def make_tap_handler(cid: str, did: str):
                def _on_tap(e: ft.ControlEvent) -> None:
                    ctrl.select_domain(cid, did)
                    render_workspace()

                return _on_tap

            def make_pan_handler(cid: str, did: str):
                def _on_pan_update(e: ft.DragUpdateEvent) -> None:
                    if e.local_delta is None:
                        return
                    for ch in ctrl.project.chains:
                        if ch.id == cid:
                            for dm in ch.domains:
                                if dm.id == did:
                                    ctrl.move_domain(
                                        cid,
                                        did,
                                        dm.x + e.local_delta.x,
                                        dm.y + e.local_delta.y,
                                    )
                                    rebuild_canvas_only()
                                    page.update()
                                    return

                return _on_pan_update

            box_col: List[ft.Control] = []
            if node.is_n_terminus and node.chain_name:
                box_col.append(
                    ft.Text(
                        node.chain_name,
                        size=10,
                        weight=ft.FontWeight.BOLD,
                        color="#0F172A",
                    )
                )
            box_col.append(
                ft.Container(
                    width=node.width,
                    height=node.height,
                    bgcolor=node.fill_color,
                    border=ft.Border.all(border_width, border_color),
                    border_radius=6,
                    alignment=ft.Alignment.CENTER,
                    content=ft.Text(
                        node.label,
                        color="#FFFFFF",
                        size=11,
                        weight=ft.FontWeight.BOLD,
                    ),
                )
            )
            if node.badge_text:
                box_col.append(
                    ft.Text(
                        node.badge_text,
                        size=9,
                        color="#B91C1C",
                        max_lines=1,
                    )
                )

            top_offset = max(0.0, node.y - (14.0 if node.is_n_terminus else 0.0))
            if dom_lookup is not None:
                controls.append(
                    ft.Container(
                        left=max(0.0, node.x),
                        top=top_offset,
                        content=ft.GestureDetector(
                            on_tap=make_tap_handler(node.chain_id, node.domain_id),
                            on_pan_update=make_pan_handler(
                                node.chain_id, node.domain_id
                            ),
                            content=ft.Column(spacing=1, controls=box_col),
                        ),
                    )
                )

        canvas_overlay.controls = controls

    def rebuild_chain_strip() -> None:
        rows: List[ft.Control] = []
        for chain in ctrl.project.chains:
            chip_controls: List[ft.Control] = []
            for idx, dom in enumerate(chain.domains):
                fill_col = DOMAIN_STYLE_MAP[dom.domain_type][2]
                is_sel = (
                    ctrl.selected_chain_id == chain.id
                    and ctrl.selected_domain_id == dom.id
                )

                def make_chip_click(cid: str, did: str):
                    def _on_click(e: ft.ControlEvent) -> None:
                        ctrl.select_domain(cid, did)
                        render_workspace()

                    return _on_click

                chip_badge = ft.Container(
                    padding=ft.Padding.symmetric(horizontal=8, vertical=4),
                    bgcolor=fill_col if is_sel else "#E2E8F0",
                    border=ft.Border.all(2.0 if is_sel else 1.2, "#FACC15" if is_sel else fill_col),
                    border_radius=6,
                    on_click=make_chip_click(chain.id, dom.id),
                    content=ft.Text(
                        dom.domain_type,
                        size=11,
                        weight=ft.FontWeight.W_600,
                        color="#FFFFFF" if is_sel else "#0F172A",
                    ),
                )
                draggable_chip = ft.Draggable(
                    group="domain_palette",
                    data=f"chip:{chain.id}:{dom.id}",
                    content=chip_badge,
                )
                slot_target = ft.DragTarget(
                    group="domain_palette",
                    on_accept=make_slot_accept(chain.id, idx),
                    content=draggable_chip,
                )
                chip_controls.append(slot_target)

            # Trailing drop target to append at chain end
            chip_controls.append(
                ft.DragTarget(
                    group="domain_palette",
                    on_accept=make_slot_accept(chain.id, None),
                    content=ft.Container(
                        padding=ft.Padding.symmetric(horizontal=8, vertical=4),
                        border=ft.Border.all(1, "#94A3B8"),
                        border_radius=6,
                        content=ft.Text("+ Drop", size=10, color="#475569"),
                    ),
                )
            )

            def make_rm_chain(cid: str):
                def _on_rm(e: ft.ControlEvent) -> None:
                    ctrl.remove_chain(cid)
                    render_workspace()

                return _on_rm

            rows.append(
                ft.Container(
                    padding=6,
                    bgcolor="#F8FAFC",
                    border=ft.Border.all(1, "#CBD5E1"),
                    border_radius=6,
                    content=ft.Row(
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        controls=[
                            ft.Row(
                                spacing=6,
                                wrap=True,
                                expand=True,
                                controls=[
                                    ft.Text(
                                        f"{chain.name}:",
                                        weight=ft.FontWeight.BOLD,
                                        size=12,
                                    ),
                                    *chip_controls,
                                ],
                            ),
                            ft.IconButton(
                                icon=ft.Icons.DELETE_OUTLINE,
                                tooltip="Remove Chain",
                                on_click=make_rm_chain(chain.id),
                            ),
                        ],
                    ),
                )
            )
        chain_strip_col.controls = rows

    def rebuild_pairings() -> None:
        endpoint_options: List[ft.dropdown.Option] = []
        for chain in ctrl.project.chains:
            for dom in chain.domains:
                key = f"{chain.id}:{dom.id}"
                endpoint_options.append(
                    ft.dropdown.Option(key, text=f"{chain.name}.{dom.domain_type}")
                )

        valid_keys = {opt.key for opt in endpoint_options}
        if pairing_state["a"] not in valid_keys:
            pairing_state["a"] = endpoint_options[0].key if endpoint_options else None
        if pairing_state["b"] not in valid_keys:
            pairing_state["b"] = (
                endpoint_options[1].key
                if len(endpoint_options) > 1
                else pairing_state["a"]
            )

        dd_a = ft.Dropdown(
            label="Domain A",
            value=pairing_state["a"],
            options=endpoint_options,
            dense=True,
            on_select=lambda e: pairing_state.__setitem__("a", e.control.value),
        )
        dd_b = ft.Dropdown(
            label="Domain B",
            value=pairing_state["b"],
            options=endpoint_options,
            dense=True,
            on_select=lambda e: pairing_state.__setitem__("b", e.control.value),
        )
        dd_ptype = ft.Dropdown(
            label="Pairing Type",
            value=pairing_state["type"],
            options=[ft.dropdown.Option(pt) for pt in VALID_PAIRING_TYPES],
            dense=True,
            on_select=lambda e: pairing_state.__setitem__("type", e.control.value),
        )

        def on_add_pairing(e: ft.ControlEvent) -> None:
            a_val = pairing_state["a"]
            b_val = pairing_state["b"]
            pt_val = pairing_state["type"] or "VH_VL"
            if not a_val or not b_val:
                notify("Select two domains to pair.", is_error=True)
                return
            try:
                ca, da = a_val.split(":", 1)
                cb, db = b_val.split(":", 1)
                ctrl.add_pairing(ca, da, cb, db, pt_val)
                render_workspace()
            except Exception as err:
                notify(str(err), is_error=True)

        pairing_list_controls: List[ft.Control] = []
        for p in ctrl.project.pairings:
            def make_del_p(pid: str):
                def _on_del(e: ft.ControlEvent) -> None:
                    ctrl.remove_pairing(pid)
                    render_workspace()

                return _on_del

            pairing_list_controls.append(
                ft.Row(
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    controls=[
                        ft.Text(
                            f"{p.pairing_type}: {p.chain_a_id}.{p.domain_a_id} ↔ {p.chain_b_id}.{p.domain_b_id}",
                            size=10,
                            expand=True,
                        ),
                        ft.IconButton(
                            icon=ft.Icons.CLOSE,
                            icon_size=14,
                            on_click=make_del_p(p.id),
                        ),
                    ],
                )
            )

        pairing_controls_col.controls = [
            dd_a,
            dd_b,
            dd_ptype,
            ft.Button("Add Pairing", on_click=on_add_pairing),
            *pairing_list_controls,
        ]

    def rebuild_inspector() -> None:
        controls: List[ft.Control] = []
        if not ctrl.selected_chain_id:
            inspector_col.controls = [
                ft.Text("Select a chain or domain to inspect.", size=12)
            ]
            return

        chain = ctrl._get_chain(ctrl.selected_chain_id)
        dom: Optional[DomainBlock] = None
        if ctrl.selected_domain_id:
            for d in chain.domains:
                if d.id == ctrl.selected_domain_id:
                    dom = d
                    break

        controls.append(
            ft.Text(
                f"Chain: {chain.name}",
                size=14,
                weight=ft.FontWeight.BOLD,
            )
        )

        if dom is not None:
            controls.append(
                ft.Row(
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    controls=[
                        ft.Text(
                            f"Selected Domain: {dom.domain_type} ({dom.id})",
                            size=12,
                            weight=ft.FontWeight.W_600,
                        ),
                        ft.OutlinedButton(
                            "Delete Domain",
                            on_click=lambda e: (
                                ctrl.remove_domain(chain.id, dom.id),
                                render_workspace(),
                            ),
                        ),
                    ],
                )
            )
            if dom.domain_type in CONSTANT_DOMAIN_EU_RANGES:
                iso_options = [
                    ft.dropdown.Option(iso_k, text=iso_k.upper())
                    for iso_k in config.SEQUENCES.keys()
                ]
                allo_map = config.SEQUENCES.get(dom.isotype, {})
                allo_options = [ft.dropdown.Option(ak) for ak in allo_map.keys()]
                current_allo = dom.allotype or (
                    next(iter(allo_map)) if allo_map else ""
                )

                def on_iso_select(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id, dom.id, isotype=e.control.value
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    render_workspace()

                def on_allo_select(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id, dom.id, allotype=e.control.value
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    render_workspace()

                controls.append(
                    ft.Row(
                        controls=[
                            ft.Dropdown(
                                label="Isotype",
                                value=dom.isotype,
                                options=iso_options,
                                dense=True,
                                width=120,
                                on_select=on_iso_select,
                            ),
                            ft.Dropdown(
                                label="Allotype",
                                value=current_allo,
                                options=allo_options,
                                dense=True,
                                expand=True,
                                on_select=on_allo_select,
                            ),
                        ]
                    )
                )
            elif dom.domain_type == "CL":
                tkeys = ALLOWED_TEMPLATE_KEYS_BY_TYPE["CL"]

                def on_cl_select(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id, dom.id, template_key=e.control.value
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    render_workspace()

                controls.append(
                    ft.Dropdown(
                        label="CL Subtype",
                        value=dom.template_key or "CL_kappa",
                        options=[ft.dropdown.Option(k) for k in tkeys],
                        dense=True,
                        on_select=on_cl_select,
                    )
                )
            elif dom.domain_type in ("VH", "VL", "G4S"):
                tkeys = ALLOWED_TEMPLATE_KEYS_BY_TYPE[dom.domain_type]
                seq_tf = ft.TextField(
                    label=f"Custom {dom.domain_type} Sequence",
                    value=dom.custom_sequence,
                    multiline=True,
                    min_lines=2,
                    max_lines=3,
                    max_length=MAX_CUSTOM_SEQ_LENGTH,
                    dense=True,
                )

                def on_tkey_select(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id, dom.id, template_key=e.control.value
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    render_workspace()

                def on_apply_cseq(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id,
                        dom.id,
                        custom_sequence=(seq_tf.value or "").strip().upper(),
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    else:
                        notify(f"Updated {dom.domain_type} sequence.")
                    render_workspace()

                controls.extend(
                    [
                        ft.Dropdown(
                            label="Default Template",
                            value=dom.template_key or DEFAULT_TEMPLATE_KEYS[dom.domain_type],
                            options=[ft.dropdown.Option(k) for k in tkeys],
                            dense=True,
                            on_select=on_tkey_select,
                        ),
                        seq_tf,
                        ft.Button("Apply Sequence", on_click=on_apply_cseq),
                    ]
                )

        # Fc Mutation Selector & Apply-to-All-Heavy Toggle
        controls.append(ft.Divider())
        controls.append(
            ft.Switch(
                label="Apply to all heavy chains",
                value=ctrl.apply_to_all_heavy,
                on_change=lambda e: setattr(
                    ctrl, "apply_to_all_heavy", bool(e.control.value)
                ),
            )
        )

        active_presets_on_chain = {
            p for d in chain.domains for p in d.applied_presets
        }
        preset_map = ctrl.get_preset_map()
        preset_list_view = ft.ListView(controls=[], spacing=2)

        def populate_preset_boxes() -> None:
            q = preset_filter_state["query"].strip().lower()
            boxes: List[ft.Control] = []
            for label, val in preset_map.items():
                if q and q not in label.lower() and q not in val.lower():
                    continue

                def make_preset_cb(lbl: str):
                    def _on_cb(e: ft.ControlEvent) -> None:
                        errs = ctrl.toggle_preset(chain.id, lbl, bool(e.control.value))
                        if errs:
                            notify("; ".join(errs), is_error=True)
                        render_workspace()

                    return _on_cb

                boxes.append(
                    ft.Checkbox(
                        label=f"{label} ({val})",
                        value=(label in active_presets_on_chain),
                        on_change=make_preset_cb(label),
                    )
                )
            preset_list_view.controls = boxes

        def on_filter_change(e: ft.ControlEvent) -> None:
            preset_filter_state["query"] = e.control.value or ""
            populate_preset_boxes()
            page.update()

        search_tf = ft.TextField(
            label="Filter Fc Presets",
            value=preset_filter_state["query"],
            dense=True,
            on_change=on_filter_change,
        )
        controls.append(search_tf)
        populate_preset_boxes()

        controls.append(
            ft.Container(
                height=160,
                border=ft.Border.all(1, "#E2E8F0"),
                border_radius=6,
                padding=4,
                content=preset_list_view,
            )
        )

        # Ad-hoc EU Mutation input
        custom_eu_tf = ft.TextField(
            label="Apply Ad-hoc EU Mutation (e.g. M252Y/S254T)",
            dense=True,
            max_length=1000,
            input_filter=ft.InputFilter(allow=True, regex_string=r"[a-zA-Z0-9/, ]"),
        )

        def on_apply_adhoc(e: ft.ControlEvent) -> None:
            errs = ctrl.apply_custom_mutation(chain.id, custom_eu_tf.value or "")
            if errs:
                notify("; ".join(errs), is_error=True)
            else:
                notify("Applied custom EU mutation.")
                render_workspace()

        controls.append(
            ft.Row(
                controls=[
                    ft.Container(expand=True, content=custom_eu_tf),
                    ft.Button("Apply", on_click=on_apply_adhoc),
                ]
            )
        )

        # File-Scoped Custom Mutation registration
        file_mut_label_tf = ft.TextField(
            label="New File Preset Label",
            dense=True,
            max_length=40,
            input_filter=ft.InputFilter(allow=False, regex_string=r"[\r\n]"),
        )
        file_mut_val_tf = ft.TextField(
            label="Preset EU Mutations (e.g. L234A/L235A)",
            dense=True,
            max_length=1000,
            input_filter=ft.InputFilter(allow=True, regex_string=r"[a-zA-Z0-9/, ]"),
        )

        def on_add_file_mut(e: ft.ControlEvent) -> None:
            errs = ctrl.register_file_custom_mutation(
                file_mut_label_tf.value or "", file_mut_val_tf.value or ""
            )
            if errs:
                notify("; ".join(errs), is_error=True)
            else:
                notify("Added custom preset to current project file.")
                render_workspace()

        controls.extend(
            [
                ft.Text(
                    "Add Custom Mutation to File",
                    size=12,
                    weight=ft.FontWeight.W_600,
                ),
                file_mut_label_tf,
                file_mut_val_tf,
                ft.OutlinedButton("Save to Project File", on_click=on_add_file_mut),
            ]
        )

        # Live Sequence, Diff & Developability Summary
        controls.append(ft.Divider())
        chain_seq, _, seq_errs = assemble_chain_sequence(chain, preset_map)
        diffs, diff_errs = diff_chain_sequences(chain, preset_map)
        liabilities, dev_errs = analyze_chain_developability(chain, preset_map)

        if seq_errs or diff_errs or dev_errs:
            all_errs = seq_errs + diff_errs + dev_errs
            controls.append(
                ft.Text(
                    "Resolution Errors: " + "; ".join(all_errs),
                    size=11,
                    color="#DC2626",
                )
            )
        else:
            diff_text = (
                ", ".join(f"{dt}:{wt}{pos}{mut}" for dt, pos, wt, mut in diffs)
                or "None (WT)"
            )
            dev_lines = [
                f"{cat}: {', '.join(str(p) for p in pos_list) or 'None'}"
                for cat, pos_list in liabilities
            ]
            controls.extend(
                [
                    ft.Text(
                        f"Chain Length: {len(chain_seq)} aa",
                        size=12,
                        weight=ft.FontWeight.BOLD,
                    ),
                    ft.Text(f"EU Diffs: {diff_text}", size=11),
                    ft.Text(
                        "Developability:\n" + "\n".join(dev_lines),
                        size=10,
                    ),
                    ft.TextField(
                        label="Assembled Chain Sequence",
                        value=chain_seq,
                        multiline=True,
                        read_only=True,
                        min_lines=3,
                        max_lines=4,
                        text_size=10,
                    ),
                ]
            )

        inspector_col.controls = controls

    def render_workspace() -> None:
        project_name_tf.value = ctrl.project.name
        format_dd.value = ctrl.project.format_type
        rebuild_canvas_only()
        rebuild_chain_strip()
        rebuild_pairings()
        rebuild_inspector()
        page.update()

    # Toolbar callbacks
    def on_project_name_commit(e: ft.ControlEvent) -> None:
        try:
            ctrl.rename_project(project_name_tf.value or "")
        except ValueError as err:
            project_name_tf.value = ctrl.project.name
            notify(str(err), is_error=True)
            page.update()

    project_name_tf.on_submit = on_project_name_commit
    project_name_tf.on_blur = on_project_name_commit

    def on_format_select(e: ft.ControlEvent) -> None:
        if e.control.value:
            ctrl.select_format(e.control.value)
            render_workspace()

    format_dd.on_select = on_format_select

    async def on_load_ab(e: ft.ControlEvent) -> None:
        files = await file_picker.pick_files(
            dialog_title="Load Antibody Project (.ab)",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["ab"],
        )
        if not files or not files[0].path:
            return
        try:
            ctrl.load_ab(files[0].path)
            notify(f"Loaded {os.path.basename(files[0].path)}")
            render_workspace()
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_save_ab(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Save Antibody Project (.ab)",
            file_name="antibody_project.ab",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["ab"],
        )
        if not dest:
            return
        try:
            ctrl.save_ab(dest)
            notify(f"Saved {os.path.basename(dest)}")
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_export_fasta(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Export Multi-FASTA",
            file_name="antibody_chains.fasta",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["fasta", "fa"],
        )
        if not dest:
            return
        try:
            ctrl.save_fasta(dest)
            notify(f"Exported {os.path.basename(dest)}")
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_export_png(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Export Transparent PNG Schematic",
            file_name="antibody_schematic.png",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["png"],
        )
        if not dest:
            return
        try:
            ctrl.save_png(dest, scale=2.0)
            notify(f"Exported {os.path.basename(dest)}")
        except Exception as err:
            notify(str(err), is_error=True)

    def on_copy_fasta(e: ft.ControlEvent) -> None:
        fasta_text, errs = ctrl.export_fasta()
        if errs:
            notify("; ".join(errs), is_error=True)
            return
        if ctrl.copy_text_secure(fasta_text, label="FASTA"):
            notify("Multi-FASTA copied! (Will auto-clear in 30s)")
        else:
            notify("Error copying to clipboard.", is_error=True)

    palette_cards: List[ft.Control] = [
        ft.Draggable(
            group="domain_palette",
            data=f"palette:{dtype}",
            content=ft.Container(
                padding=ft.Padding.symmetric(horizontal=10, vertical=6),
                bgcolor="#E2E8F0",
                border=ft.Border.all(1, "#64748B"),
                border_radius=6,
                content=ft.Text(dtype, size=12, weight=ft.FontWeight.W_600),
            ),
        )
        for dtype in VALID_DOMAIN_TYPES
    ]

    toolbar = ft.Row(
        wrap=True,
        spacing=8,
        controls=[
            project_name_tf,
            format_dd,
            ft.Button(
                "Add Chain",
                on_click=lambda e: (ctrl.add_chain(), render_workspace()),
            ),
            ft.OutlinedButton("Load .ab", on_click=on_load_ab),
            ft.OutlinedButton("Save .ab", on_click=on_save_ab),
            ft.OutlinedButton("Export FASTA", on_click=on_export_fasta),
            ft.OutlinedButton("Export PNG", on_click=on_export_png),
            ft.Button(
                "Copy FASTA (30s Auto-Clear)",
                on_click=on_copy_fasta,
            ),
        ],
    )

    left_panel = ft.Container(
        width=240,
        padding=8,
        border=ft.Border.all(1, "#CBD5E1"),
        border_radius=8,
        content=ft.Column(
            scroll=ft.ScrollMode.AUTO,
            spacing=8,
            controls=[
                ft.Text("Domain Palette", size=14, weight=ft.FontWeight.BOLD),
                ft.Row(wrap=True, spacing=6, controls=palette_cards),
                ft.Divider(),
                ft.Text("Pairings", size=14, weight=ft.FontWeight.BOLD),
                pairing_controls_col,
            ],
        ),
    )

    center_panel = ft.Container(
        expand=True,
        padding=8,
        border=ft.Border.all(1, "#CBD5E1"),
        border_radius=8,
        content=ft.Column(
            expand=True,
            spacing=8,
            controls=[
                ft.Text("2D Interactive Schematic", size=14, weight=ft.FontWeight.BOLD),
                ft.Container(
                    expand=True,
                    bgcolor="#F8FAFC",
                    border=ft.Border.all(1, "#E2E8F0"),
                    border_radius=6,
                    content=ft.Column(
                        scroll=ft.ScrollMode.AUTO,
                        controls=[
                            ft.Row(
                                scroll=ft.ScrollMode.AUTO,
                                controls=[canvas_overlay],
                            )
                        ],
                    ),
                ),
                ft.Text(
                    "Polypeptide Chains (N -> C Drop Targets)",
                    size=13,
                    weight=ft.FontWeight.BOLD,
                ),
                ft.Container(height=190, content=chain_strip_col),
            ],
        ),
    )

    right_panel = ft.Container(
        width=360,
        padding=8,
        border=ft.Border.all(1, "#CBD5E1"),
        border_radius=8,
        content=inspector_col,
    )

    page.add(
        ft.Column(
            expand=True,
            spacing=8,
            controls=[
                toolbar,
                ft.Row(
                    expand=True,
                    spacing=8,
                    vertical_alignment=ft.CrossAxisAlignment.STRETCH,
                    controls=[left_panel, center_panel, right_panel],
                ),
            ],
        )
    )

    render_workspace()
    return ctrl


def main() -> None:
    """Launch the Fc Engineering Studio Flet desktop application."""
    controller = AntibodyStudioController()
    atexit.register(controller.shutdown)
    try:
        ft.run(lambda page: build_studio_ui(page, controller))
    finally:
        try:
            controller.shutdown()
        finally:
            atexit.unregister(controller.shutdown)


if __name__ == "__main__":
    main()
