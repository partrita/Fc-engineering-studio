"""Flet local desktop GUI and state controller for Fc Engineering Studio."""

import asyncio
import atexit
import copy
import logging
import math
import os
import tempfile
import threading
import time
import webbrowser
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import flet as ft
import flet.canvas as cv
import pyperclip

from fc_engineer import config
from fc_engineer.core import (
    CONSTANT_DOMAIN_EU_RANGES,
    compute_protein_metrics,
    generate_batch,
    lookup_wt_residue,
)
from fc_engineer.designer import (
    ALLOWED_TEMPLATE_KEYS_BY_TYPE,
    DEFAULT_SCHEMATIC_THEME,
    DEFAULT_TEMPLATE_KEYS,
    DOMAIN_STYLE_MAP,
    MAX_CHAINS,
    MAX_COORD,
    MAX_CUSTOM_SEQ_LENGTH,
    MAX_DOMAINS_PER_CHAIN,
    MAX_METADATA_LENGTH,
    MAX_PAIRINGS,
    VALID_DOMAIN_TYPES,
    VALID_ORIENTATIONS,
    VALID_PAIRING_TYPES,
    VALID_SHAPES,
    AntibodyProject,
    DomainBlock,
    DomainVisualStyle,
    InterChainPairing,
    PolypeptideChain,
    SchematicScene,
    SchematicTheme,
    add_file_custom_mutation,
    analyze_chain_developability,
    apply_heterodimer_pair,
    assemble_chain_sequence,
    auto_detect_standard_pairings,
    build_format_template,
    build_preset_map,
    build_schematic_scene,
    clone_chain_in_project,
    compute_hexagon_vertices,
    delete_file_custom_mutation,
    diff_chain_sequences,
    reflow_chain_coordinates,
    reflow_project_coordinates,
    remove_chain_from_project,
    remove_custom_mutation_from_domain,
    remove_domain_from_chain,
    remove_preset_from_chain,
    reset_chain_mutations,
    resolve_domain_sequence,
    route_mutations_to_chain,
    swap_crossmab_ch1_cl,
)
from fc_engineer.exporters import (
    format_project_fasta,
    load_ab_file,
    save_ab_file,
    save_csv_report_file,
    save_fasta_file,
    save_genbank_file,
    save_schematic_png,
    serialize_ab_project,
)

logger = logging.getLogger(__name__)

_FORMAT_LABELS: List[Tuple[str, str]] = [
    ("igg", "Monospecific IgG"),
    ("crossmab", "CrossMab"),
    ("igg_scfv", "IgG-scFv"),
    ("scfv_fc", "scFv-Fc"),
    ("dvd_ig", "DVD-Ig"),
    ("bite", "BiTE"),
    ("two_plus_one", "2+1 T-Cell Engager"),
    ("trispecific", "Trispecific IgG-scFv"),
    ("fab_igg", "Fab-IgG"),
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
        self.theme_mode: str = "light"
        self.theme_presets: Dict[str, SchematicTheme] = config.load_theme_presets()
        self.theme: SchematicTheme = copy.deepcopy(
            self.theme_presets.get(DEFAULT_SCHEMATIC_THEME.name, DEFAULT_SCHEMATIC_THEME)
        )
        self.project: AntibodyProject = build_format_template(
            initial_format, isotype=isotype, allotype=allotype
        )
        reflow_project_coordinates(self.project, theme=self.theme)
        self._undo_stack: List[AntibodyProject] = []
        self._redo_stack: List[AntibodyProject] = []
        self.is_dirty: bool = False
        self.current_file_path: Optional[str] = None
        self.recent_files: List[str] = []
        self.zoom_scale: float = 1.0
        self.snap_to_grid: bool = False
        self.include_variable_dev: bool = False
        self.side_panel_width: float = 360.0
        self.wiz_state: Dict[str, Optional[str]] = {"a": None, "b": None, "mode": "kih_ss"}
        self.clipboard_expires_in: int = 0
        self._clipboard_deadline: float = 0.0
        self.auto_recovery: bool = False
        self._shutdown_flag: bool = False
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
            self._record_undo()
            self.project = new_proj
            self.current_file_path = None
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
        self._record_undo()
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
            x=40.0 + len(self.project.chains) * self.theme.chain_spacing
            if self.theme.orientation == "vertical"
            else 40.0,
            y=40.0
            if self.theme.orientation == "vertical"
            else 40.0 + len(self.project.chains) * self.theme.chain_spacing,
        )
        self._record_undo()
        self.project.chains.append(chain)
        self.selected_chain_id = chain.id
        self.selected_domain_id = None
        return chain

    def remove_chain(self, chain_id: str) -> None:
        """Remove a polypeptide chain and update current selection safely."""
        self._get_chain(chain_id)
        self._record_undo()
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
        self._record_undo()
        if index is None:
            chain.domains.append(dom)
        else:
            clamped = max(0, min(len(chain.domains), index))
            chain.domains.insert(clamped, dom)

        reflow_chain_coordinates(chain, theme=self.theme)
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
        self._record_undo()
        dom = chain.domains.pop(from_index)
        chain.domains.insert(to_index, dom)
        reflow_chain_coordinates(chain, theme=self.theme)
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
            self._record_undo()
            src_chain.domains.pop(src_idx)
            if target_index is None:
                src_chain.domains.append(dom)
            else:
                adj = target_index - 1 if target_index > src_idx else target_index
                clamped = max(0, min(len(src_chain.domains), adj))
                src_chain.domains.insert(clamped, dom)
            reflow_chain_coordinates(src_chain, theme=self.theme)
            self.selected_chain_id = src_chain.id
            self.selected_domain_id = dom.id
            return

        if len(dst_chain.domains) >= MAX_DOMAINS_PER_CHAIN:
            raise ValueError(
                f"Target chain '{target_chain_id}' reached maximum of {MAX_DOMAINS_PER_CHAIN} domains."
            )

        self._record_undo()
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

        reflow_chain_coordinates(src_chain, theme=self.theme)
        reflow_chain_coordinates(dst_chain, theme=self.theme)
        self.selected_chain_id = dst_chain.id
        self.selected_domain_id = dom.id

    def remove_domain(self, chain_id: str, domain_id: str) -> None:
        """Delete a domain from a chain and clean up any orphaned pairings."""
        chain, _, _ = self._get_domain(chain_id, domain_id)
        self._record_undo()
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
        dom.x = max(0.0, min(MAX_COORD, fx))
        dom.y = max(0.0, min(MAX_COORD, fy))
        dom.pinned = True

    def unpin_domain(self, chain_id: str, domain_id: str) -> None:
        """Unpin a domain and reflow its chain coordinates."""
        chain, dom, _ = self._get_domain(chain_id, domain_id)
        self._record_undo()
        dom.pinned = False
        reflow_chain_coordinates(chain, theme=self.theme)

    def reset_layout(self) -> None:
        """Unpin all domains across all chains and reflow the entire schematic."""
        self._record_undo()
        for chain in self.project.chains:
            for dom in chain.domains:
                dom.pinned = False
        reflow_project_coordinates(self.project, theme=self.theme)

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
        pre_project = copy.deepcopy(self.project)
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
        self._push_undo(pre_project)
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
            self._record_undo()
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

        self._record_undo()
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

        self._record_undo()
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
        record_undo: bool = True,
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
        if record_undo:
            self._record_undo()
        self.project.pairings.append(pairing)
        return pairing

    def remove_pairing(self, pairing_id: str) -> None:
        """Remove an existing pairing by ID."""
        if not any(p.id == pairing_id for p in self.project.pairings):
            raise ValueError(f"Pairing '{pairing_id}' not found.")
        self._record_undo()
        self.project.pairings = [p for p in self.project.pairings if p.id != pairing_id]

    def get_preset_map(self) -> Dict[str, str]:
        return build_preset_map(self.project)

    def export_fasta(self) -> Tuple[str, List[str]]:
        return format_project_fasta(self.project)

    def save_fasta(self, path: str) -> None:
        save_fasta_file(path, self.project)

    @property
    def can_undo(self) -> bool:
        return bool(self._undo_stack)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo_stack)

    def _push_undo(self, snapshot: AntibodyProject) -> None:
        self._undo_stack.append(snapshot)
        if len(self._undo_stack) > 50:
            self._undo_stack = self._undo_stack[-50:]
        self._redo_stack.clear()
        self.is_dirty = True

    def _record_undo(self) -> None:
        self._push_undo(copy.deepcopy(self.project))

    def undo(self) -> bool:
        if not self._undo_stack:
            return False
        self._redo_stack.append(copy.deepcopy(self.project))
        self.project = self._undo_stack.pop()
        self.is_dirty = True
        self._reset_selection()
        return True

    def redo(self) -> bool:
        if not self._redo_stack:
            return False
        self._undo_stack.append(copy.deepcopy(self.project))
        self.project = self._redo_stack.pop()
        self.is_dirty = True
        self._reset_selection()
        return True

    def rename_chain(self, chain_id: str, new_name: str) -> None:
        chain = self._get_chain(chain_id)
        if not isinstance(new_name, str) or any(ord(ch) < 0x20 or ord(ch) == 0x7F for ch in new_name):
            raise ValueError("Chain name contains disallowed control characters.")
        stripped = new_name.strip()
        if not stripped or len(stripped) > MAX_METADATA_LENGTH:
            raise ValueError("Invalid chain name length.")
        self._record_undo()
        chain.name = stripped

    def reorder_chain(self, chain_id: str, delta: int) -> None:
        idx = next((i for i, c in enumerate(self.project.chains) if c.id == chain_id), None)
        if idx is None:
            raise ValueError(f"Chain '{chain_id}' not found.")
        target = idx + delta
        if not (0 <= target < len(self.project.chains)):
            return
        self._record_undo()
        ch = self.project.chains.pop(idx)
        self.project.chains.insert(target, ch)
        reflow_project_coordinates(self.project, theme=self.theme)

    def insert_macro_block(self, chain_id: str, macro_kind: str) -> List[DomainBlock]:
        chain = self._get_chain(chain_id)
        if macro_kind == "scfv":
            specs = [("VH", "VH_A"), ("G4S", "G4S_3"), ("VL", "VL_A")]
        elif macro_kind == "fc":
            specs = [("Hinge", ""), ("CH2", ""), ("CH3", "")]
        else:
            raise ValueError(f"Unsupported macro block '{macro_kind}'.")
        if len(chain.domains) + len(specs) > MAX_DOMAINS_PER_CHAIN:
            raise ValueError(f"Chain '{chain_id}' would exceed maximum of {MAX_DOMAINS_PER_CHAIN} domains.")
        self._record_undo()
        added: List[DomainBlock] = []
        iso_map = config.SEQUENCES.get("igg1", {})
        default_allo = next(iter(iso_map)) if iso_map else ""
        for dtype, tkey in specs:
            is_const = dtype in CONSTANT_DOMAIN_EU_RANGES
            dom = DomainBlock(
                id=self._next_domain_id(chain.id, dtype),
                domain_type=dtype,
                isotype="igg1",
                allotype=default_allo if is_const else "",
                template_key=tkey if not is_const else "",
            )
            chain.domains.append(dom)
            added.append(dom)
        if macro_kind == "scfv" and len(self.project.pairings) < MAX_PAIRINGS:
            self.add_pairing(
                chain.id, added[0].id, chain.id, added[2].id, "VH_VL", record_undo=False
            )
        reflow_chain_coordinates(chain, theme=self.theme)
        self.selected_chain_id = chain.id
        self.selected_domain_id = added[0].id
        return added

    def remove_custom_mutation(self, chain_id: str, domain_id: str, token: str) -> None:
        _, dom, _ = self._get_domain(chain_id, domain_id)
        self._record_undo()
        remove_custom_mutation_from_domain(dom, token)

    def delete_file_custom_mutation(self, label: str) -> None:
        self._record_undo()
        delete_file_custom_mutation(self.project, label)

    def reset_chain_mutations(self, chain_id: str) -> None:
        chain = self._get_chain(chain_id)
        self._record_undo()
        reset_chain_mutations(chain)

    def auto_detect_pairings(self) -> int:
        self._record_undo()
        return auto_detect_standard_pairings(self.project)

    def save_ab(self, path: str) -> None:
        save_ab_file(path, self.project)
        self.current_file_path = os.path.abspath(path)
        self.is_dirty = False
        if self.current_file_path not in self.recent_files:
            self.recent_files = [self.current_file_path] + self.recent_files[:4]

    def load_ab(self, path: str) -> AntibodyProject:
        loaded = load_ab_file(path)
        try:
            self.project = loaded
            self.current_file_path = os.path.abspath(path)
            self.is_dirty = False
            self._undo_stack.clear()
            self._redo_stack.clear()
            if self.current_file_path not in self.recent_files:
                self.recent_files = [self.current_file_path] + self.recent_files[:4]
            self._reset_selection()
        finally:
            self.clear_clipboard()
        return self.project

    def save_png(self, path: str, scale: float = 2.0, bg_color: Optional[str] = None) -> None:
        save_schematic_png(path, self.project, scale=scale, theme=self.theme, bg_color=bg_color)

    def save_genbank(self, path: str) -> None:
        save_genbank_file(path, self.project)

    def save_csv(self, path: str) -> None:
        save_csv_report_file(path, self.project)

    def apply_theme_preset(self, name: str) -> None:
        if name not in self.theme_presets:
            raise ValueError(f"Unknown theme preset '{name}'.")
        self.theme = copy.deepcopy(self.theme_presets[name])
        reflow_project_coordinates(self.project, theme=self.theme)

    def set_theme(self, theme: SchematicTheme) -> None:
        self.theme = copy.deepcopy(theme)
        reflow_project_coordinates(self.project, theme=self.theme)

    def set_spacing(
        self,
        domain_spacing: Optional[float] = None,
        chain_spacing: Optional[float] = None,
    ) -> None:
        if domain_spacing is not None:
            self.theme.domain_spacing = max(4.0, min(80.0, float(domain_spacing)))
        if chain_spacing is not None:
            self.theme.chain_spacing = max(60.0, min(260.0, float(chain_spacing)))
        reflow_project_coordinates(self.project, theme=self.theme)

    def load_theme(self, path: str) -> SchematicTheme:
        loaded = config.load_theme_file(path)
        self.theme = loaded
        reflow_project_coordinates(self.project, theme=self.theme)
        return loaded

    def save_theme(self, path: str) -> None:
        config.save_theme_file(path, self.theme)

    def set_theme_mode(self, mode: str) -> None:
        if mode not in ("dark", "light"):
            raise ValueError("theme_mode must be 'dark' or 'light'.")
        self.theme_mode = mode

    def clone_chain(self, chain_id: str) -> PolypeptideChain:
        pre = copy.deepcopy(self.project)
        cloned = clone_chain_in_project(self.project, chain_id, theme=self.theme)
        self._push_undo(pre)
        self.selected_chain_id = cloned.id
        self.selected_domain_id = cloned.domains[0].id if cloned.domains else None
        return cloned

    def swap_crossmab_ch1_cl(self, chain_a_id: str, chain_b_id: str) -> List[str]:
        pre = copy.deepcopy(self.project)
        errs = swap_crossmab_ch1_cl(
            self.project, chain_a_id, chain_b_id, theme=self.theme
        )
        if not errs:
            self._push_undo(pre)
        return errs

    def apply_heterodimer_pair(
        self, chain_a_id: str, chain_b_id: str, mode: str = "kih_ss"
    ) -> List[str]:
        pre = copy.deepcopy(self.project)
        errs = apply_heterodimer_pair(
            self.project, chain_a_id, chain_b_id, mode=mode
        )
        if not errs:
            self._push_undo(pre)
        return errs

    def recovery_snapshot_path(self) -> str:
        base_dir = os.path.dirname(config.default_prefs_path())
        return os.path.join(base_dir, "recovery.ab")

    def has_recovery_snapshot(self) -> bool:
        return os.path.isfile(self.recovery_snapshot_path())

    def save_recovery_snapshot(self) -> None:
        if not self.is_dirty:
            return
        non_empty = [c for c in self.project.chains if c.domains]
        if not non_empty:
            return
        filtered_proj = copy.deepcopy(self.project)
        filtered_proj.chains = non_empty
        valid_cids = {c.id for c in non_empty}
        filtered_proj.pairings = [
            p
            for p in filtered_proj.pairings
            if p.chain_a_id in valid_cids and p.chain_b_id in valid_cids
        ]
        path = self.recovery_snapshot_path()
        rec_dir = os.path.dirname(path)
        os.makedirs(rec_dir, exist_ok=True)
        tmp_path: Optional[str] = None
        try:
            with tempfile.NamedTemporaryFile(
                "w", encoding="utf-8", dir=rec_dir, prefix=".rec_", suffix=".ab", delete=False
            ) as tmp:
                tmp_path = tmp.name
                tmp.write(serialize_ab_project(filtered_proj))
            os.replace(tmp_path, path)
        except Exception:
            if tmp_path and os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass
            raise

    def load_recovery_snapshot(self) -> AntibodyProject:
        path = self.recovery_snapshot_path()
        proj = load_ab_file(path)
        self.project = proj
        reflow_project_coordinates(self.project, theme=self.theme)
        self._undo_stack.clear()
        self._redo_stack.clear()
        self.is_dirty = True
        self._reset_selection()
        return proj

    def clear_recovery_snapshot(self) -> None:
        path = self.recovery_snapshot_path()
        if os.path.exists(path):
            try:
                os.remove(path)
            except OSError:
                pass

    def load_prefs(self, path: Optional[str] = None) -> Dict[str, Any]:
        prefs = config.load_user_prefs(path)
        self.recent_files = list(prefs.get("recent_files", []))
        self.theme_mode = prefs.get("theme_mode", self.theme_mode)
        self.snap_to_grid = bool(prefs.get("snap_to_grid", self.snap_to_grid))
        self.include_variable_dev = bool(
            prefs.get("include_variable_dev", self.include_variable_dev)
        )
        return prefs

    def save_prefs(self, path: Optional[str] = None) -> None:
        prefs = {
            "recent_files": list(self.recent_files[:5]),
            "theme_mode": self.theme_mode,
            "snap_to_grid": self.snap_to_grid,
            "include_variable_dev": self.include_variable_dev,
        }
        config.save_user_prefs(prefs, path)

    def export_batch_csv(
        self, chain_id: str, domain_id: str, batch_str: str, path: str
    ) -> Tuple[int, List[str]]:
        import csv

        _, dom, _ = self._get_domain(chain_id, domain_id)
        base_seq, res_errs = resolve_domain_sequence(dom, self.get_preset_map())
        if res_errs or not base_seq:
            return 0, res_errs or ["Domain has no sequence."]
        normalized_batch = ", ".join(
            part.strip()
            for line in batch_str.splitlines()
            for part in line.split(",")
            if part.strip()
        )
        results = generate_batch(base_seq, normalized_batch, dom.isotype)
        all_errs: List[str] = []

        def _csv_safe(val: str) -> str:
            if val and val[0] in ("=", "+", "-", "@", "\t", "\r"):
                return "'" + val
            return val

        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["Variant", "Sequence", "Length_aa", "MW_kDa", "Theoretical_pI", "Errors"]
            )
            for label, seq, errs in results:
                mw, pi = compute_protein_metrics(seq) if seq else (0.0, 0.0)
                err_txt = "; ".join(errs)
                if errs:
                    all_errs.extend(errs)
                writer.writerow(
                    [
                        _csv_safe(label),
                        _csv_safe(seq),
                        len(seq),
                        f"{mw:.2f}",
                        f"{pi:.2f}",
                        _csv_safe(err_txt),
                    ]
                )
        return len(results), all_errs

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
                self.clipboard_expires_in = 0
                self._clipboard_deadline = 0.0

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
                self.clipboard_expires_in = 30
                self._clipboard_deadline = time.monotonic() + 30.0
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
        self._shutdown_flag = True
        with self._clipboard_lock:
            self._clipboard_gen += 1
            if self._clipboard_timer is not None:
                self._clipboard_timer.cancel()
                self._clipboard_timer = None
            self.clear_clipboard()


def _build_canvas_shapes(
    scene: SchematicScene,
    theme: Optional[SchematicTheme] = None,
    selected_pair: Optional[Tuple[Optional[str], Optional[str]]] = None,
) -> List[Any]:
    active_theme = theme if theme is not None else DEFAULT_SCHEMATIC_THEME
    shapes: List[Any] = []
    for pedge in scene.pairing_edges:
        if pedge.chain_a_id == pedge.chain_b_id and pedge.pairing_type == "VH_VL":
            if active_theme.orientation == "vertical":
                ctrl_x = max(pedge.x1, pedge.x2) + 34.0
                ctrl_y = (pedge.y1 + pedge.y2) / 2.0
            else:
                ctrl_x = (pedge.x1 + pedge.x2) / 2.0
                ctrl_y = min(pedge.y1, pedge.y2) - 30.0
            shapes.append(
                cv.Path(
                    elements=[
                        cv.Path.MoveTo(pedge.x1, pedge.y1),
                        cv.Path.QuadraticTo(ctrl_x, ctrl_y, pedge.x2, pedge.y2),
                    ],
                    paint=ft.Paint(
                        color="#D97757",
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
                        paint=ft.Paint(color="#F59E0B", stroke_width=2),
                    )
                )
        elif pedge.is_disulfide:
            shapes.append(
                cv.Line(
                    pedge.x1,
                    pedge.y1,
                    pedge.x2,
                    pedge.y2,
                    paint=ft.Paint(color="#F59E0B", stroke_width=2.2),
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
                        color="#9E988C",
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
                paint=ft.Paint(color="#7D786D", stroke_width=2.5),
            )
        )

    sel_cid, sel_did = selected_pair if selected_pair is not None else (None, None)
    for node in scene.nodes:
        is_sel = sel_cid == node.chain_id and sel_did == node.domain_id
        border_color = "#FACC15" if is_sel else ("#EF4444" if node.error_flag else node.border_color)
        border_width = (
            max(2.8, node.border_width + 1.0)
            if is_sel
            else (max(3.0, node.border_width + 1.5) if node.error_flag else node.border_width)
        )
        fill_paint = ft.Paint(color=node.fill_color, style=ft.PaintingStyle.FILL)
        stroke_paint = ft.Paint(
            color=border_color,
            stroke_width=border_width,
            style=ft.PaintingStyle.STROKE,
        )

        if node.shape == "ellipse":
            shapes.append(cv.Oval(node.x, node.y, node.width, node.height, paint=fill_paint))
            shapes.append(cv.Oval(node.x, node.y, node.width, node.height, paint=stroke_paint))
        elif node.shape == "hexagon":
            verts = compute_hexagon_vertices(node.x, node.y, node.width, node.height)
            elems = [cv.Path.MoveTo(verts[0][0], verts[0][1])]
            for vx, vy in verts[1:]:
                elems.append(cv.Path.LineTo(vx, vy))
            elems.append(cv.Path.Close())
            shapes.append(cv.Path(elements=elems, paint=fill_paint))
            shapes.append(cv.Path(elements=elems, paint=stroke_paint))
        else:
            if node.shape == "rect":
                radius = 0.0
            elif node.shape == "capsule":
                radius = min(node.width, node.height) / 2.0
            else:
                radius = 8.0
            shapes.append(
                cv.Rect(
                    node.x,
                    node.y,
                    node.width,
                    node.height,
                    border_radius=radius,
                    paint=fill_paint,
                )
            )
            shapes.append(
                cv.Rect(
                    node.x,
                    node.y,
                    node.width,
                    node.height,
                    border_radius=radius,
                    paint=stroke_paint,
                )
            )
    return shapes


def _ui_palette(mode: str) -> Dict[str, str]:
    if mode == "light":
        return {
            "bg": "#F8FAFC",
            "panel": "#FFFFFF",
            "elevated": "#F1F5F9",
            "border": "#CBD5E1",
            "text": "#0F172A",
            "muted": "#475569",
            "canvas_bg": "#F8FAFC",
            "chain_text": "#0F172A",
            "badge_text": "#B91C1C",
            "chip_bg": "#E2E8F0",
            "chip_text": "#0F172A",
        }
    return {
        "bg": "#141413",
        "panel": "#1E1E1C",
        "elevated": "#262521",
        "border": "#3D3B35",
        "text": "#FAF9F5",
        "muted": "#A6A299",
        "canvas_bg": "#141413",
        "chain_text": "#FAF9F5",
        "badge_text": "#F09A7E",
        "chip_bg": "#262521",
        "chip_text": "#FAF9F5",
    }


def build_studio_ui(
    page: ft.Page, controller: Optional[AntibodyStudioController] = None
) -> AntibodyStudioController:
    """Build the minimal top-tabbed Flet Antibody & Bispecific Domain Engineering Studio."""
    ctrl = controller if controller is not None else AntibodyStudioController()
    pal = _ui_palette(ctrl.theme_mode)

    def update_page_title() -> None:
        dirty_mark = " *" if ctrl.is_dirty else ""
        page.title = f"{ctrl.project.name}{dirty_mark} — Fc Engineering Studio"

    update_page_title()
    page.padding = 10
    page.theme_mode = (
        ft.ThemeMode.LIGHT if ctrl.theme_mode == "light" else ft.ThemeMode.DARK
    )
    page.bgcolor = pal["bg"]
    page.window.prevent_close = True

    file_picker = ft.FilePicker()
    if hasattr(page, "services") and file_picker not in page.services:
        page.services.append(file_picker)

    def notify(
        msg: str,
        is_error: bool = False,
        action_label: Optional[str] = None,
        on_action: Optional[Any] = None,
    ) -> None:
        snack = ft.SnackBar(
            content=ft.Text(msg),
            bgcolor="#DC2626" if is_error else "#0F766E",
            action=action_label,
            on_action=on_action,
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(snack)

    async def on_window_event(e: ft.WindowEvent) -> None:
        if e.type != ft.WindowEventType.CLOSE:
            return
        if not ctrl.is_dirty or not hasattr(page, "show_dialog"):
            try:
                ctrl.clear_recovery_snapshot()
                ctrl.shutdown()
            finally:
                await page.window.destroy()
            return

        async def _save_and_close(ev: ft.ControlEvent) -> None:
            close_dlg.open = False
            page.update()
            try:
                if ctrl.current_file_path:
                    ctrl.save_ab(ctrl.current_file_path)
                else:
                    path = await file_picker.save_file(
                        dialog_title="Save Antibody Project (.ab)",
                        file_name=f"{ctrl.project.name.replace(' ', '_')}.ab",
                        allowed_extensions=["ab"],
                    )
                    if not path:
                        return
                    if not path.endswith(".ab"):
                        path += ".ab"
                    ctrl.save_ab(path)
                ctrl.clear_recovery_snapshot()
                ctrl.save_prefs()
                ctrl.shutdown()
                await page.window.destroy()
            except Exception as err:
                notify(f"Save failed: {err}", is_error=True)

        async def _discard_and_close(ev: ft.ControlEvent) -> None:
            close_dlg.open = False
            page.update()
            try:
                ctrl.clear_recovery_snapshot()
                ctrl.save_prefs()
                ctrl.shutdown()
            finally:
                await page.window.destroy()

        def _cancel_close(ev: ft.ControlEvent) -> None:
            close_dlg.open = False
            page.update()

        close_dlg = ft.AlertDialog(
            modal=True,
            title=ft.Text("Unsaved Changes"),
            content=ft.Text(
                f"Project '{ctrl.project.name}' has unsaved changes. Save before closing?"
            ),
            actions=[
                ft.TextButton("Cancel", on_click=_cancel_close),
                ft.TextButton("Discard & Close", on_click=_discard_and_close),
                ft.Button("Save & Close", on_click=_save_and_close),
            ],
        )
        page.show_dialog(close_dlg)

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
        expand=True,
        dense=True,
        max_length=MAX_METADATA_LENGTH,
        input_filter=ft.InputFilter(allow=False, regex_string=r"[\r\n]"),
    )

    format_dd = ft.Dropdown(
        label="Format Template",
        value=ctrl.project.format_type,
        expand=True,
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

    def open_fc_mutations_modal(chain_id: str) -> None:
        chain = ctrl._get_chain(chain_id)
        modal_pal = _ui_palette(ctrl.theme_mode)
        has_fc = any(d.domain_type in ("CH2", "CH3", "Hinge") for d in chain.domains)
        initial_hint = (
            ""
            if (has_fc or ctrl.apply_to_all_heavy)
            else "Selected chain has no CH2/CH3/Hinge domain. Enable 'Apply to all heavy chains' or select a heavy chain."
        )
        modal_status = ft.Text(
            initial_hint,
            size=11,
            color="#F59E0B" if initial_hint else modal_pal["muted"],
        )
        preset_rows_col = ft.Column(spacing=4, scroll=ft.ScrollMode.AUTO)

        cat_state = {"cat": "All"}

        file_preset_labels = {
            m.get("label")
            for m in ctrl.project.file_custom_mutations
            if isinstance(m, dict)
        }

        def _classify_preset(lbl: str, val: str) -> str:
            if lbl in file_preset_labels:
                return "Project"
            low = f"{lbl} {val}".lower()
            if any(k in low for k in ("yte", "ls", "half", "fcrn", "m252", "m428", "n434")):
                return "Half-Life"
            if any(k in low for k in ("lala", "pg", "agly", "n297", "silent", "effector", "adcc", "cdc", "gasdalie", "s239d", "k322")):
                return "Effector"
            if any(k in low for k in ("knob", "hole", "kih", "charge", "t366", "y407", "k392", "d356", "k409", "d399")):
                return "Heterodimer"
            return "Glyco/Other"

        def refresh_modal_presets() -> None:
            active_presets = {p for d in chain.domains for p in d.applied_presets}
            preset_map = ctrl.get_preset_map()
            q = preset_filter_state["query"].strip().lower()
            cat = cat_state["cat"]
            filtered = [
                (lbl, val)
                for lbl, val in preset_map.items()
                if (not q or q in lbl.lower() or q in val.lower())
                and (cat == "All" or _classify_preset(lbl, val) == cat)
            ]
            rows: List[ft.Control] = []

            def make_cb(lbl: str, val: str):
                def _on_cb(ev: ft.ControlEvent) -> None:
                    errs = ctrl.toggle_preset(chain.id, lbl, bool(ev.control.value))
                    if errs:
                        modal_status.value = "; ".join(errs)
                        modal_status.color = "#DC2626"
                    else:
                        modal_status.value = f"{'Applied' if ev.control.value else 'Removed'} {lbl} ({val})."
                        modal_status.color = "#0F766E"
                    refresh_modal_presets()
                    render_workspace()

                cb = ft.Checkbox(
                    label=f"{lbl} ({val})",
                    value=(lbl in active_presets),
                    expand=True,
                    on_change=_on_cb,
                )
                if lbl in file_preset_labels:
                    def _on_del_file_preset(ev: ft.ControlEvent) -> None:
                        ctrl.delete_file_custom_mutation(lbl)
                        modal_status.value = f"Deleted file preset '{lbl}'."
                        modal_status.color = "#0F766E"
                        refresh_modal_presets()
                        render_workspace()

                    return ft.Row(
                        spacing=2,
                        expand=True,
                        controls=[
                            cb,
                            ft.IconButton(
                                icon=ft.Icons.DELETE_OUTLINE,
                                icon_size=15,
                                tooltip=f"Delete file preset {lbl}",
                                on_click=_on_del_file_preset,
                            ),
                        ],
                    )
                return cb

            for i in range(0, len(filtered), 2):
                pair = [make_cb(*filtered[i])]
                if i + 1 < len(filtered):
                    pair.append(make_cb(*filtered[i + 1]))
                else:
                    pair.append(ft.Container(expand=True))
                rows.append(ft.Row(spacing=8, controls=pair))
            preset_rows_col.controls = rows

        cat_row = ft.Row(spacing=4, wrap=True)

        def refresh_cat_chips() -> None:
            cats = ["All", "Half-Life", "Effector", "Heterodimer", "Glyco/Other", "Project"]
            chips: List[ft.Control] = []
            for cname in cats:
                selected = cat_state["cat"] == cname

                def _make_cat(c_val: str):
                    def _on_cat(ev: ft.ControlEvent) -> None:
                        cat_state["cat"] = c_val
                        refresh_cat_chips()
                        refresh_modal_presets()
                        page.update()

                    return _on_cat

                chips.append(
                    ft.Chip(
                        label=ft.Text(cname, size=11),
                        selected=selected,
                        on_select=_make_cat(cname),
                    )
                )
            cat_row.controls = chips

        refresh_cat_chips()

        def on_modal_filter(ev: ft.ControlEvent) -> None:
            preset_filter_state["query"] = ev.control.value or ""
            refresh_modal_presets()
            page.update()

        def _target_chains_summary() -> str:
            targets = ctrl._resolve_mutation_target_chains(chain.id, None)
            return "Applies to: " + ", ".join(c.name for c in targets)

        applies_to_lbl = ft.Text(
            _target_chains_summary(),
            size=10,
            color=modal_pal["muted"],
        )

        def on_modal_all_heavy(ev: ft.ControlEvent) -> None:
            ctrl.apply_to_all_heavy = bool(ev.control.value)
            applies_to_lbl.value = _target_chains_summary()
            if has_fc or ctrl.apply_to_all_heavy:
                modal_status.value = ""
            else:
                modal_status.value = "Selected chain has no CH2/CH3/Hinge domain. Enable 'Apply to all heavy chains' or select a heavy chain."
                modal_status.color = "#F59E0B"
            render_workspace()

        search_tf = ft.TextField(
            label="Filter Fc Presets",
            value=preset_filter_state["query"],
            autofocus=True,
            dense=True,
            expand=True,
            on_change=on_modal_filter,
        )
        all_heavy_sw = ft.Switch(
            label="Apply to all heavy chains",
            value=ctrl.apply_to_all_heavy,
            on_change=on_modal_all_heavy,
        )

        wt_lookup_lbl = ft.Text(
            "EU ranges: CH1 118–215 · Hinge 216–230 · CH2 231–340 · CH3 341–447 (type e.g. 234 for WT preview).",
            size=10,
            color=modal_pal["muted"],
        )

        def on_adhoc_change(ev: ft.ControlEvent) -> None:
            import re

            raw = (ev.control.value or "").strip()
            m = re.search(r"(\d{3})", raw)
            if not m:
                wt_lookup_lbl.value = "Tip: type an EU position (e.g. 234) to preview WT residue."
                page.update()
                return
            pos = int(m.group(1))
            iso_dom = next(
                (d for d in chain.domains if d.domain_type in CONSTANT_DOMAIN_EU_RANGES),
                None,
            )
            iso = iso_dom.isotype if iso_dom is not None else "igg1"
            allo = iso_dom.allotype if iso_dom is not None else ""
            iso_seqs = config.SEQUENCES.get(iso, {})
            full_seq = iso_seqs.get(allo) or (next(iter(iso_seqs.values())) if iso_seqs else "")
            aa = lookup_wt_residue(full_seq, pos, iso) if full_seq else None
            if aa:
                wt_lookup_lbl.value = f"WT @ EU {pos} ({iso.upper()}): {aa} — e.g. {aa}{pos}A"
            else:
                wt_lookup_lbl.value = f"EU {pos} not in {iso.upper()} constant region."
            page.update()

        custom_eu_tf = ft.TextField(
            label="Apply Ad-hoc EU Mutation (e.g. M252Y/S254T)",
            dense=True,
            expand=True,
            max_length=1000,
            input_filter=ft.InputFilter(allow=True, regex_string=r"[a-zA-Z0-9/, ]"),
            on_change=on_adhoc_change,
        )

        def on_modal_adhoc(ev: ft.ControlEvent) -> None:
            errs = ctrl.apply_custom_mutation(chain.id, custom_eu_tf.value or "")
            if errs:
                modal_status.value = "; ".join(errs)
                modal_status.color = "#DC2626"
            else:
                modal_status.value = "Applied custom EU mutation."
                modal_status.color = "#0F766E"
                custom_eu_tf.value = ""
            refresh_modal_presets()
            render_workspace()

        file_mut_label_tf = ft.TextField(
            label="New File Preset Label",
            dense=True,
            expand=True,
            max_length=40,
            input_filter=ft.InputFilter(allow=False, regex_string=r"[\r\n]"),
        )
        file_mut_val_tf = ft.TextField(
            label="Preset EU Mutations (e.g. L234A/L235A)",
            dense=True,
            expand=True,
            max_length=1000,
            input_filter=ft.InputFilter(allow=True, regex_string=r"[a-zA-Z0-9/, ]"),
        )

        def on_modal_add_file_mut(ev: ft.ControlEvent) -> None:
            errs = ctrl.register_file_custom_mutation(
                file_mut_label_tf.value or "", file_mut_val_tf.value or ""
            )
            if errs:
                modal_status.value = "; ".join(errs)
                modal_status.color = "#DC2626"
            else:
                modal_status.value = "Added custom preset to current project file."
                modal_status.color = "#0F766E"
                file_mut_label_tf.value = ""
                file_mut_val_tf.value = ""
            refresh_modal_presets()
            render_workspace()

        def on_reset_chain_wt(ev: ft.ControlEvent) -> None:
            ctrl.reset_chain_mutations(chain.id)
            modal_status.value = f"Reset all mutations on {chain.name} to WT."
            modal_status.color = "#0F766E"
            refresh_modal_presets()
            render_workspace()

        def on_close_fc_modal(ev: ft.ControlEvent) -> None:
            dlg.open = False
            page.update()

        refresh_modal_presets()
        dlg = ft.AlertDialog(
            title=ft.Text(f"Edit Fc Mutations — Chain: {chain.name}"),
            content=ft.Container(
                width=740,
                content=ft.Column(
                    tight=True,
                    spacing=8,
                    controls=[
                        ft.Row(spacing=12, controls=[search_tf, all_heavy_sw]),
                        ft.Row(spacing=12, alignment=ft.MainAxisAlignment.SPACE_BETWEEN, controls=[cat_row, applies_to_lbl]),
                        ft.Container(
                            height=240,
                            border=ft.Border.all(1, modal_pal["border"]),
                            border_radius=6,
                            padding=8,
                            content=preset_rows_col,
                        ),
                        ft.Row(
                            spacing=8,
                            controls=[custom_eu_tf, ft.Button("Apply", on_click=on_modal_adhoc)],
                        ),
                        wt_lookup_lbl,
                        ft.Row(
                            spacing=8,
                            controls=[
                                file_mut_label_tf,
                                file_mut_val_tf,
                                ft.OutlinedButton("Save to Project File", on_click=on_modal_add_file_mut),
                            ],
                        ),
                        modal_status,
                    ],
                ),
            ),
            actions=[
                ft.TextButton("Reset Chain to WT", on_click=on_reset_chain_wt),
                ft.Button("Done", on_click=on_close_fc_modal),
            ],
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(dlg)

    def rebuild_canvas_only() -> None:
        scene = build_schematic_scene(
            ctrl.project, ctrl.get_preset_map(), theme=ctrl.theme
        )
        canvas_bg.shapes = _build_canvas_shapes(
            scene,
            theme=ctrl.theme,
            selected_pair=(ctrl.selected_chain_id, ctrl.selected_domain_id),
        )
        max_x, max_y = scene.bounds[2], scene.bounds[3]
        canvas_w = max(860.0, max_x + 48.0)
        canvas_h = max(420.0, max_y + 48.0)
        canvas_overlay.width = canvas_w
        canvas_overlay.height = canvas_h

        controls: List[ft.Control] = [
            ft.Container(content=canvas_bg, width=canvas_w, height=canvas_h)
        ]

        if not scene.nodes:
            active_pal = _ui_palette(ctrl.theme_mode)
            controls.append(
                ft.Container(
                    left=40,
                    top=80,
                    padding=20,
                    border=ft.Border.all(1, active_pal["border"]),
                    border_radius=10,
                    bgcolor=active_pal["panel"],
                    content=ft.Column(
                        spacing=6,
                        controls=[
                            ft.Text(
                                "Canvas is empty",
                                size=14,
                                weight=ft.FontWeight.BOLD,
                                color=active_pal["text"],
                            ),
                            ft.Text(
                                "Drag a domain block from the Left Palette onto a chain, or click '+' in the Chain Strip.",
                                size=12,
                                color=active_pal["muted"],
                            ),
                        ],
                    ),
                )
            )

        for node in scene.nodes:
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

            def make_double_tap_handler(cid: str, did: str, dtype: str):
                def _on_double_tap(e: ft.ControlEvent) -> None:
                    ctrl.select_domain(cid, did)
                    render_workspace()
                    if dtype in CONSTANT_DOMAIN_EU_RANGES:
                        open_fc_mutations_modal(cid)
                    else:
                        notify(f"Selected {dtype} — edit template or custom sequence in Inspector.")

                return _on_double_tap

            drag_state: Dict[str, Any] = {"snapshot": None, "orig_xy": None}

            def make_pan_start_handler(cid: str, did: str):
                def _on_pan_start(e: Any) -> None:
                    for ch in ctrl.project.chains:
                        if ch.id == cid:
                            for dm in ch.domains:
                                if dm.id == did:
                                    drag_state["snapshot"] = copy.deepcopy(ctrl.project)
                                    drag_state["orig_xy"] = (dm.x, dm.y)
                                    return

                return _on_pan_start

            def make_pan_handler(cid: str, did: str):
                def _on_pan_update(e: ft.DragUpdateEvent) -> None:
                    if e.local_delta is None:
                        return
                    inv_z = 1.0 / (ctrl.zoom_scale or 1.0)
                    for ch in ctrl.project.chains:
                        if ch.id == cid:
                            for dm in ch.domains:
                                if dm.id == did:
                                    if drag_state["snapshot"] is None:
                                        drag_state["snapshot"] = copy.deepcopy(ctrl.project)
                                        drag_state["orig_xy"] = (dm.x, dm.y)
                                    ctrl.move_domain(
                                        cid,
                                        did,
                                        dm.x + e.local_delta.x * inv_z,
                                        dm.y + e.local_delta.y * inv_z,
                                    )
                                    rebuild_canvas_only()
                                    page.update()
                                    return

                return _on_pan_update

            def make_pan_end_handler(cid: str, did: str):
                def _on_pan_end(e: Any) -> None:
                    snap = drag_state.get("snapshot")
                    orig_xy = drag_state.get("orig_xy")
                    drag_state["snapshot"] = None
                    drag_state["orig_xy"] = None
                    moved = False
                    for ch in ctrl.project.chains:
                        if ch.id == cid:
                            for dm in ch.domains:
                                if dm.id == did:
                                    if ctrl.snap_to_grid:
                                        ctrl.move_domain(
                                            cid,
                                            did,
                                            round(dm.x / 10.0) * 10.0,
                                            round(dm.y / 10.0) * 10.0,
                                        )
                                    moved = orig_xy is None or (dm.x, dm.y) != orig_xy
                                    break
                    if snap is not None and moved:
                        ctrl._push_undo(snap)
                    elif moved:
                        ctrl.is_dirty = True
                    render_workspace()

                return _on_pan_end

            box_col: List[ft.Control] = []
            active_pal = _ui_palette(ctrl.theme_mode)
            if node.is_n_terminus and node.chain_name:
                box_col.append(
                    ft.Text(
                        node.chain_name,
                        size=10,
                        weight=ft.FontWeight.BOLD,
                        color=active_pal["chain_text"],
                    )
                )
            mut_summary = ""
            if dom_lookup is not None:
                muts = list(dom_lookup.applied_presets) + list(dom_lookup.custom_mutations)
                if muts:
                    mut_summary = f" | Mutations: {', '.join(muts[:4])}"
            tip_text = f"{node.chain_name} · {node.domain_type}{mut_summary} (Double-click to edit)"
            box_col.append(
                ft.Container(
                    width=node.width,
                    height=node.height,
                    bgcolor=ft.Colors.TRANSPARENT,
                    border=None,
                    alignment=ft.Alignment.CENTER,
                    tooltip=tip_text,
                    content=ft.Text(
                        f"📌 {node.label}" if (dom_lookup is not None and dom_lookup.pinned) else node.label,
                        color=node.text_color,
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
                        color=active_pal["badge_text"],
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
                            on_double_tap=make_double_tap_handler(
                                node.chain_id, node.domain_id, node.domain_type
                            ),
                            on_pan_start=make_pan_start_handler(
                                node.chain_id, node.domain_id
                            ),
                            on_pan_update=make_pan_handler(
                                node.chain_id, node.domain_id
                            ),
                            on_pan_end=make_pan_end_handler(
                                node.chain_id, node.domain_id
                            ),
                            content=ft.Column(spacing=1, controls=box_col),
                        ),
                    )
                )

        for pedge in scene.pairing_edges:
            if pedge.is_disulfide:
                mx = max(0.0, (pedge.x1 + pedge.x2) / 2.0 - 18.0)
                my = max(0.0, (pedge.y1 + pedge.y2) / 2.0 - 8.0)

                def make_ss_tap(ca: str, da: str, cb: str, db: str, ptype: str):
                    def _on_ss_tap(ev: ft.ControlEvent) -> None:
                        match = next(
                            (
                                p
                                for p in ctrl.project.pairings
                                if (p.chain_a_id, p.domain_a_id, p.chain_b_id, p.domain_b_id, p.pairing_type)
                                == (ca, da, cb, db, ptype)
                            ),
                            None,
                        )
                        if match is not None:
                            ctrl.remove_pairing(match.id)
                            render_workspace()
                            notify(f"Removed pairing {match.id} ({ptype}).")

                    return _on_ss_tap

                controls.append(
                    ft.Container(
                        left=mx,
                        top=my,
                        width=36,
                        height=16,
                        bgcolor="#141413",
                        border=ft.Border.all(1, "#F59E0B"),
                        border_radius=4,
                        alignment=ft.Alignment.CENTER,
                        tooltip=f"{pedge.pairing_type} disulfide — click to remove",
                        on_click=make_ss_tap(
                            pedge.chain_a_id,
                            pedge.domain_a_id,
                            pedge.chain_b_id,
                            pedge.domain_b_id,
                            pedge.pairing_type,
                        ),
                        content=ft.Text(
                            "-S-S-",
                            size=8,
                            weight=ft.FontWeight.BOLD,
                            color="#FDE68A",
                        ),
                    )
                )

        canvas_overlay.controls = controls

    def rebuild_chain_strip() -> None:
        rows: List[ft.Control] = []
        preset_map = ctrl.get_preset_map()
        for c_idx, chain in enumerate(ctrl.project.chains):
            chip_controls: List[ft.Control] = []
            active_pal = _ui_palette(ctrl.theme_mode)
            total_aa = 0
            for idx, dom in enumerate(chain.domains):
                fill_col = DOMAIN_STYLE_MAP[dom.domain_type][2]
                is_sel = (
                    ctrl.selected_chain_id == chain.id
                    and ctrl.selected_domain_id == dom.id
                )
                dseq, derrs = resolve_domain_sequence(dom, preset_map)
                total_aa += len(dseq)
                has_mut = bool(dom.applied_presets or dom.custom_mutations)

                def make_chip_click(cid: str, did: str):
                    def _on_click(e: ft.ControlEvent) -> None:
                        ctrl.select_domain(cid, did)
                        render_workspace()

                    return _on_click

                chip_label = ft.Text(
                    dom.domain_type,
                    size=11,
                    weight=ft.FontWeight.W_600,
                    color="#FFFFFF" if is_sel else active_pal["chip_text"],
                )
                chip_dot: Optional[ft.Control] = None
                if derrs:
                    chip_dot = ft.Container(
                        width=6,
                        height=6,
                        border_radius=3,
                        bgcolor="#EF4444",
                        tooltip="; ".join(derrs),
                    )
                elif has_mut:
                    chip_dot = ft.Container(
                        width=6,
                        height=6,
                        border_radius=3,
                        bgcolor="#F59E0B",
                        tooltip="Has mutations",
                    )

                chip_badge = ft.Container(
                    padding=ft.Padding.symmetric(horizontal=8, vertical=4),
                    bgcolor=fill_col if is_sel else active_pal["chip_bg"],
                    border=ft.Border.all(2.0 if is_sel else 1.2, "#FACC15" if is_sel else fill_col),
                    border_radius=6,
                    on_click=make_chip_click(chain.id, dom.id),
                    content=ft.Row(spacing=4, tight=True, controls=[chip_label, chip_dot])
                    if chip_dot is not None
                    else chip_label,
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
                        border=ft.Border.all(1, active_pal["border"]),
                        border_radius=6,
                        content=ft.Text("+ Drop", size=10, color=active_pal["muted"]),
                    ),
                )
            )

            def make_quick_add(cid: str, dtype: str):
                def _on_add(e: ft.ControlEvent) -> None:
                    try:
                        ctrl.add_domain_to_chain(cid, dtype)
                        render_workspace()
                    except Exception as err:
                        notify(str(err), is_error=True)

                return _on_add

            def make_macro_add(cid: str, kind: str):
                def _on_macro(e: ft.ControlEvent) -> None:
                    try:
                        ctrl.insert_macro_block(cid, kind)
                        render_workspace()
                        notify(f"Inserted {kind.upper()} macro block.")
                    except Exception as err:
                        notify(str(err), is_error=True)

                return _on_macro

            quick_add_btn = ft.PopupMenuButton(
                tooltip="Quick Add Domain",
                icon=ft.Icons.ADD_CIRCLE_OUTLINE,
                icon_size=16,
                items=[
                    *[
                        ft.PopupMenuItem(content=dt, on_click=make_quick_add(chain.id, dt))
                        for dt in VALID_DOMAIN_TYPES
                    ],
                    ft.PopupMenuItem(),
                    ft.PopupMenuItem(
                        content="Macro: Insert scFv (VH-G4S-VL)",
                        on_click=make_macro_add(chain.id, "scfv"),
                    ),
                    ft.PopupMenuItem(
                        content="Macro: Insert Fc (Hinge-CH2-CH3)",
                        on_click=make_macro_add(chain.id, "fc"),
                    ),
                ],
            )

            def make_rename_chain(cid: str, current_name: str):
                def _on_rename(e: ft.ControlEvent) -> None:
                    name_tf = ft.TextField(
                        label="Chain Name",
                        value=current_name,
                        autofocus=True,
                        max_length=MAX_METADATA_LENGTH,
                    )

                    def _commit(ev: ft.ControlEvent) -> None:
                        try:
                            ctrl.rename_chain(cid, name_tf.value or "")
                            rename_dlg.open = False
                            render_workspace()
                        except Exception as err:
                            notify(str(err), is_error=True)

                    def _cancel(ev: ft.ControlEvent) -> None:
                        rename_dlg.open = False
                        page.update()

                    rename_dlg = ft.AlertDialog(
                        title=ft.Text("Rename Chain"),
                        content=name_tf,
                        actions=[
                            ft.TextButton("Cancel", on_click=_cancel),
                            ft.Button("Save", on_click=_commit),
                        ],
                    )
                    if hasattr(page, "show_dialog"):
                        page.show_dialog(rename_dlg)

                return _on_rename

            def make_reorder_chain(cid: str, delta: int):
                def _on_reord(e: ft.ControlEvent) -> None:
                    ctrl.reorder_chain(cid, delta)
                    render_workspace()

                return _on_reord

            def make_clone_chain(cid: str):
                def _on_clone(e: ft.ControlEvent) -> None:
                    try:
                        cloned = ctrl.clone_chain(cid)
                        render_workspace()
                        notify(f"Cloned chain '{cloned.name}'.")
                    except Exception as err:
                        notify(str(err), is_error=True)

                return _on_clone

            def make_rm_chain(cid: str, cname: str):
                def _on_rm(e: ft.ControlEvent) -> None:
                    ctrl.remove_chain(cid)
                    render_workspace()
                    notify(
                        f"Removed chain '{cname}'.",
                        action_label="Undo",
                        on_action=lambda _: (ctrl.undo(), render_workspace()),
                    )

                return _on_rm

            rows.append(
                ft.Container(
                    padding=6,
                    bgcolor=active_pal["elevated"],
                    border=ft.Border.all(1, active_pal["border"]),
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
                                        color=active_pal["text"],
                                    ),
                                    ft.Text(
                                        f"({len(chain.domains)} doms · {total_aa} aa)",
                                        size=10,
                                        color=active_pal["muted"],
                                    ),
                                    *chip_controls,
                                    quick_add_btn,
                                ],
                            ),
                            ft.Row(
                                spacing=0,
                                tight=True,
                                controls=[
                                    ft.IconButton(
                                        icon=ft.Icons.EDIT_OUTLINED,
                                        icon_size=15,
                                        tooltip="Rename Chain",
                                        on_click=make_rename_chain(chain.id, chain.name),
                                    ),
                                    ft.IconButton(
                                        icon=ft.Icons.ARROW_UPWARD,
                                        icon_size=15,
                                        tooltip="Move Chain Up",
                                        disabled=(c_idx == 0),
                                        on_click=make_reorder_chain(chain.id, -1),
                                    ),
                                    ft.IconButton(
                                        icon=ft.Icons.ARROW_DOWNWARD,
                                        icon_size=15,
                                        tooltip="Move Chain Down",
                                        disabled=(c_idx == len(ctrl.project.chains) - 1),
                                        on_click=make_reorder_chain(chain.id, 1),
                                    ),
                                    ft.IconButton(
                                        icon=ft.Icons.CONTENT_COPY,
                                        icon_size=15,
                                        tooltip="Clone Chain",
                                        on_click=make_clone_chain(chain.id),
                                    ),
                                    ft.IconButton(
                                        icon=ft.Icons.DELETE_OUTLINE,
                                        icon_size=16,
                                        tooltip="Remove Chain",
                                        on_click=make_rm_chain(chain.id, chain.name),
                                    ),
                                ],
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
                page.update()

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

        def on_auto_detect(e: ft.ControlEvent) -> None:
            added = ctrl.auto_detect_pairings()
            render_workspace()
            notify(f"Auto-detected {added} standard pairing(s)." if added else "No new standard pairings detected.")

        pairing_controls_col.controls = [
            dd_a,
            dd_b,
            dd_ptype,
            ft.Row(
                spacing=6,
                controls=[
                    ft.Button("Add Pairing", on_click=on_add_pairing, expand=True),
                    ft.OutlinedButton("Auto-Detect", on_click=on_auto_detect),
                ],
            ),
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
        dom_idx = -1
        if ctrl.selected_domain_id:
            for idx_d, d in enumerate(chain.domains):
                if d.id == ctrl.selected_domain_id:
                    dom = d
                    dom_idx = idx_d
                    break

        active_pal = _ui_palette(ctrl.theme_mode)
        crumb = f"{ctrl.project.name} > {chain.name}" + (f" > {dom.domain_type}" if dom else "")
        controls.append(ft.Text(crumb, size=10, color=active_pal["muted"]))

        controls.append(
            ft.Row(
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                controls=[
                    ft.Text(
                        f"Chain: {chain.name}",
                        size=14,
                        weight=ft.FontWeight.BOLD,
                    ),
                    ft.OutlinedButton(
                        "Reset Layout",
                        icon=ft.Icons.AUTO_FIX_HIGH,
                        on_click=lambda e: (
                            ctrl.reset_layout(),
                            render_workspace(),
                        ),
                    ),
                ],
            )
        )

        if dom is not None:
            dom_actions: List[ft.Control] = [
                ft.IconButton(
                    icon=ft.Icons.CHEVRON_LEFT,
                    icon_size=16,
                    tooltip="Move Domain Left",
                    disabled=(dom_idx <= 0),
                    on_click=lambda e, cid=chain.id, idx=dom_idx: (
                        ctrl.reorder_domain(cid, idx, idx - 1),
                        render_workspace(),
                    ),
                ),
                ft.IconButton(
                    icon=ft.Icons.CHEVRON_RIGHT,
                    icon_size=16,
                    tooltip="Move Domain Right",
                    disabled=(dom_idx < 0 or dom_idx >= len(chain.domains) - 1),
                    on_click=lambda e, cid=chain.id, idx=dom_idx: (
                        ctrl.reorder_domain(cid, idx, idx + 1),
                        render_workspace(),
                    ),
                ),
            ]
            if dom.pinned:
                dom_actions.append(
                    ft.OutlinedButton(
                        "Unpin",
                        on_click=lambda e, cid=chain.id, did=dom.id: (
                            ctrl.unpin_domain(cid, did),
                            render_workspace(),
                        ),
                    )
                )
            dom_actions.append(
                ft.OutlinedButton(
                    "Delete",
                    on_click=lambda e, cid=chain.id, did=dom.id, dt=dom.domain_type: (
                        ctrl.remove_domain(cid, did),
                        render_workspace(),
                        notify(
                            f"Deleted {dt} domain.",
                            action_label="Undo",
                            on_action=lambda _: (ctrl.undo(), render_workspace()),
                        ),
                    ),
                )
            )
            controls.append(
                ft.Row(
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    controls=[
                        ft.Text(
                            f"Selected Domain: {dom.domain_type} ({dom.id})",
                            size=12,
                            weight=ft.FontWeight.W_600,
                        ),
                        ft.Row(spacing=4, controls=dom_actions),
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
                seq_len_txt = ft.Text(
                    f"{len(dom.custom_sequence)} / {MAX_CUSTOM_SEQ_LENGTH}",
                    size=10,
                    color=active_pal["muted"],
                )
                seq_tf = ft.TextField(
                    label=f"Custom {dom.domain_type} Sequence",
                    value=dom.custom_sequence,
                    multiline=True,
                    min_lines=2,
                    max_lines=3,
                    max_length=MAX_CUSTOM_SEQ_LENGTH,
                    dense=True,
                    on_change=lambda e: (
                        setattr(seq_len_txt, "value", f"{len(e.control.value or '')} / {MAX_CUSTOM_SEQ_LENGTH}"),
                        page.update(),
                    ),
                )

                def on_tkey_select(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id, dom.id, template_key=e.control.value
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    render_workspace()

                def on_apply_cseq(e: ft.ControlEvent) -> None:
                    raw_lines = (seq_tf.value or "").splitlines()
                    cleaned = "".join(
                        ln.strip() for ln in raw_lines if not ln.lstrip().startswith(">")
                    ).replace(" ", "").upper()
                    errs = ctrl.update_domain_properties(
                        chain.id,
                        dom.id,
                        custom_sequence=cleaned,
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    else:
                        notify(f"Updated {dom.domain_type} sequence.")
                    render_workspace()

                def on_reset_cseq(e: ft.ControlEvent) -> None:
                    errs = ctrl.update_domain_properties(
                        chain.id, dom.id, custom_sequence=""
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    else:
                        notify(f"Reset {dom.domain_type} to template.")
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
                        ft.Row(
                            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                            controls=[
                                seq_len_txt,
                                ft.Row(
                                    spacing=6,
                                    controls=[
                                        ft.OutlinedButton("Reset to Template", on_click=on_reset_cseq),
                                        ft.Button("Apply Sequence", on_click=on_apply_cseq),
                                    ],
                                ),
                            ],
                        ),
                    ]
                )

        # Fc Mutation Summary + Modal Trigger
        controls.append(ft.Divider())
        active_presets_on_chain = sorted(
            {p for d in chain.domains for p in d.applied_presets}
        )
        preset_map = ctrl.get_preset_map()
        custom_mut_chips: List[ft.Control] = []
        for d in chain.domains:
            for tok in d.custom_mutations:
                def _make_rm_mut(cid: str, did: str, t: str):
                    def _on_del(ev: ft.ControlEvent) -> None:
                        ctrl.remove_custom_mutation(cid, did, t)
                        render_workspace()

                    return _on_del

                custom_mut_chips.append(
                    ft.Chip(
                        label=ft.Text(f"{d.domain_type}:{tok}", size=10),
                        on_delete=_make_rm_mut(chain.id, d.id, tok),
                    )
                )

        controls.extend(
            [
                ft.Text("Fc Mutations", size=13, weight=ft.FontWeight.BOLD),
                ft.Text(
                    "Active Presets: "
                    + (", ".join(active_presets_on_chain) if active_presets_on_chain else "None (WT)"),
                    size=11,
                    color=active_pal["muted"],
                ),
                ft.Text(
                    f"Apply to all heavy chains: {'ON' if ctrl.apply_to_all_heavy else 'OFF'}",
                    size=11,
                    color=active_pal["muted"],
                ),
                *(
                    [ft.Row(wrap=True, spacing=4, controls=custom_mut_chips)]
                    if custom_mut_chips
                    else []
                ),
                ft.Row(
                    spacing=6,
                    controls=[
                        ft.Button(
                            "Edit Fc Mutations...",
                            icon=ft.Icons.TUNE,
                            on_click=lambda e, cid=chain.id: open_fc_mutations_modal(cid),
                        ),
                        ft.OutlinedButton(
                            "Reset WT",
                            on_click=lambda e, cid=chain.id: (
                                ctrl.reset_chain_mutations(cid),
                                render_workspace(),
                                notify(f"Reset {chain.name} to WT."),
                            ),
                        ),
                    ],
                ),
            ]
        )

        # Live Sequence, Diff & Developability Summary
        controls.append(ft.Divider())
        chain_seq, segs, seq_errs = assemble_chain_sequence(chain, preset_map)
        diffs, diff_errs = diff_chain_sequences(chain, preset_map)
        liabilities, dev_errs = analyze_chain_developability(
            chain, preset_map, include_variable=ctrl.include_variable_dev
        )

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
            mw_kda, theo_pi = compute_protein_metrics(chain_seq)
            diff_chips = [
                ft.Chip(
                    label=ft.Text(f"{dt}:{wt}{pos}{mut}", size=10),
                    on_click=lambda e, cid=chain.id: open_fc_mutations_modal(cid),
                )
                for dt, pos, wt, mut in diffs
            ]
            dev_chips = [
                ft.Container(
                    padding=ft.Padding.symmetric(horizontal=6, vertical=2),
                    border_radius=4,
                    bgcolor="#78350F" if pos_list else active_pal["elevated"],
                    border=ft.Border.all(1, "#F59E0B" if pos_list else active_pal["border"]),
                    content=ft.Text(
                        f"{cat}: {len(pos_list)}",
                        size=10,
                        color="#FDE68A" if pos_list else active_pal["muted"],
                        tooltip=", ".join(str(p) for p in pos_list) or "None",
                    ),
                )
                for cat, pos_list in liabilities
            ]
            spans = []
            for start, end, seg_dom in segs:
                seg_color = DOMAIN_STYLE_MAP.get(seg_dom.domain_type, ("", "", "#9CA3AF"))[2]
                spans.append(
                    ft.TextSpan(
                        chain_seq[start:end],
                        style=ft.TextStyle(color=seg_color, size=10, font_family="monospace"),
                    )
                )

            def _on_toggle_var_dev(ev: ft.ControlEvent) -> None:
                ctrl.include_variable_dev = bool(ev.control.value)
                ctrl.save_prefs()
                render_workspace()

            def _on_copy_chain_seq(ev: ft.ControlEvent) -> None:
                if ctrl.copy_text_secure(chain_seq, label=f"{chain.name} sequence"):
                    render_workspace()
                    if hasattr(page, "run_task"):
                        page.run_task(_clipboard_countdown_loop)
                    notify(f"Copied {chain.name} sequence (auto-clears in 30s).")

            hist_items = [
                f"{d.domain_type}: {', '.join(d.mutation_history[-3:])}"
                for d in chain.domains
                if d.mutation_history
            ]
            controls.extend(
                [
                    ft.Row(
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        controls=[
                            ft.Text(
                                f"Length: {len(chain_seq)} aa · MW: {mw_kda:.1f} kDa · pI: {theo_pi:.2f}",
                                size=11,
                                weight=ft.FontWeight.BOLD,
                            ),
                            ft.IconButton(
                                icon=ft.Icons.COPY,
                                icon_size=15,
                                tooltip="Copy Chain Sequence (30s auto-clear)",
                                on_click=_on_copy_chain_seq,
                            ),
                        ],
                    ),
                    ft.Text("EU Diffs:", size=11, weight=ft.FontWeight.W_600),
                    ft.Row(wrap=True, spacing=4, controls=diff_chips)
                    if diff_chips
                    else ft.Text("None (WT)", size=10, color=active_pal["muted"]),
                    ft.Row(
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        controls=[
                            ft.Text("Developability:", size=11, weight=ft.FontWeight.W_600),
                            ft.Switch(
                                label="Scan VH/VL/G4S",
                                value=ctrl.include_variable_dev,
                                scale=0.8,
                                on_change=_on_toggle_var_dev,
                            ),
                        ],
                    ),
                    ft.Row(wrap=True, spacing=4, controls=dev_chips),
                    ft.Container(
                        padding=6,
                        border=ft.Border.all(1, active_pal["border"]),
                        border_radius=6,
                        bgcolor=active_pal["bg"],
                        content=ft.Text(
                            spans=spans if spans else None,
                            value="" if spans else chain_seq,
                            selectable=True,
                            size=10,
                        ),
                    ),
                    *(
                        [
                            ft.Text(
                                "History: " + " | ".join(hist_items),
                                size=9,
                                color=active_pal["muted"],
                            )
                        ]
                        if hist_items
                        else []
                    ),
                ]
            )

        inspector_col.controls = controls

    multispecific_col = ft.Column(spacing=8)

    def rebuild_multispecific_tools() -> None:
        chain_opts = [ft.dropdown.Option(c.id, text=c.name) for c in ctrl.project.chains]
        valid_ids = {c.id for c in ctrl.project.chains}
        if ctrl.wiz_state["a"] not in valid_ids:
            ctrl.wiz_state["a"] = chain_opts[0].key if chain_opts else None
        if ctrl.wiz_state["b"] not in valid_ids:
            ctrl.wiz_state["b"] = chain_opts[1].key if len(chain_opts) > 1 else ctrl.wiz_state["a"]

        dd_ca = ft.Dropdown(
            label="Chain A (Knob / Cross A)",
            value=ctrl.wiz_state["a"],
            options=chain_opts,
            dense=True,
            expand=True,
            on_select=lambda e: ctrl.wiz_state.__setitem__("a", e.control.value),
        )
        dd_cb = ft.Dropdown(
            label="Chain B (Hole / Cross B)",
            value=ctrl.wiz_state["b"],
            options=chain_opts,
            dense=True,
            expand=True,
            on_select=lambda e: ctrl.wiz_state.__setitem__("b", e.control.value),
        )
        dd_mode = ft.Dropdown(
            label="Heterodimer Mode",
            value=ctrl.wiz_state["mode"] or "kih_ss",
            options=[
                ft.dropdown.Option("kih_ss", text="Knob-into-Hole + S354C/Y349C (-S-S-)"),
                ft.dropdown.Option("kih", text="Knob-into-Hole (T366W / T366S/L368A/Y407V)"),
                ft.dropdown.Option("charge_pair", text="Charge Pair (K392D/K409D vs D356K/D399K)"),
            ],
            dense=True,
            expand=True,
            on_select=lambda e: ctrl.wiz_state.__setitem__("mode", e.control.value),
        )

        def on_apply_heterodimer(e: ft.ControlEvent) -> None:
            if not ctrl.wiz_state["a"] or not ctrl.wiz_state["b"]:
                notify("Select two heavy chains.", is_error=True)
                return
            errs = ctrl.apply_heterodimer_pair(
                ctrl.wiz_state["a"],
                ctrl.wiz_state["b"],
                mode=ctrl.wiz_state["mode"] or "kih_ss",
            )
            if errs:
                notify("; ".join(errs), is_error=True)
            else:
                notify("Applied Fc heterodimerization & disulfide pairings.")
                render_workspace()

        def on_swap_crossmab(e: ft.ControlEvent) -> None:
            if not ctrl.wiz_state["a"] or not ctrl.wiz_state["b"]:
                notify("Select two chains for CrossMab swap.", is_error=True)
                return
            errs = ctrl.swap_crossmab_ch1_cl(ctrl.wiz_state["a"], ctrl.wiz_state["b"])
            if errs:
                notify("; ".join(errs), is_error=True)
            else:
                notify("Swapped CH1 <-> CL domains (CrossMab).")
                render_workspace()

        multispecific_col.controls = [
            ft.Row(spacing=10, controls=[dd_ca, dd_cb]),
            ft.Row(controls=[dd_mode]),
            ft.Row(
                spacing=10,
                controls=[
                    ft.Button("Apply Heterodimer Fc", expand=True, on_click=on_apply_heterodimer),
                    ft.OutlinedButton("Swap CH1 <-> CL (CrossMab)", expand=True, on_click=on_swap_crossmab),
                ],
            ),
        ]

    status_bar_txt = ft.Text("", size=11, color=pal["muted"])
    clipboard_badge_txt = ft.Text("", size=11, color="#F59E0B")
    undo_btn = ft.IconButton(
        icon=ft.Icons.UNDO,
        icon_size=18,
        tooltip="Undo (Cmd/Ctrl+Z)",
        disabled=True,
    )
    redo_btn = ft.IconButton(
        icon=ft.Icons.REDO,
        icon_size=18,
        tooltip="Redo (Cmd/Ctrl+Shift+Z)",
        disabled=True,
    )
    zoom_txt = ft.Text("100%", size=11, color=pal["muted"])
    recent_files_col = ft.Column(spacing=4)

    def _open_folder_for(path: str) -> None:
        folder = os.path.dirname(os.path.abspath(path))
        try:
            webbrowser.open(Path(folder).as_uri())
        except Exception as err:
            logger.debug("Could not open folder %s: %s", folder, err)

    def _refresh_recent_files_ui() -> None:
        items: List[ft.Control] = []
        for rpath in ctrl.recent_files[:5]:
            def _make_open_recent(p: str):
                def _on_open(ev: ft.ControlEvent) -> None:
                    try:
                        ctrl.load_ab(p)
                        ctrl.save_prefs()
                        notify(f"Loaded {os.path.basename(p)}")
                        render_workspace()
                    except Exception as err:
                        notify(str(err), is_error=True)

                return _on_open

            items.append(
                ft.TextButton(
                    text=os.path.basename(rpath),
                    tooltip=rpath,
                    on_click=_make_open_recent(rpath),
                )
            )
        recent_files_col.controls = items if items else [
            ft.Text("No recent projects yet.", size=11, color=pal["muted"])
        ]

    def render_workspace() -> None:
        project_name_tf.value = ctrl.project.name
        format_dd.value = ctrl.project.format_type
        update_page_title()
        undo_btn.disabled = not ctrl.can_undo
        redo_btn.disabled = not ctrl.can_redo
        zoom_txt.value = f"{int(round(ctrl.zoom_scale * 100))}%"
        total_doms = sum(len(c.domains) for c in ctrl.project.chains)
        dirty_str = " · Unsaved*" if ctrl.is_dirty else " · Saved"
        status_bar_txt.value = (
            f"Format: {ctrl.project.format_type} · {len(ctrl.project.chains)} chains · "
            f"{total_doms} domains · {len(ctrl.project.pairings)} pairings{dirty_str}"
        )
        if ctrl.copied_fasta and ctrl._clipboard_deadline > 0:
            rem = max(0, int(round(ctrl._clipboard_deadline - time.monotonic())))
            clipboard_badge_txt.value = f"Clipboard auto-clear: {rem}s" if rem > 0 else ""
            clear_clip_btn.visible = rem > 0
        else:
            clipboard_badge_txt.value = ""
            clear_clip_btn.visible = False
        _refresh_recent_files_ui()
        rebuild_canvas_only()
        rebuild_chain_strip()
        rebuild_pairings()
        rebuild_multispecific_tools()
        rebuild_inspector()
        if ctrl.is_dirty and ctrl.auto_recovery:
            try:
                ctrl.save_recovery_snapshot()
            except Exception as err:
                logger.debug("Recovery snapshot skipped: %s", err)
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
        new_fmt = e.control.value
        if not new_fmt or new_fmt == ctrl.project.format_type:
            return
        if not ctrl.is_dirty or not hasattr(page, "show_dialog"):
            ctrl.select_format(new_fmt)
            render_workspace()
            return

        def _confirm_fmt(ev: ft.ControlEvent) -> None:
            fmt_dlg.open = False
            ctrl.select_format(new_fmt)
            render_workspace()

        def _cancel_fmt(ev: ft.ControlEvent) -> None:
            fmt_dlg.open = False
            format_dd.value = ctrl.project.format_type
            page.update()

        fmt_dlg = ft.AlertDialog(
            title=ft.Text("Switch Format Template?"),
            content=ft.Text(
                f"Replace current topology with '{new_fmt}' template? (You can still Undo.)"
            ),
            actions=[
                ft.TextButton("Cancel", on_click=_cancel_fmt),
                ft.Button("Switch Template", on_click=_confirm_fmt),
            ],
        )
        page.show_dialog(fmt_dlg)

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
            ctrl.clear_recovery_snapshot()
            ctrl.save_prefs()
            notify(f"Loaded {os.path.basename(files[0].path)}")
            render_workspace()
        except Exception as err:
            if hasattr(page, "show_dialog"):
                err_dlg = ft.AlertDialog(
                    title=ft.Text("Load .ab Failed"),
                    content=ft.Text(
                        f"Could not load '{os.path.basename(files[0].path)}':\n{err}",
                        size=12,
                    ),
                    actions=[
                        ft.Button(
                            "OK",
                            on_click=lambda ev: (setattr(err_dlg, "open", False), page.update()),
                        )
                    ],
                )
                page.show_dialog(err_dlg)
            else:
                notify(str(err), is_error=True)

    async def on_save_ab(e: ft.ControlEvent) -> None:
        try:
            serialize_ab_project(ctrl.project)
        except Exception as err:
            notify(str(err), is_error=True)
            return
        dest = ctrl.current_file_path
        if not dest:
            dest = await file_picker.save_file(
                dialog_title="Save Antibody Project (.ab)",
                file_name=f"{ctrl.project.name.replace(' ', '_')}.ab",
                file_type=ft.FilePickerFileType.CUSTOM,
                allowed_extensions=["ab"],
            )
        if not dest:
            return
        try:
            ctrl.save_ab(dest)
            ctrl.clear_recovery_snapshot()
            ctrl.save_prefs()
            render_workspace()
            notify(
                f"Saved {os.path.basename(dest)}",
                action_label="Open Folder",
                on_action=lambda ev, p=dest: _open_folder_for(p),
            )
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_save_ab_as(e: ft.ControlEvent) -> None:
        try:
            serialize_ab_project(ctrl.project)
        except Exception as err:
            notify(str(err), is_error=True)
            return
        dest = await file_picker.save_file(
            dialog_title="Save Antibody Project As (.ab)",
            file_name=f"{ctrl.project.name.replace(' ', '_')}.ab",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["ab"],
        )
        if not dest:
            return
        try:
            ctrl.save_ab(dest)
            ctrl.clear_recovery_snapshot()
            ctrl.save_prefs()
            render_workspace()
            notify(
                f"Saved {os.path.basename(dest)}",
                action_label="Open Folder",
                on_action=lambda ev, p=dest: _open_folder_for(p),
            )
        except Exception as err:
            notify(str(err), is_error=True)

    fasta_header_dd = ft.Dropdown(
        label="FASTA Header",
        value="project_chain",
        width=165,
        dense=True,
        options=[
            ft.dropdown.Option("project_chain", text="Project | Chain"),
            ft.dropdown.Option("chain_only", text="Chain Only"),
            ft.dropdown.Option("detailed", text="Detailed (Domains)"),
        ],
    )

    async def on_export_fasta(e: ft.ControlEvent) -> None:
        _, pre_errs = format_project_fasta(
            ctrl.project, header_style=fasta_header_dd.value or "project_chain"
        )
        if pre_errs:
            notify("; ".join(pre_errs), is_error=True)
            return
        dest = await file_picker.save_file(
            dialog_title="Export Multi-FASTA",
            file_name="antibody_chains.fasta",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["fasta", "fa"],
        )
        if not dest:
            return
        try:
            save_fasta_file(
                dest, ctrl.project, header_style=fasta_header_dd.value or "project_chain"
            )
            notify(
                f"Exported {os.path.basename(dest)}",
                action_label="Open Folder",
                on_action=lambda ev, p=dest: _open_folder_for(p),
            )
        except Exception as err:
            notify(str(err), is_error=True)

    png_scale_dd = ft.Dropdown(
        label="PNG Scale",
        value="2.0",
        width=90,
        dense=True,
        options=[
            ft.dropdown.Option("1.0", text="1x"),
            ft.dropdown.Option("2.0", text="2x"),
            ft.dropdown.Option("3.0", text="3x"),
            ft.dropdown.Option("4.0", text="4x"),
        ],
    )
    png_bg_dd = ft.Dropdown(
        label="PNG Background",
        value="transparent",
        width=145,
        dense=True,
        options=[
            ft.dropdown.Option("transparent", text="Transparent"),
            ft.dropdown.Option("#FFFFFF", text="White (#FFFFFF)"),
            ft.dropdown.Option("#141413", text="Dark (#141413)"),
        ],
    )

    async def on_export_png(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Export PNG Schematic",
            file_name="antibody_schematic.png",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["png"],
        )
        if not dest:
            return
        try:
            scale_val = float(png_scale_dd.value or "2.0")
            bg_val = None if (png_bg_dd.value or "transparent") == "transparent" else png_bg_dd.value
            ctrl.save_png(dest, scale=scale_val, bg_color=bg_val)
            notify(
                f"Exported {os.path.basename(dest)} ({scale_val:g}x)",
                action_label="Open Folder",
                on_action=lambda ev, p=dest: _open_folder_for(p),
            )
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_export_genbank(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Export GenBank (.gb)",
            file_name="antibody_chains.gb",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["gb", "gbk"],
        )
        if not dest:
            return
        try:
            ctrl.save_genbank(dest)
            notify(
                f"Exported {os.path.basename(dest)}",
                action_label="Open Folder",
                on_action=lambda ev, p=dest: _open_folder_for(p),
            )
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_export_csv(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Export EU Mutation Diffs (.csv)",
            file_name="antibody_diffs.csv",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["csv"],
        )
        if not dest:
            return
        try:
            ctrl.save_csv(dest)
            notify(
                f"Exported {os.path.basename(dest)}",
                action_label="Open Folder",
                on_action=lambda ev, p=dest: _open_folder_for(p),
            )
        except Exception as err:
            notify(str(err), is_error=True)

    def on_copy_fasta(e: ft.ControlEvent) -> None:
        fasta_text, errs = ctrl.export_fasta()
        if errs:
            notify("; ".join(errs), is_error=True)
            return
        if ctrl.copy_text_secure(fasta_text, label="FASTA"):
            render_workspace()
            if hasattr(page, "run_task"):
                page.run_task(_clipboard_countdown_loop)
            notify("Multi-FASTA copied! (Will auto-clear in 30s)")
        else:
            notify("Error copying to clipboard.", is_error=True)

    theme_dd = ft.Dropdown(
        label="Theme Preset",
        value=ctrl.theme.name if ctrl.theme.name in ctrl.theme_presets else None,
        width=185,
        dense=True,
        options=[ft.dropdown.Option(k) for k in ctrl.theme_presets.keys()],
    )

    def on_theme_preset_select(e: ft.ControlEvent) -> None:
        if e.control.value:
            ctrl.apply_theme_preset(e.control.value)
            render_workspace()

    theme_dd.on_select = on_theme_preset_select

    async def on_load_theme(e: ft.ControlEvent) -> None:
        files = await file_picker.pick_files(
            dialog_title="Load Schematic Theme (.yaml)",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["yaml", "yml"],
        )
        if not files or not files[0].path:
            return
        try:
            loaded = ctrl.load_theme(files[0].path)
            theme_dd.value = loaded.name if loaded.name in ctrl.theme_presets else None
            notify(f"Loaded theme {os.path.basename(files[0].path)}")
            render_workspace()
        except Exception as err:
            notify(str(err), is_error=True)

    async def on_save_theme(e: ft.ControlEvent) -> None:
        dest = await file_picker.save_file(
            dialog_title="Save Schematic Theme (.yaml)",
            file_name="schematic_theme.yaml",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["yaml", "yml"],
        )
        if not dest:
            return
        try:
            ctrl.save_theme(dest)
            notify(f"Saved theme {os.path.basename(dest)}")
        except Exception as err:
            notify(str(err), is_error=True)

    def on_open_theme_editor(e: ft.ControlEvent) -> None:
        draft = copy.deepcopy(ctrl.theme)
        orient_dd = ft.Dropdown(
            label="Orientation",
            value=draft.orientation,
            options=[ft.dropdown.Option(o) for o in VALID_ORIENTATIONS],
            dense=True,
            width=150,
        )
        domain_dd = ft.Dropdown(
            label="Domain Type",
            value=VALID_DOMAIN_TYPES[0],
            options=[ft.dropdown.Option(dt) for dt in VALID_DOMAIN_TYPES],
            dense=True,
            width=140,
        )
        cur = draft.styles[VALID_DOMAIN_TYPES[0]]
        fill_tf = ft.TextField(label="Fill (#RRGGBB)", value=cur.fill_color, dense=True, width=130)
        border_tf = ft.TextField(label="Border (#RRGGBB)", value=cur.border_color, dense=True, width=130)
        text_tf = ft.TextField(label="Text (#RRGGBB)", value=cur.text_color, dense=True, width=130)
        bw_tf = ft.TextField(label="Border Width", value=str(cur.border_width), dense=True, width=110)
        w_tf = ft.TextField(label="Width", value=str(cur.width), dense=True, width=110)
        h_tf = ft.TextField(label="Height", value=str(cur.height), dense=True, width=110)
        shape_dd = ft.Dropdown(
            label="Shape",
            value=cur.shape,
            options=[ft.dropdown.Option(s) for s in VALID_SHAPES],
            dense=True,
            width=150,
        )

        def load_domain_fields(dtype: str) -> None:
            st = draft.styles[dtype]
            fill_tf.value = st.fill_color
            border_tf.value = st.border_color
            text_tf.value = st.text_color
            bw_tf.value = str(st.border_width)
            w_tf.value = str(st.width)
            h_tf.value = str(st.height)
            shape_dd.value = st.shape

        def save_current_domain_fields(dtype: str) -> None:
            for fld in (bw_tf, w_tf, h_tf):
                fld.error_text = None
            try:
                bw_val = float(bw_tf.value or 1.5)
            except ValueError:
                bw_tf.error_text = "Invalid number"
                raise
            try:
                w_val = float(w_tf.value or 52.0)
            except ValueError:
                w_tf.error_text = "Invalid number"
                raise
            try:
                h_val = float(h_tf.value or 88.0)
            except ValueError:
                h_tf.error_text = "Invalid number"
                raise
            draft.styles[dtype] = DomainVisualStyle(
                fill_color=(fill_tf.value or "#D97757").strip(),
                border_color=(border_tf.value or "#F09A7E").strip(),
                border_width=max(0.5, min(8.0, bw_val)),
                text_color=(text_tf.value or "#FAF9F5").strip(),
                width=max(32.0, min(240.0, w_val)),
                height=max(16.0, min(120.0, h_val)),
                shape=shape_dd.value or "rounded_rect",
            )

        active_dtype = {"val": VALID_DOMAIN_TYPES[0]}

        def on_domain_switch(ev: ft.ControlEvent) -> None:
            try:
                save_current_domain_fields(active_dtype["val"])
            except ValueError:
                page.update()
                return
            active_dtype["val"] = ev.control.value or VALID_DOMAIN_TYPES[0]
            load_domain_fields(active_dtype["val"])
            page.update()

        domain_dd.on_select = on_domain_switch

        def _auto_contrast_text(ev: ft.ControlEvent) -> None:
            hx = (fill_tf.value or "#D97757").strip().lstrip("#")
            if len(hx) == 6:
                try:
                    r, g, b = int(hx[0:2], 16), int(hx[2:4], 16), int(hx[4:6], 16)
                    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
                    text_tf.value = "#141413" if lum > 150 else "#FAF9F5"
                    page.update()
                except ValueError:
                    pass

        swatch_colors = [
            "#D97757", "#2563EB", "#0F766E", "#7C3AED",
            "#DC2626", "#D97706", "#4B5563", "#141413",
        ]

        def _make_swatch(col_hex: str) -> ft.Container:
            return ft.Container(
                width=20,
                height=20,
                bgcolor=col_hex,
                border_radius=4,
                border=ft.Border.all(1, "#9CA3AF"),
                tooltip=f"Set Fill to {col_hex}",
                on_click=lambda ev, c=col_hex: (
                    setattr(fill_tf, "value", c),
                    page.update(),
                ),
            )

        def on_apply_dialog(ev: ft.ControlEvent) -> None:
            try:
                save_current_domain_fields(active_dtype["val"])
                draft.orientation = orient_dd.value or "vertical"
                ctrl.set_theme(draft)
                theme_dd.value = None
                dlg.open = False
                render_workspace()
            except Exception as err:
                notify(str(err), is_error=True)
                page.update()

        def on_close_dialog(ev: ft.ControlEvent) -> None:
            dlg.open = False
            page.update()

        dlg = ft.AlertDialog(
            title=ft.Text("Customize Schematic Theme"),
            content=ft.Column(
                tight=True,
                spacing=8,
                controls=[
                    ft.Row(controls=[orient_dd, domain_dd]),
                    ft.Row(
                        spacing=6,
                        controls=[
                            ft.Text("Swatches:", size=11),
                            *[_make_swatch(c) for c in swatch_colors],
                            ft.OutlinedButton("Auto Contrast", on_click=_auto_contrast_text),
                        ],
                    ),
                    ft.Row(controls=[fill_tf, border_tf, text_tf]),
                    ft.Row(controls=[bw_tf, w_tf, h_tf, shape_dd]),
                ],
            ),
            actions=[
                ft.TextButton("Cancel", on_click=on_close_dialog),
                ft.Button("Apply Theme", on_click=on_apply_dialog),
            ],
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(dlg)

    canvas_holder = ft.Container(
        scale=ctrl.zoom_scale,
        content=canvas_overlay,
    )

    canvas_box = ft.Container(
        expand=True,
        bgcolor=pal["canvas_bg"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=6,
        content=ft.Column(
            scroll=ft.ScrollMode.AUTO,
            controls=[
                ft.Row(
                    scroll=ft.ScrollMode.AUTO,
                    controls=[canvas_holder],
                )
            ],
        ),
    )

    def _set_zoom(new_zoom: float) -> None:
        ctrl.zoom_scale = max(0.5, min(2.0, round(new_zoom, 2)))
        canvas_holder.scale = ctrl.zoom_scale
        render_workspace()

    def _on_undo(e: Any) -> None:
        if ctrl.undo():
            render_workspace()
            notify("Undid last action.")

    def _on_redo(e: Any) -> None:
        if ctrl.redo():
            render_workspace()
            notify("Redid action.")

    undo_btn.on_click = _on_undo
    redo_btn.on_click = _on_redo

    def _on_toggle_snap(e: ft.ControlEvent) -> None:
        ctrl.snap_to_grid = not ctrl.snap_to_grid
        ctrl.save_prefs()
        e.control.text = f"Grid: {'ON' if ctrl.snap_to_grid else 'OFF'}"
        page.update()

    def _on_toggle_orient(e: ft.ControlEvent) -> None:
        new_theme = copy.deepcopy(ctrl.theme)
        new_theme.orientation = (
            "horizontal" if ctrl.theme.orientation == "vertical" else "vertical"
        )
        ctrl.set_theme(new_theme)
        e.control.text = "Orient: H" if new_theme.orientation == "horizontal" else "Orient: V"
        render_workspace()

    def open_batch_variants_modal(e: Any = None) -> None:
        if not ctrl.selected_chain_id or not ctrl.selected_domain_id:
            notify("Select a constant domain first.", is_error=True)
            return
        batch_tf = ft.TextField(
            label="Batch Variants (comma or newline separated, e.g. L234A/L235A, M252Y/S254T/T256E)",
            multiline=True,
            min_lines=4,
            max_lines=6,
            value="M252Y/S254T/T256E\nL234A/L235A",
        )
        batch_status = ft.Text("", size=11, color=pal["muted"])

        async def _export_batch(ev: ft.ControlEvent) -> None:
            dest = await file_picker.save_file(
                dialog_title="Export Batch Variants (.csv)",
                file_name="batch_variants.csv",
                file_type=ft.FilePickerFileType.CUSTOM,
                allowed_extensions=["csv"],
            )
            if not dest:
                return
            try:
                cnt, errs = ctrl.export_batch_csv(
                    ctrl.selected_chain_id or "",
                    ctrl.selected_domain_id or "",
                    batch_tf.value or "",
                    dest,
                )
                if errs:
                    batch_status.value = f"Exported {cnt} row(s) with warnings: {'; '.join(errs[:2])}"
                    batch_status.color = "#F59E0B"
                else:
                    batch_dlg.open = False
                    notify(
                        f"Exported {cnt} batch variants to {os.path.basename(dest)}",
                        action_label="Open Folder",
                        on_action=lambda _, p=dest: _open_folder_for(p),
                    )
                page.update()
            except Exception as err:
                batch_status.value = str(err)
                batch_status.color = "#DC2626"
                page.update()

        def _close_batch(ev: ft.ControlEvent) -> None:
            batch_dlg.open = False
            page.update()

        batch_dlg = ft.AlertDialog(
            title=ft.Text("Batch Variant Generator"),
            content=ft.Container(
                width=520,
                content=ft.Column(tight=True, spacing=8, controls=[batch_tf, batch_status]),
            ),
            actions=[
                ft.TextButton("Close", on_click=_close_batch),
                ft.Button("Export Batch CSV...", on_click=_export_batch),
            ],
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(batch_dlg)

    def open_help_modal(e: Any = None) -> None:
        help_lines = [
            "Cmd/Ctrl + S : Save project (.ab)",
            "Cmd/Ctrl + O : Load project (.ab)",
            "Cmd/Ctrl + Z : Undo",
            "Cmd/Ctrl + Shift + Z : Redo",
            "Cmd/Ctrl + K : Command Palette",
            "Double-click constant domain : Edit Fc Mutations",
            "Drag domain on canvas : Pin custom 2D position",
            "Click -S-S- badge : Remove disulfide pairing",
        ]
        help_dlg = ft.AlertDialog(
            title=ft.Text("Keyboard Shortcuts & Tips"),
            content=ft.Column(
                tight=True,
                spacing=4,
                controls=[ft.Text(line, size=12) for line in help_lines],
            ),
            actions=[
                ft.Button(
                    "Close",
                    on_click=lambda ev: (setattr(help_dlg, "open", False), page.update()),
                )
            ],
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(help_dlg)

    def open_command_palette(e: Any = None) -> None:
        cmd_list = ft.Column(spacing=4, scroll=ft.ScrollMode.AUTO, height=220)
        commands = [
            ("Save Project (.ab)", lambda: page.run_task(on_save_ab, None) if hasattr(page, "run_task") else None),
            ("Undo", lambda: _on_undo(None)),
            ("Redo", lambda: _on_redo(None)),
            ("Add Chain", lambda: (ctrl.add_chain(), render_workspace())),
            ("Reset 2D Layout", lambda: (ctrl.reset_layout(), render_workspace())),
            ("Auto-Detect Pairings", lambda: (ctrl.auto_detect_pairings(), render_workspace())),
            ("Edit Fc Mutations...", lambda: open_fc_mutations_modal(ctrl.selected_chain_id) if ctrl.selected_chain_id else None),
            ("Batch Variant Generator...", lambda: open_batch_variants_modal()),
            ("Keyboard Shortcuts & Help", lambda: open_help_modal()),
        ]

        def _refresh_cmds(q: str = "") -> None:
            ql = q.strip().lower()
            btns: List[ft.Control] = []
            for label, fn in commands:
                if ql and ql not in label.lower():
                    continue

                def _make_run(callback):
                    def _run(ev: ft.ControlEvent) -> None:
                        cmd_dlg.open = False
                        page.update()
                        callback()

                    return _run

                btns.append(
                    ft.TextButton(
                        text=label,
                        on_click=_make_run(fn),
                    )
                )
            cmd_list.controls = btns

        _refresh_cmds()
        cmd_tf = ft.TextField(
            label="Type a command...",
            autofocus=True,
            dense=True,
            on_change=lambda ev: (_refresh_cmds(ev.control.value or ""), page.update()),
        )
        cmd_dlg = ft.AlertDialog(
            title=ft.Text("Command Palette (Cmd/Ctrl+K)"),
            content=ft.Container(
                width=420,
                content=ft.Column(tight=True, spacing=8, controls=[cmd_tf, cmd_list]),
            ),
            actions=[
                ft.TextButton(
                    "Close",
                    on_click=lambda ev: (setattr(cmd_dlg, "open", False), page.update()),
                )
            ],
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(cmd_dlg)

    studio_toolbar = ft.Row(
        spacing=4,
        wrap=True,
        controls=[
            ft.IconButton(
                icon=ft.Icons.SAVE,
                icon_size=18,
                tooltip="Save (.ab) (Cmd/Ctrl+S)",
                on_click=on_save_ab,
            ),
            undo_btn,
            redo_btn,
            ft.VerticalDivider(width=1),
            ft.IconButton(
                icon=ft.Icons.ZOOM_OUT,
                icon_size=16,
                tooltip="Zoom Out",
                on_click=lambda e: _set_zoom(ctrl.zoom_scale - 0.15),
            ),
            zoom_txt,
            ft.IconButton(
                icon=ft.Icons.ZOOM_IN,
                icon_size=16,
                tooltip="Zoom In",
                on_click=lambda e: _set_zoom(ctrl.zoom_scale + 0.15),
            ),
            ft.OutlinedButton(
                f"Grid: {'ON' if ctrl.snap_to_grid else 'OFF'}",
                on_click=_on_toggle_snap,
            ),
            ft.OutlinedButton(
                "Orient: H" if ctrl.theme.orientation == "horizontal" else "Orient: V",
                on_click=_on_toggle_orient,
            ),
            ft.OutlinedButton("Batch...", on_click=open_batch_variants_modal),
            ft.IconButton(
                icon=ft.Icons.SEARCH,
                icon_size=18,
                tooltip="Command Palette (Cmd/Ctrl+K)",
                on_click=open_command_palette,
            ),
            ft.IconButton(
                icon=ft.Icons.HELP_OUTLINE,
                icon_size=18,
                tooltip="Shortcuts & Help",
                on_click=open_help_modal,
            ),
        ],
    )

    palette_blocks: List[ft.Control] = []
    for dt in VALID_DOMAIN_TYPES:
        p_col = DOMAIN_STYLE_MAP[dt][2]
        palette_blocks.append(
            ft.Draggable(
                group="domain_palette",
                data=f"palette:{dt}",
                content=ft.Container(
                    width=68,
                    padding=ft.Padding.symmetric(horizontal=6, vertical=6),
                    bgcolor=pal["elevated"],
                    border=ft.Border.all(1.5, p_col),
                    border_radius=6,
                    alignment=ft.Alignment.CENTER,
                    tooltip=f"Drag {dt} onto a chain",
                    content=ft.Text(dt, size=11, weight=ft.FontWeight.BOLD, color=pal["text"]),
                ),
            )
        )

    palette_panel = ft.Container(
        width=86,
        padding=8,
        bgcolor=pal["panel"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=8,
        content=ft.Column(
            spacing=6,
            scroll=ft.ScrollMode.AUTO,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            controls=[
                ft.Text("Palette", size=11, weight=ft.FontWeight.BOLD),
                *palette_blocks,
            ],
        ),
    )

    side_panel_tabs = ft.Tabs(
        selected_index=0,
        length=2,
        expand=True,
        content=ft.Column(
            expand=True,
            controls=[
                ft.TabBar(
                    tabs=[
                        ft.Tab(label="Inspector"),
                        ft.Tab(label="Pairings"),
                    ]
                ),
                ft.TabBarView(
                    expand=True,
                    controls=[
                        inspector_col,
                        ft.Column(
                            scroll=ft.ScrollMode.AUTO,
                            spacing=8,
                            controls=[pairing_controls_col],
                        ),
                    ],
                ),
            ],
        ),
    )

    studio_side_panel = ft.Container(
        width=360,
        padding=10,
        bgcolor=pal["panel"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=8,
        content=side_panel_tabs,
    )

    studio_main_col = ft.Container(
        expand=True,
        padding=8,
        bgcolor=pal["panel"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=8,
        content=ft.Column(
            expand=True,
            spacing=8,
            controls=[
                studio_toolbar,
                canvas_box,
                ft.Text(
                    "Polypeptide Chains (N -> C)",
                    size=12,
                    weight=ft.FontWeight.BOLD,
                ),
                ft.Container(height=168, content=chain_strip_col),
            ],
        ),
    )

    studio_panel = ft.Row(
        expand=True,
        spacing=8,
        vertical_alignment=ft.CrossAxisAlignment.STRETCH,
        controls=[palette_panel, studio_main_col, studio_side_panel],
    )

    project_cards: List[ft.Container] = []

    def _make_project_card(title: str, controls: List[ft.Control]) -> ft.Container:
        card = ft.Container(
            padding=16,
            bgcolor=pal["panel"],
            border=ft.Border.all(1, pal["border"]),
            border_radius=8,
            content=ft.Column(
                spacing=12,
                controls=[
                    ft.Text(title, size=15, weight=ft.FontWeight.BOLD),
                    *controls,
                ],
            ),
        )
        project_cards.append(card)
        return card

    def on_clone_selected_chain(e: ft.ControlEvent) -> None:
        if not ctrl.selected_chain_id:
            notify("Select a chain to clone.", is_error=True)
            return
        try:
            cloned = ctrl.clone_chain(ctrl.selected_chain_id)
            notify(f"Cloned chain as {cloned.name}")
            render_workspace()
        except Exception as err:
            notify(str(err), is_error=True)

    setup_card = _make_project_card(
        "Project Setup",
        [
            ft.Row(spacing=10, controls=[project_name_tf, format_dd]),
            ft.Row(
                spacing=10,
                controls=[
                    ft.Button(
                        "Add Chain",
                        icon=ft.Icons.ADD,
                        expand=True,
                        on_click=lambda e: (ctrl.add_chain(), render_workspace()),
                    ),
                    ft.OutlinedButton(
                        "Clone Selected Chain",
                        icon=ft.Icons.COPY_ALL,
                        expand=True,
                        on_click=on_clone_selected_chain,
                    ),
                ],
            ),
        ],
    )

    multispecific_card = _make_project_card(
        "Multispecific & Heterodimer Builder",
        [multispecific_col],
    )

    ab_files_card = _make_project_card(
        "Project Files (.ab)",
        [
            ft.Row(
                spacing=10,
                controls=[
                    ft.OutlinedButton(
                        "Load .ab",
                        icon=ft.Icons.UPLOAD_FILE,
                        expand=True,
                        on_click=on_load_ab,
                    ),
                    ft.Button(
                        "Save .ab",
                        icon=ft.Icons.SAVE,
                        expand=True,
                        on_click=on_save_ab,
                    ),
                    ft.OutlinedButton(
                        "Save As...",
                        icon=ft.Icons.SAVE_AS,
                        expand=True,
                        on_click=on_save_ab_as,
                    ),
                ],
            ),
            ft.Divider(),
            ft.Text("Recent Projects", size=12, weight=ft.FontWeight.W_600),
            recent_files_col,
        ],
    )

    export_card = _make_project_card(
        "Sequence & Schematic Export",
        [
            ft.Row(
                spacing=10,
                controls=[
                    ft.OutlinedButton(
                        "Export FASTA",
                        icon=ft.Icons.DOWNLOAD,
                        expand=True,
                        on_click=on_export_fasta,
                    ),
                    fasta_header_dd,
                    ft.OutlinedButton(
                        "Export GenBank (.gb)",
                        icon=ft.Icons.DESCRIPTION,
                        expand=True,
                        on_click=on_export_genbank,
                    ),
                ],
            ),
            ft.Row(
                spacing=10,
                controls=[
                    ft.OutlinedButton(
                        "Export CSV (.csv)",
                        icon=ft.Icons.TABLE_CHART,
                        expand=True,
                        on_click=on_export_csv,
                    ),
                    ft.OutlinedButton(
                        "Export PNG",
                        icon=ft.Icons.IMAGE,
                        expand=True,
                        on_click=on_export_png,
                    ),
                    png_scale_dd,
                    png_bg_dd,
                ],
            ),
            ft.Row(
                controls=[
                    ft.Button(
                        "Copy FASTA (30s Auto-Clear)",
                        icon=ft.Icons.CONTENT_COPY,
                        expand=True,
                        on_click=on_copy_fasta,
                    ),
                ],
            ),
        ],
    )

    theme_files_card = _make_project_card(
        "Schematic Theme Files (.yaml)",
        [
            ft.Text(
                "2D schematic theme files (domain colors, shapes, sizes, orientation). Saved and loaded independently from Light/Dark UI mode and .ab project files.",
                size=12,
                color=pal["muted"],
            ),
            ft.Row(
                spacing=10,
                controls=[
                    ft.OutlinedButton(
                        "Load Theme .yaml",
                        icon=ft.Icons.UPLOAD_FILE,
                        expand=True,
                        on_click=on_load_theme,
                    ),
                    ft.OutlinedButton(
                        "Save Theme .yaml",
                        icon=ft.Icons.SAVE,
                        expand=True,
                        on_click=on_save_theme,
                    ),
                ],
            ),
        ],
    )

    project_panel = ft.Container(
        expand=True,
        padding=12,
        bgcolor=pal["bg"],
        content=ft.Column(
            scroll=ft.ScrollMode.AUTO,
            controls=[
                ft.Row(
                    vertical_alignment=ft.CrossAxisAlignment.START,
                    spacing=16,
                    controls=[
                        ft.Column(
                            expand=True,
                            spacing=16,
                            controls=[setup_card, multispecific_card],
                        ),
                        ft.Column(
                            expand=True,
                            spacing=16,
                            controls=[ab_files_card, export_card],
                        ),
                    ],
                )
            ],
        ),
    )

    def on_theme_mode_toggle(e: ft.ControlEvent) -> None:
        new_mode = "light" if bool(e.control.value) else "dark"
        ctrl.set_theme_mode(new_mode)
        new_pal = _ui_palette(new_mode)
        page.theme_mode = (
            ft.ThemeMode.LIGHT if new_mode == "light" else ft.ThemeMode.DARK
        )
        page.bgcolor = new_pal["bg"]
        canvas_box.bgcolor = new_pal["canvas_bg"]
        canvas_box.border = ft.Border.all(1, new_pal["border"])
        palette_panel.bgcolor = new_pal["panel"]
        palette_panel.border = ft.Border.all(1, new_pal["border"])
        studio_main_col.bgcolor = new_pal["panel"]
        studio_main_col.border = ft.Border.all(1, new_pal["border"])
        studio_side_panel.bgcolor = new_pal["panel"]
        studio_side_panel.border = ft.Border.all(1, new_pal["border"])
        project_panel.bgcolor = new_pal["bg"]
        for card in project_cards:
            card.bgcolor = new_pal["panel"]
            card.border = ft.Border.all(1, new_pal["border"])
        settings_panel.bgcolor = new_pal["panel"]
        settings_panel.border = ft.Border.all(1, new_pal["border"])
        status_bar_row.bgcolor = new_pal["elevated"]
        status_bar_row.border = ft.Border.all(1, new_pal["border"])
        status_bar_txt.color = new_pal["muted"]
        ctrl.save_prefs()
        render_workspace()

    light_mode_switch = ft.Switch(
        label="Light Mode (uncheck for Dark Mode)",
        value=(ctrl.theme_mode == "light"),
        on_change=on_theme_mode_toggle,
    )

    style_domain_gap_slider = ft.Slider(
        min=4.0,
        max=80.0,
        divisions=38,
        value=ctrl.theme.domain_spacing,
        label="Domain Gap: {value}px",
        expand=True,
        on_change=lambda e: (
            ctrl.set_spacing(domain_spacing=float(e.control.value)),
            render_workspace(),
        ),
    )
    style_chain_gap_slider = ft.Slider(
        min=60.0,
        max=260.0,
        divisions=40,
        value=ctrl.theme.chain_spacing,
        label="Chain Gap: {value}px",
        expand=True,
        on_change=lambda e: (
            ctrl.set_spacing(chain_spacing=float(e.control.value)),
            render_workspace(),
        ),
    )

    settings_panel = ft.Container(
        expand=True,
        padding=16,
        bgcolor=pal["panel"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=8,
        content=ft.Column(
            scroll=ft.ScrollMode.AUTO,
            spacing=12,
            controls=[
                ft.Text("App Appearance & Schematic Style", size=16, weight=ft.FontWeight.BOLD),
                light_mode_switch,
                ft.Divider(),
                ft.Text("Schematic Theme Preset", size=14, weight=ft.FontWeight.W_600),
                ft.Row(
                    wrap=True,
                    spacing=8,
                    controls=[
                        theme_dd,
                        ft.OutlinedButton("Customize Theme", on_click=on_open_theme_editor),
                    ],
                ),
                ft.Text("Schematic Spacing", size=14, weight=ft.FontWeight.W_600),
                ft.Row(
                    spacing=12,
                    controls=[
                        ft.Text("Domain Gap", size=12),
                        style_domain_gap_slider,
                        ft.Text("Chain Gap", size=12),
                        style_chain_gap_slider,
                    ],
                ),
                theme_files_card,
            ],
        ),
    )

    main_tabs = ft.Tabs(
        selected_index=0,
        length=3,
        expand=True,
        content=ft.Column(
            expand=True,
            controls=[
                ft.TabBar(
                    tabs=[
                        ft.Tab(label="Studio", icon=ft.Icons.SCHEMA),
                        ft.Tab(label="Project", icon=ft.Icons.ACCOUNT_TREE),
                        ft.Tab(label="Style", icon=ft.Icons.PALETTE),
                    ]
                ),
                ft.TabBarView(
                    expand=True,
                    controls=[
                        studio_panel,
                        project_panel,
                        settings_panel,
                    ],
                ),
            ],
        ),
    )

    def _on_clear_clip_now(ev: Any) -> None:
        ctrl.clear_clipboard()
        render_workspace()
        notify("Clipboard cleared.")

    clear_clip_btn = ft.TextButton(
        "Clear Now",
        visible=False,
        on_click=_on_clear_clip_now,
    )

    async def _clipboard_countdown_loop() -> None:
        while not ctrl._shutdown_flag and ctrl.copied_fasta and ctrl._clipboard_deadline > 0:
            rem = max(0, int(round(ctrl._clipboard_deadline - time.monotonic())))
            clipboard_badge_txt.value = f"Clipboard auto-clear: {rem}s" if rem > 0 else ""
            clear_clip_btn.visible = rem > 0
            page.update()
            if rem <= 0:
                break
            await asyncio.sleep(1.0)
        if not ctrl.copied_fasta:
            clipboard_badge_txt.value = ""
            clear_clip_btn.visible = False
            page.update()

    status_bar_row = ft.Container(
        padding=ft.Padding.symmetric(horizontal=10, vertical=4),
        bgcolor=pal["elevated"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=6,
        content=ft.Row(
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            controls=[
                status_bar_txt,
                ft.Row(spacing=4, tight=True, controls=[clipboard_badge_txt, clear_clip_btn]),
            ],
        ),
    )

    def _on_keyboard_event(e: ft.KeyboardEvent) -> None:
        mod = bool(getattr(e, "meta", False) or getattr(e, "ctrl", False))
        if not mod:
            return
        key = (getattr(e, "key", "") or "").upper()
        shift = bool(getattr(e, "shift", False))
        if key == "S":
            if hasattr(page, "run_task"):
                page.run_task(on_save_ab, None)
        elif key == "O":
            if hasattr(page, "run_task"):
                page.run_task(on_load_ab, None)
        elif key == "Z" and shift:
            _on_redo(None)
        elif key == "Z":
            _on_undo(None)
        elif key == "K":
            open_command_palette()

    if hasattr(page, "on_keyboard_event"):
        page.on_keyboard_event = _on_keyboard_event

    page.add(
        ft.Column(
            expand=True,
            spacing=6,
            controls=[main_tabs, status_bar_row],
        )
    )

    render_workspace()
    if ctrl.auto_recovery and ctrl.has_recovery_snapshot() and hasattr(page, "show_dialog"):
        def _restore_rec(ev: Any) -> None:
            try:
                ctrl.load_recovery_snapshot()
                render_workspace()
                notify("Restored unsaved recovery snapshot.")
            except Exception as err:
                notify(str(err), is_error=True)
                page.update()

        notify(
            "Unsaved recovery snapshot found.",
            action_label="Restore",
            on_action=_restore_rec,
        )
    return ctrl


def main() -> None:
    """Launch the Fc Engineering Studio Flet desktop application."""
    controller = AntibodyStudioController()
    controller.load_prefs()
    controller.auto_recovery = True
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
