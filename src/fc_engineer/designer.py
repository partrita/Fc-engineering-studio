"""Modular antibody domain assembly, bispecific templates, mutation routing, and 2D schematic scene model."""

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from fc_engineer.config import (
    COMMON_MUTATIONS,
    DEFAULT_LINKERS,
    LIGHT_CONSTANT_DOMAINS,
    SEQUENCES,
    VARIABLE_DOMAINS,
)
from fc_engineer.core import (
    CONSTANT_DOMAIN_EU_RANGES,
    analyze_developability,
    apply_domain_mutations,
    diff_sequences,
    parse_mutation_list,
    slice_constant_domain,
)

MAX_CHAINS = 16
MAX_DOMAINS_PER_CHAIN = 32
MAX_PAIRINGS = 64
MAX_FILE_CUSTOM_MUTATIONS = 100
MAX_CUSTOM_SEQ_LENGTH = 2000
MAX_EXPORT_SEQUENCE_LENGTH = 10000
MAX_METADATA_LENGTH = 200
MAX_HISTORY_ENTRIES = 50
MAX_COORD = 10000.0
SCENE_PADDING: float = 24.0

VALID_DOMAIN_TYPES = ("VH", "VL", "CH1", "CL", "Hinge", "CH2", "CH3", "G4S")
NON_CONSTANT_DOMAIN_TYPES = ("VH", "VL", "CL", "G4S")
VALID_PAIRING_TYPES = ("VH_VL", "CH1_CL", "HINGE_DISULFIDE", "CH3_DIMER")
VALID_SHAPES = ("rect", "rounded_rect", "capsule", "ellipse", "hexagon")
VALID_ORIENTATIONS = ("vertical", "horizontal")
VALID_FORMAT_TYPES = (
    "igg",
    "crossmab",
    "igg_scfv",
    "scfv_fc",
    "dvd_ig",
    "bite",
    "custom",
)

ALLOWED_TEMPLATE_KEYS_BY_TYPE: Dict[str, Tuple[str, ...]] = {
    "VH": ("VH_A", "VH_B"),
    "VL": ("VL_A", "VL_B"),
    "CL": ("CL_kappa", "CL_lambda"),
    "G4S": ("G4S_3", "G4S_2", "G4S_1"),
}

DEFAULT_TEMPLATE_KEYS: Dict[str, str] = {
    "VH": "VH_A",
    "VL": "VL_A",
    "CL": "CL_kappa",
    "G4S": "G4S_3",
}

DOMAIN_STYLE_MAP: Dict[str, Tuple[float, float, str]] = {
    "VH": (96.0, 44.0, "#1D4ED8"),
    "VL": (96.0, 44.0, "#0F766E"),
    "CH1": (96.0, 44.0, "#334155"),
    "CL": (96.0, 44.0, "#475569"),
    "Hinge": (64.0, 36.0, "#6D28D9"),
    "CH2": (96.0, 44.0, "#0369A1"),
    "CH3": (96.0, 44.0, "#1E40AF"),
    "G4S": (64.0, 36.0, "#374151"),
}

_VALID_SEQ_RE = re.compile(r"[A-Z]+\Z")

_DEVELOPABILITY_CATEGORIES: Tuple[str, ...] = (
    "N-glycosylation motif (N-X-S/T)",
    "Deamidation-prone NG",
    "Isomerization-prone DG",
    "Oxidation-prone Met (M252/M428 = FcRn sites)",
)


@dataclass
class DomainBlock:
    """Single modular domain inside an expressible polypeptide chain."""

    id: str
    domain_type: str
    isotype: str = "igg1"
    allotype: str = ""
    template_key: str = ""
    custom_sequence: str = ""
    x: float = 0.0
    y: float = 0.0
    pinned: bool = False
    applied_presets: List[str] = field(default_factory=list)
    custom_mutations: List[str] = field(default_factory=list)
    mutation_history: List[str] = field(default_factory=list)


@dataclass
class PolypeptideChain:
    """Expressible N-to-C ordered polypeptide chain."""

    id: str
    name: str
    domains: List[DomainBlock] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0


@dataclass
class InterChainPairing:
    """Disulfide or non-covalent pairing interaction between two domains."""

    id: str
    chain_a_id: str
    domain_a_id: str
    chain_b_id: str
    domain_b_id: str
    pairing_type: str


@dataclass
class AntibodyProject:
    """Root container for a multi-chain antibody / bispecific engineering workspace."""

    name: str
    format_type: str = "custom"
    chains: List[PolypeptideChain] = field(default_factory=list)
    pairings: List[InterChainPairing] = field(default_factory=list)
    file_custom_mutations: List[Dict[str, str]] = field(default_factory=list)


@dataclass
class DomainVisualStyle:
    """Visual styling attributes for a single antibody domain_type."""

    fill_color: str
    border_color: str
    border_width: float
    text_color: str
    width: float
    height: float
    shape: str = "rounded_rect"


@dataclass
class SchematicTheme:
    """Complete schematic styling, orientation, and spacing configuration."""

    name: str
    orientation: str = "vertical"
    domain_spacing: float = 18.0
    chain_spacing: float = 116.0
    styles: Dict[str, DomainVisualStyle] = field(default_factory=dict)


DEFAULT_SCHEMATIC_THEME = SchematicTheme(
    name="Claude Code",
    orientation="vertical",
    domain_spacing=18.0,
    chain_spacing=116.0,
    styles={
        "VH": DomainVisualStyle("#D97757", "#F09A7E", 1.5, "#FAF9F5", 52.0, 88.0, "rounded_rect"),
        "VL": DomainVisualStyle("#E09F5A", "#F2BC82", 1.5, "#141413", 52.0, 88.0, "rounded_rect"),
        "CH1": DomainVisualStyle("#B85C3E", "#D97757", 1.5, "#FAF9F5", 52.0, 88.0, "rounded_rect"),
        "CL": DomainVisualStyle("#C48B71", "#DFB099", 1.5, "#FAF9F5", 52.0, 88.0, "rounded_rect"),
        "Hinge": DomainVisualStyle("#6E6A61", "#9E988C", 1.5, "#FAF9F5", 38.0, 52.0, "capsule"),
        "CH2": DomainVisualStyle("#5C8A72", "#84B39B", 1.5, "#FAF9F5", 52.0, 88.0, "rounded_rect"),
        "CH3": DomainVisualStyle("#3D6B58", "#639681", 1.5, "#FAF9F5", 52.0, 88.0, "rounded_rect"),
        "G4S": DomainVisualStyle("#524F48", "#7D786D", 1.5, "#ECE9E0", 38.0, 52.0, "capsule"),
    },
)


def compute_hexagon_vertices(
    x: float, y: float, width: float, height: float
) -> List[Tuple[float, float]]:
    """Compute 6 polygon vertices for a hexagon inscribed in [x, x+width] x [y, y+height]."""
    cx = x + width / 2.0
    cy = y + height / 2.0
    if height >= width:
        inset_y = height * 0.22
        return [
            (cx, y),
            (x + width, y + inset_y),
            (x + width, y + height - inset_y),
            (cx, y + height),
            (x, y + height - inset_y),
            (x, y + inset_y),
        ]
    inset_x = width * 0.22
    return [
        (x + inset_x, y),
        (x + width - inset_x, y),
        (x + width, cy),
        (x + width - inset_x, y + height),
        (x + inset_x, y + height),
        (x, cy),
    ]


@dataclass
class SchematicNode:
    """Positioned domain node consumed by the Flet canvas and Pillow PNG exporter."""

    domain_id: str
    chain_id: str
    chain_name: str
    is_n_terminus: bool
    label: str
    domain_type: str
    x: float
    y: float
    width: float
    height: float
    fill_color: str
    badge_text: str
    error_flag: bool
    border_color: str = "#0F172A"
    border_width: float = 1.5
    text_color: str = "#FFFFFF"
    shape: str = "rounded_rect"


@dataclass
class SchematicEdge:
    """Covalent N-to-C backbone segment between consecutive domains on a chain."""

    chain_id: str
    from_domain_id: str
    to_domain_id: str
    x1: float
    y1: float
    x2: float
    y2: float
    edge_kind: str


@dataclass
class SchematicPairingEdge:
    """Visual edge representing an inter-chain or intra-chain pairing."""

    pairing_id: str
    chain_a_id: str
    domain_a_id: str
    chain_b_id: str
    domain_b_id: str
    pairing_type: str
    x1: float
    y1: float
    x2: float
    y2: float


@dataclass
class SchematicScene:
    """UI-agnostic 2D scene representation with padded bounding box."""

    nodes: List[SchematicNode]
    edges: List[SchematicEdge]
    pairing_edges: List[SchematicPairingEdge]
    bounds: Tuple[float, float, float, float]


def _builtin_preset_map() -> Dict[str, str]:
    return {
        item["label"]: item["value"]
        for item in COMMON_MUTATIONS
        if isinstance(item, dict) and "label" in item and "value" in item
    }


def build_preset_map(project: Optional[AntibodyProject] = None) -> Dict[str, str]:
    """Merge built-in COMMON_MUTATIONS with project-scoped file_custom_mutations."""
    merged = _builtin_preset_map()
    if project is not None:
        for item in project.file_custom_mutations:
            if isinstance(item, dict) and "label" in item and "value" in item:
                merged[item["label"]] = item["value"]
    return merged


def get_effective_domain_mutations(
    domain: DomainBlock, preset_map: Optional[Dict[str, str]] = None
) -> Tuple[List[str], List[str]]:
    """Expand preset labels and custom point mutations for a constant domain."""
    resolved_map = preset_map if preset_map is not None else _builtin_preset_map()

    if domain.domain_type not in CONSTANT_DOMAIN_EU_RANGES:
        if domain.applied_presets or domain.custom_mutations:
            return [], [
                f"Domain '{domain.domain_type}' is outside CONSTANT_DOMAIN_EU_RANGES and cannot carry EU mutations."
            ]
        return [], []

    eu_start, eu_end = CONSTANT_DOMAIN_EU_RANGES[domain.domain_type]
    tokens: List[str] = []
    seen = set()
    errors: List[str] = []

    for label in domain.applied_presets:
        if label not in resolved_map:
            errors.append(f"Unknown mutation preset '{label}'.")
            continue
        parsed, parse_errs = parse_mutation_list(resolved_map[label])
        if parse_errs:
            errors.extend(parse_errs)
            continue
        for wt, pos, mut in parsed:
            if eu_start <= pos <= eu_end:
                norm = f"{wt}{pos}{mut}"
                if norm not in seen:
                    seen.add(norm)
                    tokens.append(norm)

    for raw_tok in domain.custom_mutations:
        parsed, parse_errs = parse_mutation_list(raw_tok)
        if parse_errs:
            errors.extend(parse_errs)
            continue
        if len(parsed) != 1:
            errors.append(f"Invalid custom mutation token '{raw_tok}'.")
            continue
        wt, pos, mut = parsed[0]
        if not (eu_start <= pos <= eu_end):
            errors.append(
                f"Custom mutation '{raw_tok}' is outside EU range {eu_start}-{eu_end} for {domain.domain_type}."
            )
            continue
        norm = f"{wt}{pos}{mut}"
        if norm not in seen:
            seen.add(norm)
            tokens.append(norm)

    if errors:
        return [], errors
    return tokens, []


def resolve_domain_base_sequence(domain: DomainBlock) -> Tuple[str, List[str]]:
    """Resolve the wild-type or custom base amino acid sequence for a DomainBlock."""
    dtype = domain.domain_type
    if dtype not in VALID_DOMAIN_TYPES:
        return "", [f"Unsupported domain type '{dtype}'."]

    if dtype in CONSTANT_DOMAIN_EU_RANGES:
        if domain.custom_sequence:
            return "", [
                f"Custom sequence override is not allowed on constant domain '{dtype}'."
            ]
        iso_map = SEQUENCES.get(domain.isotype)
        if not isinstance(iso_map, dict) or not iso_map:
            return "", [f"Unknown or empty isotype '{domain.isotype}'."]
        allotype = domain.allotype or next(iter(iso_map))
        full_seq = iso_map.get(allotype)
        if not isinstance(full_seq, str) or not full_seq:
            return "", [
                f"Unknown allotype '{allotype}' for isotype '{domain.isotype}'."
            ]
        sliced = slice_constant_domain(full_seq, dtype, domain.isotype)
        if not sliced:
            return "", [f"Failed to slice {dtype} for {domain.isotype} ({allotype})."]
        return sliced, []

    if dtype == "CL":
        if domain.custom_sequence:
            return "", ["Custom sequence override is not allowed on constant domain 'CL'."]
        tkey = domain.template_key or DEFAULT_TEMPLATE_KEYS["CL"]
        if tkey not in ALLOWED_TEMPLATE_KEYS_BY_TYPE["CL"]:
            return "", [f"Invalid CL template key '{tkey}'."]
        seq = LIGHT_CONSTANT_DOMAINS.get(tkey)
        if not isinstance(seq, str) or not seq:
            return "", [f"Missing template sequence for '{tkey}'."]
        return seq, []

    # Domains outside constant families that allow custom sequence overrides: VH, VL, G4S
    tkey = domain.template_key or DEFAULT_TEMPLATE_KEYS[dtype]
    if tkey not in ALLOWED_TEMPLATE_KEYS_BY_TYPE[dtype]:
        return "", [f"Invalid template key '{tkey}' for domain type '{dtype}'."]

    if domain.custom_sequence:
        if not isinstance(domain.custom_sequence, str):
            return "", ["Custom sequence must be a string."]
        if len(domain.custom_sequence) > MAX_CUSTOM_SEQ_LENGTH:
            return "", [
                f"Custom sequence exceeds maximum length of {MAX_CUSTOM_SEQ_LENGTH} residues."
            ]
        if not _VALID_SEQ_RE.fullmatch(domain.custom_sequence):
            return "", [
                "Custom sequence must contain only uppercase amino acid letters (A-Z)."
            ]
        return domain.custom_sequence, []

    table = VARIABLE_DOMAINS if dtype in ("VH", "VL") else DEFAULT_LINKERS
    seq = table.get(tkey)
    if not isinstance(seq, str) or not seq:
        return "", [f"Missing built-in template sequence for '{tkey}'."]
    return seq, []


def resolve_domain_sequence(
    domain: DomainBlock, preset_map: Optional[Dict[str, str]] = None
) -> Tuple[str, List[str]]:
    """Resolve the final engineered amino acid sequence for a DomainBlock."""
    base_seq, base_errs = resolve_domain_base_sequence(domain)
    muts, mut_errs = get_effective_domain_mutations(domain, preset_map)
    errors = base_errs + mut_errs
    if errors:
        return "", errors
    if not muts:
        return base_seq, []
    mut_seq, apply_errs = apply_domain_mutations(
        base_seq, "/".join(muts), domain.domain_type, domain.isotype
    )
    if apply_errs:
        return "", apply_errs
    return mut_seq, []


def assemble_chain_sequence(
    chain: PolypeptideChain, preset_map: Optional[Dict[str, str]] = None
) -> Tuple[str, List[Tuple[int, int, DomainBlock]], List[str]]:
    """Concatenate N-to-C domain sequences into an expressible chain sequence.

    Returns (full_sequence, 0-based [start, end) spans per domain, errors).
    """
    if not chain.domains:
        return "", [], [f"Chain '{chain.name or chain.id}' contains no domains."]
    if len(chain.domains) > MAX_DOMAINS_PER_CHAIN:
        return (
            "",
            [],
            [
                f"Chain '{chain.name or chain.id}' exceeds maximum of {MAX_DOMAINS_PER_CHAIN} domains."
            ],
        )

    parts: List[str] = []
    spans: List[Tuple[int, int, DomainBlock]] = []
    errors: List[str] = []
    cursor = 0

    for dom in chain.domains:
        seq, errs = resolve_domain_sequence(dom, preset_map)
        if errs:
            errors.extend(f"[{chain.name or chain.id}:{dom.id}] {e}" for e in errs)
            continue
        start = cursor
        cursor += len(seq)
        parts.append(seq)
        spans.append((start, cursor, dom))

    if errors:
        return "", [], errors

    full_seq = "".join(parts)
    if len(full_seq) > MAX_EXPORT_SEQUENCE_LENGTH:
        return (
            "",
            [],
            [
                f"Assembled chain '{chain.name or chain.id}' length ({len(full_seq)}) exceeds {MAX_EXPORT_SEQUENCE_LENGTH} residues."
            ],
        )
    return full_seq, spans, []


def _append_domain_history(domain: DomainBlock, entry: str) -> None:
    domain.mutation_history.append(entry[:MAX_METADATA_LENGTH])
    if len(domain.mutation_history) > MAX_HISTORY_ENTRIES:
        domain.mutation_history = domain.mutation_history[-MAX_HISTORY_ENTRIES:]


def route_mutations_to_chain(
    chain: PolypeptideChain,
    mutants_str: str,
    preset_label: str = "",
    preset_map: Optional[Dict[str, str]] = None,
) -> List[str]:
    """Atomically route constant-domain EU mutations (preset or ad-hoc) onto a chain."""
    resolved_map = preset_map if preset_map is not None else _builtin_preset_map()
    stripped_label = ""
    effective_mut_str = mutants_str

    if preset_label:
        stripped_label = preset_label.strip()
        if not stripped_label:
            return ["Preset label cannot be empty."]
        if stripped_label not in resolved_map:
            return [f"Unknown mutation preset '{stripped_label}'."]
        preset_tokens, preset_errs = parse_mutation_list(resolved_map[stripped_label])
        if preset_errs or not preset_tokens:
            return preset_errs or [f"Preset '{stripped_label}' contains no mutations."]
        if not isinstance(mutants_str, str) or not mutants_str.strip():
            effective_mut_str = resolved_map[stripped_label]
        else:
            parsed_input, input_errs = parse_mutation_list(mutants_str)
            if input_errs:
                return input_errs
            if parsed_input != preset_tokens:
                return [
                    f"Mutations '{mutants_str}' do not match preset '{stripped_label}'."
                ]

    parsed_tokens, parse_errs = parse_mutation_list(effective_mut_str)
    if parse_errs:
        return parse_errs
    if not parsed_tokens:
        return ["No mutations provided."]

    grouped: Dict[str, List[str]] = {}
    errors: List[str] = []

    for wt, pos, mut in parsed_tokens:
        target_type: Optional[str] = None
        for dt, (eu_start, eu_end) in CONSTANT_DOMAIN_EU_RANGES.items():
            if eu_start <= pos <= eu_end:
                target_type = dt
                break
        if target_type is None:
            errors.append(
                f"Position {pos} is outside constant domain EU range (118-447)."
            )
            continue
        grouped.setdefault(target_type, []).append(f"{wt}{pos}{mut}")

    if errors:
        return errors

    touched_pairs: List[Tuple[DomainBlock, List[str]]] = []
    for target_type, toks in grouped.items():
        matching = [d for d in chain.domains if d.domain_type == target_type]
        if not matching:
            errors.append(
                f"Chain '{chain.name or chain.id}' lacks a {target_type} domain for {'/'.join(toks)}."
            )
            continue
        for dom in matching:
            touched_pairs.append((dom, toks))

    if errors:
        return errors

    # Phase 1: Atomic dry-run across every touched domain before mutating any state
    for dom, toks in touched_pairs:
        base_seq, base_errs = resolve_domain_base_sequence(dom)
        if base_errs:
            errors.extend(base_errs)
            continue
        existing_muts, eff_errs = get_effective_domain_mutations(dom, resolved_map)
        if eff_errs:
            errors.extend(eff_errs)
            continue
        combined: List[str] = list(existing_muts)
        for tok in toks:
            if tok not in combined:
                combined.append(tok)
        _, apply_errs = apply_domain_mutations(
            base_seq, "/".join(combined), dom.domain_type, dom.isotype
        )
        if apply_errs:
            errors.extend(apply_errs)

    if errors:
        return errors

    # Phase 2: Atomic commit
    for dom, toks in touched_pairs:
        joined = "/".join(toks)
        if stripped_label:
            if stripped_label not in dom.applied_presets:
                dom.applied_presets.append(stripped_label)
            _append_domain_history(dom, f"Preset {stripped_label}: {joined}")
        else:
            for tok in toks:
                if tok not in dom.custom_mutations:
                    dom.custom_mutations.append(tok)
            _append_domain_history(dom, f"Custom: {joined}")

    return []


def remove_preset_from_chain(chain: PolypeptideChain, preset_label: str) -> None:
    """Remove preset_label from domains carrying it without affecting custom mutations or sibling presets."""
    for dom in chain.domains:
        if preset_label in dom.applied_presets:
            dom.applied_presets = [p for p in dom.applied_presets if p != preset_label]
            _append_domain_history(dom, f"Removed preset {preset_label}")


def remove_domain_from_chain(
    project: AntibodyProject, chain_id: str, domain_id: str
) -> None:
    """Remove a domain from a chain, prune orphaned pairings, and reflow unpinned coords."""
    target_chain: Optional[PolypeptideChain] = None
    for chain in project.chains:
        if chain.id == chain_id:
            target_chain = chain
            break
    if target_chain is None:
        return
    target_chain.domains = [d for d in target_chain.domains if d.id != domain_id]
    project.pairings = [
        p
        for p in project.pairings
        if not (
            (p.chain_a_id == chain_id and p.domain_a_id == domain_id)
            or (p.chain_b_id == chain_id and p.domain_b_id == domain_id)
        )
    ]
    reflow_chain_coordinates(target_chain)


def remove_chain_from_project(project: AntibodyProject, chain_id: str) -> None:
    """Remove a polypeptide chain and any InterChainPairing records referencing it."""
    project.chains = [c for c in project.chains if c.id != chain_id]
    project.pairings = [
        p
        for p in project.pairings
        if p.chain_a_id != chain_id and p.chain_b_id != chain_id
    ]


def diff_chain_sequences(
    chain: PolypeptideChain, preset_map: Optional[Dict[str, str]] = None
) -> Tuple[List[Tuple[str, int, str, str]], List[str]]:
    """Collect constant-domain EU diffs vs WT for every constant domain on a chain."""
    diffs: List[Tuple[str, int, str, str]] = []
    errors: List[str] = []

    for dom in chain.domains:
        if dom.domain_type not in CONSTANT_DOMAIN_EU_RANGES:
            continue
        wt_seq, base_errs = resolve_domain_base_sequence(dom)
        mut_seq, mut_errs = resolve_domain_sequence(dom, preset_map)
        if base_errs or mut_errs:
            errors.extend(base_errs + mut_errs)
            continue
        for pos, wt_aa, mut_aa in diff_sequences(
            wt_seq, mut_seq, dom.isotype, domain_type=dom.domain_type
        ):
            diffs.append((dom.domain_type, pos, wt_aa, mut_aa))

    return diffs, errors


def analyze_chain_developability(
    chain: PolypeptideChain, preset_map: Optional[Dict[str, str]] = None
) -> Tuple[List[Tuple[str, List[int]]], List[str]]:
    """Aggregate constant-domain PTM liabilities across a chain's constant domains."""
    by_category: Dict[str, set[int]] = {cat: set() for cat in _DEVELOPABILITY_CATEGORIES}
    errors: List[str] = []

    for dom in chain.domains:
        if dom.domain_type not in CONSTANT_DOMAIN_EU_RANGES:
            continue
        seq, errs = resolve_domain_sequence(dom, preset_map)
        if errs:
            errors.extend(errs)
            continue
        for cat, positions in analyze_developability(
            seq, dom.isotype, domain_type=dom.domain_type
        ):
            if cat in by_category:
                by_category[cat].update(positions)

    report = [(cat, sorted(by_category[cat])) for cat in _DEVELOPABILITY_CATEGORIES]
    return report, errors


def add_file_custom_mutation(
    project: AntibodyProject, label: str, mutants_str: str
) -> List[str]:
    """Validate and register a custom Fc mutation preset scoped to this AntibodyProject."""
    if not isinstance(label, str):
        return ["Label must be a string."]
    stripped = label.strip()
    if not stripped:
        return ["Label cannot be empty."]
    if len(stripped) > 40:
        return ["Label exceeds maximum length of 40 characters."]
    if not all(0x20 <= ord(ch) <= 0x7E for ch in stripped):
        return ["Label must contain only printable ASCII characters."]

    builtin_labels = set(_builtin_preset_map().keys())
    if stripped in builtin_labels:
        return [f"Label '{stripped}' shadows a built-in preset."]
    if any(
        isinstance(item, dict) and item.get("label") == stripped
        for item in project.file_custom_mutations
    ):
        return [f"Label '{stripped}' already exists in this project."]
    if len(project.file_custom_mutations) >= MAX_FILE_CUSTOM_MUTATIONS:
        return [
            f"Maximum of {MAX_FILE_CUSTOM_MUTATIONS} file custom mutations allowed."
        ]

    parsed, parse_errs = parse_mutation_list(mutants_str)
    if parse_errs:
        return parse_errs
    if not parsed:
        return ["Mutation value cannot be empty."]

    normalized = "/".join(f"{wt}{pos}{mut}" for wt, pos, mut in parsed)
    project.file_custom_mutations.append({"label": stripped, "value": normalized})
    return []


def _resolve_domain_style(
    domain_type: str, theme: Optional[SchematicTheme] = None
) -> DomainVisualStyle:
    active = theme if theme is not None else DEFAULT_SCHEMATIC_THEME
    if domain_type in active.styles:
        return active.styles[domain_type]
    if domain_type in DEFAULT_SCHEMATIC_THEME.styles:
        return DEFAULT_SCHEMATIC_THEME.styles[domain_type]
    fallback = DOMAIN_STYLE_MAP.get(domain_type)
    if fallback is None:
        raise ValueError(f"Unsupported domain type '{domain_type}'.")
    w, h, col = fallback
    return DomainVisualStyle(col, "#0F172A", 1.5, "#FFFFFF", w, h, "rounded_rect")


def reflow_chain_coordinates(
    chain: PolypeptideChain, theme: Optional[SchematicTheme] = None
) -> None:
    """Position unpinned DomainBlocks vertically (N->C) or horizontally from the chain anchor."""
    active = theme if theme is not None else DEFAULT_SCHEMATIC_THEME
    if active.orientation == "vertical":
        col_width = max(
            (_resolve_domain_style(d.domain_type, active).width for d in chain.domains),
            default=52.0,
        )
        cursor = float(chain.y)
        for dom in chain.domains:
            style = _resolve_domain_style(dom.domain_type, active)
            if not dom.pinned:
                nx = chain.x + (col_width - style.width) / 2.0
                dom.x = max(-MAX_COORD, min(MAX_COORD, float(nx)))
                dom.y = max(-MAX_COORD, min(MAX_COORD, float(cursor)))
            cursor += style.height + active.domain_spacing
    else:
        row_height = max(
            (_resolve_domain_style(d.domain_type, active).height for d in chain.domains),
            default=44.0,
        )
        cursor = float(chain.x)
        for dom in chain.domains:
            style = _resolve_domain_style(dom.domain_type, active)
            if not dom.pinned:
                ny = chain.y + (row_height - style.height) / 2.0
                dom.x = max(-MAX_COORD, min(MAX_COORD, float(cursor)))
                dom.y = max(-MAX_COORD, min(MAX_COORD, float(ny)))
            cursor += style.width + active.domain_spacing


def reflow_project_coordinates(
    project: AntibodyProject, theme: Optional[SchematicTheme] = None
) -> None:
    """Reflow all unpinned chain anchors and unpinned domains according to theme orientation and spacing."""
    active = theme if theme is not None else DEFAULT_SCHEMATIC_THEME
    for idx, chain in enumerate(project.chains):
        if active.orientation == "vertical":
            chain.x = 40.0 + idx * active.chain_spacing
            chain.y = 40.0
        else:
            chain.x = 40.0
            chain.y = 40.0 + idx * active.chain_spacing
        reflow_chain_coordinates(chain, theme=active)


def _make_domain(
    dom_id: str, domain_type: str, isotype: str, allotype: str, template_key: str = ""
) -> DomainBlock:
    is_const = domain_type in CONSTANT_DOMAIN_EU_RANGES
    return DomainBlock(
        id=dom_id,
        domain_type=domain_type,
        isotype=isotype,
        allotype=allotype if is_const else "",
        template_key="" if is_const else template_key,
    )


def _route_preset_checked(chain: PolypeptideChain, label: str) -> None:
    errs = route_mutations_to_chain(chain, "", preset_label=label)
    if errs:
        raise ValueError(
            f"Template preset '{label}' failed on chain '{chain.id}': {'; '.join(errs)}"
        )


def build_format_template(
    format_type: str, isotype: str = "igg1", allotype: Optional[str] = None
) -> AntibodyProject:
    """Construct a pre-wired AntibodyProject starter workspace for a supported format."""
    if format_type not in VALID_FORMAT_TYPES:
        raise ValueError(f"Unsupported format template '{format_type}'.")
    iso_map = SEQUENCES.get(isotype)
    if not isinstance(iso_map, dict) or not iso_map:
        raise ValueError(f"Unknown or empty isotype '{isotype}'.")
    if allotype is not None and allotype != "" and allotype not in iso_map:
        raise ValueError(f"Unknown allotype '{allotype}' for isotype '{isotype}'.")
    resolved_allotype = allotype if allotype else next(iter(iso_map))

    def mk(did: str, dtype: str, tkey: str = "") -> DomainBlock:
        return _make_domain(did, dtype, isotype, resolved_allotype, tkey)

    chains: List[PolypeptideChain] = []
    pairings: List[InterChainPairing] = []

    if format_type == "igg":
        title = "Monospecific IgG"
        chains = [
            PolypeptideChain(
                id="LC_1",
                name="LC_1",
                domains=[mk("lc1_vl", "VL", "VL_A"), mk("lc1_cl", "CL", "CL_kappa")],
            ),
            PolypeptideChain(
                id="HC_1",
                name="HC_1",
                domains=[
                    mk("hc1_vh", "VH", "VH_A"),
                    mk("hc1_ch1", "CH1"),
                    mk("hc1_hinge", "Hinge"),
                    mk("hc1_ch2", "CH2"),
                    mk("hc1_ch3", "CH3"),
                ],
            ),
            PolypeptideChain(
                id="HC_2",
                name="HC_2",
                domains=[
                    mk("hc2_vh", "VH", "VH_A"),
                    mk("hc2_ch1", "CH1"),
                    mk("hc2_hinge", "Hinge"),
                    mk("hc2_ch2", "CH2"),
                    mk("hc2_ch3", "CH3"),
                ],
            ),
            PolypeptideChain(
                id="LC_2",
                name="LC_2",
                domains=[mk("lc2_vl", "VL", "VL_A"), mk("lc2_cl", "CL", "CL_kappa")],
            ),
        ]
        pairings = [
            InterChainPairing("pair_1", "HC_1", "hc1_vh", "LC_1", "lc1_vl", "VH_VL"),
            InterChainPairing("pair_2", "HC_1", "hc1_ch1", "LC_1", "lc1_cl", "CH1_CL"),
            InterChainPairing("pair_3", "HC_2", "hc2_vh", "LC_2", "lc2_vl", "VH_VL"),
            InterChainPairing("pair_4", "HC_2", "hc2_ch1", "LC_2", "lc2_cl", "CH1_CL"),
            InterChainPairing(
                "pair_5", "HC_1", "hc1_hinge", "HC_2", "hc2_hinge", "HINGE_DISULFIDE"
            ),
            InterChainPairing("pair_6", "HC_1", "hc1_ch3", "HC_2", "hc2_ch3", "CH3_DIMER"),
        ]

    elif format_type == "crossmab":
        title = "CrossMab"
        hc_knob = PolypeptideChain(
            id="HC_Knob",
            name="HC_Knob",
            domains=[
                mk("hck_vh", "VH", "VH_A"),
                mk("hck_ch1", "CH1"),
                mk("hck_hinge", "Hinge"),
                mk("hck_ch2", "CH2"),
                mk("hck_ch3", "CH3"),
            ],
        )
        hc_hole = PolypeptideChain(
            id="HC_Hole_Cross",
            name="HC_Hole_Cross",
            domains=[
                mk("hch_vh", "VH", "VH_B"),
                mk("hch_cl", "CL", "CL_kappa"),
                mk("hch_hinge", "Hinge"),
                mk("hch_ch2", "CH2"),
                mk("hch_ch3", "CH3"),
            ],
        )
        _route_preset_checked(hc_knob, "Knob")
        _route_preset_checked(hc_hole, "Hole")

        chains = [
            PolypeptideChain(
                id="LC_A",
                name="LC_A",
                domains=[mk("lca_vl", "VL", "VL_A"), mk("lca_cl", "CL", "CL_kappa")],
            ),
            hc_knob,
            hc_hole,
            PolypeptideChain(
                id="LC_B_Cross",
                name="LC_B_Cross",
                domains=[mk("lcb_vl", "VL", "VL_B"), mk("lcb_ch1", "CH1")],
            ),
        ]
        pairings = [
            InterChainPairing("pair_1", "HC_Knob", "hck_vh", "LC_A", "lca_vl", "VH_VL"),
            InterChainPairing("pair_2", "HC_Knob", "hck_ch1", "LC_A", "lca_cl", "CH1_CL"),
            InterChainPairing(
                "pair_3", "HC_Hole_Cross", "hch_vh", "LC_B_Cross", "lcb_vl", "VH_VL"
            ),
            InterChainPairing(
                "pair_4", "HC_Hole_Cross", "hch_cl", "LC_B_Cross", "lcb_ch1", "CH1_CL"
            ),
            InterChainPairing(
                "pair_5",
                "HC_Knob",
                "hck_hinge",
                "HC_Hole_Cross",
                "hch_hinge",
                "HINGE_DISULFIDE",
            ),
            InterChainPairing(
                "pair_6", "HC_Knob", "hck_ch3", "HC_Hole_Cross", "hch_ch3", "CH3_DIMER"
            ),
        ]

    elif format_type == "igg_scfv":
        title = "IgG-scFv"
        chains = [
            PolypeptideChain(
                id="LC_1",
                name="LC_1",
                domains=[mk("lc1_vl", "VL", "VL_A"), mk("lc1_cl", "CL", "CL_kappa")],
            ),
            PolypeptideChain(
                id="HC_1",
                name="HC_1",
                domains=[
                    mk("hc1_vha", "VH", "VH_A"),
                    mk("hc1_ch1", "CH1"),
                    mk("hc1_hinge", "Hinge"),
                    mk("hc1_ch2", "CH2"),
                    mk("hc1_ch3", "CH3"),
                    mk("hc1_link1", "G4S", "G4S_3"),
                    mk("hc1_vhb", "VH", "VH_B"),
                    mk("hc1_link2", "G4S", "G4S_3"),
                    mk("hc1_vlb", "VL", "VL_B"),
                ],
            ),
            PolypeptideChain(
                id="HC_2",
                name="HC_2",
                domains=[
                    mk("hc2_vha", "VH", "VH_A"),
                    mk("hc2_ch1", "CH1"),
                    mk("hc2_hinge", "Hinge"),
                    mk("hc2_ch2", "CH2"),
                    mk("hc2_ch3", "CH3"),
                    mk("hc2_link1", "G4S", "G4S_3"),
                    mk("hc2_vhb", "VH", "VH_B"),
                    mk("hc2_link2", "G4S", "G4S_3"),
                    mk("hc2_vlb", "VL", "VL_B"),
                ],
            ),
            PolypeptideChain(
                id="LC_2",
                name="LC_2",
                domains=[mk("lc2_vl", "VL", "VL_A"), mk("lc2_cl", "CL", "CL_kappa")],
            ),
        ]
        pairings = [
            InterChainPairing("pair_1", "HC_1", "hc1_vha", "LC_1", "lc1_vl", "VH_VL"),
            InterChainPairing("pair_2", "HC_1", "hc1_ch1", "LC_1", "lc1_cl", "CH1_CL"),
            InterChainPairing("pair_3", "HC_2", "hc2_vha", "LC_2", "lc2_vl", "VH_VL"),
            InterChainPairing("pair_4", "HC_2", "hc2_ch1", "LC_2", "lc2_cl", "CH1_CL"),
            InterChainPairing(
                "pair_5", "HC_1", "hc1_hinge", "HC_2", "hc2_hinge", "HINGE_DISULFIDE"
            ),
            InterChainPairing("pair_6", "HC_1", "hc1_ch3", "HC_2", "hc2_ch3", "CH3_DIMER"),
            InterChainPairing("pair_7", "HC_1", "hc1_vhb", "HC_1", "hc1_vlb", "VH_VL"),
            InterChainPairing("pair_8", "HC_2", "hc2_vhb", "HC_2", "hc2_vlb", "VH_VL"),
        ]

    elif format_type == "scfv_fc":
        title = "scFv-Fc"
        chain_a = PolypeptideChain(
            id="Chain_A",
            name="Chain_A",
            domains=[
                mk("ca_vh", "VH", "VH_A"),
                mk("ca_link", "G4S", "G4S_3"),
                mk("ca_vl", "VL", "VL_A"),
                mk("ca_hinge", "Hinge"),
                mk("ca_ch2", "CH2"),
                mk("ca_ch3", "CH3"),
            ],
        )
        chain_b = PolypeptideChain(
            id="Chain_B",
            name="Chain_B",
            domains=[
                mk("cb_vh", "VH", "VH_B"),
                mk("cb_link", "G4S", "G4S_3"),
                mk("cb_vl", "VL", "VL_B"),
                mk("cb_hinge", "Hinge"),
                mk("cb_ch2", "CH2"),
                mk("cb_ch3", "CH3"),
            ],
        )
        _route_preset_checked(chain_a, "Knob")
        _route_preset_checked(chain_b, "Hole")
        chains = [chain_a, chain_b]
        pairings = [
            InterChainPairing("pair_1", "Chain_A", "ca_vh", "Chain_A", "ca_vl", "VH_VL"),
            InterChainPairing("pair_2", "Chain_B", "cb_vh", "Chain_B", "cb_vl", "VH_VL"),
            InterChainPairing(
                "pair_3",
                "Chain_A",
                "ca_hinge",
                "Chain_B",
                "cb_hinge",
                "HINGE_DISULFIDE",
            ),
            InterChainPairing(
                "pair_4", "Chain_A", "ca_ch3", "Chain_B", "cb_ch3", "CH3_DIMER"
            ),
        ]

    elif format_type == "dvd_ig":
        title = "DVD-Ig"
        chains = [
            PolypeptideChain(
                id="LC_1",
                name="LC_1",
                domains=[
                    mk("lc1_vlb", "VL", "VL_B"),
                    mk("lc1_link", "G4S", "G4S_2"),
                    mk("lc1_vla", "VL", "VL_A"),
                    mk("lc1_cl", "CL", "CL_kappa"),
                ],
            ),
            PolypeptideChain(
                id="HC_1",
                name="HC_1",
                domains=[
                    mk("hc1_vhb", "VH", "VH_B"),
                    mk("hc1_link", "G4S", "G4S_2"),
                    mk("hc1_vha", "VH", "VH_A"),
                    mk("hc1_ch1", "CH1"),
                    mk("hc1_hinge", "Hinge"),
                    mk("hc1_ch2", "CH2"),
                    mk("hc1_ch3", "CH3"),
                ],
            ),
            PolypeptideChain(
                id="HC_2",
                name="HC_2",
                domains=[
                    mk("hc2_vhb", "VH", "VH_B"),
                    mk("hc2_link", "G4S", "G4S_2"),
                    mk("hc2_vha", "VH", "VH_A"),
                    mk("hc2_ch1", "CH1"),
                    mk("hc2_hinge", "Hinge"),
                    mk("hc2_ch2", "CH2"),
                    mk("hc2_ch3", "CH3"),
                ],
            ),
            PolypeptideChain(
                id="LC_2",
                name="LC_2",
                domains=[
                    mk("lc2_vlb", "VL", "VL_B"),
                    mk("lc2_link", "G4S", "G4S_2"),
                    mk("lc2_vla", "VL", "VL_A"),
                    mk("lc2_cl", "CL", "CL_kappa"),
                ],
            ),
        ]
        pairings = [
            InterChainPairing("pair_1", "HC_1", "hc1_vhb", "LC_1", "lc1_vlb", "VH_VL"),
            InterChainPairing("pair_2", "HC_1", "hc1_vha", "LC_1", "lc1_vla", "VH_VL"),
            InterChainPairing("pair_3", "HC_1", "hc1_ch1", "LC_1", "lc1_cl", "CH1_CL"),
            InterChainPairing("pair_4", "HC_2", "hc2_vhb", "LC_2", "lc2_vlb", "VH_VL"),
            InterChainPairing("pair_5", "HC_2", "hc2_vha", "LC_2", "lc2_vla", "VH_VL"),
            InterChainPairing("pair_6", "HC_2", "hc2_ch1", "LC_2", "lc2_cl", "CH1_CL"),
            InterChainPairing(
                "pair_7", "HC_1", "hc1_hinge", "HC_2", "hc2_hinge", "HINGE_DISULFIDE"
            ),
            InterChainPairing("pair_8", "HC_1", "hc1_ch3", "HC_2", "hc2_ch3", "CH3_DIMER"),
        ]

    elif format_type == "bite":
        title = "BiTE"
        chains = [
            PolypeptideChain(
                id="Chain_1",
                name="Chain_1",
                domains=[
                    mk("c1_vha", "VH", "VH_A"),
                    mk("c1_l1", "G4S", "G4S_3"),
                    mk("c1_vla", "VL", "VL_A"),
                    mk("c1_l2", "G4S", "G4S_1"),
                    mk("c1_vhb", "VH", "VH_B"),
                    mk("c1_l3", "G4S", "G4S_3"),
                    mk("c1_vlb", "VL", "VL_B"),
                ],
            )
        ]
        pairings = [
            InterChainPairing("pair_1", "Chain_1", "c1_vha", "Chain_1", "c1_vla", "VH_VL"),
            InterChainPairing("pair_2", "Chain_1", "c1_vhb", "Chain_1", "c1_vlb", "VH_VL"),
        ]

    else:  # "custom"
        title = "Custom Construct"
        chains = [
            PolypeptideChain(
                id="Chain_1",
                name="Chain_1",
                domains=[mk("c1_vh", "VH", "VH_A")],
            )
        ]

    for idx, chain in enumerate(chains):
        chain.x = 40.0 + idx * DEFAULT_SCHEMATIC_THEME.chain_spacing
        chain.y = 40.0
        reflow_chain_coordinates(chain, theme=DEFAULT_SCHEMATIC_THEME)

    return AntibodyProject(
        name=title,
        format_type=format_type,
        chains=chains,
        pairings=pairings,
        file_custom_mutations=[],
    )


def _validate_coord(value: float, axis_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"Coordinate {axis_name} must be a number.")
    fval = float(value)
    if not math.isfinite(fval) or fval < -MAX_COORD or fval > MAX_COORD:
        raise ValueError(
            f"Coordinate {axis_name}={value} must be finite and within [{-MAX_COORD}, {MAX_COORD}]."
        )
    return fval


def _domain_display_label(domain: DomainBlock) -> str:
    dtype = domain.domain_type
    if dtype in ("VH", "VL"):
        if domain.custom_sequence:
            return f"{dtype}*"
        return domain.template_key or dtype
    if dtype == "CL":
        tkey = domain.template_key or "CL_kappa"
        return "CLκ" if tkey == "CL_kappa" else "CLλ"
    if dtype == "G4S":
        if domain.custom_sequence:
            return "G4S*"
        return domain.template_key or "G4S"
    return dtype


def build_schematic_scene(
    project: AntibodyProject,
    preset_map: Optional[Dict[str, str]] = None,
    theme: Optional[SchematicTheme] = None,
) -> SchematicScene:
    """Build the shared 2D scene representation for interactive Flet rendering and PNG export."""
    resolved_map = preset_map if preset_map is not None else build_preset_map(project)
    active_theme = theme if theme is not None else DEFAULT_SCHEMATIC_THEME
    nodes: List[SchematicNode] = []
    edges: List[SchematicEdge] = []
    pairing_edges: List[SchematicPairingEdge] = []
    node_lookup: Dict[Tuple[str, str], SchematicNode] = {}

    for chain in project.chains:
        chain_nodes: List[SchematicNode] = []
        for idx, dom in enumerate(chain.domains):
            nx = _validate_coord(dom.x, f"{dom.id}.x")
            ny = _validate_coord(dom.y, f"{dom.id}.y")
            style = _resolve_domain_style(dom.domain_type, active_theme)

            _, errs = resolve_domain_sequence(dom, resolved_map)
            badge_parts = list(dom.applied_presets)
            custom_preview = list(dom.custom_mutations[:2])
            if len(dom.custom_mutations) > 2:
                custom_preview.append(f"+{len(dom.custom_mutations) - 2}")
            badge_parts.extend(custom_preview)

            node = SchematicNode(
                domain_id=dom.id,
                chain_id=chain.id,
                chain_name=chain.name,
                is_n_terminus=(idx == 0),
                label=_domain_display_label(dom),
                domain_type=dom.domain_type,
                x=nx,
                y=ny,
                width=style.width,
                height=style.height,
                fill_color=style.fill_color,
                badge_text=", ".join(badge_parts),
                error_flag=bool(errs),
                border_color=style.border_color,
                border_width=style.border_width,
                text_color=style.text_color,
                shape=style.shape,
            )
            nodes.append(node)
            chain_nodes.append(node)
            node_lookup[(chain.id, dom.id)] = node

        for prev_node, next_node in zip(chain_nodes, chain_nodes[1:]):
            if active_theme.orientation == "vertical":
                x1 = prev_node.x + prev_node.width / 2.0
                y1 = prev_node.y + prev_node.height
                x2 = next_node.x + next_node.width / 2.0
                y2 = next_node.y
            else:
                x1 = prev_node.x + prev_node.width
                y1 = prev_node.y + prev_node.height / 2.0
                x2 = next_node.x
                y2 = next_node.y + next_node.height / 2.0
            edges.append(
                SchematicEdge(
                    chain_id=chain.id,
                    from_domain_id=prev_node.domain_id,
                    to_domain_id=next_node.domain_id,
                    x1=x1,
                    y1=y1,
                    x2=x2,
                    y2=y2,
                    edge_kind="backbone",
                )
            )

    for pairing in project.pairings:
        node_a = node_lookup.get((pairing.chain_a_id, pairing.domain_a_id))
        node_b = node_lookup.get((pairing.chain_b_id, pairing.domain_b_id))
        if node_a is None or node_b is None:
            continue
        pairing_edges.append(
            SchematicPairingEdge(
                pairing_id=pairing.id,
                chain_a_id=pairing.chain_a_id,
                domain_a_id=pairing.domain_a_id,
                chain_b_id=pairing.chain_b_id,
                domain_b_id=pairing.domain_b_id,
                pairing_type=pairing.pairing_type,
                x1=node_a.x + node_a.width / 2.0,
                y1=node_a.y + node_a.height / 2.0,
                x2=node_b.x + node_b.width / 2.0,
                y2=node_b.y + node_b.height / 2.0,
            )
        )

    if not nodes:
        bounds = (0.0, 0.0, 320.0, 200.0)
    else:
        min_x = min(n.x for n in nodes)
        min_y = min(n.y for n in nodes)
        max_x = max(n.x + n.width for n in nodes)
        max_y = max(n.y + n.height for n in nodes)
        bounds = (
            min_x - SCENE_PADDING,
            min_y - SCENE_PADDING,
            max_x + SCENE_PADDING,
            max_y + SCENE_PADDING,
        )

    return SchematicScene(
        nodes=nodes,
        edges=edges,
        pairing_edges=pairing_edges,
        bounds=bounds,
    )
