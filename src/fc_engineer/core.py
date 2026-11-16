import re
from typing import Dict, List, Optional, Tuple

# --- Core Logic ---

EU_START = 118

# IgG3 carries an extended hinge with tandem repeats relative to IgG1 (+47 residues
# for the bundled WT(P01860-1) entry). EU positions 223-230 fall inside the repeated
# region where unique numbering is undefined, while CH2-CH3 (EU 231-447) aligns to
# the conserved C-terminal segment.
IGG3_HINGE_EXTRA = 47
IGG3_GAP_END = 230
_IGG3_CH2_START_IDX = IGG3_GAP_END + 1 - EU_START + IGG3_HINGE_EXTRA

CONSTANT_DOMAIN_EU_RANGES: Dict[str, Tuple[int, int]] = {
    "CH1": (118, 215),
    "Hinge": (216, 230),
    "CH2": (231, 340),
    "CH3": (341, 447),
}

MAX_MUTATION_STRING_LENGTH = 1000
MAX_MUTATIONS = 50
MAX_BATCH_ITEMS = 20
MAX_BATCH_LENGTH = 2000

_VALID_SEQ_RE = re.compile(r"[A-Z]+\Z")
_VALID_MUT_LIST_CHARS_RE = re.compile(r"[^a-zA-Z0-9/, ]")
# SECURITY: Use [0-9] instead of \d and \Z instead of $ to strictly enforce ASCII digits without trailing newlines
_MUTATION_TOKEN_RE = re.compile(r"[A-Z][0-9]+[A-Z]\Z")


def get_residue_index(pos: int, isotype: str) -> Optional[int]:
    """Map an EU position to a 0-based index within a full CH1-Hinge-CH2-CH3 sequence."""
    if isotype == "igg1":
        return pos - EU_START
    elif isotype in ["igg2", "igg4"]:
        if pos <= 222:
            return pos - EU_START
        elif 223 <= pos <= 225:
            return None
        else:
            return pos - EU_START - 3
    elif isotype == "igg3":
        if pos <= 222:
            return pos - EU_START
        elif pos <= IGG3_GAP_END:
            return None
        else:
            return pos - EU_START + IGG3_HINGE_EXTRA
    return None


def get_eu_position(index: int, isotype: str = "igg1") -> Optional[int]:
    """Inverse of get_residue_index: map a 0-based sequence index back to an EU position."""
    if index < 0:
        return None
    if isotype == "igg1":
        return index + EU_START
    elif isotype in ["igg2", "igg4"]:
        if index <= 104:
            return index + EU_START
        return index + EU_START + 3
    elif isotype == "igg3":
        if index <= 104:
            return index + EU_START
        if index < _IGG3_CH2_START_IDX:
            return None
        return index + EU_START - IGG3_HINGE_EXTRA
    return None


def get_domain_slice_bounds(
    domain_type: str, isotype: str, full_seq_len: int
) -> Optional[Tuple[int, int]]:
    """Compute [start, end) slice bounds for a constant domain in a full isotype sequence."""
    if domain_type not in CONSTANT_DOMAIN_EU_RANGES:
        return None
    if not isinstance(full_seq_len, int) or isinstance(full_seq_len, bool) or full_seq_len <= 0:
        return None

    eu_start, eu_end = CONSTANT_DOMAIN_EU_RANGES[domain_type]
    start_idx = get_residue_index(eu_start, isotype)
    if start_idx is None:
        return None

    if domain_type == "CH3":
        raw_end = full_seq_len
    else:
        next_idx = get_residue_index(eu_end + 1, isotype)
        raw_end = full_seq_len if next_idx is None else next_idx

    end_idx = min(raw_end, full_seq_len)
    if start_idx < 0 or start_idx >= end_idx:
        return None
    return start_idx, end_idx


def slice_constant_domain(full_seq: str, domain_type: str, isotype: str) -> str:
    """Slice an isolated constant domain (CH1, Hinge, CH2, CH3) from a full constant sequence."""
    if not isinstance(full_seq, str) or not full_seq:
        return ""
    bounds = get_domain_slice_bounds(domain_type, isotype, len(full_seq))
    if bounds is None:
        return ""
    return full_seq[bounds[0] : bounds[1]]


def get_domain_residue_index(
    pos: int, domain_type: str, isotype: str
) -> Optional[int]:
    """Map an EU position to a 0-based index within an isolated constant domain."""
    if domain_type not in CONSTANT_DOMAIN_EU_RANGES:
        return None
    eu_start, eu_end = CONSTANT_DOMAIN_EU_RANGES[domain_type]
    if not (eu_start <= pos <= eu_end):
        return None
    start_idx = get_residue_index(eu_start, isotype)
    global_idx = get_residue_index(pos, isotype)
    if start_idx is None or global_idx is None:
        return None
    return global_idx - start_idx


def get_domain_eu_position(
    local_idx: int, domain_type: str, isotype: str
) -> Optional[int]:
    """Map a 0-based index inside an isolated constant domain back to its EU position."""
    if domain_type not in CONSTANT_DOMAIN_EU_RANGES or local_idx < 0:
        return None
    eu_start, eu_end = CONSTANT_DOMAIN_EU_RANGES[domain_type]
    start_idx = get_residue_index(eu_start, isotype)
    if start_idx is None:
        return None
    pos = get_eu_position(start_idx + local_idx, isotype)
    if pos is not None and eu_start <= pos <= eu_end:
        return pos
    return None


def _lookup_eu_position(
    idx: int, isotype: str, domain_type: Optional[str] = None
) -> Optional[int]:
    if domain_type is not None:
        return get_domain_eu_position(idx, domain_type, isotype)
    return get_eu_position(idx, isotype)


def diff_sequences(
    wt_seq: str, mut_seq: str, isotype: str, domain_type: Optional[str] = None
) -> List[Tuple[int, str, str]]:
    """Compare WT vs mutant sequence and return (EU position, WT aa, mutant aa) tuples."""
    if domain_type is not None and domain_type not in CONSTANT_DOMAIN_EU_RANGES:
        return []
    diffs: List[Tuple[int, str, str]] = []
    for idx, (wt_aa, mut_aa) in enumerate(zip(wt_seq, mut_seq)):
        if wt_aa != mut_aa:
            pos = _lookup_eu_position(idx, isotype, domain_type)
            if pos is not None:
                diffs.append((pos, wt_aa, mut_aa))
    return diffs


def parse_mutation(m_str: str) -> Tuple[str, int, str]:
    """Parse a single EU mutation token such as 'L234A' into (wt_aa, pos, mut_aa)."""
    if not isinstance(m_str, str):
        raise ValueError("Mutation token must be a string")
    if len(m_str) > 10:
        raise ValueError(f"Mutation string too long (max 10): {len(m_str)}")
    m_str = m_str.upper()
    if not _MUTATION_TOKEN_RE.fullmatch(m_str):
        raise ValueError(f"Invalid mutation format: {m_str}")
    wt_aa = m_str[0]
    pos = int(m_str[1:-1])
    mut_aa = m_str[-1]
    return wt_aa, pos, mut_aa


def parse_mutation_list(
    mutants_str: str,
) -> Tuple[List[Tuple[str, int, str]], List[str]]:
    """Validate and tokenize a delimited EU mutation string into parsed tuples and errors."""
    if not isinstance(mutants_str, str):
        return [], ["Error: Mutation string must be a string."]
    if not mutants_str or not mutants_str.strip():
        return [], []
    if len(mutants_str) > MAX_MUTATION_STRING_LENGTH:
        return [], [
            f"Error: Mutation string exceeds maximum length of {MAX_MUTATION_STRING_LENGTH} characters."
        ]
    if _VALID_MUT_LIST_CHARS_RE.search(mutants_str):
        return [], ["Error: Mutation string contains invalid characters."]

    mut_list = [m.strip() for m in mutants_str.replace(",", "/").split("/") if m.strip()]
    if len(mut_list) > MAX_MUTATIONS:
        return [], [f"Error: Maximum of {MAX_MUTATIONS} mutations allowed."]

    parsed: List[Tuple[str, int, str]] = []
    errors: List[str] = []
    for m in mut_list:
        try:
            parsed.append(parse_mutation(m))
        except ValueError:
            errors.append(f"Format error: Invalid mutation '{m}'.")
    return parsed, errors


def _apply_mutations_impl(
    sequence: str,
    mutants_str: str,
    isotype: str,
    domain_type: Optional[str] = None,
) -> Tuple[str, List[str]]:
    parsed_tokens, errors = parse_mutation_list(mutants_str)
    if not parsed_tokens and not errors:
        return sequence, []

    seq_list = list(sequence)
    applied_by_pos: Dict[int, Tuple[str, str]] = {}

    for wt_aa, pos, mut_aa in parsed_tokens:
        if domain_type is None:
            index = get_residue_index(pos, isotype)
            if index is None:
                errors.append(f"Position {pos} is a Gap in {isotype}.")
                continue
        else:
            eu_start, eu_end = CONSTANT_DOMAIN_EU_RANGES[domain_type]
            if eu_start <= pos <= eu_end and get_residue_index(pos, isotype) is None:
                errors.append(f"Position {pos} is a Gap in {isotype}.")
                continue
            index = get_domain_residue_index(pos, domain_type, isotype)
            if index is None:
                errors.append(f"Position {pos} is out of range for {domain_type}.")
                continue

        if index < 0 or index >= len(seq_list):
            errors.append(f"Position {pos} is out of range.")
            continue

        prior = applied_by_pos.get(pos)
        if prior is not None:
            if prior == (wt_aa, mut_aa):
                continue
            errors.append(
                f"Pos {pos}: Conflicting mutations '{prior[1]}' and '{mut_aa}'."
            )
            continue

        orig_wt = sequence[index]
        if orig_wt != wt_aa:
            errors.append(f"Pos {pos}: Expected '{wt_aa}', found '{orig_wt}'.")
            continue

        applied_by_pos[pos] = (wt_aa, mut_aa)
        seq_list[index] = mut_aa

    if errors:
        return sequence, errors
    return "".join(seq_list), []


def apply_mutations(
    sequence: str, mutants_str: str, isotype: str
) -> Tuple[str, List[str]]:
    """Apply EU mutations to a contiguous CH1-Hinge-CH2-CH3 constant sequence."""
    # SECURITY: Defense-in-depth validation to ensure base sequence contains only valid amino acid characters
    if not isinstance(sequence, str):
        return "", ["Error: Base sequence must be a string."]
    if sequence and not _VALID_SEQ_RE.fullmatch(sequence):
        return sequence, ["Error: Base sequence contains invalid characters."]
    if not mutants_str:
        return sequence, []
    return _apply_mutations_impl(sequence, mutants_str, isotype, domain_type=None)


def apply_domain_mutations(
    domain_seq: str, mutants_str: str, domain_type: str, isotype: str
) -> Tuple[str, List[str]]:
    """Apply EU mutations to an isolated constant domain sequence (CH1, Hinge, CH2, CH3)."""
    if not isinstance(domain_seq, str) or not domain_seq or not _VALID_SEQ_RE.fullmatch(domain_seq):
        return (
            domain_seq if isinstance(domain_seq, str) else "",
            ["Error: Base domain sequence contains invalid characters."],
        )
    if domain_type not in CONSTANT_DOMAIN_EU_RANGES:
        return domain_seq, [f"Error: Unsupported constant domain '{domain_type}'."]
    if not isinstance(mutants_str, str):
        return domain_seq, ["Error: Mutation string must be a string."]
    return _apply_mutations_impl(
        domain_seq, mutants_str, isotype, domain_type=domain_type
    )


# --- Batch Generation ---


def generate_batch(
    base_seq: str, batch_str: str, isotype: str
) -> List[Tuple[str, str, List[str]]]:
    """Generate multiple sequences from a comma-separated list of mutation sets.

    Each comma-separated entry is one combination (it may itself contain
    slash-joined mutations). Returns a list of (combo, sequence, errors) tuples.
    """
    if not isinstance(batch_str, str) or not batch_str.strip():
        return []

    # SECURITY: Defense-in-depth bounds on batch input size and character set
    if len(batch_str) > MAX_BATCH_LENGTH:
        return [
            (
                "",
                base_seq,
                [
                    f"Error: Batch input exceeds maximum length of {MAX_BATCH_LENGTH} characters."
                ],
            )
        ]
    if _VALID_MUT_LIST_CHARS_RE.search(batch_str):
        return [("", base_seq, ["Error: Batch input contains invalid characters."])]

    combos = [c.strip() for c in batch_str.split(",") if c.strip()]
    if not combos:
        return []
    if len(combos) > MAX_BATCH_ITEMS:
        return [
            (
                "",
                base_seq,
                [f"Error: Maximum of {MAX_BATCH_ITEMS} batch combinations allowed."],
            )
        ]

    results: List[Tuple[str, str, List[str]]] = []
    for combo in combos:
        seq_out, errors = apply_mutations(base_seq, combo, isotype)
        results.append((combo, seq_out, errors))
    return results


# --- Developability Analysis ---


def analyze_developability(
    sequence: str, isotype: str, domain_type: Optional[str] = None
) -> List[Tuple[str, List[int]]]:
    """Scan a sequence for developability-related liabilities (PTM hotspots).

    Returns a list of (category, EU positions) tuples. Positions that do not
    map to an EU number are skipped. Categories:
      - N-glycosylation consensus (N-X-S/T)
      - Deamidation-prone NG motifs
      - Isomerization-prone DG motifs
      - Oxidation-prone Met residues
    """
    findings: List[Tuple[str, List[int]]] = []
    if not isinstance(sequence, str) or not sequence:
        return findings
    if not isinstance(isotype, str) or not isotype:
        return findings
    if domain_type is not None and domain_type not in CONSTANT_DOMAIN_EU_RANGES:
        return findings

    n = len(sequence)

    glyco: List[int] = []
    for i in range(n - 2):
        if (
            sequence[i] == "N"
            and sequence[i + 1] != "P"
            and sequence[i + 2] in ("S", "T")
        ):
            pos = _lookup_eu_position(i, isotype, domain_type)
            if pos is not None:
                glyco.append(pos)
    findings.append(("N-glycosylation motif (N-X-S/T)", glyco))

    ng: List[int] = []
    dg: List[int] = []
    met: List[int] = []
    for i in range(n):
        pos = _lookup_eu_position(i, isotype, domain_type)
        if pos is None:
            continue
        pair = sequence[i : i + 2]
        if pair == "NG":
            ng.append(pos)
        if pair == "DG":
            dg.append(pos)
        if sequence[i] == "M":
            met.append(pos)
    findings.append(("Deamidation-prone NG", ng))
    findings.append(("Isomerization-prone DG", dg))
    findings.append(("Oxidation-prone Met (M252/M428 = FcRn sites)", met))
    return findings


_AA_WEIGHTS: Dict[str, float] = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
_PKA_POS: Dict[str, float] = {"K": 10.5, "R": 12.5, "H": 6.0}
_PKA_NEG: Dict[str, float] = {"D": 3.9, "E": 4.2, "C": 8.3, "Y": 10.1}


def compute_protein_metrics(sequence: str) -> Tuple[float, float]:
    """Compute approximate molecular weight (kDa) and theoretical pI for an amino acid sequence."""
    if not isinstance(sequence, str) or not sequence:
        return 0.0, 7.0
    mw_da = 18.01528 + sum(_AA_WEIGHTS.get(aa, 110.0) for aa in sequence)
    counts = {aa: sequence.count(aa) for aa in set(sequence)}

    def net_charge(ph: float) -> float:
        pos = 1.0 / (1.0 + 10.0 ** (ph - 8.0))  # N-terminus
        neg = 1.0 / (1.0 + 10.0 ** (3.1 - ph))  # C-terminus
        for aa, pka in _PKA_POS.items():
            if aa in counts:
                pos += counts[aa] / (1.0 + 10.0 ** (ph - pka))
        for aa, pka in _PKA_NEG.items():
            if aa in counts:
                neg += counts[aa] / (1.0 + 10.0 ** (pka - ph))
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2.0
        if net_charge(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return round(mw_da / 1000.0, 2), round((lo + hi) / 2.0, 2)


def lookup_wt_residue(sequence: str, eu_pos: int, isotype: str) -> Optional[str]:
    """Return the wild-type amino acid at EU position eu_pos in a full constant sequence."""
    if not isinstance(sequence, str) or not sequence:
        return None
    idx = get_residue_index(eu_pos, isotype)
    if idx is None or not (0 <= idx < len(sequence)):
        return None
    return sequence[idx]
