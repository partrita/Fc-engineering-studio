import re
from typing import List, Optional, Tuple

# --- Core Logic ---

EU_START = 118

# IgG3 carries an extended hinge with tandem repeats relative to IgG1 (+47 residues
# for the bundled WT(P01860-1) entry). EU positions 223-230 fall inside the repeated
# region where unique numbering is undefined, while CH2-CH3 (EU 231-447) aligns to
# the conserved C-terminal segment.
IGG3_HINGE_EXTRA = 47
IGG3_GAP_END = 230

def get_residue_index(pos: int, isotype: str) -> Optional[int]:
    if isotype == "igg1": return pos - EU_START
    elif isotype in ["igg2", "igg4"]:
        if pos <= 222: return pos - EU_START
        elif 223 <= pos <= 225: return None
        else: return pos - EU_START - 3
    elif isotype == "igg3":
        if pos <= 222: return pos - EU_START
        elif pos <= IGG3_GAP_END: return None
        else: return pos - EU_START + IGG3_HINGE_EXTRA
    return None

def get_eu_position(index: int, isotype: str) -> Optional[int]:
    """Inverse of get_residue_index: maps a 0-based sequence index back to the EU position."""
    if isotype == "igg1": return index + EU_START
    elif isotype in ["igg2", "igg4"]:
        if index <= 104: return index + EU_START
        else: return index + EU_START + 3
    elif isotype == "igg3":
        if index <= 104: return index + EU_START
        else: return index + EU_START - IGG3_HINGE_EXTRA
    return None

def diff_sequences(wt_seq: str, mut_seq: str, isotype: str) -> List[Tuple[int, str, str]]:
    """Compare WT vs mutant sequence and return (EU position, WT aa, mutant aa) tuples."""
    diffs: List[Tuple[int, str, str]] = []
    for idx, (wt_aa, mut_aa) in enumerate(zip(wt_seq, mut_seq)):
        if wt_aa != mut_aa:
            pos = get_eu_position(idx, isotype)
            if pos is not None:
                diffs.append((pos, wt_aa, mut_aa))
    return diffs

def parse_mutation(m_str: str) -> Tuple[str, int, str]:
    if len(m_str) > 10:
        raise ValueError(f"Mutation string too long (max 10): {len(m_str)}")
    m_str = m_str.upper()
    # SECURITY: Use [0-9] instead of \d to strictly enforce ASCII digits and prevent Unicode digit injection
    if not re.fullmatch(r"[A-Z][0-9]+[A-Z]", m_str):
        raise ValueError(f"Invalid mutation format: {m_str}")
    wt_aa = m_str[0]
    pos = int(m_str[1:-1])
    mut_aa = m_str[-1]
    return wt_aa, pos, mut_aa

def apply_mutations(sequence: str, mutants_str: str, isotype: str) -> Tuple[str, List[str]]:
    # SECURITY: Defense-in-depth validation to ensure base sequence contains only valid amino acid characters
    if sequence and not re.fullmatch(r"[A-Z]+", sequence):
        return sequence, ["Error: Base sequence contains invalid characters."]

    if not mutants_str: return sequence, []

    # SECURITY: Enforce strict length limits on input string to prevent DoS (memory/CPU exhaustion)
    if len(mutants_str) > 1000:
        return sequence, ["Error: Mutation string exceeds maximum length of 1000 characters."]

    # SECURITY: Defense-in-depth validation on the backend logic layer
    if re.search(r'[^a-zA-Z0-9/, ]', mutants_str):
        return sequence, ["Error: Mutation string contains invalid characters."]

    mut_list = [m.strip() for m in mutants_str.replace(',', '/').split('/') if m.strip()]

    MAX_MUTATIONS = 50
    if len(mut_list) > MAX_MUTATIONS:
        return sequence, [f"Error: Maximum of {MAX_MUTATIONS} mutations allowed."]

    seq_list = list(sequence)
    errors = []
    for m in mut_list:
        try:
            wt_aa, pos, mut_aa = parse_mutation(m)
            index = get_residue_index(pos, isotype)
            if index is None: errors.append(f"Position {pos} is a Gap in {isotype}."); continue
            if index < 0 or index >= len(seq_list): errors.append(f"Position {pos} is out of range."); continue
            if seq_list[index] != wt_aa: errors.append(f"Pos {pos}: Expected '{wt_aa}', found '{seq_list[index]}'."); continue
            seq_list[index] = mut_aa
        except ValueError: errors.append(f"Format error: Invalid mutation '{m}'.")
    return "".join(seq_list), errors
