"""Curated Fc engineering knowledge catalog, functional profiling, and EU residue browser helpers."""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from fc_engineer.core import (
    CONSTANT_DOMAIN_EU_RANGES,
    analyze_developability,
    compute_protein_metrics,
    get_domain_eu_position,
)
from fc_engineer.designer import (
    PolypeptideChain,
    assemble_chain_sequence,
    diff_chain_sequences,
    resolve_domain_base_sequence,
    resolve_domain_sequence,
)

AMINO_ACIDS: Tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")

FC_CATEGORIES: Tuple[str, ...] = (
    "All",
    "Silence Effector",
    "Enhance ADCC",
    "Enhance CDC",
    "Extend Half-Life",
    "Reduce Half-Life",
    "Heterodimer",
    "Stability / Hinge",
    "Glyco / PTM",
    "Project",
)


@dataclass(frozen=True)
class PresetMetadata:
    """Biological annotation and functional profile for an Fc mutation preset."""

    label: str
    category: str
    summary: str
    adcc: str = "neutral"  # "up" | "down" | "neutral"
    cdc: str = "neutral"  # "up" | "down" | "neutral"
    fcrn: str = "neutral"  # "up" | "down" | "neutral"
    heterodimer: bool = False
    recommended_isotypes: Tuple[str, ...] = ("igg1",)


# Curated biological metadata for built-in presets in src/fc_engineer/data/mutants.yaml
_CURATED_PRESET_METADATA: Dict[str, PresetMetadata] = {
    "S228P": PresetMetadata(
        "S228P",
        "Stability / Hinge",
        "Prevents IgG4 Fab-arm exchange by stabilizing the core hinge.",
        recommended_isotypes=("igg4",),
    ),
    "LALA": PresetMetadata(
        "LALA",
        "Silence Effector",
        "Strongly reduces FcγR and C1q binding (ADCC/CDC silencing).",
        adcc="down",
        cdc="down",
    ),
    "PGLALA": PresetMetadata(
        "PGLALA",
        "Silence Effector",
        "Complete effector silencing combining P329G with LALA (ADCC/CDC abolished).",
        adcc="down",
        cdc="down",
    ),
    "Aglycosylation": PresetMetadata(
        "Aglycosylation",
        "Glyco / PTM",
        "Eliminates N297 glycosylation; strongly reduces FcγR/C1q effector activity.",
        adcc="down",
        cdc="down",
        recommended_isotypes=("igg1", "igg2", "igg3", "igg4"),
    ),
    "D265A": PresetMetadata(
        "D265A",
        "Silence Effector",
        "Abolishes FcγR binding while retaining N297 glycosylation.",
        adcc="down",
        cdc="down",
    ),
    "K322A": PresetMetadata(
        "K322A",
        "Silence Effector",
        "Selectively reduces C1q binding and CDC activity.",
        cdc="down",
    ),
    "E318A": PresetMetadata(
        "E318A",
        "Silence Effector",
        "Disrupts C1q binding site to attenuate CDC.",
        cdc="down",
    ),
    "G236R": PresetMetadata(
        "G236R",
        "Silence Effector",
        "Lower-hinge steric block reducing FcγR engagement.",
        adcc="down",
    ),
    "L328A": PresetMetadata(
        "L328A",
        "Silence Effector",
        "Reduces FcγR binding in the CH2 FG loop.",
        adcc="down",
    ),
    "YTE": PresetMetadata(
        "YTE",
        "Extend Half-Life",
        "Increases FcRn binding at pH 6.0 for ~3–4x serum half-life extension.",
        fcrn="up",
        recommended_isotypes=("igg1", "igg2", "igg4"),
    ),
    "LS": PresetMetadata(
        "LS",
        "Extend Half-Life",
        "Enhances pH 6.0 FcRn affinity for extended half-life with minimal ADCC impact.",
        fcrn="up",
        recommended_isotypes=("igg1", "igg2", "igg4"),
    ),
    "QL": PresetMetadata(
        "QL",
        "Extend Half-Life",
        "T250Q/M428L variant extending serum half-life via FcRn.",
        fcrn="up",
        recommended_isotypes=("igg1", "igg2", "igg4"),
    ),
    "N434A": PresetMetadata(
        "N434A",
        "Extend Half-Life",
        "Single-point FcRn affinity enhancement for moderate half-life extension.",
        fcrn="up",
        recommended_isotypes=("igg1", "igg2", "igg4"),
    ),
    "FcRnKO": PresetMetadata(
        "FcRnKO",
        "Reduce Half-Life",
        "Abolishes FcRn recycling (H310A/H435A) for rapid serum clearance.",
        fcrn="down",
        recommended_isotypes=("igg1", "igg2", "igg4"),
    ),
    "Knob": PresetMetadata(
        "Knob",
        "Heterodimer",
        "T366W steric protrusion in CH3 for Knob-into-Hole bispecific assembly.",
        heterodimer=True,
        recommended_isotypes=("igg1", "igg4"),
    ),
    "Hole": PresetMetadata(
        "Hole",
        "Heterodimer",
        "T366S/L368A/Y407V CH3 cavity pairing with Knob heavy chain.",
        heterodimer=True,
        recommended_isotypes=("igg1", "igg4"),
    ),
    "S354C": PresetMetadata(
        "S354C",
        "Heterodimer",
        "Engineered CH3 cysteine forming inter-chain disulfide with Y349C.",
        heterodimer=True,
        recommended_isotypes=("igg1", "igg4"),
    ),
    "Y349C": PresetMetadata(
        "Y349C",
        "Heterodimer",
        "Engineered CH3 cysteine forming inter-chain disulfide with S354C.",
        heterodimer=True,
        recommended_isotypes=("igg1", "igg4"),
    ),
    "Knob_Cys": PresetMetadata(
        "Knob_Cys",
        "Heterodimer",
        "Knob (T366W) + S354C stabilizing disulfide for high-purity heterodimer.",
        heterodimer=True,
        recommended_isotypes=("igg1", "igg4"),
    ),
    "Hole_Cys": PresetMetadata(
        "Hole_Cys",
        "Heterodimer",
        "Hole (T366S/L368A/Y407V) + Y349C stabilizing disulfide for heterodimer.",
        heterodimer=True,
        recommended_isotypes=("igg1", "igg4"),
    ),
    "ChargePair_A": PresetMetadata(
        "ChargePair_A",
        "Heterodimer",
        "Electrostatic steering (K392D/K409D) pairing with ChargePair_B.",
        heterodimer=True,
    ),
    "ChargePair_B": PresetMetadata(
        "ChargePair_B",
        "Heterodimer",
        "Electrostatic steering (D356K/D399K) pairing with ChargePair_A.",
        heterodimer=True,
    ),
    "DE": PresetMetadata(
        "DE",
        "Enhance ADCC",
        "S239D/I332E strongly enhances FcγRIIIa binding and ADCC.",
        adcc="up",
    ),
    "DLE": PresetMetadata(
        "DLE",
        "Enhance ADCC",
        "S239D/A330L/I332E maximizes ADCC potency via FcγRIIIa.",
        adcc="up",
    ),
    "GAALIE": PresetMetadata(
        "GAALIE",
        "Enhance ADCC",
        "G236A/A330L/I332E enhances FcγRIIa/IIIa activation (ADCC & ADCP ↑).",
        adcc="up",
    ),
    "GAALIELS": PresetMetadata(
        "GAALIELS",
        "Enhance ADCC",
        "Combines GAALIE effector enhancement with LS half-life extension.",
        adcc="up",
        fcrn="up",
    ),
    "SDIEALGA": PresetMetadata(
        "SDIEALGA",
        "Enhance ADCC",
        "GASDALIE-class variant boosting both ADCC and ADCP.",
        adcc="up",
    ),
    "SDIEALGALS": PresetMetadata(
        "SDIEALGALS",
        "Enhance ADCC",
        "Combines SDIEALGA ADCC/ADCP enhancement with LS half-life extension.",
        adcc="up",
        fcrn="up",
    ),
    "SDIEGA": PresetMetadata(
        "SDIEGA",
        "Enhance ADCC",
        "G236A/S239D/I332E boosts FcγRIIa and FcγRIIIa affinity.",
        adcc="up",
    ),
    "VLPLL": PresetMetadata(
        "VLPLL",
        "Enhance ADCC",
        "Margetuximab-class Fc variant enhancing FcγRIIIa over inhibitory FcγRIIb.",
        adcc="up",
    ),
    "LPLL": PresetMetadata(
        "LPLL",
        "Enhance ADCC",
        "Core FcγRIIIa-enhancing variant (F243L/R292P/Y300L/P396L).",
        adcc="up",
    ),
    "LPL": PresetMetadata(
        "LPL",
        "Enhance ADCC",
        "F243L/R292P/Y300L variant enhancing ADCC.",
        adcc="up",
    ),
    "LET": PresetMetadata(
        "LET",
        "Enhance ADCC",
        "F243L/D270E/K392T variant modulating FcγRIIIa engagement.",
        adcc="up",
    ),
    "LLEK": PresetMetadata(
        "LLEK",
        "Enhance ADCC",
        "F243L/P247L/D270E/N421K variant enhancing ADCC.",
        adcc="up",
    ),
    "LPIL": PresetMetadata(
        "LPIL",
        "Enhance ADCC",
        "F243L/R292P/V305I/P396L variant enhancing FcγRIIIa affinity.",
        adcc="up",
    ),
    "FALE": PresetMetadata(
        "FALE",
        "Enhance ADCC",
        "G237F/S298A/A330L/I332E variant enhancing ADCC.",
        adcc="up",
    ),
    "EIYE": PresetMetadata(
        "EIYE",
        "Enhance ADCC",
        "S239E/V264I/A330Y/I332E variant enhancing FcγRIIIa binding.",
        adcc="up",
    ),
    "SDIESA": PresetMetadata(
        "SDIESA",
        "Enhance ADCC",
        "S239D/S298A/I332E variant enhancing FcγRIIIa binding.",
        adcc="up",
    ),
    "G236A": PresetMetadata(
        "G236A",
        "Enhance ADCC",
        "Selectively increases FcγRIIa binding and macrophage ADCP.",
        adcc="up",
    ),
    "I332E": PresetMetadata(
        "I332E",
        "Enhance ADCC",
        "Key single-point CH2 mutation enhancing FcγRIIIa affinity.",
        adcc="up",
    ),
    "F243L": PresetMetadata(
        "F243L",
        "Enhance ADCC",
        "Modulates glycan/CH2 conformation to increase FcγRIIIa binding.",
        adcc="up",
    ),
    "Y300L": PresetMetadata(
        "Y300L",
        "Enhance ADCC",
        "CH2 mutation supporting enhanced FcγRIIIa affinity.",
        adcc="up",
    ),
    "D270E": PresetMetadata(
        "D270E",
        "Enhance ADCC",
        "Modulates CH2 C'E loop FcγR selectivity.",
        adcc="up",
    ),
    "S267E": PresetMetadata(
        "S267E",
        "Enhance ADCC",
        "Increases inhibitory FcγRIIb binding (immune-complex / agonism tuning).",
    ),
    "KWES": PresetMetadata(
        "KWES",
        "Enhance CDC",
        "K326W/E333S increases C1q affinity and complement-dependent cytotoxicity (CDC ↑).",
        cdc="up",
    ),
    "FT": PresetMetadata(
        "FT",
        "Enhance CDC",
        "H268F/S324T enhances C1q binding and CDC activity.",
        cdc="up",
    ),
    "AEFT": PresetMetadata(
        "AEFT",
        "Enhance CDC",
        "Combines G236A/I332E with H268F/S324T for dual ADCC & CDC enhancement.",
        adcc="up",
        cdc="up",
    ),
    "ID": PresetMetadata(
        "ID",
        "Enhance ADCC",
        "L234I/L235D lower-hinge FcγR tuning variant.",
    ),
    "IQ": PresetMetadata(
        "IQ",
        "Enhance ADCC",
        "P247I/A339Q variant enhancing FcγRIIIa affinity.",
        adcc="up",
    ),
    "E258H": PresetMetadata(
        "E258H",
        "Extend Half-Life",
        "CH2 loop mutation tuning FcRn pH-dependent interaction.",
        fcrn="up",
    ),
    "T223K": PresetMetadata(
        "T223K",
        "Stability / Hinge",
        "Upper-hinge charge substitution.",
    ),
    "K246H": PresetMetadata(
        "K246H",
        "Glyco / PTM",
        "CH2 N-terminal residue tuning.",
    ),
}


# Known EU functional hotspot map: EU position -> (short_tag, biological_description)
FC_HOTSPOTS: Dict[int, Tuple[str, str]] = {
    228: ("Hinge", "IgG4 core hinge Fab-arm exchange site (S228P stabilizes)"),
    234: ("FcγR", "Lower hinge FcγR binding motif (LALA / silencing)"),
    235: ("FcγR", "Lower hinge FcγR binding motif (LALA / silencing)"),
    236: ("FcγR", "Lower hinge FcγRIIa/ADCP selectivity residue (G236A/R)"),
    237: ("FcγR", "Lower hinge FcγR contact residue"),
    238: ("FcγR", "Lower hinge FcγRIIb/IIIa contact residue"),
    239: ("FcγR", "CH2 FcγRIIIa potency hotspot (S239D)"),
    243: ("FcγR", "CH2 glycan-facing residue tuning FcγRIIIa (F243L)"),
    250: ("FcRn", "CH2 FcRn binding loop (T250Q)"),
    252: ("FcRn", "CH2 FcRn pH 6.0 hotspot & Met oxidation site (M252Y)"),
    254: ("FcRn", "CH2 FcRn pH 6.0 binding loop (S254T)"),
    256: ("FcRn", "CH2 FcRn pH 6.0 binding loop (T256E)"),
    265: ("FcγR", "CH2 BC loop FcγR silencing hotspot (D265A)"),
    267: ("FcγR", "CH2 FcγRIIb affinity residue (S267E)"),
    268: ("C1q", "CH2 C1q/CDC interface residue (H268F)"),
    270: ("FcγR/C1q", "CH2 C'E loop effector interface (D270E)"),
    297: ("N-Glyco", "Conserved Fc N-glycosylation site (N297A removes glycan)"),
    298: ("N-Glyco", "N-glycosylation motif position +1 (S298A)"),
    299: ("N-Glyco", "N-glycosylation motif position +2 (T299)"),
    310: ("FcRn", "Critical FcRn pH-sensing histidine (H310A = FcRnKO)"),
    318: ("C1q", "CH2 C1q binding residue (E318A)"),
    322: ("C1q", "CH2 C1q binding hotspot (K322A reduces CDC)"),
    324: ("C1q", "CH2 C1q binding residue (S324T)"),
    326: ("C1q", "CH2 C1q binding hotspot (K326W enhances CDC)"),
    328: ("FcγR", "CH2 FG loop FcγR interface (L328A)"),
    329: ("FcγR/C1q", "CH2 FG loop P329G effector silencing hotspot"),
    330: ("FcγR", "CH2 FG loop FcγRIIIa/C1q site (A330L)"),
    332: ("FcγR", "CH2 FcγRIIIa potency hotspot (I332E)"),
    333: ("C1q", "CH2 C1q binding hotspot (E333S enhances CDC)"),
    349: ("KiH -S-S-", "CH3 engineered disulfide partner for S354C (Y349C)"),
    354: ("KiH -S-S-", "CH3 engineered disulfide partner for Y349C (S354C)"),
    356: ("CH3 Dimer", "CH3 electrostatic steering site (D356K)"),
    366: ("CH3 Dimer", "CH3 Knob-into-Hole core site (T366W / T366S)"),
    368: ("CH3 Dimer", "CH3 Hole cavity site (L368A)"),
    392: ("CH3 Dimer", "CH3 electrostatic steering site (K392D)"),
    399: ("CH3 Dimer", "CH3 electrostatic steering site (D399K)"),
    407: ("CH3 Dimer", "CH3 Hole cavity core site (Y407V)"),
    409: ("CH3 Dimer", "CH3 dimer interface / IgG4 R409 (K409D)"),
    428: ("FcRn", "CH3 FcRn half-life hotspot & Met oxidation site (M428L)"),
    434: ("FcRn", "CH3 FcRn half-life hotspot (N434S / N434A)"),
    435: ("FcRn", "Critical FcRn pH-sensing histidine (H435A = FcRnKO)"),
}


# One-click goal-oriented Fc engineering recipes
FC_RECIPES: List[Dict[str, Any]] = [
    {
        "id": "silent_pglala",
        "title": "Silent IgG1 (PGLALA)",
        "subtitle": "Complete ADCC/CDC silencing",
        "isotype": "igg1",
        "presets": ["PGLALA"],
        "badge": "ADCC ↓ · CDC ↓",
    },
    {
        "id": "silent_lala_ls",
        "title": "Silent + Long t½ (LALA+LS)",
        "subtitle": "Effector-silent + FcRn half-life extension",
        "isotype": "igg1",
        "presets": ["LALA", "LS"],
        "badge": "Silent · t½ ↑",
    },
    {
        "id": "silent_lala_yte",
        "title": "Silent + YTE (LALA+YTE)",
        "subtitle": "Effector-silent + YTE half-life extension",
        "isotype": "igg1",
        "presets": ["LALA", "YTE"],
        "badge": "Silent · t½ ↑↑",
    },
    {
        "id": "adcc_de",
        "title": "High ADCC (DE)",
        "subtitle": "Strong FcγRIIIa enhancement (S239D/I332E)",
        "isotype": "igg1",
        "presets": ["DE"],
        "badge": "ADCC ↑↑",
    },
    {
        "id": "adcc_gaaliels",
        "title": "ADCC + Long t½ (GAALIELS)",
        "subtitle": "Enhanced ADCC/ADCP + LS half-life",
        "isotype": "igg1",
        "presets": ["GAALIELS"],
        "badge": "ADCC ↑ · t½ ↑",
    },
    {
        "id": "cdc_kwes",
        "title": "High CDC (KWES)",
        "subtitle": "Enhanced C1q complement activation",
        "isotype": "igg1",
        "presets": ["KWES"],
        "badge": "CDC ↑↑",
    },
    {
        "id": "igg4_s228p",
        "title": "Stabilized IgG4 (S228P)",
        "subtitle": "Switch to IgG4 + prevent Fab-arm exchange",
        "isotype": "igg4",
        "presets": ["S228P"],
        "badge": "IgG4 · Hinge",
    },
    {
        "id": "aglyco_n297a",
        "title": "Aglycosylated (N297A)",
        "subtitle": "Remove N297 glycan & effector binding",
        "isotype": None,
        "presets": ["Aglycosylation"],
        "badge": "Aglyco · Silent",
    },
]


def get_preset_metadata(
    label: str,
    mutation_str: str = "",
    is_project_custom: bool = False,
) -> PresetMetadata:
    """Return curated or inferred biological metadata for an Fc mutation preset."""
    if is_project_custom:
        return PresetMetadata(
            label=label,
            category="Project",
            summary=f"Project-specific custom preset ({mutation_str}).",
        )
    if label in _CURATED_PRESET_METADATA:
        return _CURATED_PRESET_METADATA[label]

    low = f"{label} {mutation_str}".lower()
    if any(k in low for k in ("h310a", "h435a", "fcrnko")):
        return PresetMetadata(label, "Reduce Half-Life", f"FcRn-attenuating variant ({mutation_str}).", fcrn="down")
    if any(k in low for k in ("yte", "m252", "m428", "n434", "t250")):
        return PresetMetadata(label, "Extend Half-Life", f"FcRn half-life variant ({mutation_str}).", fcrn="up")
    if any(k in low for k in ("knob", "hole", "t366", "y407", "k392", "d356", "k409", "d399", "s354c", "y349c")):
        return PresetMetadata(label, "Heterodimer", f"CH3 heterodimerization variant ({mutation_str}).", heterodimer=True)
    if any(k in low for k in ("lala", "p329g", "d265a", "k322a", "e318a")):
        return PresetMetadata(label, "Silence Effector", f"Effector-silencing variant ({mutation_str}).", adcc="down", cdc="down")
    if any(k in low for k in ("n297", "agly")):
        return PresetMetadata(label, "Glyco / PTM", f"Glycosylation/PTM variant ({mutation_str}).", adcc="down", cdc="down")
    if any(k in low for k in ("k326", "e333", "h268", "s324")):
        return PresetMetadata(label, "Enhance CDC", f"C1q/CDC-enhancing variant ({mutation_str}).", cdc="up")
    if any(k in low for k in ("s239d", "i332e", "a330l", "g236a", "f243l")):
        return PresetMetadata(label, "Enhance ADCC", f"FcγR/ADCC-enhancing variant ({mutation_str}).", adcc="up")
    if "228" in low:
        return PresetMetadata(label, "Stability / Hinge", f"Hinge stability variant ({mutation_str}).", recommended_isotypes=("igg4",))
    return PresetMetadata(label, "Glyco / PTM", f"Fc variant ({mutation_str}).")


def get_hotspot_info(eu_pos: Optional[int]) -> Optional[Tuple[str, str]]:
    """Return (region_tag, description) for a known functional EU position, or None."""
    if eu_pos is None:
        return None
    return FC_HOTSPOTS.get(eu_pos)


def _merge_effect(cur: str, new: str) -> str:
    if new == "neutral":
        return cur
    if cur == "neutral":
        return new
    if cur == new:
        return cur
    return "mixed"


def summarize_fc_profile(
    chain: PolypeptideChain,
    preset_map: Dict[str, str],
    file_preset_labels: Optional[set] = None,
) -> Dict[str, Any]:
    """Compute a live Fc engineering dashboard summary for a polypeptide chain."""
    file_labels = file_preset_labels or set()
    const_doms = [d for d in chain.domains if d.domain_type in CONSTANT_DOMAIN_EU_RANGES]
    fc_doms = [d for d in chain.domains if d.domain_type in ("Hinge", "CH2", "CH3")]
    has_fc = bool(fc_doms)

    primary_dom = fc_doms[0] if fc_doms else (const_doms[0] if const_doms else None)
    isotype = primary_dom.isotype if primary_dom else ""
    allotype = primary_dom.allotype if primary_dom else ""

    active_presets = sorted({p for d in chain.domains for p in d.applied_presets})
    custom_muts = [(d.domain_type, d.id, tok) for d in chain.domains for tok in d.custom_mutations]

    adcc = "neutral"
    cdc = "neutral"
    fcrn = "neutral"
    heterodimer = False

    # Base IgG2/IgG4 naturally have lower effector function than IgG1
    if isotype == "igg4":
        adcc = "low (IgG4)"
        cdc = "low (IgG4)"
    elif isotype == "igg2":
        adcc = "low (IgG2)"
        cdc = "low (IgG2)"

    preset_adcc = "neutral"
    preset_cdc = "neutral"
    preset_fcrn = "neutral"
    for p_lbl in active_presets:
        meta = get_preset_metadata(
            p_lbl,
            preset_map.get(p_lbl, ""),
            is_project_custom=(p_lbl in file_labels),
        )
        preset_adcc = _merge_effect(preset_adcc, meta.adcc)
        preset_cdc = _merge_effect(preset_cdc, meta.cdc)
        preset_fcrn = _merge_effect(preset_fcrn, meta.fcrn)
        heterodimer = heterodimer or meta.heterodimer

    if preset_adcc != "neutral":
        adcc = preset_adcc
    if preset_cdc != "neutral":
        cdc = preset_cdc
    if preset_fcrn != "neutral":
        fcrn = preset_fcrn

    mut_seq, _, seq_errs = assemble_chain_sequence(chain, preset_map)
    diffs, diff_errs = diff_chain_sequences(chain, preset_map)

    # Also infer effects from ad-hoc custom EU diffs when no preset set them
    for _, pos, _, mut_aa in diffs:
        if pos in (234, 235, 265, 297, 329) and preset_adcc == "neutral":
            adcc = _merge_effect(adcc if adcc in ("up", "down", "mixed") else "neutral", "down")
        if pos in (239, 332) and mut_aa in ("D", "E") and preset_adcc == "neutral":
            adcc = _merge_effect(adcc if adcc in ("up", "down", "mixed") else "neutral", "up")
        if pos in (252, 428, 434) and preset_fcrn == "neutral":
            fcrn = _merge_effect(fcrn, "up")
        if pos in (310, 435) and mut_aa == "A" and preset_fcrn == "neutral":
            fcrn = _merge_effect(fcrn, "down")
        if pos in (366, 368, 407, 349, 354):
            heterodimer = True

    # Compute WT sequence for delta metrics
    wt_parts: List[str] = []
    for d in chain.domains:
        bseq, berr = resolve_domain_base_sequence(d)
        if not berr:
            wt_parts.append(bseq)
    wt_seq = "".join(wt_parts)

    wt_mw, wt_pi = compute_protein_metrics(wt_seq)
    mut_mw, mut_pi = compute_protein_metrics(mut_seq) if not seq_errs else (wt_mw, wt_pi)

    return {
        "has_fc": has_fc,
        "isotype": isotype,
        "allotype": allotype,
        "active_presets": active_presets,
        "custom_mutations": custom_muts,
        "diffs": diffs,
        "errors": seq_errs + diff_errs,
        "adcc": adcc,
        "cdc": cdc,
        "fcrn": fcrn,
        "heterodimer": heterodimer,
        "length": len(mut_seq),
        "wt_mw": wt_mw,
        "mut_mw": mut_mw,
        "delta_mw": round(mut_mw - wt_mw, 3),
        "wt_pi": wt_pi,
        "mut_pi": mut_pi,
        "delta_pi": round(mut_pi - wt_pi, 2),
    }


def build_eu_residue_rows(
    chain: PolypeptideChain,
    preset_map: Dict[str, str],
    domain_types: Optional[Tuple[str, ...]] = None,
) -> List[Dict[str, Any]]:
    """Build per-residue EU-numbered entries for constant domains on a chain."""
    allowed = domain_types or ("Hinge", "CH2", "CH3", "CH1")
    entries: List[Dict[str, Any]] = []

    for dom in chain.domains:
        if dom.domain_type not in CONSTANT_DOMAIN_EU_RANGES:
            continue
        if dom.domain_type not in allowed:
            continue
        wt_seq, wt_err = resolve_domain_base_sequence(dom)
        mut_seq, mut_errs = resolve_domain_sequence(dom, preset_map)
        if wt_err or mut_errs or not wt_seq:
            continue

        # Collect PTM positions on this domain's mutated sequence
        ptm_by_pos: Dict[int, List[str]] = {}
        for cat, pos_list in analyze_developability(mut_seq, dom.isotype, domain_type=dom.domain_type):
            short_cat = cat.split("(")[0].strip()
            for p in pos_list:
                ptm_by_pos.setdefault(p, []).append(short_cat)

        for idx, (wt_aa, mut_aa) in enumerate(zip(wt_seq, mut_seq)):
            eu_pos = get_domain_eu_position(idx, dom.domain_type, dom.isotype)
            if eu_pos is None:
                continue
            hotspot = get_hotspot_info(eu_pos)
            entries.append(
                {
                    "domain_id": dom.id,
                    "domain_type": dom.domain_type,
                    "isotype": dom.isotype,
                    "eu_pos": eu_pos,
                    "wt_aa": wt_aa,
                    "mut_aa": mut_aa,
                    "is_mutated": wt_aa != mut_aa,
                    "hotspot_tag": hotspot[0] if hotspot else "",
                    "hotspot_desc": hotspot[1] if hotspot else "",
                    "ptm_tags": ptm_by_pos.get(eu_pos, []),
                }
            )
    return entries
