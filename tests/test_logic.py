import pytest
from fc_engineer.core import (
    get_residue_index,
    get_eu_position,
    parse_mutation,
    apply_mutations,
    diff_sequences,
    generate_batch,
    analyze_developability,
)
from fc_engineer.config import SEQUENCES, COMMON_MUTATIONS

# 로드맵: developability / PTM 관련 신규 프리셋
DEVELOPABILITY_PRESETS = [
    "N434A",
    "T250Q/M428L",
    "D265A",
    "K322A",
    "E318A",
    "G236R",
    "S267E",
    "H310A/H435A",
    "P329G/L234A/L235A",
]

def test_get_residue_index_igg1():
    # IgG1은 EU 118부터 연속적임
    assert get_residue_index(118, "igg1") == 0
    assert get_residue_index(119, "igg1") == 1
    assert get_residue_index(447, "igg1") == 329

def test_get_residue_index_igg2_gaps():
    # IgG2/4는 EU 223-225 부위가 결실(None)
    assert get_residue_index(222, "igg2") == 104
    assert get_residue_index(223, "igg2") is None
    assert get_residue_index(224, "igg2") is None
    assert get_residue_index(225, "igg2") is None
    assert get_residue_index(226, "igg2") == 105  # 104 다음은 105 (3개 점프)

def test_get_residue_index_igg3():
    # IgG3는 확장 힌지 보유: 222까지 정방향, 223-230은 Gap, 231+는 CH2/CH3 정렬(+47)
    assert get_residue_index(118, "igg3") == 0
    assert get_residue_index(222, "igg3") == 104
    for pos in range(223, 231):
        assert get_residue_index(pos, "igg3") is None
    assert get_residue_index(231, "igg3") == 160  # APELL 시작 위치와 일치
    assert get_residue_index(297, "igg3") == 226  # N297 글리코실화 부위
    assert get_residue_index(447, "igg3") == 376  # C-말단 (len=377, 마지막 인덱스)

def test_apply_mutations_igg3():
    seq = SEQUENCES.get("igg3", {}).get("WT(P01860-1)", "")
    assert len(seq) == 377  # IgG1(330) 대비 힌지 +47
    # N297A (Aglycosylation) 적용 검증
    result, errors = apply_mutations(seq, "N297A", "igg3")
    assert not errors
    idx = get_residue_index(297, "igg3")
    assert seq[idx] == "N"
    assert result[idx] == "A"

@pytest.mark.parametrize("value", DEVELOPABILITY_PRESETS)
def test_developability_presets_apply_cleanly(value):
    # developability 프리셋은 IgG1 WT에 오류 없이 적용되어야 함
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    result, errors = apply_mutations(seq, value, "igg1")
    assert not errors
    assert result != seq  # 실제 변이가 반영됨

def test_all_preset_values_are_parseable():
    # mutants.yaml의 모든 value는 유효한 변이 문법(들)이어야 함 (데이터 무결성)
    for item in COMMON_MUTATIONS:
        for m in item["value"].split("/"):
            wt_aa, pos, mut_aa = parse_mutation(m)  # ValueError 발생 시 실패
            assert wt_aa in "ACDEFGHIKLMNPQRSTVWY"
            assert mut_aa in "ACDEFGHIKLMNPQRSTVWY"

@pytest.mark.parametrize("iso", ["igg1", "igg2", "igg4", "igg3"])
def test_get_eu_position_roundtrip(iso):
    # 인덱스 <-> EU 위치 역변환이 각 Isotype에서 일관되어야 함
    for pos in [118, 150, 200, 222]:
        assert get_eu_position(get_residue_index(pos, iso), iso) == pos

def test_get_eu_position_after_gap():
    # Gap 이후 구간의 역변환 (IgG2: 226 -> idx 105, IgG3: 231 -> idx 160)
    assert get_eu_position(105, "igg2") == 226
    assert get_eu_position(160, "igg3") == 231
    assert get_eu_position(376, "igg3") == 447

def test_diff_sequences_igg1():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    mut, _ = apply_mutations(seq, "L234A/L235A/N297A", "igg1")
    diffs = diff_sequences(seq, mut, "igg1")
    assert diffs == [(234, "L", "A"), (235, "L", "A"), (297, "N", "A")]

def test_diff_sequences_igg3():
    seq = SEQUENCES["igg3"]["WT(P01860-1)"]
    mut, _ = apply_mutations(seq, "N297A", "igg3")
    assert diff_sequences(seq, mut, "igg3") == [(297, "N", "A")]

def test_diff_sequences_no_change():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    assert diff_sequences(seq, seq, "igg1") == []

# --- Batch generation ---

def test_generate_batch_basic():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    results = generate_batch(seq, "L234A/L235A, N297A", "igg1")
    assert len(results) == 2
    combo1, seq1, err1 = results[0]
    assert combo1 == "L234A/L235A" and not err1
    assert seq1[116:118] == "AA"
    combo2, seq2, err2 = results[1]
    assert combo2 == "N297A" and not err2
    idx297 = get_residue_index(297, "igg1")
    assert seq2[idx297] == "A"

def test_generate_batch_empty_returns_nothing():
    seq = "ASTK"
    assert generate_batch(seq, "", "igg1") == []
    assert generate_batch(seq, " , , ", "igg1") == []

def test_generate_batch_invalid_chars():
    seq = "ASTK"
    _, _, errors = generate_batch(seq, "A118X; rm -rf", "igg1")[0]
    assert any("invalid characters" in e for e in errors)

def test_generate_batch_too_many_items():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    batch_str = ",".join(["L234A"] * 21)
    results = generate_batch(seq, batch_str, "igg1")
    assert len(results) == 1
    assert "Maximum of 20 batch combinations allowed." in results[0][2][0]

def test_generate_batch_length_limit():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    results = generate_batch(seq, "A118X," * 400, "igg1")  # > 2000 chars
    assert len(results) == 1
    assert "exceeds maximum length of 2000" in results[0][2][0]

def test_generate_batch_propagates_per_item_errors():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    results = generate_batch(seq, "N297A, C118X", "igg1")
    # 두 번째 조합은 WT 불일치 오류를 포함하되 전체 배치는 계속 진행되어야 함
    assert not results[0][2]
    assert any("Expected 'C'" in e for e in results[1][2])

# --- Developability analysis ---

def test_analyze_developability_igg1_wt():
    # IgG1 WT에서 알려진 liability 위치들을 정확히 탐지해야 함
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    report = dict(analyze_developability(seq, "igg1"))
    assert report["N-glycosylation motif (N-X-S/T)"] == [297]  # Fc N-글리코실화 부위
    assert 252 in report["Oxidation-prone Met (M252/M428 = FcRn sites)"]
    assert 428 in report["Oxidation-prone Met (M252/M428 = FcRn sites)"]

def test_analyze_developability_n297a_removes_glycosylation():
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    mut, _ = apply_mutations(seq, "N297A", "igg1")
    report = dict(analyze_developability(mut, "igg1"))
    assert report["N-glycosylation motif (N-X-S/T)"] == []

def test_analyze_developability_invalid_input():
    assert analyze_developability("", "igg1") == []
    assert analyze_developability(None, "igg1") == []
    assert analyze_developability("ASTK", "") == []

def test_analyze_developability_positions_match_eu():
    # 모든 보고된 위치는 유효한 EU 번호(>=118)여야 하고, 해당 위치의 잔기와 일치해야 함
    seq = SEQUENCES["igg1"]["WT(P01857-1)"]
    for category, positions in analyze_developability(seq, "igg1"):
        for pos in positions:
            assert pos >= 118
            idx = get_residue_index(pos, "igg1")
            if "Met" in category:
                assert seq[idx] == "M"

def test_parse_mutation():
    assert parse_mutation("S228P") == ("S", 228, "P")
    assert parse_mutation("L234A") == ("L", 234, "A")
    with pytest.raises(Exception):
        parse_mutation("INVALID")

def test_apply_mutations_single():
    seq = "ABCDEFGHIJ"  # Dummy sequence (10 aa)
    # EU_START(118) 기준, 118번은 'A', 119번은 'B' ...
    # 119번 B를 X로 변경 (B119X)
    result, errors = apply_mutations(seq, "B119X", "igg1")
    assert result == "AXCDEFGHIJ"
    assert not errors

def test_apply_mutations_multiple():
    seq = "ASTKGPSVFPLAPSSK"
    # A118X / S119Y / T120Z 적용
    result, errors = apply_mutations(seq, "A118X/S119Y/T120Z", "igg1")
    assert result.startswith("XYZ")
    assert not errors

def test_apply_mutations_errors():
    seq = "ASTK" # A(118), S(119), T(120), K(121)
    
    # 1. 기대값 불일치 (기존 위치에 A가 아닌 C가 있다고 주장)
    _, errors = apply_mutations(seq, "C118X", "igg1")
    assert any("Expected 'C', found 'A'" in e for e in errors)
    
    # 2. 범위를 벗어난 위치
    _, errors = apply_mutations(seq, "A500X", "igg1")
    assert any("out of range" in e for e in errors)
    
    # 3. Gap 부위에 변이 시도 (IgG2의 223번)
    _, errors = apply_mutations("ANYSEQ", "S223P", "igg2")
    assert any("Gap" in e for e in errors)

    # 4. Too many mutations limit
    many_muts = "/".join([f"A{118+i}X" for i in range(51)])
    _, errors = apply_mutations(seq, many_muts, "igg1")
    assert "Maximum of 50 mutations allowed" in errors[0]

    # 5. Mutation string too long
    _, errors = apply_mutations(seq, "A118X" * 3, "igg1")
    assert any("Invalid mutation" in e for e in errors)

def test_apply_mutations_length_limit():
    seq = "ASTKGPSVFPLAPSSK"
    long_input = "A118X/" * 200  # length > 1000
    _, errors = apply_mutations(seq, long_input, "igg1")
    assert "Error: Mutation string exceeds maximum length of 1000 characters." in errors[0]


# --- Config domain templates & domain-aware core.py tests ---

from fc_engineer.config import (
    DEFAULT_LINKERS,
    LIGHT_CONSTANT_DOMAINS,
    VARIABLE_DOMAINS,
    load_domain_templates,
)
from fc_engineer.core import (
    CONSTANT_DOMAIN_EU_RANGES,
    apply_domain_mutations,
    get_domain_eu_position,
    get_domain_residue_index,
    get_domain_slice_bounds,
    parse_mutation_list,
    slice_constant_domain,
)
from fc_engineer.designer import (
    AntibodyProject,
    DomainBlock,
    InterChainPairing,
    PolypeptideChain,
    SCENE_PADDING,
    add_file_custom_mutation,
    analyze_chain_developability,
    assemble_chain_sequence,
    build_format_template,
    build_preset_map,
    build_schematic_scene,
    diff_chain_sequences,
    get_effective_domain_mutations,
    reflow_chain_coordinates,
    remove_chain_from_project,
    remove_domain_from_chain,
    remove_preset_from_chain,
    resolve_domain_base_sequence,
    resolve_domain_sequence,
    route_mutations_to_chain,
)


def test_load_domain_templates_non_empty():
    var_d, light_d, linkers = load_domain_templates()
    assert var_d == VARIABLE_DOMAINS
    assert light_d == LIGHT_CONSTANT_DOMAINS
    assert linkers == DEFAULT_LINKERS
    assert set(var_d.keys()) == {"VH_A", "VL_A", "VH_B", "VL_B"}
    assert set(light_d.keys()) == {"CL_kappa", "CL_lambda"}
    assert set(linkers.keys()) == {"G4S_3", "G4S_2", "G4S_1"}


@pytest.mark.parametrize("iso", ["igg1", "igg2", "igg3", "igg4"])
def test_slice_constant_domain_concatenation_matches_full_sequence(iso):
    allotype = next(iter(SEQUENCES[iso]))
    full_seq = SEQUENCES[iso][allotype]
    ch1 = slice_constant_domain(full_seq, "CH1", iso)
    hinge = slice_constant_domain(full_seq, "Hinge", iso)
    ch2 = slice_constant_domain(full_seq, "CH2", iso)
    ch3 = slice_constant_domain(full_seq, "CH3", iso)
    assert ch1 and hinge and ch2 and ch3
    assert ch1 + hinge + ch2 + ch3 == full_seq


def test_get_eu_position_negative_and_igg3_extended_hinge():
    assert get_eu_position(-1, "igg1") is None
    for idx in range(105, 160):
        assert get_eu_position(idx, "igg3") is None
    assert get_eu_position(104, "igg3") == 222
    assert get_eu_position(160, "igg3") == 231


def test_domain_eu_index_and_igg3_hinge_roundtrip():
    for local_idx in range(7):  # EU 216..222
        assert get_domain_eu_position(local_idx, "Hinge", "igg3") == 216 + local_idx
    for local_idx in range(7, 62):  # Global 105..159 repeat region
        assert get_domain_eu_position(local_idx, "Hinge", "igg3") is None
    assert get_domain_eu_position(-1, "CH2", "igg1") is None
    assert get_domain_eu_position(10, "Unknown", "igg1") is None
    assert get_domain_eu_position(999, "CH2", "igg1") is None

    assert get_domain_residue_index(234, "CH2", "igg1") == 3
    assert get_domain_residue_index(366, "CH3", "igg1") == 25
    assert get_domain_residue_index(234, "CH3", "igg1") is None
    assert get_domain_slice_bounds("Unknown", "igg1", 330) is None


def test_apply_domain_mutations_dedup_and_conflict():
    ch2 = slice_constant_domain(SEQUENCES["igg1"]["WT(P01857-1)"], "CH2", "igg1")
    mut_seq, errs = apply_domain_mutations(ch2, "L234A/L234A/L235A", "CH2", "igg1")
    assert not errs
    assert mut_seq != ch2

    conflict_seq, conflict_errs = apply_domain_mutations(
        ch2, "L234A/L234F", "CH2", "igg1"
    )
    assert conflict_seq == ch2
    assert any("Conflicting mutations" in e for e in conflict_errs)

    oor_seq, oor_errs = apply_domain_mutations(ch2, "T366W", "CH2", "igg1")
    assert oor_seq == ch2
    assert any("out of range" in e for e in oor_errs)


def test_diff_and_developability_on_isolated_domains():
    ch2 = slice_constant_domain(SEQUENCES["igg1"]["WT(P01857-1)"], "CH2", "igg1")
    ch2_mut, _ = apply_domain_mutations(ch2, "M252Y/N297A", "CH2", "igg1")
    assert diff_sequences(ch2, ch2_mut, "igg1", domain_type="CH2") == [
        (252, "M", "Y"),
        (297, "N", "A"),
    ]
    dev_wt = dict(analyze_developability(ch2, "igg1", domain_type="CH2"))
    dev_mut = dict(analyze_developability(ch2_mut, "igg1", domain_type="CH2"))
    assert dev_wt["N-glycosylation motif (N-X-S/T)"] == [297]
    assert 252 in dev_wt["Oxidation-prone Met (M252/M428 = FcRn sites)"]
    assert dev_mut["N-glycosylation motif (N-X-S/T)"] == []


# --- designer.py tests ---


@pytest.mark.parametrize(
    "fmt,expected_chains",
    [
        ("igg", {"LC_1", "HC_1", "HC_2", "LC_2"}),
        ("crossmab", {"LC_A", "HC_Knob", "HC_Hole_Cross", "LC_B_Cross"}),
        ("igg_scfv", {"LC_1", "HC_1", "HC_2", "LC_2"}),
        ("scfv_fc", {"Chain_A", "Chain_B"}),
        ("dvd_ig", {"LC_1", "HC_1", "HC_2", "LC_2"}),
        ("bite", {"Chain_1"}),
        ("two_plus_one", {"LC_A1", "LC_B", "HC_2Fab", "HC_1Fab", "LC_A2"}),
        ("trispecific", {"LC_1", "HC_1_Knob", "HC_2_Hole", "LC_2_Cross"}),
        ("fab_igg", {"LC_1", "HC_1", "HC_2", "LC_2"}),
        ("custom", {"Chain_1"}),
    ],
)
def test_build_format_template_all_formats_assemble_cleanly(fmt, expected_chains):
    proj = build_format_template(fmt)
    assert proj.format_type == fmt
    assert {c.id for c in proj.chains} == expected_chains
    preset_map = build_preset_map(proj)
    for chain in proj.chains:
        seq, spans, errs = assemble_chain_sequence(chain, preset_map)
        assert not errs
        assert seq
        assert len(spans) == len(chain.domains)

    if fmt in ("igg_scfv", "scfv_fc", "bite"):
        assert any(
            p.chain_a_id == p.chain_b_id and p.pairing_type == "VH_VL"
            for p in proj.pairings
        )
    if fmt == "crossmab":
        hck = next(c for c in proj.chains if c.id == "HC_Knob")
        hch = next(c for c in proj.chains if c.id == "HC_Hole_Cross")
        assert "Knob" in hck.domains[-1].applied_presets
        assert "Hole" in hch.domains[-1].applied_presets


def test_build_format_template_invalid_inputs():
    with pytest.raises(ValueError):
        build_format_template("unknown_fmt")
    with pytest.raises(ValueError):
        build_format_template("igg", isotype="igg9")
    with pytest.raises(ValueError):
        build_format_template("igg", isotype="igg1", allotype="missing_allotype")


def test_domain_custom_sequence_and_template_validation():
    vh_ok = DomainBlock(id="d1", domain_type="VH", custom_sequence="EVQLVESGGGLV")
    seq, errs = resolve_domain_base_sequence(vh_ok)
    assert not errs and seq == "EVQLVESGGGLV"

    # Cross-type template key rejected even when custom_sequence is present (IR3-N2)
    vh_bad_tkey = DomainBlock(
        id="d2", domain_type="VH", template_key="VL_A", custom_sequence="EVQL"
    )
    _, errs = resolve_domain_base_sequence(vh_bad_tkey)
    assert errs

    g4s_long = DomainBlock(id="d3", domain_type="G4S", custom_sequence="G" * 2001)
    _, errs = resolve_domain_base_sequence(g4s_long)
    assert any("exceeds maximum length" in e for e in errs)

    g4s_bad_chars = DomainBlock(id="d4", domain_type="G4S", custom_sequence="ggggs!")
    _, errs = resolve_domain_base_sequence(g4s_bad_chars)
    assert errs

    cl_custom = DomainBlock(id="d5", domain_type="CL", custom_sequence="RTVAAP")
    _, errs = resolve_domain_base_sequence(cl_custom)
    assert errs

    ch2_custom = DomainBlock(id="d6", domain_type="CH2", custom_sequence="APELL")
    _, errs = resolve_domain_base_sequence(ch2_custom)
    assert errs


@pytest.mark.parametrize(
    "preset_label",
    [
        "GAALIELS",
        "SDIEALGALS",
        "VLPLL",
        "LPLL",
        "LET",
        "LLEK",
        "LPIL",
        "QL",
        "FcRnKO",
        "LALA",
        "Knob",
        "Hole",
    ],
)
def test_cross_domain_preset_routing_and_atomic_removal(preset_label):
    proj = build_format_template("igg")
    hc1 = next(c for c in proj.chains if c.id == "HC_1")
    preset_map = build_preset_map(proj)
    errs = route_mutations_to_chain(
        hc1, preset_map[preset_label], preset_label=preset_label, preset_map=preset_map
    )
    assert not errs
    diffs, derrs = diff_chain_sequences(hc1, preset_map)
    assert not derrs
    assert diffs

    # Remove preset and confirm return to WT
    remove_preset_from_chain(hc1, preset_label)
    diffs_after, _ = diff_chain_sequences(hc1, preset_map)
    assert diffs_after == []


def test_overlapping_presets_and_custom_mutations_preserved_on_remove():
    proj = build_format_template("igg")
    hc1 = next(c for c in proj.chains if c.id == "HC_1")
    preset_map = build_preset_map(proj)
    assert not route_mutations_to_chain(hc1, "", preset_label="LALA", preset_map=preset_map)
    assert not route_mutations_to_chain(hc1, "", preset_label="PGLALA", preset_map=preset_map)
    assert not route_mutations_to_chain(hc1, "N297A", preset_map=preset_map)

    ch2 = next(d for d in hc1.domains if d.domain_type == "CH2")
    eff, errs = get_effective_domain_mutations(ch2, preset_map)
    assert not errs
    assert eff == ["L234A", "L235A", "P329G", "N297A"]

    # Removing PGLALA leaves LALA (L234A/L235A) and N297A intact!
    remove_preset_from_chain(hc1, "PGLALA")
    eff_after, _ = get_effective_domain_mutations(ch2, preset_map)
    assert eff_after == ["L234A", "L235A", "N297A"]


def test_route_mutations_atomic_rollback_and_edge_errors():
    # IgG3 has R435 instead of H435; FcRnKO (H310A/H435A) must fail atomically without mutating CH2!
    proj_igg3 = build_format_template("igg", isotype="igg3")
    hc3 = next(c for c in proj_igg3.chains if c.id == "HC_1")
    errs = route_mutations_to_chain(hc3, "H310A/H435A", preset_label="FcRnKO")
    assert errs
    assert all(not d.applied_presets and not d.mutation_history for d in hc3.domains)

    # Cetuximab lacks C-terminal K447
    proj_cetux = build_format_template("igg", isotype="igg1", allotype="cetuximab")
    hc_c = next(c for c in proj_cetux.chains if c.id == "HC_1")
    errs_c = route_mutations_to_chain(hc_c, "L234A/K447A")
    assert any("out of range" in e for e in errs_c)
    assert all(not d.custom_mutations for d in hc_c.domains)

    # Mixed valid Fc mutation + out-of-range EU position (>447) must fail atomically (IR3-M4)
    proj = build_format_template("igg")
    hc1 = next(c for c in proj.chains if c.id == "HC_1")
    errs_oor = route_mutations_to_chain(hc1, "L234A/A500X")
    assert any("outside constant domain EU range" in e for e in errs_oor)
    assert all(not d.custom_mutations for d in hc1.domains)

    # Mismatched mutants_str vs preset_label rejected (IR3-M4)
    errs_div = route_mutations_to_chain(hc1, "N297A", preset_label="LALA")
    assert errs_div


def test_chain_diff_and_developability_independent_of_vh_g4s():
    proj_igg = build_format_template("igg")
    hc_igg = next(c for c in proj_igg.chains if c.id == "HC_1")
    proj_dvd = build_format_template("dvd_ig")  # Extra VH + G4S at N-terminus!
    hc_dvd = next(c for c in proj_dvd.chains if c.id == "HC_1")

    assert not route_mutations_to_chain(hc_igg, "", preset_label="YTE")
    assert not route_mutations_to_chain(hc_dvd, "", preset_label="YTE")

    diffs_igg, _ = diff_chain_sequences(hc_igg)
    diffs_dvd, _ = diff_chain_sequences(hc_dvd)
    assert diffs_igg == diffs_dvd == [
        ("CH2", 252, "M", "Y"),
        ("CH2", 254, "S", "T"),
        ("CH2", 256, "T", "E"),
    ]
    dev_igg, _ = analyze_chain_developability(hc_igg)
    dev_dvd, _ = analyze_chain_developability(hc_dvd)
    assert dev_igg == dev_dvd


def test_mutation_history_fifo_and_length_bounds():
    proj = build_format_template("igg")
    hc1 = next(c for c in proj.chains if c.id == "HC_1")
    for idx in range(55):
        mut_str = "L234A" if idx % 2 == 0 else "N297A"
        assert not route_mutations_to_chain(hc1, mut_str)
    ch2 = next(d for d in hc1.domains if d.domain_type == "CH2")
    assert len(ch2.mutation_history) == 50
    assert all(len(entry) <= 200 for entry in ch2.mutation_history)


def test_reflow_pruning_and_schematic_scene():
    proj = build_format_template("igg")
    hc1 = next(c for c in proj.chains if c.id == "HC_1")
    assert not route_mutations_to_chain(hc1, "", preset_label="Knob")
    assert not route_mutations_to_chain(hc1, "K360E/Q362E/Q419K")

    ch3 = next(d for d in hc1.domains if d.domain_type == "CH3")
    ch3.x = 777.0
    ch3.y = 888.0
    ch3.pinned = True
    reflow_chain_coordinates(hc1)
    assert (ch3.x, ch3.y) == (777.0, 888.0)

    scene = build_schematic_scene(proj)
    assert scene.bounds[2] >= 777.0 + 52.0 + SCENE_PADDING
    assert scene.bounds[3] >= 888.0 + 88.0 + SCENE_PADDING
    ch3_node = next(n for n in scene.nodes if n.domain_id == ch3.id)
    assert ch3_node.badge_text == "Knob, K360E, Q362E, +1"

    # Remove domain & chain prune pairings cleanly
    remove_domain_from_chain(proj, "HC_1", "hc1_vh")
    assert not any(
        (p.chain_a_id == "HC_1" and p.domain_a_id == "hc1_vh")
        or (p.chain_b_id == "HC_1" and p.domain_b_id == "hc1_vh")
        for p in proj.pairings
    )
    remove_chain_from_project(proj, "HC_2")
    assert all(
        p.chain_a_id != "HC_2" and p.chain_b_id != "HC_2" for p in proj.pairings
    )

    # Boolean / non-finite coord rejection in scene builder
    ch3.x = True  # type: ignore[assignment]
    with pytest.raises(ValueError):
        build_schematic_scene(proj)


def test_add_file_custom_mutation_validation_and_mutants_yaml_untouched(tmp_path):
    from pathlib import Path

    mutants_yaml = Path(__file__).resolve().parents[1] / "src" / "fc_engineer" / "data" / "mutants.yaml"
    before_bytes = mutants_yaml.read_bytes()

    proj = build_format_template("igg")
    assert not add_file_custom_mutation(proj, "ProjYTE_Alt", "M252Y/T256E")
    assert proj.file_custom_mutations == [{"label": "ProjYTE_Alt", "value": "M252Y/T256E"}]

    # Duplicate & shadow rejections
    assert add_file_custom_mutation(proj, "ProjYTE_Alt", "N297A")
    assert add_file_custom_mutation(proj, "LALA", "L234A/L235A")
    assert add_file_custom_mutation(proj, "Bad\nLabel", "N297A")
    assert add_file_custom_mutation(proj, "X" * 41, "N297A")
    assert add_file_custom_mutation(proj, "ValidLabel", "")

    assert mutants_yaml.read_bytes() == before_bytes


def test_package_version_alignment():
    """Verify pyproject.toml version and fc_engineer.__version__ stay aligned."""
    import tomllib
    from pathlib import Path
    import fc_engineer

    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        project_version = tomllib.load(f)["project"]["version"]

    assert project_version == "0.0.1"
    assert fc_engineer.__version__ == project_version


def test_disulfide_detection_and_multispecific_helpers():
    from fc_engineer.designer import (
        apply_heterodimer_pair,
        clone_chain_in_project,
        swap_crossmab_ch1_cl,
    )

    proj = build_format_template("igg")
    assert not apply_heterodimer_pair(proj, "HC_1", "HC_2", mode="kih_ss")
    scene = build_schematic_scene(proj)
    ch3_edges = [e for e in scene.pairing_edges if e.pairing_type == "CH3_DIMER"]
    assert len(ch3_edges) == 1
    assert ch3_edges[0].is_disulfide is True

    hinge_edges = [e for e in scene.pairing_edges if e.pairing_type == "HINGE_DISULFIDE"]
    assert hinge_edges and all(e.is_disulfide for e in hinge_edges)

    cloned = clone_chain_in_project(proj, "LC_1")
    assert cloned.id == "Chain_1"
    assert len(cloned.domains) == 2

    # CrossMab CH1 <-> CL swap between HC_2 and LC_2
    assert not swap_crossmab_ch1_cl(proj, "HC_2", "LC_2")
    hc2 = next(c for c in proj.chains if c.id == "HC_2")
    lc2 = next(c for c in proj.chains if c.id == "LC_2")
    assert any(d.domain_type == "CL" for d in hc2.domains)
    assert any(d.domain_type == "CH1" for d in lc2.domains)


def test_metrics_wt_lookup_and_charge_pair_ux_helpers(tmp_path):
    from fc_engineer import config
    from fc_engineer.core import compute_protein_metrics, lookup_wt_residue
    from fc_engineer.designer import (
        add_file_custom_mutation,
        analyze_chain_developability,
        apply_heterodimer_pair,
        auto_detect_standard_pairings,
        delete_file_custom_mutation,
        remove_custom_mutation_from_domain,
        reset_chain_mutations,
    )

    mw, pi = compute_protein_metrics("ACDEFGHIKLMNPQRSTVWY")
    assert mw > 2.0
    assert 4.0 < pi < 9.0

    igg1_seq = next(iter(config.SEQUENCES["igg1"].values()))
    assert lookup_wt_residue(igg1_seq, 234, "igg1") == "L"
    assert lookup_wt_residue(igg1_seq, 252, "igg1") == "M"
    assert lookup_wt_residue(igg1_seq, 999, "igg1") is None

    proj = build_format_template("igg")
    assert not apply_heterodimer_pair(proj, "HC_1", "HC_2", mode="charge_pair")
    hc1 = next(c for c in proj.chains if c.id == "HC_1")
    hc2 = next(c for c in proj.chains if c.id == "HC_2")
    assert any("ChargePair_A" in d.applied_presets for d in hc1.domains)
    assert any("ChargePair_B" in d.applied_presets for d in hc2.domains)

    # Clear pairings and auto-detect
    proj.pairings.clear()
    added = auto_detect_standard_pairings(proj)
    assert added >= 2

    # Custom mutation removal & reset chain
    ch2 = next(d for d in hc1.domains if d.domain_type == "CH2")
    ch2.custom_mutations.append("L234A")
    remove_custom_mutation_from_domain(ch2, "L234A")
    assert "L234A" not in ch2.custom_mutations

    assert not add_file_custom_mutation(proj, "MyMut", "M252Y")
    delete_file_custom_mutation(proj, "MyMut")
    assert not any(m.get("label") == "MyMut" for m in proj.file_custom_mutations)

    reset_chain_mutations(hc1)
    assert all(not d.applied_presets and not d.custom_mutations for d in hc1.domains)

    # Variable developability scan toggle
    liab_var, dev_errs = analyze_chain_developability(hc1, include_variable=True)
    assert not dev_errs
    assert len(liab_var) >= 1

    # User prefs roundtrip
    prefs_path = tmp_path / "prefs.yaml"
    config.save_user_prefs({"theme_mode": "dark", "recent_files": ["/a.ab"]}, str(prefs_path))
    loaded = config.load_user_prefs(str(prefs_path))
    assert loaded["theme_mode"] == "dark"
    assert loaded["recent_files"] == ["/a.ab"]



