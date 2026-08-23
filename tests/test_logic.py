import pytest
from fc_engineer.core import get_residue_index, parse_mutation, apply_mutations
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
