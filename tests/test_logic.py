import pytest
from src.main import get_residue_index, parse_mutation, apply_mutations

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
    _, errors = apply_mutations("ANY_SEQ", "S223P", "igg2")
    assert any("Gap" in e for e in errors)
