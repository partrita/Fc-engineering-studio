import pytest

from fc_engineer.config import SEQUENCES
from fc_engineer.core import apply_mutations
from fc_engineer.exporters import (
    format_genbank,
    format_csv_report,
    sanitize_locus,
)

IGG1 = SEQUENCES["igg1"]["WT(P01857-1)"]


class TestSanitizeLocus:
    def test_replaces_invalid_characters(self):
        assert sanitize_locus('IGG1_Wt(p01857-1)_X/Y') == "IGG1_Wt_p01857-1"

    def test_truncates_to_16_chars(self):
        assert len(sanitize_locus("A" * 100)) == 16

    def test_non_string_falls_back(self):
        assert sanitize_locus(None) == "UNKNOWN"
        assert sanitize_locus("") == "UNKNOWN"


class TestFormatGenbank:
    def test_basic_structure(self):
        record = format_genbank("TEST_SEQ", IGG1, definition="Test record")
        lines = record.split("\n")
        assert lines[0].startswith("LOCUS       TEST_SEQ")
        assert f"{len(IGG1):>9} aa" in lines[0]
        assert " aa    linear   UNK " in lines[0]
        assert any(l.startswith("DEFINITION  Test record.") for l in lines)
        assert any(l.startswith("ORIGIN") for l in lines)
        assert record.endswith("//\n")
        assert record.strip().split("\n")[-1] == "//"
        assert any(f"1..{len(IGG1)}" in l for l in lines)

    def test_sequence_wrapping_60_per_line(self):
        record = format_genbank("T", IGG1)
        origin_lines = [l for l in record.split("\n") if l[:9].strip().isdigit()]
        assert len(origin_lines) == (len(IGG1) + 59) // 60
        # 첫 줄: 위치 1, 60자 서열이 10자 블록으로 분리
        assert origin_lines[0].startswith("        1 ASTK")
        joined = "".join(l[10:] for l in origin_lines).replace(" ", "")
        assert joined == IGG1

    def test_invalid_inputs_return_empty(self):
        assert format_genbank("T", "") == ""
        assert format_genbank("T", "atgc!123") == ""  # 소문자/비AA 문자 거부
        assert format_genbank("T", "A" * 20001) == ""
        # 빈 locus는 UNKNOWN으로 대체되어 생성됨 (거부 아님)
        assert format_genbank("", IGG1).startswith("LOCUS       UNKNOWN")


class TestFormatCsvReport:
    def setup_method(self):
        self.mut, _ = apply_mutations(IGG1, "L234A/L235A/N297A", "igg1")

    def test_structure_and_row_count(self):
        report = format_csv_report("igg1", "WT", "L234A/L235A/N297A", IGG1, self.mut)
        rows = report.strip().split("\n")
        # 메타 행 + 컬럼 헤더 + 서열 길이만큼 데이터 행
        assert len(rows) == 2 + len(IGG1)
        assert rows[1] == "eu_position,wt_aa,mut_aa,mutated"

    def test_mutated_rows_flagged(self):
        report = format_csv_report("igg1", "WT", "", IGG1, self.mut)
        mutated = [r for r in report.strip().split("\n") if r.endswith(",yes")]
        assert sorted(mutated) == sorted([
            "234,L,A,yes",
            "235,L,A,yes",
            "297,N,A,yes",
        ])

    def test_no_change_report(self):
        report = format_csv_report("igg1", "WT", "", IGG1, IGG1)
        assert ",yes" not in report
        assert report.count(",no") == len(IGG1)

    def test_unknown_isotype_leaves_position_empty(self):
        report = format_csv_report("unknown_iso", "WT", "", "AK", "AK")
        assert "eu_position,wt_aa,mut_aa,mutated\n,A,A,no\n,K,K,no\n" in report

    def test_invalid_inputs_return_empty(self):
        assert format_csv_report("igg1", "WT", "", "", self.mut) == ""
        assert format_csv_report("igg1", "WT", "", IGG1, None) == ""
        assert format_csv_report("igg1", "WT", "", "A" * 20001, "A") == ""
