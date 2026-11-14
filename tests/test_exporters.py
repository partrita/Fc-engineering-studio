import io
import json
import os
import pytest
from PIL import Image

from fc_engineer.config import SEQUENCES
from fc_engineer.core import apply_mutations
from fc_engineer.designer import (
    AntibodyProject,
    DomainBlock,
    InterChainPairing,
    PolypeptideChain,
    add_file_custom_mutation,
    assemble_chain_sequence,
    build_format_template,
    build_preset_map,
    route_mutations_to_chain,
)
from fc_engineer.exporters import (
    MAX_PNG_DIMENSION,
    export_schematic_png,
    format_csv_report,
    format_genbank,
    format_project_fasta,
    load_ab_file,
    parse_ab_project,
    sanitize_locus,
    save_ab_file,
    save_fasta_file,
    save_schematic_png,
    serialize_ab_project,
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


class TestProjectFastaAndAbPersistence:
    def test_format_project_fasta_and_save_fasta(self, tmp_path):
        proj = build_format_template("crossmab")
        proj.name = "CrossMab\nHeaderInjection"
        fasta, errs = format_project_fasta(proj)
        assert not errs
        assert fasta.startswith(">CrossMab HeaderInjection_")
        assert ">CrossMab HeaderInjection_HC_Knob" in fasta

        out = tmp_path / "crossmab"
        save_fasta_file(str(out), proj)
        saved_path = tmp_path / "crossmab.fasta"
        assert saved_path.is_file()
        assert saved_path.read_text(encoding="utf-8") == fasta

        # Uppercase .FA accepted case-insensitively; bad suffix rejected
        out_fa = tmp_path / "UPPER.FA"
        save_fasta_file(str(out_fa), proj)
        assert out_fa.is_file()
        with pytest.raises(ValueError):
            save_fasta_file(str(tmp_path / "bad.txt"), proj)

    def test_ab_round_trip_with_coords_history_and_file_custom_mutations(self, tmp_path):
        proj = build_format_template("crossmab")
        proj.name = "AntiHer2xCD3_CrossMab"
        assert not add_file_custom_mutation(proj, "MyHalfLife", "M252Y/N434S")
        hck = next(c for c in proj.chains if c.id == "HC_Knob")
        assert not route_mutations_to_chain(
            hck, "M252Y/N434S", preset_label="MyHalfLife", preset_map=build_preset_map(proj)
        )
        # Custom sequence & pinned coordinate preservation
        hck.domains[0].custom_sequence = "EVQLVESGGGLV"
        hck.domains[0].x = 123.5
        hck.domains[0].y = 456.25
        hck.domains[0].pinned = True

        # Generate > 50 history entries to verify serialize/parse round-trip parity (IR3-M1)
        for i in range(55):
            assert not route_mutations_to_chain(
                hck, "N297A" if i % 2 == 0 else "D265A", preset_map=build_preset_map(proj)
            )
        assert len(hck.domains[3].mutation_history) == 50

        ab_path = tmp_path / "workspace"
        save_ab_file(str(ab_path), proj)
        resolved_ab = tmp_path / "workspace.ab"
        assert resolved_ab.is_file()

        raw_text = resolved_ab.read_text(encoding="utf-8")
        assert '/mutations="T366W/N434S"' in raw_text
        assert '/note="EU_366:T->W"' in raw_text

        reloaded = load_ab_file(str(resolved_ab))
        assert reloaded.name == proj.name
        assert reloaded.format_type == "crossmab"
        assert reloaded.file_custom_mutations == proj.file_custom_mutations
        assert len(reloaded.pairings) == len(proj.pairings)
        r_hck = next(c for c in reloaded.chains if c.id == "HC_Knob")
        assert r_hck.domains[0].custom_sequence == "EVQLVESGGGLV"
        assert (r_hck.domains[0].x, r_hck.domains[0].y, r_hck.domains[0].pinned) == (
            123.5,
            456.25,
            True,
        )
        assert len(r_hck.domains[3].mutation_history) == 50

        orig_fasta, _ = format_project_fasta(proj)
        reloaded_fasta, _ = format_project_fasta(reloaded)
        assert reloaded_fasta == orig_fasta

    def test_ab_file_security_and_path_sanitization(self, tmp_path):
        proj = build_format_template("igg")
        ab_text = serialize_ab_project(proj)

        # Bad suffix / missing file / non-regular file / oversized
        with pytest.raises(ValueError):
            load_ab_file(str(tmp_path / "wrong.gb"))

        missing_path = tmp_path / "secret_dir_name" / "missing.ab"
        with pytest.raises(FileNotFoundError) as exc_info:
            load_ab_file(str(missing_path))
        assert "secret_dir_name" not in str(exc_info.value)

        subdir = tmp_path / "folder.ab"
        subdir.mkdir()
        with pytest.raises(ValueError) as dir_exc:
            load_ab_file(str(subdir))
        assert str(tmp_path) not in str(dir_exc.value)

        big_ab = tmp_path / "big.ab"
        big_ab.write_text(ab_text + (" " * (1024 * 1024 + 10)), encoding="utf-8")
        with pytest.raises(ValueError):
            load_ab_file(str(big_ab))

        if os.path.exists("/dev/zero"):
            # Wait: symlink to /dev/zero with .ab extension tests special file rejection
            zero_link = tmp_path / "zero.ab"
            zero_link.symlink_to("/dev/zero")
            with pytest.raises(ValueError):
                load_ab_file(str(zero_link))

        bad_utf8 = tmp_path / "bad_bytes.ab"
        bad_utf8.write_bytes(b"\xff\xfe\x00\x00")
        with pytest.raises(ValueError) as utf_exc:
            load_ab_file(str(bad_utf8))
        assert str(tmp_path) not in str(utf_exc.value)

    def test_ab_schema_and_integrity_rejections(self):
        proj = build_format_template("igg")
        valid_text = serialize_ab_project(proj)

        # 1. Trailing non-whitespace content after final // (IR3-N5)
        with pytest.raises(ValueError):
            parse_ab_project(valid_text + "SPOOFED TRAILER\n")

        # 2. Tampered ORIGIN sequence mismatch
        tampered_origin = valid_text.replace("astkgpsvf", "xxxxxxxxx", 1).replace(
            "ASTKGPSVF", "XXXXXXXXX", 1
        )
        # Ensure we mutate ORIGIN chunk reliably
        origin_pos = valid_text.rfind("ORIGIN")
        tampered_origin = valid_text[:origin_pos] + valid_text[origin_pos:].replace("A", "G", 1)
        with pytest.raises(ValueError):
            parse_ab_project(tampered_origin)

        # Helper to tamper manifest JSON cleanly
        def tamper_manifest(mutator_fn):
            begin_tag = "FC_STUDIO_PROJECT_V1_BEGIN"
            end_tag = "FC_STUDIO_PROJECT_V1_END"
            lines = valid_text.splitlines()
            b = next(i for i, l in enumerate(lines) if l.strip() == begin_tag)
            e = next(i for i, l in enumerate(lines) if l.strip() == end_tag)
            manifest = json.loads("\n".join(l[12:] for l in lines[b + 1 : e]))
            mutator_fn(manifest)
            new_json_lines = ["            " + jl for jl in json.dumps(manifest).splitlines()]
            return "\n".join(lines[: b + 1] + new_json_lines + lines[e:]) + "\n"

        # 3. Invalid schema / NaN / bool coords / bounds / shadowing / duplicates
        with pytest.raises(ValueError):
            parse_ab_project(tamper_manifest(lambda m: m.__setitem__("schema_version", True)))
        with pytest.raises(ValueError):
            parse_ab_project(tamper_manifest(lambda m: m.__setitem__("format_type", "bad")))
        with pytest.raises(ValueError):
            parse_ab_project(tamper_manifest(lambda m: m["chains"][0]["domains"][0].__setitem__("x", True)))
        with pytest.raises(ValueError):
            parse_ab_project(
                tamper_manifest(lambda m: m["chains"][0]["domains"][0].__setitem__("x", 99999.0))
            )
        with pytest.raises(ValueError):
            parse_ab_project(
                tamper_manifest(lambda m: m["chains"][0].__setitem__("domains", []))
            )
        with pytest.raises(ValueError):
            parse_ab_project(
                tamper_manifest(
                    lambda m: m.__setitem__(
                        "file_custom_mutations", [{"label": "LALA", "value": "L234A"}]
                    )
                )
            )
        with pytest.raises(ValueError):
            parse_ab_project(
                tamper_manifest(
                    lambda m: m["chains"][0]["domains"][0].__setitem__(
                        "custom_mutations", ["INVALID!"]
                    )
                )
            )
        with pytest.raises(ValueError):
            parse_ab_project(
                valid_text.replace('"x": 40.0', '"x": NaN', 1)
            )

        # 4. Control characters in project name blocked before save (IR3-M6)
        proj_bad_name = build_format_template("bite")
        proj_bad_name.name = "Bad\nName"
        with pytest.raises(ValueError):
            serialize_ab_project(proj_bad_name)


class TestSchematicPngExporter:
    def test_export_transparent_rgba_png_and_save(self, tmp_path):
        proj = build_format_template("igg_scfv")
        png_bytes = export_schematic_png(proj, scale=2.0)
        assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n")

        with Image.open(io.BytesIO(png_bytes)) as im:
            assert im.mode == "RGBA"
            assert im.width > 100 and im.height > 100
            # Top-left corner outside padded domain bounds must be fully transparent
            assert im.getpixel((0, 0)) == (0, 0, 0, 0)
            # Center of first domain block contains rendered opaque pixels
            pixels = im.get_flattened_data() if hasattr(im, "get_flattened_data") else im.getdata()
            non_transparent = any(px[3] > 0 for px in pixels)
            assert non_transparent

        png_path = tmp_path / "schematic"
        save_schematic_png(str(png_path), proj, scale=1.5)
        assert (tmp_path / "schematic.png").is_file()

    def test_export_empty_project_png(self):
        empty_proj = AntibodyProject(name="Empty", format_type="custom", chains=[], pairings=[])
        png_bytes = export_schematic_png(empty_proj, scale=1.0)
        assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n")
        with Image.open(io.BytesIO(png_bytes)) as im:
            assert im.mode == "RGBA"
            assert im.size == (320, 200)

    def test_export_png_rejects_invalid_scale_and_oversized_dimensions(self):
        proj = build_format_template("igg")
        for bad_scale in (0.1, 5.0, float("nan"), float("inf"), True, "2"):
            with pytest.raises(ValueError):
                export_schematic_png(proj, scale=bad_scale)  # type: ignore[arg-type]

        proj.chains[0].domains[0].x = float(MAX_PNG_DIMENSION)
        with pytest.raises(ValueError):
            export_schematic_png(proj, scale=2.0)

