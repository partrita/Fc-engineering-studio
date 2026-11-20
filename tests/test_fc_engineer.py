"""Tests for Fc engineering catalog, controller workflows, and Fc Engineer GUI tab."""

from types import SimpleNamespace

import flet as ft
import pytest

from fc_engineer import config
from fc_engineer.app import AntibodyStudioController, build_studio_ui
from fc_engineer.designer import build_format_template, build_preset_map
from fc_engineer.fc_catalog import (
    FC_CATEGORIES,
    FC_RECIPES,
    build_eu_residue_rows,
    get_hotspot_info,
    get_preset_metadata,
    summarize_fc_profile,
)


class _FakeWindow:
    def __init__(self) -> None:
        self.prevent_close = False
        self.on_event = None


class _FakePage:
    def __init__(self) -> None:
        self.title = ""
        self.padding = 0
        self.theme_mode = None
        self.window = _FakeWindow()
        self.services = []
        self.controls = []
        self.dialogs = []
        self.update_calls = 0

    def add(self, *controls) -> None:
        self.controls.extend(controls)

    def update(self) -> None:
        self.update_calls += 1

    def show_dialog(self, dlg) -> None:
        self.dialogs.append(dlg)


def _walk_controls(root):
    if root is None:
        return
    yield root
    for attr in ("controls", "tabs", "items"):
        children = getattr(root, attr, None)
        if isinstance(children, list):
            for child in children:
                yield from _walk_controls(child)
    content = getattr(root, "content", None)
    if content is not None and not isinstance(content, str):
        yield from _walk_controls(content)


class TestFcCatalog:
    def test_all_builtin_presets_have_valid_metadata(self):
        for item in config.COMMON_MUTATIONS:
            meta = get_preset_metadata(item["label"], item["value"])
            assert meta.label == item["label"]
            assert meta.category in FC_CATEGORIES
            assert meta.summary
            assert meta.adcc in ("up", "down", "neutral")
            assert meta.cdc in ("up", "down", "neutral")
            assert meta.fcrn in ("up", "down", "neutral")

    def test_project_custom_preset_metadata(self):
        meta = get_preset_metadata("MyCustomFc", "M252Y/N434S", is_project_custom=True)
        assert meta.category == "Project"
        assert "M252Y/N434S" in meta.summary

    def test_all_fc_recipes_apply_cleanly_to_igg(self):
        for rec in FC_RECIPES:
            ctrl = AntibodyStudioController("igg")
            errs = ctrl.apply_fc_recipe(
                "HC_1",
                presets=rec["presets"],
                isotype=rec.get("isotype"),
                reset_first=True,
                apply_to_all_heavy=True,
            )
            assert not errs, f"Recipe {rec['id']} failed: {errs}"
            hc1 = ctrl._get_chain("HC_1")
            hc2 = ctrl._get_chain("HC_2")
            for p_lbl in rec["presets"]:
                assert any(p_lbl in d.applied_presets for d in hc1.domains)
                assert any(p_lbl in d.applied_presets for d in hc2.domains)

    def test_summarize_fc_profile_effects_and_deltas(self):
        ctrl = AntibodyStudioController("igg")
        pmap = ctrl.get_preset_map()
        wt_prof = summarize_fc_profile(ctrl._get_chain("HC_1"), pmap)
        assert wt_prof["has_fc"] is True
        assert wt_prof["isotype"] == "igg1"
        assert wt_prof["adcc"] == "neutral"
        assert wt_prof["cdc"] == "neutral"
        assert wt_prof["fcrn"] == "neutral"
        assert wt_prof["diffs"] == []

        assert not ctrl.toggle_preset("HC_1", "PGLALA", True)
        assert not ctrl.toggle_preset("HC_1", "LS", True)
        mut_prof = summarize_fc_profile(ctrl._get_chain("HC_1"), pmap)
        assert mut_prof["adcc"] == "down"
        assert mut_prof["cdc"] == "down"
        assert mut_prof["fcrn"] == "up"
        assert len(mut_prof["diffs"]) == 5

    def test_build_eu_residue_rows_and_hotspots(self):
        proj = build_format_template("igg", isotype="igg1")
        pmap = build_preset_map(proj)
        hc1 = next(c for c in proj.chains if c.id == "HC_1")
        ch2_rows = build_eu_residue_rows(hc1, pmap, domain_types=("CH2",))
        eu_positions = [r["eu_pos"] for r in ch2_rows]
        assert 234 in eu_positions
        assert 297 in eu_positions
        assert 329 in eu_positions
        row_297 = next(r for r in ch2_rows if r["eu_pos"] == 297)
        assert row_297["wt_aa"] == "N"
        assert row_297["hotspot_tag"] == "N-Glyco"
        assert any("N-glycosylation" in t for t in row_297["ptm_tags"])
        assert get_hotspot_info(366)[0] == "CH3 Dimer"


class TestFcControllerWorkflows:
    def test_set_chain_isotype_atomic_and_incompatible_clear(self):
        ctrl = AntibodyStudioController("igg")
        assert not ctrl.toggle_preset("HC_1", "LALA", True, apply_to_all_heavy=True)

        # Switching to IgG4 fails cleanly when LALA conflicts (EU 234 is F in IgG4)
        errs = ctrl.set_chain_isotype(
            "HC_1", "igg4", apply_to_all_heavy=True, clear_incompatible_mutations=False
        )
        assert errs
        assert ctrl._get_domain("HC_1", "hc1_ch2")[1].isotype == "igg1"

        # With clear_incompatible_mutations=True, switches all constant domains on both heavy chains
        errs2 = ctrl.set_chain_isotype(
            "HC_1", "igg4", apply_to_all_heavy=True, clear_incompatible_mutations=True
        )
        assert not errs2
        for cid in ("HC_1", "HC_2"):
            for d in ctrl._get_chain(cid).domains:
                if d.domain_type in ("CH1", "Hinge", "CH2", "CH3"):
                    assert d.isotype == "igg4"

        # Single-step Undo restores IgG1 + LALA
        assert ctrl.undo()
        assert ctrl._get_domain("HC_1", "hc1_ch2")[1].isotype == "igg1"
        assert "LALA" in ctrl._get_domain("HC_1", "hc1_ch2")[1].applied_presets

    def test_mutate_residue_at_eu_and_unbundle_preset(self):
        ctrl = AntibodyStudioController("igg")
        assert not ctrl.toggle_preset("HC_1", "LALA", True)
        # Surgical edit of EU 234 from A -> F must unbundle LALA and preserve L235A
        errs = ctrl.mutate_residue_at_eu("HC_1", "hc1_ch2", 234, "F")
        assert not errs
        ch2 = ctrl._get_domain("HC_1", "hc1_ch2")[1]
        assert "LALA" not in ch2.applied_presets
        assert "L234F" in ch2.custom_mutations
        assert "L235A" in ch2.custom_mutations

        # Surgical revert of EU 234 to WT (L) leaves only L235A
        errs2 = ctrl.mutate_residue_at_eu("HC_1", "hc1_ch2", 234, "L")
        assert not errs2
        assert "L234F" not in ch2.custom_mutations
        assert "L235A" in ch2.custom_mutations

        # Invalid amino acid rejected cleanly
        assert ctrl.mutate_residue_at_eu("HC_1", "hc1_ch2", 234, "1")


class TestFcEngineerGuiTab:
    def test_fc_engineer_tab_recipes_and_eu_point_mutation_picker(self):
        fake_page = _FakePage()
        ctrl = AntibodyStudioController("igg")
        build_studio_ui(fake_page, ctrl)

        all_controls = list(_walk_controls(fake_page.controls[0]))
        tabs = [c for c in all_controls if isinstance(c, ft.Tab)]
        assert [t.label for t in tabs[:4]] == ["Studio", "Project", "Style", "Fc Engineer"]

        # Click "Silent + Long t½ (LALA+LS)" recipe button
        recipe_btn = next(
            c
            for c in all_controls
            if isinstance(c, ft.OutlinedButton)
            and getattr(c, "content", None) == "Silent + Long t½ (LALA+LS)"
        )
        recipe_btn.on_click(SimpleNamespace())
        assert "LALA" in ctrl._get_domain("HC_1", "hc1_ch2")[1].applied_presets
        assert "LS" in ctrl._get_domain("HC_1", "hc1_ch3")[1].applied_presets

        # Click EU 297 tile in the Interactive EU Browser to open Point Mutation Picker
        rebuilt = list(_walk_controls(fake_page.controls[0]))
        tile_297 = next(
            c
            for c in rebuilt
            if isinstance(c, ft.Container)
            and isinstance(c.tooltip, str)
            and c.tooltip.startswith("EU 297 (CH2):")
        )
        tile_297.on_click(SimpleNamespace())
        assert fake_page.dialogs
        picker_dlg = fake_page.dialogs[-1]
        assert isinstance(picker_dlg, ft.AlertDialog)
        assert "EU 297" in picker_dlg.title.value

        # Click 'A' button in the picker to apply N297A in 1 click
        picker_controls = list(_walk_controls(picker_dlg.content))
        aa_a_tile = next(
            c
            for c in picker_controls
            if isinstance(c, ft.Container)
            and isinstance(c.content, ft.Text)
            and c.content.value == "A"
        )
        aa_a_tile.on_click(SimpleNamespace())
        assert "N297A" in ctrl._get_domain("HC_1", "hc1_ch2")[1].custom_mutations

        # Fit zoom button in Studio toolbar updates zoom_scale cleanly
        fit_btn = next(
            c
            for c in rebuilt
            if isinstance(c, ft.OutlinedButton) and getattr(c, "content", None) == "Fit"
        )
        fit_btn.on_click(SimpleNamespace())
        assert ctrl.zoom_scale == 0.75
