from types import SimpleNamespace
from unittest.mock import AsyncMock, MagicMock

import flet as ft
import pytest

from fc_engineer.app import AntibodyStudioController, build_studio_ui
from fc_engineer.designer import build_schematic_scene


class FakeWindow:
    def __init__(self) -> None:
        self.prevent_close = False
        self.on_event = None
        self.destroy = AsyncMock()


class FakePage:
    """Lightweight headless stub matching the ft.Page surface used by build_studio_ui."""

    def __init__(self) -> None:
        self.title = ""
        self.padding = 0
        self.theme_mode = None
        self.window = FakeWindow()
        self.services = []
        self.controls = []
        self.dialogs = []
        self.update_calls = 0
        self._controls_by_id = {}

    def add(self, *controls) -> None:
        self.controls.extend(controls)

    def update(self) -> None:
        self.update_calls += 1

    def show_dialog(self, dlg) -> None:
        self.dialogs.append(dlg)

    def register_control(self, cid: int, ctrl) -> None:
        self._controls_by_id[cid] = ctrl

    def get_control(self, cid: int):
        return self._controls_by_id.get(cid)


def _walk_controls(root):
    if root is None:
        return
    yield root
    for attr in ("controls",):
        children = getattr(root, attr, None)
        if isinstance(children, list):
            for child in children:
                yield from _walk_controls(child)
    content = getattr(root, "content", None)
    if content is not None and not isinstance(content, str):
        yield from _walk_controls(content)


class TestAntibodyStudioController:
    def test_initial_defaults_and_format_switching(self):
        ctrl = AntibodyStudioController()
        assert ctrl.project.format_type == "igg"
        assert ctrl.selected_chain_id == "LC_1"
        assert ctrl.selected_domain_id == "lc1_vl"
        assert ctrl.apply_to_all_heavy is False
        assert ctrl.copied_fasta == ""

        ctrl.rename_project("My_BiSpecific_v1")
        assert ctrl.project.name == "My_BiSpecific_v1"
        with pytest.raises(ValueError):
            ctrl.rename_project("Bad\nTitle")
        with pytest.raises(ValueError):
            ctrl.rename_project("   ")

        ctrl.select_format("crossmab")
        assert ctrl.project.format_type == "crossmab"
        assert len(ctrl.project.chains) == 4

    def test_chain_and_domain_manipulation_and_pairing_cleanup(self):
        ctrl = AntibodyStudioController("igg")
        c_new = ctrl.add_chain()
        assert c_new.id == "Chain_1"
        assert ctrl.selected_chain_id == "Chain_1"
        assert ctrl.selected_domain_id is None  # IR3-M5 reset

        d1 = ctrl.add_domain_to_chain("Chain_1", "VH")
        d2 = ctrl.add_domain_to_chain("Chain_1", "CH1")
        assert [d.id for d in c_new.domains] == [d1.id, d2.id]
        ctrl.reorder_domain("Chain_1", 1, 0)
        assert [d.id for d in c_new.domains] == [d2.id, d1.id]

        # Move domain across chains prunes orphaned pairings (IR3-M2)
        pairings_before = len(ctrl.project.pairings)
        assert any(p.chain_a_id == "HC_1" and p.domain_a_id == "hc1_vh" for p in ctrl.project.pairings)
        ctrl.move_domain_across_chains("HC_1", "hc1_vh", "Chain_1", target_index=0)
        assert len(ctrl.project.pairings) < pairings_before
        assert not any(
            (p.chain_a_id == "HC_1" and p.domain_a_id == "hc1_vh")
            or (p.domain_a_id == "hc1_vh" or p.domain_b_id == "hc1_vh")
            for p in ctrl.project.pairings
        )

        # Same-chain reorder works even when chain already has 32 domains (IR3-M2)
        full_chain = ctrl.add_chain("FullChain")
        for _ in range(32):
            ctrl.add_domain_to_chain(full_chain.id, "G4S")
        assert len(full_chain.domains) == 32
        first_id = full_chain.domains[0].id
        ctrl.move_domain_across_chains(full_chain.id, first_id, full_chain.id, target_index=None)
        assert full_chain.domains[-1].id == first_id

        # Canvas manual drag pins coords; bool/NaN rejected
        ctrl.move_domain("Chain_1", d1.id, 300.0, 250.0)
        assert (d1.x, d1.y, d1.pinned) == (300.0, 250.0, True)
        with pytest.raises(ValueError):
            ctrl.move_domain("Chain_1", d1.id, True, 0.0)  # type: ignore[arg-type]

        ctrl.remove_chain("Chain_1")
        assert all(c.id != "Chain_1" for c in ctrl.project.chains)

    def test_update_domain_properties_isotype_reset_and_rollback(self):
        ctrl = AntibodyStudioController("igg")
        hc1_ch2 = "hc1_ch2"
        # Switching isotype without explicit allotype resets incompatible IgG1 allotype cleanly (IR3-M5)
        errs = ctrl.update_domain_properties("HC_1", hc1_ch2, isotype="igg4")
        assert not errs
        dom = ctrl._get_domain("HC_1", hc1_ch2)[1]
        assert dom.isotype == "igg4"
        assert dom.allotype == ""

        # Invalid custom sequence on CH2 rolls back cleanly
        errs_bad = ctrl.update_domain_properties("HC_1", hc1_ch2, custom_sequence="AAAA")
        assert errs_bad
        assert dom.custom_sequence == ""

        # Bad VH sequence rolls back without overwriting prior valid state
        vh_errs = ctrl.update_domain_properties("HC_1", "hc1_vh", custom_sequence="1234!")
        assert vh_errs
        assert ctrl._get_domain("HC_1", "hc1_vh")[1].custom_sequence == ""

    def test_apply_to_all_heavy_atomic_rollback_and_bite_fallback(self):
        ctrl = AntibodyStudioController("crossmab")
        # Make HC_Knob IgG3 (lacks H435), keep HC_Hole_Cross IgG1 -> FcRnKO must fail atomically across BOTH
        ctrl.update_domain_properties("HC_Knob", "hck_ch3", isotype="igg3")
        errs = ctrl.toggle_preset("HC_Hole_Cross", "FcRnKO", True, apply_to_all_heavy=True)
        assert errs
        assert "FcRnKO" not in ctrl._get_domain("HC_Hole_Cross", "hch_ch3")[1].applied_presets

        # Successful toggle applies to both heavy chains symmetrically
        assert not ctrl.toggle_preset("HC_Knob", "LALA", True, apply_to_all_heavy=True)
        assert "LALA" in ctrl._get_domain("HC_Knob", "hck_ch2")[1].applied_presets
        assert "LALA" in ctrl._get_domain("HC_Hole_Cross", "hch_ch2")[1].applied_presets

        # Fc-less project (BiTE) with apply_to_all_heavy=True must report error rather than silently succeed (IR3-M5)
        ctrl.select_format("bite")
        bite_errs = ctrl.toggle_preset("Chain_1", "LALA", True, apply_to_all_heavy=True)
        assert bite_errs

    def test_persistence_load_preserves_project_on_error(self, tmp_path):
        ctrl = AntibodyStudioController("scfv_fc")
        ab_file = tmp_path / "proj.ab"
        fasta_file = tmp_path / "proj.fasta"
        png_file = tmp_path / "proj.png"
        ctrl.save_ab(str(ab_file))
        ctrl.save_fasta(str(fasta_file))
        ctrl.save_png(str(png_file))
        assert ab_file.is_file() and fasta_file.is_file() and png_file.is_file()

        # Failed load must preserve current workspace project
        bad_ab = tmp_path / "corrupt.ab"
        bad_ab.write_text("NOT A VALID AB FILE", encoding="utf-8")
        with pytest.raises(ValueError):
            ctrl.load_ab(str(bad_ab))
        assert ctrl.project.format_type == "scfv_fc"

        # Valid load switches workspace
        ctrl.select_format("bite")
        assert ctrl.project.format_type == "bite"
        ctrl.load_ab(str(ab_file))
        assert ctrl.project.format_type == "scfv_fc"

    def test_clipboard_security_timer_race_and_shutdown(self, monkeypatch):
        fake_clip = {"val": ""}
        monkeypatch.setattr("fc_engineer.app.pyperclip.copy", lambda txt: fake_clip.__setitem__("val", txt))
        monkeypatch.setattr("fc_engineer.app.pyperclip.paste", lambda: fake_clip["val"])

        scheduled = []

        class StubTimer:
            def __init__(self, interval, fn, args=()):
                self.interval = interval
                self.fn = fn
                self.args = args
                self.daemon = False
                self.cancelled = False
                scheduled.append(self)

            def start(self):
                pass

            def cancel(self):
                self.cancelled = True

        monkeypatch.setattr("fc_engineer.app.threading.Timer", StubTimer)

        ctrl = AntibodyStudioController("igg")
        assert ctrl.copy_text_secure("") is False
        assert ctrl.copy_text_secure("SEQ_A") is True
        assert fake_clip["val"] == "SEQ_A"
        first_timer = scheduled[-1]
        assert first_timer.daemon is True

        # Copy SEQ_B -> cancels timer 1, wipes SEQ_A first, writes SEQ_B
        assert ctrl.copy_text_secure("SEQ_B") is True
        assert first_timer.cancelled is True
        assert fake_clip["val"] == "SEQ_B"
        second_timer = scheduled[-1]

        # Even if stale timer 1 callback runs out-of-order, SEQ_B must NOT be wiped
        first_timer.fn(*first_timer.args)
        assert fake_clip["val"] == "SEQ_B"
        assert ctrl.copied_fasta == "SEQ_B"

        # Active timer 2 fires -> wipes SEQ_B
        second_timer.fn(*second_timer.args)
        assert fake_clip["val"] == ""
        assert ctrl.copied_fasta == ""

        # Pyperclip exception returns False and preserves empty state (IR3-M5)
        monkeypatch.setattr(
            "fc_engineer.app.pyperclip.copy",
            MagicMock(side_effect=RuntimeError("no clipboard")),
        )
        assert ctrl.copy_text_secure("SEQ_C") is False
        assert ctrl.copied_fasta == ""

        # Shutdown cancels active timer and clears matching clipboard
        monkeypatch.setattr("fc_engineer.app.pyperclip.copy", lambda txt: fake_clip.__setitem__("val", txt))
        ctrl.copy_text_secure("SEQ_D")
        ctrl.shutdown()
        assert fake_clip["val"] == ""
        assert ctrl.copied_fasta == ""


@pytest.mark.asyncio
async def test_build_studio_ui_headless_dnd_pan_and_window_close():
    fake_page = FakePage()
    ctrl = AntibodyStudioController("igg")
    returned = build_studio_ui(fake_page, ctrl)
    assert returned is ctrl
    assert fake_page.title == "Fc Engineering Studio"
    assert fake_page.window.prevent_close is True
    assert len(fake_page.services) == 1
    assert isinstance(fake_page.services[0], ft.FilePicker)
    assert len(fake_page.controls) == 1

    # Locate DragTargets & GestureDetectors inside rendered tree
    all_controls = list(_walk_controls(fake_page.controls[0]))
    drag_targets = [c for c in all_controls if isinstance(c, ft.DragTarget)]
    gestures = [c for c in all_controls if isinstance(c, ft.GestureDetector)]
    draggables = [c for c in all_controls if isinstance(c, ft.Draggable)]
    assert len(drag_targets) >= 4
    assert len(gestures) >= 14  # 14 domains in Monospecific IgG
    assert any(d.data == "palette:VH" for d in draggables)

    # Guarded DnD drop: e.src_id = None must safely no-op (IR3-M6)
    hc1_domains_before = len(ctrl._get_chain("HC_1").domains)
    drag_targets[0].on_accept(SimpleNamespace(src_id=None))
    assert len(ctrl._get_chain("HC_1").domains) == hc1_domains_before

    # Valid DnD drop from palette onto LC_1 first slot inserts domain & refreshes
    vh_draggable = next(d for d in draggables if d.data == "palette:G4S")
    fake_page.register_control(999, vh_draggable)
    drag_targets[0].on_accept(SimpleNamespace(src_id=999))
    assert len(ctrl._get_chain("LC_1").domains) == 3

    # Re-collect rebuilt canvas GestureDetectors after render_workspace()
    rebuilt_gestures = [
        c for c in _walk_controls(fake_page.controls[0]) if isinstance(c, ft.GestureDetector)
    ]
    assert len(rebuilt_gestures) >= 15

    # Guarded pan update: local_delta=None is a safe no-op; valid delta pins & moves node
    rebuilt_gestures[0].on_pan_update(SimpleNamespace(local_delta=None))
    orig_x = ctrl.project.chains[0].domains[0].x
    rebuilt_gestures[0].on_pan_update(
        SimpleNamespace(local_delta=SimpleNamespace(x=25.0, y=10.0))
    )
    assert ctrl.project.chains[0].domains[0].x == orig_x + 25.0
    assert ctrl.project.chains[0].domains[0].pinned is True

    # Non-CLOSE window event does NOT shutdown or destroy window (IR3-N6)
    ctrl.copied_fasta = "SENSITIVE"
    await fake_page.window.on_event(SimpleNamespace(type=ft.WindowEventType.RESIZE))
    fake_page.window.destroy.assert_not_awaited()
    assert ctrl.copied_fasta == "SENSITIVE"

    # CLOSE window event triggers shutdown and awaits window.destroy()
    await fake_page.window.on_event(SimpleNamespace(type=ft.WindowEventType.CLOSE))
    fake_page.window.destroy.assert_awaited_once()
    assert ctrl.copied_fasta == ""
