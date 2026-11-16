# 3-Tab Studio Consolidation & Side-Panel Inspector Integration

* **Author:** fkt
* **Date:** 2026-11-16

## Overview
Fc Engineering Studio currently separates domain inspection/editing, Fc mutations, pairings, and `Domain Palette` cards into an extra `Inspector & Domain` tab while `Studio` already displays a read-only `Domain Info` side panel. This forces users to switch tabs repeatedly to edit a domain they just clicked on the 2D schematic. This feature removes the redundant `Inspector & Domain` tab and `Domain Palette` cards, simplifying the app to **3 tabs (`Studio`, `Project`, `Theme & Settings`)**, and consolidates domain editing, Fc mutations, developability, and `Pairings` directly into `Studio`'s right-hand side panel.

## Analysis

### Key Contradictions & Resolutions (7-Lens Goldfish Synthesis)
1. **Rich Side-Panel Editing vs. Minimal `Studio` Visual Clarity:**
   * **Problem:** Moving `Isotype`/`Allotype`, custom sequence editing, Fc mutation presets, ad-hoc EU mutations, file-scoped custom mutations, developability report, and `Pairings` into `Studio`'s right-hand side panel can make the side panel long and visually noisy if everything is expanded simultaneously.
   * **Resolution:** Organize `Studio`'s right-hand side panel (`width=340`) with a compact two-section or segmented toggle (`Domain` | `Pairings`, or clean `ft.ExpansionTile` / compact sub-tabs inside the side panel) so `Domain` editing (`Isotype`/`Allotype`/`Sequence`/`Fc Mutations`/`Developability`) stays immediately visible on click while `Pairings` (`VH_VL`, `CH1_CL`, `HINGE_DISULFIDE`, `CH3_DIMER`, `DISULFIDE_BOND`) sits cleanly one click away inside the same side panel without crowding the 2D canvas.
2. **Removing `Domain Palette` Cards vs. Chain Strip Reordering & Headless UI Tests (`tests/test_ui.py`):**
   * **Problem:** `tests/test_ui.py::test_build_studio_ui_headless_dnd_pan_and_window_close` currently asserts `any(d.data == "palette:VH" for d in draggables)` and simulates dropping `palette:G4S` onto a `DragTarget`, and `test_light_dark_mode_switching` checks `palette:` draggables.
   * **Resolution:** Unify domain addition onto each chain row's `+` quick-add `PopupMenuButton` (while keeping `Draggable` `chip:<chain_id>:<domain_id>` + `DragTarget` slots on the `N -> C` chain strip for drag-reordering and cross-chain moves). Update `tests/test_ui.py` to test `+` quick-add menu items and `chip:` drag-and-drop instead of `palette:` cards, and verify Light/Dark contrast on chain strip chips.

## Goals vs. Non-Goals

### Goals
* Simplify the top-level tabs from 4 to **3 (`Studio`, `Project`, `Theme & Settings`)**.
* Integrate selected domain editing (`Isotype`/`Allotype`/`CL Subtype`/custom sequence/`Delete Domain`), Fc mutations, developability summary, and `Pairings` into `Studio`'s right-hand side panel (`studio_side_panel`).
* Remove `Domain Palette` cards (`palette_row` / `rebuild_palette`) and unify domain addition via each chain row's `+` `PopupMenuButton`.
* Update `tests/test_ui.py` and `README.md` to match the 3-tab structure.

### Non-Goals
* Changing `core.py`, `designer.py`, `config.py`, or `exporters.py` domain/serialization logic.
* Removing drag-reordering (`chip:<chain_id>:<domain_id>`) from the `N -> C` chain strip.

## Alternatives Weighed

1. **Side-Panel Organization: Compact `Domain | Pairings` Sub-Tabs vs. Long Single Scrollable Column:**
   * **Chosen:** Compact `Domain | Pairings` sub-tabs (or `ft.Tabs` inside `studio_side_panel`) keeps `Domain` editing concise and `Pairings` out of the way until needed.
   * *Rejected by author: Stacking all inspector controls and pairings in one long vertical column creates excessive scrolling.*

## Key Decisions & Staging for Next Steps

| Decision Area | Agreed Specification | Target Files |
| :--- | :--- | :--- |
| **3-Tab Structure** | `ft.Tabs(length=3)`: `Studio`, `Project`, `Theme & Settings` (removing `Inspector & Domain`). | `src/fc_engineer/app.py` |
| **Unified `Studio` Side Panel** | `studio_side_panel` (`width=340`) houses `Domain` (`inspector_col`: properties, Fc presets, EU diffs, developability, sequence) and `Pairings` (`pairing_controls_col`). | `src/fc_engineer/app.py` |
| **Domain Addition** | Remove `palette_row` / `rebuild_palette`; use each chain row's `+` `PopupMenuButton` for adding domains and `chip:` `Draggable`/`DragTarget` for reordering. | `src/fc_engineer/app.py`, `tests/test_ui.py` |
