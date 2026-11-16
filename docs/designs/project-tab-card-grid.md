# Project Tab 2-Column Card Grid & Direct File Action Forms

* **Author:** fkt
* **Date:** 2026-11-16

## Overview
The `Project` tab in Fc Engineering Studio currently places a `File` `PopupMenuButton` dropdown next to fixed-width `Project Name` and `Format Template` inputs above a loose vertical stack of multispecific controls. Users must open the dropdown menu to reach `.ab`, `FASTA`, `PNG`, and `.yaml` actions, and fixed pixel widths (`220`, `185`, `190`, `290`) leave ragged right edges. This feature replaces the `File` dropdown with direct-action button cards inside a balanced 2-column card grid with uniform `expand=True` alignment and full Light/Dark theme support.

## Analysis

### Key Contradictions & Resolutions (7-Lens Goldfish Synthesis)
1. **Column Alignment (`STRETCH`) vs. Flet Scroll Safety:**
   * **Problem:** Using `CrossAxisAlignment.STRETCH` on an `ft.Row` inside a scrollable `ft.Column` can trigger unbounded-height layout errors in Flet desktop.
   * **Resolution:** Use `CrossAxisAlignment.START` with `spacing=16` on the 2-column `ft.Row` (`left_col` and `right_col`, each `expand=True`) inside a scrollable `ft.Column`, and match card padding (`16`), spacing (`12`), and `expand=True` inner rows for clean alignment.
2. **Hardcoded Control Widths vs. Flush 2-Column Alignment:**
   * **Problem:** `project_name_tf` (`width=220`), `format_dd` (`width=185`), and `dd_ca`/`dd_cb`/`dd_mode` (`width=190`/`290`) wrap unevenly inside a half-width column.
   * **Resolution:** Replace fixed widths inside `Project` tab cards with `expand=True` inside `ft.Row`s (`Project Name` + `Format Template` side-by-side, `Chain A` + `Chain B` side-by-side, and full-width `Heterodimer Mode`).
3. **Box-Within-a-Box Nesting vs. Surface Hierarchy:**
   * **Problem:** Putting bordered cards inside a bordered `project_panel` creates nested borders unlike the `Studio` tab.
   * **Resolution:** Use `pal["bg"]` (borderless) for `project_panel`, and build the 5 cards via a shared `_make_project_card()` factory (`bgcolor=pal["panel"]`, `border=ft.Border.all(1, pal["border"])`, `border_radius=8`, `padding=16`) tracked in `project_cards` so `on_theme_mode_toggle()` updates every card uniformly.

## Goals vs. Non-Goals

### Goals
* Remove `file_menu` (`ft.PopupMenuButton`) and its stale `file_menu.content.border` update in `on_theme_mode_toggle()`.
* Reorganize `project_panel` into a balanced 2-column card grid:
  * **Left Column:** `Project Setup` card (`Project Name`, `Format Template`, `Add Chain`, `Clone Selected Chain`) + `Multispecific & Heterodimer Builder` card (`Chain A`, `Chain B`, `Heterodimer Mode`, `Apply Heterodimer Fc`, `Swap CH1 <-> CL (CrossMab)`).
  * **Right Column:** `Project Files (.ab)` card (`Load .ab`, `Save .ab`) + `Sequence & Schematic Export` card (`Export FASTA`, `Copy FASTA (30s Auto-Clear)`, `Export PNG`) + `Theme Files (.yaml)` card (`Load Theme .yaml`, `Save Theme .yaml`).
* Replace hardcoded widths in `Project` controls with `expand=True` rows and add tests in `tests/test_ui.py` verifying direct-action buttons and Light/Dark card updates.

### Non-Goals
* Modifying `AntibodyStudioController` logic, file exporters, or `Studio` / `Theme & Settings` layouts.
* Refactoring `wiz_state` persistence across `render_workspace()` calls.

## Alternatives Weighed

1. **Grid Container: `ft.Row` (2 `expand=True` Columns) vs. `ft.ResponsiveRow`:**
   * **Chosen:** `ft.Row` with two `expand=True` `ft.Column`s (`spacing=16`) matches `studio_panel` and avoids sudden desktop breakpoint jumps.
   * *Rejected by author: `ft.ResponsiveRow` introduces breakpoint wrapping jumps on desktop window resizes.*
2. **Theme Files Card: Keep in `Project` Right Column vs. `Theme & Settings` Only:**
   * **Chosen:** Keep `Theme Files (.yaml)` in the right column of `Project` to preserve 100% file I/O parity with the removed `File` menu.
   * *Rejected by author: Removing `.yaml` theme I/O from `Project` splits file operations across tabs.*

## Key Decisions & Staging for Next Steps

| Decision Area | Agreed Specification | Target Files |
| :--- | :--- | :--- |
| **2-Column Card Grid** | `project_panel` (`bgcolor=pal["bg"]`, borderless) containing `ft.Row(vertical_alignment=CrossAxisAlignment.START, spacing=16)` with 2 `expand=True` columns and 5 `_make_project_card()` containers. | `src/fc_engineer/app.py` |
| **Direct File Action Cards** | Remove `file_menu` (`PopupMenuButton`); expose `Load .ab`/`Save .ab`, `Export FASTA`/`Copy FASTA (30s Auto-Clear)`/`Export PNG`, and `Load Theme .yaml`/`Save Theme .yaml` as direct icon buttons. | `src/fc_engineer/app.py` |
| **UI Tests** | Assert zero `File Menu` `PopupMenuButton`, verify all 7 file action buttons and Light/Dark updates on `project_cards`. | `tests/test_ui.py` |
