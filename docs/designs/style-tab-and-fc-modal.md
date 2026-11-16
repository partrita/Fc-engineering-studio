# Style Tab Consolidation & Wide Fc Mutation Modal

* **Author:** fkt
* **Date:** 2026-11-16

## Overview
Two UI friction points remain in Fc Engineering Studio:
1. The `Project` tab contains a `Theme Files (.yaml)` card without explaining what the theme controls, while `Theme & Settings` duplicates `Load Theme .yaml` / `Save Theme .yaml`.
2. In the `Studio` tab's 360px side panel, Fc mutation preset checkboxes live inside a narrow `160px`-tall `ListView` that is hard to scan and select, crowding out `EU Diffs` and `Developability`.

This feature renames `Theme & Settings` to **`Style`**, moves `Theme Files (.yaml)` out of `Project` into `Style` with an explicit description (*"2D 모식도(도메인 색상·모양·크기·방향) 전용 테마 파일"*), and replaces the cramped inline Fc checkbox list in `Studio`'s Inspector with an **Active Fc Mutations** summary plus an **`Edit Fc Mutations...`** button (also triggered by double-clicking any domain on the 2D canvas) that opens a wide modal dialog.

## Analysis

### Key Contradictions & Resolutions (7-Lens Goldfish Synthesis)
1. **`page.show_dialog(snack)` vs. `page.show_dialog(dlg)` Inside the Modal:**
   * **Problem:** `notify()` uses `page.show_dialog(ft.SnackBar(...))`. Calling `notify()` while the Fc Mutations `AlertDialog` is open can replace or hide behind the active dialog.
   * **Resolution:** Display validation errors and confirmation messages in an **inline status text** inside the Fc Mutations `AlertDialog`, rather than calling `notify()` while the modal is open.
2. **Double-Clicking a Non-Fc Domain (e.g., `VL`, `VH`, `G4S`) on the 2D Canvas:**
   * **Problem:** All domain nodes share the same `ft.GestureDetector`. Double-clicking a light chain (`LC_1`) when `apply_to_all_heavy` is `False` means the selected chain has no `CH2`/`CH3` domain.
   * **Resolution:** `on_double_tap` selects `(cid, did)` first (`ctrl.select_domain(cid, did)`), opens the modal titled with the active chain (`Edit Fc Mutations — Chain: <name>`), and shows an inline hint if the selected chain lacks `CH2`/`CH3` while `apply_to_all_heavy` is `False`.
3. **2-Column Preset Grid Layout vs. Long Preset Labels:**
   * **Problem:** `ft.GridView` clips long preset labels like `YTE (M252Y/S254T/T256E)` due to fixed `child_aspect_ratio`.
   * **Resolution:** Pair filtered presets into `ft.Row`s of two `expand=True` `ft.Checkbox` items inside a scrollable `ft.Column` (`width=740`), giving ample horizontal room for every label.

## Goals vs. Non-Goals

### Goals
* Rename `"Theme & Settings"` tab to `"Style"`.
* Move `theme_files_card` (`"Theme Files (.yaml)"`) from `Project` into `Style`, remove duplicate inline `Load Theme .yaml` / `Save Theme .yaml` buttons, add the explicit schematic-only description, and keep Light/Dark card styling synchronized.
* Replace the cramped `160px` Fc checkbox list + inline ad-hoc/file-custom forms in `Studio`'s Inspector with a concise **Active Fc Mutations** summary + `Edit Fc Mutations...` button.
* Add `on_double_tap` on canvas domain `ft.GestureDetector`s to select the domain and open the Fc Mutations modal.
* Provide a wide `ft.AlertDialog` (`width=740`) with `Apply to all heavy chains` switch, live search filter, 2-column preset checkboxes, ad-hoc EU mutation input + `Apply`, file-custom mutation registration (`Save to Project File`), and inline status feedback.
* Update/add tests in `tests/test_ui.py` covering the `"Style"` tab, moved `"Theme Files (.yaml)"` card, Inspector summary, and Fc Mutations modal (button + `on_double_tap`).

### Non-Goals
* Modifying `designer.py`, `config.py`, `exporters.py`, or `.ab`/`.yaml` schemas.
* Changing the `"Customize Theme"` modal (`on_open_theme_editor`).

## Alternatives Weighed

1. **Modal Preset Grid: `ft.Row` Pairs vs. `ft.GridView`:**
   * **Chosen:** `ft.Row` pairs (`expand=True`) inside a scrollable `ft.Column` avoids `GridView` aspect-ratio clipping.
   * *Rejected by author: `ft.GridView` clips long mutation strings when `child_aspect_ratio` is fixed.*
2. **Modal Error Feedback: Inline Status Text vs. `SnackBar`:**
   * **Chosen:** Inline status `ft.Text` inside the modal avoids `page.show_dialog` slot collisions.
   * *Rejected by author: `page.show_dialog(ft.SnackBar(...))` can dismiss or conflict with the open `AlertDialog`.*

## Key Decisions & Staging for Next Steps

| Decision Area | Agreed Specification | Target Files |
| :--- | :--- | :--- |
| **`Style` Tab & Card Move** | Rename `"Theme & Settings"` to `"Style"`, move `theme_files_card` from `Project` into `Style` with schematic-only description, remove duplicate buttons, keep Light/Dark sync. | `src/fc_engineer/app.py` |
| **Inspector Summary + Fc Modal** | Replace inline `160px` checkbox list with active preset summary + `Edit Fc Mutations...` button and canvas `on_double_tap`, opening a `740px`-wide `ft.AlertDialog` with 2-column presets and inline status text. | `src/fc_engineer/app.py` |
| **Headless UI Tests** | Test `"Style"` tab, moved `Theme Files (.yaml)` card, `Edit Fc Mutations...` modal & `on_double_tap`, preset toggling, and ad-hoc/file-custom mutations. | `tests/test_ui.py` |
