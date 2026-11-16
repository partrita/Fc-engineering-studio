# 100 UX Improvements Catalog — Fc Engineering Studio

* **Author:** fkt
* **Date:** 2026-11-16
* **Status:** Catalog / Partially Implemented (Canonical `[x]` / `[~]` / `[ ]` tracker; supersedes the initial draft in [design_draft.md](./design_draft.md); remaining `[ ]` and `[~]` items require dedicated follow-up specs)

## Overview
This document catalogs 100 codebase-grounded UX improvements for Fc Engineering Studio, evenly distributed across 10 functional areas (`10` improvements per area). All items map to the 3-tab Flet desktop studio (`Studio`, `Project`, `Style`) and its 5 core modules: `app.py` (Flet UI & `AntibodyStudioController`), `designer.py` (domain/chain topology & 2D scene layout), `core.py` (EU numbering, mutations & developability), `config.py` (YAML presets & user preferences), and `exporters.py` (`.ab`, FASTA, GenBank, CSV & PNG I/O).

Status legend: `[x]` = shipped in codebase, `[~]` = partially shipped (remaining sub-task noted), `[ ]` = backlog.

## Problem Statement
Earlier UX audits of Fc Engineering Studio identified three cross-cutting friction points:
1. **Freeform Canvas Dragging vs. Auto-Reflow (`dom.pinned`):** `on_pan_update` sets `dom.pinned = True` to opt a domain out of `reflow_chain_coordinates()`. While `ctrl.unpin_domain()` and `ctrl.reset_layout()` are now exposed in the Inspector, `top_offset = max(0.0, node.y - 14.0)` in `canvas_overlay` can clip N-terminus labels near `y = 0`, and `canvas_overlay` bounds need `on_pan_end` viewport padding normalization.
2. **Destructive Transitions & Undo/Redo:** Bounded 50-snapshot `_undo_stack`/`_redo_stack`, `ctrl.is_dirty` format-switch confirmation, and window-close `Save / Discard / Cancel` guards are now implemented, whereas `Remove Chain` and `Delete Domain` still execute without an Undo `SnackBar` or confirmation.
3. **Remaining UI Polish & Ergonomics:** Focus-gated `Delete` shortcuts, real-time alignment guides, canvas-based `Pairings` picking, inline EU/allotype/conflict hints, and collapsible side-panel drawers remain in the backlog below.

## Goals
* **Catalog Completeness:** Maintain 100 non-redundant, codebase-grounded UX improvements evenly organized across 10 functional areas (`10` items per area) with clear `[x]` / `[~]` / `[ ]` implementation status.
* **Technical Accuracy:** Keep every symbol reference, API signature, and security constraint aligned with `src/fc_engineer/`.

## Non-Goals
* Breaking `.ab` `schema_version: 1` or `.yaml` theme file backward compatibility.
* Introducing full multi-language i18n frameworks (standardizing UI copy in English is sufficient).

## Alternatives Considered
* **Flat uncategorized list:** *Rejected:* Uncategorized lists obscure gaps in export, error recovery, and accessibility.
* **Top-10 only shortlist:** *Rejected:* A 100-item catalog across 10 areas provides comprehensive coverage before selecting priority implementation batches.
* **Chosen:** **10x10 Categorized Catalog with Status Markers (`[x]` / `[~]` / `[ ]`)** — Preserves stable `#1`–`#100` IDs for traceability with [design_draft.md](./design_draft.md) and tests while separating shipped work from partial and backlog items.

---

## 100 UX Improvements (10 Areas x 10 Improvements)

### 1. 2D Canvas & Schematic Interaction (`Studio — app.py`, `designer.py`)
1. [x] Add an **Unpin** button (`ctrl.unpin_domain()`) and a **Reset Layout** button (`ctrl.reset_layout()`) in `rebuild_inspector()`.
2. [x] Render a pin badge (`📌`) on `SchematicNode`s where `dom.pinned == True`.
3. [ ] Normalize viewport container padding from `scene.bounds` (`min_x`, `min_y`) in `canvas_overlay` on `on_pan_end` using a shared `world_to_screen` / `screen_to_world` transform across `GestureDetector`, `cv.Canvas`, `zoom_scale`, `snap_to_grid`, and canvas picking (`#44`) to prevent hit-target desync or N-terminus label clipping near `(0, 0)`.
4. [~] Shipped: zoom controls (`0.5x–2.0x` via `ctrl.zoom_scale`); Remaining: **Fit to View** (`min(vp_w / bounds_w, vp_h / bounds_h)` clamped to `0.5–2.0`, resetting scroll offsets) and middle-drag viewport pan.
5. [~] Shipped: `10px` snap-to-grid (`ctrl.snap_to_grid`) on `on_pan_end`; Remaining: real-time alignment guide lines (`±5px` center/edge snap tolerance) on a dedicated `cv.Canvas` overlay.
6. [ ] Make chain header labels (`node.is_n_terminus`) draggable handles to translate an entire `PolypeptideChain` at once (setting `dom.pinned = True` on all chain domains until `reset_layout()`) (Optional).
7. [~] Shipped: clickable `-S-S-` disulfide badges to remove pairings; Remaining: highlight both endpoint domains and switch to `Studio > Pairings sub-tab` on edge click.
8. [x] Make `on_double_tap` context-aware: open `open_fc_mutations_modal` for `CH1`/`Hinge`/`CH2`/`CH3`, and select and show an Inspector-edit prompt for `VH`/`VL`/`CL`/`G4S`.
9. [~] Shipped: hover tooltips on canvas domain nodes showing chain name, domain type, and active mutations; Remaining: append `isotype/allotype`, `template_key`, and residue length.
10. [~] Shipped: `rebuild_canvas_only()` during `on_pan_update` and full `render_workspace()` on `on_pan_end` (updating pairing edges on drop); Remaining: `~16–33ms` drag-event throttling and recording an undo snapshot and dirty state in `move_domain()` on `on_pan_end`.

### 2. Chain & Domain Assembly (`Studio — rebuild_chain_strip`, `designer.py`)
11. [x] Add dialog-based chain renaming on each chain row in `rebuild_chain_strip()` backed by `ctrl.rename_chain()`.
12. [x] Add chain reorder controls (up/down) via `ctrl.reorder_chain()` to adjust chain order.
13. [x] Place a **Clone Chain** icon button directly on each chain row in `rebuild_chain_strip()`.
14. [x] Provide a compact `86px` Left Domain Palette (`ft.Draggable`) in `Studio` for drag-to-insert at any slot index.
15. [ ] Highlight `slot_target` borders/backgrounds via `on_will_accept` / `on_leave` during chip drag-and-drop.
16. [x] Add `◀` / `▶` domain reorder buttons in the Inspector wired to `ctrl.reorder_domain(chain_id, from_index, to_index)`.
17. [ ] Replace the fixed `height=168` chain strip container with a collapsible/resizable splitter (preserving `min_height >= 300px` for `canvas_box`) for 8+ chain multispecifics (Optional).
18. [x] Show mutation (amber) and resolution-error (red) indicator dots on domain chips in `rebuild_chain_strip()`.
19. [~] Shipped: domain count and total `aa` length next to `{chain.name}:`; Remaining: show an explicit empty-chain warning badge.
20. [x] Add multi-domain macro items to `quick_add_btn` (`+`): **Insert scFv (VH-G4S-VL)** and **Insert Fc (Hinge-CH2-CH3)**.

### 3. Fc Mutations & EU Numbering (`Studio Modal — open_fc_mutations_modal`, `core.py`)
21. [x] Group presets in `open_fc_mutations_modal` by biological category (`All`, `Half-Life`, `Effector`, `Heterodimer`, `Glyco/Other`, `Project`).
22. [ ] Show a live `WT -> Mut` preview and EU-position conflict warning before toggling overlapping presets.
23. [~] Shipped: removable `ft.Chip`s for `dom.custom_mutations` in the Inspector; Remaining: mirror removable chips inside `open_fc_mutations_modal`.
24. [ ] Display the affected chain list (e.g., `Applies to: HC_1, HC_2`) when `apply_to_all_heavy` is ON.
25. [ ] Show EU range hints (`CH1: 118-215`, `Hinge: 216-230`, `CH2: 231-340`, `CH3: 341-447`) and IgG2/3/4 hinge gap warnings inline.
26. [x] Auto-suggest the WT amino acid when typing an EU position in `custom_eu_tf` via `lookup_wt_residue()` (e.g., `234` -> `WT @ EU 234: L`).
27. [~] Shipped: delete actions for project-scoped `file_custom_mutations` inside `open_fc_mutations_modal`; Remaining: inline edit support.
28. [ ] When switching `isotype` rolls back due to WT mismatch, list the conflicting mutation and offer **Clear incompatible & switch**.
29. [~] Shipped: `dom.mutation_history` audit entries displayed in the Inspector; Remaining: wrap in a collapsible `ExpansionTile`.
30. [x] Add a **Reset Chain to WT** button in both `open_fc_mutations_modal` and Inspector via `ctrl.reset_chain_mutations()`.

### 4. Sequence, Diff & Developability Analysis (`Inspector — rebuild_inspector`, `core.py`)
31. [x] Replace plain `ft.Text("Developability: ...")` with color-coded liability count badges.
32. [x] Add a `Scan VH/VL/G4S` toggle (`ctrl.include_variable_dev`) calling `analyze_chain_developability(..., include_variable=True)` using 1-based local residue indices for non-constant domains.
33. [~] Shipped: color-code `Assembled Chain Sequence` by domain span (`ft.TextSpan` per `assemble_chain_sequence` segment); Remaining: bold/underline mutated residues.
34. [x] Add a **Copy Chain Sequence** button (via `ctrl.copy_text_secure`) next to chain metrics.
35. [x] Make `EU Diffs` clickable `ft.Chip`s (`CH2:L234A`) that open `open_fc_mutations_modal` for the chain.
36. [~] Shipped: live residue counter (`len / 2000`, matching `MAX_CUSTOM_SEQ_LENGTH = 2000`) in `seq_tf`; Remaining: auto-strip FASTA headers/whitespace on paste and flag non-`A-Z` characters.
37. [x] Add a **Reset to Template** button next to `Apply Sequence` to clear `dom.custom_sequence`.
38. [x] Display calculated molecular weight (`MW` in kDa) and theoretical isoelectric point (`pI`) via `compute_protein_metrics()`.
39. [ ] Show a concise allotype sequence difference note in `rebuild_inspector()` (e.g., `G1m1,17` vs. `nG1m1`).
40. [ ] Separate **Chain Summary** and **Selected Domain** into collapsible `ExpansionTile`s so long sequences do not push domain controls out of view.

### 5. Multispecific & Heterodimer Engineering (`Project & Studio — designer.py`)
41. [ ] Mirror `Multispecific & Heterodimer Builder` (currently in `Project`) in a compact third sub-tab (`Multispecific`) inside `Studio`'s side panel for immediate visual feedback (Optional).
42. [x] Persist `wiz_state` (`a`, `b`, `mode`) on `ctrl` so `render_workspace()` does not reset dropdown selections.
43. [ ] Warn before `swap_crossmab_ch1_cl()` if either `CH1` or `CL` carries mutations that will be cleared when domain identity changes.
44. [ ] Add a **Pick on Canvas** mode in `Studio > Pairings sub-tab` so users can click Domain A and Domain B directly on the 2D schematic.
45. [ ] Disambiguate `endpoint_options` in `Studio > Pairings sub-tab` dropdowns with domain index and `template_key` (e.g., `HC_1 [2] VH (VH_B)`).
46. [ ] Flag unpaired variable/constant Fab domains (excluding `G4S` linkers and single-domain `BiTE`/`scFv` fragments) or mismatched `VH_VL` pairings with a warning icon.
47. [ ] Add `mode="vh_vl"` (`VH <-> VL`) and `mode="fab"` (full `Fab`) crossover options alongside `CH1 <-> CL` CrossMab swap, with explicit `custom_sequence`/`template_key`/`allotype` and `InterChainPairing` validity rules.
48. [x] Add an **Auto-Detect** button in `Studio > Pairings sub-tab` calling `ctrl.auto_detect_pairings()`.
49. [x] Register `ChargePair_A` (`K392D/K409D`) and `ChargePair_B` (`D356K/D399K`) in `mutants.yaml` and wire `mode="charge_pair"` in `apply_heterodimer_pair()`.
50. [ ] Show a chain-count/topology preview summary for each format template before switching.

### 6. Project, File I/O & Export (`Project — exporters.py`, `app.py`)
51. [x] Expose **Export GenBank (.gb)** (`save_genbank_file`) and **Export CSV (.csv)** (`save_csv_report_file`) buttons in `export_card`.
52. [x] Add a PNG resolution selector (`1x`, `2x`, `3x`, `4x`) and background color selector (`Transparent`, `#FFFFFF`, `#141413`) for `export_schematic_png()`, with `MAX_FILE_SIZE` (1 MB) enforced by `save_schematic_png()`.
53. [x] Load a scalable TrueType font (`_load_scaled_font`) in `export_schematic_png()` instead of fixed bitmap `ImageFont.load_default()`.
54. [x] Track `ctrl.current_file_path` for direct **Save `.ab`** overwrite plus a separate **Save As...** action.
55. [x] Track `ctrl.is_dirty` and display `*` in `page.title` and `status_bar_row` when unsaved changes exist.
56. [x] Persist a **Recent Projects** list (last 5 `.ab` paths) in `ab_files_card` via `~/.fc_engineer/prefs.yaml`.
57. [ ] Validate `format_project_fasta` / `serialize_ab_project` *before* opening `file_picker.save_file` so errors fail fast.
58. [~] Shipped: synchronous 30-second clipboard auto-clear countdown indicator in `status_bar_row` updated on `render_workspace()`; Remaining: add an inline **Clear Now** button (live background ticks tracked in `#98`).
59. [x] Expose `core.generate_batch()` in a **Batch Variant Generator** modal (`ctrl.export_batch_csv`) with CSV formula-injection sanitization (capped at `MAX_BATCH_ITEMS = 20` combinations per export).
60. [x] Allow customizing FASTA header format (`project_chain`, `chain_only`, `detailed`) in `export_card`.

### 7. Style, Theme & Visual Polish (`Style — on_open_theme_editor`, `_ui_palette`)
61. [x] Expose `domain_spacing` (`4..80px`) and `chain_spacing` (`60..260px`) sliders in `Style` (`ctrl.set_spacing()`).
62. [ ] Add instant canvas preview updates while editing colors/shapes inside `on_open_theme_editor`.
63. [~] Shipped: clickable color swatch pickers beside `#RRGGBB` text fields in `on_open_theme_editor`; Remaining: add inline hex validation error text.
64. [x] Add an **Auto Contrast** button in `on_open_theme_editor` that picks `#FAF9F5` or `#141413` from `fill_color` luminance.
65. [ ] Repaint stale `ft.Text` colors inside `settings_panel` and `project_cards` during `on_theme_mode_toggle()`.
66. [x] Set `theme_dd.value = None` when `on_apply_dialog` applies modified custom theme settings.
67. [x] Add a compact `Orient: V / H` toggle button directly on `studio_toolbar`.
68. [ ] Allow saving customized themes into `~/.fc_engineer/themes/*.yaml` and reloading them into `ctrl.theme_presets` (Optional).
69. [ ] Standardize `theme_files_card` copy in English to match the rest of the `Style` tab labels.
70. [ ] Warn in `Style` when exported PNG label contrast against the selected `bg_color` (or `#FFFFFF` preview) falls below `WCAG 4.5:1`.

### 8. Navigation, Shortcuts & Productivity (`Studio & Global — main_tabs`, `page`)
71. [x] Add a compact `studio_toolbar` quick-action strip (Save, Undo/Redo, Zoom, Snap Grid, Orient, Batch, Command Palette, Help) inside `Studio`.
72. [~] Shipped: global keyboard shortcuts (`Cmd/Ctrl+S`, `Cmd/Ctrl+O`, `Cmd/Ctrl+Z`, `Cmd/Ctrl+Shift+Z`, `Cmd/Ctrl+K`); Remaining: focus-gated `Delete` (active only when no `TextField` has focus) and `Cmd/Ctrl+Shift+C`.
73. [ ] Support focus-gated arrow keys (aligned with active `Vertical`/`Horizontal` canvas orientation) to step domain selection along a chain and `Tab` to switch chains (Optional).
74. [x] Add a `Cmd/Ctrl+K` quick-jump Command Palette to trigger actions and dialogs.
75. [~] Shipped: `Project > Chain > Domain` breadcrumbs at the top of `inspector_col`; Remaining: make breadcrumb segments clickable.
76. [ ] Make the status indicator clickable to jump to `Project` or trigger `Save .ab`.
77. [x] Add a persistent bottom status bar (`status_bar_row`) showing active format, chain/domain/pairing counts, dirty state, and clipboard timer.
78. [ ] Wire `Enter` (`on_submit`) on modal inputs and `Escape` to close `AlertDialog`s.
79. [ ] Set `autofocus=True` on `search_tf` when `open_fc_mutations_modal` opens.
80. [ ] Auto-scroll `canvas_box` to keep the newly selected domain visible when clicked from `chain_strip_col`.

### 9. Onboarding, Guidance & Error Recovery (`Feedback — notify`, `select_format`, `load_ab`)
81. [x] Prompt confirmation on `on_format_select` when `ctrl.is_dirty` is True before replacing `ctrl.project`.
82. [ ] Prompt a confirmation dialog for bulk or destructive operations (e.g., `Reset Chain to WT` when custom mutations exist, or clearing all pairings).
83. [ ] Include actionable buttons (e.g., **Select Error Domain**) in error `SnackBar`s.
84. [ ] Show the exact `resolve_domain_sequence` error message in an Inspector banner when clicking an `error_flag` red node.
85. [x] Render a UI-only empty-canvas placeholder in `canvas_overlay` when no domains exist (excluded from `build_schematic_scene()` and `export_schematic_png()`).
86. [ ] Show a structured diagnostic dialog when `load_ab()` fails validation (`ORIGIN` mismatch, `schema_version`, `1 MB` limit).
87. [ ] Display a 1-line biological summary under `format_dd` explaining each antibody format template.
88. [x] Add an in-app `?` **Keyboard Shortcuts & Tips** dialog (`open_help_modal`).
89. [ ] Open a selectable copy-text fallback dialog if `pyperclip.copy` fails on headless/restricted systems (Optional).
90. [ ] Highlight invalid numeric inputs in red inside `on_open_theme_editor` instead of silently ignoring `ValueError`.

### 10. Undo/Redo, Desktop Ergonomics & Accessibility (`Controller — AntibodyStudioController`, `shutdown`)
91. [x] Implement bounded `_undo_stack` / `_redo_stack` (`50` deep-copy `AntibodyProject` snapshots) across `AntibodyStudioController` mutations.
92. [ ] Pair single-item deletions (`Remove Chain`, `Delete Domain`) with a `4s` non-blocking **Undo** action button in `SnackBar` (calling `ctrl.undo()` when the stack top matches).
93. [x] Prompt `"Save & Close / Discard & Close / Cancel"` in `on_window_event` (`CLOSE`) when `ctrl.is_dirty` is True.
94. [~] Shipped: write `~/.fc_engineer/recovery.ab` snapshots on dirty workspace updates when `auto_recovery` is enabled, clear on clean save/exit, and offer restore on startup; Remaining: debounced atomic writes (`tempfile` + `os.replace`) and skipping empty-chain snapshots (`parse_ab_project()` requires `>= 1` domain per chain).
95. [~] Shipped: persist `recent_files`, `theme_mode`, `snap_to_grid`, and `include_variable_dev` in `~/.fc_engineer/prefs.yaml`; Remaining: persist window dimensions and `side_panel_tabs` selection.
96. [ ] Make `studio_side_panel` width (`300px–560px`) user-resizable with a draggable splitter, clamped so `canvas_box` keeps `min_width >= 400px` alongside the `86px` Left Domain Palette (Optional).
97. [ ] Support dragging a single `.ab` regular file from Finder/Explorer onto the app window (`page.on_file_drop`) with `MAX_FILE_SIZE` (1 MB) and dirty-state validation (Optional).
98. [ ] Dispatch background `threading.Timer` countdown ticks via an `async` `page.run_task` coroutine guarded by a live-session/shutdown check so Flet controls are never mutated off-thread (Internal async sub-task of `#58`).
99. [ ] Collapse `studio_side_panel` into a toggleable drawer on narrow desktop windows (`< 1100px`) while preserving `canvas_box` `min_width >= 400px`.
100. [x] Add an **Open Folder** action button in `SnackBar` after saving `.ab`, `.fasta`, `.gb`, `.csv`, or `.png` files.

---

## Architecture & Security Constraints
* **Undo/Redo Snapshots (`_undo_stack` / `_redo_stack`):** Bounded to `50` `copy.deepcopy(AntibodyProject)` entries in `AntibodyStudioController`. View-only state (`zoom_scale`, `snap_to_grid`, `theme_mode`) does not push undo snapshots; `load_ab()` clears both stacks. Recording a single undo snapshot on `move_domain()` (`on_pan_end`) is tracked in `#10`.
* **Persistence & File Bounds:**
  * User preferences (`~/.fc_engineer/prefs.yaml`) are parsed via `config.py` using `NoAliasSafeLoader` (which disables YAML anchors/aliases) and capped at `MAX_FILE_SIZE` (1 MB).
  * Custom theme directory (`~/.fc_engineer/themes/*.yaml` in `#68`) must enforce filename/symlink sanitization, `NoAliasSafeLoader`, `1 MB` size limit, and preset-count bounds.
  * Crash-recovery snapshots (`~/.fc_engineer/recovery.ab`) use `.ab` GenBank-style serialization (`exporters.save_ab_file` / `parse_ab_project`, `1 MB` limit) and are deleted on clean save/exit; debounced atomic `tempfile` + `os.replace` writes and empty-chain guarding are tracked in `#94`.
* **Clipboard & Thread Safety:** `copy_text_secure()` wipes OS clipboard contents after 30 seconds via generation-tokenized `threading.Timer` and `shutdown()`. `#58` currently renders synchronously on `render_workspace()`; live timer ticks in `#98` must marshal via an `async` `page.run_task` coroutine with a session-alive guard.

## Glossary
* **Immunology & Sequence Terms:**
  * **Fc / CH1 / CH2 / CH3 / Hinge:** Constant regions of an `IgG` (Immunoglobulin G) heavy chain; `CH2`/`CH3` form the Fc (Fragment crystallizable) region mediating `FcRn` (neonatal Fc receptor) half-life and effector functions (`ADCC` = antibody-dependent cellular cytotoxicity, `CDC` = complement-dependent cytotoxicity).
  * **EU Numbering / Isotype / Allotype / Hinge Gaps:** Standardized residue numbering scheme (`118–447`) for human IgG constant domains across subclasses (`IgG1–4`) and polymorphic allotype sequence variants (`G1m1,17`, `nG1m1`). `IgG2/3/4 hinge gaps` refer to EU positions (`223–225` / `223–230`) absent in IgG2/3/4 relative to IgG1.
  * **VH / VL / CL / G4S / scFv / Fab / BiTE:** Variable heavy (`VH`), variable light (`VL`), constant light (`CL`), flexible glycine-serine linker (`G4S`), single-chain variable fragment (`scFv`), antigen-binding fragment (`Fab`), and bispecific T-cell engager (`BiTE`).
  * **CrossMab / ChargePair:** Bispecific engineering strategies—`CrossMab` swaps `CH1`/`CL` domains to enforce light-chain pairing; `ChargePair_A`/`ChargePair_B` (`K392D/K409D` vs. `D356K/D399K`) drive `CH3` heterodimerization.
  * **PTM / CDR / WT / MW / kDa / pI / aa / `-S-S-`:** Post-translational modification liability (`PTM`), complementarity-determining region (`CDR`), wild-type (`WT`), molecular weight (`MW` in kilodaltons `kDa`), theoretical isoelectric point (`pI`), amino acid (`aa`), and disulfide bond (`-S-S-`).
  * **`.ab` Format & Batch Variants:** GenBank-style project file (`schema_version: 1`, `1 MB` max size, `ORIGIN` sequence integrity check); `core.generate_batch()` expands comma/newline mutation sets (up to `MAX_BATCH_ITEMS = 20`) into multi-variant sequences.
* **Codebase & UI Symbols:**
  * `ctrl` (`AntibodyStudioController`), `dom` (`DomainBlock` with `template_key` for preset `VH`/`VL`/`CL`/`G4S` templates), `PolypeptideChain`, `InterChainPairing` (inter- or intra-chain pairing link), `HC` (Heavy Chain), `SchematicNode` (2D layout node with `is_n_terminus` for chain start), `cv` (`flet.canvas`), `wiz_state` (heterodimer wizard `{a, b, mode}`), `file_custom_mutations` (project-scoped custom presets stored in `.ab`), `slot_target` (chain-strip `ft.DragTarget`), `resolve_domain_sequence` / `error_flag` (sequence validation error state), and `_tf` (`ft.TextField`) / `_dd` (`ft.Dropdown`).

## Testing & Verification
* **Implementation History & Next Specs:** Initial 100-item design draft is superseded and archived in [design_draft.md](./design_draft.md); this catalog is the canonical status tracker, and remaining `[ ]` / `[~]` items require dedicated implementation specs before coding.
* **Verification Suite:** Verified via `uv run pytest` (`tests/test_logic.py`, `tests/test_exporters.py`, `tests/test_ui.py` — 112 tests passing) and `uv run bandit -r src/` (0 security issues).

## Next Steps
* **Core Ergonomics Backlog (All Remaining `[ ]` and `[~]` Items):**
  * *Canvas & Assembly:* `#3`, `#4`, `#5`, `#7`, `#9`, `#10`, `#15`, `#19`.
  * *Fc Mutations & Inspector:* `#22`, `#23`, `#24`, `#25`, `#27`, `#28`, `#29`, `#33`, `#36`, `#39`, `#40`.
  * *Multispecific & Export:* `#43`, `#44`, `#45`, `#46`, `#47`, `#50`, `#57`, `#58`.
  * *Style, Shortcuts, Recovery & Ergonomics:* `#62`, `#63`, `#65`, `#69`, `#70`, `#72`, `#75`, `#76`, `#78`, `#79`, `#80`, `#82`, `#83`, `#84`, `#86`, `#87`, `#90`, `#92`, `#94`, `#95`, `#98`, `#99`.
* **Optional / Low-Priority Backlog:** `#6`, `#17`, `#41`, `#68`, `#73`, `#89`, `#96`, `#97` remain optional and are deferred unless explicitly requested.
