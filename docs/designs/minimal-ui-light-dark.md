# Minimal Top-Tabbed Studio UI & Light/Dark Theme Support

* **Author:** fkt
* **Date:** 2026-11-16

## Overview
Fc Engineering Studio currently displays the toolbar, domain palette, pairing builder, 2D schematic canvas, chain strip, and right-hand inspector all at once on a single screen locked to dark mode. This feature reorganizes the Flet desktop interface into a minimal top-tabbed layout:
1. A **`Studio`** tab showing only the 2D interactive schematic canvas and the `N -> C` polypeptide chain strip, with domain/chain editing and inspection accessible through contextual controls or secondary tabs.
2. A top **`File`** dropdown menu (`PopupMenuButton`) consolidating `.ab` load/save, `FASTA` export/copy, `.png` export, and `.yaml` theme load/save.
3. A **`Theme & Settings`** tab dedicated to switching between Light and Dark mode (`ft.ThemeMode.LIGHT` / `ft.ThemeMode.DARK`) and customizing `SchematicTheme` presets, shapes, colors, orientation, and spacing.
4. Synchronized Light/Dark polarity across the Flet app shell, 2D canvas background, and chain header labels, while keeping transparent `.png` exports legible.

## Analysis

### Key Contradictions & Resolutions (7-Lens Goldfish Synthesis)
1. **Minimal `Studio` Tab vs. Full Authoring Workflow:**
   * **Problem:** Removing the left palette/pairing sidebar and right inspector from the `Studio` tab frees screen space for the 2D schematic, but users still need a way to add domains, configure pairings, and apply Fc mutations.
   * **Resolution:** Keep the `N -> C` chain strip directly below the 2D canvas in `Studio` for direct drag-and-drop reordering and selection, and move detailed `Domain Palette`, `Pairings`, and `Inspector / Fc Mutations` into a dedicated `Inspector & Domain` tab (or contextual modal/drawer) so they never crowd the `Studio` tab.
2. **Unified Light/Dark Mode vs. Custom `SchematicTheme` & Transparent PNG Exports:**
   * **Problem:** If Light/Dark mode is coupled to `SchematicTheme` presets, toggling Light/Dark would overwrite custom domain colors. In addition, using light-colored chain labels on a dark canvas produces illegible text if exported onto a transparent `.png` pasted into white documents.
   * **Resolution:** Keep `theme_mode` (`"light"` / `"dark"` polarity on `AntibodyStudioController` and `page.theme_mode`) orthogonal to `SchematicTheme` (domain fill colors, borders, shapes, orientation, spacing). Switching `theme_mode` updates the Flet shell surfaces, 2D canvas background (`#F8FAFC` in Light, `#141413` in Dark), and chain header text (`#0F172A` in Light, `#FAF9F5` in Dark), while `export_schematic_png()` continues to render dark-ink chain labels on a transparent RGBA background.
3. **Headless UI Testability (`tests/test_ui.py`) vs. `ft.Tabs` and `ft.PopupMenuButton`:**
   * **Problem:** `_walk_controls()` in `tests/test_ui.py` currently traverses `.controls` and `.content`, which would miss controls nested inside `ft.Tabs.tabs` or `ft.PopupMenuButton.items`.
   * **Resolution:** Keep all state mutations on `AntibodyStudioController` decoupled from Flet widgets, and update `_walk_controls()` in `tests/test_ui.py` to also traverse `.tabs` and `.items` if present.

## Goals vs. Non-Goals

### Goals
* Minimal `Studio` tab containing only the 2D schematic canvas and `N -> C` chain strip.
* Top `File` dropdown menu (`PopupMenuButton`) grouping `.ab`, `FASTA`, `PNG`, and `.yaml` theme actions.
* Dedicated `Theme & Settings` tab housing Light/Dark mode toggle, `SchematicTheme` preset selector, per-domain style editor, orientation, and spacing sliders.
* Synchronized Light/Dark mode switching across `page.theme_mode`, panel surfaces, 2D canvas background, and chain labels.
* Full preservation of headless controller and UI tests in `tests/test_ui.py`.

### Non-Goals
* Modifying core EU numbering, mutation routing, or `.ab`/`FASTA` serialization in `core.py`, `designer.py`, or `exporters.py`.
* Modifying read-only `src/fc_engineer/data/*.yaml` files at runtime.
* Introducing multi-window or web-only Flet layouts.

## Alternatives Weighed

1. **Top `ft.Tabs` vs. Left `NavigationRail`:**
   * **Chosen:** Top-tabbed layout (`Studio`, `Theme & Settings`, etc.) preserves full horizontal width for the 2D schematic canvas and `N -> C` chain strip.
   * *Rejected by author: Left `NavigationRail` consumes horizontal space needed for multi-chain bispecific schematics.*
2. **Top `File` Dropdown (`PopupMenuButton`) vs. Flat Toolbar Buttons:**
   * **Chosen:** Top `File` dropdown menu reduces toolbar clutter and groups all file I/O actions in one place.
   * *Rejected by author: A flat row of 6+ load/save/export buttons crowds the header alongside project name and format selector.*
3. **Orthogonal `theme_mode` + `SchematicTheme` vs. Coupled Light/Dark Presets:**
   * **Chosen:** Orthogonal `theme_mode` (polarity) + `SchematicTheme` (domain palette/geometry) lets any custom domain theme work in both Light and Dark mode without duplicating presets.
   * *Rejected by author: Coupling Light/Dark to `SchematicTheme` presets would overwrite user-customized domain colors when toggling Light/Dark mode.*

## Key Decisions & Staging for Next Steps

| Decision Area | Agreed Specification | Target Files |
| :--- | :--- | :--- |
| **Top-Tabbed Layout** | `Studio` tab (2D canvas + `N -> C` chain strip only), `Inspector & Domain` tab (palette, pairings, Fc mutations, sequence/developability), `Theme & Settings` tab (Light/Dark toggle + `SchematicTheme` editor + spacing). | `src/fc_engineer/app.py` |
| **Top `File` Menu** | `ft.PopupMenuButton("File")` in the top bar for `Load .ab`, `Save .ab`, `Export FASTA`, `Copy FASTA (30s Auto-Clear)`, `Export PNG`, `Load Theme .yaml`, `Save Theme .yaml`. | `src/fc_engineer/app.py` |
| **Light/Dark Polarity** | `AntibodyStudioController.theme_mode` (`"dark"` / `"light"`) drives `page.theme_mode`, panel colors, canvas background, and chain label contrast without altering `SchematicTheme` domain colors or transparent PNG ink contrast. | `src/fc_engineer/app.py` |
| **Headless UI Tests** | Extend `_walk_controls()` in `tests/test_ui.py` to traverse `.tabs` and `.items` and add tests for Light/Dark toggling. | `tests/test_ui.py` |
