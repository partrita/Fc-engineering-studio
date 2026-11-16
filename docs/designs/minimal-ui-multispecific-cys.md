# Minimal 4-Tab Studio UI, Cysteine Bond Visualization & Multispecific Design Suite

* **Author:** fkt
* **Date:** 2026-11-16

## Overview
Fc Engineering Studio currently places `Project Name`, `Format Template`, and `Add Chain` in the top bar alongside spacing sliders inside `Studio`, and `Domain Palette` cards use a hardcoded light background (`#E2E8F0`) with unreadable text in dark mode. In addition, the schematic only draws double lines for `HINGE_DISULFIDE` and lacks dedicated tools for rapid multispecific antibody assembly. This feature makes four changes:
1. It reorganizes the UI into a minimal top bar (`File` dropdown menu only) and 4 tabs: **`Studio`** (pure 2D schematic canvas + `N -> C` chain strip, with all spacing sliders/buttons removed), **`Project`** (`Project Name`, `Format Template`, `Add Chain`, `Clone Chain`, and `Heterodimer Fc Wizard`), **`Inspector & Domain`** (`Domain Palette`, `Pairings`, `CrossMab CH1<->CL Swap`, domain properties, and Fc mutations), and **`Theme & Settings`** (Light/Dark toggle, schematic theme presets, customizer, and YAML load/save).
2. It fixes dark-mode contrast across `Domain Palette` cards, chain strip chips, panels, and labels using a dynamic theme palette.
3. It adds `-S-S-` cysteine/disulfide bond visualization across both the Flet 2D canvas and Pillow `.png` exports for `HINGE_DISULFIDE`, `CH1_CL`, explicit `DISULFIDE_BOND` pairings, and engineered `CH3` `S354C`+`Y349C` pairs.
4. It adds multispecific design tools: 3 new format templates (`two_plus_one`, `trispecific`, `fab_igg`), one-click `Heterodimer Fc Wizard`, `Clone Chain`, and `CrossMab CH1<->CL Swap`.

## Analysis

### Key Contradictions & Resolutions (7-Lens Goldfish Synthesis)
1. **Static `Domain Palette` Colors vs. Dynamic Light/Dark Switching:**
   * **Problem:** `palette_cards` in `src/fc_engineer/app.py` hardcodes `bgcolor="#E2E8F0"` while `ft.Text(dtype)` inherits white text in `ThemeMode.DARK`, making labels invisible.
   * **Resolution:** Rebuild `palette_cards` inside a `rebuild_palette()` helper triggered on `render_workspace()` so `bgcolor`, `border`, and `ft.Text(color=...)` use `_ui_palette(ctrl.theme_mode)` in both Dark and Light modes.
2. **Cysteine Bond Detection vs. `.ab` `schema_version: 1` Compatibility:**
   * **Problem:** Adding a new top-level `disulfide_bonds` field to `.ab` JSON manifests would break `schema_version: 1` parsing.
   * **Resolution:** Add `"DISULFIDE_BOND"` to `VALID_PAIRING_TYPES` in `src/fc_engineer/designer.py` and append `is_disulfide: bool = False` to `SchematicPairingEdge`. `build_schematic_scene()` marks `is_disulfide = True` for `HINGE_DISULFIDE`, `CH1_CL`, `DISULFIDE_BOND`, and `CH3_DIMER` pairs carrying complementary `S354C`/`Y349C` mutations. Both `app.py` and `exporters.py` draw solid gold/amber `-S-S-` bridges with a `-S-S-` midpoint label.
3. **Multispecific Rapid Design Tools vs. Minimal `Studio` Tab:**
   * **Problem:** Adding format wizards and chain cloning tools directly to `Studio` would re-introduce visual clutter.
   * **Resolution:** Keep `Studio` strictly limited to the 2D canvas and `N -> C` chain strip. Place `Project Name`, `Format Template`, `Add Chain`, `Clone Chain`, and `Heterodimer Fc Wizard` inside `Project`, and `CrossMab CH1<->CL Swap` + `Pairings` inside `Inspector & Domain`.

## Goals vs. Non-Goals

### Goals
* Minimal top bar (`File` dropdown menu) + 4 tabs (`Studio`, `Project`, `Inspector & Domain`, `Theme & Settings`), with all spacing sliders/buttons removed.
* Theme-aware Light/Dark contrast across `Domain Palette` cards, chain strip chips, panels, canvas, and chain labels.
* `-S-S-` disulfide bond rendering on both the interactive Flet canvas and transparent `.png` exports (`HINGE_DISULFIDE`, `CH1_CL`, `DISULFIDE_BOND`, and `CH3` `S354C`+`Y349C`).
* Multispecific design suite: `two_plus_one` (`2+1 T-Cell Engager`), `trispecific` (`Trispecific IgG-scFv`), `fab_igg` (`Fab-IgG`) format templates, `Heterodimer Fc Wizard` (`kih`, `kih_ss`, `charge_pair`), `clone_chain`, and `swap_crossmab_ch1_cl`.

### Non-Goals
* Breaking `.ab` `schema_version: 1` serialization/parsing compatibility.
* Re-adding spacing sliders or sidebars to the `Studio` tab.

## Alternatives Weighed

1. **Disulfide Modeling: `DISULFIDE_BOND` Pairing Type + `is_disulfide` Flag vs. New `.ab` Manifest Section:**
   * **Chosen:** Add `"DISULFIDE_BOND"` to `VALID_PAIRING_TYPES` and `is_disulfide: bool = False` on `SchematicPairingEdge`.
   * *Rejected by author: Adding a new top-level `disulfide_bonds` list in `.ab` manifests breaks `schema_version: 1` compatibility.*
2. **Dark-Mode Palette Fix: Dynamic `rebuild_palette()` vs. Fixed Mid-Gray:**
   * **Chosen:** Theme-aware `rebuild_palette()` bound to `_ui_palette(ctrl.theme_mode)`.
   * *Rejected by author: A fixed mid-gray background lowers text contrast in both Light and Dark modes.*

## Key Decisions & Staging for Next Steps

| Decision Area | Agreed Specification | Target Files |
| :--- | :--- | :--- |
| **Minimal 4-Tab UI & Contrast** | Top bar (`File` menu only) + `Studio` (canvas + chain strip only), `Project` (name, format, add/clone chain, heterodimer wizard), `Inspector & Domain` (dynamic palette, pairings, CrossMab swap, inspector, mutations), `Theme & Settings` (Light/Dark, presets, YAML). | `src/fc_engineer/app.py`, `tests/test_ui.py` |
| **Cysteine Bond (`-S-S-`)** | Add `"DISULFIDE_BOND"` to `VALID_PAIRING_TYPES`, `is_disulfide: bool = False` on `SchematicPairingEdge`, auto-detect `HINGE_DISULFIDE` / `CH1_CL` / `S354C`+`Y349C`, render solid `#F59E0B` bridge + `-S-S-` label in Flet and Pillow PNG. | `src/fc_engineer/designer.py`, `src/fc_engineer/app.py`, `src/fc_engineer/exporters.py` |
| **Multispecific Suite** | Add `two_plus_one`, `trispecific`, `fab_igg` to `VALID_FORMAT_TYPES` & `build_format_template()`; add `apply_heterodimer_pair()`, `clone_chain()`, `swap_crossmab_ch1_cl()` to `designer.py` / `AntibodyStudioController`. | `src/fc_engineer/designer.py`, `src/fc_engineer/app.py`, `tests/test_logic.py` |
