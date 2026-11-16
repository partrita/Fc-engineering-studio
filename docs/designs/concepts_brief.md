# Cross-Platform Standalone Packaging & Customizable Vertical Schematic Theme System

* **Author:** fkt
* **Date:** 2026-11-16

## Overview
Fc Engineering Studio currently requires `uv` and Python 3.14 or newer, renders polypeptide chains horizontally left-to-right with fixed light-mode colors, and gives every domain a fixed size and rounded-rectangle shape. This feature makes four changes:
1. It lowers the minimum Python version in `pyproject.toml` to `>=3.12` and builds standalone desktop application archives for Windows, macOS, and Linux in GitHub Actions while stripping unused runtime libraries to reduce download size.
2. It rotates the default antibody schematic orientation 90 degrees so chains flow vertically top-to-bottom (N-terminus to C-terminus) and sit side-by-side left-to-right, with vertically proportioned domain blocks (`52x88` for standard domains, `38x52` for hinge and linker blocks) and theme-selectable vertical or horizontal layout.
3. It sets a Claude Code inspired warm dark theme (`"Claude Code"`) as the default preset for both the schematic canvas and the desktop application shell, and lets users customize fill color, border color, border width, text color, width, height, and 5 shapes (`rect`, `rounded_rect`, `capsule`, `ellipse`, `hexagon`) per domain type (`VH`, `VL`, `CH1`, `CL`, `Hinge`, `CH2`, `CH3`, `G4S`), saving and loading themes as standalone `.yaml` files without touching `.ab` project files.
4. It adds live adjustable inter-domain and inter-chain spacing controls (`domain_spacing` and `chain_spacing`, plus `Compact`, `Normal`, and `Wide` presets) so users can pack domains closer together or spread them farther apart on screen and in `.png` exports.

## Analysis

### Current Codebase State
* **Hardcoded Horizontal Geometry & Domain Styling (`src/fc_engineer/designer.py`):**
  * `DOMAIN_STYLE_MAP: Dict[str, Tuple[float, float, str]]` (lines 62–71) defines horizontal dimensions and blue/slate fill colors (`96.0 x 44.0` for `VH`, `VL`, `CH1`, `CL`, `CH2`, `CH3`; `64.0 x 36.0` for `Hinge` and `G4S`).
  * `SchematicNode` (lines 136–152) holds `domain_id`, `chain_id`, `chain_name`, `is_n_terminus`, `label`, `domain_type`, `x`, `y`, `width`, `height`, `fill_color`, `badge_text`, and `error_flag`.
  * `reflow_chain_coordinates()` (lines 638–646) positions unpinned domains left-to-right at `chain.x + idx * 116.0` and `chain.y`.
  * `build_format_template()` (lines 995–998) stacks chains top-to-bottom at `chain.x = 40.0` and `chain.y = 40.0 + idx * 96.0`.
  * `build_schematic_scene()` (lines 1082–1095) draws covalent backbone edges horizontally from the right edge of `prev_node` to the left edge of `next_node`.
* **2D Canvas Rendering & UI Shell (`src/fc_engineer/app.py`):**
  * `_build_canvas_shapes()` (lines 594–654) draws N->C backbone lines (`cv.Line`) and inter-chain pairing lines (`cv.Line`, upward `cv.Path.QuadraticTo` for intra-chain scFv pairs) onto `canvas_bg`.
  * `rebuild_canvas_only()` (lines 752–860) renders each `SchematicNode` as an `ft.GestureDetector` wrapping an `ft.Container` (lines 819–832) with fixed `border_radius=6`, `#FFFFFF` label text, and `#0F172A` / `#EF4444` / `#FACC15` borders.
  * `build_studio_ui()` uses `ft.ThemeMode.LIGHT` (line 664) and `#F8FAFC` for the 2D canvas background (line 1552).
* **Pillow PNG Rendering (`src/fc_engineer/exporters.py`):**
  * `export_schematic_png()` (lines 878–993) draws rounded rectangles (`radius=6.0 * scale`), upward dashed arcs (`_draw_dashed_arc`) for intra-chain `VH_VL` pairings, and white text on a transparent RGBA image.
* **Configuration & Frozen Binary Path Resolution (`src/fc_engineer/config.py`):**
  * `_data_dir()` (lines 59–60) returns `os.path.join(os.path.dirname(__file__), "data")`.
  * `NoAliasSafeLoader` (lines 34–56) and `MAX_FILE_SIZE = 1 * 1024 * 1024` (line 14) enforce safe YAML loading without aliases, deep nesting, or files larger than 1 MiB.
* **Packaging & Release Workflow (`pyproject.toml`, `.github/workflows/release.yml`):**
  * `pyproject.toml` sets `requires-python = ">=3.14"` (line 10).
  * `.github/workflows/release.yml` installs Python 3.14 on `ubuntu-latest` and builds only Python wheel/sdist packages and `SHA256SUMS.txt`.

### Key Architectural Contradictions & Resolutions

1. **90-Degree Vertical Layout vs. Intra-Chain ScFv Arcs, Chain Anchors, and `DOMAIN_STYLE_MAP` Dimensions**
   * **Problem:** Rotating polypeptide chains by 90 degrees reverses both chain progression (top-to-bottom along Y instead of left-to-right along X) and chain stacking (columns left-to-right along X instead of rows top-to-bottom along Y). In vertical mode, connecting domains at left/right edges crosses adjacent chains, and curving scFv arcs upward crosses consecutive domains along the backbone.
   * **Resolution:**
     * Add `orientation: Literal["vertical", "horizontal"] = "vertical"` to `SchematicTheme` (with `"vertical"` as the default).
     * In `"vertical"` mode, `build_format_template()`, `AntibodyStudioController.add_chain()`, and reflow helpers arrange chains side-by-side (`chain.x = 40.0 + idx * chain_spacing`, `chain.y = 40.0`). Consecutive unpinned domains progress downward (`N` at top -> `C` at bottom), centered horizontally around each chain's column axis (`dom.x = chain.x + (col_width - dom_width) / 2.0`) so narrower `Hinge` (`38.0`) and `G4S` (`38.0`) blocks align along the centerline of `52.0`-wide domains (`VH`, `VL`, `CH1`, `CL`, `CH2`, `CH3`).
     * `build_schematic_scene()` connects consecutive domains bottom-center to top-center when `orientation == "vertical"`, and right-center to left-center when `orientation == "horizontal"`.
     * Intra-chain `VH_VL` scFv arcs in both `_build_canvas_shapes()` and `export_schematic_png()` curve outward to the side (`ctrl_x = max(x1, x2) + 34.0`) in vertical mode instead of upward (`ctrl_y = min(y1, y2) - 30.0`), while horizontal mode keeps upward arcs.
     * Existing tests such as `test_reflow_pruning_and_schematic_scene()` in `tests/test_logic.py` that assert fixed horizontal dimensions will be updated to verify vertical/horizontal bounds and custom theme dimensions explicitly.
2. **Flet `ft.Container` Hit-Testing vs. Non-Rectangular Shapes (`ellipse`, `hexagon`, `capsule`)**
   * **Problem:** `ft.Container` in `src/fc_engineer/app.py` only draws rectangles with uniform corner radii (`border_radius`). It cannot render a 6-vertex hexagon or a true non-circular oval (`ellipse`), so relying on `ft.Container` alone causes the interactive canvas and `.png` exports to diverge.
   * **Resolution:** Move domain shape drawing onto `canvas_bg` (`flet.canvas.Canvas` inside `_build_canvas_shapes()`) right after the pairing and backbone lines, using `cv.Rect` (for `rect`, `rounded_rect`, `capsule`), `cv.Oval` (for `ellipse`), and `cv.Path` (for `hexagon`). Keep the `ft.GestureDetector` in `canvas_overlay` with a transparent hit-box `ft.Container` (`bgcolor=ft.Colors.TRANSPARENT`, `border=None`) containing only the centered label `ft.Text` (colored by `node.text_color`). Share a single geometry helper `compute_hexagon_vertices(x, y, width, height)` in `src/fc_engineer/designer.py` between `app.py` and `exporters.py`.
3. **Variable Domain Sizes & Live Spacing Controls vs. Fixed `116.0` Coordinate Step**
   * **Problem:** Using a fixed step (`idx * 116.0`) either leaves large empty gaps when domains are compact or causes overlap when a user increases domain `height` (in vertical mode) or `width` (in horizontal mode) in a custom theme.
   * **Resolution:** Add `domain_spacing: float = 18.0` (bounded `4.0..80.0`, edge-to-edge gap between adjacent domains on a chain) and `chain_spacing: float = 116.0` (bounded `60.0..260.0`, step between chain columns/rows) to `SchematicTheme`. `reflow_chain_coordinates(chain, theme=...)` accumulates `cursor` using each domain's resolved height (vertical) or width (horizontal) plus `theme.domain_spacing`, and a project-level helper `reflow_project_coordinates(project, theme)` updates unpinned chain anchors using `theme.chain_spacing` whenever the user changes spacing or theme orientation.
4. **Standalone Binary Size Minimization vs. Frozen Data Asset Resolution (`sys._MEIPASS`)**
   * **Problem:** `flet pack` (PyInstaller) bundles Flutter desktop binaries and Python packages; on Linux `flet-desktop` pulls in `libmpv` audio/video libraries, and frozen binaries resolve package paths relative to `sys._MEIPASS`.
   * **Resolution:**
     * Update `_data_dir()` in `src/fc_engineer/config.py` to check `os.path.join(os.path.dirname(__file__), "data")` first, and fall back to `os.path.join(getattr(sys, "_MEIPASS", ""), "fc_engineer", "data")` when running under PyInstaller (`getattr(sys, "frozen", False)`).
     * Pass `--add-data` for `src/fc_engineer/data/*.yaml` and curated `--exclude-module` flags (`tkinter`, `sqlite3`, `unittest`, `pydoc`, `doctest`, `curses`, `xmlrpc`, `pytest`, `bandit`, `PIL.ImageTk`, `PIL.ImageQt`) to `flet pack`.
     * Install `flet-desktop-light` on Linux CI runners before running `flet pack` to drop `libmpv`, and run a headless smoke check on the packaged executable before archiving.

## Goals vs. Non-Goals

### Goals
* Lower `requires-python` in `pyproject.toml` from `>=3.14` to `>=3.12` and keep all unit and headless UI tests passing on Python 3.12+.
* Extend `.github/workflows/release.yml` with a cross-platform matrix (`ubuntu-latest`, `macos-latest`, `windows-latest`) that builds standalone executables via `flet pack` (`flet-desktop-light` on Linux, curated `--exclude-module` flags), smoke-tests each build, and uploads OS archives (`fc-engineer-v<version>-windows-x86_64.zip`, `fc-engineer-v<version>-macos-arm64.tar.gz`, `fc-engineer-v<version>-linux-x86_64.tar.gz`) alongside wheel, sdist, and `SHA256SUMS.txt`.
* Default to 90-degree rotated vertical domain and chain layout (`orientation: "vertical"`) across templates, reflow, Flet canvas rendering, and Pillow PNG export, while supporting `"horizontal"` layout selection in theme settings.
* Ship `"Claude Code"` as the default built-in theme preset in `src/fc_engineer/data/themes.yaml` and style `build_studio_ui()` in `src/fc_engineer/app.py` with the warm dark Claude Code palette:
  * Studio UI shell (`ft.ThemeMode.DARK`): background `#141413`, cards `#1E1E1C`, elevated inputs `#262521`, borders `#3D3B35`, primary text `#FAF9F5`, secondary text `#A6A299`, accent `#D97757`, canvas background `#141413` (while keeping exported `.png` backgrounds transparent RGBA).
  * `"Claude Code"` domain styles (`orientation: "vertical"`, `domain_spacing: 18.0`, `chain_spacing: 116.0`):
    * `VH`: fill `#D97757`, border `#F09A7E`, text `#FAF9F5`, width `52`, height `88`, shape `rounded_rect`
    * `VL`: fill `#E09F5A`, border `#F2BC82`, text `#141413`, width `52`, height `88`, shape `rounded_rect`
    * `CH1`: fill `#B85C3E`, border `#D97757`, text `#FAF9F5`, width `52`, height `88`, shape `rounded_rect`
    * `CL`: fill `#C48B71`, border `#DFB099`, text `#FAF9F5`, width `52`, height `88`, shape `rounded_rect`
    * `Hinge`: fill `#6E6A61`, border `#9E988C`, text `#FAF9F5`, width `38`, height `52`, shape `capsule`
    * `CH2`: fill `#5C8A72`, border `#84B39B`, text `#FAF9F5`, width `52`, height `88`, shape `rounded_rect`
    * `CH3`: fill `#3D6B58`, border `#639681`, text `#FAF9F5`, width `52`, height `88`, shape `rounded_rect`
    * `G4S`: fill `#524F48`, border `#7D786D`, text `#ECE9E0`, width `38`, height `52`, shape `capsule`
* Ship additional built-in presets in `src/fc_engineer/data/themes.yaml` (`Classic Blueprint`, `Monochrome Publication`, `High Contrast`, `Pastel Candy`).
* Support full per-`domain_type` (`VH`, `VL`, `CH1`, `CL`, `Hinge`, `CH2`, `CH3`, `G4S`) customization of `fill_color`, `border_color`, `border_width` (`0.5..8.0`), `text_color`, `width` (`32.0..240.0`), `height` (`16.0..120.0`), and `shape` (`rect`, `rounded_rect`, `capsule`, `ellipse`, `hexagon`), plus `orientation`, `domain_spacing` (`4.0..80.0`), and `chain_spacing` (`60.0..260.0`).
* Add live spacing controls on the 2D schematic view (`Domain Spacing` slider, `Chain Spacing` slider, and `Compact` / `Normal` / `Wide` quick preset buttons) so users can tighten or expand domain distances immediately and export matching `.png` schematics.
* Render all 5 shapes identically in `src/fc_engineer/app.py` (`flet.canvas.Canvas`) and `src/fc_engineer/exporters.py` (`ImageDraw`), keeping error/selection border states clearly visible and printing readable chain headers and mutation badges above/beside vertical chains.
* Support standalone `.yaml` theme save and load via `src/fc_engineer/config.py` using `NoAliasSafeLoader` and `MAX_FILE_SIZE` (1 MiB), keeping `.ab` project files (`schema_version: 1`) untouched.

### Non-Goals
* Storing theme data inside `.ab` project files or incrementing `.ab` `schema_version` beyond `1`.
* Per-instance styling overrides (for example, giving two `CH2` domains in one project different colors or shapes).
* Paid Apple Developer ID notarization or Windows EV Authenticode code signing in GitHub Actions.
* Bundling custom font files (`.ttf` / `.otf`) or user-uploaded fonts.
* Rewriting chain/domain `pinned` manual drag positions when a user changes `domain_spacing` or `chain_spacing` (only unpinned chains/domains reflow unless the user resets pinning).

## Alternatives Weighed

### 1. Cross-Platform Packaging Tooling
* **Option A (Selected): `flet pack` (PyInstaller wrapper) + `flet-desktop-light` on Linux + curated `--exclude-module` flags.**
  Packages Flet desktop binaries across Windows, macOS, and Linux with minimal CI complexity while removing `libmpv` on Linux and unused Python stdlib/dev modules.
* **Option B: `flet build` (Flutter SDK compilation) in CI.**
  *Rejected by author: Installing the Flutter SDK, CMake, Ninja, GTK3 dev headers, and Xcode build steps increases CI runtime from ~3 minutes to 15–25 minutes per release without shrinking the Flutter runtime binary.*
* **Option C: Nuitka C compilation.**
  *Rejected by author: Flet desktop releases still embed the precompiled Flutter desktop engine binary; compiling Python to C with Nuitka adds long build times with negligible size savings over excluding unused Python modules.*

### 2. Flet Canvas Rendering of 5 Domain Shapes
* **Option A (Selected): Hybrid `flet.canvas.Canvas` vector underlay + transparent `ft.GestureDetector` hit-box overlay.**
  Draws all 5 vector shapes (`rect`, `rounded_rect`, `capsule`, `ellipse`, `hexagon`) in `canvas_bg` and overlays transparent `ft.GestureDetector` controls for tap selection, drag updates, and labels.
* **Option B: Use `ft.Container` `border_radius` only.**
  *Rejected by author: `ft.Container` cannot draw 6-sided hexagons or true non-circular ellipses.*
* **Option C: Rasterize per-domain Pillow images to base64 strings inside `ft.Image`.**
  *Rejected by author: Re-encoding PNGs during interactive `on_pan_update` drag events causes lag and scaling blur.*

### 3. Adjustable Spacing Mechanism
* **Option A (Selected): Edge-to-edge `domain_spacing` + `chain_spacing` inside `SchematicTheme` with live sliders and `Compact`/`Normal`/`Wide` buttons.**
  Handles arbitrary vertical or horizontal domain heights/widths without overlaps and saves spacing preferences alongside the `.yaml` theme.
* **Option B: Apply a uniform zoom/scale multiplier to domain boxes.**
  *Rejected by author: Scaling domain boxes shrinks label readability when the user only wants closer gaps (`Compact`) or wider breathing room (`Wide`) between domains.*

## Key Decisions & Staging for Next Steps

| Decision Area | Agreed Specification | Target Files |
| :--- | :--- | :--- |
| **Python Version Floor** | Lower `requires-python` to `">=3.12"` and verify tests pass on Python 3.12–3.14. | `pyproject.toml`, `.github/workflows/release.yml` |
| **Vertical Orientation & Geometry** | Default `orientation = "vertical"` (`52x88` standard blocks, `38x52` `Hinge`/`G4S` centered on chain columns), vertical N->C backbone edges, outward horizontal scFv arcs in vertical mode, and shared `compute_hexagon_vertices()`. | `src/fc_engineer/designer.py`, `src/fc_engineer/app.py`, `src/fc_engineer/exporters.py` |
| **Theme & Spacing Model** | Add `VALID_SHAPES`, `VALID_ORIENTATIONS`, `DomainVisualStyle`, `SchematicTheme` (with `orientation`, `domain_spacing`, `chain_spacing`), `DEFAULT_SCHEMATIC_THEME`, `reflow_project_coordinates()`, and optional styling fields on `SchematicNode`. | `src/fc_engineer/designer.py` |
| **Claude Code Default Theme & YAML Presets** | Add `src/fc_engineer/data/themes.yaml` with `"Claude Code"` as default plus `Classic Blueprint`, `Monochrome Publication`, `High Contrast`, `Pastel Candy`; implement `load_theme_presets()`, `parse_theme_yaml()`, `serialize_theme_yaml()`, `load_theme_file()`, `save_theme_file()` via `NoAliasSafeLoader` and `MAX_FILE_SIZE`. | `src/fc_engineer/data/themes.yaml`, `src/fc_engineer/config.py` |
| **Frozen Binary Data Lookup** | Update `_data_dir()` to fall back to `sys._MEIPASS/fc_engineer/data` when `getattr(sys, "frozen", False)` is true. | `src/fc_engineer/config.py` |
| **Warm Dark App UI & Theme/Spacing Controls** | Apply Claude Code warm dark palette across `build_studio_ui()`; render 5 shapes on `canvas_bg`; add preset selector, `Customize Theme` dialog, `Load/Save Theme .yaml`, and live `Domain Spacing` / `Chain Spacing` + `Compact`/`Normal`/`Wide` controls. | `src/fc_engineer/app.py` |
| **Transparent PNG Shape & Orientation Parity** | Render rotated orientation, spacing, all 5 shapes, custom colors/borders/text colors, and outward scFv arcs in `export_schematic_png()` while keeping RGBA transparency. | `src/fc_engineer/exporters.py` |
| **Cross-Platform CI Matrix** | Build standalone archives on `ubuntu-latest`, `macos-latest`, and `windows-latest` using `flet pack`, `flet-desktop-light` (Linux), `--exclude-module`, smoke checks, and upload with `SHA256SUMS.txt`. | `.github/workflows/release.yml`, `pyproject.toml` |
