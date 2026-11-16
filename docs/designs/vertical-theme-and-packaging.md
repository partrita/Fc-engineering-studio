# Detailed Design: Cross-Platform Packaging & Customizable Vertical Schematic Theme System

* **Author:** fkt
* **Date:** 2026-11-16
* **Source Brief:** [concepts_brief.md](./concepts_brief.md)

## 1. Packaging & Cross-Platform Release (`pyproject.toml`, `.github/workflows/release.yml`, `src/fc_engineer/config.py`)

### 1.1 Python Version Floor
* `pyproject.toml` sets `requires-python = ">=3.12"` (lowered from `>=3.14`).
* Ensure all type annotations in `src/fc_engineer/*.py` remain compatible with Python 3.12+.

### 1.2 Frozen Binary Data Asset Lookup (`src/fc_engineer/config.py`)
* Update `_data_dir() -> str`:
  ```python
  def _data_dir() -> str:
      local = os.path.join(os.path.dirname(__file__), "data")
      if os.path.isdir(local):
          return local
      if getattr(sys, "frozen", False):
          meipass = getattr(sys, "_MEIPASS", "")
          if meipass:
              frozen_data = os.path.join(meipass, "fc_engineer", "data")
              if os.path.isdir(frozen_data):
                  return frozen_data
      return local
  ```

### 1.3 Cross-Platform `flet pack` Matrix (`.github/workflows/release.yml`)
* Keep the existing Python wheel/sdist build & test job, and add a `build-desktop` matrix job (`ubuntu-latest`, `macos-latest`, `windows-latest`) on Python `3.12`:
  * On `ubuntu-latest`: install `flet-desktop-light` so `libmpv` is excluded.
  * Run `flet pack src/fc_engineer/app.py` with:
    * `--name fc-engineer`
    * `--add-data "src/fc_engineer/data/*.yaml:fc_engineer/data"` (using `;` on Windows)
    * `--exclude-module` for: `tkinter`, `sqlite3`, `unittest`, `pydoc`, `doctest`, `curses`, `xmlrpc`, `pytest`, `bandit`, `PIL.ImageTk`, `PIL.ImageQt`
  * Smoke-test the packaged binary in CI before archiving.
  * Archive per-OS bundles:
    * `fc-engineer-v<version>-windows-x86_64.zip`
    * `fc-engineer-v<version>-macos-arm64.tar.gz`
    * `fc-engineer-v<version>-linux-x86_64.tar.gz`
  * Upload all archives, wheel, sdist, and `SHA256SUMS.txt` to GitHub Releases.

---

## 2. Theme & Vertical Layout Data Model (`src/fc_engineer/designer.py`)

### 2.1 Constants & Validation
* `VALID_SHAPES: Tuple[str, ...] = ("rect", "rounded_rect", "capsule", "ellipse", "hexagon")`
* `VALID_ORIENTATIONS: Tuple[str, ...] = ("vertical", "horizontal")`
* Bounds:
  * `border_width`: `0.5 <= val <= 8.0`
  * `width`: `32.0 <= val <= 240.0`
  * `height`: `16.0 <= val <= 120.0`
  * `domain_spacing`: `4.0 <= val <= 80.0` (default `18.0`)
  * `chain_spacing`: `60.0 <= val <= 260.0` (default `116.0`)
  * Colors (`fill_color`, `border_color`, `text_color`): `#RRGGBB` hex string.

### 2.2 Dataclasses
* `DomainVisualStyle`:
  * `fill_color: str`
  * `border_color: str`
  * `border_width: float`
  * `text_color: str`
  * `width: float`
  * `height: float`
  * `shape: str` (one of `VALID_SHAPES`)
* `SchematicTheme`:
  * `name: str`
  * `orientation: str = "vertical"` (`"vertical"` or `"horizontal"`)
  * `domain_spacing: float = 18.0`
  * `chain_spacing: float = 116.0`
  * `styles: Dict[str, DomainVisualStyle]` (covering all 8 `VALID_DOMAIN_TYPES`)
* `DEFAULT_SCHEMATIC_THEME`: `"Claude Code"` preset (`orientation="vertical"`, `domain_spacing=18.0`, `chain_spacing=116.0`):
  * `VH`: `fill="#D97757"`, `border="#F09A7E"`, `border_width=1.5`, `text="#FAF9F5"`, `width=52.0`, `height=88.0`, `shape="rounded_rect"`
  * `VL`: `fill="#E09F5A"`, `border="#F2BC82"`, `border_width=1.5`, `text="#141413"`, `width=52.0`, `height=88.0`, `shape="rounded_rect"`
  * `CH1`: `fill="#B85C3E"`, `border="#D97757"`, `border_width=1.5`, `text="#FAF9F5"`, `width=52.0`, `height=88.0`, `shape="rounded_rect"`
  * `CL`: `fill="#C48B71"`, `border="#DFB099"`, `border_width=1.5`, `text="#FAF9F5"`, `width=52.0`, `height=88.0`, `shape="rounded_rect"`
  * `Hinge`: `fill="#6E6A61"`, `border="#9E988C"`, `border_width=1.5`, `text="#FAF9F5"`, `width=38.0`, `height=52.0`, `shape="capsule"`
  * `CH2`: `fill="#5C8A72"`, `border="#84B39B"`, `border_width=1.5`, `text="#FAF9F5"`, `width=52.0`, `height=88.0`, `shape="rounded_rect"`
  * `CH3`: `fill="#3D6B58"`, `border="#639681"`, `border_width=1.5`, `text="#FAF9F5"`, `width=52.0`, `height=88.0`, `shape="rounded_rect"`
  * `G4S`: `fill="#524F48"`, `border="#7D786D"`, `border_width=1.5`, `text="#ECE9E0"`, `width=38.0`, `height=52.0`, `shape="capsule"`
* `SchematicNode`: append optional fields with defaults at the end for backward compatibility:
  * `border_color: str = "#0F172A"`
  * `border_width: float = 1.5`
  * `text_color: str = "#FFFFFF"`
  * `shape: str = "rounded_rect"`

### 2.3 Geometry & Reflow Functions
* `compute_hexagon_vertices(x: float, y: float, width: float, height: float) -> List[Tuple[float, float]]`:
  * If `height >= width` (vertical-pointed): `inset_y = height * 0.22`, `cx = x + width / 2.0`, vertices: `(cx, y)`, `(x + width, y + inset_y)`, `(x + width, y + height - inset_y)`, `(cx, y + height)`, `(x, y + height - inset_y)`, `(x, y + inset_y)`.
  * Otherwise (horizontal-flat-top): `inset_x = width * 0.22`, `cy = y + height / 2.0`, vertices: `(x + inset_x, y)`, `(x + width - inset_x, y)`, `(x + width, cy)`, `(x + width - inset_x, y + height)`, `(x + inset_x, y + height)`, `(x, cy)`.
* `reflow_chain_coordinates(chain: PolypeptideChain, theme: Optional[SchematicTheme] = None) -> None`:
  * Resolves `active_theme = theme or DEFAULT_SCHEMATIC_THEME`.
  * In `"vertical"` mode: `col_width = max((style.width for d in chain.domains), default=52.0)` and `cursor = chain.y`; for each domain with `not dom.pinned`, `dom.x = chain.x + (col_width - style.width) / 2.0`, `dom.y = cursor`, then `cursor += style.height + active_theme.domain_spacing` (advancing for all domains).
  * In `"horizontal"` mode: `row_height = max((style.height for d in chain.domains), default=44.0)` and `cursor = chain.x`; for each domain with `not dom.pinned`, `dom.x = cursor`, `dom.y = chain.y + (row_height - style.height) / 2.0`, then `cursor += style.width + active_theme.domain_spacing`.
* `reflow_project_coordinates(project: AntibodyProject, theme: Optional[SchematicTheme] = None) -> None`:
  * Resolves `active_theme = theme or DEFAULT_SCHEMATIC_THEME`, updates chain anchors (`chain.x = 40.0 + idx * active_theme.chain_spacing`, `chain.y = 40.0` in vertical mode; `chain.x = 40.0`, `chain.y = 40.0 + idx * active_theme.chain_spacing` in horizontal mode), and calls `reflow_chain_coordinates(chain, theme=active_theme)` on each chain (preserving `dom.pinned == True` domain coordinates).
* `build_schematic_scene(project: AntibodyProject, preset_map: Optional[Dict[str, str]] = None, theme: Optional[SchematicTheme] = None) -> SchematicScene`:
  * Uses `active_theme = theme or DEFAULT_SCHEMATIC_THEME`.
  * Populates `SchematicNode` with `style.width`, `style.height`, `style.fill_color`, `style.border_color`, `style.border_width`, `style.text_color`, `style.shape`.
  * Backbone `SchematicEdge` connects bottom-center `(prev.x + prev.width / 2, prev.y + prev.height)` to top-center `(next.x + next.width / 2, next.y)` in `"vertical"` mode, and right-center to left-center in `"horizontal"` mode.

---

## 3. Theme YAML Persistence (`src/fc_engineer/data/themes.yaml`, `src/fc_engineer/config.py`)

* `src/fc_engineer/data/themes.yaml` uses a top-level `presets:` mapping (`presets.<name>.orientation`, `domain_spacing`, `chain_spacing`, `styles.<domain_type>`) defining 5 built-in presets:
  1. `Claude Code` (default)
  2. `Classic Blueprint`
  3. `Monochrome Publication`
  4. `High Contrast`
  5. `Pastel Candy`
* Standalone theme `.yaml` files use top-level `schema_version: 1`, `name`, `orientation`, `domain_spacing`, `chain_spacing`, and `styles.<domain_type>`. Missing domain keys in `styles` fall back to `DEFAULT_SCHEMATIC_THEME.styles[dtype]`.
* `src/fc_engineer/config.py` adds:
  * `load_theme_presets() -> Dict[str, SchematicTheme]`
  * `parse_theme_yaml(content: str) -> SchematicTheme` (validates with `NoAliasSafeLoader`, `MAX_FILE_SIZE`, hex colors, bounds, `VALID_SHAPES`, `VALID_ORIENTATIONS`)
  * `serialize_theme_yaml(theme: SchematicTheme) -> str`
  * `load_theme_file(path: str) -> SchematicTheme`
  * `save_theme_file(path: str, theme: SchematicTheme) -> None`
* `.ab` project files (`exporters.py`) remain `schema_version: 1` with no theme data stored inside them.

---

## 4. Dual Vector Rendering & Studio UI (`src/fc_engineer/app.py`, `src/fc_engineer/exporters.py`)

### 4.1 Flet 2D Canvas (`src/fc_engineer/app.py`)
* `_build_canvas_shapes(scene: SchematicScene, theme: SchematicTheme) -> List[Any]`:
  * Draws pairing edges (in `"vertical"` mode, intra-chain `VH_VL` scFv arcs curve outward horizontally: `ctrl_x = max(pedge.x1, pedge.x2) + 34.0`, `ctrl_y = (pedge.y1 + pedge.y2) / 2.0`; in `"horizontal"` mode, `ctrl_x = (pedge.x1 + pedge.x2) / 2.0`, `ctrl_y = min(pedge.y1, pedge.y2) - 30.0`).
  * Draws N->C backbone edges.
  * Draws each `SchematicNode` shape on `canvas_bg` (`cv.Rect` with `border_radius=0.0` for `rect`, `8.0` for `rounded_rect`, `min(width, height) / 2.0` for `capsule`; `cv.Oval` for `ellipse`; `cv.Path` using `compute_hexagon_vertices()` for `hexagon`) with fill and stroke (`#FACC15` if selected, `#EF4444` if `error_flag`, else `node.border_color`).
* `rebuild_canvas_only()`:
  * Places a transparent `ft.Container` (`bgcolor=ft.Colors.TRANSPARENT`, `border=None`) inside `ft.GestureDetector` in `canvas_overlay` holding only the centered `ft.Text(node.label, color=node.text_color, ...)`.

### 4.2 Claude Code Warm Dark UI & Controls (`src/fc_engineer/app.py`)
* `build_studio_ui()` uses `ft.ThemeMode.DARK`, `page.bgcolor = "#141413"`, panels `#1E1E1C`, elevated inputs `#262521`, borders `#3D3B35`, primary text `#FAF9F5`, secondary text `#A6A299`, accent `#D97757`.
* Toolbar / Canvas Header controls:
  * Theme Preset `ft.Dropdown` (`Claude Code`, `Classic Blueprint`, `Monochrome Publication`, `High Contrast`, `Pastel Candy`)
  * `Customize Theme` button opening a modal `ft.AlertDialog` to edit `orientation`, `domain_spacing`, `chain_spacing`, and per-domain `fill_color`, `border_color`, `border_width`, `text_color`, `width`, `height`, `shape`
  * `Load Theme .yaml` and `Save Theme .yaml` buttons
  * Live `Domain Spacing` slider (`4..80`), `Chain Spacing` slider (`60..260`), and `Compact` (`6 / 84`), `Normal` (`18 / 116`), `Wide` (`36 / 156`) quick buttons.

### 4.3 Pillow Transparent PNG Exporter (`src/fc_engineer/exporters.py`)
* `export_schematic_png(project: AntibodyProject, scale: float = 2.0, theme: Optional[SchematicTheme] = None) -> bytes` and `save_schematic_png(..., theme: Optional[SchematicTheme] = None)`:
  * Passes `theme` to `build_schematic_scene()`.
  * Draws outward horizontal dashed arcs for intra-chain `VH_VL` pairings when `orientation == "vertical"`.
  * Draws `rect`, `rounded_rect`, `capsule` via `draw.rounded_rectangle`, `ellipse` via `draw.ellipse`, and `hexagon` via `draw.polygon(compute_hexagon_vertices(...))`, using `node.fill_color`, `node.border_color`, `node.border_width`, and `node.text_color` on a transparent RGBA image.
