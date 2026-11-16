# Stage 8: Intermediate Zoom Out Review Report

* **Target Document:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md) (`1..294`; `FeatureRange: 1..198`, `DetailedRange: 199..294`)
* **Prior Reports Reviewed:** [zoom_out_report.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/zoom_out_report.md), [goldfish_review_report.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/goldfish_review_report.md)
* **Codebase Audited:** [core.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/core.py), [config.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/config.py), [exporters.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/exporters.py), [app.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/app.py), [sequences.yaml](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/data/sequences.yaml)
* **Summary Tally:** **0 Blockers**, **4 Majors**, **3 Nits**

---

## 1. Stage 5 Finding Status & Stage 7 Complexity Impact

### Stage 5 Tracking Summary
| Prior ID | Original Title | Current Status | Notes |
| :--- | :--- | :--- | :--- |
| **IC-1** | Unsegmented Reverse EU Coordinate Mapping | **Resolved (with new duplication)** | Stage 7 `D-1` added `get_domain_eu_position`, `diff_chain_sequences`, and `analyze_chain_developability`, but duplicates `core.py`'s diff/PTM scanner (`ZO-S8-3`). |
| **IC-2** | Dual-Projection Coordinate Contention & Duplicate Schematic Engines | **Unresolved** | Still two independent schematic renderers (`app.py` Flet canvas vs. `exporters.py` Pillow PNG), undefined strip-to-canvas reflow, plus a new YAGNI `.png` preview toggle (`ZO-S8-1`). |
| **IC-3** | Destructive Preset Flattening | **Partially Resolved / Regressed** | Stage 7 `D-2` added `remove_preset_from_chain()`, but `DomainBlock.custom_mutations` still holds flattened effective mutations, contradicting preset recomputation (`ZO-S8-2`). |
| **ZO-1** | Duplicate `sequences.yaml` Loading & Validation in `config.py` | **Unresolved** | `load_domain_templates()` still re-opens and re-validates `sequences.yaml` separately (`ZO-S8-5`). |
| **ZO-2** | Triplicated Mutation Tokenization & Validation Across `core.py` / `designer.py` | **Unresolved (Compounded by `D-2`)** | Stage 7 `D-2` added a second domain-validation + deduplication/conflict pipeline inside `route_mutations_to_chain()` duplicating `apply_domain_mutations()` (`ZO-S8-4`). |
| **ZO-3** | Re-implementing GenBank Header and `ORIGIN` Formatting in `exporters.py` | **Unresolved** | `serialize_ab_project()` still does not reuse `sanitize_locus()`, `_genbank_date()`, or a shared `ORIGIN` formatter (`ZO-S8-6`). |
| **ZO-4** | Ambiguous Dual `(x, y)` Coordinates | **Resolved** | `Comp-1` clarified `PolypeptideChain.(x, y)` (anchor) vs. `DomainBlock.(x, y)` (absolute canvas coordinates). |

### Stage 7 Additions Assessment
* **Bounded Additions:** `D-3` (`.ab` bounds and ordinal `LOCUS` record check), `D-4` (daemon `threading.Timer`, `atexit`, session-switch clipboard wipe), `D-5` (intra-chain `VH_VL` pairing rule), `ST-1`, and `SY-1`..`SY-3` are concise, trace directly to security or domain requirements, and introduce no accidental complexity.
* **Accidental Complexity / Duplication Introduced:** `D-1` introduces duplicate PTM/diff logic across `core.py` and `designer.py` plus `config.py` coupling inside `core.py` (`ZO-S8-3`). `D-2` introduces an inconsistent `DomainBlock.custom_mutations` model (`ZO-S8-2`) and duplicates `apply_domain_mutations()` validation inside `route_mutations_to_chain()` (`ZO-S8-4`).

---

## 2. Stage 8 Review Findings

### `[Major]` ZO-S8-1: Unresolved `IC-2` (Dual Schematic Engines & Missing Reflow) + YAGNI `.png` Preview Toggle
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L83-L90), [lines 141–160](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L141-L160), [lines 244–265](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L244-L265); [zoom_out_report.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/zoom_out_report.md#L34-L42)
* **Rationale:**
  * The design maintains two decoupled 2D schematic layout and styling engines: Pillow drawing in `exporters.py` (`export_schematic_png`, lines `143–148`, `255`) and Flet `ft.Stack` + `ft.canvas.Canvas` drawing in `app.py` (lines `159`, `265`). Box geometry, connector anchor points, colors, and mutation badges will diverge between workspace and export.
  * Line `90` (`PolypeptideChain`) references *"newly added or reflowed domains on this chain"*, yet no reflow helper exists in `designer.py` (`231–246`) or `app.py` (`258–267`) to update `(domain.x, domain.y)` when domains are dropped, reordered, or deleted in the N-to-C Chain Slot Strip (`160`, `265`).
  * Line `265` includes a `"+ transparent .png preview toggle"` inside the interactive canvas, creating redundant preview state and live raster re-rendering on top of a canvas that already displays the 2D schematic alongside `Export PNG`.
* **Proposed Simplification:**
  * Extract a single UI-agnostic `build_schematic_scene(project)` helper in `designer.py` (or `exporters.py`) returning domain rects, colors, labels, badges, and connector/pairing segments, consumed by both Flet's canvas (`app.py`) and Pillow's `.png` renderer (`exporters.py`).
  * Specify a concise `reflow_chain_coordinates(chain)` helper triggered on N-to-C strip insertions/reordering for unpinned domains.
  * Remove the `"transparent .png preview toggle"` from `app.py` (line `265`).

### `[Major]` ZO-S8-2: Stage 7 `D-2` Contradiction on `DomainBlock.custom_mutations` & Over-Engineered Preset Removal
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L84-L86), [lines 114–115](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L114-L115), [lines 233–241](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L233-L241)
* **Rationale:**
  * Stage 7 `D-2` added `remove_preset_from_chain(chain, preset_label)` (lines `115`, `241`) to *"recompute effective domain mutations from remaining active presets and ad-hoc custom entries without stripping shared residues still required by other presets."*
  * However, `DomainBlock.custom_mutations` (line `85`) is defined as *"Effective point mutations active on this domain"* (holding flattened preset tokens like `L234A`, `L235A`), and `resolve_domain_sequence()` (line `238`) applies `domain.custom_mutations` directly. Because no field stores ad-hoc custom entries separately from preset-expanded mutations, `remove_preset_from_chain` cannot recover ad-hoc custom entries without parsing `mutation_history` (line `86`).
  * In addition, `remove_preset_from_chain(chain: PolypeptideChain, preset_label: str)` (line `241`) does not receive `project` or `project.file_custom_mutations`, so it cannot resolve file-scoped custom presets during recomputation.
* **Proposed Simplification:**
  * Keep `DomainBlock.applied_presets: List[str]` (preset labels) and `DomainBlock.custom_mutations: List[str]` (ad-hoc point mutations only) strictly disjoint; do not flatten presets into `custom_mutations`.
  * Derive effective mutations dynamically in `get_effective_domain_mutations(domain, preset_map)` (expanding `applied_presets` from a unified `COMMON_MUTATIONS + project.file_custom_mutations` map plus `custom_mutations`, deduplicated) inside `resolve_domain_sequence()`.
  * `remove_preset_from_chain` then removes `preset_label` from `domain.applied_presets` and appends an audit entry to `mutation_history`.

### `[Major]` ZO-S8-3: Stage 7 `D-1` Duplicates Diff & PTM Scanning Across `core.py` and `designer.py` and Breaks `core.py` Layer Purity
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L101-L119), [lines 222–229](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L222-L229), [line 242](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L242); [core.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/core.py#L38-L46), [lines 128–168](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/core.py#L128-L168)
* **Rationale:**
  * `core.py` already implements `diff_sequences()` (lines `38–46`) and the 4-category PTM liability scanner `analyze_developability()` (`N-X-S/T`, `NG`, `DG`, `M`, lines `128–168`). Re-implementing parallel diff and PTM scanning inside `designer.py`'s `diff_chain_sequences()` and `analyze_chain_developability()` (lines `119`, `242`) duplicates motif rules and category strings across modules.
  * `slice_constant_domain(domain_type, isotype, allotype)` is placed in `core.py` (line `225`), forcing `core.py`—which currently has no dependency on `config.py`—to import `SEQUENCES`.
  * Line `111` states that `get_domain_residue_index(pos, domain_type, isotype)` validates `0 <= local_idx < len(domain_base_seq)` against a resolved WT sequence, even though its signature in line `226` takes no sequence or allotype argument.
* **Proposed Simplification:**
  * Add an optional `domain_type: Optional[str] = None` parameter to `diff_sequences()` and `analyze_developability()` in `core.py` (routing index mapping through `get_domain_eu_position()` when `domain_type` is set, and `get_eu_position()` otherwise). `diff_chain_sequences(chain)` and `analyze_chain_developability(chain)` in `designer.py` then become thin constant-domain aggregators over `core.py`.
  * Keep `core.py` config-free: accept `full_seq: str` in `slice_constant_domain(full_seq, domain_type, isotype)` (with `SEQUENCES` lookup in `designer.py`), and keep `0 <= local_idx < len(domain_seq)` validation in `apply_domain_mutations()` rather than `get_domain_residue_index()` (aligning line `111` with line `226`).

### `[Major]` ZO-S8-4: Unresolved `ZO-2` Compounded by Stage 7 `D-2`: Duplicated Mutation Tokenization & Domain Validation
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L114-L115), [lines 222–228](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L222-L228), [lines 240–243](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L240-L243); [core.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/core.py#L60-L92); [zoom_out_report.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/zoom_out_report.md#L61-L64)
* **Rationale:**
  * `apply_mutations()` in `core.py` (lines `60–92`) already implements the 1000-character limit, character-set guard, `/` and `,` splitting, `MAX_MUTATIONS = 50`, and `parse_mutation()`. That pipeline remains duplicated across `apply_mutations`, `apply_domain_mutations`, `route_mutations_to_chain`, and `add_file_custom_mutation`.
  * Stage 7 `D-2` compounded this by duplicating domain-level validation (`0 <= local_idx < len(domain_base_seq)`, WT amino acid match, conflict rejection, deduplication) in both `apply_domain_mutations()` (`core.py`, line `228`) and `route_mutations_to_chain()` (`designer.py`, lines `114–115`, `240`).
* **Proposed Simplification:**
  * Extract `parse_mutation_list(mutants_str: str) -> Tuple[List[Tuple[str, int, str]], List[str]]` in `core.py` and reuse it across `apply_mutations`, `apply_domain_mutations`, `route_mutations_to_chain`, and `add_file_custom_mutation`.
  * Share a common index-mapped validation/application loop between `apply_mutations()` and `apply_domain_mutations()`, and perform Phase-1 atomic validation in `route_mutations_to_chain()` by grouping parsed tokens via `CONSTANT_DOMAIN_EU_RANGES` and dry-running `apply_domain_mutations()` on each touched domain before committing.

### `[Nit]` ZO-S8-5: Unresolved `ZO-1`: Duplicate `sequences.yaml` File-Loading and Validation Pipeline in `config.py`
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L215-L218); [config.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/config.py#L13-L112); [zoom_out_report.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/zoom_out_report.md#L56-L59)
* **Rationale:**
  * `config.py` repeats `os.path.isfile()`, `os.path.getsize() <= 1MB`, `.read(MAX_FILE_SIZE + 1)`, `yaml.load(..., Loader=NoAliasSafeLoader)`, and error sanitization twice in `load_yaml_data()` (lines `50–108`) without reusing its own `log_sanitized_error()` helper (`13–17`). Adding `load_domain_templates()` (line `217`) duplicates this security wrapper a third time and parses `sequences.yaml` twice at import time.
* **Proposed Simplification:**
  * Extract a single `_load_yaml_dict(path: str)` helper (reusing `log_sanitized_error()` and promoting `MAX_FILE_SIZE = 1 * 1024 * 1024` to a module constant reusable by `exporters.py`) and parse `isotypes`, `variable_domains`, `light_constant_domains`, and `linkers` during one read of `sequences.yaml`, while keeping `load_yaml_data() -> Tuple[Dict, List]`.

### `[Nit]` ZO-S8-6: Unresolved `ZO-3`: Re-implementing GenBank Header and `ORIGIN` Formatting in `exporters.py`
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L136-L138), [lines 250–252](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L250-L252); [exporters.py](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/exporters.py#L15-L65); [zoom_out_report.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/zoom_out_report.md#L66-L69)
* **Rationale:**
  * `exporters.py` already provides `sanitize_locus()` (`22–27`), `_genbank_date()` (`30–32`), `MAX_EXPORT_SEQUENCE_LENGTH` (`15`), and the 60-aa/10-residue-block `ORIGIN` builder inside `format_genbank()` (`35–65`), yet `serialize_ab_project()` does not reference reusing them.
* **Proposed Simplification:**
  * Extract `_format_origin_block(sequence: str)` (or parameterize `format_genbank()` with optional `comment_lines` / `feature_lines`) and reuse `sanitize_locus()`, `_genbank_date()`, and `MAX_EXPORT_SEQUENCE_LENGTH` inside `serialize_ab_project()`.

### `[Nit]` ZO-S8-7: Unused `validate_project_topology()` Heuristic & Broken IgG1-Hardcoded Default `allotype` on `build_format_template()`
* **Affected Scope:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L80), [lines 244–245](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L244-L245), [lines 258–282](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L258-L282); [sequences.yaml](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/src/fc_engineer/data/sequences.yaml#L2-L15)
* **Rationale:**
  * `validate_project_topology(project)` (line `245`) defines heuristic biological warnings (*"missing Knob/Hole complementarity on asymmetric heavy chains"*, which misfires on non-KiH heterodimer designs), yet has no caller in `app.py` (`150–166`, `258–267`) or `exporters.py` (`248–256`) and no unit test in Item #8 (`269–276`).
  * `build_format_template(format_type, isotype="igg1", allotype="WT(P01857-1)")` (line `244`) hardcodes IgG1's UniProt ID (`"WT(P01857-1)"`) as the default `allotype`, contradicting line `80` (*"defaults to the first WT key for that isotype"*): calling `build_format_template("igg", isotype="igg2")` will fail lookup in `SEQUENCES["igg2"]` (`"WT(P01859-1)"`).
* **Proposed Simplification:**
  * Remove `validate_project_topology()` until a UI warning panel or export guard consumes it (or limit it strictly to structural empty-chain/dangling-pairing checks if displayed in `app.py`).
  * Change `build_format_template()`'s default to `allotype: Optional[str] = None`, resolving to `next(iter(SEQUENCES[isotype]))` when `None`.
