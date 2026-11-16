# Stage 5: Initial Zoom Out & Inventive Catalyst Review Report

* **Target Document:** [docs/designs/flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md)
* **Reviewed Scope:** Full File (`1..291`; `FeatureRange: 1..198`, `DetailedRange: 199..291`)
* **Reviewers:** `egm-zoom-out`, `egm-inventive-catalyst`

---

## Summary of Findings

| ID | Reviewer | Severity | Title |
| :--- | :--- | :--- | :--- |
| **IC-1** | `egm-inventive-catalyst` | 🟡 Major | Unsegmented Reverse EU Coordinate Mapping Breaks Diffs and Developability Analysis on Modular Chains |
| **IC-2** | `egm-inventive-catalyst` | 🟡 Major | Dual-Projection Coordinate Contention and Duplicate Schematic Layout Engines |
| **IC-3** | `egm-inventive-catalyst` | 🟢 Nit | Destructive Preset Flattening Prevents Reversible Preset Toggling and Hides Same-Residue Collisions |
| **ZO-1** | `egm-zoom-out` | 🟢 Nit | Duplicate `sequences.yaml` File-Loading and Validation Pipeline in `config.py` |
| **ZO-2** | `egm-zoom-out` | 🟢 Nit | Triplicated Mutation-String Tokenization and Bounds Validation Across `core.py` and `designer.py` |
| **ZO-3** | `egm-zoom-out` | 🟢 Nit | Re-implementing GenBank Header and `ORIGIN` Formatting in `exporters.py` |
| **ZO-4** | `egm-zoom-out` | 🟢 Nit | Ambiguous Dual `(x, y)` Coordinates on Both `PolypeptideChain` and `DomainBlock` |

---

## Part 1: Inventive Catalyst Findings (`egm-inventive-catalyst`)

### IC-1: Unsegmented Reverse EU Coordinate Mapping Breaks Diffs and Developability Analysis on Modular Chains
* **Severity:** 🟡 Major
* **Affected Scope:** Technical Plan Section 3 (`Domain-Aware EU Numbering and Cross-Domain Mutation Routing`, lines 101–120), Section 7 (`Right Panel`, lines 161–165), and Detailed Implementation Items #4, #5, #6, #7 (lines 220–264).
* **False Compromise:**
  The design segments the *forward* lookup (`EU position -> 0-based domain index` via `get_domain_residue_index` and `apply_domain_mutations`), while leaving the *reverse* lookup (`0-based index -> EU position` in `get_eu_position`, `diff_sequences`, and `analyze_developability` inside `src/fc_engineer/core.py`) bound to `EU_START = 118` at index `0`. Calling `diff_sequences` or `analyze_developability` on an assembled `PolypeptideChain` (which prepends `VH`/`VL`/`G4S` or omits `CH1` in `scFv-Fc`) or on an isolated `Hinge`/`CH2`/`CH3` domain produces shifted EU numbers and discards PTM liabilities whenever `get_eu_position` returns `None`.
* **Proposed Resolution:**
  * Add `get_domain_eu_position(local_idx: int, domain_type: str, isotype: str) -> Optional[int]` in `src/fc_engineer/core.py` as the inverse of `get_domain_residue_index`, adding the domain's slice start offset before invoking `get_eu_position()`.
  * Implement chain-aware `diff_chain_sequences(chain: PolypeptideChain)` and `analyze_chain_developability(chain: PolypeptideChain)` in `src/fc_engineer/designer.py` that iterate across `domain_spans` from `assemble_chain_sequence()`, mapping constant-domain indices through `get_domain_eu_position()`.

### IC-2: Dual-Projection Coordinate Contention and Duplicate Schematic Layout Engines
* **Severity:** 🟡 Major
* **Affected Scope:** Technical Plan Section 2 (lines 76–100), Section 6 (lines 141–149), Section 7 (lines 158–160), and Detailed Implementation Items #5, #6, #7 (lines 230–264).
* **False Compromise:**
  * **Topological Edit vs. Manual 2D Placement:** Inserting or reordering domains in the N$\rightarrow$C Chain Strip either leaves existing `(domain.x, domain.y)` coordinates untouched (causing N$\rightarrow$C backbone connectors to zig-zag across overlapping boxes) or overwrites manual 2D canvas adjustments during reflow.
  * **Dual-Renderer Drift:** Domain box dimensions, color palettes, mutation badges, N$\rightarrow$C backbone connectors, and inter-chain pairing geometries are computed separately in `app.py` (`ft.Stack` + `ft.canvas.Canvas`) and in `exporters.py` (`export_schematic_png`), creating visual divergence between the workspace and exported `.png`.
* **Proposed Resolution:**
  * Extract a UI-agnostic `build_schematic_scene(project: AntibodyProject)` helper in `designer.py` (or `exporters.py`) that computes domain bounding boxes, RGBA colors, badge labels, and connector/pairing line coordinates once for both Flet's canvas and Pillow's `.png` renderer.
  * Provide a `reflow_chain_coordinates(chain: PolypeptideChain)` helper triggered on N$\rightarrow$C domain insertion/reordering unless a domain has been manually pinned/dragged.

### IC-3: Destructive Preset Flattening Prevents Reversible Preset Toggling and Hides Same-Residue Collisions
* **Severity:** 🟢 Nit
* **Affected Scope:** Technical Plan Section 2 (`DomainBlock`, lines 78–88), Section 3 (lines 113–116), and Section 7 (line 163, line 263).
* **False Compromise:**
  Flattening multi-residue presets (e.g., `GAALIELS` or `LALA`) directly into `domain.custom_mutations` means unchecking a preset checkbox can delete shared point mutations (`L234A`, `I332E`) still required by another active preset, and checking two presets that target the same EU residue (`LS` [`N434S`] and `N434A`) causes order-dependent collisions.
* **Proposed Resolution:**
  Keep `domain.applied_presets` separate from ad-hoc `domain.custom_mutations`, derive the effective point mutations by expanding `applied_presets` + `custom_mutations`, deduplicate identical mutations, and flag conflicting mutations at the same EU residue.

---

## Part 2: Zoom Out Findings (`egm-zoom-out`)

### ZO-1: Duplicate `sequences.yaml` File-Loading and Validation Pipeline in `config.py`
* **Severity:** 🟢 Nit
* **Affected Scope:** Detailed Implementation Item #3 (`./src/fc_engineer/config.py`, lines 214–218).
* **Simplicity Rationale:** `load_yaml_data()` in `src/fc_engineer/config.py` already opens `sequences.yaml`, checks `os.path.isfile()`, enforces the 1 MB file-size cap, and parses via `NoAliasSafeLoader`. Re-opening and re-parsing `sequences.yaml` in `load_domain_templates()` duplicates the security wrapper and parses the file twice. Extract a shared internal helper or parse `variable_domains`, `light_constant_domains`, and `linkers` during the single `sequences.yaml` load while keeping `load_yaml_data()` returning `(SEQUENCES, COMMON_MUTATIONS)`.

### ZO-2: Triplicated Mutation-String Tokenization and Bounds Validation Across `core.py` and `designer.py`
* **Severity:** 🟢 Nit
* **Affected Scope:** Technical Plan Section 3 (lines 112–116) & Detailed Implementation Items #4 and #5 (lines 227, 239–240).
* **Simplicity Rationale:** `apply_mutations()` in `src/fc_engineer/core.py` already enforces the 1000-character length cap, `r'[^a-zA-Z0-9/, ]'` character check, `/` and `,` splitting, `MAX_MUTATIONS = 50` cap, and `parse_mutation()` calls. Extract a shared helper in `core.py` (`parse_mutation_list(mutants_str: str)`) reused by `apply_mutations`, `apply_domain_mutations`, `route_mutations_to_chain`, and `add_file_custom_mutation`.

### ZO-3: Re-implementing GenBank Header and `ORIGIN` Formatting in `exporters.py`
* **Severity:** 🟢 Nit
* **Affected Scope:** Technical Plan Section 5 (lines 134–138) & Detailed Implementation Item #6 (`./src/fc_engineer/exporters.py`, lines 247–250).
* **Simplicity Rationale:** `src/fc_engineer/exporters.py` already provides `sanitize_locus()`, `_genbank_date()`, and the 60-aa-per-line `ORIGIN` block builder in `format_genbank()`. Explicitly reuse `sanitize_locus()`, `_genbank_date()`, and a shared `_format_origin_block(sequence: str)` helper inside `serialize_ab_project()`.

### ZO-4: Ambiguous Dual `(x, y)` Coordinates on Both `PolypeptideChain` and `DomainBlock`
* **Severity:** 🟢 Nit
* **Affected Scope:** Technical Plan Section 2 (lines 85, 92), Section 6 (line 145), and Section 7 (line 159).
* **Simplicity Rationale:** Both `DomainBlock` and `PolypeptideChain` define `x: float, y: float`, yet both the `.png` renderer and Flet canvas position domains using `(domain.x, domain.y)`. Clarify that `PolypeptideChain.(x, y)` is the chain anchor used to auto-position newly added or reflowed domains on that chain, while `DomainBlock.(x, y)` stores each domain's rendered canvas position.
