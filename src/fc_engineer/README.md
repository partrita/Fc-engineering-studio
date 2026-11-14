# fc_engineer

Modular antibody & bispecific domain engineering package combining **EU Numbering** constant-domain slicing, atomic cross-domain mutation routing, 2D schematic scene geometry, `.ab` GenBank-style persistence, FASTA/CSV/transparent PNG exporters, validated YAML configuration loading, and a Flet desktop GUI.

## Interface

### Entry Points & Package Metadata

| Symbol | Role | Source |
| :--- | :--- | :--- |
| `__version__` | Package version string constant (`"0.0.1"`) | `[__init__.py:__version__]` |
| `main` | Re-exported `fc_engineer.app.main` executed under `if __name__ == "__main__"` (`python -m fc_engineer`) | `[__main__.py]` |
| `main()` | Console script entry point (`fc-engineer`) launching `ft.run(...)` with `atexit` shutdown cleanup | `[app.py:main]` |

### Core EU Numbering & Constant-Domain Logic (`core.py`)

| Symbol | Signature / Value | Description | Source |
| :--- | :--- | :--- | :--- |
| `EU_START` | `118` | Starting EU residue number mapped to index `0` | `[core.py:EU_START]` |
| `IGG3_HINGE_EXTRA` | `47` | Extra hinge residue offset in bundled IgG3 `WT(P01860-1)` past EU gap | `[core.py:IGG3_HINGE_EXTRA]` |
| `IGG3_GAP_END` | `230` | Last EU residue number in unmapped IgG3 hinge repeat gap (`223..230`) | `[core.py:IGG3_GAP_END]` |
| `CONSTANT_DOMAIN_EU_RANGES` | `{"CH1": (118, 215), "Hinge": (216, 230), "CH2": (231, 340), "CH3": (341, 447)}` | Canonical EU ranges for heavy constant domains | `[core.py:CONSTANT_DOMAIN_EU_RANGES]` |
| `MAX_MUTATION_STRING_LENGTH`, `MAX_MUTATIONS` | `1000`, `50` | Maximum length and token count for mutation input strings | `[core.py:MAX_MUTATION_STRING_LENGTH]` |
| `MAX_BATCH_ITEMS`, `MAX_BATCH_LENGTH` | `20`, `2000` | Maximum batch combinations and input character length | `[core.py:MAX_BATCH_ITEMS]` |
| `get_residue_index` | `(pos: int, isotype: str) -> Optional[int]` | Maps EU residue number and isotype to a 0-based index, or `None` for gaps/unknown | `[core.py:get_residue_index]` |
| `get_eu_position` | `(index: int, isotype: str = "igg1") -> Optional[int]` | Inverse-maps 0-based index back to EU number (`None` on `<0`, unknown isotype, or `igg3` repeat `105..159`) | `[core.py:get_eu_position]` |
| `get_domain_slice_bounds` | `(domain_type: str, isotype: str, full_seq_len: int) -> Optional[Tuple[int, int]]` | Returns 0-based `[start, end)` bounds for `CH1`, `Hinge`, `CH2`, or `CH3` | `[core.py:get_domain_slice_bounds]` |
| `slice_constant_domain` | `(full_seq: str, domain_type: str, isotype: str) -> str` | Extracts an isolated constant domain sequence from a full constant sequence | `[core.py:slice_constant_domain]` |
| `get_domain_residue_index` | `(pos: int, domain_type: str, isotype: str) -> Optional[int]` | Maps EU number to a 0-based index within an isolated constant domain | `[core.py:get_domain_residue_index]` |
| `get_domain_eu_position` | `(local_idx: int, domain_type: str, isotype: str) -> Optional[int]` | Maps a 0-based domain index back to an EU position within the domain's EU range | `[core.py:get_domain_eu_position]` |
| `parse_mutation` | `(m_str: str) -> Tuple[str, int, str]` | Parses `<WT><pos><Mut>` token into `(wt_aa, pos, mut_aa)` | `[core.py:parse_mutation]` |
| `parse_mutation_list` | `(mutants_str: str) -> Tuple[List[Tuple[str, int, str]], List[str]]` | Validates and tokenizes slash/comma-separated EU mutations | `[core.py:parse_mutation_list]` |
| `apply_mutations` | `(sequence: str, mutants_str: str, isotype: str) -> Tuple[str, List[str]]` | Validates and applies EU mutations onto a full `CH1-Hinge-CH2-CH3` sequence | `[core.py:apply_mutations]` |
| `apply_domain_mutations` | `(domain_seq: str, mutants_str: str, domain_type: str, isotype: str) -> Tuple[str, List[str]]` | Validates and applies EU mutations onto an isolated constant domain sequence | `[core.py:apply_domain_mutations]` |
| `diff_sequences` | `(wt_seq: str, mut_seq: str, isotype: str, domain_type: Optional[str] = None) -> List[Tuple[int, str, str]]` | Returns `(eu_pos, wt_aa, mut_aa)` diffs for full or isolated constant sequences | `[core.py:diff_sequences]` |
| `generate_batch` | `(base_seq: str, batch_str: str, isotype: str) -> List[Tuple[str, str, List[str]]]` | Applies comma-separated mutation sets in batch | `[core.py:generate_batch]` |
| `analyze_developability` | `(sequence: str, isotype: str, domain_type: Optional[str] = None) -> List[Tuple[str, List[int]]]` | Scans for EU-positioned `N-X-S/T`, `NG`, `DG`, and `M` liabilities | `[core.py:analyze_developability]` |

### Modular Domain Assembly, Templates & 2D Scene Engine (`designer.py`)

| Symbol | Role / Signature | Source |
| :--- | :--- | :--- |
| `DomainBlock`, `PolypeptideChain`, `InterChainPairing`, `AntibodyProject` | Core antibody topology and state dataclasses | `[designer.py:AntibodyProject]` |
| `SchematicNode`, `SchematicEdge`, `SchematicPairingEdge`, `SchematicScene` | UI-agnostic 2D schematic geometry dataclasses (`SCENE_PADDING = 24.0`, `DOMAIN_STYLE_MAP`) | `[designer.py:SchematicScene]` |
| `build_format_template` | `(format_type: str, isotype: str = "igg1", allotype: Optional[str] = None) -> AntibodyProject` | `[designer.py:build_format_template]` |
| `build_preset_map` | `(project: Optional[AntibodyProject] = None) -> Dict[str, str]` | `[designer.py:build_preset_map]` |
| `get_effective_domain_mutations` | `(domain: DomainBlock, preset_map=None) -> Tuple[List[str], List[str]]` | `[designer.py:get_effective_domain_mutations]` |
| `resolve_domain_base_sequence`, `resolve_domain_sequence` | Resolve WT/custom/engineered amino acid sequences for any domain block | `[designer.py:resolve_domain_sequence]` |
| `assemble_chain_sequence` | `(chain: PolypeptideChain, preset_map=None) -> Tuple[str, List[Tuple[int, int, DomainBlock]], List[str]]` | `[designer.py:assemble_chain_sequence]` |
| `route_mutations_to_chain`, `remove_preset_from_chain` | Atomic Phase 1/Phase 2 cross-domain EU mutation routing and reversible preset removal | `[designer.py:route_mutations_to_chain]` |
| `remove_domain_from_chain`, `remove_chain_from_project`, `reflow_chain_coordinates` | Topology editing with orphaned pairing pruning and coordinate reflow | `[designer.py:reflow_chain_coordinates]` |
| `diff_chain_sequences`, `analyze_chain_developability`, `add_file_custom_mutation` | Chain-level constant-domain EU diffs, PTM scans, and `.ab`-scoped custom mutation registration | `[designer.py:add_file_custom_mutation]` |
| `build_schematic_scene` | `(project: AntibodyProject, preset_map=None) -> SchematicScene` | `[designer.py:build_schematic_scene]` |

### Configuration & Safe YAML Loading (`config.py`)

| Symbol | Signature / Value | Description | Source |
| :--- | :--- | :--- | :--- |
| `MAX_FILE_SIZE` | `1048576` (`1 MiB`) | Global file size cap shared by YAML and `.ab`/FASTA/PNG handlers | `[config.py:MAX_FILE_SIZE]` |
| `logger` | `logging.Logger` | Module logger configured with `logging.NullHandler()` | `[config.py:logger]` |
| `log_sanitized_error` | `(logger_obj: logging.Logger, msg: str, e: Exception) -> None` | Redacts package base directory paths (`"."`) in error logs and tracebacks | `[config.py:log_sanitized_error]` |
| `NoAliasSafeLoader` | `yaml.SafeLoader` subclass (`MAX_DEPTH = 50`) | Blocks `AliasEvent` bombs and deep nesting | `[config.py:NoAliasSafeLoader]` |
| `load_yaml_data`, `load_domain_templates` | Side-effect-free YAML loaders | Return validated constant sequences, presets, and VH/VL/CL/G4S templates | `[config.py:load_yaml_data]` |
| `SEQUENCES`, `COMMON_MUTATIONS`, `VARIABLE_DOMAINS`, `LIGHT_CONSTANT_DOMAINS`, `DEFAULT_LINKERS` | Module constants | Populated once at import via `_load_all_config()` | `[config.py:SEQUENCES]` |

### Persistence & Export Formatters (`exporters.py`)

| Symbol | Description | Source |
| :--- | :--- | :--- |
| `MAX_EXPORT_SEQUENCE_LENGTH`, `MAX_METADATA_LENGTH`, `MAX_PNG_DIMENSION` | Bounds constants (`10000`, `200`, `8192`) | `[exporters.py:MAX_PNG_DIMENSION]` |
| `sanitize_locus`, `format_genbank`, `format_csv_report` | Locus sanitization and GenBank/CSV string builders | `[exporters.py:format_genbank]` |
| `format_project_fasta`, `save_fasta_file` | Multi-chain FASTA formatter and validated file writer (`.fasta`/`.fa`) | `[exporters.py:save_fasta_file]` |
| `serialize_ab_project`, `parse_ab_project`, `save_ab_file`, `load_ab_file` | Multi-LOCUS `.ab` project serializer, strict parser, and bounded disk I/O | `[exporters.py:load_ab_file]` |
| `export_schematic_png`, `save_schematic_png` | Pillow transparent RGBA `.png` schematic renderer and bounded file writer | `[exporters.py:export_schematic_png]` |

### Flet Desktop GUI & Workspace Controller (`app.py`)

| Symbol | Key Capabilities | Source |
| :--- | :--- | :--- |
| `AntibodyStudioController` | Manages active `AntibodyProject`, selection, format templates, chain/domain CRUD & cross-chain drag moves, atomic cross-chain mutations, `.ab`/FASTA/PNG persistence, lock-protected 30s `threading.Timer` clipboard auto-clear, and `shutdown()` | `[app.py:AntibodyStudioController]` |
| `build_studio_ui`, `main` | Builds the 3-pane Flet desktop UI (`ft.DragTarget`, `ft.GestureDetector` + `flet.canvas.Canvas`, async `ft.FilePicker`, `window.prevent_close` shutdown hook) and starts `ft.run(...)` | `[app.py:build_studio_ui]` |

### Bundled YAML Data Assets

* [sequences.yaml](./data/sequences.yaml) — Reference constant sequences (`isotypes`: `igg1`–`igg4`), fallback variable regions (`variable_domains`: `VH_A`, `VL_A`, `VH_B`, `VL_B`), light constant domains (`light_constant_domains`: `CL_kappa`, `CL_lambda`), and linkers (`linkers`: `G4S_3`, `G4S_2`, `G4S_1`).
* [mutants.yaml](./data/mutants.yaml) — Read-only list of 46 `{label, value}` preset mutations under `common_mutations`.

## Invariants

* **YAML & Project File Defense-in-Depth (`[config.py:_load_yaml_dict]`, `[exporters.py:load_ab_file]`, `[exporters.py:parse_ab_project]`):** Requires `os.path.isfile`, verifies file size and bounded reads `<= 1 MiB` (`len(raw_bytes) <= MAX_FILE_SIZE` before UTF-8 decode), blocks YAML aliases and deep JSON nesting (`RecursionError`) / `NaN`/`Infinity`, validates finite non-boolean coordinates in `[-10000.0, 10000.0]`, enforces collection caps (`1..16` chains, `1..32` domains/chain, `<= 64` pairings, `<= 100` `file_custom_mutations`, `<= 50` `mutation_history` entries `<= 200` chars), rejects trailing content after `"//"`, and requires `ORIGIN == assembled_seq`.
* **EU Alignment, Slice Boundaries & Atomic Mutation Routing (`[core.py]`, `[designer.py]`):** Maps `CH1` (`118..215`), `Hinge` (`216..230`), `CH2` (`231..340`), and `CH3` (`341..447`); returns `None` from `get_eu_position` for `<0`, unknown isotypes, and `igg3` extended-hinge repeat indices `105..159`; enforces ASCII `[0-9]` regexes with `\Z` anchors; deduplicates identical tokens while rejecting conflicts; and aborts Phase 1 before mutating any domain or chain if any token fails.
* **Path Redaction & Export Safety (`[config.py:log_sanitized_error]`, `[exporters.py]`):** Strips package paths in logs, uses `os.path.basename` in disk error messages, sanitizes FASTA headers and GenBank qualifiers deterministically, and bounds transparent RGBA PNG renders to `<= MAX_PNG_DIMENSION` (`8192`).
* **Clipboard Security Policy (`[app.py:AntibodyStudioController]`):** Coordinates clipboard state under `_clipboard_lock`, pre-wipes prior copied text before overwriting, uses monotonic `_clipboard_gen` + `content_to_clear` checks on 30-second daemon `threading.Timer` callbacks, and purges clipboard contents on `.ab` load, format change, window close, and `atexit`.

## Side Effects

* **Import-Time Config Reads (`[config.py:_load_all_config]`):** Reads `./data/sequences.yaml` once and `./data/mutants.yaml` once at import, logging sanitized diagnostics on failure.
* **OS Clipboard, Background Timers & File I/O (`[app.py]`, `[exporters.py]`):** Reads/writes OS clipboard via `pyperclip`, spawns daemon `threading.Timer` threads, writes `.ab`/`.fasta`/`.png` files, reads `.ab` files, and reads `date.today()` for GenBank dates.

## Verification

```bash
uv run pytest
```
