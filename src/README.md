# src

Source package root (`src-layout`) for **Fc Engineering Studio**, housing the installable `fc_engineer` Python package configured via `pyproject.toml`.

## Subdirectory Index

| Package | Description | Contract |
| :--- | :--- | :--- |
| `fc_engineer` | Modular antibody domain/bispecific engineering engine (`core.py`, `designer.py`), YAML loader (`config.py`), `.ab`/FASTA/GenBank/CSV/PNG exporters (`exporters.py`), and Flet desktop GUI (`app.py`) | [fc_engineer/README.md](./fc_engineer/README.md) |

## Interface

Rolled-up public surface exported from [fc_engineer](./fc_engineer/README.md):

* **Entry Points & Metadata:** `[fc_engineer/__init__.py:__version__]`, `[fc_engineer/__main__.py]`, `[fc_engineer/app.py:main]`, `[fc_engineer/app.py:AntibodyStudioController]`, `[fc_engineer/app.py:build_studio_ui]`
* **Domain Assembly & Schematic Engine (`designer.py`):** `DomainBlock`, `PolypeptideChain`, `InterChainPairing`, `AntibodyProject`, `SchematicNode`, `SchematicEdge`, `SchematicPairingEdge`, `SchematicScene`, `SCENE_PADDING`, `DOMAIN_STYLE_MAP`, `build_preset_map`, `get_effective_domain_mutations`, `resolve_domain_base_sequence`, `resolve_domain_sequence`, `assemble_chain_sequence`, `route_mutations_to_chain`, `remove_preset_from_chain`, `remove_domain_from_chain`, `remove_chain_from_project`, `diff_chain_sequences`, `analyze_chain_developability`, `add_file_custom_mutation`, `build_format_template`, `reflow_chain_coordinates`, `build_schematic_scene`
* **EU Numbering & Constant Domain Engine (`core.py`):** `EU_START`, `IGG3_HINGE_EXTRA`, `IGG3_GAP_END`, `CONSTANT_DOMAIN_EU_RANGES`, `MAX_MUTATION_STRING_LENGTH`, `MAX_MUTATIONS`, `MAX_BATCH_ITEMS`, `MAX_BATCH_LENGTH`, `get_residue_index`, `get_eu_position`, `get_domain_slice_bounds`, `slice_constant_domain`, `get_domain_residue_index`, `get_domain_eu_position`, `diff_sequences`, `parse_mutation`, `parse_mutation_list`, `apply_mutations`, `apply_domain_mutations`, `generate_batch`, `analyze_developability`
* **YAML Config & Diagnostics (`config.py` + `data/*.yaml`):** `MAX_FILE_SIZE`, `NoAliasSafeLoader`, `load_yaml_data`, `load_domain_templates`, `log_sanitized_error`, `SEQUENCES`, `COMMON_MUTATIONS`, `VARIABLE_DOMAINS`, `LIGHT_CONSTANT_DOMAINS`, `DEFAULT_LINKERS`, `logger` (reading `./fc_engineer/data/sequences.yaml` and `./fc_engineer/data/mutants.yaml`)
* **Persistence & Export Formatters (`exporters.py`):** `MAX_EXPORT_SEQUENCE_LENGTH`, `MAX_METADATA_LENGTH`, `MAX_PNG_DIMENSION`, `sanitize_locus`, `format_genbank`, `format_csv_report`, `format_project_fasta`, `save_fasta_file`, `serialize_ab_project`, `parse_ab_project`, `save_ab_file`, `load_ab_file`, `export_schematic_png`, `save_schematic_png`

## Invariants

* **Bounded Alias-Free YAML Parsing (`[fc_engineer/config.py:_load_yaml_dict]`, `[fc_engineer/config.py:NoAliasSafeLoader.compose_node]`):** Rejects non-regular/special files (directories, broken links, FIFOs, devices) via `os.path.isfile`, enforces a 1 MiB limit (`os.path.getsize` and bounded `f.read(MAX_FILE_SIZE + 1)` + UTF-8 byte length check), blocks YAML aliases (`AliasEvent`), caps nesting at `MAX_DEPTH = 50`, and validates exact schemas across `SEQUENCES`, `VARIABLE_DOMAINS`, `LIGHT_CONSTANT_DOMAINS`, `DEFAULT_LINKERS`, and `COMMON_MUTATIONS`.
* **EU Alignment, Hinge Gaps & Constant-Domain Slicing (`[fc_engineer/core.py]`):** Anchors `igg1` at `EU_START = 118`; excludes hinge gaps (`223..225` in `igg2`/`igg4`, `223..230` [`IGG3_GAP_END`] in `igg3`) in forward mapping; returns `None` from `get_eu_position` for negative indices, unrecognized isotypes, and `igg3` extended-hinge repeat indices `105..159`; and maps `CH1`, `Hinge`, `CH2`, and `CH3` slices locally without variable/linker shifts.
* **Strict ASCII Input, Atomic Routing & Export Persistence Bounds (`[fc_engineer/core.py]`, `[fc_engineer/designer.py]`, `[fc_engineer/exporters.py]`):** Enforces ASCII `[0-9]` regexes with `\Z` anchors (`len <= 10`), `[A-Z]+` base sequences, atomic Phase 1 dry-run rollback before Phase 2 mutation commits, `.ab` 1 MiB file size limits, collection caps (`1..16` chains, `1..32` domains/chain, `<= 64` pairings, `<= 100` custom muts, `<= 50` history entries of `<= 200` chars), exact `ORIGIN == assembled_seq` parity, and `MAX_PNG_DIMENSION = 8192` transparent RGBA rendering bounds.
* **Clipboard Auto-Clear Security (`[fc_engineer/app.py:AntibodyStudioController]`):** Coordinates clipboard writes, 30-second daemon `threading.Timer` callbacks, and monotonic `_clipboard_gen` tokens under `_clipboard_lock`, clearing matching clipboard data on timer expiration, format/session switch, window close, and `shutdown()`.

## Side Effects

* **Import-Time Disk I/O (`[fc_engineer/config.py:_load_all_config]`):** Imports read `./fc_engineer/data/sequences.yaml` once and `./fc_engineer/data/mutants.yaml` once, logging sanitized failures to `logger.error` and `sys.stderr`.
* **Project & Image File I/O (`[fc_engineer/exporters.py]`):** Reads `.ab` workspace files and writes `.ab`, `.fasta`/`.fa`, and `.png` files to user-selected paths, sanitizing error messages to `os.path.basename`.
* **OS Clipboard, Timers & Desktop Event Loop (`[fc_engineer/app.py]`):** Reads/writes OS clipboard via `pyperclip`, schedules daemon `threading.Timer` cleanup, and runs Flet's desktop GUI loop (`ft.run(...)`) with native `FilePicker` dialogs and `date.today()` GenBank timestamps.

## Verification

```bash
uv run pytest
```
