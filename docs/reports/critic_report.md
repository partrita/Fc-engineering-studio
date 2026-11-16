# Stage 7 Critic Review Report

* **Target Document:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md)
* **Reviewed Scope:** Feature Design (`FeatureRange`: lines 1–198)
* **Reviewers:** `egm-design`, `egm-structural`, `egm-style`
* **Summary Tally:** **1 Blocker**, **3 Majors**, **5 Nits**

---

## Design Critic Findings (`egm-design`)

### `[Blocker]` D-1: Broken EU Numbering in Per-Chain Diff and Developability Analysis for Modular Polypeptide Chains
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L101-L120), [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L164) (`src/fc_engineer/core.py`, `src/fc_engineer/designer.py`, `src/fc_engineer/app.py`)
* **Issue:** Section 3 keeps `diff_sequences()` and `analyze_developability()` intact, and Section 7 displays live per-chain EU diff vs. WT and developability (`analyze_developability`). However, `core.py`'s `get_eu_position(idx, isotype)` assumes index `0` is `EU 118` (`CH1`). Calling `diff_sequences()` or `analyze_developability()` on an assembled `PolypeptideChain` (which prepends `VH`/`VL`/`G4S` or omits/swaps `CH1`) or on an isolated `Hinge`/`CH2`/`CH3` domain shifts all EU numbers and misreports PTM liabilities such as `N297`, `M252`, and `M428`.
* **Recommendation:** Add `get_domain_eu_position(local_idx: int, domain_type: str, isotype: str) -> Optional[int]` in `src/fc_engineer/core.py` (adding the domain's slice start offset before calling `get_eu_position`), and define `diff_chain_sequences(chain: PolypeptideChain)` and `analyze_chain_developability(chain: PolypeptideChain)` in `src/fc_engineer/designer.py` so EU diffs and PTM liabilities are evaluated per constant domain (`CH1`, `Hinge`, `CH2`, `CH3`) using domain-local indices rather than whole-chain index offsets.

### `[Major]` D-2: Non-Atomic Cross-Domain Mutation Routing, Duplicate/Override Collisions, and Allotype Length Bounds
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L86-L88), [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L113-L116), [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L163) (`src/fc_engineer/core.py`, `src/fc_engineer/designer.py`)
* **Issue:**
  1. Applying tokens one-by-one mutates earlier domains (`CH2`) even if a later token (`CH3`, such as `H435A` on `igg3` which has `R435`, or a missing `CH3` domain) fails validation.
  2. Appending overlapping presets (such as `LALA` + `PGLALA`, both containing `L234A`/`L235A`) to `domain.custom_mutations` causes `apply_domain_mutations` to fail on the duplicate token because position 234 is already `'A'`; unchecking/removing presets is also unspecified.
  3. Allotypes such as `cetuximab` (`igg1`, 329 aa) omit C-terminal `K447`, making `get_domain_residue_index(447, "CH3", "igg1") == 106 == len(ch3_seq)`.
* **Recommendation:** Specify two-phase atomic validation in `route_mutations_to_chain` (validate all tokens and target domains before mutating state), enforce `0 <= local_idx < len(domain_base_seq)` against the resolved domain WT sequence, deduplicate identical mutation tokens across overlapping presets, flag conflicting amino acid substitutions at the same EU position, and define preset removal/toggling.

### `[Major]` D-3: `.ab` Multi-LOCUS Record Count Validation and Empty-Chain Edge Case
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L132-L140) (`src/fc_engineer/exporters.py`)
* **Issue:** Line 139 specifies max bounds (`max 16 chains`, `max 32 domains per chain`) without lower bounds (`0` chains or `0` domains yields no valid `1..N` `LOCUS` sequence record). Additionally, truncating `LOCUS` names can cause name collisions across chains, and a modified `.ab` file with fewer `LOCUS ... //` records than `len(project.chains)` could skip full `ORIGIN` verification if not checked against record count.
* **Recommendation:** Enforce `1 <= len(project.chains) <= 16` and `1 <= len(chain.domains) <= 32`, parse `LOCUS ... //` records in ordinal sequence, require `len(locus_records) == len(project.chains)`, and verify record `i`'s `ORIGIN` sequence against `project.chains[i]`.

### `[Major]` D-4: Non-Daemon `threading.Timer` Process Hang and Cross-Session Clipboard Teardown
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L165) (`src/fc_engineer/app.py`)
* **Issue:** Python's `threading.Timer` defaults to `daemon = False`, which keeps the Python process or test runner alive for up to 30 seconds if window-close teardown is bypassed. In addition, session switches (loading an `.ab` file or switching format templates) should wipe previously copied sequence data from the OS clipboard.
* **Recommendation:** Set `self._clipboard_timer.daemon = True`, register `atexit.register(controller.shutdown)` alongside the Flet window-close teardown hook, and execute the `try...finally` clipboard wipe when loading `.ab` files or switching format templates.

### `[Nit]` D-5: Intra-Chain Pairing Rule in `InterChainPairing` and Symmetric Heavy-Chain Mutation UX
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L93-L94), [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L125-L131)
* **Issue & Recommendation:** Explicitly note in Section 2 that `chain_a_id == chain_b_id` is valid for intra-chain `VH_VL` scFv pairings (`IgG-scFv`, `scFv-Fc`, `BiTE`) provided `(chain_a_id, domain_a_id) != (chain_b_id, domain_b_id)`, and clarify whether Fc preset selection in the Inspector applies to the currently selected chain or optionally syncs across symmetric heavy chains.

---

## Structural Findings (`egm-structural`)

### `[Nit]` ST-1: Inconsistent Directional Notation (`N-to-C` vs `N$\rightarrow$C`)
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L147), [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L159-L160) vs. [L17](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L17), [L24](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L24), [L91](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L91), [L177](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L177)
* **Issue & Recommendation:** Lines 147, 159, and 160 use LaTeX math (`N$\rightarrow$C`), whereas lines 17, 24, 91, and 177 use `N-to-C` (and line 53 uses `N->C`). Standardize prose instances on `N-to-C` so Markdown renderers without MathJax do not display raw `$` tokens.

---

## Style Findings (`egm-style`)

### `[Nit]` SY-1: Uncalibrated Claim on Existing Clipboard Implementation
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L165)
* **Issue & Recommendation:** Line 165 says "Retains the 30-second auto-clear timer (`threading.Timer`) and window-close teardown (`clear_clipboard`)", which implies the existing Textual app already uses `threading.Timer`. Clarify that the design retains the 30-second auto-clear policy and `clear_clipboard` teardown behavior while re-implementing them using `threading.Timer` and Flet window-close hooks instead of Textual's `set_timer` and `on_unmount`.

### `[Nit]` SY-2: Cross-Section Repetition Between `Overview` and `Goals`
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L15-L29)
* **Issue & Recommendation:** `## Overview` repeats the full enumerated list of domains, formats, and `.ab` fields nearly verbatim from `## Goals`. Condense `## Overview` into a brief high-level summary and let `## Goals` carry the enumerated capability list.

### `[Nit]` SY-3: Synonym Redundancy and Filler Modifiers
* **Affected Scope:** [flet-antibody-studio.md](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L8), [L10](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L10), [L99](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L99), [L182](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L182), [L197](file:///Users/fkt/Downloads/PROJECT/Fc-engineering-studio/docs/designs/flet-antibody-studio.md#L197)
* **Issue & Recommendation:** Simplify redundant phrasing: `"single, monolithic"` → `"single"` (L8); `"multiple distinct ... four unique ... discrete"` → `"multiple ... four ... modular"` (L10); drop trailing `"in one self-contained text file"` (L182); and `"scoped strictly to"` → `"scoped to"` (L99, L197).
