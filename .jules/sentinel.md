## 2023-10-24 - [Rich Markup Injection in Textual/Rich Log Widget]
**Vulnerability:** User-provided inputs (such as custom mutations like `S239[bold]X`) could be passed directly to a `Log` widget's `write` method, which by default parses strings as Rich markup. This could lead to UI rendering issues or crashing the application entirely with unclosed tags (e.g. `[unclosed`).
**Learning:** In Textual and Rich, if you accept user-provided strings and display them using a `write` or print method that evaluates markup by default, you must explicitly disable markup parsing or escape the string first. The `escape` function from `rich.markup` is necessary when interpolating unvalidated user inputs within other markup strings (e.g., `f"[red]{user_input}[/red]"` will crash if `user_input` contains unclosed brackets).
**Prevention:** Always use `escape()` from `rich.markup` when displaying dynamically generated text containing user input via Rich components, or pass `markup=False` when printing plain text.

## 2026-04-17 - [Missing Input Validation & Length Limits]
**Vulnerability:** The application parsed user-provided mutation strings via basic index slicing and lacked bounds/length constraints in the UI. This allowed invalid formats to crash the application, leak internal Python exception traces, and exposed a potential DoS risk through unbounded input length processing.
**Learning:** In Textual UIs, custom input fields should leverage built-in properties like `max_length` and `restrict` regex to limit size and characters early. Similarly, logic functions parsing strings must enforce strict regex structural validation instead of relying on optimistic string slicing, to ensure predictable fail-states.
**Prevention:** Always implement defense-in-depth: constrain input length and character sets at the UI layer using Textual's input properties, and rigorously validate structural format at the core logic layer using `re.match` before parsing.

## 2024-05-24 - [Unhandled UI State Exceptions Exposing Stack Traces]
**Vulnerability:** TUI application components (e.g. OptionList selections, sequence dictionary lookups) lacked validation for missing states or invalid keys. When a user advanced without making a selection or if corrupted sequence data was provided, unhandled `TypeError` or `KeyError` exceptions bypassed application logic, crashing the UI and dumping raw Python stack traces into the terminal. This exposed internal codebase paths and architecture details.
**Learning:** In TUI frameworks like Textual, user actions can occur when UI state is not fully populated (e.g., pressing Enter on an empty list). Dictionary lookups for critical data (like base sequences) must not assume keys exist. Failing to handle these states gracefully leaks sensitive technical information and degrades security through poor error handling.
**Prevention:** Always validate UI state (e.g., check `highlighted is not None`) before accessing properties. Use `.get()` with defaults for dictionary lookups instead of direct bracket access. Wrap critical processing blocks in `try...except` and log errors securely while presenting a generic, non-exposing message to the UI via user notifications or log widgets.

## 2024-06-15 - [Regex Validation Bypass via Trailing Newline]
**Vulnerability:** The application used `$` as an end-of-string anchor in `re.match` for input validation (e.g., parsing mutations). In Python, `$` matches either the end of the string OR just before a trailing newline. This allowed inputs with trailing newlines (e.g., `A118X\n`) to pass validation but cause unhandled exceptions later when parsing (e.g., `int()` on sliced strings).
**Learning:** Python's `re` module behavior for `$` can lead to validation bypasses if inputs contain trailing newlines. The `\Z` anchor should be used for strict end-of-string matching.
**Prevention:** Always use `\Z` instead of `$` in Python regular expressions when strict input format validation is required, and ensure UI input fields restrict whitespace characters to explicitly permitted ones (like a space) rather than the broad `\s` which includes newlines.

## 2024-06-25 - [Missing Input Limits leading to DoS Risk]
**Vulnerability:** The mutation parser (`parse_mutation`) and processor (`apply_mutations`) lacked bounds checking on the size and quantity of inputs. This exposed the application to potential resource exhaustion (Denial of Service) attacks if a user provided an excessively long mutation string or an enormous list of mutations.
**Learning:** Even internal or UI-driven string parsing functions need constraints. Regular expressions and loops processing user input without bounds can be abused to consume excessive CPU or memory.
**Prevention:** Always implement hard limits on input lengths (e.g., max string length) and processing bounds (e.g., maximum number of items in a list) at the core logic layer, regardless of UI-level restrictions.

## 2026-04-26 - [Missing YAML File Size Validation]
**Vulnerability:** The application loaded configuration `yaml` files without checking their sizes. This exposed a potential Denial of Service (DoS) vulnerability via memory exhaustion, where processing an exceptionally large file would consume excessive resources.
**Learning:** File input streams must be strictly bounded before parsing content. Relying entirely on safe loaders without checking the incoming buffer size leaves applications vulnerable to resource exhaustion.
**Prevention:** Always verify file sizes against a defined maximum limit (e.g., 1MB) using `os.path.getsize` before opening and parsing data, particularly with external configuration files or user-provided files.
## 2024-08-16 - [Missing File Size Limitations for YAML Loading]
**Vulnerability:** The application used `yaml.safe_load()` without checking the size of the underlying files (`sequences.yaml` and `mutants.yaml`). This allowed the potential for Denial of Service (DoS) attacks via memory exhaustion if a user provided an excessively large file.
**Learning:** `yaml.safe_load()` prevents arbitrary code execution but does not protect against memory exhaustion from very large files. File sizes should always be validated before attempting to read and parse them into memory.
**Prevention:** Implemented a file size check (`os.path.getsize(path) <= MAX_FILE_SIZE`) before opening and parsing YAML files. Added explicit error logging when the limit is exceeded.
## 2026-04-24 - [Missing Sequence Type Validation leading to Stack Trace Leak]
**Vulnerability:** The application retrieved the base sequence from a YAML-sourced dictionary but did not validate its type. If the sequence data was malformed (e.g., an integer or list instead of a string), passing it to the `apply_mutations` function would cause an unhandled `TypeError` (e.g., `'int' object is not iterable`) when attempting to convert it to a list. This would crash the TUI and leak internal Python exception traces.
**Learning:** In addition to validating the top-level structure of loaded data files (e.g., ensuring the root is a dictionary), deeply nested values that are passed to critical processing logic must also be explicitly type-checked before use to prevent unexpected fail-states and architecture leaks.
**Prevention:** Always validate the type of data retrieved from configuration files (e.g., using `isinstance(val, str)`) before processing it, especially when the logic assumes a specific iterable or string behavior. Log formatting errors securely without leaking system internals.
## 2026-04-27 - [Unicode Digit Injection in Regex Validation]
**Vulnerability:** The application used  in the regex for mutation parsing (). In Python 3,  matches all Unicode digits (e.g., Arabic numerals `١٢٣`), not just ASCII digits `0-9`. This could allow an attacker to bypass basic input validation expectations, leading to unexpected behavior when  processes these non-ASCII digits.
**Learning:** Python's  character class is broader than expected and matches any Unicode digit. When strict input format validation is required, especially for numeric inputs that map to specific logic, ASCII digits should be explicitly enforced.
**Prevention:** Always use  instead of  in Python regular expressions when strict validation of ASCII-only numeric input is required.
## 2024-04-27 - [Unicode Digit Injection in Regex Validation]
**Vulnerability:** The application used `\d+` in the regex for mutation parsing (`parse_mutation`). In Python 3, `\d` matches all Unicode digits (e.g., Arabic numerals `١٢٣`), not just ASCII digits `0-9`. This could allow an attacker to bypass basic input validation expectations, leading to unexpected behavior when `int()` processes these non-ASCII digits.
**Learning:** Python's `\d` character class is broader than expected and matches any Unicode digit. When strict input format validation is required, especially for numeric inputs that map to specific logic, ASCII digits should be explicitly enforced.
**Prevention:** Always use `[0-9]` instead of `\d` in Python regular expressions when strict validation of ASCII-only numeric input is required.

## 2024-10-27 - [Type Coercion Regressions vs Explicit Validation]
**Vulnerability:** While trying to secure YAML data structures against type-based crashes (e.g., encountering ints instead of strings), modifying `main.py` to mutate loaded dictionaries (`isotypes[k] = ...`) or silently casting arbitrary keys/values to strings without bounds caused severe logical regressions and potential `NameError` crashes due to breaking explicit fallback logic.
**Learning:** Security input validation should never silently mutate or truncate valid configuration structures out of strict schema guessing. Relying on strict type assertions via `isinstance` while gracefully falling back to empty/safe defaults is the correct way to validate input without hallucinating a rigid schema that breaks functionality.
**Prevention:** Use defensive type checking (`isinstance()`) combined with conditional ternary fallbacks (e.g. `val if isinstance(val, dict) and all(...) else {}`) to secure input without introducing side effects or losing application state.
## $(date +%Y-%m-%d) - [Missing Sequence Type Validation leading to Stack Trace Leak]
**Vulnerability:** The application retrieved the base sequence from a YAML-sourced dictionary but did not validate the intermediate dictionary's type (e.g., `isotype_data`). If the `yaml` file contained invalid data (e.g., `isotypes: { igg1: "not a dict" }`), calling `isotype_data.get(allotype)` would raise an `AttributeError` on a string, crashing the TUI and leaking internal Python exception traces.
**Learning:** In addition to validating the top-level structure of loaded data files, deeply nested values that are passed to critical processing logic (or used as dictionaries) must also be explicitly type-checked before use to prevent unexpected fail-states and architecture leaks.
**Prevention:** Always validate the type of data retrieved from configuration files (e.g., using `isinstance(val, dict)`) before calling dictionary methods like `.get()` on it.

## 2024-05-24 - [Unnecessary Markup Escaping on Trusted Internal Data]
**Vulnerability:** The application was proposed to have Rich markup injection vulnerabilities from strings loaded via static internal configuration YAML files (e.g. `sequences.yaml`).
**Learning:** Escaping markup on trusted internal data files that are not modified by end users is a form of "security theater". If an attacker can modify internal application source or config files, they already have a higher level of compromise.
**Prevention:** Focus input escaping specifically on actual external user input paths, rather than treating trusted internal static configurations as a threat model.

## 2024-05-24 - [Missing Type Validation for List Labels in TUI]
**Vulnerability:** The application loaded `common_mutations` presets from a YAML file, checking that `value` was a string but omitting the type check for `label`. The UI component `SelectionList` expects a string for the prompt label. Providing a malformed YAML containing list types or integers as labels would crash the application with a `TypeError` and leak stack traces to the terminal.
**Learning:** All properties loaded from YAML that interact with specific UI expectations (e.g., Textual labels) must be strictly type-checked.
**Prevention:** Ensured `isinstance(item.get("label"), str)` is verified before returning the parsed presets from `load_yaml_data`.

## 2026-05-02 - [YAML Bomb (Billion Laughs) DoS Vulnerability despite File Size Limits]
**Vulnerability:** The application used `yaml.safe_load()` combined with a 1MB file size limit to load configurations. However, `yaml.safe_load()` still evaluates YAML aliases and anchors. An attacker could provide a very small YAML file (well under 1MB) containing heavily nested aliases (a "YAML Bomb" or "Billion Laughs" attack) that expand exponentially in memory, causing a Denial of Service (DoS) via memory exhaustion.
**Learning:** Checking file size is insufficient to prevent memory exhaustion when parsing formats that support data expansion features like aliases. `yaml.safe_load` protects against arbitrary object instantiation but does not block alias expansion by default.
**Prevention:** Always implement a custom SafeLoader that explicitly raises an error on `yaml.events.AliasEvent` (e.g., overriding `compose_node` and using `self.check_event()`) when loading untrusted YAML, to fully mitigate exponential expansion attacks.

## 2026-05-03 - [Infinite Stream DoS Vulnerability bypassing File Size Limits]
**Vulnerability:** The application used `os.path.exists()` and `os.path.getsize()` to validate file size before passing the file object to `yaml.load()`. However, special device files (like `/dev/zero`) or named pipes can return a size of `0`, bypassing the size check. Because `yaml.load(f)` reads the stream until EOF, reading an infinite stream would lead to CPU and memory exhaustion (Denial of Service).
**Learning:** `os.path.getsize()` is unreliable for special files. Furthermore, passing an unconstrained file-like object directly to a parser allows the parser to consume unbounded memory if the stream does not end.
**Prevention:** Always verify that a path points to a regular file using `os.path.isfile()` rather than just `os.path.exists()`. Additionally, apply defense-in-depth by explicitly reading the file contents with a bounded size limit (e.g., `content = f.read(MAX_FILE_SIZE + 1)`) and checking the length before passing the data to the parser.
## $(date +%Y-%m-%d) - [Silent Exception Handling Hiding Security/Config Issues]
**Vulnerability:** The application used `try...except Exception as e: pass` blocks when attempting to load critical configuration files (`sequences.yaml` and `mutants.yaml`). This pattern silently swallowed errors, effectively blinding the application to potential security issues or misconfigurations (e.g., malformed data, permissions errors).
**Learning:** Silent failures violate the "fail securely" and "maintain security visibility" principles. If configuration parsing fails, the system might operate in an unpredictable default state without any diagnostic trace for developers or security monitors.
**Prevention:** Always log exceptions explicitly, even if the application can recover or fallback. Use mechanisms like `print(..., file=sys.stderr)` or standard logging libraries to record the `Exception` details when handling external inputs or critical configurations.
