## 2023-10-24 - [Rich Markup Injection in Textual/Rich Log Widget]
**Vulnerability:** User-provided inputs (such as custom mutations like `S239[bold]X`) could be passed directly to a `Log` widget's `write` method, which by default parses strings as Rich markup. This could lead to UI rendering issues or crashing the application entirely with unclosed tags (e.g. `[unclosed`).
**Learning:** In Textual and Rich, if you accept user-provided strings and display them using a `write` or print method that evaluates markup by default, you must explicitly disable markup parsing or escape the string first. The `escape` function from `rich.markup` is necessary when interpolating unvalidated user inputs within other markup strings (e.g., `f"[red]{user_input}[/red]"` will crash if `user_input` contains unclosed brackets).
**Prevention:** Always use `escape()` from `rich.markup` when displaying dynamically generated text containing user input via Rich components, or pass `markup=False` when printing plain text.

## 2026-04-18 - [Exception Details Leakage & Missing Input Limits]
**Vulnerability:** Raw Python exceptions were being exposed to the user UI, potentially leaking internal application logic. The custom mutations input field lacked a length limit, posing a DoS risk.
**Learning:** In user-facing inputs and error messages, always provide generic errors to avoid info-leakage and enforce reasonable boundaries (like max_length) on components that accept user text.
**Prevention:** Apply input validation/limits on all interactive inputs and sanitize exception messages before they are rendered.
