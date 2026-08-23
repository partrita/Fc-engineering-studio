import re
import functools
from typing import Tuple
from textual import on
from textual.app import App, ComposeResult
from textual.containers import Vertical, Horizontal
from textual.widgets import (
    Button,
    Footer,
    Header,
    Input,
    Label,
    Log,
    Static,
    SelectionList,
    OptionList,
    Pretty,
)

from textual.widgets.selection_list import Selection
from textual.screen import Screen
import pyperclip
from rich.markup import escape

from fc_engineer.config import SEQUENCES, COMMON_MUTATIONS, log_sanitized_error
from fc_engineer.core import apply_mutations, diff_sequences, generate_batch, analyze_developability
from fc_engineer.exporters import format_genbank, format_csv_report

ANTIBODY_ASCII = r"""
  _____                                                 
 |  ___|__                                              
 | |_ / __|                                             
 |  _| (__         _                      _             
 |_|__\___|   __ _(_)_ __   ___  ___ _ __(_)_ __   __ _ 
  / _ \ '_ \ / _` | | '_ \ / _ \/ _ \ '__| | '_ \ / _` |
 |  __/ | | | (_| | | | | |  __/  __/ |  | | | | | (_| |
  \___|_|\__, |_|_| |_|\___|\___|_|  |_|_| |_|\__, |
  ___| |_ _  |___/_| (_) ___                      |___/ 
 / __| __| | | |/ _` | |/ _ \                           
 \__ \ |_| |_| | (_| | | (_) |                          
 |___/\__|\__,_|\__,_|_|\___/                           
                                                         
"""

# --- Screens ---

class WelcomeScreen(Screen):
    BINDINGS = [("enter", "next", "Next")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container", id="welcome-container"):
            yield Static(ANTIBODY_ASCII, id="antibody-ascii")
            yield Label("Antibody Fc Sequence Designer", classes="subtitle")
            yield Static("[dim]Design human IgG Fc mutants based on EU Numbering.[/]", id="overview")
            yield Label("Press ENTER to start", id="press-enter")
        yield Footer()

    def action_next(self) -> None:
        self.app.push_screen(IsotypeScreen())

class IsotypeScreen(Screen):
    BINDINGS = [("enter", "next", "Next"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("STEP 1", classes="title")
            yield Label("Select IgG Isotype", classes="subtitle")
            yield OptionList(
                *[escape(iso.upper()) for iso in SEQUENCES.keys()],
                id="iso-list"
            )
        yield Footer()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.action_next()

    def action_next(self) -> None:
        self.app.selected_isotype = ""
        try:
            iso_list = self.query_one("#iso-list", OptionList)
            if iso_list.highlighted is None:
                self.notify("Please select an Isotype.", severity="warning")
                return
            selected_iso = list(SEQUENCES.keys())[iso_list.highlighted]
            self.app.selected_isotype = selected_iso
            self.app.push_screen(AllotypeScreen())
        except Exception as e:
            log_sanitized_error(self.log, f"Error in IsotypeScreen.action_next", e)
            self.notify("An error occurred. Please check configuration.", severity="error")

    def action_back(self) -> None:
        self.app.selected_isotype = ""
        self.app.pop_screen()

class AllotypeScreen(Screen):
    BINDINGS = [("enter", "next", "Next"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("STEP 2", classes="title")
            yield Label(f"Select {escape(self.app.selected_isotype.upper())} Allotype", classes="subtitle")
            allotypes = SEQUENCES.get(self.app.selected_isotype, {})
            if not isinstance(allotypes, dict):
                allotypes = {}
            yield OptionList(
                *[escape(allo.capitalize()) for allo in allotypes.keys()],
                id="allo-list"
            )
        yield Footer()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.action_next()

    def action_next(self) -> None:
        self.app.selected_allotype = ""
        try:
            allo_list = self.query_one("#allo-list", OptionList)
            if allo_list.highlighted is None:
                self.notify("Please select an Allotype.", severity="warning")
                return
            allotypes = SEQUENCES.get(self.app.selected_isotype, {})
            if not isinstance(allotypes, dict):
                allotypes = {}
            selected_allo = list(allotypes.keys())[allo_list.highlighted]
            self.app.selected_allotype = selected_allo
            self.app.push_screen(MutationScreen())
        except Exception as e:
            log_sanitized_error(self.log, f"Error in AllotypeScreen.action_next", e)
            self.notify("An error occurred. Please check configuration.", severity="error")

    def action_back(self) -> None:
        self.app.selected_allotype = ""
        self.app.pop_screen()

class MutationScreen(Screen):
    BINDINGS = [("enter", "generate", "Generate"), ("escape", "back", "Back"), ("b", "batch", "Batch")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("STEP 3", classes="title")
            yield Label("Mutations (Optional)", classes="subtitle")
            yield Label("Presets:", classes="subtitle")
            with Horizontal(id="selection-area"):
                yield SelectionList[str](
                    *[Selection(escape(item.get("label", "Unknown")), item.get("value", ""), False) for item in COMMON_MUTATIONS if isinstance(item, dict)],
                    id="list-common"
                )
                yield Pretty([], id="selected-preview")
            yield Label("Custom (e.g. S239D/I332E):", classes="subtitle")
            yield Input(placeholder="None", id="input-custom", max_length=100, restrict=r"^[a-zA-Z0-9/, ]*\Z")
            yield Button("Generate FASTA", variant="primary", id="btn-gen")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#list-common").border_title = "Common Mutants"
        self.query_one("#selected-preview").border_title = "Selected"

    @on(SelectionList.SelectedChanged)
    def update_selected_view(self) -> None:
        self.query_one("#selected-preview", Pretty).update(self.query_one("#list-common", SelectionList).selected)

    def action_generate(self) -> None:
        self.app.all_mutants = ""
        try:
            selected_presets = self.query_one("#list-common", SelectionList).selected
            preset_str = "/".join(selected_presets)
            custom_str = self.query_one("#input-custom", Input).value

            # Security Enhancement: Strictly sanitize custom input to remove any invalid characters
            # acting as a secondary defense layer behind the Input widget's restrict regex.
            custom_str = re.sub(r'[^a-zA-Z0-9/, ]', '', custom_str)

            self.app.all_mutants = "/".join(filter(None, [preset_str, custom_str]))
            self.app.push_screen(ResultScreen())
        except Exception as e:
            log_sanitized_error(self.log, f"Error in MutationScreen.action_generate", e)
            self.notify("An error occurred. Please check your inputs.", severity="error")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "btn-gen": self.action_generate()

    def action_back(self) -> None:
        self.app.all_mutants = ""
        self.app.pop_screen()

    def action_batch(self) -> None:
        self.app.push_screen(BatchScreen())

class BatchScreen(Screen):
    BINDINGS = [
        ("escape", "back", "Back"),
        ("enter", "generate", "Generate"),
        ("ctrl+y", "copy_all", "Copy All"),
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("BATCH MODE", classes="title")
            yield Label("Comma-separated mutation sets", classes="subtitle")
            yield Input(
                placeholder="e.g. L234A/L235A, M428L/N434S, N297A",
                id="input-batch",
                max_length=500,
                restrict=r"^[a-zA-Z0-9/, ]*\Z",
            )
            yield Button("Generate All", variant="primary", id="btn-batch-gen")
            yield Log(id="batch-box")
            yield Static("[dim]Enter: Generate | Ctrl+Y: Copy All | Esc: Back[/]", id="batch-help")
        yield Footer()

    def action_generate(self) -> None:
        self.app.batch_results = []  # SECURITY: Reset previous batch state first

        isotype = self.app.selected_isotype
        allotype = self.app.selected_allotype
        batch_box = self.query_one("#batch-box", Log)
        batch_box.clear()

        try:
            iso_data = SEQUENCES.get(isotype, {})
            if not isinstance(iso_data, dict):
                iso_data = {}
            base_seq = iso_data.get(allotype, "")
            if not isinstance(base_seq, str) or not base_seq:
                batch_box.write(f"[bold red]Error: Base sequence for {escape(isotype)} {escape(allotype)} not found.[/]")
                return

            raw_value = self.query_one("#input-batch", Input).value
            # Security Enhancement: Secondary sanitization behind the Input widget's restrict regex
            raw_value = re.sub(r"[^a-zA-Z0-9/, ]", "", raw_value)

            results = generate_batch(base_seq, raw_value, isotype)
            if not results:
                batch_box.write("[yellow]No batch combinations provided.[/]")
                return
            if len(results) == 1 and results[0][0] == "" and results[0][2]:
                # Bounds violation reported by core logic
                for err in results[0][2]:
                    batch_box.write(f"[red]• {escape(err)}[/]")
                return

            success_count = 0
            for combo, seq_out, errors in results:
                display_muts = combo.replace("/", "_")
                header = f"{isotype.upper()}_{allotype.capitalize()}_{display_muts}"
                batch_box.write(escape(f">{header}\n{seq_out}"))
                if errors:
                    for err in errors:
                        batch_box.write(f"[red]• {escape(combo)}: {escape(err)}[/]")
                else:
                    success_count += 1
                    self.app.batch_results.append((header, seq_out))
                batch_box.write("")
            self.log.info(f"Audit: Batch generated {success_count} sequences for {isotype} {allotype}.")
        except Exception as e:
            batch_box.write("[bold red]An unexpected error occurred during batch generation.[/]")
            log_sanitized_error(self.log, "Error in BatchScreen.action_generate", e)

    def action_copy_all(self) -> None:
        if not getattr(self.app, "batch_results", []):
            self.notify("Nothing to copy. Generate a batch first.", severity="warning")
            return
        joined = "\n\n".join(f">{h}\n{s}" for h, s in self.app.batch_results)
        self.app.copy_text_secure(joined, f"{len(self.app.batch_results)} FASTA sequences")

    def action_back(self) -> None:
        # SECURITY: Wipe sensitive batch state on backward navigation
        try:
            if hasattr(self.app, "_clipboard_timer") and self.app._clipboard_timer is not None:
                self.app._clipboard_timer.stop()
                self.app._clipboard_timer = None
            if hasattr(self.app, "copied_fasta") and self.app.copied_fasta:
                self.app.clear_clipboard(self.app.copied_fasta)
        finally:
            self.app.batch_results = []
            self.app.pop_screen()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "btn-batch-gen": self.action_generate()


class ResultScreen(Screen):
    BINDINGS = [
        ("escape", "back", "Back"),
        ("ctrl+y", "copy_to_clipboard", "Copy"),
        ("ctrl+g", "copy_genbank", "Copy GB"),
        ("ctrl+r", "copy_csv", "Copy CSV"),
        ("q", "quit_to_main", "Main Menu")
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("RESULT", classes="title")
            yield Label("FASTA Sequence", classes="subtitle")
            yield Log(id="result-box")
            yield Static("[dim]Ctrl+Y: FASTA | Ctrl+G: GenBank | Ctrl+R: CSV | Esc: Back | Q: Menu[/]", id="result-help")
        yield Footer()

    def on_mount(self) -> None:
        self.generate_fasta()

    def generate_fasta(self) -> None:
        self.app.last_fasta = ""  # SECURITY: Clear previous state to prevent stale data leakage

        isotype = self.app.selected_isotype
        allotype = self.app.selected_allotype
        all_mutants = self.app.all_mutants

        result_box = self.query_one("#result-box", Log)
        result_box.clear()

        try:
            isotype_data = SEQUENCES.get(isotype, {})
            if not isinstance(isotype_data, dict):
                isotype_data = {}
            base_seq = isotype_data.get(allotype)
            if not base_seq:
                result_box.write(f"[bold red]Error: Base sequence for {escape(isotype)} {escape(allotype)} not found.[/]")
                return
            if not isinstance(base_seq, str):
                result_box.write(f"[bold red]Error: Invalid sequence data format for {escape(isotype)} {escape(allotype)}.[/]")
                return

            mutant_seq, errors = apply_mutations(base_seq, all_mutants, isotype)

            if errors:
                result_box.write("[bold red]Errors during processing:[/]")
                # SECURITY: Escape user-provided input in errors to prevent Rich markup injection
                for err in errors: result_box.write(f"[red]• {escape(err)}[/]")
                return

            display_muts = all_mutants.replace("/", "_") if all_mutants else "WT"
            header = f"{isotype.upper()}_{allotype.capitalize()}_{display_muts}"
            fasta = f">{header}\n{mutant_seq}"

            # SECURITY: Escape user-provided sequence/header data to prevent Rich markup injection
            result_box.write(escape(fasta))
            self.app.last_fasta = fasta

            # Sequence comparison view: per-position diff against WT (EU Numbering)
            diffs = diff_sequences(base_seq, mutant_seq, isotype)
            if diffs:
                result_box.write("")
                result_box.write(f"[bold]Mutations vs WT ({len(diffs)}):[/]")
                for pos, wt_aa, mut_aa in diffs:
                    result_box.write(f"[yellow]• {escape(wt_aa)}{pos}{escape(mut_aa)}[/]")

            # Developability summary: sequence liabilities annotated with EU positions
            result_box.write("")
            result_box.write("[bold]Developability summary:[/]")
            for category, positions in analyze_developability(mutant_seq, isotype):
                if positions:
                    positions_str = ", ".join(str(p) for p in positions)
                    result_box.write(f"• {escape(category)}: {positions_str}")
                else:
                    result_box.write(f"• {escape(category)}: none")

            # SECURITY: Audit log for sensitive intellectual property operation (Sequence Generation)
            self.log.info(f"Audit: Generated FASTA sequence for {isotype} {allotype} with mutations {display_muts}")
        except Exception as e:
            result_box.write("[bold red]An unexpected error occurred during sequence generation.[/]")
            log_sanitized_error(self.log, f"Error in generate_fasta", e)

    def action_copy_to_clipboard(self) -> None:
        last_fasta = self.app.last_fasta if hasattr(self.app, "last_fasta") else ""
        self.app.copy_text_secure(last_fasta, "FASTA sequence")

    def _current_export_state(self) -> Tuple[str, str, str, str]:
        """Return (header, base_seq, mut_seq, mutations) derived from current app state."""
        header, _, mut_seq = (self.app.last_fasta or "").partition("\n")
        header = header.lstrip(">")
        isotype = getattr(self.app, "selected_isotype", "")
        allotype = getattr(self.app, "selected_allotype", "")
        mutations = getattr(self.app, "all_mutants", "")
        base_seq = ""
        iso_data = SEQUENCES.get(isotype, {})
        if isinstance(iso_data, dict):
            candidate = iso_data.get(allotype, "")
            base_seq = candidate if isinstance(candidate, str) else ""
        return header, base_seq, mut_seq, mutations

    def action_copy_genbank(self) -> None:
        header, _, mut_seq, _ = self._current_export_state()
        if not mut_seq:
            self.notify("Nothing to export.", severity="warning")
            return
        record = format_genbank(header, mut_seq, definition=header)
        if not record:
            self.notify("Failed to build GenBank record.", severity="error")
            return
        self.app.copy_text_secure(record, "GenBank record")

    def action_copy_csv(self) -> None:
        header, base_seq, mut_seq, mutations = self._current_export_state()
        if not mut_seq:
            self.notify("Nothing to export.", severity="warning")
            return
        report = format_csv_report(
            getattr(self.app, "selected_isotype", ""),
            getattr(self.app, "selected_allotype", ""),
            mutations, base_seq, mut_seq,
        )
        if not report:
            self.notify("Failed to build CSV report.", severity="error")
            return
        self.app.copy_text_secure(report, f"CSV report ({header})")

    def action_quit_to_main(self) -> None:
        # SECURITY: Wipe sensitive state upon returning to main menu
        try:
            if hasattr(self.app, "_clipboard_timer") and self.app._clipboard_timer is not None:
                self.app._clipboard_timer.stop()
                self.app._clipboard_timer = None
            if hasattr(self.app, "copied_fasta") and self.app.copied_fasta:
                self.app.clear_clipboard(self.app.copied_fasta)
        finally:
            self.app.selected_isotype = ""
            self.app.selected_allotype = ""
            self.app.all_mutants = ""
            self.app.last_fasta = ""
            self.app.copied_fasta = ""
            self.app.batch_results = []
            while len(self.app.screen_stack) > 1:
                self.app.pop_screen()

    def action_back(self) -> None:
        # SECURITY: Wipe external state on backward navigation
        try:
            if hasattr(self.app, "_clipboard_timer") and self.app._clipboard_timer is not None:
                self.app._clipboard_timer.stop()
                self.app._clipboard_timer = None
            if hasattr(self.app, "copied_fasta") and self.app.copied_fasta:
                self.app.clear_clipboard(self.app.copied_fasta)
        finally:
            self.app.last_fasta = ""
            self.app.copied_fasta = ""
            self.app.pop_screen()

# --- Main App ---

class MutantApp(App):
    TITLE = "Fc Engineering Studio"
    CSS = """
    Screen {
        align: center middle;
    }

    .minimal-container {
        width: 70;
        height: auto;
        padding: 1 2;
        background: $panel;
        color: $text;
        border: $secondary tall;
    }

    .title {
        text-align: center;
        width: 100%;
        text-style: bold;
        margin-bottom: 1;
    }

    .subtitle {
        text-align: center;
        width: 100%;
        margin-bottom: 1;
    }

    SelectionList {
        padding: 1;
        border: solid $accent;
        width: 80%;
        height: 80%;
    }

    SelectionList .selection-list--button {
        background: transparent;
        color: transparent;
    }

    SelectionList .selection-list--button-selected {
        background: $accent;
        color: white;
        text-style: bold;
    }

    SelectionList .selection-list--button-highlighted {
        background: transparent;
        color: transparent;
    }

    SelectionList .selection-list--button-selected-highlighted {
        background: $accent;
        color: white;
        text-style: bold;
    }

    #antibody-ascii {
        width: 100%;
        text-align: center;
        color: $accent;
        margin-bottom: 1;
    }

    #iso-list, #allo-list {
        height: 10;
        border: solid $accent;
    }

    #list-common {
        /* Overridden by SelectionList styles if needed */
    }

    #selected-preview {
        width: 1fr;
        height: 10;
        border: solid $accent;
    }

    #selection-area {
        height: 12;
    }

    #input-custom {
        margin-bottom: 1;
    }

    #input-batch {
        margin-bottom: 1;
    }

    #btn-batch-gen {
        width: 100%;
        margin-bottom: 1;
    }

    #batch-box {
        height: 16;
        border: solid $accent;
        margin-bottom: 1;
    }

    #batch-help {
        text-align: center;
        width: 100%;
    }

    #btn-gen {
        width: 100%;
    }

    #result-box {
        height: 10;
        border: solid $accent;
        margin-bottom: 1;
    }

    #result-help {
        text-align: center;
        width: 100%;
    }
    """

    def on_mount(self) -> None:
        self.theme = "nord"
        self.selected_isotype = ""
        self.selected_allotype = ""
        self.all_mutants = ""
        self.last_fasta = ""
        self.copied_fasta = ""
        self.batch_results = []
        self.push_screen(WelcomeScreen())

    def on_unmount(self) -> None:
        """Security: Clear clipboard on application exit if it contains sensitive data."""
        try:
            if hasattr(self, "_clipboard_timer") and self._clipboard_timer is not None:
                self._clipboard_timer.stop()
                self._clipboard_timer = None
            if hasattr(self, "copied_fasta") and self.copied_fasta:
                self.clear_clipboard(self.copied_fasta)
        finally:
            if hasattr(self, "copied_fasta"):
                self.copied_fasta = ""

    def clear_clipboard(self, content_to_clear: str) -> None:
        """Security: Clears the clipboard to prevent sensitive data exposure."""
        try:
            if pyperclip.paste() == content_to_clear:
                pyperclip.copy("")
                self.log.info("Clipboard automatically cleared for security.")
        except Exception as e:
            log_sanitized_error(self.log, f"Error clearing clipboard", e)
        finally:
            if hasattr(self, "copied_fasta"):
                self.copied_fasta = ""

    def copy_text_secure(self, text: str, notify_label: str = "Sequence") -> None:
        """Copy text to the OS clipboard under the auto-clear security policy.

        Lives on the persistent App (not on ephemeral screens) so the auto-clear
        timer and teardown hooks always remain in control of the copied data.
        """
        # SECURITY: Clear any previous clipboard state before copying new data
        if hasattr(self, "copied_fasta") and self.copied_fasta:
            self.clear_clipboard(self.copied_fasta)
        self.copied_fasta = ""  # SECURITY: Clear state before copying
        if not text:
            return
        try:
            pyperclip.copy(text)
            self.copied_fasta = text
            # SECURITY: Audit log for sensitive intellectual property operation (Data Export to OS)
            self.log.info("Audit: Copied proprietary sequence data to OS clipboard.")
            self.notify(f"{notify_label} copied! (Will auto-clear in 30s)")
            # Security: Auto-clear clipboard after 30 seconds
            # Ensure overlapping timers are cancelled so the timer doesn't prematurely clear a newly copied item
            if hasattr(self, "_clipboard_timer") and self._clipboard_timer is not None:
                self._clipboard_timer.stop()
            self._clipboard_timer = self.set_timer(30, functools.partial(self.clear_clipboard, self.copied_fasta))
        except Exception as e:
            log_sanitized_error(self.log, "Error copying to clipboard", e)
            self.notify("Error copying to clipboard. See logs.", severity="error")


def main():
    MutantApp().run()

if __name__ == "__main__":
    main()
