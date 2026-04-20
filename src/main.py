import sys
import os
import re
from typing import Dict, List, Optional, Tuple, Set
import yaml

from textual import on
from textual.app import App, ComposeResult
from textual.containers import Container, Vertical, Center, Middle, Horizontal
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
from textual.binding import Binding
from textual.screen import Screen
import pyperclip
from rich.markup import escape

# --- Configuration & Data Loading ---

def load_yaml_data():
    base_path = os.path.dirname(__file__)
    seq_path = os.path.join(base_path, "sequences.yaml")
    mut_path = os.path.join(base_path, "mutants.yaml")
    
    isotypes = {}
    common_muts = []
    
    try:
        if os.path.exists(seq_path):
            with open(seq_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                isotypes = data.get("isotypes", {})
    except Exception as e:
        pass

    try:
        if os.path.exists(mut_path):
            with open(mut_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                common_muts = data.get("common_mutations", [])
    except Exception as e:
        pass
        
    return isotypes, common_muts

SEQUENCES, COMMON_MUTATIONS = load_yaml_data()

# --- Core Logic ---

EU_START = 118

def get_residue_index(pos: int, isotype: str) -> Optional[int]:
    if isotype == "igg1": return pos - EU_START
    elif isotype in ["igg2", "igg4"]:
        if pos <= 222: return pos - EU_START
        elif 223 <= pos <= 225: return None
        else: return pos - EU_START - 3
    return None

def parse_mutation(m_str: str) -> Tuple[str, int, str]:
    m_str = m_str.upper()
    if not re.match(r"^[A-Z]\d+[A-Z]$", m_str):
        raise ValueError(f"Invalid mutation format: {m_str}")
    wt_aa = m_str[0]
    pos = int(m_str[1:-1])
    mut_aa = m_str[-1]
    return wt_aa, pos, mut_aa

def apply_mutations(sequence: str, mutants_str: str, isotype: str) -> Tuple[str, List[str]]:
    if not mutants_str: return sequence, []
    mut_list = [m.strip() for m in mutants_str.replace(',', '/').split('/') if m.strip()]
    seq_list = list(sequence)
    errors = []
    for m in mut_list:
        try:
            wt_aa, pos, mut_aa = parse_mutation(m)
            index = get_residue_index(pos, isotype)
            if index is None: errors.append(f"Position {pos} is a Gap in {isotype}."); continue
            if index < 0 or index >= len(seq_list): errors.append(f"Position {pos} is out of range."); continue
            if seq_list[index] != wt_aa: errors.append(f"Pos {pos}: Expected '{wt_aa}', found '{seq_list[index]}'."); continue
            seq_list[index] = mut_aa
        except Exception as e: errors.append(f"Format error ({m}): {str(e)}")
    return "".join(seq_list), errors

ANTIBODY_ASCII = r"""
  _____                                                 
 |  ___|__                                              
 | |_ / __|                                             
 |  _| (__         _                      _             
 |_|__\___|   __ _(_)_ __   ___  ___ _ __(_)_ __   __ _ 
  / _ \ '_ \ / _` | | '_ \ / _ \/ _ \ '__| | '_ \ / _` |
 |  __/ | | | (_| | | | | |  __/  __/ |  | | | | | (_| |
  \___|_| |_|\__, |_|_| |_|\___|\___|_|  |_|_| |_|\__, |
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
                *[iso.upper() for iso in SEQUENCES.keys()],
                id="iso-list"
            )
        yield Footer()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.action_next()

    def action_next(self) -> None:
        iso_list = self.query_one("#iso-list", OptionList)
        if iso_list.highlighted is None:
            self.notify("Please select an Isotype.", severity="warning")
            return
        selected_iso = list(SEQUENCES.keys())[iso_list.highlighted]
        self.app.selected_isotype = selected_iso
        self.app.push_screen(AllotypeScreen())

    def action_back(self) -> None:
        self.app.pop_screen()

class AllotypeScreen(Screen):
    BINDINGS = [("enter", "next", "Next"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("STEP 2", classes="title")
            yield Label(f"Select {self.app.selected_isotype.upper()} Allotype", classes="subtitle")
            allotypes = SEQUENCES.get(self.app.selected_isotype, {})
            yield OptionList(
                *[allo.capitalize() for allo in allotypes.keys()],
                id="allo-list"
            )
        yield Footer()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.action_next()

    def action_next(self) -> None:
        allo_list = self.query_one("#allo-list", OptionList)
        if allo_list.highlighted is None:
            self.notify("Please select an Allotype.", severity="warning")
            return
        selected_allo = list(SEQUENCES.get(self.app.selected_isotype, {}).keys())[allo_list.highlighted]
        self.app.selected_allotype = selected_allo
        self.app.push_screen(MutationScreen())

    def action_back(self) -> None:
        self.app.pop_screen()

class MutationScreen(Screen):
    BINDINGS = [("enter", "generate", "Generate"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("STEP 3", classes="title")
            yield Label("Mutations (Optional)", classes="subtitle")
            yield Label("Presets:", classes="subtitle")
            with Horizontal(id="selection-area"):
                yield SelectionList[str](
                    *[Selection(item["label"], item["value"], False) for item in COMMON_MUTATIONS],
                    id="list-common"
                )
                yield Pretty([], id="selected-preview")
            yield Label("Custom (e.g. S239D/I332E):", classes="subtitle")
            yield Input(placeholder="None", id="input-custom", max_length=100, restrict=r"^[a-zA-Z0-9/,\s]*$")
            yield Button("Generate FASTA", variant="primary", id="btn-gen")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#list-common").border_title = "Common Mutants"
        self.query_one("#selected-preview").border_title = "Selected"

    @on(SelectionList.SelectedChanged)
    def update_selected_view(self) -> None:
        self.query_one("#selected-preview", Pretty).update(self.query_one("#list-common", SelectionList).selected)

    def action_generate(self) -> None:
        selected_presets = self.query_one("#list-common", SelectionList).selected
        preset_str = "/".join(selected_presets)
        custom_str = self.query_one("#input-custom", Input).value
        self.app.all_mutants = "/".join(filter(None, [preset_str, custom_str]))
        self.app.push_screen(ResultScreen())

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "btn-gen": self.action_generate()

    def action_back(self) -> None:
        self.app.pop_screen()

class ResultScreen(Screen):
    BINDINGS = [
        ("escape", "back", "Back"),
        ("ctrl+y", "copy_to_clipboard", "Copy"),
        ("q", "quit_to_main", "Main Menu")
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(classes="minimal-container"):
            yield Label("RESULT", classes="title")
            yield Label("FASTA Sequence", classes="subtitle")
            yield Log(id="result-box")
            yield Static("[dim]Ctrl+Y: Copy | Esc: Back | Q: Menu[/]", id="result-help")
        yield Footer()

    def on_mount(self) -> None:
        self.generate_fasta()

    def generate_fasta(self) -> None:
        isotype = self.app.selected_isotype
        allotype = self.app.selected_allotype
        all_mutants = self.app.all_mutants
        
        result_box = self.query_one("#result-box", Log)
        result_box.clear()

        try:
            isotype_data = SEQUENCES.get(isotype, {})
            base_seq = isotype_data.get(allotype)
            if not base_seq:
                result_box.write(f"[bold red]Error: Base sequence for {escape(isotype)} {escape(allotype)} not found.[/]")
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
        except Exception as e:
            result_box.write("[bold red]An unexpected error occurred during sequence generation.[/]")
            self.log.error(f"Error in generate_fasta: {e}", exc_info=True)

    def action_copy_to_clipboard(self) -> None:
        if hasattr(self.app, "last_fasta"):
            pyperclip.copy(self.app.last_fasta)
            self.notify("FASTA sequence copied!")

    def action_quit_to_main(self) -> None:
        while len(self.app.screen_stack) > 1:
            self.app.pop_screen()

    def action_back(self) -> None:
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
        self.push_screen(WelcomeScreen())


def main():
    MutantApp().run()

if __name__ == "__main__":
    main()
