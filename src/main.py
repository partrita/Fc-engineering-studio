import sys
import os
from typing import Dict, List, Optional, Tuple, Set
import yaml

from textual.app import App, ComposeResult
from textual.containers import Container, Vertical, Center, Middle
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
)
from textual.widgets.selection_list import Selection
from textual.binding import Binding
from textual.screen import Screen
import pyperclip

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

# --- Catppuccin Latte Colors & Retro Styling ---

CATPPUCCIN_LATTE = """
$rosewater: #dc8a78;
$flamingo: #dd7878;
$pink: #ea76cb;
$mauve: #8839ef;
$red: #d20f39;
$maroon: #e64553;
$peach: #fe640b;
$yellow: #df8e1d;
$green: #40a02b;
$teal: #179287;
$sky: #04a5e5;
$sapphire: #209fb5;
$blue: #1e66f5;
$lavender: #7287fd;
$text: #4c4f69;
$subtext1: #5c5f77;
$subtext0: #6c6f85;
$overlay2: #7c7f93;
$overlay1: #8c8fa1;
$overlay0: #9ca0b0;
$surface2: #acb0be;
$surface1: #bcc0cc;
$surface0: #ccd0da;
$base: #eff1f5;
$mantle: #e6e9ef;
$crust: #dce0e8;

$primary: $blue;
$secondary: $mauve;
$accent: $peach;
$success: $green;
$error: $red;
$warning: $yellow;
"""

ANTIBODY_ASCII = r"""
     \ /     \ /
      V       V
      |       |
      |---H---|
      |       |
      |   C   |
      |       |
      \_______/
"""

RETRO_STYLE = """
Screen {
    background: $base;
    color: $text;
}

Header {
    background: $mantle;
    color: $blue;
    text-style: bold;
    border-bottom: heavy $blue;
}

Footer {
    background: $mantle;
    color: $subtext0;
}

.title {
    color: $mauve;
    text-style: bold italic underline;
    margin-bottom: 1;
}

.subtitle {
    color: $subtext1;
    margin-bottom: 1;
}

.ascii-art {
    color: $peach;
    margin: 1 0;
}

.retro-box {
    border: double $blue;
    background: $surface0;
    padding: 1 2;
    margin: 1;
    width: 80;
    align: center middle;
    /* transition: opacity 500ms in_out_cubic; */
}

Button {
    background: $surface1;
    color: $text;
    border: tall $blue;
    margin-top: 1;
}

Button:hover {
    background: $blue;
    color: $base;
}

OptionList {
    background: $surface0;
    color: $text;
    border: double $blue;
    height: 6;
}

OptionList > .option-list--highlighted {
    background: $blue;
    color: $base;
}

Log {
    background: $crust;
    color: $green;
    border: double $green;
    height: 10;
}

Input {
    background: $surface0;
    color: $text;
    border: double $lavender;
}

SelectionList {
    background: $surface0;
    color: $text;
    border: double $blue;
    height: 8;
}

#press-enter {
    color: $red;
    text-style: bold;
}
"""

# --- Screens ---

class WelcomeScreen(Screen):
    BINDINGS = [("enter", "next", "Next")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Center():
            with Middle():
                with Vertical(classes="retro-box", id="welcome-container"):
                    yield Label("Fc Engineering Studio Pro", classes="title")
                    yield Label("Professional Antibody Fc Sequence Designer", classes="subtitle")
                    yield Static(ANTIBODY_ASCII, classes="ascii-art")
                    yield Static("[bold]Overview:[/]\nThis tool helps you design mutant sequences of the human IgG Fc region based on [b]EU Numbering[/b].", id="overview")
                    yield Label("\nPress ENTER to start...", id="press-enter")
        yield Footer()

    def action_next(self) -> None:
        self.app.push_screen(IsotypeScreen())

class IsotypeScreen(Screen):
    BINDINGS = [("enter", "next", "Next"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Center():
            with Middle():
                with Vertical(classes="retro-box"):
                    yield Label("STEP 1: Select IgG Isotype", classes="subtitle")
                    yield OptionList(
                        *[iso.upper() for iso in SEQUENCES.keys()],
                        id="iso-list"
                    )
                    yield Label("\n[dim]Use arrows to select, ENTER to continue[/]")
        yield Footer()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.action_next()

    def action_next(self) -> None:
        selected_iso = str(self.query_one("#iso-list").get_option_at_index(self.query_one("#iso-list").highlighted).prompt).lower()
        self.app.selected_isotype = selected_iso
        self.app.push_screen(AllotypeScreen())

    def action_back(self) -> None:
        self.app.pop_screen()

class AllotypeScreen(Screen):
    BINDINGS = [("enter", "next", "Next"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Center():
            with Middle():
                with Vertical(classes="retro-box"):
                    yield Label(f"STEP 2: Select {self.app.selected_isotype.upper()} Allotype", classes="subtitle")
                    allotypes = SEQUENCES.get(self.app.selected_isotype, {})
                    yield OptionList(
                        *[allo.capitalize() for allo in allotypes.keys()],
                        id="allo-list"
                    )
                    yield Label("\n[dim]Use arrows to select, ENTER to continue[/]")
        yield Footer()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.action_next()

    def action_next(self) -> None:
        selected_allo = str(self.query_one("#allo-list").get_option_at_index(self.query_one("#allo-list").highlighted).prompt).lower()
        self.app.selected_allotype = selected_allo
        self.app.push_screen(MutationScreen())

    def action_back(self) -> None:
        self.app.pop_screen()

class MutationScreen(Screen):
    BINDINGS = [("enter", "generate", "Generate"), ("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Center():
            with Middle():
                with Vertical(classes="retro-box"):
                    yield Label("STEP 3: Mutations (Optional)", classes="subtitle")
                    yield Label("Common Mutations:")
                    yield SelectionList[str](
                        *[Selection(item["label"], item["value"], False) for item in COMMON_MUTATIONS],
                        id="list-common"
                    )
                    yield Label("Custom Mutations (e.g. S239D/I332E):")
                    yield Input(placeholder="None", id="input-custom")
                    yield Button("Generate FASTA", variant="primary", id="btn-gen")
        yield Footer()

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
        with Center():
            with Middle():
                with Vertical(classes="retro-box"):
                    yield Label("FINAL RESULT: FASTA Sequence", classes="subtitle")
                    yield Log(id="result-box")
                    yield Static("\n[bold]Shortcuts:[/]\n[b]Ctrl+Y[/]: Copy | [b]Esc[/]: Back | [b]Q[/]: Menu", id="result-help")
        yield Footer()

    def on_mount(self) -> None:
        self.generate_fasta()

    def generate_fasta(self) -> None:
        isotype = self.app.selected_isotype
        allotype = self.app.selected_allotype
        all_mutants = self.app.all_mutants
        
        base_seq = SEQUENCES[isotype][allotype]
        mutant_seq, errors = apply_mutations(base_seq, all_mutants, isotype)
        
        result_box = self.query_one("#result-box", Log)
        result_box.clear()
        
        if errors:
            result_box.write("[bold red]Errors during processing:[/]")
            for err in errors: result_box.write(f"[red]• {err}[/]")
            return

        display_muts = all_mutants.replace("/", "_") if all_mutants else "WT"
        header = f"{isotype.upper()}_{allotype.capitalize()}_{display_muts}"
        fasta = f">{header}\n{mutant_seq}"
        
        result_box.write(fasta)
        self.app.last_fasta = fasta

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
    TITLE = "Fc Engineering Studio Pro"
    CSS = CATPPUCCIN_LATTE + RETRO_STYLE
    
    def on_mount(self) -> None:
        self.selected_isotype = ""
        self.selected_allotype = ""
        self.all_mutants = ""
        self.push_screen(WelcomeScreen())

def main():
    MutantApp().run()

if __name__ == "__main__":
    main()
