import sys
import os
from typing import Dict, List, Optional, Tuple, Set
import yaml

from textual.app import App, ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import (
    Button,
    Footer,
    Header,
    Input,
    Label,
    Log,
    Select,
    Static,
    SelectionList,
)
from textual.widgets.selection_list import Selection
from textual.binding import Binding
import pyperclip

# --- Configuration & Data Loading ---

def load_yaml_data():
    base_path = os.path.dirname(__file__)
    seq_path = os.path.join(base_path, "sequences.yaml")
    mut_path = os.path.join(base_path, "mutants.yaml")
    
    isotypes = {}
    common_muts = []
    
    # Load Isotypes/Sequences
    try:
        if os.path.exists(seq_path):
            with open(seq_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                isotypes = data.get("isotypes", {})
    except Exception as e:
        print(f"Error loading sequences.yaml: {e}")

    # Load Mutations
    try:
        if os.path.exists(mut_path):
            with open(mut_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                common_muts = data.get("common_mutations", [])
    except Exception as e:
        print(f"Error loading mutants.yaml: {e}")
        
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

# --- Textual TUI ---

USAGE_GUIDE = """
[bold cyan]Welcome to Fc Engineering Studio Pro![/]

이 도구는 [b]EU Numbering[/b] 체계를 기준으로 인간 IgG Fc 변이 서열을 생성합니다.

[b]How to use:[/b]
1. [cyan]Isotype & Allotype[/]을 선택하여 베이스 서열을 결정합니다.
2. [cyan]Common Mutations[/] 목록에서 널리 알려진 변이를 체크합니다.
3. [cyan]Custom Mutants[/] 입력창에 추가하고 싶은 변이를 입력합니다.
4. [cyan]Generate[/] 버튼을 누르거나 [b]Enter[/b] 키를 입력하여 결과를 확인합니다.
5. [b]Ctrl + Y[/b]를 눌러 FASTA 서열을 클립보드에 복사합니다.

[dim]* 서열 데이터는 [b]sequences.yaml[/b], 프리셋 변이는 [b]mutants.yaml[/b]에서 관리됩니다.[/]
"""

class MutantApp(App):
    TITLE = "Fc Engineering Studio Pro"
    CSS = """
    Container { padding: 1; }
    .sidebar { width: 40%; border-right: solid $primary; padding-right: 1; }
    .main-content { width: 60%; padding-left: 2; }
    Label { margin-top: 1; color: $accent; font-weight: bold; }
    #guide-box { 
        background: $surface; 
        border: solid $primary-darken-2; 
        padding: 1; 
        margin-bottom: 1;
        color: $text;
    }
    #result-box { 
        height: 1fr; 
        border: panel $primary; 
        padding: 1; 
        background: $surface; 
        color: $text; 
    }
    SelectionList { height: 10; border: sunken $panel; margin-top: 1; }
    .help { margin-top: 1; color: $text-disabled; font-size: 80%; }
    """

    BINDINGS = [
        Binding("ctrl+c", "quit", "Quit"),
        Binding("enter", "generate", "Generate"),
        Binding("ctrl+y", "copy_to_clipboard", "Copy"),
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with Container():
            with Horizontal():
                with Vertical(classes="sidebar"):
                    yield Label("1. Isotype & Allotype")
                    yield Select([(iso.upper(), iso) for iso in SEQUENCES.keys()], value="igg1" if "igg1" in SEQUENCES else None, id="sel-iso")
                    yield Select([], id="sel-allo")
                    
                    yield Label("2. Common Mutations")
                    yield SelectionList[str](
                        [Selection(item["label"], item["value"], False) for item in COMMON_MUTATIONS],
                        id="list-common"
                    )
                    
                    yield Label("3. Custom Mutants")
                    yield Input(placeholder="e.g. S239D/I332E", id="input-custom")
                    
                    yield Button("Generate Sequence", variant="primary", id="btn-gen")
                    yield Static("[b]Shortcuts:[/b] [Enter] Gen | [Ctrl+Y] Copy | [Ctrl+C] Quit", classes="help")
                
                with Vertical(classes="main-content"):
                    yield Static(USAGE_GUIDE, id="guide-box")
                    yield Label("Generated FASTA Sequence")
                    yield Log(id="result-box")
        yield Footer()

    def on_mount(self) -> None:
        if "igg1" in SEQUENCES:
            self.update_allotypes("igg1")

    def update_allotypes(self, isotype: str) -> None:
        allotypes = SEQUENCES.get(isotype, {})
        select = self.query_one("#sel-allo", Select)
        select.set_options([(a.capitalize(), a) for a in allotypes.keys()])
        if "wt" in allotypes: select.value = "wt"
        elif allotypes: select.value = list(allotypes.keys())[0]

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "sel-iso" and event.value:
            self.update_allotypes(str(event.value))

    def action_generate(self) -> None:
        iso_val = self.query_one("#sel-iso", Select).value
        allo_val = self.query_one("#sel-allo", Select).value
        
        if not iso_val or not allo_val:
            self.notify("Please select Isotype and Allotype.", severity="error")
            return
            
        isotype = str(iso_val)
        allotype = str(allo_val)
        selected_presets = self.query_one("#list-common", SelectionList).selected
        preset_str = "/".join(selected_presets)
        custom_str = self.query_one("#input-custom", Input).value
        all_mutants = "/".join(filter(None, [preset_str, custom_str]))
        
        base_seq = SEQUENCES[isotype][allotype]
        mutant_seq, errors = apply_mutations(base_seq, all_mutants, isotype)
        
        result_box = self.query_one("#result-box", Log)
        result_box.clear()
        
        if errors:
            result_box.write("[bold red]Errors found during processing:[/]")
            for err in errors: result_box.write(f"[red]• {err}[/]")
            return

        display_muts = all_mutants.replace("/", "_") if all_mutants else "WT"
        header = f"{isotype.upper()}_{allotype.capitalize()}_{display_muts}"
        fasta = f">{header}\n{mutant_seq}"
        
        result_box.write(fasta)
        self.last_fasta = fasta

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "btn-gen": self.action_generate()

    def action_copy_to_clipboard(self) -> None:
        if hasattr(self, "last_fasta"):
            pyperclip.copy(self.last_fasta)
            self.notify("FASTA sequence copied to clipboard!")
        else:
            self.notify("Nothing to copy yet.", severity="error")

if __name__ == "__main__":
    MutantApp().run()
