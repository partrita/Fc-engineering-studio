"""Dedicated Fc Engineering Workbench tab for Fc Engineering Studio."""

import re
from typing import Any, Callable, Dict, List, Optional, Tuple

import flet as ft

from fc_engineer import config
from fc_engineer.core import CONSTANT_DOMAIN_EU_RANGES, lookup_wt_residue
from fc_engineer.designer import (
    PolypeptideChain,
    analyze_chain_developability,
    assemble_chain_sequence,
)
from fc_engineer.fc_catalog import (
    AMINO_ACIDS,
    FC_CATEGORIES,
    FC_RECIPES,
    build_eu_residue_rows,
    get_hotspot_info,
    get_preset_metadata,
    summarize_fc_profile,
)


def _effect_pill(label: str, state: str, pal: Dict[str, str]) -> ft.Container:
    """Render a compact status badge for ADCC / CDC / FcRn / Heterodimer state."""
    low = state.lower()
    if low == "up":
        bg, border, fg, text = "#064E3B", "#10B981", "#A7F3D0", f"{label}: ↑ Enhanced"
    elif low == "down":
        bg, border, fg, text = "#7C2D12", "#F97316", "#FED7AA", f"{label}: ↓ Reduced"
    elif low == "mixed":
        bg, border, fg, text = "#581C87", "#A855F7", "#E9D5FF", f"{label}: ⇅ Mixed"
    elif low.startswith("low"):
        bg, border, fg, text = "#1E3A5F", "#38BDF8", "#BAE6FD", f"{label}: {state}"
    elif low == "yes":
        bg, border, fg, text = "#312E81", "#818CF8", "#E0E7FF", f"{label}: Active"
    else:
        bg, border, fg, text = pal["elevated"], pal["border"], pal["muted"], f"{label}: WT"

    return ft.Container(
        padding=ft.Padding.symmetric(horizontal=8, vertical=4),
        border_radius=6,
        bgcolor=bg,
        border=ft.Border.all(1, border),
        content=ft.Text(text, size=11, weight=ft.FontWeight.W_600, color=fg),
    )


def _mini_effect_badges(meta: Any) -> List[ft.Control]:
    """Return compact inline effect pills for a preset card."""
    badges: List[ft.Control] = []
    specs = [
        ("ADCC", meta.adcc),
        ("CDC", meta.cdc),
        ("t½", meta.fcrn),
    ]
    for tag, st in specs:
        if st == "up":
            badges.append(
                ft.Container(
                    padding=ft.Padding.symmetric(horizontal=5, vertical=1),
                    border_radius=4,
                    bgcolor="#065F46",
                    content=ft.Text(f"{tag} ↑", size=9, weight=ft.FontWeight.BOLD, color="#D1FAE5"),
                )
            )
        elif st == "down":
            badges.append(
                ft.Container(
                    padding=ft.Padding.symmetric(horizontal=5, vertical=1),
                    border_radius=4,
                    bgcolor="#9A3412",
                    content=ft.Text(f"{tag} ↓", size=9, weight=ft.FontWeight.BOLD, color="#FFEDD5"),
                )
            )
    if meta.heterodimer:
        badges.append(
            ft.Container(
                padding=ft.Padding.symmetric(horizontal=5, vertical=1),
                border_radius=4,
                bgcolor="#3730A3",
                content=ft.Text("Heterodimer", size=9, weight=ft.FontWeight.BOLD, color="#E0E7FF"),
            )
        )
    return badges


def build_fc_engineer_tab(
    page: ft.Page,
    ctrl: Any,
    palette_fn: Callable[[str], Dict[str, str]],
    notify: Callable[..., None],
    render_workspace: Callable[[], None],
    clipboard_countdown_coro: Callable[[], Any],
) -> Tuple[ft.Container, Callable[[], None], Callable[[str], None]]:
    """Build the dedicated Fc Engineering Workbench tab and its refresh/theme hooks."""
    pal = palette_fn(ctrl.theme_mode)

    wb_state: Dict[str, Any] = {
        "chain_id": None,
        "category": "All",
        "query": "",
        "eu_domain": "CH2",
        "hotspots_only": True,
    }

    # Dynamic containers refreshed on every render_workspace()
    header_row = ft.Row(spacing=10, wrap=True, vertical_alignment=ft.CrossAxisAlignment.CENTER)
    recipes_row = ft.Row(spacing=6, wrap=True, run_spacing=6)
    cat_chips_row = ft.Row(spacing=4, wrap=True, run_spacing=4)
    presets_list_col = ft.Column(spacing=6, scroll=ft.ScrollMode.AUTO)
    profile_col = ft.Column(spacing=8)
    eu_domain_row = ft.Row(spacing=6, wrap=True, vertical_alignment=ft.CrossAxisAlignment.CENTER)
    eu_grid_row = ft.Row(spacing=4, wrap=True, run_spacing=4)
    diffs_table_col = ft.Column(spacing=4)
    adhoc_wt_lbl = ft.Text(
        "Tip: type an EU position (e.g. 234, 252, 297, 366, 428) to preview WT residue.",
        size=10,
        color=pal["muted"],
    )

    def _pick_active_chain() -> Optional[PolypeptideChain]:
        chains = ctrl.project.chains
        if not chains:
            wb_state["chain_id"] = None
            return None
        by_id = {c.id: c for c in chains}
        # Prefer explicit workbench selection if still valid
        if wb_state["chain_id"] in by_id:
            return by_id[wb_state["chain_id"]]
        # Next prefer currently selected chain if it has Fc domains
        if ctrl.selected_chain_id in by_id:
            sel = by_id[ctrl.selected_chain_id]
            if any(d.domain_type in ("Hinge", "CH2", "CH3") for d in sel.domains):
                wb_state["chain_id"] = sel.id
                return sel
        # Next prefer first heavy chain with Fc domains
        for c in chains:
            if any(d.domain_type in ("Hinge", "CH2", "CH3") for d in c.domains):
                wb_state["chain_id"] = c.id
                return c
        wb_state["chain_id"] = chains[0].id
        return chains[0]

    def _open_point_mutation_picker(chain: PolypeptideChain, row: Dict[str, Any]) -> None:
        active_pal = palette_fn(ctrl.theme_mode)
        eu_pos = row["eu_pos"]
        dom_id = row["domain_id"]
        dom_type = row["domain_type"]
        wt_aa = row["wt_aa"]
        mut_aa = row["mut_aa"]
        hotspot_desc = row["hotspot_desc"]
        ptm_tags = row["ptm_tags"]

        def _apply_aa(aa: str) -> None:
            errs = ctrl.mutate_residue_at_eu(chain.id, dom_id, eu_pos, aa)
            picker_dlg.open = False
            if errs:
                notify("; ".join(errs), is_error=True)
                page.update()
            else:
                if aa == wt_aa:
                    notify(f"Reverted EU {eu_pos} ({dom_type}) to WT ({wt_aa}).")
                else:
                    scope = "all heavy chains" if ctrl.apply_to_all_heavy else chain.name
                    notify(f"Applied {wt_aa}{eu_pos}{aa} on {scope}.")
                render_workspace()

        aa_buttons: List[ft.Control] = []
        for aa in AMINO_ACIDS:
            is_cur = aa == mut_aa
            is_wt = aa == wt_aa
            label_txt = f"{aa}*" if is_wt else aa
            btn_bg = "#D97706" if (is_cur and not is_wt) else ("#0F766E" if is_cur else active_pal["elevated"])
            btn_fg = "#FFFFFF" if is_cur else active_pal["text"]
            aa_buttons.append(
                ft.Container(
                    width=44,
                    height=34,
                    border_radius=6,
                    bgcolor=btn_bg,
                    border=ft.Border.all(1.5 if is_cur else 1.0, "#F59E0B" if is_cur else active_pal["border"]),
                    alignment=ft.Alignment.CENTER,
                    tooltip=f"{'WT residue' if is_wt else f'Mutate to {wt_aa}{eu_pos}{aa}'}",
                    on_click=lambda ev, target_aa=aa: _apply_aa(target_aa),
                    content=ft.Text(label_txt, size=12, weight=ft.FontWeight.BOLD, color=btn_fg),
                )
            )

        info_lines: List[ft.Control] = [
            ft.Text(
                f"Chain: {chain.name} · Domain: {dom_type} ({row['isotype'].upper()}) · WT: {wt_aa} · Current: {mut_aa}",
                size=12,
                weight=ft.FontWeight.W_600,
            ),
        ]
        if hotspot_desc:
            info_lines.append(
                ft.Text(
                    f"Functional Site [{row['hotspot_tag']}]: {hotspot_desc}",
                    size=11,
                    color="#0EA5E9",
                )
            )
        if ptm_tags:
            info_lines.append(
                ft.Text(
                    f"PTM Hotspot: {', '.join(ptm_tags)}",
                    size=11,
                    color="#F59E0B",
                )
            )
        info_lines.append(
            ft.Text(
                "Click an amino acid below (* = WT). "
                + ("Applies to all heavy chains." if ctrl.apply_to_all_heavy else f"Applies to {chain.name}."),
                size=10,
                color=active_pal["muted"],
            )
        )

        picker_dlg = ft.AlertDialog(
            title=ft.Text(f"Point Mutation — EU {eu_pos} ({wt_aa} → {mut_aa})"),
            content=ft.Container(
                width=420,
                content=ft.Column(
                    tight=True,
                    spacing=8,
                    controls=[
                        *info_lines,
                        ft.Row(wrap=True, spacing=6, run_spacing=6, controls=aa_buttons),
                    ],
                ),
            ),
            actions=[
                *(
                    [
                        ft.TextButton(
                            f"Revert to WT ({wt_aa})",
                            on_click=lambda ev: _apply_aa(wt_aa),
                        )
                    ]
                    if mut_aa != wt_aa
                    else []
                ),
                ft.Button(
                    "Cancel",
                    on_click=lambda ev: (setattr(picker_dlg, "open", False), page.update()),
                ),
            ],
        )
        if hasattr(page, "show_dialog"):
            page.show_dialog(picker_dlg)

    def _on_search_change(ev: ft.ControlEvent) -> None:
        wb_state["query"] = ev.control.value or ""
        _refresh_preset_catalog()
        page.update()

    search_tf = ft.TextField(
        label="Search Fc Presets (name, EU position, or effect)",
        value="",
        dense=True,
        expand=True,
        prefix_icon=ft.Icons.SEARCH,
        on_change=_on_search_change,
    )

    def _on_adhoc_input_change(ev: ft.ControlEvent) -> None:
        chain = _pick_active_chain()
        raw = (ev.control.value or "").strip()
        m = re.search(r"(\d{3})", raw)
        if not m or chain is None:
            adhoc_wt_lbl.value = "Tip: type an EU position (e.g. 234, 252, 297, 366, 428) to preview WT residue."
            page.update()
            return
        pos = int(m.group(1))
        iso_dom = next(
            (d for d in chain.domains if d.domain_type in CONSTANT_DOMAIN_EU_RANGES),
            None,
        )
        iso = iso_dom.isotype if iso_dom is not None else "igg1"
        allo = iso_dom.allotype if iso_dom is not None else ""
        iso_seqs = config.SEQUENCES.get(iso, {})
        full_seq = iso_seqs.get(allo) or (next(iter(iso_seqs.values())) if iso_seqs else "")
        aa = lookup_wt_residue(full_seq, pos, iso) if full_seq else None
        hs = get_hotspot_info(pos)
        if aa:
            hs_note = f" [{hs[0]}: {hs[1]}]" if hs else ""
            adhoc_wt_lbl.value = f"WT @ EU {pos} ({iso.upper()}): {aa} — e.g. {aa}{pos}A{hs_note}"
        else:
            adhoc_wt_lbl.value = f"EU {pos} not present in {iso.upper()} constant region."
        page.update()

    adhoc_tf = ft.TextField(
        label="Ad-hoc EU Mutation (e.g. L234A/L235A or M252Y)",
        dense=True,
        expand=True,
        max_length=1000,
        input_filter=ft.InputFilter(allow=True, regex_string=r"[a-zA-Z0-9/, ]"),
        on_change=_on_adhoc_input_change,
    )

    def _on_apply_adhoc(ev: ft.ControlEvent) -> None:
        chain = _pick_active_chain()
        if chain is None:
            notify("No chain available.", is_error=True)
            return
        errs = ctrl.apply_custom_mutation(chain.id, adhoc_tf.value or "")
        if errs:
            notify("; ".join(errs), is_error=True)
        else:
            adhoc_tf.value = ""
            notify("Applied custom EU mutation.")
            render_workspace()

    save_preset_label_tf = ft.TextField(
        label="New Project Preset Name",
        dense=True,
        expand=True,
        max_length=40,
        input_filter=ft.InputFilter(allow=False, regex_string=r"[\r\n]"),
    )
    save_preset_val_tf = ft.TextField(
        label="EU Mutations (e.g. L234A/L235A)",
        dense=True,
        expand=True,
        max_length=1000,
        input_filter=ft.InputFilter(allow=True, regex_string=r"[a-zA-Z0-9/, ]"),
    )

    def _on_save_project_preset(ev: ft.ControlEvent) -> None:
        errs = ctrl.register_file_custom_mutation(
            save_preset_label_tf.value or "",
            save_preset_val_tf.value or "",
        )
        if errs:
            notify("; ".join(errs), is_error=True)
        else:
            save_preset_label_tf.value = ""
            save_preset_val_tf.value = ""
            notify("Saved custom preset to current .ab project.")
            render_workspace()

    def _refresh_preset_catalog() -> None:
        active_pal = palette_fn(ctrl.theme_mode)
        chain = _pick_active_chain()
        if chain is None:
            presets_list_col.controls = [
                ft.Text("Add a heavy chain to engineer Fc mutations.", size=12, color=active_pal["muted"])
            ]
            return

        file_labels = {
            m.get("label")
            for m in ctrl.project.file_custom_mutations
            if isinstance(m, dict)
        }
        active_presets = {p for d in chain.domains for p in d.applied_presets}
        preset_map = ctrl.get_preset_map()
        q = wb_state["query"].strip().lower()
        cat = wb_state["category"]

        # Category filter chips
        chips: List[ft.Control] = []
        for cname in FC_CATEGORIES:
            is_sel = cat == cname

            def _make_cat_cb(c_val: str):
                def _on_cat(ev: ft.ControlEvent) -> None:
                    wb_state["category"] = c_val
                    _refresh_preset_catalog()
                    page.update()

                return _on_cat

            chips.append(
                ft.Chip(
                    label=ft.Text(cname, size=10),
                    selected=is_sel,
                    on_select=_make_cat_cb(cname),
                )
            )
        cat_chips_row.controls = chips

        cards: List[ft.Control] = []
        for lbl, val in preset_map.items():
            is_proj = lbl in file_labels
            meta = get_preset_metadata(lbl, val, is_project_custom=is_proj)
            if cat != "All" and meta.category != cat:
                continue
            if q and not (
                q in lbl.lower()
                or q in val.lower()
                or q in meta.summary.lower()
                or q in meta.category.lower()
            ):
                continue

            is_active = lbl in active_presets

            def _make_toggle(preset_lbl: str, preset_val: str):
                def _on_toggle(ev: ft.ControlEvent) -> None:
                    errs = ctrl.toggle_preset(chain.id, preset_lbl, bool(ev.control.value))
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    else:
                        action = "Applied" if ev.control.value else "Removed"
                        notify(f"{action} {preset_lbl} ({preset_val}).")
                    render_workspace()

                return _on_toggle

            right_controls: List[ft.Control] = _mini_effect_badges(meta)
            if is_proj:
                def _make_del_proj(plbl: str):
                    def _on_del(ev: ft.ControlEvent) -> None:
                        ctrl.delete_file_custom_mutation(plbl)
                        notify(f"Deleted project preset '{plbl}'.")
                        render_workspace()

                    return _on_del

                right_controls.append(
                    ft.IconButton(
                        icon=ft.Icons.DELETE_OUTLINE,
                        icon_size=15,
                        tooltip=f"Delete project preset {lbl}",
                        on_click=_make_del_proj(lbl),
                    )
                )

            cards.append(
                ft.Container(
                    padding=ft.Padding.symmetric(horizontal=8, vertical=6),
                    border_radius=6,
                    bgcolor=active_pal["elevated"],
                    border=ft.Border.all(
                        1.5 if is_active else 1.0,
                        "#0F766E" if is_active else active_pal["border"],
                    ),
                    content=ft.Row(
                        spacing=6,
                        vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        controls=[
                            ft.Checkbox(
                                value=is_active,
                                on_change=_make_toggle(lbl, val),
                            ),
                            ft.Column(
                                spacing=2,
                                expand=True,
                                controls=[
                                    ft.Row(
                                        spacing=6,
                                        wrap=True,
                                        controls=[
                                            ft.Text(
                                                lbl,
                                                size=12,
                                                weight=ft.FontWeight.BOLD,
                                                color=active_pal["text"],
                                            ),
                                            ft.Text(
                                                val,
                                                size=10,
                                                font_family="monospace",
                                                color="#0EA5E9",
                                            ),
                                            *right_controls,
                                        ],
                                    ),
                                    ft.Text(
                                        meta.summary,
                                        size=10,
                                        color=active_pal["muted"],
                                    ),
                                ],
                            ),
                        ],
                    ),
                )
            )

        if not cards:
            cards = [
                ft.Text(
                    "No Fc presets match the current filter.",
                    size=11,
                    color=active_pal["muted"],
                )
            ]
        presets_list_col.controls = cards

    def _refresh_eu_browser(chain: PolypeptideChain) -> None:
        active_pal = palette_fn(ctrl.theme_mode)
        preset_map = ctrl.get_preset_map()
        avail_doms = [
            d.domain_type
            for d in chain.domains
            if d.domain_type in CONSTANT_DOMAIN_EU_RANGES
        ]
        # Preserve canonical Fc-first ordering
        ordered_doms = [dt for dt in ("Hinge", "CH2", "CH3", "CH1") if dt in avail_doms]
        if not ordered_doms:
            eu_domain_row.controls = []
            eu_grid_row.controls = [
                ft.Text(
                    "Selected chain has no constant domains (CH1/Hinge/CH2/CH3).",
                    size=11,
                    color=active_pal["muted"],
                )
            ]
            return

        if wb_state["eu_domain"] not in ordered_doms:
            wb_state["eu_domain"] = "CH2" if "CH2" in ordered_doms else ordered_doms[0]

        dom_chips: List[ft.Control] = []
        for dt in ordered_doms:
            is_sel = wb_state["eu_domain"] == dt
            lo, hi = CONSTANT_DOMAIN_EU_RANGES[dt]

            def _make_dom_sel(d_val: str):
                def _on_sel(ev: ft.ControlEvent) -> None:
                    wb_state["eu_domain"] = d_val
                    _refresh_eu_browser(chain)
                    page.update()

                return _on_sel

            dom_chips.append(
                ft.Chip(
                    label=ft.Text(f"{dt} ({lo}–{hi})", size=10),
                    selected=is_sel,
                    on_select=_make_dom_sel(dt),
                )
            )

        def _on_hotspots_toggle(ev: ft.ControlEvent) -> None:
            wb_state["hotspots_only"] = bool(ev.control.value)
            _refresh_eu_browser(chain)
            page.update()

        eu_domain_row.controls = [
            *dom_chips,
            ft.Switch(
                label="Hotspots Only",
                value=bool(wb_state["hotspots_only"]),
                scale=0.82,
                on_change=_on_hotspots_toggle,
            ),
        ]

        rows = build_eu_residue_rows(chain, preset_map, domain_types=(wb_state["eu_domain"],))
        if wb_state["hotspots_only"]:
            filtered_rows = [
                r for r in rows if r["is_mutated"] or r["hotspot_tag"] or r["ptm_tags"]
            ]
            if not filtered_rows:
                filtered_rows = rows
        else:
            filtered_rows = rows

        tiles: List[ft.Control] = []
        for r in filtered_rows:
            is_mut = r["is_mutated"]
            hs_tag = r["hotspot_tag"]
            has_ptm = bool(r["ptm_tags"])

            if is_mut:
                bg = "#9A3412"
                bdr = "#F59E0B"
                aa_col = "#FEF3C7"
                pos_col = "#FDE68A"
                aa_text = f"{r['wt_aa']}→{r['mut_aa']}"
            elif hs_tag:
                bg = active_pal["elevated"]
                bdr = "#0EA5E9"
                aa_col = active_pal["text"]
                pos_col = "#0EA5E9"
                aa_text = r["wt_aa"]
            elif has_ptm:
                bg = active_pal["elevated"]
                bdr = "#F59E0B"
                aa_col = active_pal["text"]
                pos_col = "#D97706"
                aa_text = r["wt_aa"]
            else:
                bg = active_pal["bg"]
                bdr = active_pal["border"]
                aa_col = active_pal["text"]
                pos_col = active_pal["muted"]
                aa_text = r["wt_aa"]

            tip_parts = [f"EU {r['eu_pos']} ({r['domain_type']}): WT={r['wt_aa']}, Current={r['mut_aa']}"]
            if r["hotspot_desc"]:
                tip_parts.append(f"[{hs_tag}] {r['hotspot_desc']}")
            if r["ptm_tags"]:
                tip_parts.append(f"PTM: {', '.join(r['ptm_tags'])}")
            tip_parts.append("Click to mutate")

            tile_controls: List[ft.Control] = [
                ft.Text(str(r["eu_pos"]), size=8, color=pos_col, weight=ft.FontWeight.BOLD),
                ft.Text(aa_text, size=11, color=aa_col, weight=ft.FontWeight.BOLD, font_family="monospace"),
            ]
            if hs_tag and wb_state["hotspots_only"]:
                tile_controls.append(
                    ft.Text(hs_tag, size=7, color=pos_col, max_lines=1)
                )

            tiles.append(
                ft.Container(
                    width=52 if (is_mut or (hs_tag and wb_state["hotspots_only"])) else 38,
                    padding=ft.Padding.symmetric(horizontal=3, vertical=3),
                    border_radius=5,
                    bgcolor=bg,
                    border=ft.Border.all(1.5 if is_mut else 1.0, bdr),
                    alignment=ft.Alignment.CENTER,
                    tooltip=" | ".join(tip_parts),
                    on_click=lambda ev, ch=chain, row_data=r: _open_point_mutation_picker(ch, row_data),
                    content=ft.Column(
                        spacing=0,
                        tight=True,
                        horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                        controls=tile_controls,
                    ),
                )
            )
        eu_grid_row.controls = tiles

    def rebuild_fc_engineer() -> None:
        active_pal = palette_fn(ctrl.theme_mode)
        chain = _pick_active_chain()
        if chain is None:
            header_row.controls = [
                ft.Text("No chains in project. Add a chain in Studio or Project tab.", size=12)
            ]
            recipes_row.controls = []
            profile_col.controls = []
            eu_domain_row.controls = []
            eu_grid_row.controls = []
            diffs_table_col.controls = []
            _refresh_preset_catalog()
            return

        preset_map = ctrl.get_preset_map()
        file_labels = {
            m.get("label")
            for m in ctrl.project.file_custom_mutations
            if isinstance(m, dict)
        }
        profile = summarize_fc_profile(chain, preset_map, file_preset_labels=file_labels)

        # 1. Header bar controls
        chain_opts = []
        for c in ctrl.project.chains:
            has_fc_tag = " (Fc)" if any(d.domain_type in ("Hinge", "CH2", "CH3") for d in c.domains) else ""
            chain_opts.append(ft.dropdown.Option(c.id, text=f"{c.name}{has_fc_tag}"))

        def _on_chain_select(ev: ft.ControlEvent) -> None:
            wb_state["chain_id"] = ev.control.value
            if ev.control.value:
                ctrl.select_domain(ev.control.value)
            render_workspace()

        def _on_all_heavy_toggle(ev: ft.ControlEvent) -> None:
            ctrl.apply_to_all_heavy = bool(ev.control.value)
            render_workspace()

        cur_iso = profile["isotype"] or "igg1"
        allo_map = config.SEQUENCES.get(cur_iso, {})
        cur_allo = profile["allotype"] or (next(iter(allo_map)) if allo_map else "")

        def _on_iso_change(ev: ft.ControlEvent) -> None:
            new_iso = ev.control.value
            if not new_iso or new_iso == cur_iso:
                return
            errs = ctrl.set_chain_isotype(chain.id, new_iso, clear_incompatible_mutations=False)
            if not errs:
                hint = " Tip: Apply S228P for IgG4 hinge stability." if new_iso == "igg4" else ""
                notify(f"Switched constant domains to {new_iso.upper()}.{hint}")
                render_workspace()
                return

            if not hasattr(page, "show_dialog"):
                notify("; ".join(errs), is_error=True)
                return

            def _force_switch(ev2: ft.ControlEvent) -> None:
                iso_dlg.open = False
                errs2 = ctrl.set_chain_isotype(chain.id, new_iso, clear_incompatible_mutations=True)
                if errs2:
                    notify("; ".join(errs2), is_error=True)
                    page.update()
                else:
                    notify(f"Cleared incompatible mutations and switched to {new_iso.upper()}.")
                    render_workspace()

            iso_dlg = ft.AlertDialog(
                title=ft.Text(f"Switch to {new_iso.upper()}?"),
                content=ft.Text(
                    f"Existing mutations conflict with {new_iso.upper()} WT residues:\n"
                    + "\n".join(errs[:4])
                    + "\n\nClear incompatible mutations and switch?"
                ),
                actions=[
                    ft.TextButton(
                        "Cancel",
                        on_click=lambda ev2: (setattr(iso_dlg, "open", False), render_workspace()),
                    ),
                    ft.Button("Clear & Switch", on_click=_force_switch),
                ],
            )
            page.show_dialog(iso_dlg)

        def _on_allo_change(ev: ft.ControlEvent) -> None:
            new_allo = ev.control.value
            if not new_allo:
                return
            errs = ctrl.set_chain_isotype(chain.id, cur_iso, allotype=new_allo)
            if errs:
                notify("; ".join(errs), is_error=True)
            else:
                notify(f"Updated allotype to {new_allo}.")
            render_workspace()

        def _on_reset_wt(ev: ft.ControlEvent) -> None:
            targets = ctrl._resolve_mutation_target_chains(chain.id, None)
            for tc in targets:
                ctrl.reset_chain_mutations(tc.id)
            render_workspace()
            notify("Reset Fc mutations to WT.")

        def _on_copy_chain(ev: ft.ControlEvent) -> None:
            cseq, _, cerrs = assemble_chain_sequence(chain, preset_map)
            if cerrs:
                notify("; ".join(cerrs), is_error=True)
                return
            if ctrl.copy_text_secure(cseq, label=f"{chain.name} sequence"):
                render_workspace()
                if hasattr(page, "run_task"):
                    page.run_task(clipboard_countdown_coro)
                notify(f"Copied {chain.name} sequence (auto-clears in 30s).")

        header_row.controls = [
            ft.Dropdown(
                label="Target Chain",
                value=chain.id,
                options=chain_opts,
                dense=True,
                width=180,
                on_select=_on_chain_select,
            ),
            ft.Dropdown(
                label="Whole-Chain Isotype",
                value=cur_iso,
                options=[ft.dropdown.Option(k, text=k.upper()) for k in config.SEQUENCES.keys()],
                dense=True,
                width=145,
                disabled=not bool(profile["isotype"]),
                on_select=_on_iso_change,
            ),
            ft.Dropdown(
                label="Allotype",
                value=cur_allo,
                options=[ft.dropdown.Option(ak) for ak in allo_map.keys()],
                dense=True,
                width=190,
                disabled=not bool(allo_map),
                on_select=_on_allo_change,
            ),
            ft.Switch(
                label="All Heavy Chains",
                value=ctrl.apply_to_all_heavy,
                on_change=_on_all_heavy_toggle,
            ),
            ft.OutlinedButton("Reset WT", icon=ft.Icons.RESTART_ALT, on_click=_on_reset_wt),
            ft.OutlinedButton("Copy Seq", icon=ft.Icons.COPY, on_click=_on_copy_chain),
        ]

        # 2. Quick Recipes strip
        recipe_btns: List[ft.Control] = []
        for rec in FC_RECIPES:
            def _make_recipe_cb(recipe: Dict[str, Any]):
                def _on_rec(ev: ft.ControlEvent) -> None:
                    errs = ctrl.apply_fc_recipe(
                        chain.id,
                        presets=recipe["presets"],
                        isotype=recipe.get("isotype"),
                        reset_first=True,
                    )
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    else:
                        notify(f"Applied recipe: {recipe['title']}")
                        render_workspace()

                return _on_rec

            recipe_btns.append(
                ft.OutlinedButton(
                    content=f"{rec['title']}",
                    tooltip=f"{rec['subtitle']} ({rec['badge']})",
                    on_click=_make_recipe_cb(rec),
                )
            )
        recipes_row.controls = recipe_btns

        # 3. Preset catalog
        _refresh_preset_catalog()

        # 4. Live Fc Profile Dashboard
        dmw = profile["delta_mw"]
        dpi = profile["delta_pi"]
        dmw_str = f" ({dmw:+.2f})" if abs(dmw) >= 0.005 else ""
        dpi_str = f" ({dpi:+.2f})" if abs(dpi) >= 0.01 else ""

        liabilities, _ = analyze_chain_developability(
            chain, preset_map, include_variable=ctrl.include_variable_dev
        )
        ptm_badges = [
            ft.Container(
                padding=ft.Padding.symmetric(horizontal=6, vertical=2),
                border_radius=4,
                bgcolor="#78350F" if pos_list else active_pal["elevated"],
                border=ft.Border.all(1, "#F59E0B" if pos_list else active_pal["border"]),
                tooltip=", ".join(str(p) for p in pos_list) or "None",
                content=ft.Text(
                    f"{cat.split('(')[0].strip()}: {len(pos_list)}",
                    size=10,
                    color="#FDE68A" if pos_list else active_pal["muted"],
                ),
            )
            for cat, pos_list in liabilities
        ]

        profile_col.controls = [
            ft.Row(
                spacing=6,
                wrap=True,
                controls=[
                    _effect_pill("ADCC", profile["adcc"], active_pal),
                    _effect_pill("CDC", profile["cdc"], active_pal),
                    _effect_pill("FcRn t½", profile["fcrn"], active_pal),
                    _effect_pill("Heterodimer", "yes" if profile["heterodimer"] else "neutral", active_pal),
                ],
            ),
            ft.Text(
                f"Isotype: {(profile['isotype'] or 'N/A').upper()} · "
                f"Length: {profile['length']} aa · "
                f"MW: {profile['mut_mw']:.2f} kDa{dmw_str} · "
                f"pI: {profile['mut_pi']:.2f}{dpi_str} · "
                f"Mutations: {len(profile['diffs'])}",
                size=11,
                weight=ft.FontWeight.W_600,
                color=active_pal["text"],
            ),
            ft.Row(spacing=4, wrap=True, controls=ptm_badges),
        ]

        # 5. Interactive EU Residue Browser
        _refresh_eu_browser(chain)

        # 6. WT vs. Mutant Diff Table
        diff_rows: List[ft.Control] = []
        dom_id_by_type = {d.domain_type: d.id for d in chain.domains}
        for dt, eu_pos, wt_aa, mut_aa in profile["diffs"]:
            hs = get_hotspot_info(eu_pos)
            hs_text = f"[{hs[0]}] {hs[1]}" if hs else "Constant region"
            did = dom_id_by_type.get(dt, "")

            def _make_revert(cid: str, domain_id: str, pos: int, orig_wt: str):
                def _on_rev(ev: ft.ControlEvent) -> None:
                    errs = ctrl.mutate_residue_at_eu(cid, domain_id, pos, orig_wt)
                    if errs:
                        notify("; ".join(errs), is_error=True)
                    else:
                        notify(f"Reverted EU {pos} to WT ({orig_wt}).")
                        render_workspace()

                return _on_rev

            diff_rows.append(
                ft.Container(
                    padding=ft.Padding.symmetric(horizontal=8, vertical=4),
                    border_radius=5,
                    bgcolor=active_pal["elevated"],
                    border=ft.Border.all(1, active_pal["border"]),
                    content=ft.Row(
                        spacing=8,
                        vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        controls=[
                            ft.Text(dt, size=10, weight=ft.FontWeight.BOLD, width=42, color="#0EA5E9"),
                            ft.Text(
                                f"{wt_aa}{eu_pos}{mut_aa}",
                                size=11,
                                weight=ft.FontWeight.BOLD,
                                font_family="monospace",
                                width=64,
                                color="#F59E0B",
                            ),
                            ft.Text(
                                hs_text,
                                size=10,
                                color=active_pal["muted"],
                                expand=True,
                                overflow=ft.TextOverflow.ELLIPSIS,
                            ),
                            ft.IconButton(
                                icon=ft.Icons.UNDO,
                                icon_size=14,
                                tooltip=f"Revert EU {eu_pos} to WT ({wt_aa})",
                                disabled=not bool(did),
                                on_click=_make_revert(chain.id, did, eu_pos, wt_aa),
                            ),
                        ],
                    ),
                )
            )

        if not diff_rows:
            diff_rows = [
                ft.Text(
                    "Wild-Type (no Fc mutations applied). Choose a Quick Recipe, check a Preset, or click any EU residue above.",
                    size=11,
                    color=active_pal["muted"],
                )
            ]
        diffs_table_col.controls = diff_rows

    # Static card shells (themed dynamically via apply_fc_theme)
    themed_cards: List[ft.Container] = []

    def _make_card(title: str, subtitle: str, body: List[ft.Control], expand: bool = False) -> ft.Container:
        card = ft.Container(
            expand=expand,
            padding=12,
            bgcolor=pal["panel"],
            border=ft.Border.all(1, pal["border"]),
            border_radius=8,
            content=ft.Column(
                spacing=8,
                expand=expand,
                controls=[
                    ft.Column(
                        spacing=2,
                        controls=[
                            ft.Text(title, size=14, weight=ft.FontWeight.BOLD),
                            *(
                                [ft.Text(subtitle, size=10, color=pal["muted"])]
                                if subtitle
                                else []
                            ),
                        ],
                    ),
                    *body,
                ],
            ),
        )
        themed_cards.append(card)
        return card

    top_bar_card = ft.Container(
        padding=10,
        bgcolor=pal["panel"],
        border=ft.Border.all(1, pal["border"]),
        border_radius=8,
        content=header_row,
    )
    themed_cards.append(top_bar_card)

    recipes_card = _make_card(
        "Quick Fc Engineering Recipes",
        "One-click goal-oriented presets (resets existing mutations on target chain(s) atomically).",
        [recipes_row],
    )

    catalog_card = _make_card(
        "Fc Mutation Presets & Custom EU Mutations",
        "Filter by engineering goal, combine presets, or apply custom EU point mutations.",
        [
            search_tf,
            cat_chips_row,
            ft.Container(
                height=260,
                border=ft.Border.all(1, pal["border"]),
                border_radius=6,
                padding=6,
                content=presets_list_col,
            ),
            ft.Row(
                spacing=6,
                controls=[
                    adhoc_tf,
                    ft.Button("Apply", on_click=_on_apply_adhoc),
                ],
            ),
            adhoc_wt_lbl,
            ft.Row(
                spacing=6,
                controls=[
                    save_preset_label_tf,
                    save_preset_val_tf,
                    ft.OutlinedButton("Save Preset", on_click=_on_save_project_preset),
                ],
            ),
        ],
    )

    profile_card = _make_card(
        "Live Fc Functional Profile & Developability",
        "Real-time effector/PK summary, MW/pI delta vs WT, and PTM hotspot counts.",
        [profile_col],
    )

    eu_browser_card = _make_card(
        "Interactive EU-Numbered Sequence Browser",
        "Click any EU residue tile to mutate or revert in 1 click. Blue = functional hotspot, Amber = mutated.",
        [
            eu_domain_row,
            ft.Container(
                padding=8,
                border=ft.Border.all(1, pal["border"]),
                border_radius=6,
                content=eu_grid_row,
            ),
        ],
    )

    diffs_card = _make_card(
        "Engineered EU Substitutions (WT → Mutant)",
        "Active substitutions mapped to EU positions and Fc functional sites.",
        [diffs_table_col],
    )

    fc_engineer_panel = ft.Container(
        expand=True,
        padding=10,
        bgcolor=pal["bg"],
        content=ft.Column(
            expand=True,
            spacing=10,
            scroll=ft.ScrollMode.AUTO,
            controls=[
                top_bar_card,
                ft.Row(
                    spacing=12,
                    vertical_alignment=ft.CrossAxisAlignment.START,
                    controls=[
                        ft.Column(
                            expand=5,
                            spacing=10,
                            controls=[recipes_card, catalog_card],
                        ),
                        ft.Column(
                            expand=6,
                            spacing=10,
                            controls=[profile_card, eu_browser_card, diffs_card],
                        ),
                    ],
                ),
            ],
        ),
    )

    def apply_fc_theme(new_mode: str) -> None:
        new_pal = palette_fn(new_mode)
        fc_engineer_panel.bgcolor = new_pal["bg"]
        for c in themed_cards:
            c.bgcolor = new_pal["panel"]
            c.border = ft.Border.all(1, new_pal["border"])
        adhoc_wt_lbl.color = new_pal["muted"]

    return fc_engineer_panel, rebuild_fc_engineer, apply_fc_theme
