"""Export and persistence formatters (.ab project files, FASTA, GenBank, CSV, transparent PNG)."""

import csv
import io
import json
import math
import os
import re
from datetime import date
from typing import Any, Dict, List, Optional, Tuple

from PIL import Image, ImageDraw, ImageFont

from fc_engineer import config
from fc_engineer.core import (
    CONSTANT_DOMAIN_EU_RANGES,
    diff_sequences,
    get_domain_eu_position,
    get_domain_residue_index,
    get_eu_position,
    parse_mutation_list,
)
from fc_engineer.designer import (
    DEFAULT_SCHEMATIC_THEME,
    MAX_CHAINS,
    MAX_COORD,
    MAX_CUSTOM_SEQ_LENGTH,
    MAX_DOMAINS_PER_CHAIN,
    MAX_EXPORT_SEQUENCE_LENGTH,
    MAX_FILE_CUSTOM_MUTATIONS,
    MAX_HISTORY_ENTRIES,
    MAX_METADATA_LENGTH,
    MAX_PAIRINGS,
    NON_CONSTANT_DOMAIN_TYPES,
    VALID_DOMAIN_TYPES,
    VALID_FORMAT_TYPES,
    VALID_PAIRING_TYPES,
    ALLOWED_TEMPLATE_KEYS_BY_TYPE,
    AntibodyProject,
    DomainBlock,
    InterChainPairing,
    PolypeptideChain,
    SchematicTheme,
    assemble_chain_sequence,
    build_preset_map,
    build_schematic_scene,
    compute_hexagon_vertices,
    get_effective_domain_mutations,
    resolve_domain_base_sequence,
    resolve_domain_sequence,
)

MAX_PNG_DIMENSION = 8192

_MANIFEST_BEGIN = "FC_STUDIO_PROJECT_V1_BEGIN"
_MANIFEST_END = "FC_STUDIO_PROJECT_V1_END"
_VALID_ID_RE = re.compile(r"[A-Za-z0-9_]+\Z")
_VALID_SEQ_RE = re.compile(r"[A-Z]+\Z")
_VALID_OPT_SEQ_RE = re.compile(r"[A-Z]*\Z")
_MUTATION_TOKEN_RE = re.compile(r"[A-Z][0-9]+[A-Z]\Z")

_MONTHS = [
    "JAN",
    "FEB",
    "MAR",
    "APR",
    "MAY",
    "JUN",
    "JUL",
    "AUG",
    "SEP",
    "OCT",
    "NOV",
    "DEC",
]


def sanitize_locus(locus: str) -> str:
    """Reduce a locus name to GenBank-safe characters ([A-Za-z0-9_.-]), max 16 chars."""
    if not isinstance(locus, str):
        return "UNKNOWN"
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", locus)
    return safe[:16] if safe else "UNKNOWN"


def _genbank_date() -> str:
    today = date.today()
    return f"{today.day:02d}-{_MONTHS[today.month - 1]}-{today.year}"


def _format_origin_block(sequence: str) -> List[str]:
    """Format an uppercase amino acid sequence into GenBank ORIGIN lines through '//'."""
    lines: List[str] = ["ORIGIN"]
    # 60 residues per line in 10-residue blocks, position right-justified (width 9)
    for start in range(0, len(sequence), 60):
        chunk = sequence[start : start + 60]
        blocks = " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10))
        lines.append(f"{start + 1:>9} {blocks}")
    lines.append("//")
    return lines


def _sanitize_header_text(val: Any, fallback: str) -> str:
    if not isinstance(val, str):
        return fallback
    cleaned = re.sub(r"[\x00-\x1f\x7f]+", " ", val).strip()
    return cleaned[:MAX_METADATA_LENGTH] if cleaned else fallback


def _validate_human_name(val: Any, field_name: str) -> str:
    if not isinstance(val, str):
        raise ValueError(f"{field_name} must be a string.")
    if any(ord(ch) < 0x20 or ord(ch) == 0x7F for ch in val):
        raise ValueError(f"{field_name} contains disallowed control characters.")
    stripped = val.strip()
    if not stripped:
        raise ValueError(f"{field_name} cannot be empty.")
    if len(stripped) > MAX_METADATA_LENGTH:
        raise ValueError(
            f"{field_name} exceeds maximum length of {MAX_METADATA_LENGTH} characters."
        )
    return stripped


def _sanitize_qualifier(val: str) -> str:
    """Strip quotes and ASCII control characters deterministically for GenBank qualifiers."""
    cleaned = re.sub(r'["\x00-\x1f\x7f]+', "", val)
    return cleaned[:MAX_METADATA_LENGTH]


def _validate_file_path(
    path: str, allowed_exts: Tuple[str, ...], auto_append: bool = True
) -> str:
    if not isinstance(path, str) or not path.strip():
        raise ValueError("File path cannot be empty.")
    _, ext = os.path.splitext(path)
    resolved = f"{path}{allowed_exts[0]}" if (not ext and auto_append) else path
    _, resolved_ext = os.path.splitext(resolved)
    if resolved_ext.lower() not in allowed_exts:
        raise ValueError(
            f"Unsupported extension for '{os.path.basename(resolved)}' (expected {', '.join(allowed_exts)})."
        )
    return resolved


def format_genbank(locus: str, sequence: str, definition: str = "") -> str:
    """Build a minimal GenBank-style protein record. Returns '' on invalid input."""
    if not isinstance(sequence, str) or not sequence:
        return ""
    if len(sequence) > MAX_EXPORT_SEQUENCE_LENGTH:
        return ""
    if not _VALID_SEQ_RE.fullmatch(sequence):
        return ""
    definition = (definition or "")[:MAX_METADATA_LENGTH]
    safe_locus = sanitize_locus(locus)

    lines: List[str] = [
        f"LOCUS       {safe_locus:<16}{len(sequence):>9} aa    linear   UNK {_genbank_date()}",
        f"DEFINITION  {definition}.",
        f"ACCESSION   {safe_locus}",
        "KEYWORDS    .",
        "SOURCE      .",
        "  ORGANISM  .",
        "COMMENT     Generated by Fc Engineering Studio.",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(sequence)}",
        '                     /organism=["Homo sapiens"]',
    ]
    lines.extend(_format_origin_block(sequence))
    return "\n".join(lines) + "\n"


def format_csv_report(
    isotype: str, allotype: str, mutations: str, wt_seq: str, mut_seq: str
) -> str:
    """Build a per-position CSV report comparing WT and mutant sequences."""
    if not isinstance(wt_seq, str) or not wt_seq or not isinstance(mut_seq, str):
        return ""
    if max(len(wt_seq), len(mut_seq)) > MAX_EXPORT_SEQUENCE_LENGTH:
        return ""

    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["# isotype", isotype or "", allotype or "", mutations or ""])
    writer.writerow(["eu_position", "wt_aa", "mut_aa", "mutated"])
    for idx in range(max(len(wt_seq), len(mut_seq))):
        wt = wt_seq[idx] if idx < len(wt_seq) else ""
        mut = mut_seq[idx] if idx < len(mut_seq) else ""
        pos = get_eu_position(idx, isotype or "")
        writer.writerow(
            [pos if pos is not None else "", wt, mut, "yes" if wt != mut else "no"]
        )
    return buf.getvalue()


def format_project_fasta(
    project: AntibodyProject, header_style: str = "project_chain"
) -> Tuple[str, List[str]]:
    """Format expressible per-chain polypeptide sequences into multi-FASTA output."""
    if not isinstance(project, AntibodyProject):
        return "", ["Project object is required."]
    if not (1 <= len(project.chains) <= MAX_CHAINS):
        return "", [
            f"Project must contain between 1 and {MAX_CHAINS} polypeptide chains."
        ]

    preset_map = build_preset_map(project)
    safe_project = _sanitize_header_text(project.name, "Project")
    records: List[str] = []
    errors: List[str] = []

    for chain in project.chains:
        seq, _, chain_errs = assemble_chain_sequence(chain, preset_map)
        if chain_errs:
            errors.extend(chain_errs)
            continue
        safe_chain = _sanitize_header_text(chain.name, chain.id)
        if header_style == "chain_only":
            hdr = safe_chain
        elif header_style == "detailed":
            isos = ",".join(sorted({d.isotype for d in chain.domains if d.isotype}))
            muts = ",".join(
                p for d in chain.domains for p in (d.applied_presets + d.custom_mutations)
            )
            hdr = _sanitize_header_text(
                f"{safe_project}_{safe_chain}_{isos}_{muts or 'WT'}", safe_chain
            )
        else:
            hdr = f"{safe_project}_{safe_chain}"
        records.append(f">{hdr}\n{seq}")

    if errors:
        return "", errors
    return "\n".join(records) + "\n", []


def save_fasta_file(
    path: str, project: AntibodyProject, header_style: str = "project_chain"
) -> None:
    """Export project multi-FASTA sequences safely to disk."""
    resolved = _validate_file_path(path, (".fasta", ".fa"), auto_append=True)
    fasta_text, errors = format_project_fasta(project, header_style=header_style)
    if errors:
        raise ValueError("; ".join(errors))
    if len(fasta_text.encode("utf-8")) > config.MAX_FILE_SIZE:
        raise ValueError("FASTA output exceeds maximum allowed 1MB size.")
    try:
        with open(resolved, "w", encoding="utf-8") as f:
            f.write(fasta_text)
    except OSError:
        raise ValueError(f"Failed to write '{os.path.basename(resolved)}'.") from None


def _validate_ascii_id(val: Any, label: str) -> str:
    if not isinstance(val, str) or len(val) > 64 or not _VALID_ID_RE.fullmatch(val):
        raise ValueError(
            f"{label} must match [A-Za-z0-9_]+ and be at most 64 characters."
        )
    return val


def _validate_coord_num(val: Any, label: str) -> float:
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        raise ValueError(f"Coordinate {label} must be numeric.")
    fval = float(val)
    if not math.isfinite(fval) or fval < -MAX_COORD or fval > MAX_COORD:
        raise ValueError(f"Coordinate {label} out of allowed bounds.")
    return fval


def _validate_custom_mutation_defs(raw_defs: Any) -> List[Dict[str, str]]:
    if not isinstance(raw_defs, list):
        raise ValueError("file_custom_mutations must be a list.")
    if len(raw_defs) > MAX_FILE_CUSTOM_MUTATIONS:
        raise ValueError(
            f"Maximum of {MAX_FILE_CUSTOM_MUTATIONS} file custom mutations allowed."
        )
    builtin_labels = {
        item["label"]
        for item in config.COMMON_MUTATIONS
        if isinstance(item, dict) and "label" in item
    }
    seen_labels = set()
    validated: List[Dict[str, str]] = []
    for item in raw_defs:
        if not isinstance(item, dict):
            raise ValueError("Custom mutation entry must be a dictionary.")
        label = item.get("label")
        value = item.get("value")
        if not isinstance(label, str) or not isinstance(value, str):
            raise ValueError("Custom mutation label and value must be strings.")
        stripped = label.strip()
        if (
            not stripped
            or len(stripped) > 40
            or not all(0x20 <= ord(ch) <= 0x7E for ch in stripped)
        ):
            raise ValueError(f"Invalid custom mutation label '{label}'.")
        if stripped in builtin_labels:
            raise ValueError(
                f"Custom mutation label '{stripped}' shadows a built-in preset."
            )
        if stripped in seen_labels:
            raise ValueError(f"Duplicate custom mutation label '{stripped}'.")
        seen_labels.add(stripped)
        parsed, errs = parse_mutation_list(value)
        if errs or not parsed:
            raise ValueError(
                f"Invalid mutation value for '{stripped}': {'; '.join(errs) or 'empty'}"
            )
        validated.append(
            {
                "label": stripped,
                "value": "/".join(f"{wt}{pos}{mut}" for wt, pos, mut in parsed),
            }
        )
    return validated


def serialize_ab_project(project: AntibodyProject) -> str:
    """Serialize an AntibodyProject into a multi-LOCUS GenBank-style .ab string."""
    if not isinstance(project, AntibodyProject):
        raise ValueError("Expected AntibodyProject instance.")
    safe_proj_name = _validate_human_name(project.name, "Project name")
    if project.format_type not in VALID_FORMAT_TYPES:
        raise ValueError(f"Invalid format type '{project.format_type}'.")
    if not (1 <= len(project.chains) <= MAX_CHAINS):
        raise ValueError(f"Project must contain between 1 and {MAX_CHAINS} chains.")
    if len(project.pairings) > MAX_PAIRINGS:
        raise ValueError(f"Project exceeds maximum of {MAX_PAIRINGS} pairings.")
    custom_defs = _validate_custom_mutation_defs(project.file_custom_mutations)

    seen_chain_ids = set()
    seen_domain_ids = set()
    domain_endpoints = set()
    serialized_chains: List[Dict[str, Any]] = []
    assembled_records: List[Tuple[PolypeptideChain, str, List[Tuple[int, int, DomainBlock]]]] = []

    preset_map = build_preset_map(project)

    for chain in project.chains:
        cid = _validate_ascii_id(chain.id, "Chain id")
        if cid in seen_chain_ids:
            raise ValueError(f"Duplicate chain id '{cid}'.")
        seen_chain_ids.add(cid)
        cname = _validate_human_name(chain.name, f"Chain name ({cid})")
        cx = _validate_coord_num(chain.x, f"{cid}.x")
        cy = _validate_coord_num(chain.y, f"{cid}.y")
        if not (1 <= len(chain.domains) <= MAX_DOMAINS_PER_CHAIN):
            raise ValueError(
                f"Chain '{cid}' must contain between 1 and {MAX_DOMAINS_PER_CHAIN} domains."
            )

        serialized_domains: List[Dict[str, Any]] = []
        for dom in chain.domains:
            did = _validate_ascii_id(dom.id, "Domain id")
            if did in seen_domain_ids:
                raise ValueError(f"Duplicate domain id '{did}'.")
            seen_domain_ids.add(did)
            domain_endpoints.add((cid, did))
            dx = _validate_coord_num(dom.x, f"{did}.x")
            dy = _validate_coord_num(dom.y, f"{did}.y")
            if not isinstance(dom.pinned, bool):
                raise ValueError(f"Domain '{did}' pinned flag must be a boolean.")

            history_capped = [
                str(entry)[:MAX_METADATA_LENGTH]
                for entry in dom.mutation_history[-MAX_HISTORY_ENTRIES:]
            ]
            serialized_domains.append(
                {
                    "id": did,
                    "domain_type": dom.domain_type,
                    "isotype": dom.isotype,
                    "allotype": dom.allotype,
                    "template_key": dom.template_key,
                    "custom_sequence": dom.custom_sequence,
                    "x": dx,
                    "y": dy,
                    "pinned": dom.pinned,
                    "applied_presets": list(dom.applied_presets),
                    "custom_mutations": list(dom.custom_mutations),
                    "mutation_history": history_capped,
                }
            )

        chain_seq, spans, errs = assemble_chain_sequence(chain, preset_map)
        if errs:
            raise ValueError(f"Cannot serialize chain '{cid}': {'; '.join(errs)}")
        assembled_records.append((chain, chain_seq, spans))
        serialized_chains.append(
            {
                "id": cid,
                "name": cname,
                "x": cx,
                "y": cy,
                "domains": serialized_domains,
            }
        )

    seen_pairing_ids = set()
    serialized_pairings: List[Dict[str, str]] = []
    for p in project.pairings:
        pid = _validate_ascii_id(p.id, "Pairing id")
        if pid in seen_pairing_ids:
            raise ValueError(f"Duplicate pairing id '{pid}'.")
        seen_pairing_ids.add(pid)
        if p.pairing_type not in VALID_PAIRING_TYPES:
            raise ValueError(f"Invalid pairing type '{p.pairing_type}'.")
        ep_a = (p.chain_a_id, p.domain_a_id)
        ep_b = (p.chain_b_id, p.domain_b_id)
        if ep_a not in domain_endpoints or ep_b not in domain_endpoints or ep_a == ep_b:
            raise ValueError(f"Pairing '{pid}' references invalid or identical endpoints.")
        serialized_pairings.append(
            {
                "id": pid,
                "chain_a_id": p.chain_a_id,
                "domain_a_id": p.domain_a_id,
                "chain_b_id": p.chain_b_id,
                "domain_b_id": p.domain_b_id,
                "pairing_type": p.pairing_type,
            }
        )

    manifest = {
        "schema_version": 1,
        "name": safe_proj_name,
        "format_type": project.format_type,
        "file_custom_mutations": custom_defs,
        "chains": serialized_chains,
        "pairings": serialized_pairings,
    }
    manifest_lines = json.dumps(manifest, ensure_ascii=True, indent=2).splitlines()

    genbank_date = _genbank_date()
    output_lines: List[str] = []

    for idx, (chain, chain_seq, spans) in enumerate(assembled_records):
        safe_locus = sanitize_locus(chain.id)
        def_str = _sanitize_header_text(
            f"{safe_proj_name} - {chain.name}", safe_locus
        )
        output_lines.append(
            f"LOCUS       {safe_locus:<16}{len(chain_seq):>9} aa    linear   UNK {genbank_date}"
        )
        output_lines.append(f"DEFINITION  {def_str}.")
        output_lines.append(f"ACCESSION   {safe_locus}")
        output_lines.append("KEYWORDS    .")
        output_lines.append("SOURCE      .")
        output_lines.append("  ORGANISM  .")
        if idx == 0:
            output_lines.append("COMMENT     Generated by Fc Engineering Studio.")
            output_lines.append(f"            {_MANIFEST_BEGIN}")
            for line in manifest_lines:
                output_lines.append(f"            {line}")
            output_lines.append(f"            {_MANIFEST_END}")
        else:
            output_lines.append("COMMENT     Generated by Fc Engineering Studio.")
        output_lines.append("FEATURES             Location/Qualifiers")
        output_lines.append(f"     source          1..{len(chain_seq)}")
        output_lines.append('                     /organism=["Homo sapiens"]')

        for start, end, dom in spans:
            output_lines.append(f"     region          {start + 1}..{end}")
            output_lines.append(
                f'                     /domain_type="{_sanitize_qualifier(dom.domain_type)}"'
            )
            detail = dom.template_key or dom.allotype or dom.isotype
            output_lines.append(
                f'                     /label="{_sanitize_qualifier(f"{dom.domain_type}:{detail}")}"'
            )
            output_lines.append(
                f'                     /coords="{dom.x:.1f},{dom.y:.1f}"'
            )
            eff_muts, _ = get_effective_domain_mutations(dom, preset_map)
            if eff_muts:
                output_lines.append(
                    f'                     /mutations="{_sanitize_qualifier("/".join(eff_muts))}"'
                )
            if eff_muts and dom.domain_type in CONSTANT_DOMAIN_EU_RANGES:
                wt_dom, _ = resolve_domain_base_sequence(dom)
                mut_dom, _ = resolve_domain_sequence(dom, preset_map)
                for eu_pos, wt_aa, mut_aa in diff_sequences(
                    wt_dom, mut_dom, dom.isotype, domain_type=dom.domain_type
                ):
                    local_idx = get_domain_residue_index(
                        eu_pos, dom.domain_type, dom.isotype
                    )
                    if local_idx is not None:
                        chain_pos_1 = start + local_idx + 1
                        output_lines.append(f"     variation       {chain_pos_1}")
                        note = _sanitize_qualifier(f"EU_{eu_pos}:{wt_aa}->{mut_aa}")
                        output_lines.append(f'                     /note="{note}"')

        output_lines.extend(_format_origin_block(chain_seq))

    payload = "\n".join(output_lines) + "\n"
    if len(payload.encode("utf-8")) > config.MAX_FILE_SIZE:
        raise ValueError("Serialized .ab project exceeds maximum allowed 1MB size.")
    return payload


def _split_genbank_records(content: str) -> List[List[str]]:
    lines = content.splitlines()
    records: List[List[str]] = []
    current: List[str] = []
    saw_first_record = False

    for line in lines:
        if not current and not line.strip():
            continue
        current.append(line)
        if line.strip() == "//":
            records.append(current)
            current = []
            saw_first_record = True

    if current:
        if any(item.strip() for item in current):
            if saw_first_record:
                raise ValueError("Unexpected trailing content after final '//' delimiter.")
            raise ValueError("Incomplete GenBank record missing '//' delimiter.")

    if not records:
        raise ValueError("No GenBank LOCUS records found in .ab payload.")
    return records


def _extract_manifest_json(first_record: List[str]) -> str:
    begin_indices = [
        i
        for i, line in enumerate(first_record)
        if line.strip() == _MANIFEST_BEGIN and i > 0
    ]
    end_indices = [
        i
        for i, line in enumerate(first_record)
        if line.strip() == _MANIFEST_END and i > 0
    ]
    if len(begin_indices) != 1 or len(end_indices) != 1:
        raise ValueError("Missing or ambiguous project manifest markers in first record.")
    start_idx, end_idx = begin_indices[0], end_indices[0]
    if end_idx <= start_idx:
        raise ValueError("Invalid project manifest marker ordering.")

    # Verify markers reside between COMMENT and FEATURES/ORIGIN sections
    has_comment_before = any(
        line.startswith("COMMENT") for line in first_record[:start_idx]
    )
    has_features_before = any(
        line.startswith("FEATURES") for line in first_record[:start_idx]
    )
    if not has_comment_before or has_features_before:
        raise ValueError("Project manifest must reside inside COMMENT block.")

    json_lines: List[str] = []
    for raw in first_record[start_idx + 1 : end_idx]:
        json_lines.append(raw[12:] if raw.startswith("            ") else raw.lstrip())
    return "\n".join(json_lines)


def _extract_origin_sequence(record_lines: List[str]) -> str:
    origin_idx: Optional[int] = None
    for idx, line in enumerate(record_lines):
        if line.strip() == "ORIGIN":
            origin_idx = idx
            break
    if origin_idx is None:
        raise ValueError("LOCUS record missing ORIGIN section.")
    body_lines = record_lines[origin_idx + 1 : -1]  # exclude trailing //
    letters = re.sub(r"[^A-Za-z]", "", "".join(body_lines)).upper()
    if not letters or not _VALID_SEQ_RE.fullmatch(letters):
        raise ValueError("Invalid or empty ORIGIN amino acid sequence.")
    return letters


def _validate_str_list(
    raw: Any, field_name: str, max_len: int, is_mut_token: bool = False
) -> List[str]:
    if not isinstance(raw, list):
        raise ValueError(f"{field_name} must be a list.")
    if len(raw) > max_len:
        raise ValueError(f"{field_name} exceeds maximum of {max_len} entries.")
    result: List[str] = []
    for item in raw:
        if not isinstance(item, str):
            raise ValueError(f"{field_name} items must be strings.")
        if is_mut_token:
            if len(item) > 10 or not _MUTATION_TOKEN_RE.fullmatch(item):
                raise ValueError(f"Invalid mutation token '{item}' in {field_name}.")
        else:
            if len(item) > MAX_METADATA_LENGTH or not all(
                0x20 <= ord(ch) <= 0x7E for ch in item
            ):
                raise ValueError(f"Invalid metadata string in {field_name}.")
        result.append(item)
    return result


def _reject_json_constant(const_str: str) -> Any:
    raise ValueError(f"Invalid JSON numeric constant '{const_str}'.")


def parse_ab_project(content: str) -> AntibodyProject:
    """Parse and strictly validate a GenBank-style .ab string back into an AntibodyProject."""
    if not isinstance(content, str):
        raise ValueError("Project content must be a string.")
    if len(content.encode("utf-8")) > config.MAX_FILE_SIZE:
        raise ValueError("Project content exceeds maximum 1MB size.")

    records = _split_genbank_records(content)
    raw_manifest = _extract_manifest_json(records[0])

    try:
        data = json.loads(raw_manifest, parse_constant=_reject_json_constant)
    except (json.JSONDecodeError, RecursionError):
        raise ValueError("Malformed or deeply nested JSON manifest in .ab file.") from None

    if not isinstance(data, dict):
        raise ValueError("Manifest payload must be a JSON object.")

    sv = data.get("schema_version")
    if isinstance(sv, bool) or not isinstance(sv, int) or sv != 1:
        raise ValueError(f"Unsupported .ab schema version: {sv!r}.")

    proj_name = _validate_human_name(data.get("name"), "Project name")
    format_type = data.get("format_type")
    if format_type not in VALID_FORMAT_TYPES:
        raise ValueError(f"Unsupported format_type '{format_type}'.")

    custom_defs = _validate_custom_mutation_defs(data.get("file_custom_mutations", []))

    raw_chains = data.get("chains")
    if not isinstance(raw_chains, list) or not (1 <= len(raw_chains) <= MAX_CHAINS):
        raise ValueError(f"Project must contain between 1 and {MAX_CHAINS} chains.")

    seen_chain_ids = set()
    seen_domain_ids = set()
    domain_endpoints = set()
    chains: List[PolypeptideChain] = []

    for c_item in raw_chains:
        if not isinstance(c_item, dict):
            raise ValueError("Chain entry must be a JSON object.")
        cid = _validate_ascii_id(c_item.get("id"), "Chain id")
        if cid in seen_chain_ids:
            raise ValueError(f"Duplicate chain id '{cid}'.")
        seen_chain_ids.add(cid)
        cname = _validate_human_name(c_item.get("name"), f"Chain name ({cid})")
        cx = _validate_coord_num(c_item.get("x"), f"{cid}.x")
        cy = _validate_coord_num(c_item.get("y"), f"{cid}.y")

        raw_domains = c_item.get("domains")
        if not isinstance(raw_domains, list) or not (
            1 <= len(raw_domains) <= MAX_DOMAINS_PER_CHAIN
        ):
            raise ValueError(
                f"Chain '{cid}' must contain between 1 and {MAX_DOMAINS_PER_CHAIN} domains."
            )

        domains: List[DomainBlock] = []
        for d_item in raw_domains:
            if not isinstance(d_item, dict):
                raise ValueError("Domain entry must be a JSON object.")
            did = _validate_ascii_id(d_item.get("id"), "Domain id")
            if did in seen_domain_ids:
                raise ValueError(f"Duplicate domain id '{did}'.")
            seen_domain_ids.add(did)
            domain_endpoints.add((cid, did))

            dtype = d_item.get("domain_type")
            if dtype not in VALID_DOMAIN_TYPES:
                raise ValueError(f"Invalid domain_type '{dtype}'.")

            isotype = d_item.get("isotype", "igg1")
            if not isinstance(isotype, str) or isotype not in config.SEQUENCES:
                raise ValueError(f"Invalid isotype '{isotype}' on domain '{did}'.")

            allotype = d_item.get("allotype", "")
            if not isinstance(allotype, str) or (
                allotype != "" and allotype not in config.SEQUENCES[isotype]
            ):
                raise ValueError(
                    f"Invalid allotype '{allotype}' for isotype '{isotype}'."
                )

            tkey = d_item.get("template_key", "")
            if not isinstance(tkey, str):
                raise ValueError(f"template_key must be a string on domain '{did}'.")
            if dtype in ALLOWED_TEMPLATE_KEYS_BY_TYPE and tkey != "":
                if tkey not in ALLOWED_TEMPLATE_KEYS_BY_TYPE[dtype]:
                    raise ValueError(
                        f"Invalid template_key '{tkey}' for domain_type '{dtype}'."
                    )

            cseq = d_item.get("custom_sequence", "")
            if not isinstance(cseq, str) or len(cseq) > MAX_CUSTOM_SEQ_LENGTH or not _VALID_OPT_SEQ_RE.fullmatch(cseq):
                raise ValueError(f"Invalid custom_sequence on domain '{did}'.")
            if dtype in CONSTANT_DOMAIN_EU_RANGES or dtype == "CL":
                if cseq:
                    raise ValueError(
                        f"custom_sequence is forbidden on constant domain '{dtype}'."
                    )

            dx = _validate_coord_num(d_item.get("x"), f"{did}.x")
            dy = _validate_coord_num(d_item.get("y"), f"{did}.y")
            pinned = d_item.get("pinned", False)
            if not isinstance(pinned, bool):
                raise ValueError(f"pinned must be boolean on domain '{did}'.")

            applied_presets = _validate_str_list(
                d_item.get("applied_presets", []), "applied_presets", 50
            )
            custom_mutations = _validate_str_list(
                d_item.get("custom_mutations", []),
                "custom_mutations",
                50,
                is_mut_token=True,
            )
            mutation_history = _validate_str_list(
                d_item.get("mutation_history", []), "mutation_history", 50
            )

            domains.append(
                DomainBlock(
                    id=did,
                    domain_type=dtype,
                    isotype=isotype,
                    allotype=allotype,
                    template_key=tkey,
                    custom_sequence=cseq,
                    x=dx,
                    y=dy,
                    pinned=pinned,
                    applied_presets=applied_presets,
                    custom_mutations=custom_mutations,
                    mutation_history=mutation_history,
                )
            )

        chains.append(PolypeptideChain(id=cid, name=cname, domains=domains, x=cx, y=cy))

    raw_pairings = data.get("pairings", [])
    if not isinstance(raw_pairings, list) or len(raw_pairings) > MAX_PAIRINGS:
        raise ValueError(f"pairings exceeds maximum of {MAX_PAIRINGS} entries.")

    pairings: List[InterChainPairing] = []
    seen_pairing_ids = set()
    for p_item in raw_pairings:
        if not isinstance(p_item, dict):
            raise ValueError("Pairing entry must be a JSON object.")
        pid = _validate_ascii_id(p_item.get("id"), "Pairing id")
        if pid in seen_pairing_ids:
            raise ValueError(f"Duplicate pairing id '{pid}'.")
        seen_pairing_ids.add(pid)
        ptype = p_item.get("pairing_type")
        if ptype not in VALID_PAIRING_TYPES:
            raise ValueError(f"Invalid pairing_type '{ptype}'.")
        ca = p_item.get("chain_a_id")
        da = p_item.get("domain_a_id")
        cb = p_item.get("chain_b_id")
        db = p_item.get("domain_b_id")
        if (ca, da) not in domain_endpoints or (cb, db) not in domain_endpoints:
            raise ValueError(f"Pairing '{pid}' references missing chain/domain.")
        if (ca, da) == (cb, db):
            raise ValueError(f"Pairing '{pid}' cannot connect a domain to itself.")
        pairings.append(
            InterChainPairing(
                id=pid,
                chain_a_id=ca,
                domain_a_id=da,
                chain_b_id=cb,
                domain_b_id=db,
                pairing_type=ptype,
            )
        )

    project = AntibodyProject(
        name=proj_name,
        format_type=format_type,
        chains=chains,
        pairings=pairings,
        file_custom_mutations=custom_defs,
    )

    if len(records) != len(project.chains):
        raise ValueError(
            f"Record count mismatch: found {len(records)} LOCUS records for {len(project.chains)} chains."
        )

    preset_map = build_preset_map(project)
    for idx, (chain, rec_lines) in enumerate(zip(project.chains, records)):
        origin_seq = _extract_origin_sequence(rec_lines)
        if len(origin_seq) > MAX_EXPORT_SEQUENCE_LENGTH:
            raise ValueError(f"ORIGIN sequence in record {idx + 1} exceeds 10000 aa.")
        assembled_seq, _, errs = assemble_chain_sequence(chain, preset_map)
        if errs:
            raise ValueError(
                f"Domain resolution failed for chain '{chain.id}': {'; '.join(errs)}"
            )
        if assembled_seq != origin_seq:
            raise ValueError(
                f"ORIGIN sequence mismatch for chain '{chain.id}' in record {idx + 1}."
            )

    return project


def save_ab_file(path: str, project: AntibodyProject) -> None:
    """Serialize and save an AntibodyProject as a GenBank-style .ab file."""
    resolved = _validate_file_path(path, (".ab",), auto_append=True)
    payload = serialize_ab_project(project)
    if len(payload.encode("utf-8")) > config.MAX_FILE_SIZE:
        raise ValueError("Serialized project exceeds 1MB file limit.")
    try:
        with open(resolved, "w", encoding="utf-8") as f:
            f.write(payload)
    except OSError:
        raise ValueError(f"Failed to write '{os.path.basename(resolved)}'.") from None


def save_genbank_file(path: str, project: AntibodyProject) -> None:
    """Export all assembled chains as a multi-record GenBank (.gb) file."""
    resolved = _validate_file_path(path, (".gb", ".gbk"), auto_append=True)
    preset_map = build_preset_map(project)
    records: List[str] = []
    for chain in project.chains:
        seq, _, errs = assemble_chain_sequence(chain, preset_map)
        if errs:
            raise ValueError("; ".join(errs))
        rec = format_genbank(
            f"{project.name}_{chain.id}",
            seq,
            definition=f"{project.name} - {chain.name}",
        )
        if not rec:
            raise ValueError(f"Failed to format GenBank record for chain '{chain.name}'.")
        records.append(rec)
    payload = "".join(records)
    if len(payload.encode("utf-8")) > config.MAX_FILE_SIZE:
        raise ValueError("Serialized GenBank file exceeds 1MB limit.")
    try:
        with open(resolved, "w", encoding="utf-8") as f:
            f.write(payload)
    except OSError:
        raise ValueError(f"Failed to write '{os.path.basename(resolved)}'.") from None


def save_csv_report_file(path: str, project: AntibodyProject) -> None:
    """Export constant-domain EU diffs across all chains to a .csv report."""
    resolved = _validate_file_path(path, (".csv",), auto_append=True)
    preset_map = build_preset_map(project)
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["chain_id", "chain_name", "domain_type", "isotype", "eu_position", "wt_aa", "mut_aa", "mutated"])
    for chain in project.chains:
        for dom in chain.domains:
            if dom.domain_type not in CONSTANT_DOMAIN_EU_RANGES:
                continue
            wt_seq, base_errs = resolve_domain_base_sequence(dom)
            mut_seq, mut_errs = resolve_domain_sequence(dom, preset_map)
            if base_errs or mut_errs:
                raise ValueError("; ".join(base_errs + mut_errs))
            for idx in range(max(len(wt_seq), len(mut_seq))):
                wt = wt_seq[idx] if idx < len(wt_seq) else ""
                mut = mut_seq[idx] if idx < len(mut_seq) else ""
                eu_pos = get_domain_eu_position(idx, dom.domain_type, dom.isotype)
                writer.writerow(
                    [
                        chain.id,
                        chain.name,
                        dom.domain_type,
                        dom.isotype,
                        eu_pos if eu_pos is not None else "",
                        wt,
                        mut,
                        "yes" if wt != mut else "no",
                    ]
                )
    payload = buf.getvalue()
    if len(payload.encode("utf-8")) > config.MAX_FILE_SIZE:
        raise ValueError("Serialized CSV report exceeds 1MB limit.")
    try:
        with open(resolved, "w", encoding="utf-8") as f:
            f.write(payload)
    except OSError:
        raise ValueError(f"Failed to write '{os.path.basename(resolved)}'.") from None


def load_ab_file(path: str) -> AntibodyProject:
    """Read and validate an .ab project file safely within 1MB size bounds."""
    resolved = _validate_file_path(path, (".ab",), auto_append=False)
    base_name = os.path.basename(resolved)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"File '{base_name}' does not exist.")
    if not os.path.isfile(resolved):
        raise ValueError(f"'{base_name}' is not a regular file.")
    try:
        if os.path.getsize(resolved) > config.MAX_FILE_SIZE:
            raise ValueError(f"File '{base_name}' exceeds maximum 1MB size.")
        with open(resolved, "rb") as f:
            raw_bytes = f.read(config.MAX_FILE_SIZE + 1)
        if len(raw_bytes) > config.MAX_FILE_SIZE:
            raise ValueError(f"File '{base_name}' exceeds maximum 1MB size.")
        content = raw_bytes.decode("utf-8")
    except (OSError, UnicodeDecodeError):
        raise ValueError(f"Unable to read '{base_name}' as UTF-8.") from None
    return parse_ab_project(content)


def _draw_dashed_line(
    draw: ImageDraw.ImageDraw,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    fill: Tuple[int, int, int, int],
    width: int = 2,
    dash_len: float = 8.0,
    gap_len: float = 5.0,
) -> None:
    dx = x2 - x1
    dy = y2 - y1
    dist = math.hypot(dx, dy)
    if dist < 1e-6:
        return
    ux, uy = dx / dist, dy / dist
    pos = 0.0
    while pos < dist:
        seg_end = min(dist, pos + dash_len)
        draw.line(
            [
                (x1 + ux * pos, y1 + uy * pos),
                (x1 + ux * seg_end, y1 + uy * seg_end),
            ],
            fill=fill,
            width=width,
        )
        pos += dash_len + gap_len


def _draw_dashed_arc(
    draw: ImageDraw.ImageDraw,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    fill: Tuple[int, int, int, int],
    width: int = 2,
    lift: float = 34.0,
    orientation: str = "vertical",
) -> None:
    if orientation == "vertical":
        ctrl_x = max(x1, x2) + lift
        ctrl_y = (y1 + y2) / 2.0
    else:
        ctrl_x = (x1 + x2) / 2.0
        ctrl_y = min(y1, y2) - lift
    steps = 18
    pts: List[Tuple[float, float]] = []
    for i in range(steps + 1):
        t = i / steps
        inv = 1.0 - t
        px = inv * inv * x1 + 2.0 * inv * t * ctrl_x + t * t * x2
        py = inv * inv * y1 + 2.0 * inv * t * ctrl_y + t * t * y2
        pts.append((px, py))
    for i in range(0, len(pts) - 1, 2):
        draw.line([pts[i], pts[i + 1]], fill=fill, width=width)


def _load_scaled_font(scale: float) -> ImageFont.ImageFont | ImageFont.FreeTypeFont:
    font_px = max(10, int(round(10.0 * scale)))
    for candidate in (
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "DejaVuSans.ttf",
        "Arial.ttf",
    ):
        try:
            return ImageFont.truetype(candidate, size=font_px)
        except OSError:
            continue
    return ImageFont.load_default()


def export_schematic_png(
    project: AntibodyProject,
    scale: float = 2.0,
    theme: Optional[SchematicTheme] = None,
    bg_color: Optional[str] = None,
) -> bytes:
    """Render a transparent-background (or solid bg_color) RGBA .png 2D schematic."""
    if (
        isinstance(scale, bool)
        or not isinstance(scale, (int, float))
        or not math.isfinite(scale)
        or scale < 0.5
        or scale > 4.0
    ):
        raise ValueError("Scale must be a finite number between 0.5 and 4.0.")

    active_theme = theme if theme is not None else DEFAULT_SCHEMATIC_THEME
    scene = build_schematic_scene(project, build_preset_map(project), theme=active_theme)
    min_x, min_y, max_x, max_y = scene.bounds
    width_px = max(1, int(math.ceil((max_x - min_x) * scale)))
    height_px = max(1, int(math.ceil((max_y - min_y) * scale)))
    if width_px > MAX_PNG_DIMENSION or height_px > MAX_PNG_DIMENSION:
        raise ValueError(
            f"Rendered PNG dimensions ({width_px}x{height_px}) exceed MAX_PNG_DIMENSION={MAX_PNG_DIMENSION}."
        )

    img = Image.new("RGBA", (width_px, height_px), bg_color if bg_color else (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    font = _load_scaled_font(scale)

    def tx(x: float) -> float:
        return (x - min_x) * scale

    def ty(y: float) -> float:
        return (y - min_y) * scale

    stroke_w = max(1, int(round(1.5 * scale)))

    # 1. Draw inter-chain and intra-chain pairings behind domain blocks
    disulfide_midpoints: List[Tuple[float, float]] = []
    for pedge in scene.pairing_edges:
        x1, y1, x2, y2 = tx(pedge.x1), ty(pedge.y1), tx(pedge.x2), ty(pedge.y2)
        if pedge.chain_a_id == pedge.chain_b_id and pedge.pairing_type == "VH_VL":
            _draw_dashed_arc(
                draw,
                x1,
                y1,
                x2,
                y2,
                fill=(217, 119, 87, 230),
                width=stroke_w,
                lift=34.0 * scale,
                orientation=active_theme.orientation,
            )
        elif pedge.pairing_type == "HINGE_DISULFIDE":
            dx = x2 - x1
            dy = y2 - y1
            dist = math.hypot(dx, dy) or 1.0
            nx, ny = (-dy / dist) * (3.0 * scale), (dx / dist) * (3.0 * scale)
            for sign in (-1.0, 1.0):
                draw.line(
                    [(x1 + nx * sign, y1 + ny * sign), (x2 + nx * sign, y2 + ny * sign)],
                    fill=(245, 158, 11, 245),
                    width=stroke_w,
                )
            disulfide_midpoints.append(((x1 + x2) / 2.0, (y1 + y2) / 2.0))
        elif pedge.is_disulfide:
            draw.line(
                [(x1, y1), (x2, y2)],
                fill=(245, 158, 11, 245),
                width=max(2, int(round(2.0 * scale))),
            )
            disulfide_midpoints.append(((x1 + x2) / 2.0, (y1 + y2) / 2.0))
        else:
            _draw_dashed_line(
                draw,
                x1,
                y1,
                x2,
                y2,
                fill=(158, 152, 140, 230),
                width=stroke_w,
                dash_len=6.0 * scale,
                gap_len=4.0 * scale,
            )

    # 2. Draw solid N-to-C backbone connectors
    for edge in scene.edges:
        draw.line(
            [(tx(edge.x1), ty(edge.y1)), (tx(edge.x2), ty(edge.y2))],
            fill=(125, 120, 109, 255),
            width=max(2, int(round(2.0 * scale))),
        )

    # 3. Draw domain shapes, chain names, domain labels, and mutation badges
    for node in scene.nodes:
        x0 = tx(node.x)
        y0 = ty(node.y)
        w_px = node.width * scale
        h_px = node.height * scale
        x1 = x0 + w_px
        y1 = y0 + h_px
        outline_col = "#EF4444" if node.error_flag else node.border_color
        outline_w = max(1, int(round((max(3.0, node.border_width + 1.5) if node.error_flag else node.border_width) * scale)))

        if node.shape == "ellipse":
            draw.ellipse([x0, y0, x1, y1], fill=node.fill_color, outline=outline_col, width=outline_w)
        elif node.shape == "hexagon":
            verts = compute_hexagon_vertices(x0, y0, w_px, h_px)
            draw.polygon(verts, fill=node.fill_color, outline=outline_col, width=outline_w)
        else:
            if node.shape == "rect":
                radius = 0.0
            elif node.shape == "capsule":
                radius = min(w_px, h_px) / 2.0
            else:
                radius = 8.0 * scale
            draw.rounded_rectangle(
                [x0, y0, x1, y1],
                radius=radius,
                fill=node.fill_color,
                outline=outline_col,
                width=outline_w,
            )

        cx = (x0 + x1) / 2.0
        cy = (y0 + y1) / 2.0
        draw.text((cx, cy), node.label, fill=node.text_color, font=font, anchor="mm")

        if node.is_n_terminus and node.chain_name:
            draw.text(
                (cx, max(2.0, y0 - 11.0 * scale)),
                node.chain_name,
                fill=(100, 116, 139, 255),
                font=font,
                anchor="mm",
            )
        if node.badge_text:
            draw.text(
                (cx, y1 + 8.0 * scale),
                node.badge_text,
                fill=(220, 38, 38, 255),
                font=font,
                anchor="mm",
            )

    # 4. Draw -S-S- cysteine disulfide bond badges on top of bridges
    for mx, my in disulfide_midpoints:
        pw, ph = 18.0 * scale, 8.0 * scale
        draw.rounded_rectangle(
            [mx - pw, my - ph, mx + pw, my + ph],
            radius=4.0 * scale,
            fill=(20, 20, 19, 235),
            outline=(245, 158, 11, 255),
            width=max(1, int(round(1.0 * scale))),
        )
        draw.text((mx, my), "-S-S-", fill=(253, 230, 138, 255), font=font, anchor="mm")

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def save_schematic_png(
    path: str,
    project: AntibodyProject,
    scale: float = 2.0,
    theme: Optional[SchematicTheme] = None,
    bg_color: Optional[str] = None,
) -> None:
    """Save a transparent-background (or solid bg_color) RGBA .png schematic to disk safely."""
    resolved = _validate_file_path(path, (".png",), auto_append=True)
    png_bytes = export_schematic_png(project, scale=scale, theme=theme, bg_color=bg_color)
    if len(png_bytes) > config.MAX_FILE_SIZE:
        raise ValueError("Rendered PNG exceeds maximum allowed 1MB file size.")
    try:
        with open(resolved, "wb") as f:
            f.write(png_bytes)
    except OSError:
        raise ValueError(f"Failed to write '{os.path.basename(resolved)}'.") from None
