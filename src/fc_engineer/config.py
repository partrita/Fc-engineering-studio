import logging
import os
import re
import sys
import traceback
from typing import Any, Dict, List, Optional, Tuple

import yaml

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# SECURITY: Max file size of 1MB to prevent DoS via memory exhaustion
MAX_FILE_SIZE = 1 * 1024 * 1024

_VALID_SEQ_RE = re.compile(r"[A-Z]+\Z")


def log_sanitized_error(logger_obj: logging.Logger, msg: str, e: Exception) -> None:
    """Log exception details with package base directory paths redacted."""
    base_path = os.path.dirname(__file__)
    safe_err = str(e).replace(base_path, ".") if base_path else str(e)
    safe_traceback = (
        traceback.format_exc().replace(base_path, ".")
        if base_path
        else traceback.format_exc()
    )
    logger_obj.error(f"{msg}: {safe_err}\n{safe_traceback}")


# --- Configuration & Data Loading ---


class NoAliasSafeLoader(yaml.SafeLoader):
    """Safe YAML loader blocking alias expansion and deep nesting bombs."""

    MAX_DEPTH = 50

    def compose_node(self, parent: Any, index: Any) -> Any:
        if self.check_event(yaml.events.AliasEvent):
            raise yaml.constructor.ConstructorError(
                "Aliases are not allowed to prevent YAML bomb (DoS) attacks."
            )

        if not hasattr(self, "depth"):
            self.depth = 0
        self.depth += 1

        if self.depth > self.MAX_DEPTH:
            raise yaml.constructor.ConstructorError(
                "Maximum YAML nesting depth exceeded to prevent DoS attacks."
            )

        node = super().compose_node(parent, index)
        self.depth -= 1
        return node


def _data_dir() -> str:
    local = os.path.join(os.path.dirname(__file__), "data")
    if os.path.isdir(local):
        return local
    if getattr(sys, "frozen", False):
        meipass = getattr(sys, "_MEIPASS", "")
        if meipass:
            frozen_data = os.path.join(meipass, "fc_engineer", "data")
            if os.path.isdir(frozen_data):
                return frozen_data
    return local


def _load_yaml_dict(path: str) -> Dict[str, Any]:
    """Load a YAML file safely with regular-file, 1MB size, and alias guards."""
    file_name = os.path.basename(path)
    try:
        if not os.path.isfile(path):
            print(f"Error: Missing configuration file {file_name}.", file=sys.stderr)
            return {}
        if os.path.getsize(path) > MAX_FILE_SIZE:
            print(f"Error: {file_name} exceeds 1MB limit.", file=sys.stderr)
            raise ValueError(f"File {file_name} exceeds maximum size of 1MB")
        with open(path, "r", encoding="utf-8") as f:
            content = f.read(MAX_FILE_SIZE + 1)
            if len(content.encode("utf-8")) > MAX_FILE_SIZE:
                print(f"Error: {file_name} content exceeds 1MB limit.", file=sys.stderr)
                raise ValueError(
                    f"File {file_name} content exceeds maximum size of 1MB"
                )
            data = yaml.load(content, Loader=NoAliasSafeLoader)  # nosec B506
            if not isinstance(data, dict):
                print(
                    f"Error: Parsed data from {file_name} is not a dictionary.",
                    file=sys.stderr,
                )
                return {}
            return data
    except Exception as e:
        log_sanitized_error(logger, f"Error loading {file_name}", e)
        print(
            f"Error loading {file_name}: An unexpected error occurred.",
            file=sys.stderr,
        )
        return {}


def _is_valid_seq_map(val: Any) -> bool:
    return isinstance(val, dict) and all(
        isinstance(k, str)
        and bool(k)
        and isinstance(v, str)
        and bool(_VALID_SEQ_RE.fullmatch(v))
        for k, v in val.items()
    )


def _load_sequences_bundle() -> (
    Tuple[Dict[str, Dict[str, str]], Dict[str, str], Dict[str, str], Dict[str, str]]
):
    """Read sequences.yaml once and return validated constant and template tables."""
    seq_path = os.path.join(_data_dir(), "sequences.yaml")
    data = _load_yaml_dict(seq_path)
    if not data:
        return {}, {}, {}, {}

    raw_iso = data.get("isotypes")
    isotypes: Dict[str, Dict[str, str]] = (
        raw_iso
        if (
            isinstance(raw_iso, dict)
            and all(
                isinstance(iso_k, str) and bool(iso_k) and _is_valid_seq_map(iso_v)
                for iso_k, iso_v in raw_iso.items()
            )
        )
        else {}
    )

    raw_var = data.get("variable_domains")
    var_domains: Dict[str, str] = raw_var if _is_valid_seq_map(raw_var) else {}

    raw_light = data.get("light_constant_domains")
    light_domains: Dict[str, str] = raw_light if _is_valid_seq_map(raw_light) else {}

    raw_linkers = data.get("linkers")
    linkers: Dict[str, str] = raw_linkers if _is_valid_seq_map(raw_linkers) else {}

    return isotypes, var_domains, light_domains, linkers


def _load_common_mutations() -> List[Dict[str, str]]:
    """Load and schema-validate common Fc mutation presets from mutants.yaml."""
    mut_path = os.path.join(_data_dir(), "mutants.yaml")
    data = _load_yaml_dict(mut_path)
    if not data:
        return []
    val = data.get("common_mutations")
    if isinstance(val, list) and all(
        isinstance(item, dict)
        and isinstance(item.get("value"), str)
        and isinstance(item.get("label"), str)
        for item in val
    ):
        return val
    return []


def load_yaml_data() -> Tuple[Dict[str, Dict[str, str]], List[Dict[str, str]]]:
    """Load isotype sequences and common mutation presets without side effects."""
    isotypes, _, _, _ = _load_sequences_bundle()
    common_muts = _load_common_mutations()
    return isotypes, common_muts


def load_domain_templates() -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Load built-in variable domain, light constant, and linker templates."""
    _, var_domains, light_domains, linkers = _load_sequences_bundle()
    return var_domains, light_domains, linkers


SEQUENCES: Dict[str, Dict[str, str]] = {}
COMMON_MUTATIONS: List[Dict[str, str]] = []
VARIABLE_DOMAINS: Dict[str, str] = {}
LIGHT_CONSTANT_DOMAINS: Dict[str, str] = {}
DEFAULT_LINKERS: Dict[str, str] = {}


def _load_all_config() -> None:
    """Populate module configuration constants from a single pass over data YAML files."""
    global SEQUENCES, COMMON_MUTATIONS, VARIABLE_DOMAINS, LIGHT_CONSTANT_DOMAINS, DEFAULT_LINKERS
    SEQUENCES, VARIABLE_DOMAINS, LIGHT_CONSTANT_DOMAINS, DEFAULT_LINKERS = (
        _load_sequences_bundle()
    )
    COMMON_MUTATIONS = _load_common_mutations()


_load_all_config()


_HEX_COLOR_RE = re.compile(r"#[0-9A-Fa-f]{6}\Z")


def _validate_hex_color(val: Any, field_name: str) -> str:
    if not isinstance(val, str) or not _HEX_COLOR_RE.fullmatch(val):
        raise ValueError(f"{field_name} must be a #RRGGBB hex color string.")
    return val.upper()


def _validate_bounded_float(
    val: Any, field_name: str, low: float, high: float
) -> float:
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        raise ValueError(f"{field_name} must be a number.")
    fval = float(val)
    if fval < low or fval > high:
        raise ValueError(f"{field_name} must be between {low} and {high}.")
    return fval


def _parse_theme_dict(name: str, raw: Any) -> Any:
    from fc_engineer.designer import (
        DEFAULT_SCHEMATIC_THEME,
        DomainVisualStyle,
        SchematicTheme,
        VALID_DOMAIN_TYPES,
        VALID_ORIENTATIONS,
        VALID_SHAPES,
    )

    if not isinstance(raw, dict):
        raise ValueError("Theme entry must be a mapping.")
    orientation = raw.get("orientation", DEFAULT_SCHEMATIC_THEME.orientation)
    if orientation not in VALID_ORIENTATIONS:
        raise ValueError(f"Invalid theme orientation '{orientation}'.")
    domain_spacing = _validate_bounded_float(
        raw.get("domain_spacing", DEFAULT_SCHEMATIC_THEME.domain_spacing),
        "domain_spacing",
        4.0,
        80.0,
    )
    chain_spacing = _validate_bounded_float(
        raw.get("chain_spacing", DEFAULT_SCHEMATIC_THEME.chain_spacing),
        "chain_spacing",
        60.0,
        260.0,
    )
    raw_styles = raw.get("styles")
    if not isinstance(raw_styles, dict):
        raise ValueError("Theme 'styles' must be a mapping.")

    styles: Dict[str, DomainVisualStyle] = {}
    for dtype in VALID_DOMAIN_TYPES:
        dflt = DEFAULT_SCHEMATIC_THEME.styles[dtype]
        entry = raw_styles.get(dtype)
        if entry is None:
            styles[dtype] = DomainVisualStyle(
                dflt.fill_color,
                dflt.border_color,
                dflt.border_width,
                dflt.text_color,
                dflt.width,
                dflt.height,
                dflt.shape,
            )
            continue
        if not isinstance(entry, dict):
            raise ValueError(f"Style for '{dtype}' must be a mapping.")
        shape = entry.get("shape", dflt.shape)
        if shape not in VALID_SHAPES:
            raise ValueError(f"Invalid shape '{shape}' for '{dtype}'.")
        styles[dtype] = DomainVisualStyle(
            fill_color=_validate_hex_color(
                entry.get("fill_color", dflt.fill_color), f"{dtype}.fill_color"
            ),
            border_color=_validate_hex_color(
                entry.get("border_color", dflt.border_color), f"{dtype}.border_color"
            ),
            border_width=_validate_bounded_float(
                entry.get("border_width", dflt.border_width),
                f"{dtype}.border_width",
                0.5,
                8.0,
            ),
            text_color=_validate_hex_color(
                entry.get("text_color", dflt.text_color), f"{dtype}.text_color"
            ),
            width=_validate_bounded_float(
                entry.get("width", dflt.width), f"{dtype}.width", 32.0, 240.0
            ),
            height=_validate_bounded_float(
                entry.get("height", dflt.height), f"{dtype}.height", 16.0, 120.0
            ),
            shape=shape,
        )
    return SchematicTheme(
        name=name,
        orientation=orientation,
        domain_spacing=domain_spacing,
        chain_spacing=chain_spacing,
        styles=styles,
    )


def load_theme_presets() -> Dict[str, Any]:
    """Load built-in schematic theme presets from themes.yaml."""
    from fc_engineer.designer import DEFAULT_SCHEMATIC_THEME

    themes_path = os.path.join(_data_dir(), "themes.yaml")
    data = _load_yaml_dict(themes_path)
    presets: Dict[str, Any] = {}
    raw_presets = data.get("presets") if isinstance(data, dict) else None
    if isinstance(raw_presets, dict):
        for name, raw in raw_presets.items():
            if isinstance(name, str) and name.strip():
                try:
                    presets[name] = _parse_theme_dict(name, raw)
                except Exception as e:
                    log_sanitized_error(logger, f"Invalid built-in theme '{name}'", e)
    if DEFAULT_SCHEMATIC_THEME.name not in presets:
        presets[DEFAULT_SCHEMATIC_THEME.name] = DEFAULT_SCHEMATIC_THEME
    return presets


def parse_theme_yaml(content: str) -> Any:
    """Parse and validate a standalone .yaml schematic theme string."""
    if not isinstance(content, str):
        raise ValueError("Theme content must be a string.")
    if len(content.encode("utf-8")) > MAX_FILE_SIZE:
        raise ValueError("Theme content exceeds maximum 1MB size.")
    try:
        data = yaml.load(content, Loader=NoAliasSafeLoader)  # nosec B506
    except yaml.YAMLError as e:
        raise ValueError(f"Malformed theme YAML: {e}") from None
    if not isinstance(data, dict):
        raise ValueError("Theme YAML root must be a mapping.")
    name = str(data.get("name", "Custom Theme")).strip() or "Custom Theme"
    return _parse_theme_dict(name, data)


def serialize_theme_yaml(theme: Any) -> str:
    """Serialize a SchematicTheme to a standalone YAML string."""
    payload = {
        "schema_version": 1,
        "name": theme.name,
        "orientation": theme.orientation,
        "domain_spacing": float(theme.domain_spacing),
        "chain_spacing": float(theme.chain_spacing),
        "styles": {
            dtype: {
                "fill_color": st.fill_color,
                "border_color": st.border_color,
                "border_width": float(st.border_width),
                "text_color": st.text_color,
                "width": float(st.width),
                "height": float(st.height),
                "shape": st.shape,
            }
            for dtype, st in theme.styles.items()
        },
    }
    return yaml.safe_dump(payload, sort_keys=False, allow_unicode=True)


def load_theme_file(path: str) -> Any:
    """Load and validate a standalone .yaml/.yml theme file within 1MB bounds."""
    if not isinstance(path, str) or not path.strip():
        raise ValueError("Theme file path cannot be empty.")
    resolved = os.path.abspath(path)
    base_name = os.path.basename(resolved)
    if not resolved.lower().endswith((".yaml", ".yml")):
        raise ValueError("Theme file must have a .yaml or .yml extension.")
    if not os.path.isfile(resolved):
        raise ValueError(f"Theme file '{base_name}' does not exist.")
    if os.path.getsize(resolved) > MAX_FILE_SIZE:
        raise ValueError(f"Theme file '{base_name}' exceeds maximum 1MB size.")
    with open(resolved, "r", encoding="utf-8") as f:
        content = f.read(MAX_FILE_SIZE + 1)
    if len(content.encode("utf-8")) > MAX_FILE_SIZE:
        raise ValueError(f"Theme file '{base_name}' exceeds maximum 1MB size.")
    return parse_theme_yaml(content)


def save_theme_file(path: str, theme: Any) -> None:
    """Save a SchematicTheme to a standalone .yaml file."""
    if not isinstance(path, str) or not path.strip():
        raise ValueError("Theme file path cannot be empty.")
    resolved = os.path.abspath(path)
    if not resolved.lower().endswith((".yaml", ".yml")):
        resolved += ".yaml"
    payload = serialize_theme_yaml(theme)
    if len(payload.encode("utf-8")) > MAX_FILE_SIZE:
        raise ValueError("Serialized theme exceeds 1MB file limit.")
    with open(resolved, "w", encoding="utf-8") as f:
        f.write(payload)


def default_prefs_path() -> str:
    return os.path.join(os.path.expanduser("~"), ".fc_engineer", "prefs.yaml")


def load_user_prefs(path: Optional[str] = None) -> Dict[str, Any]:
    """Load user desktop preferences safely within 1MB YAML bounds."""
    target = os.path.abspath(path) if path else default_prefs_path()
    if not os.path.isfile(target):
        return {}
    try:
        if os.path.getsize(target) > MAX_FILE_SIZE:
            return {}
        with open(target, "r", encoding="utf-8") as f:
            raw = f.read(MAX_FILE_SIZE + 1)
        if len(raw.encode("utf-8")) > MAX_FILE_SIZE:
            return {}
        data = yaml.load(raw, Loader=NoAliasSafeLoader)  # nosec B506
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def save_user_prefs(prefs: Dict[str, Any], path: Optional[str] = None) -> None:
    """Persist user desktop preferences safely to YAML."""
    if not isinstance(prefs, dict):
        raise ValueError("Preferences must be a dictionary.")
    target = os.path.abspath(path) if path else default_prefs_path()
    parent = os.path.dirname(target)
    if parent:
        os.makedirs(parent, exist_ok=True)
    payload = yaml.safe_dump(prefs, sort_keys=False, allow_unicode=True)
    if len(payload.encode("utf-8")) > MAX_FILE_SIZE:
        raise ValueError("Preferences payload exceeds 1MB limit.")
    with open(target, "w", encoding="utf-8") as f:
        f.write(payload)
