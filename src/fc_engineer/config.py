import logging
import os
import re
import sys
import traceback
from typing import Any, Dict, List, Tuple

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
    return os.path.join(os.path.dirname(__file__), "data")


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
