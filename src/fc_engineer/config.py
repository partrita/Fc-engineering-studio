import sys
import os
import logging
import traceback
from typing import Dict, List, Tuple

import yaml

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def log_sanitized_error(logger_obj, msg: str, e: Exception):
    base_path = os.path.dirname(__file__)
    safe_err = str(e).replace(base_path, '.') if base_path else str(e)
    safe_traceback = traceback.format_exc().replace(base_path, '.') if base_path else traceback.format_exc()
    logger_obj.error(f"{msg}: {safe_err}\n{safe_traceback}")

# --- Configuration & Data Loading ---

class NoAliasSafeLoader(yaml.SafeLoader):
    MAX_DEPTH = 50

    def compose_node(self, parent, index):
        if self.check_event(yaml.events.AliasEvent):
            raise yaml.constructor.ConstructorError("Aliases are not allowed to prevent YAML bomb (DoS) attacks.")

        if not hasattr(self, "depth"):
            self.depth = 0
        self.depth += 1

        if self.depth > self.MAX_DEPTH:
            raise yaml.constructor.ConstructorError("Maximum YAML nesting depth exceeded to prevent DoS attacks.")

        node = super().compose_node(parent, index)
        self.depth -= 1
        return node

def load_yaml_data():
    base_path = os.path.join(os.path.dirname(__file__), "data")
    seq_path = os.path.join(base_path, "sequences.yaml")
    mut_path = os.path.join(base_path, "mutants.yaml")

    isotypes = {}
    common_muts = []

    # SECURITY: Max file size of 1MB to prevent DoS via memory exhaustion
    MAX_FILE_SIZE = 1 * 1024 * 1024

    try:
        seq_name = os.path.basename(seq_path)
        if os.path.isfile(seq_path):
            if os.path.getsize(seq_path) > MAX_FILE_SIZE:
                print(f"Error: {seq_name} exceeds 1MB limit.", file=sys.stderr)
                raise ValueError(f"File {seq_name} exceeds maximum size of 1MB")
            with open(seq_path, "r", encoding="utf-8") as f:
                content = f.read(MAX_FILE_SIZE + 1)
                if len(content) > MAX_FILE_SIZE:
                    print(f"Error: {seq_name} content exceeds 1MB limit.", file=sys.stderr)
                    raise ValueError(f"File {seq_name} content exceeds maximum size of 1MB")
                data = yaml.load(content, Loader=NoAliasSafeLoader)  # nosec B506
                if isinstance(data, dict):
                    val = data.get("isotypes")
                    isotypes = val if (
                        isinstance(val, dict) and all(
                            isinstance(k, str) and isinstance(v, dict) and all(isinstance(ik, str) and isinstance(iv, str) for ik, iv in v.items())
                            for k, v in val.items()
                        )
                    ) else {}
                else:
                    print(f"Error: Parsed data from {seq_name} is not a dictionary.", file=sys.stderr)
        else:
            print(f"Error: Missing configuration file {seq_name}.", file=sys.stderr)
    except Exception as e:
        safe_err = str(e).replace(base_path, '.') if base_path else str(e)
        safe_traceback = traceback.format_exc().replace(base_path, '.') if base_path else traceback.format_exc()
        logger.error(f"Error loading {seq_name}: {safe_err}\n{safe_traceback}")
        print(f"Error loading {seq_name}: An unexpected error occurred.", file=sys.stderr)

    try:
        mut_name = os.path.basename(mut_path)
        if os.path.isfile(mut_path):
            if os.path.getsize(mut_path) > MAX_FILE_SIZE:
                print(f"Error: {mut_name} exceeds 1MB limit.", file=sys.stderr)
                raise ValueError(f"File {mut_name} exceeds maximum size of 1MB")
            with open(mut_path, "r", encoding="utf-8") as f:
                content = f.read(MAX_FILE_SIZE + 1)
                if len(content) > MAX_FILE_SIZE:
                    print(f"Error: {mut_name} content exceeds 1MB limit.", file=sys.stderr)
                    raise ValueError(f"File {mut_name} content exceeds maximum size of 1MB")
                data = yaml.load(content, Loader=NoAliasSafeLoader)  # nosec B506
                if isinstance(data, dict):
                    val = data.get("common_mutations")
                    common_muts = val if (
                        isinstance(val, list) and all(
                            isinstance(item, dict) and isinstance(item.get("value"), str) and isinstance(item.get("label"), str)
                            for item in val
                        )
                    ) else []
                else:
                    print(f"Error: Parsed data from {mut_name} is not a dictionary.", file=sys.stderr)
        else:
            print(f"Error: Missing configuration file {mut_name}.", file=sys.stderr)
    except Exception as e:
        safe_err = str(e).replace(base_path, '.') if base_path else str(e)
        safe_traceback = traceback.format_exc().replace(base_path, '.') if base_path else traceback.format_exc()
        logger.error(f"Error loading {mut_name}: {safe_err}\n{safe_traceback}")
        print(f"Error loading {mut_name}: An unexpected error occurred.", file=sys.stderr)

    return isotypes, common_muts

SEQUENCES, COMMON_MUTATIONS = load_yaml_data()
