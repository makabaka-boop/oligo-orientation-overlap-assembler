"""Command line interface: JSON in (stdin or file), JSON out (stdout)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .core import AssemblyError, assemble, assembly_to_data, reads_from_data


def _fail(message: str, code: int) -> int:
    print(json.dumps({"error": message}), file=sys.stderr)
    return code


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="assembler",
        description=(
            "Assemble short DNA reads into the shortest string by choosing "
            "read order and strand, joining adjacent reads on their maximal "
            "exact suffix/prefix overlap."
        ),
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="JSON input file (default: read from stdin)",
    )
    parser.add_argument(
        "--pretty", action="store_true", help="pretty-print the JSON output"
    )
    args = parser.parse_args(argv)

    try:
        if args.input == "-":
            text = sys.stdin.read()
        else:
            text = Path(args.input).read_text(encoding="utf-8")
    except OSError as exc:
        return _fail(f"cannot read input: {exc}", 2)

    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        return _fail(f"invalid JSON: {exc}", 1)

    try:
        reads = reads_from_data(data)
    except AssemblyError as exc:
        return _fail(str(exc), 1)

    result = assembly_to_data(assemble(reads))
    json.dump(result, sys.stdout, indent=2 if args.pretty else None)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
