"""Dispatch ``python -m base_edit_safety.evaluation {tuning,held-out}``."""

from __future__ import annotations

import argparse
import sys
from typing import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] in {"-h", "--help"}:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("mode", choices=("tuning", "held-out"))
        parser.print_help()
        return 0
    mode, remainder = arguments[0], arguments[1:]
    if mode == "tuning":
        from .tuning import main as selected
    elif mode == "held-out":
        from .held_out import main as selected
    else:
        raise SystemExit(f"unknown evaluation mode {mode!r}; choose 'tuning' or 'held-out'")
    return selected(remainder)


if __name__ == "__main__":
    raise SystemExit(main())
