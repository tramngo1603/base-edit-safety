"""CLI entry point for repeatable tuning-set evaluation."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .runner import (
    DEFAULT_FIXTURE_ROOT,
    DEFAULT_HAN_OVERLAP_MAP,
    DEFAULT_PROTOCOL,
    DEFAULT_REGISTRY,
    REPO_ROOT,
    run_tuning,
)


def _tool_versions(values: list[str]) -> dict[str, str]:
    result = {}
    for value in values:
        if "=" not in value or not all(part.strip() for part in value.split("=", 1)):
            raise argparse.ArgumentTypeError("--tool-version must be NAME=VERSION")
        name, version = value.split("=", 1)
        result[name.strip()] = version.strip()
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--measurement", action="append", required=True, help="normalized ingested TSV/CSV/JSON; repeatable")
    value.add_argument("--pipeline-output", action="append", required=True, help="tool output JSON; repeatable")
    value.add_argument("--golden-set", required=True)
    value.add_argument("--output-dir", required=True)
    value.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    value.add_argument("--protocol", default=str(DEFAULT_PROTOCOL))
    value.add_argument("--han-overlap-map", default=str(DEFAULT_HAN_OVERLAP_MAP))
    value.add_argument("--repo-root", default=str(REPO_ROOT))
    value.add_argument("--fixture-mode", action="store_true")
    value.add_argument("--fixture-root", default=str(DEFAULT_FIXTURE_ROOT))
    value.add_argument("--tool-version", action="append", default=[], metavar="NAME=VERSION")
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    artifacts = run_tuning(
        measurement_paths=args.measurement,
        pipeline_output_paths=args.pipeline_output,
        golden_set_path=args.golden_set,
        output_dir=args.output_dir,
        registry_path=args.registry,
        protocol_path=args.protocol,
        han_overlap_map_path=args.han_overlap_map,
        fixture_mode=args.fixture_mode,
        fixture_root=args.fixture_root,
        repo_root=args.repo_root,
        tool_versions=_tool_versions(args.tool_version),
    )
    print(Path(artifacts.report_markdown))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
