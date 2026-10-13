#!/usr/bin/env python3
"""Assemble a local handoff from a named, clean private source revision."""
import argparse
from datetime import datetime, timezone
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile

ALLOWLIST = ["README.md", "LICENSE", "SECURITY.md", "pyproject.toml", "MANIFEST.in",
             ".gitignore", ".github",
             ".env.example", "requirements-external.txt", "data/README.md", "src", "tests", "examples",
             "docs", "tools", "handoff"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            h.update(block)
    return h.hexdigest()


def checksums(root):
    (root / "SHA256SUMS").write_text("".join(
        f"{sha(p)}  {p.relative_to(root)}\n" for p in sorted(root.rglob("*"))
        if p.is_file() and p.name != "SHA256SUMS" and not any(part.startswith(".") for part in p.relative_to(root).parts)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, assets, output = [p.resolve() for p in (args.source, args.assets, args.output)]
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=source).strip():
        raise ValueError("Commit the reviewed source before building a handoff.")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    output.mkdir(parents=True, exist_ok=False)
    pins = json.loads((source / "handoff/dependencies.json").read_text())
    vendor = output / "vendor"
    vendor.mkdir()
    (vendor / "wheels").mkdir()
    for name, record, path in [
        ("python-runtime.tar.gz", pins["python"], assets / "python-runtime.tar.gz"),
        ("cas-offinder-source.tar.gz", pins["cas_offinder"], assets / "cas-offinder-source.tar.gz"),
    ]:
        if sha(path) != record["sha256"]:
            raise ValueError("Asset checksum mismatch: " + name)
        shutil.copyfile(path, vendor / name)
    for record in pins["wheels"].values():
        path = assets / "wheels" / record["filename"]
        if sha(path) != record["sha256"]:
            raise ValueError("Dependency checksum mismatch: " + path.name)
        shutil.copyfile(path, vendor / "wheels" / path.name)
    if sha(source / "tools/cas-offinder-2.4.1-compat.patch") != pins["cas_offinder"]["patch_sha256"]:
        raise ValueError("Backend patch differs from its pin")
    with (output / "source.tar.gz").open("wb") as handle:
        subprocess.run(["git", "archive", "--format=tar.gz", "--prefix=source/", revision, *ALLOWLIST],
                       cwd=source, stdout=handle, check=True)
    with tempfile.TemporaryDirectory(prefix="bes-source-build-") as tmp:
        with tarfile.open(output / "source.tar.gz") as archive:
            archive.extractall(tmp, filter="data")
        subprocess.run([sys.executable, "-I", "-m", "pip", "--isolated", "wheel", "--no-deps",
                        "--no-build-isolation", "--no-index", "--wheel-dir", str(vendor / "wheels"),
                        str(Path(tmp) / "source")], check=True)
    wheel, = (vendor / "wheels").glob("base_edit_safety-*.whl")
    with zipfile.ZipFile(wheel) as archive:
        metadata_path, = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        version = BytesParser().parsebytes(archive.read(metadata_path))["Version"]
    for name in ("Start.command", "build-backend.sh", "QUICKSTART.html", "ADVANCED.md", "PUBLIC_INPUTS.md",
                 "dependencies.json", "build-requirements.lock", "expected-hg003.json"):
        shutil.copy2(source / "handoff" / name, output / name)
    shutil.copyfile(source / "LICENSE", output / "LICENSE")
    safety_sentence = (source / "README.md").read_text().splitlines()[2]
    (output / "README.md").write_text("# Base Edit Safety\n\n" + safety_sentence +
        "\n\nDouble-click **Start.command**. Follow [the quickstart](QUICKSTART.html) for installation, "
        "examples and your own inputs. [Advanced instructions](ADVANCED.md) cover command-line use "
        "and rebuilding. The complete source and documentation are in source.tar.gz.\n")
    shutil.copyfile(source / "tools/cas-offinder-2.4.1-compat.patch", output / "cas-offinder-2.4.1-compat.patch")
    (output / "requirements.lock").write_text((source / "handoff/runtime-requirements.lock").read_text() +
        f"base-edit-safety=={version} --hash=sha256:{sha(wheel)}\n")
    (output / "backend").mkdir()
    binary = output / "backend/cas-offinder"
    shutil.copy2(assets / "backend-build/cas-offinder", binary)
    binary.chmod(0o755)
    licenses = output / "licenses"
    licenses.mkdir()
    shutil.copyfile(source / "tools/Cas-OFFinder-LICENSE.txt", licenses / "Cas-OFFinder-LICENSE.txt")
    shutil.copytree(assets / "python-notices/python/licenses", licenses / "python")
    for archive in (vendor / "wheels").glob("*.whl"):
        with zipfile.ZipFile(archive) as z:
            for name in z.namelist():
                if not name.endswith("/") and (".dist-info/licenses/" in name or Path(name).name.upper().startswith(("LICENSE", "COPYING", "NOTICE", "AUTHORS"))):
                    target = licenses / archive.name / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(z.read(name))
    release = {
        "release_version": version, "release_revision": revision,
        "cleanup_revision": revision,
        "previous_release_revision": "72e50f74d5f604b800a9891b594a98648989894e",
        "source_clean": True, "created_utc": datetime.now(timezone.utc).isoformat(),
        "allowlist": ALLOWLIST, "source_archive_sha256": sha(output / "source.tar.gz"),
        "wheel": {"path": str(wheel.relative_to(output)), "sha256": sha(wheel), "bytes": wheel.stat().st_size},
        "backend": {"path": "backend/cas-offinder", "sha256": sha(binary), **pins["cas_offinder"],
                    "compiler": "Apple clang 15.0.0 (clang-1500.3.9.4)", "cmake": "4.0.3",
                    "build_command": "./build-backend.sh NEW_BUILD_FOLDER", "device": "Apple OpenCL GPU G0"},
        "validation": {"status": "pending"}, "public_repository_updated": False,
        "private_data_in_bundle": False,
    }
    # Source archive and executable have distinct hashes.
    release["backend"]["source_archive_sha256"] = pins["cas_offinder"]["sha256"]
    release["backend"]["sha256"] = sha(binary)
    (output / "RELEASE.json").write_text(json.dumps(release, indent=2, sort_keys=True) + "\n")
    checksums(output)
    print(output)


if __name__ == "__main__":
    main()
