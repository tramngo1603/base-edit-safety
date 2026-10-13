"""Remember the chosen results folder and check that it can accept a run."""
import json
import os
from pathlib import Path
import tempfile


def preferences_path():
    return Path.home() / "Library/Application Support/Base Edit Safety/preferences.json"


def read_location(path):
    path = Path(path)
    try:
        value = json.loads(path.read_text())["results_folder"]
        if not isinstance(value, str) or not value or not Path(value).is_absolute():
            raise ValueError("Invalid folder record")
        return value
    except FileNotFoundError:
        return ""
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise ValueError("The remembered save location could not be read. Choose a results folder again.") from exc


def check_location(folder):
    if not str(folder).strip():
        raise ValueError("Choose a save location before running.")
    path = Path(folder).expanduser()
    if not path.is_dir():
        raise ValueError(f"Save location is unavailable: {path}. Reconnect it or choose another folder.")
    if not os.access(path, os.W_OK | os.X_OK):
        raise ValueError(f"Cannot write to {path}. Check its permissions or choose another folder.")
    try:
        with tempfile.TemporaryFile(dir=path) as probe:
            probe.write(b"BES folder check\n")
            probe.flush()
    except OSError as exc:
        raise ValueError(f"Cannot write to {path}: {exc.strerror or exc}. Choose another folder.") from exc
    return path.resolve()


def remember_location(folder, preferences):
    path = Path(preferences)
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".preferences-", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump({"results_folder": str(Path(folder).expanduser().resolve())}, handle, indent=2)
            handle.write("\n")
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
