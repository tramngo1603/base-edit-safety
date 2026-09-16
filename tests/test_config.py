"""Local external-tool configuration stays explicit and portable."""
from __future__ import annotations

import pytest

from base_edit_safety import config


def test_bedict_config_can_fall_back_to_ignored_project_env(tmp_path, monkeypatch):
    repo = tmp_path / "bedict"
    repo.mkdir()
    python = tmp_path / "bedict-python"
    python.write_text("")
    dotenv = tmp_path / ".env"
    dotenv.write_text(
        f'export BE_SAFETY_BEDICT_REPO="{repo}"\n'
        f'BE_SAFETY_BEDICT_PYTHON="{python}"\n'
    )
    monkeypatch.delenv("BE_SAFETY_BEDICT_REPO", raising=False)
    monkeypatch.delenv("BE_SAFETY_BEDICT_PYTHON", raising=False)
    monkeypatch.setattr(config, "_DOTENV_PATH", dotenv)

    loaded = config.load_bedict_config()

    assert loaded.repo_root == str(repo)
    assert loaded.python_executable == str(python)
    assert config.bedict_configured() is True


def test_process_environment_takes_precedence_over_project_env(tmp_path, monkeypatch):
    dotenv_repo = tmp_path / "dotenv-repo"
    dotenv_repo.mkdir()
    dotenv_python = tmp_path / "dotenv-python"
    dotenv_python.write_text("")
    env_repo = tmp_path / "env-repo"
    env_repo.mkdir()
    env_python = tmp_path / "env-python"
    env_python.write_text("")
    dotenv = tmp_path / ".env"
    dotenv.write_text(
        f"BE_SAFETY_BEDICT_REPO={dotenv_repo}\n"
        f"BE_SAFETY_BEDICT_PYTHON={dotenv_python}\n"
    )
    monkeypatch.setattr(config, "_DOTENV_PATH", dotenv)
    monkeypatch.setenv("BE_SAFETY_BEDICT_REPO", str(env_repo))
    monkeypatch.setenv("BE_SAFETY_BEDICT_PYTHON", str(env_python))

    loaded = config.load_bedict_config()

    assert loaded.repo_root == str(env_repo)
    assert loaded.python_executable == str(env_python)


def test_crispron_abe_config_can_fall_back_to_ignored_project_env(tmp_path, monkeypatch):
    repo = tmp_path / "crispron-BE"
    (repo / "bin").mkdir(parents=True)
    (repo / "bin" / "CRISPRonBE.sh").write_text("")
    python = tmp_path / "crispron-python"
    python.write_text("")
    dotenv = tmp_path / ".env"
    dotenv.write_text(
        f"BE_SAFETY_CRISPRON_ABE_REPO={repo}\n"
        f"BE_SAFETY_CRISPRON_ABE_PYTHON={python}\n"
    )
    monkeypatch.delenv("BE_SAFETY_CRISPRON_ABE_REPO", raising=False)
    monkeypatch.delenv("BE_SAFETY_CRISPRON_ABE_PYTHON", raising=False)
    monkeypatch.setattr(config, "_DOTENV_PATH", dotenv)

    loaded = config.load_crispron_abe_config()

    assert loaded.repo_root == str(repo)
    assert loaded.python_executable == str(python)
    assert config.crispron_abe_configured() is True


def test_crispor_config_can_fall_back_to_ignored_project_env(tmp_path, monkeypatch):
    repo = tmp_path / "crispor"
    repo.mkdir()
    (repo / "crispor.py").write_text("")
    python = tmp_path / "crispor-python"
    python.write_text("")
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    dotenv = tmp_path / ".env"
    dotenv.write_text(
        f"BE_SAFETY_CRISPOR_REPO={repo}\n"
        f"BE_SAFETY_CRISPOR_PYTHON={python}\n"
        "BE_SAFETY_CRISPOR_GENOME=hg38\n"
        f"BE_SAFETY_CRISPOR_GENOMEDIR={genome_dir}\n"
    )
    for name in (
        "BE_SAFETY_CRISPOR_REPO",
        "BE_SAFETY_CRISPOR_PYTHON",
        "BE_SAFETY_CRISPOR_GENOME",
        "BE_SAFETY_CRISPOR_GENOMEDIR",
    ):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(config, "_DOTENV_PATH", dotenv)

    loaded = config.load_crispor_config()

    assert loaded.repo_root == str(repo)
    assert loaded.python_executable == str(python)
    assert loaded.genome == "hg38"
    assert loaded.genome_dir == str(genome_dir)
    assert loaded.pam == "NGG"
    assert config.crispor_configured() is True


@pytest.mark.parametrize("pam", ["NG", "NGA", "NGAG", "NGN"])
def test_crispor_pam_is_loaded_without_narrowing_supported_patterns(tmp_path, monkeypatch, pam):
    repo = tmp_path / "crispor"
    repo.mkdir()
    (repo / "crispor.py").write_text("")
    python = tmp_path / "crispor-python"
    python.write_text("")
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    monkeypatch.setenv("BE_SAFETY_CRISPOR_REPO", str(repo))
    monkeypatch.setenv("BE_SAFETY_CRISPOR_PYTHON", str(python))
    monkeypatch.setenv("BE_SAFETY_CRISPOR_GENOME", "hg38")
    monkeypatch.setenv("BE_SAFETY_CRISPOR_GENOMEDIR", str(genome_dir))
    monkeypatch.setenv("BE_SAFETY_CRISPOR_PAM", pam.lower())

    loaded = config.load_crispor_config()

    assert loaded.pam == pam
