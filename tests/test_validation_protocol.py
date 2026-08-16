"""Check the required sections and tables in the versioned validation protocol."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "validation" / "validation-protocol.md"


def test_validation_protocol_exists():
    assert PROTOCOL.is_file()


def test_validation_protocol_has_required_section_anchors():
    text = PROTOCOL.read_text()

    for heading in (
        "## 2. Dataset registry",
        "## 3. Tuning / held-out split (pre-registered)",
        "## 6. Freeze procedure",
        "## 11. Human sign-off",
    ):
        assert heading in text


def test_validation_protocol_section_2_retains_registry_column_spec():
    text = PROTOCOL.read_text()

    assert (
        "| ID | Source (citation) | Accession / table | Assay | Cell type | "
        "Nuclease/editor | Build | Guide(s) | Validates layer | Role | QC status | "
        "Known issues |"
    ) in text
