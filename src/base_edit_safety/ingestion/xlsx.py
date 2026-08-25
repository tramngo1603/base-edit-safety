"""Small, dependency-free reader for tabular cells in OOXML ``.xlsx`` workbooks.

The validation inputs are publisher-supplied workbooks, while the core project intentionally has
no Python dependencies.  This reader covers the cell encodings used by those workbooks (shared
strings, inline strings, booleans, cached formula values, and numeric/string cells).  Formatting,
charts, macros, and formula evaluation are outside the ingestion contract.
"""
from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile


_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
_REL_NS = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}
_RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
_CELL_RE = re.compile(r"([A-Z]+)([0-9]+)$")


def _column_index(reference: str) -> int:
    match = _CELL_RE.fullmatch(reference)
    if not match:
        raise ValueError(f"invalid XLSX cell reference: {reference!r}")
    value = 0
    for char in match.group(1):
        value = value * 26 + ord(char) - ord("A") + 1
    return value - 1


def _text(node: ET.Element) -> str:
    """Join plain and rich-text runs in an OOXML string node."""
    return "".join(part.text or "" for part in node.findall(".//m:t", _NS))


def _shared_strings(book: ZipFile) -> list[str]:
    try:
        root = ET.fromstring(book.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return [_text(item) for item in root.findall("m:si", _NS)]


def _sheet_paths(book: ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(book.read("xl/workbook.xml"))
    relationships = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
    targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in relationships.findall("r:Relationship", _REL_NS)
    }
    paths: dict[str, str] = {}
    for sheet in workbook.findall("m:sheets/m:sheet", _NS):
        target = targets[sheet.attrib[_RID]].lstrip("/")
        if not target.startswith("xl/"):
            target = f"xl/{target}"
        paths[sheet.attrib["name"]] = target
    return paths


def sheet_names(path: str | Path) -> list[str]:
    """Return workbook sheet names in publisher order."""
    with ZipFile(path) as book:
        return list(_sheet_paths(book))


def iter_rows(path: str | Path, sheet_name: str) -> Iterator[list[str]]:
    """Yield rows with blank cells preserved up to each row's rightmost populated cell."""
    with ZipFile(path) as book:
        paths = _sheet_paths(book)
        if sheet_name not in paths:
            raise KeyError(f"sheet {sheet_name!r} not found; available: {list(paths)}")
        shared = _shared_strings(book)
        root = ET.fromstring(book.read(paths[sheet_name]))
        for row_node in root.findall("m:sheetData/m:row", _NS):
            values: dict[int, str] = {}
            for cell in row_node.findall("m:c", _NS):
                column = _column_index(cell.attrib["r"])
                kind = cell.attrib.get("t")
                value_node = cell.find("m:v", _NS)
                raw = "" if value_node is None or value_node.text is None else value_node.text
                if kind == "s":
                    value = shared[int(raw)]
                elif kind == "inlineStr":
                    inline = cell.find("m:is", _NS)
                    value = "" if inline is None else _text(inline)
                elif kind == "b":
                    value = "true" if raw == "1" else "false"
                else:
                    value = raw
                values[column] = value
            if values:
                width = max(values) + 1
                yield [values.get(column, "") for column in range(width)]
            else:
                yield []
