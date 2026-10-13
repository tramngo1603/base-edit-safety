"""Write the readable HTML report for a completed search."""
import html
from pathlib import Path

from .lab_view import STATUS_LABELS, coverage_facts, donor_scope, result_status

STYLE = """
:root{color:#202c38;background:#fff;font-family:'DM Sans',-apple-system,'Helvetica Neue',Arial,sans-serif;
font-size:16px;font-variant-numeric:tabular-nums}*{box-sizing:border-box}body{margin:0;line-height:1.6}
.shell{max-width:1160px;margin:auto;padding:0 36px}header{background:#f2f0ec;border-bottom:1px solid #dcdad4;padding:26px 0}
.brand{font-size:14px;font-weight:600;letter-spacing:.5px;color:#66717c;margin-bottom:12px}
h1,h2,h3{color:#172b43;line-height:1.25;margin:0;font-weight:600}h1{font-size:34px;letter-spacing:-.8px}
h2{font-size:22px;margin:0 0 16px}h3{font-size:17px;margin:16px 0 8px}p{margin:8px 0 16px}
main.shell{padding-top:24px;padding-bottom:44px}.heading{display:flex;align-items:center;justify-content:space-between;gap:18px}
.muted{color:#64717e}.meta{margin:10px 0 0}.status,.tag{display:inline-block;color:#27465d;background:#e8eff6;
font-size:14px;padding:4px 10px;border-radius:4px}.warn{color:#76591f;background:#f6f0e1;padding:12px 16px;border-left:3px solid #b87a2f}
.metrics{display:grid;grid-template-columns:repeat(3,1fr);gap:30px;padding:0 0 22px;margin-bottom:20px;border-bottom:1px solid #e2e5e8}
.metric strong{display:block;font-size:34px;line-height:1.25;font-weight:600;letter-spacing:-.7px;color:#172b43}
.metric span{display:block;color:#5e6872;font-size:15px;margin-top:6px}.scope{background:#eef3f9;border-left:3px solid #2b4162;
color:#33465b;border-radius:3px;padding:13px 16px;font-size:15px;margin:0 0 24px}.scope a{float:right;margin-left:14px}
.columns{display:grid;grid-template-columns:minmax(0,1.1fr) minmax(0,1fr);gap:48px;margin-bottom:26px}
table{border-collapse:collapse;width:100%;font-size:14px;margin:0 0 14px}th{background:#f5f6f7;color:#526274;text-align:left;font-weight:600;
padding:11px 14px;white-space:nowrap}td{padding:10px 14px;border-bottom:1px solid #e9edf0;vertical-align:top}
.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}.table-scroll{overflow:auto;border:1px solid #e0e5e9;border-radius:5px;margin-top:16px}
.table-scroll table{margin:0}.table-scroll th{white-space:normal}.table-scroll td{white-space:nowrap}code,.seq{font-family:Menlo,monospace;font-size:13px}
a{color:#2b4162;text-underline-offset:3px}a:focus-visible,summary:focus-visible{outline:2px solid #315679;outline-offset:4px}
.files{list-style:none;padding:0;margin:0}.files li{border-bottom:1px solid #e9edf0;padding:10px 0}.files a{font-weight:500}
.files small{display:block;color:#64717e;font-size:13px}details{border-top:1px solid #e0e5e9;padding:12px 0}
summary{color:#2b4162;font-size:16px;font-weight:600;cursor:pointer;padding:4px 0}details>div{padding-top:12px}
dl{display:grid;grid-template-columns:240px minmax(0,1fr);gap:10px 24px;font-size:14px}dt{color:#64717e}dd{margin:0;overflow-wrap:anywhere}
.note{font-size:14px;color:#64717e}.rule{border-top:1px solid #e0e5e9;padding-top:18px;margin-top:26px}
@media(max-width:750px){.shell{padding-left:20px;padding-right:20px}.columns{grid-template-columns:1fr;gap:26px}.heading{align-items:flex-start;flex-direction:column;gap:10px}.metrics{gap:16px}.metric strong{font-size:27px}.metric span{font-size:13px}.scope a{float:none;display:block;margin:6px 0 0}dl{grid-template-columns:1fr;gap:4px}dd{margin-bottom:10px}}
@media print{header{background:none}.shell{max-width:none;padding:0}details{break-inside:avoid}a{color:inherit}.scope a{display:none}}
"""


def render_coverage_table(rows):
    """Show exact denominators beside comparable coverage and copy counts."""
    parts = ['<div class="table-scroll"><table><thead><tr><th>Chromosome</th><th class=num>Donor modeled / reference bases</th>'
             '<th class=num>Comparable / possible windows</th><th class=num>Donor copy windows</th>'
             '<th class=num>Applied / stars resolved / unsupported records</th></tr></thead><tbody>']
    def count(value):
        return f"{value:,}" if isinstance(value, int) else "Not assessed"
    for row in rows:
        cells = [row["contig"],
                 count(row["modeled_reference_bases"]) + " / " + count(row["reference_bases"]),
                 count(row["comparable_oriented_windows_per_guide"]) + " / " + count(row["reference_oriented_windows_per_guide"]),
                 count(row["sample_evaluated_copy_oriented_windows_per_guide"]),
                 " / ".join(count(row.get(key, 0)) for key in
                            ("applied_variant_records", "resolved_spanning_deletion_variant_records", "unsupported_variant_records"))]
        parts.append('<tr><td class=seq>' + '</td><td class=num>'.join(html.escape(str(value)) for value in cells) + '</td></tr>')
    parts.append('</tbody></table></div><p class=note>Windows include both strands, per guide, before PAM or mismatch filtering. '
                 'Donor copy windows include indel paths and have a different denominator. '
                 'Record dispositions may overlap between copies. Benchmark-region membership is separate from confident genotypes; see coverage.tsv.</p>')
    return "\n".join(parts)


def write_report(output, state, coverage, config):
    output = Path(output)
    esc = lambda value: html.escape(str(value))
    sample = config.get("sample", {}).get("id", "Reference only")
    verification = state.get("example_verification", {})
    check = ""
    if verification.get("status") == "matched":
        check = '<p class=note>Example verification: all expected scientific results matched.</p>'
    elif verification.get("status") == "mismatched":
        check = ('<p class=warn><strong>Example verification failed.</strong> ' + esc(verification.get("message", "")) +
                 ' Keep this result folder and check the release and input hashes before using the example results.</p>')
    metrics = "".join(f'<div class=metric><strong>{value}</strong><span>{label}</span></div>' for label, value in (
        ("Reference loci", f"{state['reference_loci']:,}"),
        ("Sample paths", f"{state['sample_loci']:,}" if config.get("sample") else "—"),
        ("Sequence configurations", f"{state['candidate_configurations']:,}")))
    changes = "".join(f'<tr><td>{esc(STATUS_LABELS.get(key, key.replace("_", " ")))}</td><td class=num>{value:,}</td></tr>'
                      for key, value in state["locus_changes"].items())
    main_files = {
        "candidates.tsv": "Candidate sequences & coordinates",
        "locus-changes.tsv": "Sample locus comparisons",
        "sequence-changes.tsv": "Changed target sequences",
        "exclusions.tsv": "Excluded regions & phase notes",
    }
    detail_files = {
        "coverage.tsv": "Per-chromosome coverage and fixed denominators",
        "coverage-regions.tsv": "Evaluated and excluded reference bases",
        "coverage-windows.tsv": "Comparable and excluded reference-window starts",
        "coverage.json": "Coverage definitions and complete accounting",
        "request.json": "Exact replay request and input identities",
        "run-manifest.json": "Completion, file checksums and measured resources",
    }
    if (output / "evidence.tsv").exists():
        main_files["evidence.tsv"] = "Linked experimental evidence"
    def links(files):
        return "".join(f'<li><a href="{name}">{esc(label)}</a><small>{name}</small></li>' for name, label in files.items())
    exclusions = "; ".join(f"{row['contig']}: {row['reason']}" for row in coverage.get("sample_scope_exclusions", [])) or "None"
    unknown = ", ".join(row["contig"] for row in coverage.get("sample_undeclared_ploidy_intervals", [])) or "None"
    facts = "".join(f'<dt>{esc(label)}</dt><dd>{esc(value)}</dd>' for label, value in coverage_facts(coverage, config))
    comparison = (f'<table><thead><tr><th>Locus status</th><th class=num>Count</th></tr></thead><tbody>{changes}</tbody></table>'
                  f'<p><strong>{state["changed_sequence_loci"]:,}</strong> retained loci have changed target sequences.</p>'
                  '<p class=note>Changed sequences stay distinct from locus gains and losses. Unresolved and unpaired paths retain their status.</p>') if config.get("sample") else (
                  '<p>Reference-only analysis. Candidate review and export are available without a VCF.</p>'
                  '<p class=note>To add a sample, start a new analysis from the app. The original results remain intact.</p>')
    status_class = "status" if state.get("status") == "complete" else "warn"
    (output / "report.html").write_text(f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Base Edit Safety results</title><style>{STYLE}</style></head><body>
<header><div class=shell><div class=brand>BASE EDIT SAFETY</div><div class=heading><h1>Review candidates</h1>
<span class="{status_class}">{esc(result_status(state))}</span></div>
<p class="meta muted">{esc(sample)} &nbsp;·&nbsp; {coverage['guide_count']} guides &nbsp;·&nbsp; {esc(config['assembly'])}
&nbsp;·&nbsp; {len(coverage['searched_contigs'])} chromosomes</p></div></header>
<main class=shell><div class=metrics>{metrics}</div>
<div class=scope>{esc(donor_scope(coverage, config))}<a href="#coverage" onclick="document.getElementById('coverage').open=true">Coverage details ↓</a></div>{check}
<div class=columns><section><h2>{'Sample differences' if config.get('sample') else 'Reference candidates'}</h2>{comparison}</section>
<section><h2>Candidates for follow-up</h2><ul class=files>{links(main_files)}</ul>
<p class=note>These links open complete result tables. For a follow-up export with source details and a field dictionary, open this results folder in the app and choose Export…. All and matching exports include rows beyond the display limit.</p></section></div>
<details id=coverage><summary>Coverage & exclusions</summary><div><dl>{facts}</dl>
<p>Processing a chromosome does not mean every region was resolved.</p>
{render_coverage_table(coverage.get('chromosomes', [])) if coverage.get('chromosomes') and config.get('sample') else ''}
<p>Sample accounting: {esc(coverage.get('sample_status', 'Not recorded').replace('_', ' '))}.</p>
<p>Sample chromosome exclusions: {esc(exclusions)}. Undeclared sample ploidy: {esc(unknown)}.</p>
<p>Other excluded scope: {esc(', '.join(coverage.get('unsupported_scope', [])))}.</p>
<p class=note>Sample sequence at absent VCF records is conditional on the retained callset. Local copy labels do not establish long-range phase. Missing calls and incompatible phase remain unresolved.</p></div></details>
<details><summary>Search settings & coordinate conventions</summary><div>
<p>20-base guides · NGG PAM · both strands · up to {config.get('max_mismatches', 6)} mismatches · no gaps.</p>
<p>Coordinates start at zero and exclude the end position. Sequences include the PAM in guide orientation.</p></div></details>
<details><summary>Files, input identities & reproduction</summary><div><ul class=files>{links(detail_files)}</ul>
<p class=note>Keep this folder and its input files together. Timestamps and memory/runtime measurements belong to this run and are separate from the expected scientific results.</p></div></details>
<p class="note rule">Base Edit Safety · local sequence analysis</p></main></body></html>''', encoding="utf-8")
