#!/usr/bin/env python3
"""
fungar_report.py — Generate a self-contained HTML report from FUNGAR output.

Usage (called automatically by FUNGAR.sh, or standalone):
    python3 fungar_report.py \\
        --sample-id SAMPLE01 \\
        --species Aspergillus_fumigatus \\
        --results path/to/SAMPLE01_results.csv \\
        --summary path/to/SAMPLE01_summary.csv \\
        --outdir  path/to/output/ \\
        [--version 2.0] \\
        [--runtime 42s] \\
        [--params "threads=4, min-support=1"]
"""

import argparse
import os
import re
import sys
from datetime import datetime

try:
    import pandas as pd
except ImportError:
    sys.exit("❌ pandas is required: pip install pandas")

# ---------------------------------------------------------------------------
# Drug classification
# ---------------------------------------------------------------------------
DRUG_CLASSES = {
    "azole": {
        "keywords": [
            "itraconazole", "voriconazole", "posaconazole", "isavuconazole",
            "fluconazole", "ravuconazole", "clotrimazole", "miconazole",
            "ketoconazole", "bromuconazole", "difenoconazole", "epoxiconazole",
            "imazalil", "prochloraz", "propiconazole", "tebuconazole",
            "azole",
        ],
        "label": "Azole",
        "color": "#f97316",   # orange
        "bg": "rgba(249,115,22,0.12)",
    },
    "echinocandin": {
        "keywords": ["caspofungin", "micafungin", "anidulafungin", "echinocandin"],
        "label": "Echinocandin",
        "color": "#06b6d4",   # cyan
        "bg": "rgba(6,182,212,0.12)",
    },
    "polyene": {
        "keywords": ["amphotericin_b", "amphotericin b", "amphotericin", "polyene", "nystatin"],
        "label": "Polyene",
        "color": "#ec4899",   # pink
        "bg": "rgba(236,72,153,0.12)",
    },
    "qoi": {
        "keywords": ["azoxystrobin", "strobilurin", "qoi"],
        "label": "QoI Fungicide",
        "color": "#84cc16",   # lime
        "bg": "rgba(132,204,22,0.12)",
    },
    "sdhi": {
        "keywords": ["boscalid", "fluxapyroxad", "isopyrazam", "penthiopyrad", "sdhi"],
        "label": "SDHI",
        "color": "#a855f7",   # purple
        "bg": "rgba(168,85,247,0.12)",
    },
}

SOURCE_BADGE = {
    "clinical":     ('<span class="badge badge-clinical">Clinical</span>', "#22c55e"),
    "agricultural": ('<span class="badge badge-agri">Agricultural</span>', "#eab308"),
    "both":         ('<span class="badge badge-both">Both</span>', "#6366f1"),
    "unknown":      ('<span class="badge badge-unknown">Unknown</span>', "#6b7280"),
}


def classify_drug(fungicide_str: str):
    raw = str(fungicide_str).lower()
    for cls, info in DRUG_CLASSES.items():
        if any(kw in raw for kw in info["keywords"]):
            return info
    return {"label": "Other", "color": "#6b7280", "bg": "rgba(107,114,128,0.12)"}


def mutation_notation(row: dict) -> str:
    """p.Ref{Pos}Mut  e.g. p.L98H"""
    return f"p.{row.get('Reference','?')}{row.get('Position','?')}{row.get('Mutation','?')}"


def gene_short(gene: str) -> str:
    """Strip accession suffix: Cyp51A_AGH55425.1 → Cyp51A"""
    parts = gene.split("_")
    if len(parts) > 1:
        return "_".join(parts[:-1])
    return gene


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------
def build_html(
    sample_id: str,
    species: str,
    results_df,
    summary_df,
    version: str,
    runtime: str,
    params: str,
    generated_at: str,
) -> str:

    # ---- aggregate stats ----
    n_mutations = len(summary_df)
    n_genes = summary_df["Gene"].nunique() if not summary_df.empty else 0
    total_support = int(summary_df["Support_Reads"].sum()) if not summary_df.empty else 0
    total_pairs = int(summary_df["Support_Pairs"].sum()) if not summary_df.empty and "Support_Pairs" in summary_df.columns else total_support
    n_reads = len(results_df) if not results_df.empty else 0

    # ---- build summary table rows ----
    summary_rows_html = ""
    chart_labels = []
    chart_data = []
    chart_colors = []

    if not summary_df.empty:
        for _, row in summary_df.sort_values("Support_Reads", ascending=False).iterrows():
            drug_info = classify_drug(str(row.get("Fungicide", "")))
            src = str(row.get("Source", "unknown")).lower()
            badge_html, _ = SOURCE_BADGE.get(src, SOURCE_BADGE["unknown"])
            gene_s = gene_short(str(row["Gene"]))
            notation = mutation_notation(row)
            support_reads = int(row.get("Support_Reads", 0))
            support_pairs = int(row.get("Support_Pairs", support_reads)) if "Support_Pairs" in summary_df.columns else support_reads

            species_val = str(row.get("Species", species)).replace("_", " ")

            if "Support_Pairs" in summary_df.columns:
                display_support = f"{support_pairs} pairs / {support_reads} reads"
                chart_val = support_pairs
            else:
                display_support = f"{support_reads} reads"
                chart_val = support_reads

            summary_rows_html += f"""
            <tr>
              <td><strong>{gene_s}</strong></td>
              <td><code>{notation}</code></td>
              <td>{int(row.get('Position', 0))}</td>
              <td><span class="drug-badge" style="background:{drug_info['bg']};color:{drug_info['color']};border:1px solid {drug_info['color']}33">{drug_info['label']}</span></td>
              <td style="max-width:220px;word-break:break-word;font-size:0.78rem">{row.get('Fungicide','')}</td>
              <td>{badge_html}</td>
              <td style="font-size:0.78rem"><em>{species_val}</em></td>
              <td><span class="support-bar-wrap"><span class="support-bar" style="width:min(100%,{min(support_reads*4,100)}px);background:{drug_info['color']}"></span><span class="support-num">{display_support}</span></span></td>
            </tr>"""

            chart_labels.append(f"{gene_s} {notation}")
            chart_data.append(chart_val)
            chart_colors.append(drug_info["color"])
    else:
        summary_rows_html = '<tr><td colspan="8" class="empty-row">No mutations detected above the minimum read-support threshold.</td></tr>'

    # ---- build detailed results table rows (first 200 rows) ----
    detail_rows_html = ""
    if not results_df.empty:
        for _, row in results_df.head(200).iterrows():
            gene_s = gene_short(str(row.get("Gene", "")))
            notation = mutation_notation(row)
            accession = str(row.get("Subject_Accession", ""))
            strand = str(row.get("Strand", ""))
            aln_file = str(row.get("Alignment", ""))
            src = str(row.get("Source", "unknown")).lower()
            badge_html, _ = SOURCE_BADGE.get(src, SOURCE_BADGE["unknown"])
            species_val = str(row.get("Species", species)).replace("_", " ")
            detail_rows_html += f"""
            <tr>
              <td style="font-size:0.78rem;word-break:break-all">{row.get('Sample','')}</td>
              <td><strong>{gene_s}</strong></td>
              <td><code>{notation}</code></td>
              <td style="font-size:0.78rem">{accession}</td>
              <td>{strand}</td>
              <td>{badge_html}</td>
              <td style="font-size:0.78rem"><em>{species_val}</em></td>
              <td style="font-size:0.78rem">{aln_file}</td>
            </tr>"""
    else:
        detail_rows_html = '<tr><td colspan="8" class="empty-row">No raw alignments to display.</td></tr>'

    # ---- chart.js datasets ----
    import json
    labels_js = json.dumps(chart_labels)
    data_js = json.dumps(chart_data)
    colors_js = json.dumps(chart_colors)

    # ---- assemble drug class legend ----
    legend_items = ""
    for info in DRUG_CLASSES.values():
        legend_items += f'<span class="legend-dot" style="background:{info["color"]}"></span>{info["label"]} &nbsp; '
    legend_items += '<span class="legend-dot" style="background:#6b7280"></span>Other'

    species_display = species.replace("_", " ")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>FUNGAR Report — {sample_id}</title>
<meta name="description" content="FUNGAR antifungal resistance report for sample {sample_id} ({species_display})"/>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet"/>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
  :root {{
    --bg:        #0d1117;
    --surface:   #161b22;
    --surface2:  #21262d;
    --border:    #30363d;
    --text:      #e6edf3;
    --text-muted:#7d8590;
    --accent:    #3fb950;
    --accent2:   #58a6ff;
    --header-h:  64px;
  }}
  *,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
  html{{scroll-behavior:smooth}}
  body{{
    font-family:'Inter',sans-serif;
    background:var(--bg);
    color:var(--text);
    min-height:100vh;
    line-height:1.6;
  }}
  /* ---- header ---- */
  header{{
    position:sticky;top:0;z-index:100;
    height:var(--header-h);
    background:rgba(13,17,23,0.85);
    backdrop-filter:blur(12px);
    border-bottom:1px solid var(--border);
    display:flex;align-items:center;justify-content:space-between;
    padding:0 2rem;
  }}
  .logo{{display:flex;align-items:center;gap:.6rem;font-weight:700;font-size:1.15rem;letter-spacing:-.02em}}
  .logo .mushroom{{font-size:1.4rem}}
  .logo .ver{{font-size:.7rem;font-weight:500;color:var(--accent);background:rgba(63,185,80,.12);padding:2px 7px;border-radius:20px;border:1px solid rgba(63,185,80,.3)}}
  .header-meta{{font-size:.8rem;color:var(--text-muted)}}
  .btn-print{{
    cursor:pointer;background:transparent;border:1px solid var(--border);
    color:var(--text-muted);padding:6px 14px;border-radius:8px;font-size:.8rem;
    transition:all .2s;
  }}
  .btn-print:hover{{border-color:var(--accent2);color:var(--accent2)}}

  /* ---- layout ---- */
  main{{max-width:1200px;margin:0 auto;padding:2rem 1.5rem 4rem}}

  /* ---- hero ---- */
  .hero{{
    background:linear-gradient(135deg,rgba(63,185,80,.06) 0%,rgba(88,166,255,.06) 100%);
    border:1px solid var(--border);border-radius:16px;
    padding:2rem 2.5rem;margin-bottom:2rem;
  }}
  .hero h1{{font-size:1.6rem;font-weight:700;margin-bottom:.4rem}}
  .hero h1 em{{font-style:normal;color:var(--accent2)}}
  .hero .subtitle{{color:var(--text-muted);font-size:.9rem;display:flex;flex-wrap:wrap;gap:.5rem 1.5rem}}
  .hero .subtitle span{{display:flex;align-items:center;gap:.35rem}}
  .hero .subtitle .icon{{color:var(--accent)}}

  /* ---- stat cards ---- */
  .cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem;margin-bottom:2rem}}
  .card{{
    background:var(--surface);border:1px solid var(--border);border-radius:12px;
    padding:1.25rem 1.5rem;transition:border-color .2s;
  }}
  .card:hover{{border-color:var(--accent2)}}
  .card .label{{font-size:.75rem;text-transform:uppercase;letter-spacing:.06em;color:var(--text-muted);margin-bottom:.4rem}}
  .card .value{{font-size:2rem;font-weight:700;line-height:1.1}}
  .card .sub{{font-size:.75rem;color:var(--text-muted);margin-top:.3rem}}

  /* ---- sections ---- */
  section{{margin-bottom:2.5rem}}
  .section-title{{
    font-size:1.05rem;font-weight:600;margin-bottom:1rem;
    display:flex;align-items:center;gap:.5rem;
  }}
  .section-title .icon{{color:var(--accent)}}

  /* ---- tables ---- */
  .table-wrap{{overflow-x:auto;border:1px solid var(--border);border-radius:12px}}
  table{{width:100%;border-collapse:collapse;font-size:.85rem}}
  thead{{background:var(--surface2)}}
  th{{padding:.75rem 1rem;text-align:left;font-weight:600;font-size:.75rem;text-transform:uppercase;letter-spacing:.05em;color:var(--text-muted);white-space:nowrap}}
  td{{padding:.7rem 1rem;border-top:1px solid var(--border);vertical-align:middle}}
  tr:hover td{{background:rgba(255,255,255,.02)}}
  code{{font-family:'JetBrains Mono',monospace;font-size:.82em;background:var(--surface2);padding:2px 5px;border-radius:4px;color:var(--accent2)}}
  .empty-row{{text-align:center;padding:2rem!important;color:var(--text-muted)}}

  /* ---- support bar ---- */
  .support-bar-wrap{{display:flex;align-items:center;gap:.5rem}}
  .support-bar{{display:inline-block;height:8px;border-radius:4px;min-width:4px;transition:width .4s}}
  .support-num{{font-weight:600;font-size:.85rem;white-space:nowrap}}

  /* ---- badges ---- */
  .drug-badge{{display:inline-block;padding:2px 8px;border-radius:20px;font-size:.72rem;font-weight:600;white-space:nowrap}}
  .badge{{display:inline-block;padding:2px 8px;border-radius:20px;font-size:.72rem;font-weight:600}}
  .badge-clinical   {{background:rgba(34,197,94,.15); color:#22c55e;border:1px solid rgba(34,197,94,.3)}}
  .badge-agri       {{background:rgba(234,179,8,.15);  color:#eab308;border:1px solid rgba(234,179,8,.3)}}
  .badge-both       {{background:rgba(99,102,241,.15); color:#818cf8;border:1px solid rgba(99,102,241,.3)}}
  .badge-unknown    {{background:rgba(107,114,128,.15);color:#9ca3af;border:1px solid rgba(107,114,128,.3)}}

  /* ---- chart ---- */
  .chart-wrap{{
    background:var(--surface);border:1px solid var(--border);border-radius:12px;
    padding:1.5rem;position:relative;height:320px;
  }}

  /* ---- legend ---- */
  .legend{{display:flex;flex-wrap:wrap;gap:.4rem 1rem;font-size:.78rem;color:var(--text-muted);margin-bottom:1rem}}
  .legend-dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:.25rem;vertical-align:middle}}

  /* ---- params box ---- */
  .params-box{{
    background:var(--surface);border:1px solid var(--border);border-radius:12px;
    padding:1.25rem 1.5rem;font-size:.82rem;color:var(--text-muted);
    display:flex;flex-wrap:wrap;gap:.5rem 2rem;
  }}
  .params-box span{{display:flex;flex-direction:column;gap:.1rem}}
  .params-box strong{{color:var(--text);font-size:.88rem}}

  /* ---- warnings box ---- */
  .warning-box{{
    background:rgba(234,179,8,.07);border:1px solid rgba(234,179,8,.3);border-radius:12px;
    padding:1rem 1.25rem;font-size:.85rem;
  }}
  .warning-box h4{{color:#eab308;margin-bottom:.4rem;font-size:.88rem}}
  .warning-box ul{{padding-left:1.2rem;color:var(--text-muted);line-height:1.8}}

  /* ---- footer ---- */
  footer{{
    text-align:center;padding:1.5rem;border-top:1px solid var(--border);
    color:var(--text-muted);font-size:.78rem;
  }}
  footer a{{color:var(--accent2);text-decoration:none}}
  footer a:hover{{text-decoration:underline}}

  /* ---- details/summary ---- */
  details summary{{cursor:pointer;user-select:none;color:var(--accent2);font-size:.85rem;margin-bottom:.75rem}}
  details summary:hover{{text-decoration:underline}}

  /* ---- print ---- */
  @media print{{
    header{{position:relative;background:#fff;color:#000;border-bottom:1px solid #ddd}}
    body{{background:#fff;color:#000}}
    .btn-print{{display:none}}
    .card,.hero,.chart-wrap,.table-wrap,.params-box{{border:1px solid #ccc!important;background:#f9f9f9!important}}
  }}
</style>
</head>
<body>

<header>
  <div class="logo">
    <span class="mushroom">🍄</span>
    <span>FUNGAR</span>
    <span class="ver">v{version}</span>
  </div>
  <div class="header-meta">Report generated {generated_at}</div>
  <button class="btn-print" onclick="window.print()">🖨 Print / PDF</button>
</header>

<main>

  <!-- hero -->
  <div class="hero">
    <h1>🔬 Sample: <em>{sample_id}</em></h1>
    <div class="subtitle">
      <span><span class="icon">🍄</span><em>{species_display}</em></span>
      <span><span class="icon">📅</span>{generated_at}</span>
      <span><span class="icon">⏱</span>Pipeline runtime: {runtime}</span>
    </div>
  </div>

  <!-- stat cards -->
  <div class="cards">
    <div class="card">
      <div class="label">Mutations Detected</div>
      <div class="value" style="color:var(--accent)">{n_mutations}</div>
      <div class="sub">above min-support threshold</div>
    </div>
    <div class="card">
      <div class="label">Genes Affected</div>
      <div class="value" style="color:var(--accent2)">{n_genes}</div>
      <div class="sub">unique target genes</div>
    </div>
    <div class="card">
      <div class="label">Supporting Pairs</div>
      <div class="value">{total_pairs}</div>
      <div class="sub">unique sequence fragments</div>
    </div>
    <div class="card">
      <div class="label">Supporting Reads</div>
      <div class="value">{total_support}</div>
      <div class="sub">individual aligned mates</div>
    </div>
  </div>

  <!-- Summary table -->
  <section>
    <div class="section-title"><span class="icon">📊</span> Resistance Mutation Summary</div>
    <div class="legend">{legend_items}</div>
    <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>Gene</th>
            <th>Mutation</th>
            <th>Position</th>
            <th>Drug Class</th>
            <th>Fungicide(s)</th>
            <th>Source</th>
            <th>Species</th>
            <th>Read Support (Pairs / Reads)</th>
          </tr>
        </thead>
        <tbody>
          {summary_rows_html}
        </tbody>
      </table>
    </div>
  </section>

  <!-- Chart -->
  {"" if not chart_data else f'''
  <section>
    <div class="section-title"><span class="icon">📈</span> Supporting Pairs per Mutation</div>
    <div class="chart-wrap">
      <canvas id="supportChart"></canvas>
    </div>
  </section>
  '''}

  <!-- Caveats -->
  <section>
    <div class="section-title"><span class="icon">⚠️</span> Interpretation Notes</div>
    <div class="warning-box">
      <h4>Important caveats</h4>
      <ul>
        <li>Mutations are classified as <strong>computationally detected</strong> based on short-read alignment. Biological confirmation (e.g. culture-based MIC, whole-gene sequencing) is required before clinical use.</li>
        <li>Reads where the diagnostic residue falls <strong>outside the aligned region</strong> are not counted — check raw results for partial alignments.</li>
        <li>Hits from <strong>closely related species</strong> sharing conserved protein regions may generate false positives. Verify the Subject Accession in the detailed results below.</li>
        <li>Multiple ORFs and frameshifts within a single read are handled by DIAMOND; see the FUNGAR log for details.</li>
      </ul>
    </div>
  </section>

  <!-- Detailed results (collapsible) -->
  <section>
    <details>
      <summary>▶ Show raw per-read alignments (up to 200 rows)</summary>
      <div class="table-wrap">
        <table>
          <thead>
            <tr>
              <th>Read ID</th>
              <th>Gene</th>
              <th>Mutation</th>
              <th>Subject Accession</th>
              <th>Strand</th>
              <th>Source</th>
              <th>Species</th>
              <th>Alignment file</th>
            </tr>
          </thead>
          <tbody>
            {detail_rows_html}
          </tbody>
        </table>
      </div>
    </details>
  </section>

  <!-- Parameters -->
  <section>
    <div class="section-title"><span class="icon">⚙️</span> Run Parameters</div>
    <div class="params-box">
      <span><span>Sample</span><strong>{sample_id}</strong></span>
      <span><span>Species</span><strong>{species_display}</strong></span>
      <span><span>FUNGAR version</span><strong>v{version}</strong></span>
      <span><span>Runtime</span><strong>{runtime}</strong></span>
      <span><span>Parameters</span><strong>{params}</strong></span>
    </div>
  </section>

</main>

<footer>
  <p>
    Generated by <a href="https://github.com/resgen-br/fungar" target="_blank">FUNGAR v{version}</a> &mdash;
    antiFUNGAl gene Resistance detection &mdash;
    Universidade Federal do Rio Grande do Sul
  </p>
</footer>

{"" if not chart_data else f"""
<script>
const ctx = document.getElementById('supportChart').getContext('2d');
new Chart(ctx, {{
  type: 'bar',
  data: {{
    labels: {labels_js},
    datasets: [{{
      label: 'Read Support',
      data: {data_js},
      backgroundColor: {colors_js}.map(c => c + '99'),
      borderColor: {colors_js},
      borderWidth: 2,
      borderRadius: 6,
    }}]
  }},
  options: {{
    responsive: true,
    maintainAspectRatio: false,
    plugins: {{
      legend: {{ display: false }},
      tooltip: {{
        backgroundColor: '#161b22',
        borderColor: '#30363d',
        borderWidth: 1,
        titleColor: '#e6edf3',
        bodyColor: '#7d8590',
        callbacks: {{
          label: ctx => ' ' + ctx.parsed.y + ' supporting reads',
        }}
      }}
    }},
    scales: {{
      x: {{
        ticks: {{ color: '#7d8590', font: {{ size: 11, family: 'JetBrains Mono' }} }},
        grid: {{ color: '#21262d' }}
      }},
      y: {{
        beginAtZero: true,
        ticks: {{ color: '#7d8590', stepSize: 1 }},
        grid: {{ color: '#21262d' }},
        title: {{ display: true, text: 'Supporting Reads', color: '#7d8590' }}
      }}
    }}
  }}
}});
</script>
"""}

</body>
</html>"""

    return html


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate a self-contained HTML report from FUNGAR output."
    )
    parser.add_argument("--sample-id",  required=True, help="Sample identifier")
    parser.add_argument("--species",    required=True, help="Species name")
    parser.add_argument("--results",    required=True, help="Path to *_results.csv")
    parser.add_argument("--summary",    required=True, help="Path to *_summary.csv")
    parser.add_argument("--outdir",     required=True, help="Output directory")
    parser.add_argument("--version",    default="2.0", help="FUNGAR version string")
    parser.add_argument("--runtime",    default="N/A",  help="Pipeline runtime string")
    parser.add_argument("--params",     default="",     help="Run parameters string")
    args = parser.parse_args()

    # Load CSVs
    if os.path.isfile(args.results) and os.path.getsize(args.results) > 0:
        results_df = pd.read_csv(args.results)
    else:
        results_df = pd.DataFrame()

    if os.path.isfile(args.summary) and os.path.getsize(args.summary) > 0:
        summary_df = pd.read_csv(args.summary)
    else:
        summary_df = pd.DataFrame(
            columns=["Gene", "Position", "Reference", "Mutation",
                     "Fungicide", "Source", "Species", "Support_Reads", "Support_Pairs"]
        )

    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")

    html = build_html(
        sample_id=args.sample_id,
        species=args.species,
        results_df=results_df,
        summary_df=summary_df,
        version=args.version,
        runtime=args.runtime,
        params=args.params,
        generated_at=generated_at,
    )

    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, f"{args.sample_id}_report.html")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(f"✅ HTML report written to: {out_path}")


if __name__ == "__main__":
    main()
