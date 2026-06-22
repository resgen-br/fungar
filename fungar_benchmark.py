#!/usr/bin/env python3
"""
fungar_benchmark.py — Quantitative benchmarking for FUNGAR.

Generates synthetic reads from protein database sequences with known mutations
planted at defined coverage depths, runs FUNGAR, and reports sensitivity,
specificity, false-positive rate, precision, F1, runtime, and memory.

Usage:
    python3 fungar_benchmark.py \\
        -sp Aspergillus_fumigatus \\
        -d  database/ \\
        -o  benchmark_out/ \\
        --fungar FUNGAR.sh \\
        [--depths 5 10 20 50] \\
        [--read-length 150] \\
        [--threads 4] \\
        [--min-support 1] \\
        [--seed 42]
"""

import argparse
import csv
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path


def fmt_depth(d: float) -> str:
    """Format a depth nicely, e.g. 5.0 -> 5, 0.5 -> 0.5."""
    return f"{int(d)}" if float(d).is_integer() else f"{d}"

try:
    import pandas as pd
except ImportError:
    sys.exit("❌ pandas is required: conda install pandas")

# ---------------------------------------------------------------------------
# Standard genetic code (table 1) — preferred codon per amino acid
# ---------------------------------------------------------------------------
_CODONS = {
    'A': 'GCC', 'R': 'CGG', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
    'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTC', 'P': 'CCC',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTG',
    '*': 'TGA', '-': 'NNN',
}


def rev_translate(protein: str) -> str:
    """Reverse-translate a protein string to a DNA string."""
    return "".join(_CODONS.get(aa, 'NNN') for aa in protein.upper())


def complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


# ---------------------------------------------------------------------------
# FASTA reader
# ---------------------------------------------------------------------------
def read_fasta(path: str) -> dict:
    seqs = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


# ---------------------------------------------------------------------------
# Synthetic read generation
# ---------------------------------------------------------------------------
def make_reads(protein: str, mut_pos_0: int, mut_aa: str,
               read_length: int, depth: float, rng: random.Random,
               plant_mutation: bool = True, variable_lengths: bool = False,
               min_read_length: int = 20) -> list:
    """
    Return a list of (read_id, dna_seq, qual) tuples.

    Strategy:
    - Introduce (or keep) the amino acid at mut_pos_0 (0-indexed).
    - Reverse-translate the protein to DNA.
    - Sample the actual read count covering the mutation stochastically from a
      Poisson distribution with mean equal to 'depth'.
    - If variable_lengths is True, sample read lengths dynamically using a mixture:
      70% full length, 30% uniform between min_read_length and full length.
    - Slide a window covering mut_pos_0 across the region.
    """
    prot = list(protein)
    if plant_mutation:
        prot[mut_pos_0] = mut_aa
    dna = rev_translate("".join(prot))

    mut_nt_centre = mut_pos_0 * 3 + 1   # centre nucleotide of mutated codon
    reads = []

    # Sample read count stochastically from a Poisson distribution of mean depth
    if depth > 30.0:
        # Gaussian approximation to prevent float exponent underflow
        u1 = rng.random()
        u2 = rng.random()
        z = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
        num_reads = max(0, int(round(depth + z * math.sqrt(depth))))
    else:
        # Knuth's algorithm
        L = math.exp(-depth)
        k = 0
        p = 1.0
        while p > L:
            k += 1
            p *= rng.random()
        num_reads = k - 1

    for i in range(num_reads):
        if variable_lengths:
            # 70% probability of full read_length, 30% uniform down to min_read_length
            if rng.random() < 0.7:
                curr_read_length = read_length
            else:
                curr_read_length = rng.randint(min_read_length, read_length)
        else:
            curr_read_length = read_length

        qual = 'I' * curr_read_length             # Phred 40

        # jitter the window start so reads don't all start at the same position
        half = curr_read_length // 2
        jitter = rng.randint(-10, 10)
        start = max(0, mut_nt_centre - half + jitter)
        end = start + curr_read_length
        if end > len(dna):
            start = max(0, len(dna) - curr_read_length)
            end = len(dna)
        fragment = dna[start:end]
        if len(fragment) < min_read_length:          # too short to be useful
            continue
        # pad quality to fragment length
        q = qual[:len(fragment)]
        strand = rng.choice(['fwd', 'rev'])
        if strand == 'rev':
            fragment = reverse_complement(fragment)
        reads.append((f"synth_{i}_{strand}", fragment, q))

    return reads


def write_fastq(reads: list, path: str) -> None:
    with open(path, 'w') as fh:
        for rid, seq, qual in reads:
            fh.write(f"@{rid}\n{seq}\n+\n{qual}\n")


# ---------------------------------------------------------------------------
# Run FUNGAR
# ---------------------------------------------------------------------------
def run_fungar(fungar_script: str, sample_id: str, fastq: str,
               species: str, db_dir: str, outdir: str,
               threads: int, min_support: int, min_orf: int, min_pident: float) -> tuple:
    """
    Returns (summary_df, wall_s, peak_rss_kb).
    """
    cmd = [
        'bash', fungar_script,
        '-id', sample_id,
        '-sp', species,
        '-s',  fastq,
        '-d',  db_dir,
        '-o',  outdir,
        '-c',  str(threads),
        '--min-support', str(min_support),
        '-orf', str(min_orf),
        '--min-pident', str(min_pident),
    ]
    t0 = time.perf_counter()
    my_env = os.environ.copy()
    my_env["PATH"] = os.path.dirname(sys.executable) + os.path.pathsep + my_env.get("PATH", "")
    result = subprocess.run(
        cmd, capture_output=True, text=True, env=my_env
    )
    wall_s = time.perf_counter() - t0

    # Try to read peak RSS from runtime file
    peak_rss_kb = None
    runtime_file = os.path.join(outdir, sample_id, f"{sample_id}_runtime.txt")
    if os.path.isfile(runtime_file):
        with open(runtime_file) as fh:
            for line in fh:
                m = re.search(r'Maximum resident set size.*?(\d+)', line)
                if m:
                    peak_rss_kb = int(m.group(1))
                    break

    summary_path = os.path.join(outdir, sample_id, f"{sample_id}_summary.csv")
    if os.path.isfile(summary_path) and os.path.getsize(summary_path) > 0:
        try:
            summary_df = pd.read_csv(summary_path)
        except Exception:
            summary_df = pd.DataFrame()
    else:
        summary_df = pd.DataFrame()

    if result.returncode != 0:
        print(f"    ⚠ FUNGAR exited with code {result.returncode}")
        print(result.stderr[-800:] if result.stderr else "")

    return summary_df, wall_s, peak_rss_kb


# ---------------------------------------------------------------------------
# Evaluate one run
# ---------------------------------------------------------------------------
def evaluate(planted: list, detected_df) -> dict:
    """
    planted: list of (gene, position, mutation) tuples that were planted.
    detected_df: summary DataFrame from FUNGAR.
    Returns dict with TP, FP, FN, sensitivity, precision, f1.
    """
    if detected_df.empty:
        detected_set = set()
    else:
        detected_set = set(
            zip(detected_df['Gene'].astype(str),
                detected_df['Position'].astype(int),
                detected_df['Mutation'].astype(str))
        )
    planted_set = set(planted)

    TP = len(planted_set & detected_set)
    FN = len(planted_set - detected_set)
    FP = len(detected_set - planted_set)
    TN = 0   # undefined in mutation-only mode; use WT run for specificity

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else float('nan')
    precision   = TP / (TP + FP) if (TP + FP) > 0 else float('nan')
    f1 = (2 * sensitivity * precision / (sensitivity + precision)
          if (sensitivity + precision) > 0 else float('nan'))

    return dict(TP=TP, FP=FP, FN=FN,
                sensitivity=sensitivity, precision=precision, f1=f1)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
def _fmt(v, decimals=3):
    if v != v:   # nan
        return "—"
    return f"{v:.{decimals}f}"


def build_benchmark_html(rows: list, species: str, depths: list,
                         version: str, generated_at: str) -> str:
    table_rows = ""
    for r in rows:
        sens_color = ("#22c55e" if (r['sensitivity'] == r['sensitivity'] and r['sensitivity'] >= 0.8)
                      else "#ef4444" if r['sensitivity'] == r['sensitivity'] else "#6b7280")
        fpr_color  = ("#22c55e" if (r['fpr'] == r['fpr'] and r['fpr'] <= 0.1)
                      else "#ef4444" if r['fpr'] == r['fpr'] else "#6b7280")
        d_fmt = fmt_depth(r['depth'])
        table_rows += f"""
        <tr>
          <td>{r['gene_short']}</td>
          <td><code>{r['notation']}</code></td>
          <td>{d_fmt}×</td>
          <td style="color:{sens_color};font-weight:600">{_fmt(r['sensitivity'])}</td>
          <td>{_fmt(r['precision'])}</td>
          <td>{_fmt(r['f1'])}</td>
          <td style="color:{fpr_color};font-weight:600">{_fmt(r['fpr'])}</td>
          <td>{r['TP']}</td><td>{r['FP']}</td><td>{r['FN']}</td>
          <td>{r['wall_s']:.1f}s</td>
          <td>{'N/A' if r['peak_rss_kb'] is None else str(r['peak_rss_kb'])}</td>
        </tr>"""

    # aggregate across all depths
    if rows:
        all_sens = [r['sensitivity'] for r in rows if r['sensitivity'] == r['sensitivity']]
        avg_sens = sum(all_sens) / len(all_sens) if all_sens else float('nan')
        all_fpr  = [r['fpr'] for r in rows if r['fpr'] == r['fpr']]
        avg_fpr  = sum(all_fpr) / len(all_fpr) if all_fpr else float('nan')
        avg_wall = sum(r['wall_s'] for r in rows) / len(rows)

        # Calculate optimum coverage: minimum depth where avg sensitivity >= 95% and avg FDR <= 5%
        optimum_depth = None
        depths_sorted = sorted(list(set(r['depth'] for r in rows)))
        for d in depths_sorted:
            d_rows = [r for r in rows if r['depth'] == d]
            s_vals = [r['sensitivity'] for r in d_rows if r['sensitivity'] == r['sensitivity']]
            avg_sens_d = sum(s_vals) / len(s_vals) if s_vals else 0.0

            fdr_vals = []
            for r in d_rows:
                tp_fp = r['TP'] + r['FP']
                fdr = r['FP'] / tp_fp if tp_fp > 0 else 0.0
                fdr_vals.append(fdr)
            avg_fdr_d = sum(fdr_vals) / len(fdr_vals) if fdr_vals else 0.0

            if avg_sens_d >= 0.90 and avg_fdr_d <= 0.10:
                optimum_depth = d
                break
        opt_depth_str = f"{fmt_depth(optimum_depth)}×" if optimum_depth is not None else "N/A"
    else:
        avg_sens = avg_fpr = avg_wall = float('nan')
        opt_depth_str = "N/A"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>FUNGAR Benchmark — {species.replace('_',' ')}</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet"/>
<style>
:root{{--bg:#0d1117;--surf:#161b22;--surf2:#21262d;--bdr:#30363d;--txt:#e6edf3;--muted:#7d8590;--green:#22c55e;--red:#ef4444;--blue:#58a6ff;}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'Inter',sans-serif;background:var(--bg);color:var(--txt);min-height:100vh;}}
header{{background:rgba(13,17,23,.9);backdrop-filter:blur(10px);border-bottom:1px solid var(--bdr);padding:.9rem 2rem;display:flex;align-items:center;gap:.8rem;position:sticky;top:0;z-index:10}}
header h1{{font-size:1.1rem;font-weight:700}} header .ver{{font-size:.7rem;background:rgba(63,185,80,.15);color:var(--green);padding:2px 8px;border-radius:20px;border:1px solid rgba(63,185,80,.3)}}
main{{max-width:1300px;margin:0 auto;padding:2rem 1.5rem 4rem}}
.hero{{background:linear-gradient(135deg,rgba(63,185,80,.05),rgba(88,166,255,.05));border:1px solid var(--bdr);border-radius:14px;padding:1.8rem 2rem;margin-bottom:2rem}}
.hero h2{{font-size:1.4rem;margin-bottom:.4rem}} .hero p{{color:var(--muted);font-size:.88rem}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem;margin-bottom:2rem}}
.card{{background:var(--surf);border:1px solid var(--bdr);border-radius:12px;padding:1.1rem 1.3rem}}
.card .lbl{{font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);margin-bottom:.3rem}}
.card .val{{font-size:1.9rem;font-weight:700;line-height:1.1}}
.section-title{{font-size:1rem;font-weight:600;margin-bottom:.9rem;display:flex;align-items:center;gap:.4rem}}
.section-title .ic{{color:var(--green)}}
.table-wrap{{overflow-x:auto;border:1px solid var(--bdr);border-radius:12px;margin-bottom:2rem}}
table{{width:100%;border-collapse:collapse;font-size:.83rem}}
thead{{background:var(--surf2)}}
th{{padding:.65rem .9rem;text-align:left;font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);white-space:nowrap}}
td{{padding:.65rem .9rem;border-top:1px solid var(--bdr);vertical-align:middle}}
tr:hover td{{background:rgba(255,255,255,.02)}}
code{{font-family:'JetBrains Mono',monospace;font-size:.8em;background:var(--surf2);padding:2px 5px;border-radius:4px;color:var(--blue)}}
.warn{{background:rgba(234,179,8,.07);border:1px solid rgba(234,179,8,.28);border-radius:10px;padding:.9rem 1.1rem;font-size:.83rem;color:var(--muted);margin-bottom:2rem}}
.warn strong{{color:#eab308}}
footer{{text-align:center;padding:1.2rem;border-top:1px solid var(--bdr);color:var(--muted);font-size:.75rem}}
footer a{{color:var(--blue);text-decoration:none}}
</style>
</head>
<body>
<header>
  <span style="font-size:1.4rem">🍄</span>
  <h1>FUNGAR Benchmark</h1>
  <span class="ver">v{version}</span>
  <span style="margin-left:auto;font-size:.78rem;color:var(--muted)">{generated_at}</span>
</header>
<main>
<div class="hero">
  <h2>Species: <em style="font-style:normal;color:var(--blue)">{species.replace('_',' ')}</em></h2>
  <p>Synthetic read benchmark across depths: {', '.join(fmt_depth(d)+'×' for d in depths)}.
  Reads generated by reverse-translating database protein sequences with mutations planted at target positions.</p>
</div>
<div class="cards">
  <div class="card"><div class="lbl">Avg Sensitivity</div><div class="val" style="color:var(--green)">{_fmt(avg_sens)}</div></div>
  <div class="card"><div class="lbl">Avg FPR</div><div class="val" style="color:{'var(--green)' if avg_fpr==avg_fpr and avg_fpr<=0.1 else 'var(--red)'}">{_fmt(avg_fpr)}</div></div>
  <div class="card">
    <div class="lbl">Optimum Coverage</div>
    <div class="val" style="color:var(--blue)">{opt_depth_str}</div>
    <div class="lbl" style="font-size:.62rem;margin-top:0.25rem;text-transform:none;color:var(--muted);line-height:1.2">(Sens &ge; 90%, FDR &le; 10%)</div>
  </div>
  <div class="card"><div class="lbl">Avg Runtime</div><div class="val">{avg_wall:.1f}s</div></div>
  <div class="card"><div class="lbl">Mutations Tested</div><div class="val">{len(set((r['gene'],r['pos'],r['mut']) for r in rows))}</div></div>
</div>
<div class="warn">
  <strong>⚠ Note:</strong>
  Benchmarks use synthetic reads produced by reverse-translating database sequences.
  Real-world sensitivity may differ due to sequencing error, read length variation, indels,
  frameshifts, and homologous species cross-mapping.
  FPR is computed from wild-type (no mutation planted) runs.
</div>
<div class="section-title"><span class="ic">📊</span> Per-Mutation Metrics</div>
<div class="table-wrap">
<table>
  <thead>
    <tr>
      <th>Gene</th><th>Mutation</th><th>Depth</th>
      <th>Sensitivity</th><th>Precision</th><th>F1</th><th>FPR</th>
      <th>TP</th><th>FP</th><th>FN</th>
      <th>Wall Time</th><th>Peak RSS (KB)</th>
    </tr>
  </thead>
  <tbody>{table_rows}</tbody>
</table>
</div>
</main>
<footer>
  Generated by <a href="https://github.com/resgen-br/fungar">FUNGAR v{version}</a> benchmark tool &mdash; {generated_at}
</footer>
</body>
</html>"""
    return html


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="FUNGAR quantitative benchmark — synthetic read simulation."
    )
    parser.add_argument('-sp', '--species',  required=True)
    parser.add_argument('-d',  '--db-dir',   required=True)
    parser.add_argument('-o',  '--outdir',   required=True)
    parser.add_argument('--fungar',          default=os.path.join(os.path.dirname(__file__), 'FUNGAR.sh'),
                        help="Path to FUNGAR.sh (default: same directory as this script)")
    parser.add_argument('--depths',  nargs='+', type=float, default=[0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0],
                        metavar='D', help="Coverage depths to simulate (default: 0.1 0.5 1.0 5.0 10.0 20.0 50.0)")
    parser.add_argument('--read-length', type=int, default=150)
    parser.add_argument('--no-variable-read-lengths', dest='variable_read_lengths', action='store_false',
                        help="Disable simulating variable read lengths (i.e. all reads will be full length)")
    parser.set_defaults(variable_read_lengths=True)
    parser.add_argument('--min-read-length', type=int, default=20,
                        help="Minimum read length to simulate when using variable read lengths (default: 20)")
    parser.add_argument('--threads',     type=int, default=4)
    parser.add_argument('--min-support', type=int, default=1)
    parser.add_argument('--min-orf',     type=int, default=50)
    parser.add_argument('--min-pident',  type=float, default=0.0)
    parser.add_argument('--max-mutations', type=int, default=None,
                        help="Maximum number of mutations to test (randomly sampled)")
    parser.add_argument('--seed',        type=int, default=42)
    parser.add_argument('--version',     default='2.0')
    args = parser.parse_args()

    rng = random.Random(args.seed)

    # ---- locate DB files ----
    faa_path = os.path.join(args.db_dir, f"{args.species}_db.faa")
    csv_path = os.path.join(args.db_dir, f"{args.species}_gene_mutations.csv")
    for p in [faa_path, csv_path, args.fungar]:
        if not os.path.isfile(p):
            sys.exit(f"❌ File not found: {p}")

    proteins = read_fasta(faa_path)
    mut_df = (
    pd.read_csv(csv_path)
    .drop_duplicates(subset=['gene', 'position', 'reference', 'mutation']))

    if args.max_mutations and len(mut_df) > args.max_mutations:
        print(f"   Downsampling mutations from {len(mut_df)} to {args.max_mutations}...")
        mut_df = mut_df.sample(n=args.max_mutations, random_state=args.seed)

    os.makedirs(args.outdir, exist_ok=True)
    fastq_dir = os.path.join(args.outdir, '_fastq')
    os.makedirs(fastq_dir, exist_ok=True)

    print(f"🍄 FUNGAR Benchmark — {args.species}")
    depth_strs = [fmt_depth(d) for d in args.depths]
    len_str = f"varying ({args.min_read_length}-{args.read_length} bp)" if args.variable_read_lengths else f"{args.read_length} bp"
    print(f"   Depths: {depth_strs}  |  Read length: {len_str}  |  Min support: {args.min_support}")
    print(f"   Mutations in CSV: {len(mut_df)}")
    print()

    rows = []

    for _, mrow in mut_df.iterrows():
        gene     = str(mrow['gene'])
        pos_1    = int(mrow['position'])   # 1-indexed
        ref_aa   = str(mrow['reference'])
        mut_aa   = str(mrow['mutation'])

        if gene not in proteins:
            print(f"   ⚠ Gene {gene} not in FASTA — skipping")
            continue

        prot   = proteins[gene]
        pos_0  = pos_1 - 1                # convert to 0-indexed

        if pos_0 >= len(prot):
            print(f"   ⚠ Position {pos_1} out of range for {gene} (len={len(prot)}) — skipping")
            continue

        if prot[pos_0] != ref_aa:
            print(f"   ℹ {gene} pos {pos_1}: expected ref {ref_aa}, found {prot[pos_0]} — proceeding anyway")

        gene_parts = gene.split("_")
        gene_short = "_".join(gene_parts[:-1]) if len(gene_parts) > 1 else gene
        notation   = f"p.{ref_aa}{pos_1}{mut_aa}"

        # skip gap/deletion mutations (we can't easily synthesise them)
        if mut_aa in ('-', '*'):
            print(f"   ⚠ Skipping deletion/stop mutation {gene} {notation}")
            continue

        for depth in args.depths:
            d_formatted = fmt_depth(depth)
            safe_id = re.sub(r'[^A-Za-z0-9_]', '_', f"{gene_short}_{pos_1}{mut_aa}_d{d_formatted}")
            print(f"   ▶ {gene_short} {notation}  depth={d_formatted}×", end="  ", flush=True)

            # --- positive run (mutation planted) ---
            pos_reads = make_reads(
                prot, pos_0, mut_aa, args.read_length, depth, rng,
                plant_mutation=True, variable_lengths=args.variable_read_lengths,
                min_read_length=args.min_read_length
            )
            pos_fq    = os.path.join(fastq_dir, f"{safe_id}_MUT.fastq")
            write_fastq(pos_reads, pos_fq)

            pos_outdir = os.path.join(args.outdir, 'runs', f"{safe_id}_MUT")
            pos_df, wall_mut, rss_mut = run_fungar(
                args.fungar, safe_id, pos_fq, args.species, args.db_dir,
                pos_outdir, args.threads, args.min_support, args.min_orf, args.min_pident
            )
            planted = [(gene, pos_1, mut_aa)]
            metrics_pos = evaluate(planted, pos_df)

            # --- wild-type run (no mutation) to measure FPR ---
            wt_reads  = make_reads(
                prot, pos_0, mut_aa, args.read_length, depth, rng,
                plant_mutation=False, variable_lengths=args.variable_read_lengths,
                min_read_length=args.min_read_length
            )
            wt_fq     = os.path.join(fastq_dir, f"{safe_id}_WT.fastq")
            write_fastq(wt_reads, wt_fq)

            wt_outdir = os.path.join(args.outdir, 'runs', f"{safe_id}_WT")
            wt_df, wall_wt, rss_wt = run_fungar(
                args.fungar, f"{safe_id}_wt", wt_fq, args.species, args.db_dir,
                wt_outdir, args.threads, args.min_support, args.min_orf, args.min_pident
            )
            # FPR from WT run: proportion of mutations flagged when none planted
            n_wt_fp = 0 if wt_df.empty else len(wt_df)
            n_possible = len(mut_df)
            fpr = n_wt_fp / n_possible if n_possible > 0 else float('nan')

            total_wall = wall_mut + wall_wt
            peak_rss   = rss_mut or rss_wt   # whichever is available

            row = dict(
                gene=gene, gene_short=gene_short, pos=pos_1, mut=mut_aa,
                notation=notation, depth=depth,
                wall_s=total_wall, peak_rss_kb=peak_rss, fpr=fpr,
                **metrics_pos
            )
            rows.append(row)

            status = ("✅" if metrics_pos['sensitivity'] >= 0.8 else
                      "⚠" if metrics_pos['sensitivity'] > 0 else "❌")
            print(f"Sens={_fmt(metrics_pos['sensitivity'],2)}  FPR={_fmt(fpr,2)}  {status}")

    # ---- write CSVs ----
    csv_path_out = os.path.join(args.outdir, f"benchmark_{args.species}.csv")
    if rows:
        fieldnames = [
            'gene_short', 'notation', 'gene', 'pos', 'mut', 'depth',
            'TP', 'FP', 'FN', 'sensitivity', 'precision', 'f1', 'fpr',
            'wall_s', 'peak_rss_kb'
        ]
        with open(csv_path_out, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n📄 Results CSV: {csv_path_out}")

    # ---- write HTML ----
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")
    html = build_benchmark_html(rows, args.species, args.depths, args.version, generated_at)
    html_path = os.path.join(args.outdir, f"benchmark_{args.species}_report.html")
    with open(html_path, 'w', encoding='utf-8') as fh:
        fh.write(html)
    print(f"📊 HTML report: {html_path}")

    # ---- summary ----
    if rows:
        sens_vals = [r['sensitivity'] for r in rows if r['sensitivity'] == r['sensitivity']]
        fpr_vals  = [r['fpr'] for r in rows if r['fpr'] == r['fpr']]

        # Calculate optimum coverage: minimum depth where avg sensitivity >= 95% and avg FDR <= 5%
        optimum_depth = None
        depths_sorted = sorted(list(set(r['depth'] for r in rows)))
        for d in depths_sorted:
            d_rows = [r for r in rows if r['depth'] == d]
            s_vals = [r['sensitivity'] for r in d_rows if r['sensitivity'] == r['sensitivity']]
            avg_sens_d = sum(s_vals) / len(s_vals) if s_vals else 0.0

            fdr_vals = []
            for r in d_rows:
                tp_fp = r['TP'] + r['FP']
                fdr = r['FP'] / tp_fp if tp_fp > 0 else 0.0
                fdr_vals.append(fdr)
            avg_fdr_d = sum(fdr_vals) / len(fdr_vals) if fdr_vals else 0.0

            if avg_sens_d >= 0.90 and avg_fdr_d <= 0.10:
                optimum_depth = d
                break
        opt_depth_str = f"{fmt_depth(optimum_depth)}×" if optimum_depth is not None else "Not reached under tested depths"

        print("\n========================================")
        print(f"  Average sensitivity : {sum(sens_vals)/len(sens_vals):.3f}" if sens_vals else "  No sensitivity data")
        print(f"  Average FPR         : {sum(fpr_vals)/len(fpr_vals):.3f}"  if fpr_vals  else "  No FPR data")
        print(f"  Optimum Coverage    : {opt_depth_str}")
        print(f"  Total mutations tested (unique): {len(set((r['gene'],r['pos'],r['mut']) for r in rows))}")
        print("========================================")
    else:
        print("\n⚠ No mutations could be benchmarked.")

    print("\n✅ Benchmark complete.")


def _fmt(v, decimals=3):
    if v != v:
        return "—"
    return f"{v:.{decimals}f}"


if __name__ == '__main__':
    main()
