#!/usr/bin/env python3
"""
MEGARes Cross-Database Identity Report

Compares a MEGARes FASTA against other AMR database FASTA files using
CD-HIT-EST at 100% identity to find which MEGARes genes have identical
sequences in CARD, ResFinder, PointFinder, and NCBI AMRFinderPlus.

NOTE: BacMet is intentionally excluded from this report. BacMet files are
skipped during auto-discovery and there is no --bacmet option; any BacMet
sequences that somehow appear in a reused .clstr will be reported under
the 'Other/Unknown' column rather than silently folded into a DB column.

Output: a TSV with one row per MEGARes accession and one column per
database, showing matching accession(s) from each.

Usage:
  # Using the pipeline's input directory (auto-discovers DB files):
  python megares_identity_report.py --megares new_megares.fasta --input-dir ./input

  # Specifying individual database files:
  python megares_identity_report.py --megares new_megares.fasta \\
      --card input/card/nucleotide_fasta_protein_homolog_model.fasta \\
      --card input/card/nucleotide_fasta_protein_variant_model.fasta \\
      --resfinder input/resfinder/.../all.fsa \\
      --ncbi input/ncbi/AMR_CDS.fa

  # Skip CD-HIT (reuse existing .clstr + source_map.tsv from a prior run):
  python megares_identity_report.py --parse-only \\
      --cluster-file workdir/dedup.fsa.clstr \\
      --source-map input/source_map.tsv
"""

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


# Database labels (column headers in output). BacMet deliberately omitted.
DB_COLUMNS = ["CARD", "ResFinder", "PointFinder", "NCBI_AMRFinderPlus"]

# Sources that are excluded from this report entirely.
EXCLUDED_SOURCES = {"BacMet"}


# =============================================================================
# Utility
# =============================================================================

def _seq_id(header):
    """Extract the sequence ID (first whitespace-delimited token) from a FASTA header."""
    return header.split()[0]


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def ensure_trailing_newline(filepath):
    with open(filepath, 'rb+') as f:
        f.seek(-1, 2)
        if f.read(1) != b'\n':
            f.write(b'\n')


# =============================================================================
# Source map
# =============================================================================

def build_source_map_from_files(tagged_files):
    """
    Build {seq_id: source_label} from a list of (filepath_or_list, label) tuples.
    """
    source_map = {}
    for entry, label in tagged_files:
        paths = entry if isinstance(entry, list) else [entry]
        for path in paths:
            if not path or not os.path.exists(path):
                continue
            with open(path) as f:
                for line in f:
                    if line.startswith('>'):
                        source_map[_seq_id(line[1:].strip())] = label
    return source_map


def load_source_map(path):
    """Load source map TSV (supports both Sequence_ID and Gene_Header columns)."""
    source_map = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = row.get('Sequence_ID') or row.get('Gene_Header', '')
            source_map[key] = row['Source_Database']
    return source_map


def save_source_map(source_map, path):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Sequence_ID', 'Source_Database'])
        for k, v in source_map.items():
            w.writerow([k, v])


# =============================================================================
# Auto-discover database files from pipeline input directory
# =============================================================================

def discover_db_files(input_dir):
    """
    Look inside the pipeline's input/ directory for downloaded database files.
    Returns a list of (filepath_or_list, label) tuples.

    BacMet is deliberately NOT discovered here.
    """
    tagged = []

    # CARD
    card_dir = os.path.join(input_dir, "card")
    if os.path.isdir(card_dir):
        card_files = sorted(glob.glob(os.path.join(card_dir, "nucleotide_fasta_*.fasta")))
        if card_files:
            tagged.append((card_files, "CARD"))
            print(f"    CARD:        {len(card_files)} files")

    # ResFinder
    rf_dir = os.path.join(input_dir, "resfinder")
    if os.path.isdir(rf_dir):
        rf_subdirs = glob.glob(os.path.join(rf_dir, "genomicepidemiology-resfinder_db-*"))
        rf_subdirs = [d for d in rf_subdirs if os.path.isdir(d)]
        if rf_subdirs:
            all_fsa = os.path.join(rf_subdirs[0], "all.fsa")
            if os.path.exists(all_fsa):
                tagged.append((all_fsa, "ResFinder"))
                print(f"    ResFinder:   {all_fsa}")

    # PointFinder
    pf_dir = os.path.join(input_dir, "pointfinder")
    if os.path.isdir(pf_dir):
        pf_subdirs = glob.glob(os.path.join(pf_dir, "genomicepidemiology-pointfinder_db-*"))
        pf_subdirs = [d for d in pf_subdirs if os.path.isdir(d)]
        if pf_subdirs:
            pf_fastas = sorted(glob.glob(os.path.join(pf_subdirs[0], "*", "*.fsa")))
            if pf_fastas:
                tagged.append((pf_fastas, "PointFinder"))
                print(f"    PointFinder: {len(pf_fastas)} species files")

    # NCBI
    ncbi_fa = os.path.join(input_dir, "ncbi", "AMR_CDS.fa")
    if os.path.exists(ncbi_fa):
        tagged.append((ncbi_fa, "NCBI_AMRFinderPlus"))
        print(f"    NCBI:        {ncbi_fa}")

    # BacMet: intentionally skipped
    bacmet_fa = os.path.join(input_dir, "bacmet", "BacMet2_EXP_database.fasta")
    if os.path.exists(bacmet_fa):
        print(f"    BacMet:      found but EXCLUDED from this report")

    return tagged


# =============================================================================
# Concatenation + CD-HIT
# =============================================================================

def concatenate_fastas(megares_path, tagged_files, output_path):
    """Concatenate MEGARes (first) + all other DB files into one FASTA."""
    print(f"  Concatenating (MEGARes first)...")
    with open(output_path, 'wb') as out:
        # MEGARes first
        print(f"    + {os.path.basename(megares_path)} (MEGARes)")
        with open(megares_path, 'rb') as f:
            shutil.copyfileobj(f, out)

        for entry, label in tagged_files:
            paths = entry if isinstance(entry, list) else [entry]
            for p in paths:
                if p and os.path.exists(p):
                    ensure_trailing_newline(p)
                    print(f"    + {os.path.basename(p)} ({label})")
                    with open(p, 'rb') as f:
                        shutil.copyfileobj(f, out)
    print(f"  ✓ Combined: {output_path}")


def run_cdhit_dedup(input_fasta, output_fasta, threads=6, cdhit_path="cd-hit-est"):
    """Run CD-HIT-EST at 100% identity."""
    cmd = [
        cdhit_path,
        "-i", input_fasta,
        "-o", output_fasta,
        "-G", "0",
        "-c", "1.0",
        "-AS", "0",
        "-AL", "0",
        "-aL", "1.0",
        "-aS", "1.0",
        "-d", "200",
        "-T", str(threads),
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ✗ CD-HIT failed (exit {result.returncode})")
        if result.stderr:
            print(f"    {result.stderr[:500]}")
        return False
    print(f"  ✓ Deduplication complete")
    return True


# =============================================================================
# Cluster parsing
# =============================================================================

def parse_cluster_file(filepath):
    """Parse .clstr file into list of cluster dicts."""
    clusters = []
    current = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>Cluster'):
                if current is not None:
                    clusters.append(current)
                current = {'cluster_id': int(line.split()[1]), 'sequences': []}
            else:
                seq = _parse_seq_line(line)
                if seq and current is not None:
                    current['sequences'].append(seq)
    if current is not None:
        clusters.append(current)
    return clusters


def _parse_seq_line(line):
    match = re.match(r'(\d+)\s+(\d+)nt,\s+>(.+?)\.\.\.\s*(.*)', line)
    if not match:
        return None
    gene_name = match.group(3)
    length = int(match.group(2))
    sim_part = match.group(4).strip()
    is_rep = (sim_part == '*')
    return {
        'gene_name': gene_name,
        'length': length,
        'is_representative': is_rep,
        'is_meg': gene_name.startswith('MEG_'),
    }


# =============================================================================
# Report generation
# =============================================================================

def build_identity_report(clusters, source_map):
    """
    Build per-MEGARes-accession identity report.

    Genes whose source is in EXCLUDED_SOURCES (BacMet) are dropped entirely --
    they contribute to no column and are not counted in Total_Matches. This
    matters when reusing a .clstr from a run that did include BacMet.

    Returns list of dicts, one per MEGARes gene, with columns:
      MEGARes_Accession, Length, CARD, ResFinder, PointFinder,
      NCBI_AMRFinderPlus, Other/Unknown, Total_Matches
    """
    records = []
    excluded_hits = 0

    for cluster in clusters:
        meg_genes = [s for s in cluster['sequences'] if s['is_meg']]
        non_meg_genes = [s for s in cluster['sequences'] if not s['is_meg']]

        if not meg_genes:
            continue

        # Group non-MEG genes by source database, skipping excluded sources
        by_db = defaultdict(list)
        for ext in non_meg_genes:
            source = source_map.get(_seq_id(ext['gene_name']), 'Unknown')
            if source in EXCLUDED_SOURCES:
                excluded_hits += 1
                continue
            by_db[source].append(ext['gene_name'])

        for meg in meg_genes:
            record = {
                'MEGARes_Accession': meg['gene_name'],
                'Length': meg['length'],
            }
            total = 0
            for db_col in DB_COLUMNS:
                matches = by_db.get(db_col, [])
                record[db_col] = ' ; '.join(matches) if matches else ''
                total += len(matches)

            # Capture any Unknown-source matches too
            unknown = by_db.get('Unknown', [])
            record['Other/Unknown'] = ' ; '.join(unknown) if unknown else ''
            total += len(unknown)

            record['Total_Matches'] = total
            records.append(record)

    if excluded_hits:
        print(f"  Excluded {excluded_hits} match(es) from: {', '.join(sorted(EXCLUDED_SOURCES))}")

    return records


def write_report(records, output_path):
    """Write identity report to TSV."""
    columns = [
        'MEGARes_Accession', 'Length',
        *DB_COLUMNS,
        'Other/Unknown', 'Total_Matches',
    ]
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def print_summary(records):
    total_meg = len(records)
    with_any = sum(1 for r in records if r['Total_Matches'] > 0)
    unique = total_meg - with_any
    print(f"  MEGARes genes:               {total_meg}")
    print(f"    With ≥1 identical match:    {with_any}")
    print(f"    Unique to MEGARes:          {unique}")
    print(f"  Per-database match counts:")
    for db in DB_COLUMNS:
        count = sum(1 for r in records if r[db])
        print(f"    {db:25s} {count}")
    unk_count = sum(1 for r in records if r['Other/Unknown'])
    if unk_count:
        print(f"    {'Other/Unknown':25s} {unk_count}")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compare a MEGARes FASTA against other AMR databases at 100%% identity "
                    "(BacMet excluded)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-discover databases from pipeline's input/ directory:
  %(prog)s --megares new_megares_v4.fasta --input-dir ./input

  # Specify individual database files:
  %(prog)s --megares megares.fasta \\
      --card card/nucleotide_fasta_protein_homolog_model.fasta \\
      --card card/nucleotide_fasta_protein_variant_model.fasta \\
      --resfinder resfinder/all.fsa \\
      --ncbi ncbi/AMR_CDS.fa

  # Reuse existing dedup results (skip CD-HIT):
  %(prog)s --parse-only \\
      --cluster-file workdir/dedup.fsa.clstr \\
      --source-map input/source_map.tsv

Note:
  BacMet is excluded. It is skipped during --input-dir auto-discovery, and
  BacMet-sourced hits in a reused .clstr are dropped from the report.

Output:
  TSV with one row per MEGARes accession, columns for each database
  showing 100%%-identical accessions from that database.
        """
    )

    parser.add_argument("--megares", type=str,
                        help="Path to MEGARes FASTA file")
    parser.add_argument("--input-dir", type=str,
                        help="Pipeline input/ directory (auto-discovers DB files)")
    parser.add_argument("--card", action="append", default=[],
                        help="CARD FASTA file(s) (can specify multiple)")
    parser.add_argument("--resfinder", action="append", default=[],
                        help="ResFinder FASTA file(s)")
    parser.add_argument("--pointfinder", action="append", default=[],
                        help="PointFinder FASTA file(s)")
    parser.add_argument("--ncbi", action="append", default=[],
                        help="NCBI AMRFinderPlus FASTA file(s)")

    parser.add_argument("--output", "-o", default="megares_identity_report.tsv",
                        help="Output TSV path (default: megares_identity_report.tsv)")
    parser.add_argument("--workdir", "-w", default="./identity_workdir",
                        help="Working directory for intermediate files (default: ./identity_workdir)")
    parser.add_argument("--cdhit-path", default="cd-hit-est",
                        help="Path to cd-hit-est (default: cd-hit-est)")
    parser.add_argument("--threads", "-T", type=int, default=6,
                        help="Threads for CD-HIT (default: 6)")

    # Parse-only mode
    parser.add_argument("--parse-only", action="store_true",
                        help="Skip CD-HIT; parse existing .clstr + source_map")
    parser.add_argument("--cluster-file", type=str,
                        help="Path to .clstr file (for --parse-only)")
    parser.add_argument("--source-map", type=str,
                        help="Path to source_map.tsv (for --parse-only)")

    args = parser.parse_args()

    print("=" * 60)
    print("MEGARes Cross-Database Identity Report")
    print("=" * 60)
    print(f"  Databases compared: {', '.join(DB_COLUMNS)}")
    print(f"  Excluded:           {', '.join(sorted(EXCLUDED_SOURCES))}")

    # ------------------------------------------------------------------
    # Parse-only mode
    # ------------------------------------------------------------------
    if args.parse_only:
        if not args.cluster_file:
            print("Error: --parse-only requires --cluster-file")
            sys.exit(1)
        if not args.source_map:
            print("Error: --parse-only requires --source-map")
            sys.exit(1)

        print(f"\n  Cluster file: {args.cluster_file}")
        print(f"  Source map:   {args.source_map}")

        source_map = load_source_map(args.source_map)
        print(f"  Loaded {len(source_map)} source map entries")

        clusters = parse_cluster_file(args.cluster_file)
        print(f"  Parsed {len(clusters)} clusters")

        records = build_identity_report(clusters, source_map)
        print_summary(records)
        write_report(records, args.output)
        print(f"\n  Output: {args.output}")
        return

    # ------------------------------------------------------------------
    # Full run: need MEGARes + at least one database
    # ------------------------------------------------------------------
    if not args.megares:
        print("Error: --megares is required (path to MEGARes FASTA)")
        sys.exit(1)
    if not os.path.exists(args.megares):
        print(f"Error: MEGARes file not found: {args.megares}")
        sys.exit(1)

    ensure_dir(args.workdir)

    # Collect database files
    tagged_files = []

    if args.input_dir:
        print(f"\n  Auto-discovering databases in: {args.input_dir}")
        tagged_files = discover_db_files(args.input_dir)

    # CLI-specified files override / supplement auto-discovery
    if args.card:
        tagged_files.append((args.card, "CARD"))
    if args.resfinder:
        tagged_files.append((args.resfinder, "ResFinder"))
    if args.pointfinder:
        tagged_files.append((args.pointfinder, "PointFinder"))
    if args.ncbi:
        tagged_files.append((args.ncbi, "NCBI_AMRFinderPlus"))

    if not tagged_files:
        print("Error: no database files found. Use --input-dir or specify with --card, --resfinder, etc.")
        sys.exit(1)

    # Build source map
    print("\n  Building source map...")
    # Include MEGARes in source map for completeness
    all_tagged = [(args.megares, "MEGARes")] + tagged_files
    source_map = build_source_map_from_files(all_tagged)
    source_map_path = os.path.join(args.workdir, "source_map.tsv")
    save_source_map(source_map, source_map_path)
    print(f"  Source map: {len(source_map)} entries → {source_map_path}")

    # Concatenate (MEGARes first)
    print()
    combined = os.path.join(args.workdir, "combined.fsa")
    ensure_trailing_newline(args.megares)
    concatenate_fastas(args.megares, tagged_files, combined)

    # Run CD-HIT
    print()
    dedup_out = os.path.join(args.workdir, "dedup.fsa")
    if not run_cdhit_dedup(combined, dedup_out,
                           threads=args.threads, cdhit_path=args.cdhit_path):
        sys.exit(1)

    clstr_path = dedup_out + ".clstr"

    # Parse clusters and build report
    print()
    clusters = parse_cluster_file(clstr_path)
    print(f"  Parsed {len(clusters)} clusters")

    records = build_identity_report(clusters, source_map)
    print_summary(records)

    write_report(records, args.output)
    print(f"\n  Output: {args.output}")
    print(f"  Working files: {args.workdir}/")


if __name__ == '__main__':
    main()
