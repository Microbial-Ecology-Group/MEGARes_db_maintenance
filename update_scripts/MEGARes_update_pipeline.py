#!/usr/bin/env python3
"""
AMR Database Aggregation, Clustering, and Annotation Pipeline

Downloads AMR gene sequences from multiple public databases, concatenates them
(MEGARes first to guarantee representative status), deduplicates, clusters,
and parses clusters to inherit MEGARes annotations for novel genes.

Workflow:
  Step 1:  Download all databases → concatenate into input/all.fsa
           Also builds input/source_map.tsv (header → source database)
  Step 2:  Deduplicate at 100% identity (cd-hit-est) → input/dedup.fsa
  Step 2b: (Optional) Parse dedup clusters → output/megares_dedup_report.tsv
           Shows which external genes are 100% identical to each MEGARes gene
  Step 3:  Cluster at 80% similarity (cd-hit-est) → output/megares-update-clustering.fsa
  Step 4:  Parse clusters for annotation inheritance → output/..._parsed.tsv

Databases (all nucleotide sequences):
  - MEGARes v3.00
  - CARD v4.0.1           (nucleotide_fasta_*.fasta files)
  - ResFinder 4.7.2       (nucleotide .fsa files)
  - PointFinder 4.7.1     (nucleotide .fsa files)
  - NCBI AMRFinderPlus 4.2 (AMR_CDS.fa — coding DNA sequences)

Note: BacMet is excluded because it only provides protein sequences
(BacMet2_EXP_database.fasta). cd-hit-est requires nucleotide input.
BacMet-derived genes already in MEGARes were incorporated via manual
UniProt cross-referencing by the MEGARes team (see MEGARes 2.0 paper).
"""

import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path
from urllib.request import urlretrieve, urlopen, Request
from urllib.error import URLError


# =============================================================================
# Database Version Configuration
# =============================================================================

MEGARES_VERSION = "v3.00"
CARD_VERSION = "v4.0.1"
RESFINDER_VERSION = "4.7.2"
POINTFINDER_VERSION = "4.7.1"
NCBI_VERSION = "4.2/2026-03-24.1"

# Source database labels used in source_map and dedup report
SOURCE_MEGARES = "MEGARes"
SOURCE_CARD = "CARD"
SOURCE_RESFINDER = "ResFinder"
SOURCE_POINTFINDER = "PointFinder"
SOURCE_NCBI = "NCBI_AMRFinderPlus"

POINTFINDER_SPECIES = [
    "campylobacter",
    "enterococcus_faecalis",
    "enterococcus_faecium",
    "escherichia_coli",
    "helicobacter_pylori",
    "klebsiella",
    "mycobacterium_tuberculosis",
    "neisseria_gonorrhoeae",
    "plasmodium_falciparum",
    "salmonella",
    "staphylococcus_aureus",
]


# =============================================================================
# Utility Functions
# =============================================================================

def ensure_dir(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def ensure_trailing_newline(filepath):
    """Ensure file ends with a newline character."""
    with open(filepath, 'rb+') as f:
        f.seek(-1, 2)
        if f.read(1) != b'\n':
            f.write(b'\n')


def download_file(url, dest_path, description="file"):
    """Download a file with progress indication."""
    print(f"  Downloading {description}...")
    print(f"    URL: {url}")
    ensure_dir(os.path.dirname(dest_path))
    try:
        urlretrieve(url, dest_path)
        print(f"    ✓ Downloaded")
        return True
    except URLError as e:
        print(f"    ✗ Failed: {e}")
        return False


def extract_tar_bz2(archive_path, dest_dir):
    """Extract a .tar.bz2 archive."""
    print(f"  Extracting {os.path.basename(archive_path)}...")
    ensure_dir(dest_dir)
    with tarfile.open(archive_path, 'r:bz2') as tar:
        tar.extractall(dest_dir)


def extract_zip(archive_path, dest_dir):
    """Extract a .zip archive."""
    print(f"  Extracting {os.path.basename(archive_path)}...")
    ensure_dir(dest_dir)
    with zipfile.ZipFile(archive_path, 'r') as zip_ref:
        zip_ref.extractall(dest_dir)


def run_command(cmd, description="command"):
    """Run a shell command and handle errors."""
    print(f"  Running: {description}")
    print(f"    {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    ✗ Failed (return code {result.returncode})")
        if result.stderr:
            print(f"    stderr: {result.stderr[:500]}")
        return False
    print(f"    ✓ Complete")
    return True


def find_extracted_subdir(extract_dir, pattern):
    """
    Find the actual extracted subdirectory after unpacking a Bitbucket tarball.
    Discovers directory by glob rather than hardcoding the commit hash.
    """
    matches = glob.glob(os.path.join(extract_dir, pattern))
    matches = [m for m in matches if os.path.isdir(m)]
    if not matches:
        return None
    if len(matches) > 1:
        matches.sort(key=os.path.getmtime, reverse=True)
        print(f"    Warning: multiple directories match '{pattern}', using: {os.path.basename(matches[0])}")
    return matches[0]


# =============================================================================
# Source Mapping — tracks which database each gene header came from
# =============================================================================

def _seq_id(header):
    """
    Extract the sequence ID (first whitespace-delimited token) from a FASTA header.
    This matches what CD-HIT preserves in .clstr files regardless of -d setting.
    """
    return header.split()[0]


def scan_fasta_headers(fasta_path):
    """Extract all FASTA headers (without '>') from a file."""
    headers = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:].strip())
    return headers


def build_source_map(tagged_fasta_list):
    """
    Build a mapping of sequence_id → source_database.
    Keys on _seq_id() to match CD-HIT cluster file gene names.
    """
    source_map = {}
    for fasta_entry, source in tagged_fasta_list:
        paths = fasta_entry if isinstance(fasta_entry, list) else [fasta_entry]
        for path in paths:
            if path and os.path.exists(path):
                for header in scan_fasta_headers(path):
                    source_map[_seq_id(header)] = source
    return source_map


def save_source_map(source_map, output_path):
    """Save source map as a TSV file."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Sequence_ID', 'Source_Database'])
        for seq_id, source in source_map.items():
            writer.writerow([seq_id, source])
    print(f"  Source map: {output_path} ({len(source_map)} entries)")


def load_source_map(source_map_path):
    """Load source map from TSV file."""
    source_map = {}
    with open(source_map_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = row.get('Sequence_ID') or row.get('Gene_Header', '')
            source_map[key] = row['Source_Database']
    return source_map


# =============================================================================
# Database Download Functions
# =============================================================================

def download_megares(input_dir):
    """Download and extract MEGARes database (nucleotide)."""
    archive_path = os.path.join(input_dir, f"megares_{MEGARES_VERSION}.zip")
    extract_dir = os.path.join(input_dir, "megares")
    fasta_path = os.path.join(
        extract_dir, f"megares_{MEGARES_VERSION}",
        f"megares_database_{MEGARES_VERSION}.fasta"
    )
    if os.path.exists(fasta_path):
        print(f"  Already exists: {fasta_path}")
        return fasta_path
    url = f"https://www.meglab.org/downloads/megares_{MEGARES_VERSION}.zip"
    if not os.path.exists(archive_path):
        if not download_file(url, archive_path, "MEGARes"):
            return None
    extract_zip(archive_path, extract_dir)
    ensure_trailing_newline(fasta_path)
    return fasta_path


def download_card(input_dir):
    """Download and extract CARD database (nucleotide FASTA files only)."""
    archive_path = os.path.join(input_dir, f"broadstreet-{CARD_VERSION}.tar.bz2")
    extract_dir = os.path.join(input_dir, "card")
    sentinel_path = os.path.join(extract_dir, ".sentinel")
    # These are the NUCLEOTIDE fasta files — CARD also ships protein_fasta_*
    # files which we intentionally exclude.
    card_files = [
        "nucleotide_fasta_protein_homolog_model.fasta",
        "nucleotide_fasta_protein_overexpression_model.fasta",
        "nucleotide_fasta_protein_variant_model.fasta",
        "nucleotide_fasta_rRNA_gene_variant_model.fasta",
    ]
    if os.path.exists(sentinel_path):
        print(f"  Already extracted")
        return [os.path.join(extract_dir, f) for f in card_files
                if os.path.exists(os.path.join(extract_dir, f))]
    url = f"https://card.mcmaster.ca/download/0/broadstreet-{CARD_VERSION}.tar.bz2"
    if not os.path.exists(archive_path):
        if not download_file(url, archive_path, "CARD"):
            return None
    extract_tar_bz2(archive_path, extract_dir)
    fasta_paths = []
    for fname in card_files:
        fpath = os.path.join(extract_dir, fname)
        if os.path.exists(fpath):
            ensure_trailing_newline(fpath)
            fasta_paths.append(fpath)
    Path(sentinel_path).touch()
    return fasta_paths


def download_resfinder(input_dir):
    """Download and extract ResFinder database (nucleotide)."""
    archive_path = os.path.join(input_dir, f"resfinder-{RESFINDER_VERSION}.tar.bz2")
    extract_dir = os.path.join(input_dir, "resfinder")
    sentinel_path = os.path.join(extract_dir, ".sentinel")

    if os.path.exists(sentinel_path):
        print(f"  Already extracted")
        return _find_resfinder_fasta(extract_dir)

    url = (f"https://bitbucket.org/genomicepidemiology/resfinder_db/get/"
           f"resfinder-{RESFINDER_VERSION}.tar.bz2")
    if not os.path.exists(archive_path):
        if not download_file(url, archive_path, "ResFinder"):
            return None

    extract_tar_bz2(archive_path, extract_dir)
    fasta_path = _find_resfinder_fasta(extract_dir)
    if fasta_path is None:
        print(f"    ✗ Could not find ResFinder FASTA files after extraction")
        return None
    Path(sentinel_path).touch()
    return fasta_path


def _find_resfinder_fasta(extract_dir):
    """Locate or build the ResFinder combined FASTA from an extracted directory."""
    subdir = find_extracted_subdir(extract_dir, "genomicepidemiology-resfinder_db-*")
    if subdir is None:
        print(f"    ✗ No directory matching 'genomicepidemiology-resfinder_db-*' in {extract_dir}")
        return None
    print(f"    Found ResFinder directory: {os.path.basename(subdir)}")

    all_fsa = os.path.join(subdir, "all.fsa")
    if os.path.exists(all_fsa) and os.path.getsize(all_fsa) > 0:
        ensure_trailing_newline(all_fsa)
        return all_fsa

    print(f"    all.fsa not found, concatenating individual .fsa files...")
    individual_fsa = sorted(glob.glob(os.path.join(subdir, "*.fsa")))
    if not individual_fsa:
        individual_fsa = sorted(glob.glob(os.path.join(subdir, "*.fasta")))
    if not individual_fsa:
        print(f"    ✗ No .fsa or .fasta files found in {subdir}")
        return None

    with open(all_fsa, 'wb') as outfile:
        for fsa in individual_fsa:
            ensure_trailing_newline(fsa)
            print(f"      + {os.path.basename(fsa)}")
            with open(fsa, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
    print(f"    ✓ Built all.fsa from {len(individual_fsa)} files")
    return all_fsa


def download_pointfinder(input_dir):
    """Download and extract PointFinder database (nucleotide)."""
    archive_path = os.path.join(input_dir, f"pointfinder-{POINTFINDER_VERSION}.tar.bz2")
    extract_dir = os.path.join(input_dir, "pointfinder")
    sentinel_path = os.path.join(extract_dir, ".sentinel")

    if os.path.exists(sentinel_path):
        print(f"  Already extracted")
        return _find_pointfinder_fastas(extract_dir)

    url = (f"https://bitbucket.org/genomicepidemiology/pointfinder_db/get/"
           f"resfinder-{POINTFINDER_VERSION}.tar.bz2")
    if not os.path.exists(archive_path):
        if not download_file(url, archive_path, "PointFinder"):
            return None

    extract_tar_bz2(archive_path, extract_dir)
    fasta_paths = _find_pointfinder_fastas(extract_dir)
    if not fasta_paths:
        print(f"    ✗ Could not find PointFinder FASTA files after extraction")
        return None
    Path(sentinel_path).touch()
    return fasta_paths


def _find_pointfinder_fastas(extract_dir):
    """Locate PointFinder species FASTA files from an extracted directory."""
    subdir = find_extracted_subdir(extract_dir, "genomicepidemiology-pointfinder_db-*")
    if subdir is None:
        print(f"    ✗ No directory matching 'genomicepidemiology-pointfinder_db-*' in {extract_dir}")
        return None
    print(f"    Found PointFinder directory: {os.path.basename(subdir)}")

    fasta_paths = []
    for species in POINTFINDER_SPECIES:
        fpath = os.path.join(subdir, species, f"{species}.fsa")
        if os.path.exists(fpath):
            ensure_trailing_newline(fpath)
            fasta_paths.append(fpath)
        else:
            alt_files = glob.glob(os.path.join(subdir, species, "*.fsa"))
            if alt_files:
                for af in alt_files:
                    ensure_trailing_newline(af)
                    fasta_paths.append(af)
            else:
                print(f"    Warning: no FASTA for species '{species}'")
    print(f"    Found {len(fasta_paths)} species FASTA files")
    return fasta_paths


def download_ncbi_amrfinder(input_dir):
    """Download NCBI AMRFinderPlus database (AMR_CDS.fa — coding DNA sequences, nucleotide)."""
    ncbi_dir = os.path.join(input_dir, "ncbi")
    fasta_path = os.path.join(ncbi_dir, "AMR_CDS.fa")
    if os.path.exists(fasta_path):
        print(f"  Already exists: {fasta_path}")
        return fasta_path
    ensure_dir(ncbi_dir)
    url = (f"ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/"
           f"AMRFinderPlus/database/{NCBI_VERSION}/AMR_CDS.fa")
    if not download_file(url, fasta_path, "NCBI AMRFinderPlus"):
        return None
    ensure_trailing_newline(fasta_path)
    return fasta_path


# =============================================================================
# Processing Functions
# =============================================================================

def concatenate_fasta_files(fasta_files, output_path):
    """
    Concatenate multiple FASTA files into one.
    MEGARes must be first to ensure representative status in CD-HIT.
    """
    print(f"  Output: {output_path}")
    print(f"  Input files ({len(fasta_files)}):")
    with open(output_path, 'wb') as outfile:
        for fasta in fasta_files:
            if fasta and os.path.exists(fasta):
                print(f"    + {os.path.basename(fasta)}")
                with open(fasta, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)
    print(f"  ✓ Concatenation complete")
    return output_path


def run_cdhit_dedup(input_fasta, output_fasta, threads=6, cdhit_path="cd-hit-est"):
    """Deduplicate at 100% identity with cd-hit-est (nucleotide only)."""
    cmd = [
        cdhit_path,
        "-i", input_fasta,
        "-o", output_fasta,
        "-G", "0", "-c", "1.0",
        "-AS", "0", "-AL", "0", "-aL", "1.0", "-aS", "1.0",
        "-d", "200", "-T", str(threads),
    ]
    return run_command(cmd, "CD-HIT-EST deduplication (100% identity)")


def run_cdhit_cluster(input_fasta, output_fasta, identity=0.8, threads=6,
                      cdhit_path="cd-hit-est"):
    """Cluster at specified identity threshold with cd-hit-est (nucleotide only)."""
    cmd = [
        cdhit_path,
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(identity),
        "-d", "200", "-T", str(threads),
    ]
    return run_command(cmd, f"CD-HIT-EST clustering ({identity*100:.0f}% identity)")


# =============================================================================
# Cluster file parsing (shared between Steps 2b and 4)
# =============================================================================

def parse_cluster_file(filepath):
    """Parse CD-HIT .clstr file into list of cluster dicts."""
    clusters = []
    current = None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>Cluster'):
                if current is not None:
                    clusters.append(current)
                current = {'cluster_id': int(line.split()[1]), 'sequences': []}
            else:
                seq = _parse_sequence_line(line)
                if seq and current is not None:
                    current['sequences'].append(seq)
    if current is not None:
        clusters.append(current)
    return clusters


def _parse_sequence_line(line):
    """Parse a sequence line from a CD-HIT .clstr file."""
    match = re.match(r'(\d+)\s+(\d+)nt,\s+>(.+?)\.\.\.\s*(.*)', line)
    if not match:
        return None
    length = int(match.group(2))
    gene_name = match.group(3)
    sim_part = match.group(4).strip()
    if sim_part == '*':
        is_representative = True
        similarity = 100.0
    else:
        is_representative = False
        sim_match = re.search(r'at\s+[+-]/(\d+\.?\d*)%', sim_part)
        similarity = float(sim_match.group(1)) if sim_match else None
    return {
        'gene_name': gene_name,
        'length': length,
        'similarity': similarity,
        'is_representative': is_representative,
        'is_meg': gene_name.startswith('MEG_'),
    }


# =============================================================================
# Dedup Report — Step 2b
# =============================================================================

def generate_dedup_report(dedup_clstr_path, source_map, output_path):
    """Parse the 100% dedup .clstr to report MEGARes gene redundancy across databases."""
    clusters = parse_cluster_file(dedup_clstr_path)
    records = []
    for cluster in clusters:
        seqs = cluster['sequences']
        meg_genes = [s for s in seqs if s['is_meg']]
        non_meg_genes = [s for s in seqs if not s['is_meg']]
        if not meg_genes:
            continue
        for meg in meg_genes:
            if not non_meg_genes:
                records.append({
                    'MEGARes_Accession': meg['gene_name'],
                    'MEGARes_Length': meg['length'],
                    'Num_Identical_Genes': 0,
                    'Identical_Genes': '', 'Source_Databases': '',
                    'Databases_Represented': '', 'Cluster_ID': cluster['cluster_id'],
                })
                continue
            match_genes, match_sources, source_set = [], [], set()
            for ext in non_meg_genes:
                match_genes.append(ext['gene_name'])
                source = source_map.get(_seq_id(ext['gene_name']), 'Unknown')
                match_sources.append(source)
                source_set.add(source)
            records.append({
                'MEGARes_Accession': meg['gene_name'],
                'MEGARes_Length': meg['length'],
                'Num_Identical_Genes': len(match_genes),
                'Identical_Genes': ' ; '.join(match_genes),
                'Source_Databases': ' ; '.join(match_sources),
                'Databases_Represented': ', '.join(sorted(source_set)),
                'Cluster_ID': cluster['cluster_id'],
            })

    columns = ['MEGARes_Accession', 'MEGARes_Length', 'Num_Identical_Genes',
               'Identical_Genes', 'Source_Databases', 'Databases_Represented', 'Cluster_ID']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)

    with_matches = sum(1 for r in records if r['Num_Identical_Genes'] > 0)
    without_matches = sum(1 for r in records if r['Num_Identical_Genes'] == 0)
    print(f"  MEGARes genes total:        {len(records)}")
    print(f"    With external duplicates:  {with_matches}")
    print(f"    Unique to MEGARes:         {without_matches}")
    return records


# =============================================================================
# Cluster Annotation Parsing — Step 4
# =============================================================================

def parse_meg_annotation(gene_name):
    """Parse MEGARes annotation: MEG_ID|Type|Class|Mechanism|Group|SNP"""
    if gene_name.startswith('>'):
        gene_name = gene_name[1:]
    parts = gene_name.split('|')
    return {
        'Type': parts[1] if len(parts) > 1 else '',
        'Class': parts[2] if len(parts) > 2 else '',
        'Mechanism': parts[3] if len(parts) > 3 else '',
        'Group': parts[4] if len(parts) > 4 else '',
        'SNP': parts[5] if len(parts) > 5 else '',
    }


def check_annotation_consistency(meg_genes):
    """Check annotation consistency among MEG genes in a cluster."""
    if not meg_genes:
        return 'no_meg', None, None
    annotations = [parse_meg_annotation(g['gene_name']) for g in meg_genes]
    first = annotations[0]
    if all(a == first for a in annotations[1:]):
        return 'consistent', first, None
    non_snp_keys = ['Type', 'Class', 'Mechanism', 'Group']
    first_non_snp = {k: first[k] for k in non_snp_keys}
    only_snp_differs = all(
        {k: a[k] for k in non_snp_keys} == first_non_snp for a in annotations[1:]
    )
    if only_snp_differs:
        result = first.copy()
        result['SNP'] = ''
        return 'snp_only_difference', result, None
    unique_annotations = []
    seen = set()
    for a in annotations:
        key = (a['Type'], a['Class'], a['Mechanism'], a['Group'])
        if key not in seen:
            seen.add(key)
            unique_annotations.append(a)
    return 'mismatch', None, unique_annotations


def process_clusters_for_annotations(clusters, source_map=None):
    """Process clusters into annotation records for non-MEG genes."""
    records = []
    for cluster in clusters:
        seqs = cluster['sequences']
        cluster_id = cluster['cluster_id']
        meg_genes = [s for s in seqs if s['is_meg']]
        non_meg_genes = [s for s in seqs if not s['is_meg']]
        if not non_meg_genes:
            continue
        if meg_genes:
            status, annotation, mismatch_detail = check_annotation_consistency(meg_genes)
            for gene in non_meg_genes:
                record = {
                    'Gene_Accession': gene['gene_name'], 'Length': gene['length'],
                    'Similarity': gene['similarity'], 'Cluster_ID': cluster_id,
                }
                if source_map:
                    record['Source_Database'] = source_map.get(_seq_id(gene['gene_name']), '')
                if status == 'consistent':
                    record.update({
                        'Type': annotation['Type'], 'Class': annotation['Class'],
                        'Mechanism': annotation['Mechanism'], 'Group': annotation['Group'],
                        'SNP': annotation['SNP'],
                        'Annotation_Status': 'Inherited annotations', 'Mismatch_Groups': '',
                    })
                elif status == 'snp_only_difference':
                    record.update({
                        'Type': annotation['Type'], 'Class': annotation['Class'],
                        'Mechanism': annotation['Mechanism'], 'Group': annotation['Group'],
                        'SNP': '',
                        'Annotation_Status': 'Confirm SNP difference', 'Mismatch_Groups': '',
                    })
                else:
                    conflict_strs = [
                        f"{a['Type']}|{a['Class']}|{a['Mechanism']}|{a['Group']}"
                        for a in mismatch_detail
                    ]
                    record.update({
                        'Type': '', 'Class': '', 'Mechanism': '', 'Group': '', 'SNP': '',
                        'Annotation_Status': 'Mismatch in cluster',
                        'Mismatch_Groups': ' ; '.join(conflict_strs),
                    })
                records.append(record)
        else:
            for gene in non_meg_genes:
                record = {
                    'Gene_Accession': gene['gene_name'], 'Length': gene['length'],
                    'Similarity': gene['similarity'],
                    'Type': '', 'Class': '', 'Mechanism': '', 'Group': '', 'SNP': '',
                    'Annotation_Status': 'Needs new annotations', 'Mismatch_Groups': '',
                    'Cluster_ID': cluster_id,
                }
                if source_map:
                    record['Source_Database'] = source_map.get(_seq_id(gene['gene_name']), '')
                records.append(record)
    return records


def write_parsed_output(records, output_path, include_source=False):
    """Write parsed cluster records to TSV."""
    columns = [
        'Gene_Accession', 'Length', 'Similarity',
        'Type', 'Class', 'Mechanism', 'Group', 'SNP',
        'Annotation_Status', 'Mismatch_Groups', 'Cluster_ID',
    ]
    if include_source:
        columns.insert(3, 'Source_Database')
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


# =============================================================================
# Version Checking
# =============================================================================

def check_latest_versions():
    """Check for latest database versions via API queries."""
    print("=" * 60)
    print("Checking for database updates")
    print("=" * 60)
    print()

    print("  CARD:")
    print(f"    Current: {CARD_VERSION}")
    try:
        req = Request("https://card.mcmaster.ca/latest/data", headers={"User-Agent": "Mozilla/5.0"})
        resp = urlopen(req, timeout=15)
        ver_match = re.search(r'broadstreet-(v[\d.]+)', resp.geturl())
        if ver_match:
            latest = ver_match.group(1)
            flag = " ← UPDATE AVAILABLE" if latest != CARD_VERSION else " (up to date)"
            print(f"    Latest:  {latest}{flag}")
        else:
            print(f"    Check:   https://card.mcmaster.ca/download")
    except Exception as e:
        print(f"    Check failed: {e}")

    print("  ResFinder:")
    print(f"    Current: {RESFINDER_VERSION}")
    try:
        req = Request(
            "https://api.bitbucket.org/2.0/repositories/genomicepidemiology/"
            "resfinder_db/refs/tags?sort=-name&pagelen=1",
            headers={"User-Agent": "Mozilla/5.0"})
        data = json.loads(urlopen(req, timeout=15).read())
        if data.get('values'):
            latest_tag = data['values'][0]['name']
            flag = " ← UPDATE AVAILABLE" if latest_tag != f"resfinder-{RESFINDER_VERSION}" else " (up to date)"
            print(f"    Latest tag: {latest_tag}{flag}")
    except Exception as e:
        print(f"    Check failed: {e}")

    print("  PointFinder:")
    print(f"    Current: {POINTFINDER_VERSION}")
    try:
        req = Request(
            "https://api.bitbucket.org/2.0/repositories/genomicepidemiology/"
            "pointfinder_db/refs/tags?sort=-name&pagelen=1",
            headers={"User-Agent": "Mozilla/5.0"})
        data = json.loads(urlopen(req, timeout=15).read())
        if data.get('values'):
            latest_tag = data['values'][0]['name']
            flag = " ← UPDATE AVAILABLE" if latest_tag != f"resfinder-{POINTFINDER_VERSION}" else " (up to date)"
            print(f"    Latest tag: {latest_tag}{flag}")
    except Exception as e:
        print(f"    Check failed: {e}")

    print("  NCBI AMRFinderPlus:")
    print(f"    Current: {NCBI_VERSION}")
    print(f"    Manual: ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/")

    print("  MEGARes:")
    print(f"    Current: {MEGARES_VERSION}")
    print(f"    Manual:  https://www.meglab.org/megares")

    print()
    print("  To update: edit the version constants at the top of this script.")


# =============================================================================
# Main Pipeline
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="AMR Database Aggregation, Clustering, and Annotation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
  1.  Download databases → concatenate (MEGARes first) → input/all.fsa
      Also produces input/source_map.tsv
  2.  Deduplicate at 100%% identity → input/dedup.fsa
  2b. (Optional) Parse dedup clusters → output/megares_dedup_report.tsv
  3.  Cluster at 80%% similarity → output/megares-update-clustering.fsa
  4.  Parse clusters for annotation → output/..._parsed.tsv

Note: All databases are nucleotide-only. BacMet is excluded because it
only provides protein sequences (incompatible with cd-hit-est).

Step flags:
  --download-only    Stop after Step 1
  --skip-download    Skip Step 1 (requires existing input/all.fsa + source_map.tsv)
  --skip-dedup       Skip Step 2 (requires existing input/dedup.fsa)
  --cluster-only     Skip Steps 1-2, run from existing dedup.fsa
  --parse-only       Skip Steps 1-3, parse existing .clstr file
  --no-dedup-report  Skip Step 2b
  --check-versions   Check for database updates and exit

Examples:
  %(prog)s                                    # Full pipeline
  %(prog)s --download-only                    # Just download + concatenate
  %(prog)s --skip-download                    # Steps 2-4 from existing data
  %(prog)s --parse-only --cluster-file path/to/file.clstr
  %(prog)s --check-versions                   # Check for updates
  %(prog)s --identity 0.9 --threads 12        # Custom params
        """
    )

    parser.add_argument("--input-dir", "-i", default="./input")
    parser.add_argument("--output-dir", "-o", default="./output")
    parser.add_argument("--cdhit-path", default="cd-hit-est")
    parser.add_argument("--threads", "-T", type=int, default=6)
    parser.add_argument("--identity", type=float, default=0.8)
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-dedup", action="store_true")
    parser.add_argument("--cluster-only", action="store_true")
    parser.add_argument("--parse-only", action="store_true")
    parser.add_argument("--cluster-file", default=None)
    parser.add_argument("--no-dedup-report", action="store_true")
    parser.add_argument("--check-versions", action="store_true")

    args = parser.parse_args()

    if args.check_versions:
        check_latest_versions()
        return

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    ensure_dir(input_dir)
    ensure_dir(output_dir)

    combined_fasta = os.path.join(input_dir, "all.fsa")
    source_map_path = os.path.join(input_dir, "source_map.tsv")
    dedup_fasta = os.path.join(input_dir, "dedup.fsa")
    dedup_clstr = dedup_fasta + ".clstr"
    dedup_report = os.path.join(output_dir, "megares_dedup_report.tsv")
    cluster_fasta = os.path.join(output_dir, "megares-update-clustering.fsa")
    cluster_file = cluster_fasta + ".clstr"
    parsed_tsv = os.path.join(output_dir, "megares-update-clustering_parsed.tsv")

    print("=" * 60)
    print("AMR Database Aggregation & Clustering Pipeline")
    print("=" * 60)
    print(f"  Input dir:  {input_dir}")
    print(f"  Output dir: {output_dir}")
    print()

    # --parse-only shortcut
    if args.parse_only:
        cf = args.cluster_file or cluster_file
        if not os.path.exists(cf):
            print(f"Error: cluster file not found: {cf}"); sys.exit(1)
        source_map = None
        if os.path.exists(source_map_path):
            source_map = load_source_map(source_map_path)
            print(f"  Loaded source map ({len(source_map)} entries)")
        print(f"[Step 4] Parsing cluster file: {cf}")
        clusters = parse_cluster_file(cf)
        records = process_clusters_for_annotations(clusters, source_map)
        _print_parse_summary(clusters, records)
        write_parsed_output(records, parsed_tsv, include_source=(source_map is not None))
        print(f"  Output: {parsed_tsv}")
        return

    # Step 1: Download and concatenate
    source_map = {}
    if not args.skip_download and not args.cluster_only and not args.skip_dedup:
        print("[Step 1] Downloading databases and concatenating")
        print("-" * 60)

        tagged_fastas = []

        print("\n  --- MEGARes (first — ensures representative status) ---")
        megares_fasta = download_megares(input_dir)
        if megares_fasta:
            tagged_fastas.append((megares_fasta, SOURCE_MEGARES))

        print("\n  --- CARD (nucleotide) ---")
        card_fastas = download_card(input_dir)
        if card_fastas:
            tagged_fastas.append((card_fastas, SOURCE_CARD))

        print("\n  --- ResFinder (nucleotide) ---")
        resfinder_fasta = download_resfinder(input_dir)
        if resfinder_fasta:
            tagged_fastas.append((resfinder_fasta, SOURCE_RESFINDER))

        print("\n  --- PointFinder (nucleotide) ---")
        pointfinder_fastas = download_pointfinder(input_dir)
        if pointfinder_fastas:
            tagged_fastas.append((pointfinder_fastas, SOURCE_POINTFINDER))

        print("\n  --- NCBI AMRFinderPlus (nucleotide CDS) ---")
        ncbi_fasta = download_ncbi_amrfinder(input_dir)
        if ncbi_fasta:
            tagged_fastas.append((ncbi_fasta, SOURCE_NCBI))

        if not tagged_fastas:
            print("  Error: No FASTA files downloaded"); sys.exit(1)

        print("\n  --- Building source map ---")
        source_map = build_source_map(tagged_fastas)
        save_source_map(source_map, source_map_path)

        flat_fasta_files = []
        for entry, _label in tagged_fastas:
            if isinstance(entry, list):
                flat_fasta_files.extend(entry)
            else:
                flat_fasta_files.append(entry)

        print("\n  --- Concatenating (MEGARes first) ---")
        concatenate_fasta_files(flat_fasta_files, combined_fasta)

        if args.download_only:
            print(f"\n✓ Download complete.")
            print(f"  Combined FASTA: {combined_fasta}")
            print(f"  Source map:     {source_map_path}")
            return
    else:
        if os.path.exists(source_map_path):
            source_map = load_source_map(source_map_path)
            print(f"  Loaded source map ({len(source_map)} entries)")
        if not os.path.exists(combined_fasta) and not args.cluster_only and not args.skip_dedup:
            print(f"Error: {combined_fasta} not found"); sys.exit(1)

    # Step 2: Deduplicate
    if not args.skip_dedup and not args.cluster_only:
        print()
        print("[Step 2] Deduplication (100% identity)")
        print("-" * 60)
        if not run_cdhit_dedup(combined_fasta, dedup_fasta,
                               threads=args.threads, cdhit_path=args.cdhit_path):
            print("Error: Deduplication failed"); sys.exit(1)
        print(f"  Deduplicated FASTA: {dedup_fasta}")

        if not args.no_dedup_report and os.path.exists(dedup_clstr):
            print()
            print("[Step 2b] MEGARes dedup redundancy report")
            print("-" * 60)
            if not source_map and os.path.exists(source_map_path):
                source_map = load_source_map(source_map_path)
            generate_dedup_report(dedup_clstr, source_map, dedup_report)
            print(f"  Report: {dedup_report}")
    else:
        if not os.path.exists(dedup_fasta):
            print(f"Error: {dedup_fasta} not found"); sys.exit(1)

    # Step 3: Cluster
    print()
    print(f"[Step 3] Clustering ({args.identity*100:.0f}% identity)")
    print("-" * 60)
    if not run_cdhit_cluster(dedup_fasta, cluster_fasta,
                             identity=args.identity, threads=args.threads,
                             cdhit_path=args.cdhit_path):
        print("Error: Clustering failed"); sys.exit(1)
    print(f"  Clustered FASTA: {cluster_fasta}")

    # Step 4: Parse clusters
    print()
    print("[Step 4] Parsing clusters for annotation inheritance")
    print("-" * 60)
    if not source_map and os.path.exists(source_map_path):
        source_map = load_source_map(source_map_path)
    clusters = parse_cluster_file(cluster_file)
    records = process_clusters_for_annotations(clusters, source_map or None)
    _print_parse_summary(clusters, records)
    write_parsed_output(records, parsed_tsv, include_source=bool(source_map))
    print(f"  Output: {parsed_tsv}")

    print()
    print("=" * 60)
    print("Pipeline Complete")
    print("=" * 60)
    print(f"  Combined FASTA:    {combined_fasta}")
    print(f"  Source map:        {source_map_path}")
    print(f"  Deduplicated:      {dedup_fasta}")
    if not args.no_dedup_report:
        print(f"  Dedup report:      {dedup_report}")
    print(f"  Clustered:         {cluster_fasta}")
    print(f"  Parsed output:     {parsed_tsv}")


def _print_parse_summary(clusters, records):
    print(f"  Total clusters: {len(clusters)}")
    print(f"  Non-MEG gene records: {len(records)}")
    for label, status in [("Inherited annotations", "Inherited annotations"),
                          ("Confirm SNP difference", "Confirm SNP difference"),
                          ("Mismatch in cluster", "Mismatch in cluster"),
                          ("Needs new annotations", "Needs new annotations")]:
        count = sum(1 for r in records if r['Annotation_Status'] == status)
        print(f"    {label:25s} {count}")


if __name__ == '__main__':
    main()