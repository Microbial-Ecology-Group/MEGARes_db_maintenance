#!/usr/bin/env python3
"""
Parse CD-HIT cluster file to identify non-MEGARes genes and inherit/flag annotations.

This script:
1. Identifies clusters containing both MEG_ and non-MEG_ genes
   - If MEG annotations match: inherits annotations ("Inherited annotations")
   - If only SNP differs: inherits non-SNP annotations ("Confirm SNP difference")
   - If other fields differ: flags as "Mismatch in cluster" with conflicting values
2. Identifies clusters with only non-MEG_ genes as "Needs new annotations"
3. Outputs a TSV file with gene information and annotation status

Output columns:
  Gene_Accession     Non-MEGARes gene identifier
  Length             Sequence length (nt)
  Similarity         Percent similarity to cluster representative
  Type               AMR type (inherited or empty)
  Class              AMR class (inherited or empty)
  Mechanism          Resistance mechanism (inherited or empty)
  Group              Gene group (inherited or empty)
  SNP                SNP requirement (inherited or empty)
  Annotation_Status  Inheritance outcome (see above)
  Mismatch_Groups    For mismatches: semicolon-separated Type|Class|Mechanism|Group
                     combos found among MEG genes in the cluster
  Cluster_ID         CD-HIT cluster number
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def parse_meg_annotation(gene_name):
    """
    Parse MEGARes gene annotation from gene name.
    Format: >MEG_ID|Type|Class|Mechanism|Group|SNP (SNP is optional)
    """
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


def parse_sequence_line(line):
    """
    Parse a sequence line from the cluster file.
    Returns dict with gene_name, length, similarity, is_representative, is_meg.
    """
    match = re.match(r'(\d+)\s+(\d+)nt,\s+>(.+?)\.\.\.\s*(.*)', line)
    if not match:
        return None

    length = int(match.group(2))
    gene_name = match.group(3)
    similarity_part = match.group(4).strip()

    if similarity_part == '*':
        is_representative = True
        similarity = 100.0
    else:
        is_representative = False
        sim_match = re.search(r'at\s+[+-]/(\d+\.?\d*)%', similarity_part)
        similarity = float(sim_match.group(1)) if sim_match else None

    return {
        'gene_name': gene_name,
        'length': length,
        'similarity': similarity,
        'is_representative': is_representative,
        'is_meg': gene_name.startswith('MEG_'),
    }


def parse_clusters(filepath):
    """Parse the entire cluster file and return list of clusters."""
    clusters = []
    current_cluster = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>Cluster'):
                if current_cluster is not None:
                    clusters.append(current_cluster)
                current_cluster = {
                    'cluster_id': int(line.split()[1]),
                    'sequences': [],
                }
            else:
                seq_info = parse_sequence_line(line)
                if seq_info and current_cluster is not None:
                    current_cluster['sequences'].append(seq_info)

    if current_cluster is not None:
        clusters.append(current_cluster)

    return clusters


def check_annotation_consistency(meg_genes):
    """
    Check if all MEG genes in a cluster have consistent annotations.

    Returns (status, annotation, mismatch_detail):
      - 'consistent', annotation_dict, None
      - 'snp_only_difference', annotation_dict (SNP=''), None
      - 'mismatch', None, list of unique annotation dicts
      - 'no_meg', None, None
    """
    if not meg_genes:
        return 'no_meg', None, None

    annotations = [parse_meg_annotation(g['gene_name']) for g in meg_genes]
    first = annotations[0]

    # Full match
    if all(a == first for a in annotations[1:]):
        return 'consistent', first, None

    # Check if only SNP differs
    non_snp_keys = ['Type', 'Class', 'Mechanism', 'Group']
    first_non_snp = {k: first[k] for k in non_snp_keys}
    only_snp_differs = all(
        {k: a[k] for k in non_snp_keys} == first_non_snp
        for a in annotations[1:]
    )

    if only_snp_differs:
        result = first.copy()
        result['SNP'] = ''
        return 'snp_only_difference', result, None

    # Real mismatch — collect unique annotation combos for reporting
    unique_annotations = []
    seen = set()
    for a in annotations:
        key = (a['Type'], a['Class'], a['Mechanism'], a['Group'])
        if key not in seen:
            seen.add(key)
            unique_annotations.append(a)

    return 'mismatch', None, unique_annotations


def process_clusters(clusters):
    """
    Process clusters and generate output records for non-MEG genes.

    Annotation_Status values:
      'Inherited annotations'  — MEG annotations consistent, inherited
      'Confirm SNP difference' — only SNP field differs among MEG genes
      'Mismatch in cluster'    — non-SNP annotations differ; see Mismatch_Groups
      'Needs new annotations'  — no MEG genes in cluster
    """
    output_records = []

    for cluster in clusters:
        sequences = cluster['sequences']
        cluster_id = cluster['cluster_id']

        meg_genes = [s for s in sequences if s['is_meg']]
        non_meg_genes = [s for s in sequences if not s['is_meg']]

        # Skip MEG-only clusters (not of interest)
        if not non_meg_genes:
            continue

        if meg_genes:
            status, annotation, mismatch_detail = check_annotation_consistency(meg_genes)

            for gene in non_meg_genes:
                record = {
                    'Gene_Accession': gene['gene_name'],
                    'Length': gene['length'],
                    'Similarity': gene['similarity'],
                    'Cluster_ID': cluster_id,
                }

                if status == 'consistent':
                    record.update({
                        'Type': annotation['Type'],
                        'Class': annotation['Class'],
                        'Mechanism': annotation['Mechanism'],
                        'Group': annotation['Group'],
                        'SNP': annotation['SNP'],
                        'Annotation_Status': 'Inherited annotations',
                        'Mismatch_Groups': '',
                    })

                elif status == 'snp_only_difference':
                    record.update({
                        'Type': annotation['Type'],
                        'Class': annotation['Class'],
                        'Mechanism': annotation['Mechanism'],
                        'Group': annotation['Group'],
                        'SNP': '',
                        'Annotation_Status': 'Confirm SNP difference',
                        'Mismatch_Groups': '',
                    })

                else:  # mismatch
                    conflict_strs = [
                        f"{a['Type']}|{a['Class']}|{a['Mechanism']}|{a['Group']}"
                        for a in mismatch_detail
                    ]
                    record.update({
                        'Type': '',
                        'Class': '',
                        'Mechanism': '',
                        'Group': '',
                        'SNP': '',
                        'Annotation_Status': 'Mismatch in cluster',
                        'Mismatch_Groups': ' ; '.join(conflict_strs),
                    })

                output_records.append(record)

        else:
            # No MEG genes — needs new annotations
            for gene in non_meg_genes:
                output_records.append({
                    'Gene_Accession': gene['gene_name'],
                    'Length': gene['length'],
                    'Similarity': gene['similarity'],
                    'Type': '',
                    'Class': '',
                    'Mechanism': '',
                    'Group': '',
                    'SNP': '',
                    'Annotation_Status': 'Needs new annotations',
                    'Mismatch_Groups': '',
                    'Cluster_ID': cluster_id,
                })

    return output_records


def write_output(records, output_path):
    """Write records to TSV file."""
    columns = [
        'Gene_Accession', 'Length', 'Similarity',
        'Type', 'Class', 'Mechanism', 'Group', 'SNP',
        'Annotation_Status', 'Mismatch_Groups', 'Cluster_ID',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(
        description='Parse CD-HIT cluster file to extract non-MEGARes genes with annotations.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s megares-update-clustering.fsa.clstr
  %(prog)s megares-update-clustering.fsa.clstr -o custom_output.tsv

Output columns:
  Gene_Accession     The gene identifier (without '>' prefix)
  Length             Sequence length in nucleotides
  Similarity         Percent similarity to cluster representative
  Type               AMR type (e.g., Drugs, Metals, Biocides)
  Class              AMR class (e.g., Rifampin, Glycopeptides)
  Mechanism          Resistance mechanism
  Group              Gene group identifier
  SNP                SNP confirmation requirement (if applicable)
  Annotation_Status  One of:
                       'Inherited annotations' - MEG genes match, annotations inherited
                       'Confirm SNP difference' - Only SNP field differs among MEG genes
                       'Mismatch in cluster' - Non-SNP annotations differ among MEG genes
                       'Needs new annotations' - No MEG genes in cluster
  Mismatch_Groups    For mismatches: semicolon-separated Type|Class|Mechanism|Group
                     values found among the MEG genes in that cluster
  Cluster_ID         CD-HIT cluster number
        """
    )

    parser.add_argument('cluster_file', type=str,
                        help='Path to the CD-HIT .clstr file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output TSV file path (default: <input_basename>_parsed.tsv)')

    args = parser.parse_args()

    input_path = Path(args.cluster_file)
    if not input_path.exists():
        print(f"Error: Input file '{args.cluster_file}' not found.", file=sys.stderr)
        sys.exit(1)

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.parent / f"{input_path.stem}_parsed.tsv"

    print(f"Parsing cluster file: {input_path}")

    clusters = parse_clusters(input_path)
    print(f"Found {len(clusters)} clusters")

    records = process_clusters(clusters)
    print(f"Generated {len(records)} output records")

    inherited = sum(1 for r in records if r['Annotation_Status'] == 'Inherited annotations')
    snp_diff = sum(1 for r in records if r['Annotation_Status'] == 'Confirm SNP difference')
    mismatch = sum(1 for r in records if r['Annotation_Status'] == 'Mismatch in cluster')
    needs_new = sum(1 for r in records if r['Annotation_Status'] == 'Needs new annotations')

    print(f"  - Inherited annotations:  {inherited}")
    print(f"  - Confirm SNP difference: {snp_diff}")
    print(f"  - Mismatch in cluster:    {mismatch}")
    print(f"  - Needs new annotations:  {needs_new}")

    write_output(records, output_path)
    print(f"Output written to: {output_path}")


if __name__ == '__main__':
    main()
