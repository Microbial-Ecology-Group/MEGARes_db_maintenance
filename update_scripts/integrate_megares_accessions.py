#!/usr/bin/env python3
"""
MEGARes Accession Integration & Header Reformatting

Integrates curated new AMR gene sequences (from the database aggregation /
clustering pipeline) into the MEGARes database, assigning them MEGARes-format
MEG_ headers that continue the existing numbering.

Inputs
------
  --original-fasta        Original MEGARes database (FASTA; single- or multi-line)
  --original-annotations  Original MEGARes annotations
                          (columns: header, type, class, mechanism, group, snp
                           [, "Can be confirmed"])
  --clustering-fasta      megares-update-clustering.fsa (contains the new sequences;
                          may be single- or multi-line)
  --new-annotations       Curated annotations for the new sequences. Must contain a
                          'Gene_Accession' column, a 'Source_Database' column, an
                          'ED_confirmation' column, and Type / Class / Mechanism /
                          Group / SNP columns.
  --original-mappings     megares_to_external_header_mappings CSV
                          (columns: Database, Source_header, MEGARes_header,
                           UpdatedHeader)
  --snpinfo-fasta         (optional) AmrPlusPlus_SNP data/SNPinfo.fasta. When given,
                          it is authoritative for the "Can be confirmed" column and
                          any disagreement with the annotation file is reported.

Behavior
--------
  1. Determine the highest existing MEG_<n> id (scanning both the original FASTA
     headers and the original annotation headers).
  2. Read the curated new-gene annotations; drop rows whose ED_confirmation == 'remove'.
  3. Pull each surviving new gene's sequence out of the clustering FASTA, matched on
     the first whitespace-delimited token of the accession / header (this mirrors the
     _seq_id() logic in the aggregation pipeline, so it is robust to CD-HIT's -d
     header truncation and to trailing descriptions).
  4. Assign new MEG_<n> ids sequentially, starting after the current maximum.
  5. Build MEGARes-format headers:
         MEG_<n>|Type|Class|Mechanism|Group[|SNP]
     spaces -> underscores; the trailing |SNP is omitted when SNP is blank.
     The annotation-file columns keep their original spaces (headers do not).
  6. Populate the "Can be confirmed" column:
       - original rows pass their existing value through unchanged
       - new rows carrying the SNP label are set to "No" (no SNPinfo entry exists
         for them yet, so AmrPlusPlus_SNP cannot verify them)
       - new rows without the SNP label get --default-can-be-confirmed (blank by
         default: there is nothing to confirm for a non-variant gene)
       - if --snpinfo-fasta is supplied, membership in that file decides the value
         for every row and disagreements are reported
  7. Write a single-line FASTA containing every original sequence followed by the
     new sequences with their new headers.
  8. Write an updated annotation file in the same columns / delimiter as the original,
     with "Can be confirmed" appended.
  9. Write an updated external-header mappings file: every original row verbatim,
     plus one row per newly integrated gene recording its source database, its
     original accession, and its new MEG_ header.
"""

import argparse
import csv
import re
import sys


# =============================================================================
# FASTA helpers
# =============================================================================

def read_fasta(path):
    """
    Yield (header, sequence) tuples from a FASTA file.

    header is returned without the leading '>'. Multi-line sequences are joined
    into a single string, so downstream output can be written single-line.
    """
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n').rstrip('\r')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_parts)


def seq_id(header):
    """
    Return the first whitespace-delimited token of a header.

    Matches the aggregation pipeline's _seq_id(): keying on the sequence ID makes
    matching reliable regardless of trailing descriptions or CD-HIT truncation.
    """
    header = header.lstrip('>').strip()
    parts = header.split()
    return parts[0] if parts else header


def meg_id_of(header):
    """Return the MEG_<n> token at the start of a MEGARes header, or ''."""
    if not header:
        return ''
    return header.split('|')[0].strip()


def load_snpinfo_ids(path):
    """
    Return the set of MEG ids that have a curated mutation entry in SNPinfo.fasta.

    SNPinfo headers look like:
        MEG_1|Drugs|Aminoglycosides|...|A16S|Nuc:A1192G_GGATG_CGTCA|Nuc:...
    so the MEG id is simply the first pipe-delimited field. Presence in this file
    is exactly what makes a gene verifiable by AmrPlusPlus_SNP.
    """
    ids = set()
    for header, _seq in read_fasta(path):
        mid = meg_id_of(header)
        if mid:
            ids.add(mid)
    return ids


# =============================================================================
# Table (CSV/TSV) helpers
# =============================================================================

def detect_delimiter(path):
    """Guess the delimiter from the header line: prefer tab, then comma."""
    with open(path, newline='') as fh:
        first = fh.readline()
    if '\t' in first:
        return '\t'
    if ',' in first:
        return ','
    return '\t'


def load_table(path):
    """Load a delimited table into a list of dict rows. Returns (rows, fieldnames, delim)."""
    delim = detect_delimiter(path)
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return rows, fieldnames, delim


def get_ci(row, name, default=''):
    """Case-insensitive column lookup on a dict row."""
    if name in row and row[name] is not None:
        return row[name]
    lower = {k.lower(): k for k in row.keys() if k is not None}
    key = lower.get(name.lower())
    if key is not None and row[key] is not None:
        return row[key]
    return default


def has_column(fieldnames, name):
    """Case-insensitive membership test against a fieldname list."""
    return name.lower() in {f.lower() for f in fieldnames if f}


# =============================================================================
# Formatting helpers
# =============================================================================

def format_field(value):
    """Strip and convert internal whitespace runs to single underscores."""
    if value is None:
        return ''
    return re.sub(r'\s+', '_', value.strip())


def find_max_meg_id(headers, prefix="MEG_"):
    """Return the highest integer <n> appearing as <prefix><n> across the headers."""
    pat = re.compile(re.escape(prefix) + r'(\d+)')
    mx = 0
    for h in headers:
        if not h:
            continue
        for m in pat.finditer(h):
            mx = max(mx, int(m.group(1)))
    return mx


def derive_can_be_confirmed(meg_id, snp_value, snpinfo_ids, snp_label, default_blank):
    """
    Decide the "Can be confirmed" value for one accession.

    A gene is confirmable only if AmrPlusPlus_SNP has a curated mutation entry for
    it in SNPinfo.fasta. Genes with no SNP label have nothing to confirm, so they
    get the caller's default rather than a misleading "No".
    """
    is_snp = snp_value.strip().lower() == snp_label.lower()
    if not is_snp:
        return default_blank
    if snpinfo_ids is None:
        return 'No'
    return 'Yes' if meg_id in snpinfo_ids else 'No'


# =============================================================================
# Main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Integrate curated new accessions into MEGARes with MEG_-format headers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--original-fasta", "-f", required=True,
                   help="Original MEGARes database FASTA")
    p.add_argument("--original-annotations", "-a", required=True,
                   help="Original MEGARes annotations (header,type,class,mechanism,group,snp)")
    p.add_argument("--clustering-fasta", "-c", required=True,
                   help="megares-update-clustering.fsa (holds the new sequences)")
    p.add_argument("--new-annotations", "-n", required=True,
                   help="Curated annotations for the new sequences (with ED_confirmation)")
    p.add_argument("--original-mappings", "-m", default=None,
                   help="Original megares_to_external_header_mappings CSV "
                        "(Database,Source_header,MEGARes_header,UpdatedHeader)")
    p.add_argument("--snpinfo-fasta", "-s", default=None,
                   help="AmrPlusPlus_SNP data/SNPinfo.fasta. When supplied it is "
                        "authoritative for the 'Can be confirmed' column.")
    p.add_argument("--output-fasta", "-of", default="megares_updated.fasta",
                   help="Output combined single-line FASTA (default: megares_updated.fasta)")
    p.add_argument("--output-annotations", "-oa", default="megares_updated_annotations.csv",
                   help="Output combined annotation file (default: megares_updated_annotations.csv)")
    p.add_argument("--output-mappings", "-om", default="megares_to_external_header_mappings_updated.csv",
                   help="Output external-header mappings file "
                        "(default: megares_to_external_header_mappings_updated.csv)")
    p.add_argument("--meg-prefix", default="MEG_",
                   help="MEG id prefix (default: MEG_)")
    p.add_argument("--start-id", type=int, default=None,
                   help="Override the first new MEG id (default: max existing + 1)")
    p.add_argument("--accession-column", default="Gene_Accession",
                   help="Column in new annotations holding the accession (default: Gene_Accession)")
    p.add_argument("--source-column", default="Source_Database",
                   help="Column in new annotations holding the source database "
                        "(default: Source_Database)")
    p.add_argument("--confirmation-column", default="ED_confirmation",
                   help="Column marking genes to drop (default: ED_confirmation)")
    p.add_argument("--remove-value", default="remove",
                   help="Value (case-insensitive) that flags a gene for removal (default: remove)")
    p.add_argument("--can-confirm-column", default="Can be confirmed",
                   help='Name of the confirmability column (default: "Can be confirmed")')
    p.add_argument("--snp-label", default="RequiresSNPConfirmation",
                   help="SNP field value marking variant-dependent genes "
                        "(default: RequiresSNPConfirmation)")
    p.add_argument("--default-can-be-confirmed", default="",
                   help='Value for new genes with no SNP label (default: blank). '
                        'Use "No" if you prefer an explicit value on every row.')
    return p.parse_args()


def main():
    args = parse_args()

    print("=" * 60)
    print("MEGARes Accession Integration")
    print("=" * 60)

    cc_col = args.can_confirm_column

    # --- 1. Original FASTA -----------------------------------------------------
    original_records = list(read_fasta(args.original_fasta))
    original_headers = [h for h, _ in original_records]
    print(f"  Original sequences:        {len(original_records)}")

    # --- 2. Original annotations ----------------------------------------------
    orig_ann_rows, orig_fields, orig_delim = load_table(args.original_annotations)
    orig_ann_headers = [get_ci(r, 'header') for r in orig_ann_rows]
    print(f"  Original annotations:      {len(orig_ann_rows)}")

    orig_has_cc = has_column(orig_fields, cc_col)
    if not orig_has_cc:
        print(f'    Note: original annotations have no "{cc_col}" column; '
              f"it will be created")

    # --- 3. SNPinfo (optional, authoritative for confirmability) --------------
    snpinfo_ids = None
    if args.snpinfo_fasta:
        snpinfo_ids = load_snpinfo_ids(args.snpinfo_fasta)
        print(f"  SNPinfo entries:           {len(snpinfo_ids)}")

    # --- 4. Highest existing MEG id -------------------------------------------
    max_id = max(find_max_meg_id(original_headers, args.meg_prefix),
                 find_max_meg_id(orig_ann_headers, args.meg_prefix))
    next_id = (args.start_id - 1) if args.start_id is not None else max_id
    print(f"  Highest existing MEG id:   {args.meg_prefix}{max_id}")
    print(f"  New ids start at:          {args.meg_prefix}{next_id + 1}")

    # --- 5. Original mappings --------------------------------------------------
    orig_map_rows, orig_map_fields, map_delim = [], [], ','
    if args.original_mappings:
        orig_map_rows, orig_map_fields, map_delim = load_table(args.original_mappings)
        print(f"  Original mapping rows:     {len(orig_map_rows)}")
    else:
        orig_map_fields = ['Database', 'Source_header', 'MEGARes_header', 'UpdatedHeader']
        print("  Original mapping rows:     (none supplied)")

    # --- 6. Clustering FASTA -> sequence lookup by seq id ----------------------
    clustering_seqs = {}
    dup_ids = set()
    for h, s in read_fasta(args.clustering_fasta):
        sid = seq_id(h)
        if sid in clustering_seqs:
            dup_ids.add(sid)
            continue
        clustering_seqs[sid] = s
    print(f"  Clustering FASTA sequences: {len(clustering_seqs)}")
    if dup_ids:
        print(f"    Warning: {len(dup_ids)} duplicate seq-id(s) in clustering FASTA "
              f"(kept first occurrence of each)")

    # --- 7. New annotations ----------------------------------------------------
    new_rows, new_fields, _new_delim = load_table(args.new_annotations)
    print(f"  New-annotation rows:       {len(new_rows)}")

    if not has_column(new_fields, args.source_column):
        print(f'    Warning: new annotations have no "{args.source_column}" column; '
              f"mapping rows will have a blank Database")

    new_meg_records = []   # (new_header, sequence)
    new_ann_rows = []      # dicts keyed header/type/class/mechanism/group/snp/<cc_col>
    new_map_rows = []      # dicts for the mappings output
    removed = 0
    missing_seq = []
    incomplete_ann = []
    seen_ids = set()
    cc_counts = {}

    for row in new_rows:
        acc = get_ci(row, args.accession_column).strip()
        if not acc:
            continue

        conf = get_ci(row, args.confirmation_column).strip().lower()
        if conf == args.remove_value.lower():
            removed += 1
            continue

        sid = seq_id(acc)
        if sid in seen_ids:
            # already integrated this accession; skip the duplicate row
            continue

        if sid not in clustering_seqs:
            missing_seq.append(acc)
            continue

        seen_ids.add(sid)
        next_id += 1
        meg_id = f"{args.meg_prefix}{next_id}"

        # Raw values keep internal spaces (for the annotation columns / plotting);
        # the *_h versions replace spaces with underscores (for the FASTA header only).
        type_raw = get_ci(row, 'Type').strip()
        class_raw = get_ci(row, 'Class').strip()
        mech_raw = get_ci(row, 'Mechanism').strip()
        group_raw = get_ci(row, 'Group').strip()
        snp_raw = get_ci(row, 'SNP').strip()
        source_raw = get_ci(row, args.source_column).strip()

        type_h = format_field(type_raw)
        class_h = format_field(class_raw)
        mech_h = format_field(mech_raw)
        group_h = format_field(group_raw)
        snp_h = format_field(snp_raw)

        if not (type_raw and class_raw and mech_raw and group_raw):
            incomplete_ann.append(f"{meg_id}  <- {sid}")

        # Type|Class|Mechanism|Group are positional and always kept (even if blank)
        # so header positions stay aligned; SNP is optional and dropped when blank.
        header_fields = [meg_id, type_h, class_h, mech_h, group_h]
        if snp_h:
            header_fields.append(snp_h)
        new_header = '|'.join(header_fields)

        # New genes have no SNPinfo entry yet, so anything carrying the SNP label
        # is not confirmable by AmrPlusPlus_SNP.
        cc_value = derive_can_be_confirmed(
            meg_id, snp_raw, snpinfo_ids, args.snp_label,
            args.default_can_be_confirmed)
        cc_counts[cc_value or '(blank)'] = cc_counts.get(cc_value or '(blank)', 0) + 1

        new_meg_records.append((new_header, clustering_seqs[sid]))
        new_ann_rows.append({
            'header': new_header,      # underscored, matches the FASTA
            'type': type_raw,          # spaces preserved for plotting
            'class': class_raw,
            'mechanism': mech_raw,
            'group': group_raw,
            'snp': snp_raw,
            cc_col: cc_value,
        })
        new_map_rows.append({
            'Database': source_raw,
            'Source_header': acc,      # the original accession, verbatim
            'MEGARes_header': new_header,
            'UpdatedHeader': '',
        })

    # --- 8. Resolve "Can be confirmed" for the original rows ------------------
    cc_conflicts = []
    orig_cc_values = []
    for r in orig_ann_rows:
        hdr = get_ci(r, 'header')
        existing = get_ci(r, cc_col).strip() if orig_has_cc else ''
        if snpinfo_ids is not None:
            derived = derive_can_be_confirmed(
                meg_id_of(hdr), get_ci(r, 'snp'), snpinfo_ids, args.snp_label,
                args.default_can_be_confirmed)
            if orig_has_cc and existing and existing.lower() != derived.lower():
                cc_conflicts.append(f"{meg_id_of(hdr)}: file={existing!r} SNPinfo={derived!r}")
            orig_cc_values.append(derived)
        else:
            orig_cc_values.append(existing)

    # --- 9. Write combined single-line FASTA ----------------------------------
    with open(args.output_fasta, 'w') as out:
        for h, s in original_records:
            out.write(f">{h}\n{s}\n")
        for h, s in new_meg_records:
            out.write(f">{h}\n{s}\n")

    # --- 10. Write combined annotation file -----------------------------------
    out_cols = ['header', 'type', 'class', 'mechanism', 'group', 'snp', cc_col]
    with open(args.output_annotations, 'w', newline='') as out:
        writer = csv.writer(out, delimiter=orig_delim, lineterminator='\n')
        writer.writerow(out_cols)
        for r, cc in zip(orig_ann_rows, orig_cc_values):
            writer.writerow([get_ci(r, c) for c in out_cols[:-1]] + [cc])
        for r in new_ann_rows:
            writer.writerow([r[c] for c in out_cols])

    # --- 11. Write updated mappings file --------------------------------------
    map_cols = list(orig_map_fields) if orig_map_fields else \
        ['Database', 'Source_header', 'MEGARes_header', 'UpdatedHeader']
    with open(args.output_mappings, 'w', newline='') as out:
        writer = csv.writer(out, delimiter=map_delim, lineterminator='\n')
        writer.writerow(map_cols)
        for r in orig_map_rows:
            writer.writerow([get_ci(r, c) for c in map_cols])
        for r in new_map_rows:
            writer.writerow([r.get(c, '') for c in map_cols])

    # --- 12. Mapping coverage QC ----------------------------------------------
    all_out_headers = [h for h, _ in original_records] + [h for h, _ in new_meg_records]
    all_out_megs = {meg_id_of(h) for h in all_out_headers}
    mapped_megs = {meg_id_of(get_ci(r, 'MEGARes_header')) for r in orig_map_rows}
    mapped_megs |= {meg_id_of(r['MEGARes_header']) for r in new_map_rows}
    unmapped = sorted(all_out_megs - mapped_megs,
                      key=lambda s: int(s.split('_')[1]) if '_' in s and
                      s.split('_')[1].isdigit() else 0)
    orphan_maps = sorted(mapped_megs - all_out_megs)

    # --- Summary ---------------------------------------------------------------
    print()
    print("-" * 60)
    print(f"  New genes integrated:      {len(new_meg_records)}")
    print(f"    New id range:            "
          f"{args.meg_prefix}{max_id + 1} .. {args.meg_prefix}{next_id}"
          if new_meg_records else "    New id range:            (none)")
    print(f"  Dropped ('remove'):        {removed}")
    if missing_seq:
        print(f"  Skipped (no sequence found in clustering FASTA): {len(missing_seq)}")
        for a in missing_seq[:10]:
            print(f"    - {a}")
        if len(missing_seq) > 10:
            print(f"    ... and {len(missing_seq) - 10} more")
    if incomplete_ann:
        print(f"  Warning: {len(incomplete_ann)} new gene(s) missing one or more of "
              f"Type/Class/Mechanism/Group:")
        for a in incomplete_ann[:10]:
            print(f"    - {a}")
        if len(incomplete_ann) > 10:
            print(f"    ... and {len(incomplete_ann) - 10} more")

    print()
    print(f'  "{cc_col}" for new genes:')
    for k in sorted(cc_counts):
        print(f"    {k:<10} {cc_counts[k]}")
    if snpinfo_ids is not None:
        n_yes = sum(1 for v in orig_cc_values if v.strip().lower() == 'yes')
        n_no = sum(1 for v in orig_cc_values if v.strip().lower() == 'no')
        print(f'  "{cc_col}" for original genes (from SNPinfo): '
              f"Yes={n_yes}  No={n_no}")
        if cc_conflicts:
            print(f"  Warning: {len(cc_conflicts)} original row(s) disagreed with "
                  f"SNPinfo (SNPinfo used):")
            for c in cc_conflicts[:10]:
                print(f"    - {c}")
            if len(cc_conflicts) > 10:
                print(f"    ... and {len(cc_conflicts) - 10} more")

    print()
    print(f"  Mapping rows out:          {len(orig_map_rows) + len(new_map_rows)}"
          f"  (+{len(new_map_rows)} new)")
    if unmapped:
        print(f"  Warning: {len(unmapped)} accession(s) in the database have no "
              f"mapping row:")
        for a in unmapped[:10]:
            print(f"    - {a}")
        if len(unmapped) > 10:
            print(f"    ... and {len(unmapped) - 10} more")
    if orphan_maps:
        print(f"  Warning: {len(orphan_maps)} mapping row(s) point to a MEG id that "
              f"is not in the output database (removed/renumbered accessions):")
        for a in orphan_maps[:10]:
            print(f"    - {a}")
        if len(orphan_maps) > 10:
            print(f"    ... and {len(orphan_maps) - 10} more")

    print()
    print(f"  Total sequences out:       {len(original_records) + len(new_meg_records)}")
    print(f"  Total annotations out:     {len(orig_ann_rows) + len(new_ann_rows)}")
    print()
    print(f"  FASTA:       {args.output_fasta}")
    print(f"  Annotations: {args.output_annotations}")
    print(f"  Mappings:    {args.output_mappings}")


if __name__ == '__main__':
    main()