#!/usr/bin/env python3
"""
Parse a CD-HIT .clstr file and re-assess MEGARes group annotations cluster by cluster.

Usage:
    python parse_clstr.py megares_v4_clust80.clstr [-o cluster_qc.tsv]
        [-a megares_v4_annotations.csv] [-oa megares_v4_annotations_with_clusters.csv]

Outputs
-------
  <output>                 one row per sequence: cluster, representative, similarity,
                           parsed annotation fields, and an Annotation_Status derived
                           from how many distinct groups share the cluster
  <output>.summary.tsv     one row per cluster
  <annotation-output>      one row per MEG accession in MEGARes annotation format,
                           with three cluster columns appended:
                             Cluster_ID            the CD-HIT cluster the gene fell in
                             Cluster_other_groups  every OTHER group in that cluster,
                                                   ';'-delimited (blank when the gene
                                                   is the only group present)
                             Cluster_classes       every class present in that cluster,
                                                   ';'-delimited, including the gene's own

When -a/--annotations is supplied, its columns are preserved verbatim (so extra
fields such as "Can be confirmed" survive) and the cluster columns are appended,
joined on the MEG id. Without it, the annotation file is rebuilt from the .clstr
headers alone.
"""
import argparse
import re
from collections import defaultdict

import pandas as pd

# >MEG_1|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S...|A16S|RequiresSNPConfirmation
MEG_RE = re.compile(r"^MEG_\d+\|")
# 0       843nt, >HEADER... *
# 1       927nt, >HEADER... at +/94.50%
MEMBER_RE = re.compile(
    r"^\s*(?P<idx>\d+)\s+(?P<len>\d+)(?:nt|aa),\s+>(?P<name>.*?)\.\.\.\s*"
    r"(?P<rest>\*|at\s+.*)$"
)
PCT_RE = re.compile(r"([\d.]+)%\s*$")


def parse(path):
    """Yield (cluster_id, header, length, pct_identity, is_representative)."""
    cluster_id = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
                continue
            m = MEMBER_RE.match(line)
            if not m:
                continue
            rest = m.group("rest")
            is_rep = rest.strip() == "*"
            pct = None
            if not is_rep:
                p = PCT_RE.search(rest)
                if p:
                    pct = float(p.group(1))
            yield cluster_id, m.group("name"), int(m.group("len")), pct, is_rep


def split_header(header):
    """Return (type, class, mechanism, group, snp) for MEGARes headers, else Nones."""
    if not MEG_RE.match(header):
        return (None, None, None, None, None)
    f = header.split("|")
    # MEG_id | type | class | mechanism | group [| snp]
    if len(f) < 5:
        return (None, None, None, None, None)
    return (f[1], f[2], f[3], f[4], f[5] if len(f) > 5 else "")


def meg_id_of(header):
    """Return the leading MEG_<n> token of a header, or '' if there is none."""
    if not header or not MEG_RE.match(header):
        return ""
    return header.split("|")[0].strip()


def join_unique(values, sep="; "):
    """Sorted, de-duplicated, delimiter-joined string; '' when empty."""
    return sep.join(sorted({v for v in values if v}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clstr")
    ap.add_argument("-o", "--output", default="cluster_qc.tsv",
                    help="per-sequence QC table (default: cluster_qc.tsv)")
    ap.add_argument("-a", "--annotations", default=None,
                    help="existing MEGARes annotation CSV/TSV to append cluster "
                         "columns to; all of its columns are preserved")
    ap.add_argument("-oa", "--annotation-output",
                    default="annotations_with_clusters.csv",
                    help="annotation file with cluster columns "
                         "(default: annotations_with_clusters.csv)")
    ap.add_argument("--sep", default="; ",
                    help="delimiter for the multi-value cluster columns "
                         "(default: '; ')")
    args = ap.parse_args()

    rows = []
    for cid, header, length, pct, is_rep in parse(args.clstr):
        t, c, mech, grp, snp = split_header(header)
        rows.append(
            dict(
                Cluster_ID=cid,
                Header=header,
                Length=length,
                Similarity=100.0 if is_rep else pct,
                Is_representative=is_rep,
                Type=t,
                Class=c,
                Mechanism=mech,
                Group=grp,
                SNP=snp,
                MEG_ID=meg_id_of(header),
                Is_existing_MEG=bool(MEG_RE.match(header)),
            )
        )
    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("No cluster members parsed - is this a .clstr file?")

    # --- cluster-level verdicts -------------------------------------------------
    groups = defaultdict(set)
    mechs = defaultdict(set)
    classes = defaultdict(set)
    n_annot = defaultdict(int)
    n_new = defaultdict(int)
    for r in rows:
        cid = r["Cluster_ID"]
        if r["Group"]:
            groups[cid].add(r["Group"])
            mechs[cid].add(r["Mechanism"])
            classes[cid].add(r["Class"])
            n_annot[cid] += 1
        else:
            n_new[cid] += 1

    def verdict(cid):
        g, a, n = groups[cid], n_annot[cid], n_new[cid]
        if a == 0:
            return "Needs new annotations"
        if len(g) > 1:
            return "Mismatch in cluster"
        if n > 0:
            return "Inherited annotations"
        return "Existing annotation"

    df["Annotation_Status"] = df.Cluster_ID.map(verdict)
    df["Mismatch_Groups"] = df.Cluster_ID.map(
        lambda c: " ; ".join(sorted(groups[c])) if len(groups[c]) > 1 else ""
    )
    # single unambiguous parent to inherit, when there is one
    df["Inherit_Group"] = df.Cluster_ID.map(
        lambda c: next(iter(groups[c])) if len(groups[c]) == 1 else ""
    )
    df["Inherit_Mechanism"] = df.Cluster_ID.map(
        lambda c: next(iter(mechs[c])) if len(mechs[c]) == 1 else ""
    )
    df.to_csv(args.output, sep="\t", index=False)

    summary = (
        df.groupby("Cluster_ID")
        .agg(
            n_members=("Header", "size"),
            n_existing=("Is_existing_MEG", "sum"),
            min_similarity=("Similarity", "min"),
            min_length=("Length", "min"),
            max_length=("Length", "max"),
            n_groups=("Group", lambda s: s.dropna().nunique()),
            n_mechanisms=("Mechanism", lambda s: s.dropna().nunique()),
            n_classes=("Class", lambda s: s.dropna().nunique()),
            status=("Annotation_Status", "first"),
            groups=("Group", lambda s: " ; ".join(sorted(set(s.dropna())))),
        )
        .reset_index()
    )
    # length spread flags the coverage artefacts (a 525nt SAT beside an 843nt OXA)
    summary["length_ratio"] = (summary.min_length / summary.max_length).round(3)
    summary["length_artefact_suspect"] = summary.length_ratio < 0.8
    summary["crosses_class"] = summary.n_classes > 1

    out2 = args.output.rsplit(".", 1)[0] + ".summary.tsv"
    summary.to_csv(out2, sep="\t", index=False)

    # --- annotation file with cluster columns ---------------------------------
    # Cluster_other_groups excludes the gene's own group (that is what "other"
    # means); Cluster_classes lists every class present, including its own, so a
    # single-class cluster still shows what it is rather than being blank.
    sep = args.sep
    cluster_groups = {c: groups[c] for c in groups}
    cluster_classes = {c: classes[c] for c in classes}

    ann_rows = []
    for r in rows:
        mid = r["MEG_ID"]
        if not mid:
            continue
        cid = r["Cluster_ID"]
        others = cluster_groups.get(cid, set()) - {r["Group"]}
        ann_rows.append(
            {
                "MEG_ID": mid,
                "header": r["Header"],
                "type": r["Type"],
                "class": r["Class"],
                "mechanism": r["Mechanism"],
                "group": r["Group"],
                "snp": r["SNP"],
                "Cluster_ID": cid,
                "Cluster_other_groups": join_unique(others, sep),
                "Cluster_classes": join_unique(cluster_classes.get(cid, set()), sep),
            }
        )
    ann = pd.DataFrame(ann_rows)

    dup = ann.MEG_ID.duplicated().sum()
    if dup:
        print(f"warning: {dup} MEG id(s) appear in more than one cluster; "
              f"keeping the first occurrence of each")
        ann = ann.drop_duplicates("MEG_ID", keep="first")

    cluster_cols = ["Cluster_ID", "Cluster_other_groups", "Cluster_classes"]

    if args.annotations:
        # Preserve the existing annotation file exactly and append the cluster
        # columns, joined on MEG id. Extra columns (e.g. "Can be confirmed")
        # survive untouched.
        sep_in = "\t" if args.annotations.endswith((".tsv", ".tab")) else ","
        base = pd.read_csv(args.annotations, sep=sep_in, dtype=str).fillna("")
        hdr_col = next((c for c in base.columns if c.lower() == "header"), None)
        if hdr_col is None:
            raise SystemExit(
                f"{args.annotations} has no 'header' column to join on"
            )
        base["MEG_ID"] = base[hdr_col].map(meg_id_of)
        merged = base.merge(
            ann[["MEG_ID"] + cluster_cols], on="MEG_ID", how="left"
        )
        missing = merged.Cluster_ID.isna().sum()
        merged[cluster_cols] = merged[cluster_cols].fillna("")
        out_cols = [c for c in base.columns if c != "MEG_ID"] + cluster_cols
        merged[out_cols].to_csv(args.annotation_output, index=False)
        extra = len(ann) - (len(merged) - missing)
        if missing:
            print(f"warning: {missing} annotation row(s) had no cluster "
                  f"(not present in the .clstr); cluster columns left blank")
        if extra > 0:
            print(f"warning: {extra} clustered MEG id(s) are absent from "
                  f"{args.annotations}")
        n_ann_out = len(merged)
    else:
        out_cols = [
            "header", "type", "class", "mechanism", "group", "snp"
        ] + cluster_cols
        ann[out_cols].to_csv(args.annotation_output, index=False)
        n_ann_out = len(ann)

    print(f"{len(df)} sequences in {df.Cluster_ID.nunique()} clusters")
    print(df.Annotation_Status.value_counts().to_string())
    print(f"\nclusters spanning >1 group : {(summary.n_groups > 1).sum()}")
    print(f"clusters spanning >1 class : {summary.crosses_class.sum()}")
    print(f"length-artefact suspects   : {summary.length_artefact_suspect.sum()}")
    n_shared = int((ann.Cluster_other_groups != "").sum())
    print(f"accessions sharing a cluster with another group : {n_shared}")
    print(f"\nwrote {args.output}, {out2} and {args.annotation_output} "
          f"({n_ann_out} rows)")


if __name__ == "__main__":
    main()