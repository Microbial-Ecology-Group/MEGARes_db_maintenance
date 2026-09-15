# MEGARes Update Workflow

How to produce a new MEGARes release: aggregate sequences from the public AMR databases, identify what is novel, propagate annotations where they can be inherited, curate the remainder by hand, and integrate the result into a new version of the database.

This document covers the full process. For the database itself and the repository layout, see [README.md](README.md). For what changed in each release, see [CHANGELOG.md](CHANGELOG.md).

## Contents

- [Requirements](#requirements)
- [Workflow at a glance](#workflow-at-a-glance)
- [Stage 1: aggregate, cluster and propose annotations](#stage-1-aggregate-cluster-and-propose-annotations)
- [Stage 2: manual curation](#stage-2-manual-curation)
- [Stage 3: integrate into MEGARes](#stage-3-integrate-into-megares)
- [Stage 4: cluster QC](#stage-4-cluster-qc)
- [Stage 5: cross-database identity report](#stage-5-cross-database-identity-report)
- [Finalizing a release](#finalizing-a-release)
- [Pipeline internals](#pipeline-internals)

---

## Requirements

- **Python 3.7+**
- **pandas**, the only third-party Python package required
- **CD-HIT** (`cd-hit-est`)

```bash
conda create -n megares_update python=3.11 pandas cd-hit
conda activate megares_update
```

On an HPC system with modules:

```bash
module load GCC/13.2.0 CD-HIT/4.8.1
```

Stage 1 downloads databases over the network. Where compute nodes have no outbound access, run the download on a login node with `--download-only` and the remaining steps under the scheduler.

### Set up a working directory

Each update gets its own directory under `update_records/`, named `<YYYYMonth>_megares_<version>_update`:

```bash
mkdir -p update_records/2028January_megares_v5_update
cd update_records/2028January_megares_v5_update
```

All commands below are run from inside that directory, with scripts referenced as `../../update_scripts/`.

---

## Workflow at a glance

| Stage | Script | Key output |
|---|---|---|
| 1 | `MEGARes_update_pipeline.py` | `output/megares-update-clustering_parsed.tsv` |
| 2 | manual curation | `draft_annotations_<date>.csv` |
| 3 | `integrate_megares_accessions.py` | new FASTA, annotations, header mappings |
| 4 | `cd-hit-est` + `parse_clstr_qc.py` | annotations with cluster columns, QC tables |
| 5 | `megares_identity_report.py` | `megares_identity_report.tsv` |

Stage 2 is the only step that cannot be automated, and it is where most of the effort goes.

---

## Stage 1: aggregate, cluster and propose annotations

### Check for new source releases first

```bash
python3 ../../update_scripts/MEGARes_update_pipeline.py --check-versions
```

This queries CARD's download redirect, the Bitbucket API for ResFinder and PointFinder tags, and NCBI's FTP, comparing against the version constants near the top of `MEGARes_update_pipeline.py`. It reports which are out of date and which variables to change.

To apply an update, edit the constants and delete the relevant `input/` subdirectory to force a re-download.

### Download

```bash
python3 ../../update_scripts/MEGARes_update_pipeline.py --download-only -i ./input
```

Downloads are cached in `./input/` and reused on subsequent runs.

### Cluster and parse

```bash
python3 ../../update_scripts/MEGARes_update_pipeline.py --skip-download \
    -i ./input -o ./output --threads 12
```

Or run everything in one pass when the network is available from the compute node:

```bash
python3 ../../update_scripts/MEGARes_update_pipeline.py -i ./input -o ./output -T 12
```

Useful options:

| Option | Effect |
|---|---|
| `--cdhit-path PATH` | Point at a specific `cd-hit-est` binary |
| `--identity 0.9` | Change the clustering threshold from the default 0.8 |
| `--no-dedup-report` | Skip the MEGARes redundancy report (stage 2b) |
| `--parse-only --cluster-file FILE` | Re-parse an existing `.clstr` without re-clustering |
| `--cluster-only` | Stop after clustering |

### What stage 1 produces

```
input/
├── all.fsa           Concatenated FASTA, MEGARes first
├── dedup.fsa         After 100% deduplication
├── dedup.fsa.clstr
├── source_map.tsv    Header -> source database
└── <database>/       One directory per downloaded source
output/
├── megares_dedup_report.tsv              One row per MEGARes gene, with 100% matches
├── megares-update-clustering.fsa
├── megares-update-clustering.fsa.clstr
└── megares-update-clustering_parsed.tsv  THE FILE FOR STAGE 2
```

`megares-update-clustering_parsed.tsv` has one row per non-MEGARes gene, with an `Annotation_Status`:

| Status | Meaning | Action |
|---|---|---|
| `Inherited annotations` | Every MEGARes gene in the cluster agrees | Accept |
| `Confirm SNP difference` | They agree except on the SNP field | Review the SNP field only |
| `Mismatch in cluster` | They disagree on Type, Class, Mechanism or Group | Curate; see `Mismatch_Groups` |
| `Needs new annotations` | No MEGARes gene in the cluster | Annotate de novo |

---

## Stage 2: manual curation

Work through `megares-update-clustering_parsed.tsv` and produce a curated CSV. This is expert work: deciding whether a novel sequence belongs in an existing gene group, warrants a new one, or should be excluded.

The curated file must contain:

| Column | Purpose |
|---|---|
| `Gene_Accession` | Matches the accession in the parsed TSV |
| `Source_Database` | Which database the sequence came from |
| `ED_confirmation` | Curator's verdict; rows marked `remove` are dropped at stage 3 |
| `Type`, `Class`, `Mechanism`, `Group`, `SNP` | The assigned annotation |

Name it by date, for example `draft_annotations_20260818.csv`, and keep every draft. The difference between successive drafts is the record of what curation changed, and it is worth having when a reviewer asks why an accession is annotated as it is.

Two things worth checking while curating. `megares_dedup_report.tsv` shows which external genes are already 100% identical to an existing MEGARes gene, so those need no new accession. And the `Mismatch_Groups` column names the conflicting annotations directly, which usually makes the correct call obvious.

---

## Stage 3: integrate into MEGARes

`integrate_megares_accessions.py` assigns MEG_ accessions to the curated sequences, continuing the existing numbering, and writes the new database files.

```bash
python3 ../../update_scripts/integrate_megares_accessions.py \
    -f  megares_database_v3.00_updated.fasta \
    -a  megares_annotations_v3.00_updated.csv \
    -c  input/dedup.fsa \
    -n  draft_annotations_20260818.csv \
    -m  megares_to_external_header_mappings_v3.00_updated.csv \
    -s  /path/to/AmrPlusPlus_SNP/data/SNPinfo.fasta \
    -of megares_v4_draft4.fasta \
    -oa megares_v4_draft4_annotations.csv \
    -om megares_v4_draft4_to_external_header_mappings.csv
```

| Flag | Input |
|---|---|
| `-f` / `--original-fasta` | Previous release FASTA |
| `-a` / `--original-annotations` | Previous release annotations |
| `-c` / `--clustering-fasta` | `input/dedup.fsa`, holding the new sequences |
| `-n` / `--new-annotations` | The curated file from stage 2 |
| `-m` / `--original-mappings` | Previous external-header mapping table |
| `-s` / `--snpinfo-fasta` | Optional; see below |
| `-of`, `-oa`, `-om` | Output FASTA, annotations, mappings |

### What it does

It finds the highest existing `MEG_<n>`, drops curated rows marked `remove`, pulls each surviving sequence from the clustering FASTA, and assigns new ids sequentially. Headers are built as `MEG_<n>|Type|Class|Mechanism|Group[|SNP]` with spaces converted to underscores; the annotation file keeps the spaces. The trailing `|SNP` is omitted when the SNP field is blank.

The mapping output contains every original row verbatim plus one row per newly integrated gene, recording its source database, its original accession, and its new MEG_ header.

### The `-s` flag is worth using

`--snpinfo-fasta` points at `AmrPlusPlus_SNP/data/SNPinfo.fasta` and makes that file **authoritative** for the `Can be confirmed` column, reporting any disagreement with the annotation file.

Without it, new accessions carrying the SNP label default to `No`, on the reasonable grounds that no SNPinfo entry exists for them yet, and existing rows pass their value through unchanged. Supplying `-s` instead decides the value for every row from actual membership in the variant database, which is what makes the "how many SNP-confirmation genes can actually be confirmed" figure defensible. Read the reported disagreements rather than ignoring them: each one is either a stale annotation or a missing SNPinfo entry.

### Iterating

Expect several drafts. Re-run stage 3 with each new curated file, incrementing the draft number in the output names, and keep them all:

```bash
# an earlier draft, before mappings and SNPinfo were wired in
python3 ../../update_scripts/integrate_megares_accessions.py \
    -f megares_database_v3.00_updated.fasta \
    -a megares_annotations_v3.00_updated.csv \
    -c input/dedup.fsa \
    -n draft_annotations_20260817.csv \
    -of megares_v4_draft3.fasta \
    -oa megares_v4_draft3_annotations.csv
```

---

## Stage 4: cluster QC

Cluster the candidate release against itself to find annotation problems the per-cluster inheritance could not catch: gene groups that now split across clusters, clusters spanning more than one class, and length artefacts.

```bash
module load GCC/13.2.0 CD-HIT/4.8.1

cd-hit-est -i megares_v4_draft4.fasta -o megares_v4_clust80 \
    -c 0.8 -n 5 -G 0 -aS 0.8 -aL 0.6 -g 1 -sc 1 -d 0 \
    -M 16000 -T 8
```

Then parse the result, passing the annotation file so the cluster columns are appended to it:

```bash
python3 ../../update_scripts/parse_clstr_qc.py megares_v4_clust80.clstr \
    -a  megares_v4_draft4_annotations.csv \
    -oa megares_v4_draft4_annotations_with_clusters.csv \
    -o  v4_cluster_qc.tsv
```

Running it without `-a` rebuilds the annotation file from the cluster headers alone, which is useful for a quick look but discards any extra columns:

```bash
python3 ../../update_scripts/parse_clstr_qc.py megares_v4_clust80.clstr -o v4_cluster_qc.tsv
```

### Outputs

| File | Contents |
|---|---|
| `v4_cluster_qc.tsv` | One row per sequence: cluster, representative, similarity, parsed annotation, status |
| `v4_cluster_qc.summary.tsv` | One row per cluster, with `crosses_class` and `length_artefact_suspect` flags |
| `megares_v4_draft4_annotations_with_clusters.csv` | The annotation file with three columns appended |

The three appended columns are `Cluster_ID`, `Cluster_other_groups` (every *other* group sharing the cluster, blank when the gene is alone) and `Cluster_classes` (every class present, including the gene's own). Together they make the redundancy structure of the release explicit, which is what lets a downstream user see which accessions will compete for the same reads during alignment.

### What to look at

The console output reports clusters spanning more than one group, clusters spanning more than one class, and length-artefact suspects. Cross-class clusters are the highest-value check: a cluster containing both a betalactamase and a macrolide esterase usually means one of them is misannotated. Several of the v4.00 reclassifications listed in the CHANGELOG were found this way.

`length_artefact_suspect` flags clusters where the shortest member is under 80% the length of the longest, which catches cases where a short fragment clustered with a full-length gene through partial coverage.

### Optional: exact-duplicate check

To confirm no two accessions are identical, cluster at 100% instead:

```bash
cd-hit-est -i megares_v4_draft4.fasta -o megares_v4_derep \
    -G 0 -c 1.0 -AS 0 -AL 0 -AI 1.0 -aS 1.0 \
    -d 0 -M 16000 -T 8
```

---

## Stage 5: cross-database identity report

Reports which MEGARes genes have identical sequences in each source database, by clustering the candidate release against the downloaded databases at 100% identity.

```bash
python3 ../../update_scripts/megares_identity_report.py \
    --megares megares_v4_draft4.fasta \
    --input-dir ./input \
    -o megares_identity_report.tsv \
    -T 6
```

`--input-dir` auto-discovers the database files downloaded in stage 1. Individual files can be given instead:

```bash
python3 ../../update_scripts/megares_identity_report.py \
    --megares megares_v4_draft4.fasta \
    --card input/card/nucleotide_fasta_protein_homolog_model.fasta \
    --card input/card/nucleotide_fasta_protein_variant_model.fasta \
    --resfinder input/resfinder/all.fsa \
    --ncbi input/ncbi/AMR_CDS.fa \
    -o megares_identity_report.tsv
```

To re-parse without re-running CD-HIT:

```bash
python3 ../../update_scripts/megares_identity_report.py --parse-only \
    --cluster-file workdir/dedup.fsa.clstr \
    --source-map input/source_map.tsv
```

The output has one row per MEGARes accession and one column per database, listing matching accessions. **BacMet is deliberately excluded**; any BacMet sequences appearing in a reused `.clstr` are reported under `Other/Unknown` rather than being folded silently into a database column.

---

## Finalizing a release

1. Rename the accepted draft to the release name, for example `megares_v4_draft4.fasta` to `megares_database_v4.00.fasta`, and likewise for the annotation and mapping files.
2. Copy them into `database_files/v4/`, leaving the drafts in `update_records/`.
3. Build the representative database if one is being distributed.
4. Update [CHANGELOG.md](CHANGELOG.md) with the new counts and every annotation change, including the basis for each reclassification.
5. Update the version constants at the top of `MEGARes_update_pipeline.py` so the next update starts from this release.
6. Update the source-database version table in [README.md](README.md).

Step 5 is easy to forget and causes the next update to aggregate against the previous release.

---

## Pipeline internals

Background on why stage 1 works the way it does. Not needed to run an update.

### Database sources

| Database | Description |
|---|---|
| MEGARes | Reference AMR database, always first in concatenation |
| CARD | Comprehensive Antibiotic Resistance Database |
| ResFinder | Acquired antimicrobial resistance genes |
| PointFinder | Chromosomal point mutations |
| NCBI AMRFinderPlus | NCBI curated AMR reference genes |
| BacMet | Biocide and metal resistance genes |

The specific versions used for each release are recorded in [README.md](README.md); the constants currently set in `MEGARes_update_pipeline.py` are the defaults for the *next* update and may already be ahead of the last release.

### MEGARes-first ordering

MEGARes sequences are always placed first in the concatenated FASTA. CD-HIT processes sequences in input order and selects the first qualifying sequence as the cluster representative, so whenever a MEGARes gene is present in a cluster it becomes the representative. That preserves its header, which encodes the full annotation, as the reference for the cluster, and is what makes annotation inheritance possible in the parsing step.

### `cd-hit-est`, not `cd-hit-est-2d`

The two-database tool outputs sequences from DB2 with no match in DB1; it does not produce mixed clusters. Running single-database `cd-hit-est` on the combined set puts MEGARes and external genes in the same cluster whenever they share at least 80% similarity, which is exactly what annotation inheritance requires.

### Source tracking

During concatenation the pipeline scans every header and records its source database in `input/source_map.tsv`. Both the dedup report and the parsed output use this lookup for their `Source_Database` column rather than guessing from header patterns.

### CD-HIT parameters

Deduplication, step 2:

```bash
cd-hit-est -i input/all.fsa -o input/dedup.fsa \
  -G 0 -c 1.0 -AS 0 -AL 0 -aL 1.0 -aS 1.0 -d 200 -T 6
```

| Flag | Value | Purpose |
|---|---|---|
| `-G 0` | local | Local rather than global sequence identity |
| `-c 1.0` | 100% | Cluster only truly identical sequences |
| `-AS 0`, `-AL 0` | 0 | No length tolerance |
| `-aL 1.0`, `-aS 1.0` | 100% | Full-length alignment required |
| `-d 200` | 200 | Keep enough header to preserve MEGARes annotations |

Clustering, step 3:

```bash
cd-hit-est -i input/dedup.fsa -o output/megares-update-clustering.fsa \
  -c 0.8 -d 200 -T 6
```

`-d 200` matters in both: CD-HIT truncates headers in the `.clstr` file by default, which would cut off the pipe-delimited annotation the parsing step depends on.

### `parse_cdhit_clusters.py` versus `parse_clstr_qc.py`

Both read a CD-HIT `.clstr` file, but they answer different questions and are not interchangeable.

`parse_cdhit_clusters.py` runs inside stage 1. It asks, for each **non-MEGARes** gene, whether an annotation can be inherited, and emits one row per candidate with an `Annotation_Status`. It can also be run standalone on any `.clstr` file, though without `source_map.tsv` it omits the `Source_Database` column.

`parse_clstr_qc.py` runs at stage 4, after integration. It asks whether the **finished release** is internally consistent, emitting per-sequence and per-cluster QC tables and appending cluster columns to the annotation file.
