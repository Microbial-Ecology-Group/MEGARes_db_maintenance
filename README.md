# MEGARes

MEGARes is a curated database of antimicrobial resistance gene (ARG) sequences, organized under an acyclic hierarchical ontology that allows alignment counts to be aggregated to any level without multiple counting. It is designed for high-throughput analysis of unassembled short-read metagenomic data, and is the reference database used by [AMR++](https://github.com/Microbial-Ecology-Group/AMRplusplus).

Current release: **v4.00**

Database files, annotations and documentation are also available at [meglab.org](https://www.meglab.org).

---

## Repository layout

```
.
├── README.md                      This file
├── CHANGELOG.md                   Version history and all annotation changes
├── megares_update_workflow.md     How to run a database update
│
├── database_files/                Current and previous releases
│   ├── v1/
│   ├── v2/
│   ├── v3/
│   └── v4/
│
├── update_scripts/                Tooling used to build each release
│   ├── MEGARes_update_pipeline.py
│   ├── integrate_megares_accessions.py
│   ├── parse_cdhit_clusters.py
│   ├── parse_clstr_qc.py
│   └── megares_identity_report.py
│
└── update_records/                Working files from each update
    └── 2026July_megares_v4_update/
```

Three documents sit in the main directory. **README.md** describes the database and repository. **[CHANGELOG.md](CHANGELOG.md)** records what changed between releases. **[megares_update_workflow.md](megares_update_workflow.md)** is the step-by-step guide to producing a new release.

### `database_files/`

Each version directory holds the sequence file and its annotation file. Previous releases are retained so that published analyses remain reproducible: results generated against v3.00 cannot be reproduced with v4.00, because some accessions were reannotated between releases.

The v4.00 directory also contains the representative database, a dereplicated version retaining cluster representatives for each gene group, and the cross-database mapping table linking MEGARes accessions to their source-database identifiers.

### `update_scripts/`

| Script | Role |
|---|---|
| `MEGARes_update_pipeline.py` | Downloads the source databases, deduplicates, clusters, and proposes annotations by inheritance |
| `integrate_megares_accessions.py` | Assigns MEG_ accessions to curated new sequences and writes the new FASTA, annotation and mapping files |
| `parse_cdhit_clusters.py` | Standalone parser for a CD-HIT `.clstr` file, used inside the pipeline |
| `parse_clstr_qc.py` | Cluster QC on a finished release, and appends cluster columns to the annotation file |
| `megares_identity_report.py` | Reports which MEGARes genes have identical sequences in each source database |

### `update_records/`

Working files from each update, kept so a release can be audited or reproduced. One subdirectory per update, named `<YYYYMonth>_megares_<version>_update`:

```
update_records/2026July_megares_v4_update/
├── draft_annotations_20260817.csv                  Curated annotations, draft 3
├── draft_annotations_20260818.csv                  Curated annotations, final
├── megares_v4_draft3.fasta
├── megares_v4_draft3_annotations.csv
├── megares_v4_draft4.fasta                         Became v4.00
├── megares_v4_draft4_annotations.csv
├── megares_v4_draft4_annotations_with_clusters.csv
├── megares_v4_draft4_to_external_header_mappings.csv
├── v4_cluster_qc.tsv                               Per-sequence cluster QC
├── v4_cluster_qc.summary.tsv                       Per-cluster summary
└── megares_identity_report.tsv                     Cross-database identity
```

Subsequent updates add a parallel directory, for example `2028January_megares_v5_update`, leaving earlier ones untouched. Draft files are retained deliberately: the difference between successive drafts is the record of what manual curation changed.

---

## Using the database

For resistome analysis, see [AMR++](https://github.com/Microbial-Ecology-Group/AMRplusplus). To use the sequences directly, index the FASTA with your aligner of choice:

```bash
bwa index database_files/v4/megares_database_v4.00.fasta
samtools faidx database_files/v4/megares_database_v4.00.fasta
```

### Header format

```
MEG_<id>|<Type>|<Class>|<Mechanism>|<Group>[|RequiresSNPConfirmation]
```

for example:

```
MEG_3244|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation
```

The five pipe-delimited fields are the ontology levels, so annotation can be parsed from the reference name alone without a separate lookup file. Spaces are replaced by underscores in headers; the annotation file retains the spaces.

The optional sixth field marks accessions whose resistance depends on a specific sequence variant rather than on gene presence. Coverage of these genes is **not** evidence of resistance without variant confirmation. The annotation file's `Can_be_confirmed` column records whether an accession carrying this label has sufficient variant information for computational confirmation.

A detailed description of the columns found in the annotation files can be found [here](database_files/v4/megares_annotation_fields_v4.00.md).

---

## Updating the database

Full instructions are in **[megares_update_workflow.md](megares_update_workflow.md)**. In outline:

| Stage | Script | Produces |
|---|---|---|
| 1. Aggregate and cluster | `MEGARes_update_pipeline.py` | `dedup.fsa`, cluster file, parsed annotation candidates |
| 2. Manual curation | (by hand) | `draft_annotations_<date>.csv` |
| 3. Integrate | `integrate_megares_accessions.py` | New FASTA, annotation file and mapping file with MEG_ accessions assigned |
| 4. Cluster QC | `parse_clstr_qc.py` | Cluster columns appended to annotations, plus QC tables |
| 5. Identity report | `megares_identity_report.py` | Cross-database identity table |

Stage 2 is the part that cannot be automated. The pipeline partitions candidate sequences into those whose annotation can be inherited from an unambiguous cluster and those requiring expert judgement; the second group has to be annotated by a curator.

### Requirements

- **Python 3.7+**
- **pandas** — the only third-party Python package these scripts require
- **CD-HIT** (`cd-hit-est`), via conda, an HPC module, or compiled from source

```bash
conda create -n megares_update python=3.11 pandas cd-hit
conda activate megares_update
```

On an HPC system with modules:

```bash
module load GCC/13.2.0 CD-HIT/4.8.1
```

Everything else is Python standard library. Stage 1 downloads source databases over the network, so it needs outbound access; where compute nodes are offline, run the download separately on a login node with `--download-only`.

### Quick start

```bash
cd update_records/2026July_megares_v4_update

# Stage 1a: download source databases (needs network)
python3 ../../update_scripts/MEGARes_update_pipeline.py --download-only -i ./input

# Stage 1b: deduplicate, cluster and parse (compute-heavy; submit to a scheduler)
python3 ../../update_scripts/MEGARes_update_pipeline.py --skip-download \
    -i ./input -o ./output --threads 12
```

Stage 1 produces `output/megares-update-clustering_parsed.tsv`, one row per non-MEGARes gene with an `Annotation_Status` of `Inherited annotations`, `Confirm SNP difference`, `Mismatch in cluster`, or `Needs new annotations`. The last two require curation; the curated result becomes the `draft_annotations_<date>.csv` input to stage 3.

Before starting, check whether any source database has published a new release:

```bash
python3 ../../update_scripts/MEGARes_update_pipeline.py --check-versions
```

---

## Source databases

MEGARes aggregates sequences from:

| Database | Used in v4.00 |
|---|---|
| CARD | v4.0.0 |
| ResFinder | 4.7.1 |
| PointFinder | 4.7.1 |
| NCBI AMRFinderPlus | 4.0 |
| BacMet | 2.0 |

BacMet has not been updated since its sequences were incorporated into an earlier release and contributed no new accessions to v4.00.

Users redistributing MEGARes or a derivative should check the terms of each source database, since an aggregated resource cannot be licensed more permissively than its most restrictive source allows.

---

## Citation

If you use MEGARes, please cite the release you used. Version history is in [CHANGELOG.md](CHANGELOG.md).

[V4 CITATION TO BE ADDED]

Previous releases:

- **v3.00** — Bonin N, et al. MEGARes and AMR++, v3.0: an updated comprehensive database of antimicrobial resistance determinants and an improved software pipeline for classification using high-throughput sequencing. *Nucleic Acids Research* 51:D744–D752 (2023).
- **v2.00** — Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. *Nucleic Acids Research* 48:D561–D569 (2020).
- **v1.01** — Lakin SM, et al. MEGARes: an antimicrobial resistance database for high throughput sequencing. *Nucleic Acids Research* 45:D574–D580 (2017).

