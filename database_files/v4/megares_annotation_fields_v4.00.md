# MEGARes 4.0 annotation and provenance files

MEGARes 4.0 is distributed as a sequence file, two annotation files and three provenance files. Every table has a header row and one row per accession (the mapping table also keeps rows for removed accessions).

| File | Format | Used by | Contents |
|---|---|---|---|
| `megares_database_v4.00.fasta` | FASTA | AMR++ (reference database) | Nucleotide sequences, one per accession |
| `megares_annotations_v4.00.csv` | CSV | AMR++ (default annotation file) | 7 columns: accession header, resistance hierarchy, SNP confirmation fields |
| `megares_extended_annotations_v4.00.csv` | CSV | Downstream analysis of AMR++ output | 14 columns: the 7 standard columns, 3 sequence-cluster columns, sequence length and 3 supplementary context columns |
| `megares_external_header_mappings_v4.00.csv` | CSV | Provenance | The source database each sequence was taken from, and its header there |
| `megares_identity_report_v4.00.tsv` | TSV | Provenance | Which current source-database releases hold an identical sequence |
| `megares_ncbi_provenance_v4.00.tsv` | TSV | Provenance | NCBI nucleotide records carrying an identical sequence, with their dates, taxa, replicon types and cited publications |

The three provenance files answer different questions. The mapping table records where a sequence was originally obtained, which is one source per accession. The identity report records which AMR databases hold an identical sequence today, which can be several or none. The NCBI provenance table records where the identical sequence occurs in the NCBI nucleotide collection, and when it was first deposited.

The extended file contains every column of the standard file, so it can be used in its place for any downstream analysis. Its rows are ordered by sequence cluster rather than by accession number. It can also be passed to AMR++ as the annotation file.

## Accession headers

Each accession header has the form:

```
MEG_ID|Type|Class|Mechanism|Group[|RequiresSNPConfirmation]
```

For example:

```
MEG_1|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation
```

Spaces in class and mechanism names are written as underscores in the header. The `class` and `mechanism` columns of the annotation files use spaces.

**Joining files.** The `MEG_ID` (the part of the header before the first `|`) is unique and never changes between releases, whereas the rest of the header changes whenever an accession is re-annotated. Joining on `MEG_ID` is therefore the most robust way to combine these files with each other or with AMR++ output.

## The resistance hierarchy

Every accession belongs to exactly one type, one class, one mechanism and one group. The hierarchy is acyclic, so counts can be summed to any level without counting an alignment twice. Group is the lowest level at which counts are generally recommended; accession-level counts imply more precision than short reads usually support.

## Standard annotation file: `megares_annotations_v4.00.csv`

| Column | Description | Values |
|---|---|---|
| `header` | Full accession header, exactly as it appears in the database FASTA and in AMR++ output. Use this column, or the `MEG_ID` within it, to join annotations to AMR++ results. | e.g. `MEG_1\|Drugs\|Aminoglycosides\|...` |
| `type` | Broadest level: the kind of compound resistance is conferred against. | `Drugs`, `Metals`, `Biocides`, `Multi-compound` |
| `class` | Resistance class. | e.g. `Aminoglycosides`, `betalactams` |
| `mechanism` | Resistance mechanism within a class. | e.g. `Class A betalactamases` |
| `group` | Gene group: accessions representing the same gene or closely related variants of it. | e.g. `CTX`, `SHV`, `TETA` |
| `snp` | Marks accessions for which alignment alone is not evidence of resistance, because resistance is conferred by specific mutations in a gene that is otherwise not a resistance gene (for example *gyrA* or 16S rRNA). With SNP confirmation enabled, AMR++ counts these only when a resistance-conferring variant is found in the reads. | `RequiresSNPConfirmation`; empty for all others |
| `Can_be_confirmed` | For accessions marked `RequiresSNPConfirmation`, whether the variant database holds enough information about their resistance-conferring variants to confirm them computationally. | `Yes`, `No`; empty for all unflagged accessions |

**Accessions that cannot be confirmed.** The variant database was not expanded in release 4.0, so no accession newly flagged in this release can be confirmed. An alignment to an accession whose `Can_be_confirmed` value is `No` is not, on its own, evidence of resistance.

## Extended annotation file: `megares_extended_annotations_v4.00.csv`

Columns 1 to 7 are the same as in the standard annotation file. Seven columns are added: three describe the sequence-similarity structure of the database, one gives the sequence length, and three give supplementary biological context.

### Sequence clusters

| Column | Description | Values |
|---|---|---|
| `Cluster_ID` | Identifier of the sequence cluster to which the accession was assigned. All accessions in the database were clustered together at 80% nucleotide sequence identity using CD-HIT-EST v4.8.1 (see "Clustering parameters" below). Accessions that share a `Cluster_ID` are similar enough that short reads may not reliably tell them apart. | Integer, numbered from the largest cluster (0) to the smallest |
| `Cluster_other_groups` | Gene groups, other than the accession's own, that are represented in the same cluster. | Group names separated by `; `, e.g. `AAK; LEN; OHIO; OKP; OXA`; empty when every member of the cluster belongs to the same group |
| `Cluster_classes` | Resistance classes represented in the cluster. Usually a single class; a cluster whose groups come from more than one class lists each of them. | Class names separated by `; `, e.g. `Aminocoumarins; Fluoroquinolones` |

#### Clustering parameters

```
cd-hit-est -i megares_v4_draft4.fasta -o megares_v4_clust80 \
    -c 0.8 -n 5 -G 0 -aS 0.8 -aL 0.6 -g 1 -sc 1 -d 0 \
    -M 16000 -T 8
```

| Option | Meaning |
|---|---|
| `-c 0.8` | 80% nucleotide sequence identity threshold |
| `-n 5` | Word length 5, the recommended value for thresholds of 0.80 to 0.85 |
| `-G 0` | Local identity: identity is calculated over the aligned region rather than the full length of the shorter sequence |
| `-aS 0.8` | The alignment must cover at least 80% of the shorter sequence |
| `-aL 0.6` | The alignment must cover at least 60% of the longer sequence |
| `-g 1` | Each sequence joins the most similar cluster that meets the threshold, rather than the first one found |
| `-sc 1` | Clusters are numbered by size, so cluster 0 is the largest |
| `-d 0` | Sequence names are written to the cluster file up to the first space; MEGARes headers contain no spaces, so full headers are kept |
| `-M 16000 -T 8` | 16,000 MB of memory and 8 threads (no effect on the result) |

Because `-G 0` and the coverage limits work together, two accessions can share a cluster when they are at least 80% identical over most of the shorter sequence, even if they differ in length.

#### Interpreting clusters

- **Single-group cluster** (`Cluster_other_groups` empty): reads that align to any member support the accession's gene group.
- **Multi-group cluster** (`Cluster_other_groups` not empty): short reads may not separate the groups in the cluster, so counts are more reliably reported at the cluster level than at the group level. The largest cluster, for example, combines the SHV, OKP and LEN beta-lactamases of the *Klebsiella pneumoniae* species complex.
- Cluster membership shows only that short reads may not distinguish two accessions. It does not show that two detections come from the same locus, or that the enzymes they encode are functionally similar.

### Sequence length

| Column | Description | Values |
|---|---|---|
| `Accession_length` | Length of the accession's sequence. | Nucleotides |

### Supplementary context

These three columns give additional information about each accession. They were generated largely automatically, mainly from a search of the NCBI nucleotide collection (see the NCBI provenance table below), and are provided as context for users. They are not a definitive description of the biology of any accession, and they do not affect AMR++ counts.

| Column | Description | Values |
|---|---|---|
| `Origin` | Whether the sequence, as held in MEGARes, appears to be mobile or native to a species, judged from the NCBI records that carry an identical copy. An identical sequence found in two or more genera, or mostly on plasmids, has moved between hosts and is called `acquired`. A sequence found on chromosomes of a single genus, with no mobile-element gene consistently nearby, is called `intrinsic`. Conflicting evidence (a mix of plasmid and chromosomal records, a mobile element near a chromosomal copy, or disagreement between the accession's cluster and group) gives `mixed`. Accessions with too few records took the label agreed by most of their sequence cluster or gene group where there was one, and are otherwise `undetermined`. The label describes the allele MEGARes holds: a chromosomal gene of one species that has since been mobilized (*bla*CMY-2 from *Citrobacter ampC*, for example) is `acquired` as that mobile allele. | `acquired`, `intrinsic`, `mixed`, `undetermined` |
| `Gene_role` | Whether the accession is a resistance determinant itself or a gene that regulates resistance. | `determinant`, `regulator` (e.g. *vanR*/*vanS*, *tetR*, *acrR*, *mexR*) |
| `SNP_requirement_evidence` | For accessions marked `RequiresSNPConfirmation`, whether the label is documented for the accession or inferred, for example from its sequence cluster or gene group. | `documented`, `inferred`; empty for all unflagged accessions |

**Scope.** MEGARes and AMR++ are designed to facilitate compositional analysis of the resistome in metagenomic data across all types of resistance elements. These labels can help interpret results, but partitioning the resistome into separate compartments, such as acquired and intrinsic determinants, would require additional analysis beyond them.

## Cross-database mapping table: `megares_external_header_mappings_v4.00.csv`

Comma-separated. Records, for each accession, the source database from which its sequence was taken and the header the sequence carried there, so every accession can be traced back to its original record. Rows are also kept for accessions removed from the database, marked in `UpdatedHeader`.

| Column | Description | Values |
|---|---|---|
| `Database` | Source database from which the sequence was taken. Labels reflect the release in which the accession was added, so the same database can appear under more than one label. | `CARD`; `AMRFinder` and `NCBI_AMRFinderPlus` (the NCBI Reference Gene Catalog); `NCBI` (legacy GenBank records); `BacMet`; `ResFinder`; `ARG-ANNOT`; `PointFinder`; `Lahey`; `NotFoundDuringUpdate` (legacy accessions whose source record could not be located) |
| `Source_header` | The sequence's header in the source database, as distributed by that database. It usually contains a GenBank or RefSeq accession and, for most sources, the coordinates of the gene in that record. BacMet headers give a protein accession. Legacy `NCBI` rows may list more than one record, separated by spaces. | Format depends on the source, e.g. CARD `CARD\|gb\|ACCESSION\|+\|start-stop\|ARO:...`; AMRFinder `...\|gene\|product NG_ACCESSION:start-stop`; ResFinder `gene_N_ACCESSION`; ARG-ANNOT `(Class)gene:ACCESSION:start-stop:length`; legacy NCBI `gi\|N\|gb\|ACCESSION\|start-stop`; empty for `NotFoundDuringUpdate` rows |
| `MEGARes_header` | The accession's MEGARes header when it was added to the database. | Accession header |
| `UpdatedHeader` | The accession's header in release 4.00, where it differs from `MEGARes_header` (for example, after a gene group was reassigned or a class collapsed), or `REMOVED` for an accession no longer in the database. | Accession header; `REMOVED`; empty when the header is unchanged |

To check a sequence against NCBI directly, use the `verified_accession` column of the NCBI provenance table, which gives a record and coordinates confirmed to hold the exact sequence.

## Identity report: `megares_identity_report_v4.00.tsv`

Tab-separated. For each accession, records which source databases contain an identical nucleotide sequence in their current release (CARD v4.0.1, ResFinder v4.7.2, PointFinder v4.7.1, NCBI AMRFinderPlus v4.2.0). An accession can be identical to sequences in several databases, or in none: for example, a legacy accession whose sequence the source database has since revised.

| Column | Description | Values |
|---|---|---|
| `MEGARes_Accession` | Accession header. | Accession header |
| `Length` | Sequence length. | Nucleotides |
| `CARD` | The header of each CARD sequence identical to the accession. | Headers separated by ` ; `; empty when there is no identical sequence |
| `ResFinder` | As above, for ResFinder. | As above |
| `PointFinder` | As above, for PointFinder. The same header can be listed twice when PointFinder distributes the sequence twice. | As above |
| `NCBI_AMRFinderPlus` | As above, for the NCBI Reference Gene Catalog distributed with AMRFinderPlus. | As above |
| `BacMet` | As above, for BacMet. Empty for all accessions in 4.00, because BacMet distributes protein sequences and was not updated. | Empty |
| `Other/Unknown` | Reserved for identical sequences in other sources. Empty for all accessions in 4.00. | Empty |
| `Total_Matches` | The number of identical source-database entries across all columns. An accession identical to two CARD entries and one ResFinder entry has 3. | Integer; `0` when no source database holds an identical sequence |

## NCBI provenance table: `megares_ncbi_provenance_v4.00.tsv`

Tab-separated, one row per accession. For every accession, the NCBI nucleotide records that carry an identical sequence were identified, counted and summarized. Each sequence was searched against the NCBI `core_nt` database (built October 3, 2026) with BLASTn (megablast), keeping only full-length, gap-free, 100% identical matches on either strand. Record metadata (creation date, organism, replicon type, host, isolation source, country) and cited PubMed articles were then retrieved through NCBI E-utilities.

Three limits apply to every column. `core_nt` holds GenBank, ENA, DDBJ and RefSeq records but not WGS draft contigs, so counts describe that collection on the day it was built rather than every genome sequenced with the gene. The creation date is when the record became public, not when the isolate was collected. PubMed IDs are the articles cited in a record, and a genome record usually cites the genome paper, so an article is a lead to check rather than proof that it first described the gene.

**Accession**

| Column | Description |
|---|---|
| `meg_id`, `header`, `length` | MEGARes ID, accession header and sequence length (nucleotides). |
| `type`, `class`, `mechanism`, `group` | Resistance hierarchy, as in the annotation files. |
| `requires_snp` | `1` if the accession is marked `RequiresSNPConfirmation`, otherwise `0`. |

**Identical NCBI records**

| Column | Description |
|---|---|
| `n_identical_accessions` | Distinct NCBI records with a full-length, 100% identical match. |
| `n_identical_insdc`, `n_identical_refseq` | Of those, GenBank/ENA/DDBJ records and RefSeq records. RefSeq records are later copies of INSDC submissions. |
| `n_identical_loci` | Matches counted per locus; a record can carry the gene more than once. |
| `count_is_lower_bound` | `1` if the search reached its limit on matching records, so the counts are minimums. `0` for every accession in 4.00. |
| `oldest_accession`, `oldest_date`, `oldest_organism`, `oldest_title` | The matching record with the earliest NCBI creation date (INSDC records preferred), with its organism and title. |
| `oldest_is_vector_or_synthetic` | `1` if that record is a cloning vector or synthetic construct, which often carry common resistance markers and are often their earliest records. |
| `n_accessions_with_date` | Matching records with a creation date. |
| `earliest_accessions` | The first ten matching records in date order, as `ACCESSION (date)`. |

**Source-database provenance**

| Column | Description |
|---|---|
| `recorded_source_dbs` | The source database recorded for the accession in the mapping table. |
| `identical_in_source_dbs`, `n_source_dbs_identical` | The source databases holding an identical sequence today (from the identity report), separated by `;`, and their number. |
| `verified_accession` | A source record and coordinates confirmed to hold the exact sequence, as `ACCESSION:start-stop` (`c` before the coordinates marks the complementary strand). Empty when none was confirmed. |
| `verify_status` | Result of checking the source record: `exact_at_coords` (exact match at the recorded coordinates), `exact_in_record` (found elsewhere in the recorded record), `mismatch` (same length, different sequence), `coord_length_differs` or `length_differs` (the recorded coordinates give a different length), `needs_blast` (no coordinates recorded, so only the sequence search was used), `fetch_failed` (the record could not be retrieved), `not_run` (no nucleotide accession in the source header). |

**Taxa and context** (INSDC records only, so RefSeq copies are not counted twice)

| Column | Description |
|---|---|
| `n_artificial_or_vector` | Matching records from vectors or synthetic constructs. |
| `n_species`, `n_genera`, `n_families` | Distinct taxa among the matching records. |
| `top_genera` | Up to five most frequent genera, as `Genus (records)`. |
| `n_plasmid`, `n_chromosome` | Matching records declared as a plasmid or a chromosome. |
| `top_hosts`, `top_isolation_sources`, `top_countries` | Most frequent values of these record fields, as `value (records)`. |
| `earliest_collection_date` | Earliest isolate collection date among the matching records. |

**Cited publications**

| Column | Description |
|---|---|
| `oldest_pmids`, `oldest_citations` | PubMed IDs cited by the oldest record, and `First author (year) Journal` for each. |
| `n_accessions_pmid_checked`, `n_accessions_with_pmid` | Matching records checked for citations, and those citing at least one article. |
| `n_distinct_pmids` | Distinct articles cited across the matching records. |
| `earliest_pub_pmid`, `earliest_pub_date`, `earliest_pub_citation`, `earliest_pub_accessions` | The earliest-published article cited by any matching record, and the records citing it. |

**Gene-level evidence**

A record no longer than the gene plus 2 kb is a gene-level record (a submission of the gene itself); one up to 20 kb is a locus-level record (an operon, cassette, integron or transposon); anything longer is a plasmid or genome.

| Column | Description |
|---|---|
| `n_gene_level_records`, `n_locus_level_records` | Matching records of each size. |
| `oldest_gene_level_accession`, `oldest_gene_level_date`, `oldest_gene_level_title` | The oldest gene-level record. |
| `candidate_pub_pmid`, `candidate_pub_date`, `candidate_pub_citation`, `candidate_pub_accession` | The article most likely to describe the gene: the earliest-published one with gene-level evidence, and the record citing it. |
| `candidate_pub_evidence` | The evidence for that article: `feature_citation` (the annotated gene feature itself cites it), `gene_record` (cited by a gene-level record), `locus_record` (cited by a locus-level record). Empty when no article had such evidence. |
| `n_candidate_pmids` | Articles with gene- or locus-level evidence. |

**Annotation at the match**

For a sample of matching records per accession (up to the ten earliest and ten shortest), the annotation at the matched coordinates was retrieved.

| Column | Description |
|---|---|
| `n_loci_annotation_checked` | Matched loci whose annotation was retrieved. |
| `n_loci_annotated_exact` | Of those, loci where the record annotates a gene, CDS or RNA feature with exactly the same start and end. |
| `annotation_summary` | Tally of what the records annotate at the match: `exact`, `covers` (a feature that extends beyond the match), `overlap` (a feature covering only part of it) or `unannotated`, e.g. `CDS exact (12); rRNA covers (3)`. |
| `annotated_gene_names`, `annotated_products` | Most frequent gene names and product descriptions the records give the matched feature. |

**Other**

| Column | Description |
|---|---|
| `Origin`, `Gene_role` | [Empty in this file. The labels are in the extended annotation file; remove these columns or fill them before release.] |
| `Notes` | Curator notes, e.g. where a gene was re-extracted from different coordinates. |

## Using the annotations with AMR++ output

Join the annotations to an AMR++ count matrix on the accession header (or its `MEG_ID`), then sum counts at the level you want to report. A minimal example in Python, assuming the count matrix has one row per accession header and one column per sample [confirm the AMR++ output file name and layout]:

```python
import pandas as pd

counts = pd.read_csv("AMR_analytic_matrix.csv", index_col=0)  # rows: accession headers
ann = pd.read_csv("megares_extended_annotations_v4.00.csv", index_col="header")
ann = ann.loc[counts.index]  # annotations in the same row order as the counts

# Counts by gene group and by sequence cluster
by_group = counts.groupby(ann["group"]).sum()
by_cluster = counts.groupby(ann["Cluster_ID"]).sum()

# Leave out accessions that require SNP confirmation
not_snp = counts[ann["snp"].isna()]
```

Empty cells are read as missing values (`NaN`) by pandas and as `NA` by R.

## Version

These descriptions apply to release 4.00. If a file name changes, update it here and in the Data Availability statement. Changes between releases are recorded in `CHANGELOG.md`.