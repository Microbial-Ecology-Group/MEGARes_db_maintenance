# Changelog

All notable changes to the MEGARes database are documented here.

---

## [v4.00] - 2026

`megares_database_v4.00.fasta`

### Summary

| | v3.00 | v4.00 | New in v4 |
|---|---|---|---|
| Gene accessions | 8,733 | 11,506 | 2,773 |
| Accessions requiring SNP confirmation | 572 | 628 | 54 |
| Gene groups | 1,448 | 1,609 | 162 |
| Mechanisms | 233 | 275 | 43 |
| Classes | 59 | 63 | 4 |
| Compound types | 5 | 4 | 0 |

The "New in v4" column counts only what the added sequences introduced: for accessions, the number of new sequences; for the hierarchical levels, the number of categories that did not exist in v3.00. Most new accessions were assigned to categories that already existed, so the number of new categories is far smaller than the number of new sequences. **No v3.00 accessions were removed.**

The compound-type count falls from five to four because a single misspelled label (`Multi-compount`) was corrected, not because a category was removed. See "Annotation corrections" below.

Growth is concentrated in the most redundant part of the database: beta-lactamases account for **60.7%** of v4.00 accessions (6,989 of 11,506) and **86.3%** of the sequences added in this release (2,392 of 2,773).

---

### Added

**New resistance classes (4)**

- Ionophores
- Nitrofuran
- Nitroimidazoles
- Organoarsenicals

**New annotation columns**

The annotation file gains four columns relative to v3.00:

- `Can be confirmed` — whether an accession labelled `RequiresSNPConfirmation` has sufficient variant information in the accompanying variant database for its resistance-conferring variants to be verified computationally. Of the 628 accessions carrying the label, 339 can be confirmed and 289 cannot. This distinction was not previously exposed to users, so an accession that required confirmation but could not receive it was indistinguishable from one that could.
- `Cluster_ID` — the cluster each accession was assigned to during dereplication, making the redundancy structure of the database explicit.
- `Cluster_other_groups` — other gene groups represented in the same cluster, identifying accessions that are expected to compete for the same reads during alignment.
- `Cluster_classes` — the resistance classes represented in the same cluster.

**Cross-database accession mapping**

A mapping table linking each MEGARes accession to the accession identifiers of the source databases from which it was derived, supporting comparison of annotations for the same underlying sequence across resources.

**Representative database**

A dereplicated version of the database (2,376 sequences) retaining cluster representatives for each gene group. Distributed alongside the complete database for analyses where alignment ambiguity from near-identical accessions is a concern.

**Semi-automated update script**

A Python utility that retrieves current releases of the source repositories, removes sequences already present in MEGARes, and clusters the remainder against existing accessions. Sequences clustering exclusively with a single gene group have that group's annotation proposed for inheritance; sequences clustering with more than one group, or with none, are flagged for manual annotation. The script partitions each update rather than replacing manual curation.

---

### Changed

BacMet has not been updated since its sequences were incorporated into earlier releases and therefore contributed no new accessions to v4.00.

#### Annotation corrections

Thirty-three accessions carried over from v3.00 received corrected annotations. Twenty-one are spelling, formatting or consolidation fixes that do not change how an accession is binned. Twelve are reclassifications that do.

**Reclassifications**

These move an accession to a different type, class, mechanism or group, and change which category its counts fall into. All were confirmed by BLAST.

| Accession | v3.00 | v4.00 | Basis |
|---|---|---|---|
| MEG_5254 | betalactams / Class D betalactamases / OXA | **MLS / Macrolide esterases / ESTX** | BLAST |
| MEG_5255 | betalactams / Class D betalactamases / OXA | **MLS / Macrolide esterases / ESTX** | BLAST |
| MEG_9979 | betalactams / ACC | **Aminoglycosides / AAC3** | BLAST |
| MEG_8624 | betalactams / Class A betalactamases / PENA | **betalactams / Penicillin binding protein / PBP2** | BLAST |
| MEG_3648 | Multi-drug resistance / Multi-drug ABC efflux pumps / LSAE | **MLS / Lincosamide nucleotidyltransferases / LNUB** | BLAST |
| MEG_3516 | Multi-compound / Drug and biocide SMR efflux pumps / KPNO | **Multi-compound / Drug and biocide mutant porin protein / OMPC**, plus `RequiresSNPConfirmation` | BLAST; clusters with OMPK |
| MEG_10076 | Multi-drug resistance | **Trimethoprim** | |
| MEG_11357, MEG_11358, MEG_11359 | Type: Drugs | Type: **Biocides** | |
| MEG_11344 | Group: OKPSHV | Group: **SHV** | |
| MEG_11390 | Group: EMRD | Group: **EMRDD** | |
| MEG_8674 | | gained `RequiresSNPConfirmation` | |

Two notes on specific cases. MEG_3648 was annotated as *lsa(E)*, an ABC-F protein of roughly 1,479 nt, but is *lnu(B)* at roughly 804 nt; the two are typically adjacent on the same cassette, which is the likely source of the original confusion. MEG_3516 was reclassified as a mutant porin and therefore also gained the `RequiresSNPConfirmation` label, since porin-mediated resistance depends on specific mutations rather than gene presence.

**Spelling corrections**

| Accessions | Correction |
|---|---|
| MEG_10100, MEG_10120, MEG_10121 | `Pyrazinamide_resistant_polyketide_systhase` → `synthase` |
| MEG_10101 | `Multi-drug_resitant_outer_membrane_protein` → `resistant` |
| MEG_10141 | `Sulfonamide_resistant_ECF_transporter_comonent` → `component` |
| MEG_10128 | `Rifampin_monooyxgenase` → `Monooxygenase` |
| MEG_11399 | `Macrolide_eflux_MFS_transporter` → `MLS_resistance_MFS_efflux_pumps` |
| MEG_8700 | Type `Multi-compount` → `Multi-compound` |

**Formatting and capitalization**

- MEG_3075: trailing underscore removed from `Aluminum_ATPase_`
- MEG_9982, MEG_11500, MEG_11505: `nim_nitroimidazole_reductase` → `NIM_nitroimidazole_reductase`
- MEG_10045: `Daptomycin_resistant_mutant` → `Daptomycin-resistant_mutant`
- Trailing whitespace removed from group and mechanism fields throughout
- All mechanism names now begin with a capital letter

**Consolidation of near-duplicate mechanism names**

| Accessions | Consolidation |
|---|---|
| MEG_10066, MEG_10130, MEG_10133, MEG_10135, MEG_10136 | `Pyrazinamide-resistant_gene_mutation` → `Pyrazinamide-resistant_mutant` |
| MEG_10082, MEG_10153 | `Daptomycin_resistant_glycerolphosphodiesterase` and `Daptomycin_resistant_mutant_YvlB` → `Daptomycin-resistant_mutant` |
| MEG_8677, MEG_8678 | `Aminoglycoside-resistant_16S_ribosomal_subunit_protein` → `Aminoglycoside-resistant_30S_ribosomal_protein_S12` |

The last is a correction rather than a consolidation: MEG_8677 and MEG_8678 are *rpsL*-type ribosomal proteins, not 16S rRNA, confirmed by BLAST. Both retain the RPSL group and the `RequiresSNPConfirmation` label.

---

### Notes for users upgrading from v3.00

All v3.00 accessions remain present, so v3.00 identifiers resolve in v4.00. However, twelve accessions were reclassified, so **counts aggregated by class, mechanism or group are not directly comparable between versions** for those determinants. Analyses that compare results generated under different database versions should state which version was used.

Because clustering is performed within gene groups, the seven accessions that changed group also alter the clustering input for their source and destination groups. Users of the representative database should regenerate it rather than mixing representatives derived from different releases.

---

## [v3.00] - 2023

See Bonin et al., *Nucleic Acids Research* 51:D744–D752 (2023).