# RNA-seq Reanalysis: Plant Immunity and Temperature

## Overview

This repository contains an independent reanalysis of publicly available RNA-seq data from the study:

> **Kim et al. (2022)** – *Increasing the resilience of plant immunity to a warming climate*
> Nature (2022) | [https://doi.org/10.1038/s41586-022-04902-y](https://doi.org/10.1038/s41586-022-04902-y)

The purpose of this project is methodological rather than biological novelty. The analysis was carried out as a skills-building exercise to consolidate practical experience with RNA-seq differential expression workflows, data handling, and visualization in **R**.

The results are not intended to reproduce the original paper exactly, but to demonstrate a transparent, end‑to‑end analytical workflow using standard tools in the RNA‑seq ecosystem.

---

## Analysis Scope

The workflow covers:

* Import of gene-level count data
* Differential expression analysis with **DESeq2**
* Exploratory data analysis (PCA, heatmaps)
* Visualization of biologically relevant gene sets
* Gene Ontology (GO) enrichment analysis

Upstream processing steps (read QC, alignment, and gene-level quantification) were performed prior to this analysis using **Galaxy**. This repository focuses on the **downstream statistical and exploratory analysis in R**.

---

## Upstream Galaxy workflow (details)

The upstream RNA-seq processing was carried out in **Galaxy** using a standard, transparent pipeline. The included Galaxy workflow export (`.ga`) documents the exact tools and connections. In brief, the workflow comprised:

* **Quality control of raw reads**: *Falco* (FastQC-compatible summaries)
* **Adapter and quality trimming**: *Trim Galore*
* **Reference preparation / annotation handling**: *gffread*
* **Read alignment**: *HISAT2* against the *Arabidopsis thaliana* reference genome
* **Post-alignment processing**: *samtools* (sorting) and alignment statistics (*flagstat*)
* **Alignment-level QC**: *Qualimap BAM QC*
* **Gene-level quantification**: *featureCounts* to produce count tables used downstream
* **Aggregate QC reporting**: *MultiQC*

Only the resulting **gene-level count matrices** and sample metadata are used for downstream analysis in R.

---

## Experimental Design (Simplified)

Samples are grouped by:

* **Treatment**: mock vs *Pseudomonas syringae* pv. *tomato* DC3000
* **Temperature**: 23 °C vs 30 °C

Differential expression contrasts include:

* Infection effects at each temperature
* Temperature effects within each treatment

The design is intentionally kept simple to highlight core DESeq2 usage rather than complex modeling.

---

## Repository Structure

```
├── featureCounts_collection/   # Gene-level count files (input)
├── tabularmeta.txt             # Sample metadata
├── Plots_Generated/            # All generated figures
│   ├── PCA/
│   ├── heatmaps/
│   ├── boxplots/
│   └── GO/
├── Galaxy-Workflow-Workflow_constructed_from_history__RNA-seq_Analysis_.ga  # Galaxy upstream workflow
├── RNA_Seq.R                   # Main downstream DESeq2 analysis script
├── RNA_Seq.Rproj               # RStudio project file
├── .gitignore
└── README.md
```

├── featureCounts_collection/   # Gene-level count files (input)
├── tabularmeta.txt             # Sample metadata
├── Plots_Generated/            # All generated figures
│   ├── PCA/
│   ├── heatmaps/
│   ├── boxplots/
│   └── GO/
├── RNAseq_DESeq2_analysis.R    # Main analysis script
└── README.md

```

All figures are automatically written to disk to ensure full reproducibility.

---

## Key Methods

- **Differential expression**: DESeq2
- **Normalization / transformation**: Variance stabilizing transformation (VST)
- **Visualization**: ggplot2, pheatmap
- **Annotation**: org.At.tair.db (Arabidopsis thaliana)
- **Functional enrichment**: clusterProfiler (GO Biological Process)

Only genes with minimal expression are retained, and standard multiple‑testing correction is applied.

---

## Notes on Interpretation

- This analysis is a *reanalysis* and *learning exercise*, not a formal reproduction study.
- Some parameter choices (e.g. thresholds, contrasts) may differ from those in the original publication.
- Biological conclusions should be interpreted cautiously and in the context of the original study.

---

## Motivation

This project was undertaken to:
- Strengthen practical skills in RNA-seq data analysis
- Gain familiarity with real, complex experimental designs
- Practice producing publication-quality figures in R

This repository is intended to demonstrate **analytical workflow competence rather than biological discovery**.

---

## Software Environment

Analysis performed in **R** using Bioconductor packages, including:
- DESeq2
- clusterProfiler
- AnnotationDbi
- tidyverse

Exact versions may vary; the script is intended to be readable and adaptable.

---

## Data Availability

Raw sequencing data and experimental details are available via the original publication and associated public repositories. This repository does not host raw FASTQ files.

Upstream processing (quality control, alignment, and gene-level quantification) was performed using **Galaxy**. The Galaxy workflow used to generate the gene-level count tables is included in this repository as a workflow export (`.ga`) for transparency and reproducibility.

---

## Acknowledgements

All credit for experimental design, data generation, and primary biological insights belongs to the original authors:
Kim et al., Nature (2022).

This repository contains only independent downstream analysis code.

```
