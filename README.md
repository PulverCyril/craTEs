# *craTEs* (*cis*-regulatory activities of Transposable Element subfamilies)
*craTEs* is a statistical framework that estimates the *cis*-regulatory activity of transposable element (TE) subfamilies from standardly processed RNA-seq data. It is provided as an R package.

## Why should I care?
Are you looking for *cis*-regulatory sequences that may explain your transcriptomic data? Then read on. 

Transposable elements (TEs) are virus-like selfish genetic parasites of eukaryote genomes. More than half of the human genome is TE-derived, owing to the copy-and-paste mechanisms that the majority of human TEs leverage to replicate. To do so, TEs must first attract the host transcriptional machinery and have therefore evolved embedded cis-regulatory sequences to attract host transcription factors and epigenetic modifiers. In short, as TEs spread, they litter genomes with hundreds to thousands copies of ready-to-use *cis*-regulatory platforms. It is now widely believed that waves of TE invasion contribute to enhancer turnover. 

Here is what you will benefit from using *craTEs*:
- Interpretable quantitative estimates of TE-mediated *cis*-regulatory activities from standardly processed RNA-seq counts
- Leveraging replicates to maximize statistical power
- A lightweight method that easily scales to screening thousands of RNA-seq experiments for TE-mediated *cis*-regulation

Here are the pain points you will avoid thanks to *craTEs*: 
- Remapping your RNA-seq data to TE sequences, which are masked in standard mapping pipelines. In addition, we have shown that TE-derived transcription is not necessarily a good proxy for TE-mediated *cis*-regulation.
- Having to generate epigenomics data - e.g. ChIP-seq, ATAC-seq or DNase-seq - to get insights into the role of specific TE subfamilies in gene regulation.

## Why this acronym?
TEs were long discarded as non-functional and uninteresting "junk DNA". *craTEs* thus alludes to the crates found in vinyl record shops, in particular "dollar bins" in second hand shops where "diggers" go looking for forgotten gems.

## Try it!

Kickstart your explorations with *craTEs* using [this repository](https://renkulab.io/projects/cyril.pulver/crates-basics) on the reproducible data science platform `renku`, powered by the Swiss Federal Instituted of Technology.

## Usage

```
library(craTEs)

# genes as rows, samples as columns
countTable = readRDS('path_to_countTable')

# cis-regulatory susceptibility matrix
N = craTEs::N_from_tsv('path_to_N_matrix') #e.g. downloaded from https://doi.org/10.5281/zenodo.8117257

preprocessed = craTEs::preprocess_E_N_for_activities(countTable, N, log_tpm_plot_path = 'path_to_qc_file.pdf')

# estimating TE-dependent cis-regulatory activities from RNA-seq data
activities = craTEs::getSignifActivities(preprocessed$E_centered, preprocessed$N, treatment_group = c("treatment_s1", "treatment_s2"), control_group = c("control_s1", "control_s2"))

# plotting cis-regulatory subfamilies (with example output)
craTEs::plot_signif_subfams(activities, 0.05, 4, "Epigenetic repression of LTR5-Hs/SVA g#1")
```
![image](https://github.com/bopekno/craTEs/assets/44056089/02c2017a-819a-4dd1-a9fd-3706b70d7538)

See [this jupyter notebook](https://renkulab.io/gitlab/crates/klf4-znf611-sva-crispri/-/blob/master/notebooks/TE_subfamily_diff_activity_poc.md) for examples of more complex use cases.


## Installation
Depends on the following Bioconductor packages: `Biobase`, `GenomicFeatures`, `GenomicRanges`, `RMariaDB`, `ensembldb` and `EnsDb.Hsapiens.v86` that must be installed before `craTEs`, as per the following R commands:

```
install.packages("devtools", quietly = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c("Biobase", "GenomicFeatures", "GenomicRanges", "RMariaDB", "ensembldb", "EnsDb.Hsapiens.v86"))
library("devtools")
devtools::install_github("pulvercyril/craTEs")
```

## Publications
[Pulver et al., Statistical learning quantifies transposable element-mediated cis-regulation](https://www.biorxiv.org/content/10.1101/2022.09.23.509180v1)

[Pontis et al., Primate-specific transposable elements shape transcriptional networks during human development, Nature Communications 2022](https://www.nature.com/articles/s41467-022-34800-w)

[Martins et al., KRAB zinc finger proteins ZNF587/ZNF417 protect lymphoma cells from replicative stress-induced inflammation](https://www.biorxiv.org/content/10.1101/2023.03.08.531722v1)

## Attributions
Concept: Cyril Pulver, Didier Trono

Method development: Cyril Pulver, Raphaël de Fondeville, Julien Pontis

Implementation: Cyril Pulver

Funding: Didier Trono
