# craTEs
Estimates the cis-regulatory activity of transposable element (TEs) subfamilies using RNA-seq data.

# Installation
Depends on the following Bioconductor packages: `Biobase`, `GenomicFeatures`, `GenomicRanges` and `RMariaDB` that must be installed before `craTEs`.

```
install.packages("devtools", quietly = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c("Biobase", "GenomicFeatures", "GenomicRanges", "RMariaDB"))
library("devtools")
devtools::install_github("pulvercyril/craTEs")
```
# Usage

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


# Publications
[Pulver et al., Statistical learning quantifies transposable element-mediated cis-regulation](https://www.biorxiv.org/content/10.1101/2022.09.23.509180v1)

[Pontis et al., Primate-specific transposable elements shape transcriptional networks during human development, Nature Communications 2022](https://www.nature.com/articles/s41467-022-34800-w)

[Martins et al., KRAB zinc finger proteins ZNF587/ZNF417 protect lymphoma cells from replicative stress-induced inflammation](https://www.biorxiv.org/content/10.1101/2023.03.08.531722v1)

# Attributions
Concept: Cyril Pulver, Didier Trono

Method development: Cyril Pulver, Raphaël de Fondeville, Julien Pontis

Implementation: Cyril Pulver

Funding: Didier Trono
