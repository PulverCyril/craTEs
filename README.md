# craTEs
Estimate the cis-regulatory activity of transposable element (TEs) subfamilies

# Installation
Depends on the following Bioconductor packages: `Biobase`, `GenomicFeatures` and `GenomicRanges` that must be installed before `craTEs`.

```
install.packages("devtools", quietly = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c("Biobase", "GenomicFeatures", "GenomicRanges"))
library("devtools")
devtools::install_github("bopekno/craTEs")
library("craTEs")
```

# Publications
[Pulver et al., Statistical learning quantifies transposable element-mediated cis-regulation](https://www.biorxiv.org/content/10.1101/2022.09.23.509180v1)

[Pontis et al., Primate-specific transposable elements shape transcriptional networks during human development, Nature Communications 2022](https://www.nature.com/articles/s41467-022-34800-w)

[Martins et al., KRAB zinc finger proteins ZNF587/ZNF417 protect lymphoma cells from replicative stress-induced inflammation](https://www.biorxiv.org/content/10.1101/2023.03.08.531722v1)
