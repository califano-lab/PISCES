---
title: "PISCES Tutorial"
authors: 
 - Lukas Vlahos
 - Aleksandar Obradovic
 - Pasquale Laise
 - Andrea Califano
---
**Authors:** Lukas Vlahos, Aleksandar Obradovic, Pasquale Laise, Andrea Califano  
**Contacts:**

* Lukas Vlahos: lv2395@cumc.columbia.edu
* Aleksandar Obradovic: azo2104@cumc.columbia.edu

### Overview

The pipeline for Protein Activity Inference in Single Cells (PISCES) is a regulatory-network-based methdology for the analysis of single cell gene expression profiles.

Currently, the PISCES manuscript is available on bioRxiv: https://doi.org/10.1101/2021.05.20.445002

NOTE: This version of the pipeline is a newer iteration, implementing new algorithms developed by the Califano lab. An updated manuscript and more robust set of vignettes for the newer iteration of the pipeline is forthcoming.

### Installation

Here's how you can install the PISCES package:

```
install.packages(c("cluster", "ggplot2", "devtools", "Seurat", 
                   "pheatmap", "BiocManager", "RColorBrewer"))
BiocManager::install("viper")
BiocManager::install("biomaRt")
devtools::install_github("JEFworks/MUDAN")
devtools::install_github(repo = "califano-lab/PISCES", force = TRUE, build_vignettes = TRUE)
```
You can then learn about how to use PISCES with our vignettes:

```
library(PISCES)
browseVignettes(package = "PISCES")
```

Some other features we're working on right now:
* Vignette demonstating the functionality of MWKMeans for analyzing trajectories
* RCPP ARACNe for easier network generation

=======
### References

1.	Lachmann, A., et al., *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information*. Bioinformatics, 2016. 32(14): p. 2233-5.  
2.	Califano, H.D.a.A., *iterClust: Iterative Clustering*. R package version 1.4.0. 2018: https://github.com/hd2326/iterClust.  
3.	Ding, H., et al., *Quantitative assessment of protein activity in orphan tissues and single cells using the metaVIPER algorithm*. Nat Commun, 2018. 9(1): p. 1471.  
4.  Rosseeuw, P.J., *Journal of Computational and Applied Mathematics* 20 (1987) 53-65  
5.  Izenman, A.J., *Modern Multivariate Statistical Techniques. Regression, Classification, and Manifold Learning*. Springer text in statistics, 2008 (Chapter 12)

#### Acknowledgements

Jeremy Dooley - for his advice and expertise in single cell sequencing experiments.  
Hongxu Ding - whose work in the Califano laid the groundwork for the development of this pipeline.  
Evan Paull - for help with software and tutorial development and testing.

