# *bixverse* roadmap

*Last update: 19.03.2025* </br>

The idea is to over time implement more and more functionalities that are deemed
useful and asked by the users of this package. We wish to avoid implementing
niche algorithms that have not been battle-tested, but provide simple and fast
core function/methods.

## Usability of the package

- Various vignettes for workflows and (more) detailed documentation in terms of 
method implementation.

## General methods

- Gene set enrichment analysis on top of continuous data based on 
[blitzGSEA](https://academic.oup.com/bioinformatics/article/38/8/2356/6526383) 
and/or [fgsea](https://www.biorxiv.org/content/10.1101/060012v3).
- Wrapper class/functions for differential gene expression analyses with 
[limma/voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)
to generate structured objects.
- Add mouse gene ontology data for gene ontology elimination enrichment on top
of mouse data.

## Gene module detection

- Eigengene calculations for correlation-based methods (especially the single
correlation based one), akin to [WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)
- Hierarchical clustering-based gene module detection on top of correlation-based
methods.
- Interpretation layers on top of gene modules, i.e., upstream regulators, 
annotations of modules, etc.
- [Reciprocal best hit graphs based on correlation](https://academic.oup.com/bioinformatics/article/35/21/4307/5426054) 
structure of gene loadings for matrix factorisation-based co-expression module detection.
- Implement [sparse dictionary learning](https://pubmed.ncbi.nlm.nih.gov/35085500/).

## Single cell class and methods 

This is a mid/longer term project to leverage Rust to keep the count matrices
on disk (via h5?) and leverage Rust for fast on-disk computations of the data 
and avoid loading unnecessary data into memory where avoidable. Think 
[BPCell](https://bnprks.github.io/BPCells/index.html) with Rust instead of C++. 
Core functionality to be implemented:

- On disk normalisation and rapid retrieval of count data via CSC and CSR format
- On disk HVG detection
- On disk scaling and PCA for dimension reduction
- Rapid kNN (sNN?) graph generation for community detection and 2D visualisations
- Wrappers to quickly do (pseudo-bulked) DGE via appropriate methods and Wilcox-
based DGEs.
