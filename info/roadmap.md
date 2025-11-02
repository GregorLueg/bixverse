# *bixverse* roadmap

*Last update: 01.11.2025* </br>

The idea is to over time implement more and more functionalities that are deemed
useful and asked by the users of this package. We wish to avoid implementing
niche algorithms that have not been battle-tested, but provide simple and fast
core function/methods that are broadly used.

## Stability of the package

- ~~Include proper tests via [tinytest](https://github.com/markvanderloo/tinytest/tree/master).~~ 
(Test coverage is slowly but surely increasing. Bugs you find are features ;P )

### Documentation of the package

- ~~Various vignettes for workflows and (more) detailed documentation in terms of 
method implementation.~~ Missing ones:
    * Usage of the ICA class.
    * Correlation-based co-expression module detection.
    * How to run the DGRDL method.

## General methods

- ~~Gene set enrichment analysis on top of continuous data based on 
[blitzGSEA](https://academic.oup.com/bioinformatics/article/38/8/2356/6526383) 
and/or [fgsea](https://www.biorxiv.org/content/10.1101/060012v3).~~ (The fgsea
multi-level and simple method have been implemented in Rust.)
- ~~Gene ontology enrichment with elimination/pruning for continuous values.~~ 
- ~~Gene set variation analysis in [Rust](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)~~
- Add mouse gene ontology data for gene ontology elimination enrichment on top
of mouse data.
- ~~A CisTarget implementation for identification of upstream TFs.~~

#### Biomedical ontologies

- ~~Semantic similarities for ontologies.~~
- ~~Wang similarities for ontologies.~~
- ~~Simplification methods for gene ontology results with ancestry~~

#### Graph methods

- ~~The workflow to diffuse genetic information through networks, see 
[Barrio-Hernandez et al.](https://www.nature.com/articles/s41588-023-01327-9)
and identify 'privileged' sub communities in the graph.~~
- ~~A version of tied diffusion, see 
[Paull et al.](https://academic.oup.com/bioinformatics/article/29/21/2757/195824). 
(Also with implementation for 'privileged' community detection.)~~
- ~~Rapid permutation-based methods for both to generate topology-aware 
background diffusion scores for statistical testing.~~
- ~~Constrained PageRank implementation for heterogenous graphs with sink nodes
and/or edges, see [Ruiz et al.](https://www.nature.com/articles/s41467-021-21770-8).~~
- ~~Similarity network fusion, see [Wang, et al.](https://www.nature.com/articles/nmeth.2810)~~

#### bulkRNAseq methods

- ~~Wrapper class/functions for differential gene expression analyses with 
[limma/voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)
to generate structured objects.~~
- ~~Semantic similarities for ontologies.~~
- ~~Gene ontology enrichment with elimination/pruning for continuous values.~~ 
(First version implemented with fgsea simple.)
- ~~Gene set variation analysis in 
[Rust](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)~~
- ~~Mitch-like multi contrast analysis in 
[Rust](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06856-9)~~

## Rust

- ~~Make usage of borrowing/lifetimes where possible to avoid expensive 
in-memory copying.~~
- ~~Write tests for key functions within Rust.~~ Key functions are also tested
in Rust now.

## Gene module detection

- NMF implementations as a different way to do matrix factorisations. TBD in 
terms of algorithm.
- ~~TOM (topological overlap measure) for correlation-based methods.~~
- ~~Eigengene calculations for correlation-based methods (especially the single
correlation based one), akin to 
[WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)~~
- ~~Hierarchical clustering-based gene module detection on top of correlation-based
methods, inspired from [Srivastava et al.](https://www.nature.com/articles/s41467-018-06008-4).~~
- ~~Implement [sparse dictionary learning](https://pubmed.ncbi.nlm.nih.gov/35085500/)
as a bi-clustering type of version.~~
- ~~Versions of fastICA in Rust with stabilised ICA and RBH graphs based on
correlations (and set similarities), see 
[Cantini et al.](https://academic.oup.com/bioinformatics/article/35/21/4307/5426054)~~
- ~~Correlation and differential correlation based methods with subsequent 
sparsification/pruning of the correlation-based graph and community detection.~~
- NMF implementations as a different way to do matrix factorisations. TBD in 
terms of algorithm.
- Interpretation layers on top of gene modules, i.e., upstream regulators, 
annotations of modules (wrappers over pathway enrichment function), etc.
- ~~[Reciprocal best hit graphs based on correlation](https://academic.oup.com/bioinformatics/article/35/21/4307/5426054) 
structure of gene loadings for matrix factorisation-based co-expression module detection.~~
- ~~Implement [sparse dictionary learning](https://pubmed.ncbi.nlm.nih.gov/35085500/)
as a bi-clustering type of version.~~

## Single cell class and methods 

This is a mid/longer term project to leverage Rust to keep the (sparse) count 
matrices on disk (likely via some serialised binary format) and leverage Rust 
for fast on-disk computations of the data and avoid loading unnecessary data 
into memory where avoidable for single cell and spatial datasets. Think 
[BPCell](https://bnprks.github.io/BPCells/index.html) with Rust instead of C++. 
This would allow analyses of much larger datasets on local infrastructure. Core
functionality to be implemented:

- On disk normalisation and rapid retrieval of count data via CSC and CSR format
- On disk HVG detection
- On disk scaling and PCA for dimension reduction
- Rapid kNN (sNN?) graph generation for community detection and 2D visualisations
- Wrappers to quickly do (pseudo-bulked) DGE via appropriate methods and Wilcox-
based DGEs.
- Wrappers to provide guidelines on co-expression module detection and upstream
regulator prediction (e.g. [hdWGCNA](https://smorabit.github.io/hdWGCNA)).
