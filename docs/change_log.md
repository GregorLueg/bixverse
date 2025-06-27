# *bixverse* change log

*Last update: 27.06.2025* </br>

### Version **0.0.1.5**

#### New features

- Rework of the ontology class. Additionally, added Wang similarity as an 
additional measure.

#### Bug fixes, documentation updates

...

#### Breaking changes

- Functions and methods related to the ontology class will be breaking.

### Version **0.0.1.4**

#### New features

- Implementation of TOM in Rust for usage for correlation-based module detection.
- CoReMo-based gene module detection has been added.

#### Bug fixes, documentation updates

...

#### Breaking changes

- Renaming of `cor_module_check_res()` to `cor_module_graph_check_res()` and
`cor_module_final_modules()` to `cor_module_graph_final_modules()`.

### Version **0.0.1.3**

**Beta release v3** 

#### New features

- Implementation for fgsea multi level method. Additionally, n_more_extreme
is now being returned by the GSEA functions.
- Implementations of the RBF function for R matrices.

#### Bug fixes, documentation updates

...

### Version **0.0.1.2**

**Beta release v2** 

#### New features

- Implementation for fgsea simple on the gene ontology with elimination method.
- Improvement in the gene ontology elimination methods to reduce unneeded 
copying.

#### Bug fixes, documentation updates

- Minor documentation updates and fixes.

### Version **0.0.1.1**

**Beta release** Improved stability and bug fixes.

#### New features

- Semantic similarities for ontologies added.
- Speed improvements in various (Rust) functions due to less unnecessary copying
of data in memory in the Rust back-end and leverage of lifetimes. 
- Stability of the package massively improved with various test, powered by
[tinytest](https://github.com/markvanderloo/tinytest).
- Set similarities rust functions exposed.
- DGE class for leveraging Limma Voom to do differential gene expression 
calculations.
- Wrapper functions into h5ad objects to load data into R memory.
- Traditional gene set enrichment analysis and the simple version from the fgsea
package ported into Rust.

#### Bug fixes, documentation updates

- Various documentation updates, fixes in spelling etc.
- Bug fix for the `future::plan()` for iterating over different resolutions in
the reciprocal best hit graph generation.
- Bug fix in the hypergeometric calculations and RBH graph.

#### Breaking changes

- The `community_detection()` function now uses `params_community_detection()`
to retrieve parameters. Prior versions would not work anymore.

### Version **0.0.1.0**

**Alpha release** Package initialisation and first functionalities added.

#### New features

- Hypergeometric tests for gene set analysis for single target gene sets or list
of target gene sets.
- Implementation of a gene ontology aware hypergeometric test for gene ontology
analysis (using hypergeometric tests under the hood).
- Network diffusion methods based on personalised page-rank (single diffusion
and tied diffusion).
- Community detection algorithms on top of diffusion scores over network.
- Reciprocal best hit graphs using set similarities between different gene modules
from different data sets/methods.
- Contrastive PCA implementation for gene module detection.
- (Differential) Correlation-based methods for gene module detection, using 
network-based detection of graph communities.
- Independent component analysis (with the versions to run stabilised versions).
- Smaller helper function to calculate the Hedge's G effect size or using the
OpenTargets method to summarise scores over different evidences.

#### Bug fixes, documentation updates

*NA*

#### Breaking changes

*NA*
