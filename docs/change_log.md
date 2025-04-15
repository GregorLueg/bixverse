# *bixverse* change log

*Last update: 13.04.2025* </br>

### Version **0.0.1.1**

**Beta release** Package initialisation and first functionalities added.

#### New features

- Semantic similarities for ontologies added.
- Speed improvements in various functions due to less unnecessary copying of 
data in memory in the Rust back-end.
- First set of tests for package stability and quick identification of 
unexpected behaviour via [tinytest](https://github.com/markvanderloo/tinytest).

#### Bug fixes, documentation updates

- Various documentation updates, fixes in spelling etc.
- Bug fix for the `future::plan()` for iterating over different resolutions in
the reciprocal best hit graph generation.

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
