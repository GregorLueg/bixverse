# *bixverse package*

![r_package](https://img.shields.io/badge/R_package-0.0.2.0-orange) 

</br>

<img src="/misc/pics/bixverse_logo.png" width="128" height="128" alt="bixverse logo">

</br>

## *Description* 

This package contains various bioinformatics and computational biology workflows
that are being routinely used, ranging from gene set enrichment analyses in all
colours and forms, various analysis methods for networks to different methods
for co-expression module detection in bulk RNAseq. The package provides simple, 
bare bone versions of quite a few functions from your favourite 
[BioConductor](https://www.bioconductor.org/) packages and leverages Rust to 
make any computational bottlenecks in these go **brrrrrrr** (i.e., fast).

## *Release notes*

This is now officially the version **0.0.2.0** release. The package has now 
reached high degrees of maturity with more and more tests in place, vignettes 
to explain various workflows and further bug and documentation fixes. With this 
update the following has been updated/changed:

- Improved ontology class and more methods, and the addition of the Wang
similarity measure.
- Permutation-based tests for the genetic diffusion approaches that allow 
selection of significanlty enriched areas of the network (topology-aware).
- Updates to the ICA detection approaches and correlation-based reciprocal best
hits.
- Further improvements in speed in various Rust functions (less unnecessary
copying and changes to the used HashMaps and HashSets to faster versions).
- Addition of a `simplify()`-type method to reduce Gene Ontology results 
to the most relevant ones.
- Vignettes explaining different methods.
- GSVA and ssGSEA implemented leveraging Rust with in parts substantial speed
increases.
- Update of the `furrr::future_map()` type functions to the faster 
`mirai::mirai_map()` for concurrent/parallel operations in R.
- Mutual information calculations between continuous variables and also 
normalised pointwise mutual information calculations between boolean type
variables.
- Improvements to the synthetic data with better plotting functionality.

**Warnings**: Some of the previous interfaces to the functions were changed
in this release and might break compared to prior releases!

## *Installation*

You will need Rust on your system to have the package working. An installation
guide is provided [here](https://www.rust-lang.org/tools/install). There is a 
bunch of further help written [here](https://extendr.github.io/rextendr/index.html)
by the rextendr guys in terms of Rust set up. (bixverse uses rextendr to 
interface with Rust.)

Steps for installation (assuming MacOS or Linux-based systems): 

1. In the terminal, install [Rust](https://www.rust-lang.org/tools/install) 

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
   
2. In R, install [rextendr](https://extendr.github.io/rextendr/index.html):

```
install.packages("rextendr")
```

3. Finally install bixverse:

```
devtools::install_github("https://github.com/GregorLueg/bixverse")
```

## *Docs*

- [Roadmap](/docs/roadmap.md)
- [Change log](/docs/change_log.md)
- [Why Rust](/docs/why_rust.md)
- [Code Style](/docs/code_style.md) (If you want to contribute).

## *Aim*

<img src="/misc/pics/but_why.png" width="418" height="218" alt="but why">

Good question, why this package? Basically, there are three reasons for this:
1. Rust makes a lot of things in R so much faster and who does not like faster
code? The aim is to share the joy of speedy code. Additionally, it's a fun 
language to write.
2. BioConductor is great, but you end up with quite large dependency graphs, 
lots of packages that have to be downloaded and quite a few algorithms can 
benefit from the speed and memory control that Rust offers. The idea of 
`bixverse` is to accelerate key functions used in CompBio and bioinformatics and
reduce them to the barebone methods with simple interfaces (which can be
easily used in your workflows).
3. Having worked in biotech and pharma, one realises that a surprising amount
of time is spent on rewriting and reimplementing published methods for internal
usage. Better to make it fast/good once, and open it up to the public via open 
source.

*Last update to the read-me: 06.08.2025*
