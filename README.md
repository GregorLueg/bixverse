# *bixverse package*

![r_package](https://img.shields.io/badge/R_package-0.0.2.0-orange) 

</br>

<img src="/misc/pics/bixverse_logo.png" width="128" height="128" alt="bixverse logo">

</br>

## *Description* 

This package contains various bioinformatics and computational biology workflows
that are being routinely used, ranging from gene set enrichment analyses, to 
network-based approaches for module detection in bulk RNAseq. The package 
provides useful, bare bone versions of most bioinformatics functionalities and
leverages Rust to make any computational bottlenecks go *brrrrrrr* (i.e., fast).

## *Release notes*

We have now officially released version **0.0.2.0**. The package has now reached
high degrees of maturity with more and more tests in place, vignettes to explain
various workflows and further bug and documentation fixes. With this update the
following has been updated/changed:

- Improved ontology class and more methods, and the addition of the Wang
similarity measure.
- Permutation-based tests for the genetic diffusion approaches that allow selection
of significanlty enriched areas of the network (topology-aware).
- Updates to the ICA detection approaches and correlation-based reciprocal best
hits.
- Further improvements in speed in various Rust functions (less unnecessary
copying and changes to the used HashMaps and HashSets).
- Addition of the simplify() type method to reduce Gene Ontology results 
to the most relevant ones.
- Vignettes explaining different methods.

**Warnings**: Some of the previous interfaces to the functions were changed
in this release and might break compared to prior releases!

## *Installation*

You will need Rust on your system to have the package working. An installation
guide is provided [here](https://www.rust-lang.org/tools/install). There is a 
bunch of further help written [here](https://extendr.github.io/rextendr/index.html)
by the rextendr guys in terms of Rust set up. (bixverse uses rextendr to interface
with Rust.)

Steps for installation: 

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
1. Initially, the package was born out of the desire to harmonise in a single 
package lots and lots of different packages. (One package to rule them all.) 
And while doing so, slim the functionality down to the bare bone essentials and
remove barely used options.
2. After having reimplemented for the 3rd time some version of a 
hypergeometric test for some project, one is just tired of it. The same applies 
to other established methods in the field that need some internal reimplementation
(with maybe some tweaks here and there) for reasons... ? Better to do it once
properly and make it public.
3. Rust makes everything in R so much faster, that we just wished to share the 
joy.

*Last update to the read-me: 13.07.2025*
