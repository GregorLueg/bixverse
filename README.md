# bixverse package

![r_package](https://img.shields.io/badge/R_package-0.0.2.2-orange)

</br>

<img src="man/figures/bixverse_logo.png" width="128" height="128" alt="bixverse logo">

</br>

## Description

This package contains various bioinformatics and computational biology workflows
that are being routinely used, ranging from gene set enrichment analyses, to
network-based approaches for module detection in bulk RNAseq. The package
provides useful, bare bone versions of most bioinformatics functionalities and
leverages Rust to make any computational bottlenecks go _brrrrrrr_ (i.e., fast).

## Release notes

This is now officially the version **0.0.2.2** release. With this
update the following has been updated/changed:

- Mitch multi-contrast enrichment, see
  [Kaspi and Ziemann](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06856-9).

## Installation

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

## Docs

- [Roadmap](/info/roadmap.md)
- [Change log](/info/change_log.md)
- [Why Rust](/info/why_rust.md)
- [Code Style](/info/code_style.md) (If you want to contribute).

## Aim

<img src="man/figures/but_why.png" width="418" height="218" alt="but why">

Good question, why this package? Basically, there are three reasons for this:

1. Rust makes everything in R so much faster, that we just wished to share the
   joy. Also, it's a fun language write.
2. BioConductor is great, but you end up with quite large dependency graphs,
   lots of packages that have to be downloaded and quite a few algorithms can
   benefit from the speed that Rust offers. The idea of `bixverse` is to accelerate
   key functions used in CompBio and bioinformatics and reduce them to the
   barebone methods with simple interfaces.
3. Having worked in biotech and pharma, one realises that a surprising amount
   of time is spent on rewriting and reimplementing published methods for internal
   usage. Better to make it fast/good once, and open it up to the public via open
   source.

_Last update to the read-me: 07.09.2025_
