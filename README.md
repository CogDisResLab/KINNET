
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KINNET

<!-- badges: start -->
<!-- [![DOI](https://zenodo.org/badge/301042598.svg)](https://zenodo.org/badge/latestdoi/301042598) -->

[![GitHub
issues](https://img.shields.io/github/issues/AliSajid/KINNET)](https://github.com/AliSajid/KINNET/issues)
[![GitHub
forks](https://img.shields.io/github/forks/AliSajid/KINNET)](https://github.com/AliSajid/KINNET/network)
[![GitHub
stars](https://img.shields.io/github/stars/AliSajid/KINNET)](https://github.com/AliSajid/KINNET/stargazers)
[![Travis-CI](https://img.shields.io/travis/com/AliSajid/KINNET)]()
[![codecov](https://codecov.io/gh/AliSajid/KINNET/branch/main/graph/badge.svg?token=HP1VFCRZ0M)]()

![GitHub release (latest SemVer including
pre-releases)](https://img.shields.io/github/v/release/AliSajid/KINNET?include_prereleases&label=latest-release)
![GitHub release (latest
SemVer)](https://img.shields.io/github/v/release/AliSajid/KINNET?label=latest-stable)
[![GitHub
license](https://img.shields.io/github/license/AliSajid/KINNET)](https://github.com/AliSajid/KINNET/blob/main/LICENSE)

![GitHub language
count](https://img.shields.io/github/languages/count/AliSajid/KINNET)
![GitHub top
language](https://img.shields.io/github/languages/top/AliSajid/KINNET)
![Lines of
code](https://img.shields.io/tokei/lines/github/AliSajid/KINNET)
![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/AliSajid/KINNET)
![GitHub repo
size](https://img.shields.io/github/repo-size/AliSajid/KINNET) [![Travis
build
status](https://travis-ci.com/AliSajid/KINNET.svg?branch=master)](https://travis-ci.com/AliSajid/KINNET)
<!-- badges: end -->

This package provides the functionality to process PamGene’s PamChip
Data Output and generate kinase interaction networks from that.

This project uses a bayesian algorithm to generate bayesian networks for
defining dependence relationships between peptide sequences in the
PamChip data. It then uses a novel kinase assignment method to assign
upstream kinases to each peptide which is then output as a graph.

## Installation

You can install the released version of KINNET from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("KINNET")
```

You can install the latest development version of KINNET from [Github]()
with:

``` r
remotes::install_github("AliSajid/KINNET")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(KINNET)
## basic example code
```
