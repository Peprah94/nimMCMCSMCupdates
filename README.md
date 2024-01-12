# nimMCMCSMCupdates

# Introduction

This package comes with a sequential Monte Carlo algorithm for data
assimilation problems in ecology. This package is editted from the
*nimbleSMC* package Michaud et al. (2021), so the users of this package
need to be familiar with the NIMBLE software through its R-package
*nimble* (de Valpine et al. 2022) and then its R-package for sequential
Monte Carlo problems *nimbleSMC* (the reader is referred to Chapter 8 of
de Valpine et al. (2022) and Michaud et al. (2021) for details on how to
fit SSMs using SMC approach in NIMBLE).

## Installation

You can install the development version of *nimMCMCSMCupdates* from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("Peprah94/nimMCMCSMCupdates")
```

## Functions in the R-package

*nimMCMCSMCupdates* includes a set of functions that are used to update
posterior distribution of latent states and model parameters using the
bootstap and auxiliary particle filters. Refer to the \[main
paper\]\[https://peprah94.github.io/#Publications\] and vignette on
\[method\]\[https://github.com/Peprah94/nimMCMCSMCupdates/blob/main/vignette/methodUnderstanding.pdf\]
for further details on the methodology.

## Updating models with the R-package

## Other Materials

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-nimblepackage" class="csl-entry">

de Valpine, Perry, Christopher Paciorek, Daniel Turek, Nick Michaud,
Cliff Anderson-Bergman, Fritz Obermeyer, Claudia Wehrhahn Cortes, Abel
Rodrìguez, Duncan Temple Lang, and Sally Paganin. 2022. *NIMBLE User
Manual* (version 0.13.1). <https://doi.org/10.5281/zenodo.1211190>.

</div>

<div id="ref-michaud2021sequential" class="csl-entry">

Michaud, Nicholas, Perry de Valpine, Daniel Turek, Christopher J
Paciorek, and Dao Nguyen. 2021. “Sequential Monte Carlo Methods in the
Nimble and nimbleSMC r Packages.” *Journal of Statistical Software* 100:
1–39.

</div>

</div>
