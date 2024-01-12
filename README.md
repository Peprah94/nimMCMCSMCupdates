# nimMCMCSMCupdates

# Introduction

This package comes with a sequential Monte Carlo algorithm for data
assiilation problems in ecology. This document seeks to demonstrate how
once can use MCMC models already fitted to the SSMs and update them
using the SMC approach. The bootstrap and auxiliary PFs were discussed
in the main paper, but this document will focus on bootstrap particle
filter. The reader is expected to be familiar with the `nimbleSMC`
package and its functionalities (the reader is referred to Chapter 8 of
de Valpine et al. (2022) and Michaud et al. (2021) for details on how to
fit SSMs using SMC approach in NIMBLE). The first part provides a brief
introduction to the state space models (SSMs), sequential Monte Carlo
(SMC) approaches (specifically the bootstrap particle filter) and the
particle MCMC with the necessary changes made to accomodate the updating
process when new streams of data are obtained. The next part shows how
to fit the model to the simulated data used in Example 1 in the main
paper.

## Functions in the R-package

## Updating models with the R-package

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
