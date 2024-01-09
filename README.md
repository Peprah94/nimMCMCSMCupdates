# nimMCMCSMCupdates: A package with a sequential Monte Carlo algorithm for data assiilation problems in ecology.

## Introduction

This document seeks to demonstrate how once can use MCMC models already fitted to the SSMs and update them using the SMC approach. The bootstrap and auxiliary PFs were discussed in the main paper, but this document will focus on bootstrap particle filter. The reader is expected to be familiar with the `nimbleSMC` package and its functionalities (the reader is referred to Chapter 8 of @nimblepackage and @michaud2021sequential for details on how to fit SSMs using SMC approach in NIMBLE). The first part provides a brief introduction to the state space models (SSMs), sequential Monte Carlo (SMC) approaches (specifically the bootstrap particle filter) and the particle MCMC with the necessary changes made to accomodate the updating process when new streams of data are obtained. The next part shows how to fit the model to the simulated data used in Example 1 in the main paper.

## Functions in the R-package


## Updating models with the R-package

