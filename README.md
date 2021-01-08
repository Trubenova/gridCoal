# gridCoal Tutorial

## Introduction

gridCoal is an efficient spatially explicit coalescent simulation tool (Szep, Trubenova, Csillery 2021) and is a wrapper for the software msprime (Kellher et al. ). It implements a flexible two-dimensional stepping stone model, where the size of demes can change arbitrarily across space and in time, and the migration rate between individual demes can be specified. Under Kimura’s two dimensional stepping stone model (), summary statistics of genetic diversity and divergence between demes can be approximated from coalescence times (Slatkin 1985). Bypassing the simulation of genetic data in a spatially explicit context, combined with the efficiency of msprime. Each gridCoal simulation returns coalescent times between a set of sample pairs, that can be used to calculate summary statistics such as FST or generate isolation by distance patterns. Many simulation repeats (1000-3000) of the same scenario (parameter set) are required to obtain reliable estimates of the statistics.

Note that approximating summary statistics of genetic diversity and FST from coalescence times holds only when the mutation rate is low, and migration is possible to neighboring demes only (Slatkin 1985). The assumption about the migration to neighboring cells may be a limitation for simulating data for some organisms, such as migrating animals or wind-dispersed plants. Thus, even though long-distance dispersal events could be easily simulated by directly modifying the migration matrix, it is advisable to simulate genetic data in these cases for calculating summary statistics.


## How to use:

Main.py is the main python file that is necessary for running the simulations. Python package msprime needs to be installed.

## Table of contents

GridCoal Tutorial.ipynb is a jupyter notebook file that explains input and output files and how to use them.

AnalysingOutput.ipynb is a jupyter notebook file that allows to analyse outputs of simulation repeats (coalescent times) and calculte Fst as well as to create isolation by distance patterns.

InputFileGenerator.ipynb is a jupyter notebook file that allows to prepare simple inputs for the simulator. Note that this functionality is limited, and is mainly useful for testing purposes.

Examples is a directory, in which complete examples are worked out, together with input files.

## References

* Kelleher J, Etheridge AM, McVean G (2016) Efficient coalescent simulation and genealogical analysis for large sample sizes.PLoS Computational Biology,12, e1004842

* Slatkin M (1985) Gene #ow in natural populations. Annual Review of Ecology and Systematics,16: 393–430

* Szep, Trubenova, Csillery 2021
