# Readme file for gridCoal simulation tool

## Introduction

gridCoal  is an efficient spatially explicit coalescent simulation tool, a wrapper for the software msprime. It implements a flexible two-dimensional stepping stone model, where the size of demes can change arbitrarily across space and in time, and the migration rate between individual demes can be specified. Under this model, summary statistics of genetic diversity and divergence between demes can be approximated from coalescence times. Bypassing the simulation of genetic data in a spatially explicit context, combined with the efficiency of msprime.

Each gridCoal simulation returns coalescent times between a set of sample pairs, that can be used to calculate summary statistics such as Fst or generate isolation by distance patterns.


## How to use:

Main.py is the main python file that is necessary for running the simulations. Python package msprime needs to be installed. Many simulation repeats (1000-3000) of the same scenario (parameter set) are required in order to achieve reasonable statistics.

GridCoal Tutorial.ipynb is a jupyter notebook file that explains input and output files and how to use them.

AnalysingOutput.ipynb is a jupyter notebook file that allows to analyse outputs of simulation repeats (coalescent times) and calculte Fst as well as to create isolation by distance patterns.

InputFileGenerator.ipynb is a jupyter notebook file that allows to prepare simple inputs for the simulator. Note that this functionality is limited, and is mainly useful for testing purposes.

Examples is a directory, in which complete examples are worked out, together with input files.
