# gridCoal Tutorial

## Table of contents

<!-- TOC START min:2 max:3 link:true asterisk:false update:true -->
- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
  - [References](#references)
- [How to use](#how-to-use)
- [Input files and parameters](#input-files-and-parameters)
  - [Demographic history input files [txt file, required]](#demographic-history-input-files-txt-file-required)
  - [Row number [int, required]](#row-number-int-required)
  - [Migration list [txt, optional]](#migration-list-txt-optional)
  - [List of sampled demes   [txt, optional]](#list-of-sampled-demes---txt-optional)
  - [Time periods and generation times [int, optional]](#time-periods-and-generation-times-int-optional)
  - [Replicate ID numer [int, optional]](#replicate-id-numer-int-optional)
  - [Ancestral populations [txt, optional]](#ancestral-populations-txt-optional)
  - [Output directory [string, optional]](#output-directory-string-optional)
  - [Printing demography file [bool, optional]](#printing-demography-file-bool-optional)
  - [Setting random seed number [int, optional]](#setting-random-seed-number-int-optional)
- [Running the simulation](#running-the-simulation)
- [Output](#output)
  - [Analysing output files](#analysing-output-files)
<!-- TOC END -->



## Introduction

gridCoal is an efficient spatially explicit coalescent simulation tool (Szep, Trubenova, Csillery 2021) and is a wrapper for the software msprime (Kellher et al. ). It implements a flexible two-dimensional stepping stone model, where the size of demes can change arbitrarily across space and in time, and the migration rate between individual demes can be specified. Under Kimura’s two dimensional stepping stone model (), summary statistics of genetic diversity and divergence between demes can be approximated from coalescence times (Slatkin 1985). Bypassing the simulation of genetic data in a spatially explicit context, combined with the efficiency of msprime. Each gridCoal simulation returns coalescent times between a set of sample pairs, that can be used to calculate summary statistics such as FST or generate isolation by distance patterns. Many simulation repeats (1000-3000) of the same scenario (parameter set) are required to obtain reliable estimates of the statistics.

Note that approximating summary statistics of genetic diversity and FST from coalescence times holds only when the mutation rate is low, and migration is possible to neighboring demes only (Slatkin 1985). The assumption about the migration to neighboring cells may be a limitation for simulating data for some organisms, such as migrating animals or wind-dispersed plants. Thus, even though long-distance dispersal events could be easily simulated by directly modifying the migration matrix, it is advisable to simulate genetic data in these cases for calculating summary statistics.

### References

* Kelleher J, Etheridge AM, McVean G (2016) Efficient coalescent simulation and genealogical analysis for large sample sizes.PLoS Computational Biology,12, e1004842

* Slatkin M (1985) Gene #ow in natural populations. Annual Review of Ecology and Systematics,16: 393–430

* Szep, Trubenova, Csillery 2021

## How to use

Main.py is the main python file that is necessary for running the simulations. Python package msprime needs to be installed.

Below is a detailed description of required and optional inputs for gridCoal. A jupyter notebook file, InputFileGenerator.ipynb is present for your convenience, that allows  to prepare simple inputs for the simulator. Note that its functionality is limited, and is mainly useful for testing and tutorial purposes.

Description of output files follows. Another jupyter notebook file, AnalysingOutput.ipynb allows to analyse outputs of simulation repeats (coalescent times) and calculte Fst as well as to create isolation by distance patterns.

Example is a directory, in which a complete example are worked out, together with input files.




## Input files and parameters

GridCoal requires following inputs:

### Demographic history input files [txt file, required]


The demographic history of the collection of demes distributed on a grid is represented by matrix of size $T \times n$ (T rows of n numbers),  $T$ being the number of time points one wishes to define the population sizes at, and $n$ being the number of grid cells. The matrix contains the population sizes of the grid cells at given time points. Each row is the flattened two-dimensional grid, indexed from $0$ to $n-1$, defining the sizes of the demes. The first line corresponds to the oldest time point.

Empty demes, i.e. zero population sizes, are not allowed in msprime, thereby in gridCoal. To allow for the possibility that populations can go extinct and get recolonized, even repeatedly, upon loading the demography input file we replace all empty demes with population size  $10^{-10}$.

#### Example
InputData.txt:

>1	2	3	4	2	3	4	5	3	4	5	5	4	5	6	7 <br>
>1	1	1	1	2	2	2	2	3	3	3	3	4	4	4	4 <br>
>1	2	3	4	2	3	4	5	3	4	5	5	4	5	6	7 <br>

Enter as: <br>
<code>--pop_sizes InputData.txt</code>
or
<code>-pop InputData.txt</code>



### Row number [int, required]
The number of rows in the grid. The number of columns is calculated from the number of columns in the Demographic history input file. The total number of demes must be divisible by the row number.
The row number together with the demographic history file defines the shape of the spatial map.

Enter as: <br>
 <code>--row_number 5 </code> or <code>-row 5 </code>

### Migration list [txt, optional]

A migration file specifies the migration rate between pairs of connected cells in a list, and to calculate backward migration during a simulation. Even if two cells are connected by an edge, if the migration is not specified in this file, it is considered 0. Format as lines of triples: Source cell   (integer), target cell   (integer), migration rate  (positive float between 0 and 1). Note that if migration exchange migrants both ways, both directions need to be specified, by two lines of triples.

If no file, but a single number m (float) is specified, a migration matrix is automatically generated, assuming that migration occurs between all adjacent cells symmetrically with rate m (i.e. the classical stepping stone model). By default migration rate is 0.1.![image.png](attachment:image.png)

#### Example
MigrationList.txt:

>0 1 1.00e-06<br>
>0 5 1.00e-06<br>
>1 0 1.00e-06<br>
>1 2 1.00e-06<br>
>1 6 1.00e-06<br>
>2 1 1.00e-06<br>
>2 3 1.00e-06<br>
>2 7 1.00e-06<br>
>3 2 1.00e-06<br>

Enter as:

<code>--migration_matrix MigrationList.txt</code>
or
<code>--migration_matrix 0.000001</code>
or
<code>-mig MigrationList.txt</code>
or
<code>-mig 0.000001</code> <br>


### List of sampled demes   [txt, optional]
A  list of demes (indexes, starting at 0) from which the samples are taken can be specified.
These cells must not be empty at the final time point (presence), but could be empty in the past.
Two samples are taken from each indicated deme.

If no file is supplied, all samples that are not empty at presence are sampled, and a Sample_list.txt file is created in the output directory. 

#### Example
SampleList.txt:

>0 1 2 3 4 7 8 33 34 35<br>


Enter as: <br>
<code>--sample_coords SampleList.txt</code> or <code>-sam SampleList.txt</code>


### Time periods and generation times [int, optional]

Time between two time steps, denoted <i>dt</i>, is given in arbitrary time units (year, months, days, minutes). Time is measured in generations in <i>msprime</i>, therefore we need to specify the generation time of the population at hand. We define the generation time as the time at which the species comes to a reproductive age in same units, as dt.  Timing of demographic events (into the same units) is re-calculated by dividing the time of events specified in demography input files by the generation time of the simulated organism in same units.   

Be default, generation time is set to one, and time between two supplied data points is defined as ten.
To supply different values, enter as: <br>
<code>-dt 100 -gen 25</code>
or
<code>--delta_t 100 --generation_time 25</code>


### Replicate ID numer [int, optional]
Replicate ID number can be suplied specifying the simulation run and the output file. Default value is 1. If value is 1, log file with all input files and their values is created. When running multiple simulations in parallel, different replicate ID numbers need to be specified to avoid overwriting the output files.

Enter as: <br>
<code>--replicate 7</code> or <code>-rep 7</code>

### Ancestral populations [txt, optional]

At the time point, beyond which the demography is unknown or ignored, all lineages are merged into one or more spatially non-explicit panmictic ancestral populations that follow the standard coalescent process. Multiple ancestral populations allow users to define ancestral populations with specified sizes and migration between them ($10^{-8}$). Furthermore, users can specify which of cells originate from which ancestral population. A list determining the origin of each cell can be supplied as a txt file with n (number of all cells) lines. If no file is supplied, all cells are expected to originate from a single ancestral population. The default size of the ancestral populations is one.


### Output directory [string, optional]

The output directory into which outputs (log file with input parameters, demography debugger, random seed numbers and coalescent times) are saved. Default value is OUTPUT.!
Enter as: <br>
<code>--output_dir MY_OUTPUT_DIR</code>
or
<code>-odir MY_OUTPUT_DIR</code>

### Printing demography file [bool, optional]

Option that allows users printing out the detailed demography debugger file supplied by msprime. Replicate number must be set to one to use this option. This allows to start simultaneous runs of many simulations, with only one debugger file (identical for all simulations) created.

Enter as: <br>
<code>--print_debugger BOOL</code>
or
<code>-pdeb BOOL</code>

### Setting random seed number [int, optional]

This option allows reproducing the same simulation.
Enter as: <br>
<code>--set_seed 19</code>
or
<code>-seed 19</code>

## Running the simulation


Simulation is run by specifying all the necessary input files, in terminal (on Mac) for instance as

<code> python3 Main.py -pop DEMOGRAPHY_INPUT [txt file name with txt extension] -row ROW_NUMBER [integer]
</code>

other optional parameters:

<code>-mig MIGRATION_FILE[txt file name without txt extension]
-sam SAMPLE_LIST [txt file name without txt extension ,otherwise all non-empty cells]
-apl ANCESTRAL_POPULATION_LIST [txt file name without txt extension ,otherwise 1 for all]
-pdeb PRINTIN_DEMOGRAPHY [True or False, default False]
-odir OUTPUT_DIR_NAME [string ]
-d TIME_PERIOD [integer, default 10]
-gen GENERATION_TIME [integer, default 1]
-pdeb True [boolean]
-rep REPLICATE [integer]
</code>

For multiple runs:

<code>for i in {1..100} ; do
  python3 Main.py \
    -pop Test1_dataT10Row4Col5N10 \
    -sam Test1_sample_list \
    -mig Test1_mig_list \
    -row 4 \
    -odir XXX \
    -pdeb True \
    -rep $i &
done
</code>

## Output

For <code>REPLICATE = 1</code>, All the input files are collected and saved into a created <code>OUTPUT_DIR_NAME</code> directory as <code>OUTPUT_DIR_NAME/Output.txt</code>. Random seed number is also saved in the same file.

If --print_debugger is set to True, also detailed demographic history with all population and migration rate changes is produced by msprime and printer into OUTPUT_DIR_NAME directory as <code>OUTPUT_DIR_NAME/DemographyDebugger.txt</code>.


The result of the simulation itself is a square matrix of coalescence times of samples from all sampled demes, saved as <code>OUTPUT_DIR_NAME/CoalTimes${REPLICATE}.txt</code>

### Analysing output files

Jupyter notebook file AnalysingOutput allows to calculate summary statistics from the simulation outputs.
The original input files used for the simulations, as well as output files with coalescence times are specified in the first cell. The following cells define the function and run them, to calculate and print out:
* mean coalescence time for samples taken from the same deme ('within deme coalescence time'), as well as from different demes ('between deme coalescence time'), and the total mean coalescence time ('mean total coalescence time').  
* population wide Fst
* F* for samples taken from various (Manhattan) distance classes.
