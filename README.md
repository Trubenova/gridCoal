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

gridCoal is an efficient, spatially explicit coalescent simulation tool -- a wrapper for the software msprime (Kellher et al. 2016). It implements a flexible two-dimensional stepping stone model, where the size of demes can change arbitrarily across space and in time, and the migration rate between individual demes can be specified. Under Kimura’s two dimensional stepping stone model, summary statistics of genetic diversity and divergence between demes can be approximated from coalescence times (Slatkin 1985). Bypassing the simulation of genetic data in a spatially explicit context, combined with the efficiency of msprime. Each gridCoal simulation returns coalescent times between a set of sample pairs, that can be used to calculate summary statistics such as Fst or generate isolation by distance patterns. Many simulation repeats (1000-3000) of the same scenario (parameter set) are required to obtain reliable estimates of the statistics.

Note that approximating summary statistics of genetic diversity and FST from coalescence times holds only when the mutation rate is low, and migration is possible to neighboring demes only (Slatkin 1985). The assumption about the migration to neighboring cells may be a limitation for simulating data for some organisms, such as migrating animals or wind-dispersed plants. Thus, even though long-distance dispersal events could be easily simulated by directly modifying the migration matrix, it is advisable to simulate genetic data in these cases for calculating summary statistics.

### References

* Kelleher J, Etheridge AM, McVean G (2016) Efficient coalescent simulation and genealogical analysis for large sample sizes.PLoS Computational Biology,12, e1004842

* Slatkin M (1985) Gene #ow in natural populations. Annual Review of Ecology and Systematics,16: 393–430

* Szep, Trubenova, Csillery. Under review in Molecular Ecology Resources. 

## How to use

Main.py is the main python file that is necessary for running the simulations. Python package msprime needs to be installed.

Below is a detailed description of required and optional inputs for gridCoal. A jupyter notebook file, InputFileGenerator.ipynb is present for your convenience, that allows you to prepare simple inputs for the simulator. Note that its functionality is limited, and is mainly useful for testing and tutorial purposes.

Description of output files follows. Another jupyter notebook file, AnalysingOutput.ipynb allows to analyse outputs of simulation repeats (coalescent times) and calculte Fst as well as to create isolation by distance patterns.

There are two completely worked out examples provided. The first one - Simulating increase in population sizes - generates gridCoal inputs from defined parameters. Then gridCoal simulations are performed and the results analysed. The second example uses input data obtained from the LPX-Bern
dynamic global vegetation model. 

## Input files and parameters

In order to run the simulations, it is necessary to define the following input files and parameters.

### Demographic history input files [txt file, required]

The demographic history of the collection of demes distributed on a grid is represented by a matrix of size $T \times n$, where $T$ is the number of time points at which one wishes to define the population sizes and $n$ is the number of grid cells. The matrix contains the population sizes of the grid cells at given time points. Each row is the flattened two dimensional grid, indexed from $0$ to $n-1$, defining the sizes of the subpopulations. The first line is the oldest time point. 

In \textit{msprime}, a population is not allowed to have size $0$. In our case, however, we do not want to exclude the possibility that populations become extinct and the demes are subsequently recolonised, even repeatedly. We therefore set populations with size 0 as $10^{-10}$. This is done automatically -- before the simulations start, the program replaces any $0$ in the input data with $10^{-10}$. 

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
This number, together with the demographic history file, defines the shape of the spatial map. It must be an integer. Additionally, the size of the grid (number of all cells -- row length of the demographic history input file) must be divisible by the row number. 

Enter as: <br>
 <code>--row_number 5 </code> or <code>-row 5 </code>

### Migration list [txt, optional]
A migration list indicating the migration rate between two connected cells is used to build a migration matrix, and to calculate backward migration during the simulation. Even if two cells are connected by an edge, if the migration is not specified in the list, it is considered 0. The migration list must be formatted as lines of three values: source cell i (integer), target cell j (integer), and migration rate m_{ij} from i to j (positive float between 0 and 1). Note that if migration involves the exchange of migrants, both directions need to be specified.

If a file is not specified but rather a single number m (float) is supplied, a migration matrix is generated in which migration is assumed to occur between adjacent cells symmetrically with rate m.
If no value is suplied, the default migration rate is 0.1.

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
hese cells must not be empty at the final time point (present), but could be empty in the past. 
For efficiency, two samples are taken from each sampled deme. 

If only a number (float, <1) is supplied, random, non-empty cells will be sampled, with the number representing the sampled fraction.  
If no file is supplied, all demes that are not empty at presence are sampled, and a Sample_list.txt file is created in the output directory. 



#### Example
SampleList.txt:

>0 1 2 3 4 7 8 33 34 35<br>


Enter as: <br>
<code>--sample_coords SampleList.txt</code> or <code>-sam SampleList.txt</code>


### Time periods and generation times [int, optional]

The amount of time between two time points, denoted <i>dt</i>, is given in arbitrary time units (years, months, days, minutes). 

Time is measured in generations in <i>msprime</i>, and we therefore need to specify the generation time of the population at hand. We define the generation time, dt, as the time it takes for a species to reach a reproductive age, expressed in the same units as other supplied times. The timing of demographic events (expressed in the same units) is re-calculated by dividing the time point of events specified in demography input files by the generation time of the simulated organism expressed in the same units. Therefore, it is possible to run the simulations for any organism with an arbitrary generation time, from bacterial populations to trees. 

By default, generation time is set to 1, and the time between two supplied data points is set at 10.

To supply different values, enter as: <br>
<code>-dt 100 -gen 25</code>
or
<code>--delta_t 100 --generation_time 25</code>

### Ancestral populations [txt, optional]

At the point in time beyond which the demography is unknown, all lineages are merged into spatially non-explicit ancestral populations where they follow the standard coalescence process. We assume either a single or multiple panmictic ancestral populations with specified sizes and a very low rate of migration between them ($10^{-8}$). Furthermore, it is necessary to specify which of the cells originate in each ancestral population.

A list determining the origin of each cell can be supplied as a txt file with n (number of all cells) lines. If no file is supplied, all cells are expected to originate in a single spatially non-explicit population. 
The size of the ancestral populations can be set, with the default as 1. 

### Replicate ID numer [int, optional]
Replicate ID number can be suplied specifying the simulation run and the output file. Default value is 1. If value is 1, log file with all input files and their values is created. When running multiple simulations in parallel, different replicate ID numbers need to be specified to avoid overwriting the output files.

Enter as: <br>
<code>--replicate 7</code> or <code>-rep 7</code>


### Output directory [string, optional]

The output directory into which outputs (log file with input parameters, demography debugger, random seed numbers and coalescent times) are saved. Default value is OUTPUT.!
Enter as: <br>
<code>--output_dir MY_OUTPUT_DIR</code>
or
<code>-odir MY_OUTPUT_DIR</code>

### Printing demography file [bool, optional]


This option can be used to print a detailed demography debugger file, supplied by <i>msprime</i>. The replicate number must be set to 1. This makes it possible to simultaneously run many simulations with only one debugger file (identical for all).

Enter as: <br>
<code>--print_debugger BOOL</code>
or
<code>-pdeb BOOL</code>

### Setting random seed number [int, optional]

Finally, it is possible to set a <i>random seed number</i>, which makes it possible to reproduce a given simulation. 

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
Furthermore, Spearman Rank correlation between mean coalescence times and population demography in the past can be calculated, indicating which time points (e.d. near or distant past) are correlated with the current diversity the most. Note that in order to do this analysis, the whole map (all demes that are not empty in the final point) must be sampled. 
