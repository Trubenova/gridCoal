# Coalescence simulations on LPX Bern model data

## Input data


The demographic history of silver fir in the past 22,000 years was obtained from the LPX-Bern dynamic global vegetation model with a resolution of 1 x 1 degree.  The output of LPX-Bern is the foliar projective cover (FPC), which is the fraction of a grid cell that is covered by silver fir. 

We estimated the number of trees (N) from FPC assuming that a mature tree occupies 40 m x m and effective population size is Ne= c x N where c is a scaling factor (0.0001, 0.001,0.01,0.1). 
Furthermore, we created 3 types of data: a) static and homogeneous in space (in the predicted fir range), b) static in time but based on the final distribution of silver fir as predicted by the LPX Bern model,  and c) dynamic, based on the population sizes and distributions predicted by the LPX Bern model output data. 
Input data is represented in s form of 221 time points spaced 100 years (4 generations) apart. Each time point captures N in 1272 demes (53 x 24 grid). 
Input files are called TreeData[scenario]Div[1/scaling].txt. 


## Running the simulations
Run the simulation after specifying all the necessary input files, in terminal (on Mac) for instance as
```
python3 Main.py -pop LPXBernSims/TreeDataStep10Div1000.txt -row 24 -mig 0.0001 -d 100 -gen 25 -odir LPXBernSims/Scenario3Div1000/OUTPUT_DIR -rep 50

```
or, when run on the scenario, where we assumed 2 ancestral populations: 
```
python3 Main.py -pop LPXBernSims/TreeDataStep10Div1000.txt -row 24 -mig 0.0001 -d 100 -gen 25 --ancpop_list LPXBernSims/ancestralPop.txt -odir LPXBernSims/Scenario4Div1000/OUTPUT_DIR -rep 50
```

All output files are stored in specified output directory, or default 'OUTPUT' directory is created. File 'Output.txt' in the output directory contains all the input information used for the simulations, that are identical. Files 'CoalTimesI.txt' contain coalescence times for I-th simulation replicate.   


## Analysing outputs


Open jupyter notebook AnalysingOutput.ipynb file. Specify the input files, the name of the files containing coalescence times (without the replicate number and txt extension) and their number.

Run the second cell to import relevant packages and define functions.

Run the third cell to load the files and create a matrix with mean coalescence times. This step may take a while, if many replicates must be loaded. A file named with the name of the time-containing files with extension MEAN.txt is created.  

Finally, run the fourth cell to calculate mean total coalescence time, mean within deme and between deme coalescence times, Fst and mean coalescence times for samples taken from different (Manhattan) distance classes. A heatmap with diversity (approximated by within deme coalescence times) and isolation by distance plots are created.
