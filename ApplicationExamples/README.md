# Coalescence simulations on LPX Bern model data

## Input data


## Running the simulations
Run the simulation after specifying all the necessary input files, in terminal (on Mac) for instance as
```
python3 Main.py -pop example_simple_exp_lin_inc_N_10.txt -row 5 -mig example_simple_exp_mig_list0.1.txt -sam example_simple_exp_sample_list.txt -d 20 -gen 2 -rep 1
```
or, when run on the specific example provided on github
```
python3 Main.py -pop Examples/Increasing_population/example_simple_exp_lin_inc_N_10.txt -row 5 -mig Examples/Increasing_population/example_simple_exp_mig_list0.1.txt -sam Examples/Increasing_population/example_simple_exp_sample_list.txt -d 20 -gen 2 -odir Examples/Increasing_population/OUTPUT_DIR -rep 1
```

You can also specify additional parameters, such as printing out demography debugger by msprime by setting
```
-pdeb TRUE
```
sett a seed number by adding
```
-seed INT
```
or specify output file by
```
-odir OUTPUT_DIR_NAME [string ]
```

Multiple parallel simulations can be run in parallel by:

```
for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize1/example_static_N_100.txt -row 5 -mig Examples/GridSize1/example_mig_list0.1.txt -sam Examples/GridSize1/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize1/OUTPUT_DIR  -rep 100; done

```


or in series, for timing:
```
time for i in {1..100} ; do python3 Main.py -pop Examples/Increasing_population/example_simple_exp_lin_inc_N_10.txt -row 5 -mig Examples/Increasing_population/example_simple_exp_mig_list0.1.txt -sam Examples/Increasing_population/example_simple_exp_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Increasing_population/OUTPUT_DIR  -rep $i done
```

All output files are stored in specified output directory, or default 'OUTPUT' directory is created. File 'Output.txt' in the output directory contains all the input information used for the simulations, that are identical. Files 'CoalTimesI.txt' contain coalescence times for I-th simulation replicate.   


## Analysing outputs




Open jupyter notebook AnalysingOutput.ipynb file. Specify the input files, the name of the files containing coalescence times (without the replicate number and txt extension) and their number.

Run the second cell to import relevant packages and define functions.

Run the third cell to load the files and create a matrix with mean coalescence times. This step may take a while, if many replicates must be loaded. A file named with the name of the time-containing files with extension MEAN.txt is created.  

Finally, run the fourth cell to calculate mean total coalescence time, mean within deme and between deme coalescence times, Fst and mean coalescence times for samples taken from different (Manhattan) distance classes. A heatmap with diversity (approximated by within deme coalescence times) and isolation by distance plots are created.
