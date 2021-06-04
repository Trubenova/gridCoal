# Simulating increase in population sizes

In this worked out example, we create simple input data representing an expansion of a fictitious species over T=10 time points spaced dt=20 years apart, with generation time of gt=2 years, on a rectangular grid of 5 rows x 4 columns . At the final time point, all cells are occupied, with  population sizes drawn from Poisson distribution with mean size of N= 10.
Migration rate between neighbouring cells is 0.1.  

We sample half of the cells (s=0.5).

Then we run 100 coalescence simulations using the main gridCoal simulator (Main.py), generating 100 output files with coalescence times, one with a summary of inputs, as well as very detailed demography debugger produced by msprime.

Then we analyse the outputs, calculating global Fst anf F*.

## Generating inputs
Open jupyter notebook GeneratingInputs.ipynb and define the parameters in the second cell.
```python
rows = 5    # number of rows
cols = 4    # number of columns
T = 10      # number of time steps
gt = 2      # generation time
dt = 20     # time between defined time steps
N = 10      # average population size
mu = 0.1    # migration rate between neighbouring demes
s = 0.5     # s is coverage - fraction of sampled grid cells. Note that using this function may include
            # some that are empty in the input data, which will cause error. Check before submitting for simulations.
seed = 10   #initialize random seed if you want
batch_name = 'example_simple_exp_'  #prefix used for input data associated with this simulation
```

then run the following cell to define the functions.

The next cell generated input files, that are saved as txt files with prefix specified in 'batch_name':

```python
np.random.seed(seed)
final_map = generate_final_map(rows, cols, N)   #generate final map
my_demography = generate_lin_increase_data(final_map, T, batch_name)  #generate linear increasing pop sizes leading to final map
mig_list = generate_migration_list(rows, cols, mu, batch_name)
sample_list = generate_sample_list(rows, cols,0.5, batch_name)
ancpop_list = generate_ancestral_pop(rows, cols, batch_name=batch_name )
```
Finally, the last cell print out the inputs on the screen.


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

Multiple parallel simulations can be run by:
```   
for i in {1..100} ; do python3 Main.py -pop example_simple_exp_lin_inc_N_10.txt -row 5
-mig example_simple_exp_mig_list0.1.txt -sam example_simple_exp_sample_list.txt -d 20
-gen 2 -seed 1 -rep $i & done
```

```
for i in {1..100} ; do python3 Main.py -pop Examples/Increasing_population/example_simple_exp_lin_inc_N_10.txt -row 5 -mig Examples/Increasing_population/example_simple_exp_mig_list0.1.txt -sam Examples/Increasing_population/example_simple_exp_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Increasing_population/OUTPUT_DIR  -rep $i & done
```


All output files are stored in specified output directory, or default 'OUTPUT' directory is created. File 'Output.txt' in the output directory contains all the input information used for the simulations, that are identical. Files 'CoalTimesI.txt' contain coalescence times for I-th simulation replicate.   

## Analysing outputs

Open jupyter notebook AnalysingOutput.ipynb file. Specify the input files, the name of the files containing coalescence times (without the replicate number and txt extension) and their number.

Run the second cell to import relevant packages and define functions.

Run the third cell to load the files and create a matrix with mean coalescence times. This step may take a while, if many replicates must be loaded. A file named with the name of the time-containing files with extension MEAN.txt is created.  

Finally, run the fourth cell to calculate mean total coalescence time, mean within deme and between deme coalescence times, Fst and mean coalescence times for samples taken from different (Manhattan) distance classes. A heatmap with diversity (approximated by within deme coalescence times) and isolation by distance plots are created.
