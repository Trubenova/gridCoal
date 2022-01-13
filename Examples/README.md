SampledFraction07# Simulating increase in population sizes

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
T  = 10      # number of time steps
gt = 2      # generation time
dt = 20     # time between defined time steps
N  = 10      # average population size
mu = 0.1    # migration rate between neighbouring demes
s  = 0.5     # s is coverage - fraction of sampled grid cells. Note that using this function may include
            # some that are empty in the input data, which will cause error. Check before submitting for simulations.
seed = 10   #initialize random seed if you want
dir_name='Examples/Static_population/' #set name of the directory
if not os.path.exists(dir_name):
    os.makedirs(dir_name)
    print("Directory " , dir_name ,  " Created ")
batch_name = dir_name+'example_simple_exp_'  #prefix used for input data associated with this simulation
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

Multiple parallel simulations can be run in parallel by:

```
for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize1/example_static_N_100.txt -row 5 -mig Examples/GridSize1/example_mig_list0.1.txt -sam Examples/GridSize1/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize1/OUTPUT_DIR  -rep 100; done

```
for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize1/example_static_N_100.txt -row 5 -mig Examples/GridSize1/example_mig_list0.1.txt -sam Examples/GridSize1/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize1/OUTPUT_DIR  -rep 100; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize2/example_static_N_100.txt -row 10 -mig Examples/GridSize2/example_mig_list0.1.txt -sam Examples/GridSize2/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize2/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize3/example_static_N_100.txt -row 15 -mig Examples/GridSize3/example_mig_list0.1.txt -sam Examples/GridSize3/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize3/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize4/example_static_N_100.txt -row 20 -mig Examples/GridSize4/example_mig_list0.1.txt -sam Examples/GridSize4/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize4/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/GridSize5/example_static_N_100.txt -row 25 -mig Examples/GridSize5/example_mig_list0.1.txt -sam Examples/GridSize5/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GridSize5/OUTPUT_DIR  -rep 100 -ser $i; done

---------------------
for i in {1..5} ; do  time python3 Main.py -pop Examples/SampledFraction01/example_static_N_100.txt -row 15 -mig Examples/SampledFraction01/example_mig_list0.1.txt -sam Examples/SampledFraction01/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/SampledFraction01/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/SampledFraction03/example_static_N_100.txt -row 15 -mig Examples/SampledFraction03/example_mig_list0.1.txt -sam Examples/SampledFraction03/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/SampledFraction03/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/SampledFraction05/example_static_N_100.txt -row 15 -mig Examples/SampledFraction05/example_mig_list0.1.txt -sam Examples/SampledFraction05/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/SampledFraction05/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/SampledFraction07/example_static_N_100.txt -row 15 -mig Examples/SampledFraction07/example_mig_list0.1.txt -sam Examples/SampledFraction07/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/SampledFraction07/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/SampledFraction1/example_static_N_100.txt -row 15 -mig Examples/SampledFraction1/example_mig_list0.1.txt -sam Examples/SampledFraction1/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/SampledFraction1/OUTPUT_DIR  -rep 100 -ser $i; done
------------

for i in {1..5} ; do  time python3 Main.py -pop Examples/MeanPopSize10/example_static_N_10.txt -row 15 -mig Examples/MeanPopSize10/example_mig_list0.1.txt -sam Examples/MeanPopSize10/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MeanPopSize10/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MeanPopSize100/example_static_N_100.txt -row 15 -mig Examples/MeanPopSize100/example_mig_list0.1.txt -sam Examples/MeanPopSize100/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MeanPopSize100/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MeanPopSize500/example_static_N_500.txt -row 15 -mig Examples/MeanPopSize500/example_mig_list0.1.txt -sam Examples/MeanPopSize500/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MeanPopSize500/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MeanPopSize1000/example_static_N_1000.txt -row 15 -mig Examples/MeanPopSize1000/example_mig_list0.1.txt -sam Examples/MeanPopSize1000/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MeanPopSize1000/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MeanPopSize2000/example_static_N_2000.txt -row 15 -mig Examples/MeanPopSize2000/example_mig_list0.1.txt -sam Examples/MeanPopSize2000/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MeanPopSize2000/OUTPUT_DIR  -rep 100 -ser $i; done

----------

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimePoints5/example_static_N_100.txt -row 15 -mig Examples/TimePoints5/example_mig_list0.1.txt -sam Examples/TimePoints5/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/TimePoints5/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimePoints10/example_static_N_100.txt -row 15 -mig Examples/TimePoints10/example_mig_list0.1.txt -sam Examples/TimePoints10/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/TimePoints10/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimePoints20/example_static_N_100.txt -row 15 -mig Examples/TimePoints20/example_mig_list0.1.txt -sam Examples/TimePoints20/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/TimePoints20/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimePoints50/example_static_N_100.txt -row 15 -mig Examples/TimePoints50/example_mig_list0.1.txt -sam Examples/TimePoints50/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/TimePoints50/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimePoints100/example_static_N_100.txt -row 15 -mig Examples/TimePoints100/example_mig_list0.1.txt -sam Examples/TimePoints100/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/TimePoints100/OUTPUT_DIR  -rep 100 -ser $i; done
--------

for i in {1..5} ; do  time python3 Main.py -pop Examples/GenTime1/example_static_N_100.txt -row 15 -mig Examples/GenTime1/example_mig_list0.1.txt -sam Examples/GenTime1/example_sample_list.txt -d 200 -gen 1 -seed 1 -odir Examples/GenTime1/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/GenTime5/example_static_N_100.txt -row 15 -mig Examples/GenTime5/example_mig_list0.1.txt -sam Examples/GenTime5/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/GenTime5/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/GenTime10/example_static_N_100.txt -row 15 -mig Examples/GenTime10/example_mig_list0.1.txt -sam Examples/GenTime10/example_sample_list.txt -d 200 -gen 10 -seed 1 -odir Examples/GenTime10/OUTPUT_DIR  -rep 100 -ser $i; done


for i in {1..5} ; do  time python3 Main.py -pop Examples/GenTime50/example_static_N_100.txt -row 15 -mig Examples/GenTime50/example_mig_list0.1.txt -sam Examples/GenTime50/example_sample_list.txt -d 200 -gen 50 -seed 1 -odir Examples/GenTime50/OUTPUT_DIR  -rep 100 -ser $i; done


for i in {1..5} ; do  time python3 Main.py -pop Examples/GenTime100/example_static_N_100.txt -row 15 -mig Examples/GenTime100/example_mig_list0.1.txt -sam Examples/GenTime100/example_sample_list.txt -d 200 -gen 100 -seed 1 -odir Examples/GenTime100/OUTPUT_DIR  -rep 100 -ser $i; done



-----------


for i in {1..5} ; do  time python3 Main.py -pop Examples/TimeStep10/example_static_N_100.txt -row 15 -mig Examples/TimeStep10/example_mig_list0.1.txt -sam Examples/TimeStep10/example_sample_list.txt -d 10 -gen 5 -seed 1 -odir Examples/TimeStep10/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimeStep50/example_static_N_100.txt -row 15 -mig Examples/TimeStep50/example_mig_list0.1.txt -sam Examples/TimeStep50/example_sample_list.txt -d 50 -gen 5 -seed 1 -odir Examples/TimeStep50/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimeStep100/example_static_N_100.txt -row 15 -mig Examples/TimeStep100/example_mig_list0.1.txt -sam Examples/TimeStep100/example_sample_list.txt -d 100 -gen 5 -seed 1 -odir Examples/TimeStep100/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimeStep200/example_static_N_100.txt -row 15 -mig Examples/TimeStep200/example_mig_list0.1.txt -sam Examples/TimeStep200/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/TimeStep200/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/TimeStep500/example_static_N_100.txt -row 15 -mig Examples/TimeStep500/example_mig_list0.1.txt -sam Examples/TimeStep500/example_sample_list.txt -d 500 -gen 5 -seed 1 -odir Examples/TimeStep500/OUTPUT_DIR  -rep 100 -ser $i; done

----------

for i in {1..5} ; do  time python3 Main.py -pop Examples/MigRate0/example_static_N_100.txt -row 15 -mig Examples/MigRate0/example_mig_list0.0.txt -sam Examples/MigRate0/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MigRate0/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MigRate0001/example_static_N_100.txt -row 15 -mig Examples/MigRate0001/example_mig_list0.001.txt -sam Examples/MigRate0001/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MigRate0001/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MigRate001/example_static_N_100.txt -row 15 -mig Examples/MigRate001/example_mig_list0.01.txt -sam Examples/MigRate001/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MigRate001/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MigRate01/example_static_N_100.txt -row 15 -mig Examples/MigRate01/example_mig_list0.1.txt -sam Examples/MigRate01/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MigRate01/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/MigRate1/example_static_N_100.txt -row 15 -mig Examples/MigRate1/example_mig_list1.txt -sam Examples/MigRate1/example_sample_list.txt -d 200 -gen 5 -seed 1 -odir Examples/MigRate1/OUTPUT_DIR  -rep 100 -ser $i; done


```   
for i in {1..100} ; do python3 Main.py -pop example_simple_exp_lin_inc_N_10.txt -row 5
-mig example_simple_exp_mig_list0.1.txt -sam example_simple_exp_sample_list.txt -d 20
-gen 2 -seed 1 -rep $i & done


for i in {1..5} ; do  time python3 Main.py -pop Examples/Splatche2/example_static_N_2.txt -row 50 -mig Examples/Splatche2/example_mig_list0.1.txt -sam Examples/Splatche2/example_sample_list.txt -d 200 -gen 100 -seed 1 -odir Examples/Splatche2/OUTPUT_DIR  -rep 100 -ser $i; done

for i in {1..5} ; do  time python3 Main.py -pop Examples/Splatche3/example_static_N_0.txt -row 50 -mig Examples/Splatche3/example_mig_list0.1.txt -d 200 -gen 2 -seed 1 -odir Examples/Splatche3/OUTPUT_DIR  -rep 10 -ser $i; done




```
or in series, for timing:
```
time for i in {1..100} ; do python3 Main.py -pop Examples/Increasing_population/example_simple_exp_lin_inc_N_10.txt -row 5 -mig Examples/Increasing_population/example_simple_exp_mig_list0.1.txt -sam Examples/Increasing_population/example_simple_exp_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Increasing_population/OUTPUT_DIR  -rep $i done
```

time for i in {1..100} ; do python3 Main.py -pop Examples/Static_population/example_static_static_N_10.txt -row 5 -mig Examples/Static_population/example_static_mig_list0.1.txt -sam Examples/Static_population/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Static_population/OUTPUT_DIR  -rep $i; done
```
time for i in {1..100} ; do python3 Main.py -pop Examples/Static_populationL3x/example_static_static_N_100.txt -row 15 -mig Examples/Static_populationL3x/example_static_mig_list0.1.txt -sam Examples/Static_populationL3x/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Static_populationL3x/OUTPUT_DIR  -rep $i; done
```

for j in {1..5} ; do time for i in {1..100} ; do python3 Main.py -pop Examples/Static_populationMig001/example_static_static_N_100.txt -row 15 -mig Examples/SpeedTest/example_static_mig_list0.001.txt -sam Examples/SpeedTest/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/SpeedTest/OUTPUT_DIR  -rep $i; done; done



for j in {1..5} ; do time for i in {1..100} ; do python3 Main.py -pop Examples/SpeedTest/example_static_static_N_100.txt -row 15 -mig Examples/Static_populationL3xDt200/example_static_mig_list0.1.txt -sam Examples/Static_populationL3xDt200/example_static_sample_list.txt -d 200 -gen 2 -seed 1 -odir Examples/Static_populationL3xDt200/OUTPUT_DIR  -rep $i; done; done

for j in {1..5} ; do time for i in {1..100} ; do python3 Main.py -pop Examples/Static_populationL5x/example_static_static_N_100.txt -row 25 -mig Examples/Static_populationL5x/example_static_mig_list0.1.txt -sam Examples/Static_populationL5x/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Static_populationL5x/OUTPUT_DIR  -rep $i; done; done


for j in {1..5} ; do time for i in {1..1} ; do python3 Main.py -pop Examples/Splatche1/example_static_static_N_0.23529.txt -row 25 -mig Examples/Splatche1/example_static_mig_list0.1.txt -sam Examples/Splatche1/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Splatche1/OUTPUT_DIR  -rep $i; done; done

for j in {1..5} ; do time for i in {1..1} ; do  Applications/pypy3.7-v7.3.5-osx64/bin/pypy3 Main.py -pop Examples/Splatche1/example_static_static_N_0.23529.txt -row 25 -mig Examples/Splatche1/example_static_mig_list0.1.txt -sam Examples/Splatche1/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Splatche1/OUTPUT_DIR  -rep $i; done; done



for j in {1..5} ; do time for i in {1..1} ; do ./Main  -pop Examples/Splatche1/example_static_static_N_0.23529.txt -row 25 -mig Examples/Splatche1/example_static_mig_list0.1.txt -sam Examples/Splatche1/example_static_sample_list.txt -d 20 -gen 2 -seed 1 -odir Examples/Splatche1/OUTPUT_DIR  -rep $i; done; done

All output files are stored in specified output directory, or default 'OUTPUT' directory is created. File 'Output.txt' in the output directory contains all the input information used for the simulations, that are identical. Files 'CoalTimesI.txt' contain coalescence times for I-th simulation replicate.   

## Analysing outputs

Open jupyter notebook AnalysingOutput.ipynb file. Specify the input files, the name of the files containing coalescence times (without the replicate number and txt extension) and their number.

Run the second cell to import relevant packages and define functions.

Run the third cell to load the files and create a matrix with mean coalescence times. This step may take a while, if many replicates must be loaded. A file named with the name of the time-containing files with extension MEAN.txt is created.  

Finally, run the fourth cell to calculate mean total coalescence time, mean within deme and between deme coalescence times, Fst and mean coalescence times for samples taken from different (Manhattan) distance classes. A heatmap with diversity (approximated by within deme coalescence times) and isolation by distance plots are created.
