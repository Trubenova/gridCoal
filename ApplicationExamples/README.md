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


Open jupyter notebook AnalysingOutput.ipynb file. 
Run the first two cells to import relevant packages and define functions.

Specify the input files, the name of the files containing coalescence times (without the replicate number and txt extension) and their number.
```
demography_file='LPXBernSims/TreeDataStep10Div1000.txt'
prefix='LPXBernSims/Scenario3Div1000/'
time_files=prefix+'OUTPUT_DIR/CoalTimes1Rep'
sample_file=prefix+'OUTPUT_DIR/Sample_list.txt' #this we don't have yet
rows=24                     # specify the number of rows that were used to run the simulation
file_number=50              # specify the number of output files that you want to analyse. Files from 0 to the specified number will be analyzed. 
t=np.arange(-22000,1,100)   # time to be used for plotting
```

Run the third cell to load the files and create a matrix with mean coalescence times. This step may take a while, if many replicates must be loaded. A file named with the name of the time-containing files with extension MEAN.txt is created.  

```
my_demography=np.loadtxt(demography_file)
my_samples=np.loadtxt(sample_file) 
my_samples=my_samples.astype(int)
mean_coal_times=get_mean_times(time_files, file_number)
```

The following code snippets are used to create various figures. 

### Plotting tree distribution (final - present - population sizes)
```
plt.figure(figsize=(10.6,4.8), dpi=80)
[T,map_size]=np.shape(my_demography)
cols=int(map_size/rows)
final_map_2D=np.reshape(my_demography[-1,:], [rows, cols])
plt.pcolor(final_map_2D)
plt.colorbar()
plt.savefig(prefix+'Demography.png')
plt.show()
```
### Plotting coalescent times
```
[mean_total, mean_within, mean_between] = get_mean_partial_times(mean_coal_times) 
my_fst=calculate_fst(mean_total, mean_within)
plt.figure(figsize=(10.6,4.8), dpi=80)
plot_within_ctime(my_samples,cols, rows, mean_coal_times)
plt.savefig(prefix+'CoalTimesMap.png')
```
### Plotting isolation by distance patterns



```
name=prefix
plt.figure(figsize=(10.6,4.8), dpi=80)
[T,map_size]=np.shape(my_demography)
cols=int(map_size/rows)

final_map_2D=np.reshape(my_demography[-1,:], [rows, cols])
plt.pcolor(final_map_2D)
plt.colorbar()
plt.savefig(name+'Demography.png')
plt.show()

[mean_total, mean_within, mean_between] = get_mean_partial_times(mean_coal_times) 
print (mean_total, mean_within, mean_between)
my_fst=calculate_fst(mean_total, mean_within)
plt.figure(figsize=(10.6,4.8), dpi=80)
plot_within_ctime(my_samples,cols, rows, mean_coal_times)
plt.savefig(name+'CoalTimesMap.png')

plt.figure(figsize=(10.6,4.8), dpi=80)
box_means=make_ibd_plots(my_samples, mean_coal_times, rows, cols)
plt.savefig(name+'IBD.png')
print ('Mean total coalescence time', mean_total)
print ('Mean within deme coalescence time', mean_within)
print ('Mean between deme coalescence time', mean_between)
print ('Fst', my_fst)
print ('Mean times in distance classes', box_means)
```



Finally, run the fourth cell to calculate mean total coalescence time, mean within deme and between deme coalescence times, Fst and mean coalescence times for samples taken from different (Manhattan) distance classes. A heatmap with diversity (approximated by within deme coalescence times) and isolation by distance plots are created.
