#!/usr/bin/env python3
description = '''Simulate a coalescent within the history of Abies expansion.'''


import msprime
import numpy as np
import os
import timeit

import argparse

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--tree_nums", "-t", type=str, dest="tree_nums_file", required=True,
                    help="name of file to load tree number history from")
parser.add_argument("--sample_coords", "-s", type=str, dest="sample_coords_file", required=True,
                    help="name of file containing locations of sampling coordinates")
parser.add_argument("--migr_mat", "-m", type=str, dest="migr_mat_file", required=True,
                    help="name of file containing migration matrix")
parser.add_argument("--adj_mat", "-am", type=str, dest="adj_mat_file", required=True,
                    help="name of file containing adjacency matrix")
parser.add_argument("--basedir", "-o", type=str, dest="basedir", required=True,
                    help="name of directory to save output files to.")
parser.add_argument("--dt", "-d", type=float, dest="dt", required=True,
                    help="time interval between steps")
parser.add_argument("--recomb", "-r", type=float, dest="recomb", required=True,
                    help="total recombination rate")
parser.add_argument("--mut", "-u", type=float, dest="mut", required=True,
                    help="total mutation rate")
parser.add_argument("--num_reps", "-n", type=float, dest="num_reps", default=1,
                    help="number of replicates")
parser.add_argument("--logfile", "-g", type=float, dest="logfile",
                    help="name of log file")
parser.add_argument("--generation_time", "-gt", type=int, dest="gen_time",
                    help="generation time")
parser.add_argument("--row_number", "-rn", type=int, dest="row_num",
                    help="number of rows in the originical grid")
parser.add_argument("--ancpop_size", "-aps",  type=int,nargs='+',dest="ancpop_size", required=True,
                    help="list of the sizes of the ancestral populations")
parser.add_argument("--ancpop_list", "-apl", type=str, dest="ancpop_list", required=True,
                    help="name of file to load a list with the destination ancestral population for each gridcell")
#parser.add_argument("--ancpop_loc", "-apl",type=int,nargs='+', dest="ancpop_loc", required=True,
#                    help="list of the locations of the ancestral populations")
#parser.add_argument("--ancpop_time", "-apt", type=int, dest="ancpop_time", required=True,
#                    help="the time while several ancestral populations exist")
parser.add_argument("--scale", "-sc", type=int, dest="sc",
                    help="scaling factor")
parser.add_argument("--serial", "-serial", type=int, dest="serial",
                    help="serial number")





############################### PARAMETERS, FUNCTIONS #################################
#######################################################################################

args = parser.parse_args()
m=args.migr_mat_file
am=args.adj_mat_file
ps=args.tree_nums_file
serial=args.serial
scale=args.sc


a_s= []
for x in args.ancpop_size:
    a_s.append(int(x/scale))

anc_num=len(args.ancpop_size)
ancpop_list=np.loadtxt(args.ancpop_list)

# timestep
dt=args.dt
#dt=100
print(dt)
# dealing with empty cells
MIN_SIZE = 1e-10



def read_tree_nums():
    # e.g. tree_num_epsg3035.tsv
    trees = np.loadtxt("{}.tsv".format(ps))
    return trees

def read_migration_matrix(fname):
    f = open(fname, "r")
    header = f.readline().split()
    f.close()
    i = np.loadtxt(fname, skiprows=1, usecols=header.index('"i"'), dtype='int16')
    j = np.loadtxt(fname, skiprows=1, usecols=header.index('"j"'), dtype='int16')
    x = np.loadtxt(fname, skiprows=1, usecols=header.index('"x"'))
    if (x < 0.0).any():
        raise ValueError("Migration rates must be positive.")
    ids = list(set(list(set(i)) + list(set(j))))
    ids.sort()
    elem=ids[-1]+1
    M = np.zeros((elem,elem))#np.zeros((len(ids), len(ids)))
    for ii, jj, xx in zip(i, j, x):
        M[ii, jj] = xx
    # msprime wants diag elements to be zero
    for ii in ids:
        M[ii, ii] = 0.0
    return M

def read_coords():
    # e.g. coords_epsg3035.tsv
    locfile = open(args.sample_coords_file, "r")
    header = locfile.readline().split()
    assert(header == ['"X"', '"Y"', '"pop"', '"id"', '"cell"', '"msp_id"'])
    poplist = []
    idlist = []
    celllist = []
    for line in locfile:
        x, y, pop, id, cell, msp = [u.strip('"\n') for u in line.split()]
        poplist.append(pop)
        idlist.append(id)
        celllist.append(int(cell))
    return poplist, idlist, celllist


#######################################################################################
#######################################################################################

# read the migration rates and the population sizes
tree_num = read_tree_nums()
M = read_migration_matrix("M{}.tsv".format(m))
ADJ = read_migration_matrix("{}.tsv".format(am))



n = M.shape[0]


# the extended migration includes the ancestral populations, but the migration rate
# to and from these cells are 0. This is only used for the simulate() function
# in every other case the original migration matrix is used together with the population sizes
# without the ancestral ones.


M_ext=M
for i in range(anc_num):
    M_ext=np.vstack([M_ext,np.zeros(n)])

for i in range(anc_num):
    M_ext=np.hstack((M_ext,np.zeros((M_ext.shape[0],1),dtype=M_ext.dtype)))

n_ext=M_ext.shape[0]

# n is the number of demes in total on the grid
ngens = tree_num.shape[0]


file1 = "ET_pairwise.mig{}.ancpopnum{}.ancpopsize{}.dem_{}.dt{}.serial{}.tsv".format(m,anc_num,a_s[0],ps,dt,args.serial)





###################### POPULATION CONFIGURATIONS, SAMPLES ###############################
#########################################################################################




# to include the extra ancestral populations, we need to add extra 0s to the vector
# ndip and and nsamples, however these will not affect the iteration or any other calculations


poplist, idlist, celllist = read_coords()

nsamples = [0 for _ in range(n_ext)]
ndip = [0 for _ in range(n_ext)]
for k in range(n):
	nsamples[k] = 2*sum([u == k for u in celllist])
	ndip[k]=sum([u == k for u in celllist])
number_of_samples=sum([nsamples[k] for k in range(n_ext)])




sampled_demes=list(set(celllist))
sampled_demes.sort()
samp=list(set(tree_num[-1][sampled_demes]))
samp.sort()


assert tree_num.shape[1]%args.row_num == 0, "Oh no! Row numbers are not correct!"
assert n == M.shape[1], "Oh no! The migration matrix is not a square matrix!"
assert n == tree_num.shape[1], "Oh no, the number of demes is different in M and tree_num!"
assert n_ext == M_ext.shape[1], "Oh no! The extended migration matrix is not a square matrix!"
assert n_ext == n+anc_num, "Oh no! The extended population is not the sum of n and anc pop number!"
assert samp[0]>0, "Oh-oooo.... We sample from cells that are empty!"



sd=len(sampled_demes)



subpops=[]
counter=0
for i in range(len(ndip)):
    if ndip[i]!=0:
        elem=[]
        for i in range(counter, counter+ndip[i]):
            elem.append(i)
        subpops.append(elem)
        counter=elem[-1]+1

population_configurations=[]
demographic_events = []


tree_origi=[0 for _ in range(n_ext)]
tree_origi_last=[0 for _ in range(n_ext)]

for k in range(n):
	tree_origi[k] = tree_num[ngens-1,k]
	tree_origi_last[k] = tree_num[ngens-2,k]






population_configurations = [
    msprime.PopulationConfiguration(sample_size=ndip[k],
        initial_size = max(MIN_SIZE, tree_origi[k]),
        growth_rate=0) for k in range(n_ext) ]


#########################################################################################
#########################################################################################



######### DEFINE DEMOGRAPHIC EVENTS THROUGH THE POPULATION HISTORY ######################
#########################################################################################



##################################### FIRST TIMESTEP ####################################
#########################################################################################


# *_ext are only to account for the ancestral populations, have no role until the final timestep
# N is a matrix: how many individuals in pop. i have parents in j
# if a deme has individuals now but in the previous timestep was empty, we need a mass migration (backward in time):
# mass migration from this cell to the neighbors, proportions are weighted by the sizes of the population sizes in the previous generation
# pT, popT is to decide about these ratios when mass migration is needed
# when mass migration happens, a population goes empty, so the order of events (may?) matter
# if there is mass migration from i to j, we first set the new size of j, mass migrate the lineages, set the size of i to MIN_SIZE
# migration matrix (i,j): how many individuals come from cell j to i, divided by the size of deme i

N=(M.T*tree_num[ngens-2])
mig=[]
N_ext=(M_ext.T*tree_origi_last)
mig_ext=[]


for i in range(n):
    mig.append([tree_num[ngens-1,i]])

migrate=M[0,1]
if migrate==0:
    pT=(M.T*tree_num[ngens-2])
else:
    pT=(M.T*tree_num[ngens-2]/migrate)

popT=pT.sum(axis=1)
#popT[np.isnan(popT)] = 0.0
#print(popT)

for i in range(n):
    if tree_num[ngens-2,i]==0 and tree_num[ngens-1,i]!=0:
	    for j in range(n):
	        if ADJ[i,j]!=0:
	            demographic_events.append(msprime.PopulationParametersChange(
			        time=dt/args.gen_time, initial_size=max(MIN_SIZE, tree_num[ngens-2,j]),
			            population_id=j, growth_rate=0))

	            if popT[i]==0.0:
	                prop=0
	            else:
	                prop=tree_num[ngens-2,j]/popT[i]

	            demographic_events.append(msprime.MassMigration(
                    time=dt/args.gen_time, source=i, destination=j, proportion=min(1,prop)))

	            demographic_events.append(msprime.PopulationParametersChange(
			        time=dt/args.gen_time, initial_size=max(MIN_SIZE, tree_num[ngens-2,i]),
			            population_id=i, growth_rate=0))
    else:
	    demographic_events.append(msprime.PopulationParametersChange(
			        time=dt/args.gen_time, initial_size=max(MIN_SIZE, tree_num[ngens-2,i]),
			            population_id=i, growth_rate=0))


BM=N/mig
BM[np.isnan(BM)] = 0.0
from numpy import inf
BM[BM == inf] = 0.0


#same for _ext

for i in range(n_ext):
    mig_ext.append([tree_origi[i]])

BM_ext=N_ext/mig_ext
BM_ext[np.isnan(BM_ext)] = 0.0

from numpy import inf
BM_ext[BM_ext == inf] = 0.0


#this is the starting migration matrix, if in the previous generation a cell was empty, the value should be 0
bwM_ext=BM_ext

for i in range(n):
    for j in range(n):
        if tree_num[ngens-2,j]==0:
            bwM_ext[i,j]=0


################################### ITERATION ###########################################
#########################################################################################


# we need to go through all the timesteps, the first (above) and the last are treated separately,
# between there is a while loop:
# x is the tree numbers for all the cells now, x_next is the previous generation
# migration matrix is made the same way as before
# mass migrations and population size updates work the same way
# population size may get updated twice...

t= ngens-2
x_last=tree_num[t]


BM_last=BM
x=tree_num[t]


while t>=1:
    print(t)
    x_next=tree_num[t-1]
    t_ago = dt * (ngens - t)   #dt * (ngens -1- t)


    N=(M.T*x_next)
    mig=[]
    for i in range(n):
        mig.append([x[i]])


    migrate=M[0,1]

    if migrate==0:
        pT=(M.T*x_next)
    else:
        pT=(M.T*x_next/migrate)

    popT=pT.sum(axis=1)


    for i in range(n):
        if x_next[i]==0 and x[i]!=0:
	        for j in range(n):
	            if ADJ[i,j]!=0:
	                demographic_events.append(msprime.PopulationParametersChange(
			            time=t_ago/args.gen_time, initial_size=max(MIN_SIZE, x_next[j]),
			            population_id=j, growth_rate=0))

	                if popT[i]==0.0:
	                    prop=0
	                else:
	                    prop=x_next[j]/popT[i]

	                demographic_events.append(msprime.MassMigration(
                        time=t_ago/args.gen_time, source=i, destination=j, proportion=min(1, prop)))

	                demographic_events.append(msprime.PopulationParametersChange(
			            time=t_ago/args.gen_time, initial_size=max(MIN_SIZE, x_next[i]),
			            population_id=i, growth_rate=0))
        else:
	        demographic_events.append(msprime.PopulationParametersChange(
			    time=t_ago/args.gen_time, initial_size=max(MIN_SIZE, x_next[i]),
			    population_id=i, growth_rate=0))


    BM=N/mig
    BM[np.isnan(BM)] = 0.0
    BM[BM == inf] = 0.0


    for i in range(n):
        for j in range(n):
            if x[j]!=0 and BM[i,j]!=BM_last[i,j]:
                demographic_events.append(msprime.MigrationRateChange(
                     time=t_ago/args.gen_time, rate=BM[i,j], matrix_index=(i, j)))
            if x[j]==0 and BM[i,j]!=BM_last[i,j]:
                demographic_events.append(msprime.MigrationRateChange(
                     time=t_ago/args.gen_time, rate=0, matrix_index=(i, j)))

    x=x_next
    BM_last=BM
    t-=1

################################# LAST TIMESTEP #########################################
#########################################################################################



N=(M.T*x)
mig=[]
for i in range(n):
    mig.append([x[i]])

if migrate==0:
    pT=(M.T*x_next)
else:
    pT=(M.T*x_next/migrate)

popT=pT.sum(axis=1)
BM=N/mig
BM[np.isnan(BM)] = 0.0
BM[BM == inf] = 0.0
for i in range(n):
    for j in range(n):
        if BM[i,j]!=BM_last[i,j]:
            demographic_events.append(msprime.MigrationRateChange(
                    time=dt * (ngens -1- t)/args.gen_time, rate=BM[i,j], matrix_index=(i, j)))




###################### MASS MIGRATION TO ANCESTRAL POPULATIONS ##########################
#########################################################################################


# first grow the population sizes of the ancestral ones to the desired number
# find the closest ancestral population for each lineage and then mass-migrate the lineage there


for i in range(anc_num):
    demographic_events.append(msprime.PopulationParametersChange(
    time=dt * (ngens- t)/args.gen_time, initial_size=a_s[i],
	population_id=n+i, growth_rate=0))


for i in range(n):
    demographic_events.append(msprime.MassMigration(
        time=dt * (ngens - t)/args.gen_time, source=i, destination=n+int(ancpop_list[[i]])-1, proportion=1.0))


## migration rates are 0s for the grid, and for the ancestral population it's something small

for i in range(n):
    for j in range(n):
        if BM[i,j]!=0.0:
            demographic_events.append(msprime.MigrationRateChange(
                time=dt * (ngens - t)/args.gen_time, rate=0.0, matrix_index=(i, j)))


for i in range(anc_num):
    for j in range(anc_num):
        if i!=j:
            demographic_events.append(msprime.MigrationRateChange(
                    time=dt * (ngens - t)/args.gen_time, rate=1e-8, matrix_index=(n+i, n+j)))


##################################### SIMULATION ########################################
#########################################################################################

print("sim start")
#dd = msprime.DemographyDebugger(
#    population_configurations=population_configurations,
#    migration_matrix=bwM_ext,#
#    demographic_events=demographic_events)

#dd.print_history()



replicates = msprime.simulate(
	population_configurations=population_configurations,
	demographic_events=demographic_events,
	migration_matrix=bwM_ext,
	recombination_rate=args.recomb,
	mutation_rate=args.mut)


########################CALCULATE EXPECTED PAIRWISE COALESCENCE TIMES####################
#########################################################################################

tree = replicates.first()
aa=replicates.simplify().first()
lP=len(subpops)




result=[]
result=np.zeros((lP,lP))

for i in range(lP):
    for j in range(lP):
        if i != j:
            result[i,j]=(aa.tmrca(subpops[i][1],subpops[j][1])+aa.tmrca(subpops[i][0],subpops[j][0])+aa.tmrca(subpops[i][0],subpops[j][1])+aa.tmrca(subpops[i][1],subpops[j][0]))/4
        else:
            result[i,j]=aa.tmrca(subpops[i][0],subpops[j][1])


np.savetxt(file1, result)
