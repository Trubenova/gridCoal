#!/usr/bin/env python3
description = '''Simulate a coalescent within the history of Abies expansion.'''


import msprime
import numpy as np
import os
import argparse
#import timeit


parser = argparse.ArgumentParser(description=description)
parser.add_argument("--tree_nums", "-t", type=str, dest="tree_nums_file", required=True,
                    help="name of file to load tree number history from")
#this loads the file with the demographic history, all cells population sizes in time

parser.add_argument("--migr_mat", "-m", type=str, dest="migr_mat_file", required=True,
                    help="name of file containing migration matrix")
#this loads migration matrix

parser.add_argument("--adj_mat", "-am", type=str, dest="adj_mat_file", required=True,
                    help="name of file containing adjacency matrix")
#this loads adjacency matrix for when we need migration = 0

parser.add_argument("--dt", "-d", type=float, dest="dt", required=True,
                    help="time interval between steps")
#This gives dt- time between two known points


parser.add_argument("--ancpop_list", "-apl", type=str, dest="ancpop_list", required=True,
                    help="name of file to load a list with the destination ancestral population for each gridcell")
#this defines ancestral population for each cell. 

parser.add_argument("--ancpop_size", "-aps",  type=int,nargs='+',dest="ancpop_size", required=True,
                    help="list of the sizes of the ancestral populations")
#this should give sizes of ancestral populations, should be of dimension of the number of anc pop


#############################################  SYS stuff.  check if this is needed ########################################## 
######################################################################################################

parser.add_argument("--serial", "-serial", type=int, dest="serial",
                    help="serial number")
#I think this only assignes simulation serial number so it is saved into a different output file. 

parser.add_argument("--basedir", "-o", type=str, dest="basedir", required=True,
                    help="name of directory to save output files to.")
#this is important!!!

parser.add_argument("--logfile", "-g", type=float, dest="logfile",
                    help="name of log file")
#define log file

#############################  I don't know what this is for or what does it do####################### 
######################################################################################################

parser.add_argument("--generation_time", "-gt", type=int, dest="gen_time",
                    help="generation time")
#generation time, in our case 25. What does it say? can we decrease it?

parser.add_argument("--row_number", "-rn", type=int, dest="row_num",
                    help="number of rows in the originical grid")
#number of rows in the original grid. Just for the first check


###################################  I am not sure if we need this stuff ############################## 
#######################################################################################################

parser.add_argument("--scale", "-sc", type=int, dest="sc",
                    help="scaling factor")
#scaling factor, not sure what for - to use the same numbers? will it round up stuff? 


parser.add_argument("--sample_coords", "-s", type=str, dest="sample_coords_file", required=True,
                    help="name of file containing locations of sampling coordinates")
#this loads sample coordinates -  i am not sure what is this for. 

parser.add_argument("--recomb", "-r", type=float, dest="recomb", required=True,
                    help="total recombination rate")  #this is 0
#this gives recombination rate - must be 0, maybe we should just remove it?!

parser.add_argument("--mut", "-u", type=float, dest="mut", required=True,
                    help="total mutation rate")  #this is 0
#this defines mutation rate, I would suggest removing it. 

parser.add_argument("--num_reps", "-n", type=float, dest="num_reps", default=1,
                    help="number of replicates")  # it is just 1 in our case
#number of replicates - what for? we allways use 1


#parser.add_argument("--ancpop_loc", "-apl",type=int,nargs='+', dest="ancpop_loc", required=True,
#                    help="list of the locations of the ancestral populations")
#parser.add_argument("--ancpop_time", "-apt", type=int, dest="ancpop_time", required=True,
#                    help="the time while several ancestral populations exist")


############################### DEFINING PARAMETERS #################################
#######################################################################################

args = parser.parse_args()
m = args.migr_mat_file
am = args.adj_mat_file
ps = args.tree_nums_file
serial = args.serial
scale = args.sc  #possibly remove this
dt = args.dt  #timestep
genTime=args.gen_time
#dt=100
minSize = 1e-10 # dealing with empty cells, this could be done in the load function?


print('I have: mig file', m)
print('Adj matrix', am)
print('TreeInputfile', ps)
print( 'serial', serial)
print( 'scaling', scale)
print('dt is',  dt)
print( 'min size is', minSize)

##########################################  FUNCTIONS #################################
#######################################################################################




def read_inputs():
    # e.g. tree_num_epsg3035.tsv
    trees = np.loadtxt("{}.tsv".format(ps))
    celllist =np.loadtxt(args.sample_coords_file, dtype='int')
    
    return trees, celllist

def checkParameters():
    assert tree_num.shape[1]%args.row_num == 0, "Oh no! Row numbers are not correct!"  #checks if n could be divided by row numb
    assert n == M.shape[1], "Oh no! The migration matrix is not a square matrix!"
    assert n == tree_num.shape[1], "Oh no, the number of demes is different in M and tree_num!"
    assert n_ext == M_ext.shape[1], "Oh no! The extended migration matrix is not a square matrix!"
    assert n_ext == n+anc_num, "Oh no! The extended population is not the sum of n and anc pop number!"
    assert samp[0]>0, "Oh-no! We sample from cells that are empty!"
    print ('All is fine')




def read_migration_matrix(fname):
    f = open(fname, "r")
    header = f.readline().split()
    f.close()
    # why is skiprows here?!
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




def generateBackwardMigMatrix(M, PastNe, RecentNe):  
#this function 'generateBackwardMigMatrix' should generate backward migration matrix from any M and any two past points. 
    n=np.size(M,0)
    print (n)
    BM=np.zeros([n,n])
    print (BM)
    for i in range (n):
        for j in range (n):
            BM[i,j]=M[j,i]*PastNe[j]/RecentNe[i]  # this may cause problems when RecentNe[i] is zero. 

    BM[np.isnan(BM)] = 0.0
    return BM

def extendMigMatrix(M):  #this is probably not necessary, only should be used at the end? It creates extended mig matrix
    M_ext=M
    for i in range(anc_num):
        M_ext=np.vstack([M_ext,np.zeros(n)])
    for i in range(anc_num):
        M_ext=np.hstack((M_ext,np.zeros((M_ext.shape[0],1),dtype=M_ext.dtype)))
    n_ext=M_ext.shape[0]
    return n_ext, M_ext






################# I think this compares two time points. If the past one is empty and the never one is not, it mass migrate
#if it is not empty, it does param change. What is BM used for?

def DefinePopParamForStep1(M, ADJ):  #this is time step between presence and the recent past
    if M[0,1]==0:    #if migration is 0... this is only used for mass migration
        pT=(M.T*tree_num[ngens-2])
    else:
        #pT=(M.T*tree_num[ngens-2]/M[0,1])
        pT=ADJ.T*tree_num[ngens-2]

    popT=pT.sum(axis=1)  #this generate a vector - for each cell it contains sum of all neighbours in the penultimate step. 
    
    print ('starting step 1')
    
    #popT[np.isnan(popT)] = 0.0
    print ('DefinePopParamForStep1 running')
    for i in range(n):
        if tree_num[ngens-2,i]==0 and tree_num[ngens-1,i]!=0:
            print ('colonisation encountered')  #this means that a cell previously empty is now colonised. 
            for j in range(n):
                if ADJ[i,j]!=0:  # ok this goes through all neighbouring cells, then makes popu. param change, which is identical to the one ourside of the colonisation check. 
                    demographic_events.append(msprime.PopulationParametersChange(
                        time=dt/genTime, initial_size=max(minSize, tree_num[ngens-2,j]),
                            population_id=j, growth_rate=0))

                    if popT[i]==0.0:  # and then in compares proportion 
                        prop=0  #this is where sad stuff happens. Means that a lineage gets stuck somewhere until the end of time. 
                        print ('sad stuff hapenning')
                    
                    else:
                        prop=tree_num[ngens-2,j]/popT[i]
#what would happen here if we exchanged the order of stuff hapenning? 
                    demographic_events.append(msprime.MassMigration(
                        time=dt/genTime, source=i, destination=j, proportion=min(1,prop)))

                    demographic_events.append(msprime.PopulationParametersChange(
                        time=dt/genTime, initial_size=max(minSize, tree_num[ngens-2,i]),
                            population_id=i, growth_rate=0))
                    #what is the difference bweteen this one and the one on line 217?! maybe this one should be inside?
        
        else:
            demographic_events.append(msprime.PopulationParametersChange(
                        time=dt/genTime, initial_size=max(minSize, tree_num[ngens-2,i]),
                            population_id=i, growth_rate=0))
    print ('DefinePopParamForStep1 ended')





def DefinePopParamForMiddleSteps(BM):
    t= ngens-2  
    x_last=tree_num[t]
    BM_last=BM
    x=tree_num[t]
    while t>=1:
        print(t)
        x_next=tree_num[t-1]
        t_ago = dt * (ngens - t)   #dt * (ngens -1- t)
        N=(M.T*x_next)  #why is this needed?
        mig=[]
        for i in range(n):
            mig.append([x[i]])
        migrate=M[0,1]   #why is this a special case??
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
                            time=t_ago/args.gen_time, initial_size=max(minSize, x_next[j]),
                            population_id=j, growth_rate=0))
                        if popT[i]==0.0:
                            prop=0
                        else:
                            prop=x_next[j]/popT[i]
                        demographic_events.append(msprime.MassMigration(
                            time=t_ago/args.gen_time, source=i, destination=j, proportion=min(1, prop)))
                        demographic_events.append(msprime.PopulationParametersChange(
                            time=t_ago/args.gen_time, initial_size=max(minSize, x_next[i]),
                            population_id=i, growth_rate=0))
            else:
                demographic_events.append(msprime.PopulationParametersChange(
                    time=t_ago/args.gen_time, initial_size=max(minSize, x_next[i]),
                    population_id=i, growth_rate=0))
        BM=N/mig  #why? 
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
        x=x_next.copy()  #this should be .copy() ???
        BM_last=BM
        t-=1
        print ('this is t', t)
        return x, BM_last

    
def populationMerge(BM):  #this is recalculating BM so it is ready for the last step. 
    t=0  
    N=(M.T*x)
    mig=[]
    migrate=M[0,1]

    for i in range(n):
        mig.append([x[i]])

    if migrate==0:
        pT=(M.T*x_next)
    else:
        pT=(M.T*x_next/migrate)
    popT=pT.sum(axis=1)
    BM=N/mig
    print ('new BM', BM)
    BM[np.isnan(BM)] = 0.0
    BM[BM == inf] = 0.0


    for i in range(n):
        for j in range(n):
            #print ('for i, j', i, j)
            #print (BM[i,j])
            #print (BM_last[i,j])
            if BM[i,j]!=BM_last[i,j]:
                
                demographic_events.append(msprime.MigrationRateChange(
                        time=dt * (ngens -1- t)/args.gen_time, rate=BM[i,j], matrix_index=(i, j)))  #for some reason this is not happening! 

    return (BM)


def MassMigrationToAncPop(BM):

    kk=0
    print ('kk was zeroed')
    t=0  #again this was added but maybe it does not have to be there?!


    for i in range(anc_num):
        demographic_events.append(msprime.PopulationParametersChange(
        time=dt * (ngens- t)/args.gen_time, initial_size=ancPopSizes[i],
        population_id=n+i, growth_rate=0))  #I think this puts initial pop sizes (scaled) to all new ancestral populations


    for i in range(n):   #and this mass migrates all the other lineages into these ancestral populations
        demographic_events.append(msprime.MassMigration(
            time=dt * (ngens - t)/args.gen_time, source=i, destination=n+int(ancpop_list[[i]])-1, proportion=1.0))


    ## migration rates are 0s for the grid, and for the ancestral population it's something small
    print ('this is BM', BM)  #Here BM is different from Eniko's
    for i in range(n):
        for j in range(n):
            #print (BM[i,j])
            if BM[i,j]!=0.0:
                demographic_events.append(msprime.MigrationRateChange(
                    time=dt * (ngens - t)/args.gen_time, rate=0.0, matrix_index=(i, j)))

                
                kk+=1
    print ('third',kk)
    for i in range(anc_num):
        for j in range(anc_num):
            if i!=j:
                demographic_events.append(msprime.MigrationRateChange(
                        time=dt * (ngens - t)/args.gen_time, rate=1e-8, matrix_index=(n+i, n+j)))
                print ('fourth',kk)
                kk+=1
  
    
####################################### STARING TO READ IN DATA ################################################
#######################################################################################

# read the migration rates and the population sizes
tree_num, celllist = read_inputs()
M = read_migration_matrix("M{}.tsv".format(m))
ADJ = read_migration_matrix("{}.tsv".format(am))
anc_num=len(args.ancpop_size)   #what is this for?
ancpop_list=np.loadtxt(args.ancpop_list)   #this is a list assigning ancestral pop to each cell of the grid.  
ngens = tree_num.shape[0]
ancPopSizes=args.ancpop_size.copy()
n = M.shape[0]   #this gives the dimension of migration matrix - basically the whole number of cells, LxL = why? to check?


file1 = "NEW_pairwise.mig{}.ancpopnum{}.dem_{}.dt{}.serial{}.tsv".format(m,anc_num,ps,dt,args.serial)  


print ('My Adj matrix:', ADJ)
print ('migration matrix:', M)
print ('treee numbers are:', tree_num)
print ('number of known steps:',ngens)
print ('ancestral populations:', ancpop_list)
print ('output file name is:', file1)
print ('ancestral population sizes:', ancPopSizes)
print ('total number of cells is:', n)

######I think we should remove this######################################################

#def ancPopSize(scale=1):
#    a_s=[]
#    for x in args.ancpop_size:
#        a_s.append(int(x))
#    print (a_s)   
#    print (a_s[0])            
#    return a_s



#exit()
########################################################################






# the extended migration includes the ancestral populations, but the migration rate
# to and from these cells are 0. This is only used for the simulate() function
# in every other case the original migration matrix is used together with the population sizes
# without the ancestral ones.




###################### POPULATION CONFIGURATIONS, SAMPLES ###############################
#########################################################################################


[n_ext, M_ext] = extendMigMatrix(M)
print ('extenden matrix size', n_ext)  #this is all original cells plus all non-spatially explicit ones ancestral



# to include the extra ancestral populations, we need to add extra 0s to the vector
# ndip and and nsamples, however these will not affect the iteration or any other calculations

#celllist = read_coords()  


########################No idea what this is and why we need it###########################

#########################################################################################

#nsamples = [0 for _ in range(n_ext)]  #this is actuallt never used?!!!
ndip = [0 for _ in range(n_ext)]

#print (nsamples)
print (ndip)  # - just plenty of zeros

for k in range(n):
#	nsamples[k] = 2*sum([u == k for u in celllist])
	ndip[k]=sum([u == k for u in celllist])
#number_of_samples=sum([nsamples[k] for k in range(n_ext)])


#print (nsamples)
print (ndip)  # and here it is the same 

#########################################################################################

########################No idea what this is and why we need it###########################

print ('this is cellist', celllist)  #this is coordinates (cell number) of each sample

sampled_demes=list(set(celllist))
print ('this is sampled demes', sampled_demes)  #this is which cells got sampled
sampled_demes.sort()
print ('this is sorter sampled demes', sampled_demes)


samp=list(set(tree_num[-1][sampled_demes]))  #why is list and set here?   IT IS NOT AS BELOW. WHY. 
print ('unsorted samp', samp)
samp2=(tree_num[-1][sampled_demes])
print ('samples 2',samp2)


print ('test cell', tree_num[-1][4])
samp.sort()
print ('sorted samp', samp)

sd=len(sampled_demes)   #how many demes got sampled

#########################################################################################



# all this below could be in the function that reads in the file


checkParameters()



#########No idea what this is for###############

subpops=[]
counter=0
for i in range(len(ndip)):  #this generates list of subpopulations, with indexed samples (0-59)
    if ndip[i]!=0:
        elem=[]
        for i in range(counter, counter+ndip[i]):
            elem.append(i)
        subpops.append(elem)
        counter=elem[-1]+1
        print ('counter',counter)
        print ('subpop',subpops)




tree_origi=[0 for _ in range(n_ext)]
tree_origi_last=[0 for _ in range(n_ext)]

print (tree_origi)



for k in range(n):
	tree_origi[k] = tree_num[ngens-1,k]
	tree_origi_last[k] = tree_num[ngens-2,k]

    #how is all this stuff dfferent from just taking the corresponing line?!


print (tree_origi)
print (tree_origi_last)




#########################################################################################
#########################################################################################



######### DEFINE DEMOGRAPHIC EVENTS THROUGH THE POPULATION HISTORY ######################
#########################################################################################

population_configurations=[]
demographic_events = []



##################################### FIRST TIMESTEP ####################################
#########################################################################################


# *_ext are only to account for the ancestral populations, have no role until the final timestep
# N is a matrix: how many individuals in pop. i have parents in j
# if a deme has individuals now but in the previous timestep was empty, we need a mass migration (backward in time):
# mass migration from this cell to the neighbors, proportions are weighted by the sizes of the population sizes in the previous generation
# pT, popT is to decide about these ratios when mass migration is needed
# when mass migration happens, a population goes empty, so the order of events (may?) matter
# if there is mass migration from i to j, we first set the new size of j, mass migrate the lineages, set the size of i to minSize
# migration matrix (i,j): how many individuals come from cell j to i, divided by the size of deme i


    

backwardMigM=generateBackwardMigMatrix(M, tree_num[ngens-2], tree_num[ngens-1])
print ('This is BM', backwardMigM)

N=(M.T*tree_num[ngens-2])  #M.T is transposed M... 
mig=[]
N_ext=(M_ext.T*tree_origi_last)
mig_ext=[]

print ('this is N', N)
print ('this is M.t', M.T)
print ('this is M', M)
print (type(M))


for i in range(n):
    mig.append([tree_num[ngens-1,i]])

    




DefinePopParamForStep1(M, ADJ)



            
BM=N/mig   #why this   #not yet used for anything, it will be used to    
BM[np.isnan(BM)] = 0.0
from numpy import inf
BM[BM == inf] = 0.0


#same for _ext

for i in range(n_ext):
    mig_ext.append([tree_origi[i]])



            
#why is there no migration change such as as in the next steps?!

################################### ITERATION ###########################################
#########################################################################################

#why is this different?


# we need to go through all the timesteps, the first (above) and the last are treated separately,
# between there is a while loop:
# x is the tree numbers for all the cells now, x_next is the previous generation
# migration matrix is made the same way as before
# mass migrations and population size updates work the same way
# population size may get updated twice...





        
[x, BM_last]= DefinePopParamForMiddleSteps(BM)
x_next=x.copy()
print ('BM', BM)
print ('BM last', BM_last)

################################# LAST TIMESTEP #########################################
#########################################################################################

print ('this is BM before merge', BM)
BM=populationMerge(BM)
print ('this is BM after merge', BM)


BM_ext=N_ext/mig_ext
BM_ext[np.isnan(BM_ext)] = 0.0

from numpy import inf
BM_ext[BM_ext == inf] = 0.0

bwM_ext=BM_ext



###################### MASS MIGRATION TO ANCESTRAL POPULATIONS ##########################
#########################################################################################


# first grow the population sizes of the ancestral ones to the desired number
# find the closest ancestral population for each lineage and then mass-migrate the lineage there

              
MassMigrationToAncPop(BM)



##################################### Initial pop config ########################################           

population_configurations = [
    msprime.PopulationConfiguration(sample_size=ndip[k],
        initial_size = max(minSize, tree_origi[k]),
        growth_rate=0) for k in range(n_ext) ]

#makes a list of length n_ext and for every cell size it sets sample size and cell size, only time when you do pop configuration This is a starting point for the msprime simulations. 
 
            
##################################### SIMULATION ########################################
#########################################################################################

#this is the starting migration matrix, if in the previous generation a cell was empty, the value should be 0
BM_ext=N_ext/mig_ext
BM_ext[np.isnan(BM_ext)] = 0.0

from numpy import inf
BM_ext[BM_ext == inf] = 0.0

bwM_ext=BM_ext

for i in range(n):
    for j in range(n):
        if tree_num[ngens-2,j]==0:
            bwM_ext[i,j]=0


BM_ext=N_ext/mig_ext   #only needed for starting the simulation
BM_ext[np.isnan(BM_ext)] = 0.0


print("sim start")
dd = msprime.DemographyDebugger(
    population_configurations=population_configurations,
    migration_matrix=bwM_ext,#
    demographic_events=demographic_events)

log_file = open("DebuggerNew.log","w")
dd.print_history(output=log_file)


replicates = msprime.simulate(
	population_configurations=population_configurations,
	demographic_events=demographic_events,
	migration_matrix=bwM_ext,  #this is weird why is it here if it is not necessary, check
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
print ('all finished, huraaay!')