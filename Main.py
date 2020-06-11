#!/usr/bin/env python3
description = '''Simulate a coalescent within the history of Abies expansion.'''
import time
startT = time.time()

import msprime
import numpy as np
import os
import argparse
#import timeit

##########################################  FUNCTIONS #################################
#######################################################################################


def parseInputFiles(parser):
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
    return args, m, am, ps, serial, scale, dt, genTime, minSize




def read_inputs():
    # e.g. tree_num_epsg3035.tsv
    trees = np.loadtxt("{}.tsv".format(ps))
    celllist = np.loadtxt(args.sample_coords_file, dtype='int')
    migList = np.loadtxt("M1Simple.tsv", dtype='float')
    
    return trees, celllist, migList

def checkParameters(celllist, tree_num):
    sampled_demes=list(set(celllist))
    sampled_demes.sort()
    samp=list(set(tree_num[-1][sampled_demes]))  #all this fun is just to ceck if we don't sample an empty cell. Why 'set'?
    
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



def generateExtendedBMMatrix(n_ext, bMigList):
    BMM=np.zeros([n_ext, n_ext])
    for i in range(len(bMigList[:,0])):
        sourceCell=int(bMigList[i,0])  #source cell when thinking forward in time
        targetCell=int(bMigList[i,1])
        bmigRate=(bMigList[i,2])
        BMM[sourceCell,targetCell]= bmigRate
    print (BMM)
    return BMM

def generateExtPopN(tree_num, anc_num, n_ext, ngens):
    presentPop=np.zeros(n_ext)
    for i in range(n):
        presentPop[i]=tree_num[ngens-1,i]
    return presentPop

def generateBackwardMigMatrix(M, PastNe, RecentNe):  #I think this is not needed anymore
#this function 'generateBackwardMigMatrix' should generate backward migration matrix from any M and any two past points. 
    n = np.size(M,0)
    print (n)
    BM = np.zeros([n,n])
    print (BM)
    for i in range (n):
        for j in range (n):
            BM[i,j]=M[j,i]*PastNe[j]/RecentNe[i]  # this may cause problems when RecentNe[i] is zero. 

    BM[np.isnan(BM)] = 0.0
    return BM



#I think I am not using this
def extendMigMatrix(M):  #this is probably not necessary, only should be used at the end? It creates extended mig matrix
    M_ext=M
    for i in range(anc_num):
        M_ext=np.vstack([M_ext,np.zeros(n)])
    for i in range(anc_num):
        M_ext=np.hstack((M_ext,np.zeros((M_ext.shape[0],1),dtype=M_ext.dtype)))
    n_ext=M_ext.shape[0]
    return n_ext, M_ext




def generateBackwardMigList(migList, RecentNe, PastNe):  
#this function 'generateBackwardMigMatrix' should generate backward migration matrix from any M and any two past points. 
    bMigList = []
    for i in range(len(migList[:,0])):
        sourceCell=int(migList[i,0])  #source cell when thinking forward in time
        targetCell=int(migList[i,1])  #target cell when thinking forward in time
        if RecentNe[targetCell]>0:
            bwMigRate=migList[i,2]*PastNe[sourceCell]/RecentNe[targetCell]  ### it take -1 because 
        else: 
            bwMigRate=0 #
        bMig=[targetCell,sourceCell,bwMigRate ]  #note the opposite direction
        bMigList.append(bMig)
    bMigList=np.array(bMigList)
    bMigList=bMigList[np.argsort(bMigList[:,0])]
    return bMigList




def MassMigrationToAncPop():

    # Creates ancestral populations that were zeros until now
    for i in range(anc_num):
        demographic_events.append(msprime.PopulationParametersChange(
        time=dt * (ngens)/genTime, initial_size=ancPopSizes[i],
        population_id=n+i, growth_rate=0))  #I think this puts initial pop sizes (scaled) to all new ancestral populations

    #### And this mass migrates all the other lineages into these ancestral populations
    for i in range(n):   
        demographic_events.append(msprime.MassMigration(
            time=dt * (ngens)/genTime, source=i, destination=n+int(ancpop_list[[i]])-1, proportion=1.0))
        demographic_events.append(msprime.PopulationParametersChange(
                            time=dt * (ngens)/genTime, initial_size=minSize,
                                population_id=i, growth_rate=0))

    #### Migration rates are 0s for the grid, and for the ancestral population it's something small
    for i in range(n):
        for j in range(n):
            if i != j: 
                demographic_events.append(msprime.MigrationRateChange(
                    time=dt * (ngens)/genTime, rate=0.0, matrix_index=(i, j)))

    for i in range(anc_num):
        for j in range(anc_num):
            if i != j:
                demographic_events.append(msprime.MigrationRateChange(
                        time=dt * (ngens )/args.gen_time, rate=1e-8, matrix_index=(n+i, n+j)))
    




##################################### FIRST TIMESTEP ####################################
# if a deme has individuals now but in the previous timestep was empty, we need a mass migration (backward in time):
# mass migration from this cell to the neighbors, proportions are weighted by the sizes of the population sizes in the previous generation
# pT, popT is to decide about these ratios when mass migration is needed
# when mass migration happens, a population goes empty, so the order of events (may?) matter
# if there is mass migration from i to j, we first set the new size of j, mass migrate the lineages, set the size of i to minSize
# migration matrix (i,j): how many individuals come from cell j to i, divided by the size of deme i

def DefinePopParamInFirstStep(migList):  #this is time step between presence and the recent past
    popT=np.zeros(np.shape(tree_num[ngens-1]))
    if migList[0,2]!=0:  #this generate a vector - for each cell it contains sum of all neighbours in the penultimate step. 
        for i in range(len(migList[:,0])):
            cellOrigin=int(migList[i,0])
            cellTarget=int(migList[i,1])
            popT[cellOrigin]=popT[cellOrigin]+tree_num[ngens-2,cellTarget]
        print (popT)
    
    for i in range(n):
        if tree_num[ngens-2,i]==0 and tree_num[ngens-1,i]!=0: #this means that a cell previously empty is now colonised. 
            cellColonised=i
            for j in range(len(migList[:,0])):
                if cellColonised==int(migList[j,0]):
                    cellSource=int(migList[j,1])
                    demographic_events.append(msprime.PopulationParametersChange(  
                        time=dt/genTime, initial_size=max(minSize, tree_num[ngens-2,cellSource]),
                            population_id=cellSource, growth_rate=0))
                     #this will increase source pop size, so it is filled before mass migrating the stuff into it
                    if popT[cellColonised]==0.0:  # this means that sum od neighbours is 0 - sad stuff
                        prop=0  #this is where sad stuff happens. Means that a lineage gets stuck somewhere until the end of time. 
                        print ('sad stuff hapenning: some lineages may get stuck, as cell no ', cellColonised, 'has no neighbours in the past that could colonize it. ')
                    else:
                        prop=tree_num[ngens-2,cellSource]/popT[cellColonised]

                    demographic_events.append(msprime.MassMigration(
                        time=dt/genTime, source=cellColonised, destination=cellSource, proportion=min(1,prop)))
 #maybe the problem here is that it is hapenning to early. They should be allowed to move freely first then to mass migrate. Maybe this could be corrected by changing the time for mass migration further back? But this seems fine, they should be able to move around for time dt/genTime... Why do they get stuck? 

                    demographic_events.append(msprime.PopulationParametersChange(
                        time=dt/genTime, initial_size=max(minSize, tree_num[ngens-2,cellColonised]),
                            population_id=cellColonised, growth_rate=0))
        else:  #some of these mau be already altered, but it should not matter. 
            demographic_events.append(msprime.PopulationParametersChange(
                        time=dt/genTime,
                        initial_size=max(minSize, tree_num[ngens-2,i]),
                        population_id=i,
                        growth_rate=0))







def DefinePopParamForMiddleSteps2(migList):

    t= ngens-2  
    for t in range(ngens-2, 0, -1):
        print('this is time', t)
        thisN=tree_num[t]
        prevN=tree_num[t-1]  #this is more in the 

       ######################## This modifies migrations #################
        bMigList=generateBackwardMigList(migList.copy(), thisN.copy(), prevN.copy())  
        #here add stuff for checking the lists! only updating when necessary. 
        
        for i in range(len(bMigList[:,0])):
            sourceCell=int(bMigList[i,0])  #source cell when thinking forward in time -  PROBABLY NOT?!!!
            targetCell=int(bMigList[i,1])
            bmigRate=(bMigList[i,2])
            demographic_events.append(msprime.MigrationRateChange(
                         time=dt*(ngens - t-1)/genTime, rate=bmigRate, matrix_index=(sourceCell, targetCell)))
            
            print ('in gen ',dt*(ngens - t-1)/genTime,'modified for cell', sourceCell, 'to cell',  targetCell,' rate ',bmigRate)
        
        t_ago = dt * (ngens - t)   #dt * (ngens -1- t)
        popT=np.zeros(np.shape(tree_num[ngens-1]))   #######maybe this whole chunk could be a special function?
        if migList[0,2]!=0:  # this checks if migration is zero, only if not, this generates a vector - for each cell it contains sum of all neighbours in the previous step. 
            for i in range(len(migList[:,0])):
                cellOrigin=int(migList[i,0])
                cellTarget=int(migList[i,1])
                popT[cellOrigin]=popT[cellOrigin]+prevN[cellTarget]  #CHECK THIS


##############
        for i in range(n):
            if prevN[i]==0 and thisN[i]!=0:#this means that a cell previously empty is now colonised. 
                cellColonised=i
                for j in range(len(migList[:,0])):
                    if cellColonised==int(migList[j,0]):
                        cellSource=int(migList[j,1])
                        demographic_events.append(msprime.PopulationParametersChange(  
                            time=t_ago/genTime, initial_size=max(minSize, prevN[cellSource]),
                                population_id=cellSource, growth_rate=0))
                         #this will increase source pop size, so it is filled before mass migrating the stuff into it
                        if popT[cellColonised]==0.0:  # this means that sum od neighbours is 0 - sad stuff
                            prop=0  #this is where sad stuff happens. Means that a lineage gets stuck somewhere until the end of time.                
                        else:
                            prop=prevN[cellSource]/popT[cellColonised]

                        demographic_events.append(msprime.MassMigration(
                            time=t_ago/genTime, source=cellColonised, destination=cellSource, proportion=min(1,prop)))
                        demographic_events.append(msprime.PopulationParametersChange(
                            time=t_ago/genTime, initial_size=max(minSize, prevN[cellColonised]),
                                population_id=cellColonised, growth_rate=0))
            else:  #some of these mau be already altered, but it should not matter. 
                demographic_events.append(msprime.PopulationParametersChange(
                            time=t_ago/genTime, initial_size=max(minSize, prevN[i]),
                                population_id=i, growth_rate=0))

             
    #exit()
    
    return bMigList

def generateSampleList(celllist, n_ext, n):
    ndip = np.zeros(n_ext)
    for i in range(n):
        ndip[i]=int(sum([u == i for u in celllist]))
    ndip=ndip.astype(int)
    return ndip

def defineSubpopulations(ndip):  #these subpopulations are needed for final results, getting T between and T within. 
    subpops=[]
    counter=0
    for i in range(len(ndip)):  #this generates list of subpopulations, with indexed samples (0-59)
        if ndip[i]!=0:
            elem=[]
            for i in range(counter, counter+ndip[i]):
                elem.append(i)
            subpops.append(elem)
            counter=elem[-1]+1
    return subpops 
            
            
def populationMerge( migList, initialN, genTime):  #this is recalculating BM so it is ready for the last step. 
    #this cheats with the backward migration - assumes that the population is fixed, takes the same time point (the last one) and calculates bwmig between those. 
    bMigList=generateBackwardMigList(migList.copy(), initialN, initialN)  #twice the same thing
    t_ago=dt*(ngens-1)
    #print ('new way in generation', t_ago/genTime)
    
    for i in range(len(bMigList[:,0])):
            sourceCell=int(bMigList[i,0])  #source cell when thinking forward in time
            targetCell=int(bMigList[i,1])
            bmigRate=(bMigList[i,2])
            demographic_events.append(msprime.MigrationRateChange(
                         time=t_ago/genTime, rate=bmigRate, matrix_index=(sourceCell, targetCell)))
            
            print ('in gen ',t_ago/genTime,'modified for cell', sourceCell, 'to cell',  targetCell,' rate ',bmigRate)
 
    
####################################### STARING TO READ IN DATA ################################################
#######################################################################################

# read the migration rates and the population sizes



parser = argparse.ArgumentParser(description = description)
[args, m, am, ps, serial, scale, dt, genTime, minSize] = parseInputFiles(parser)



print('I have: mig file', m)
print('Adj matrix', am)
print('TreeInputfile', ps)
print('serial', serial)
print('scaling', scale)
print('dt is',  dt)
print('min size is', minSize)



tree_num, celllist, migList = read_inputs()
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





# the extended migration includes the ancestral populations, but the migration rate
# to and from these cells are 0. This is only used for the simulate() function
# in every other case the original migration matrix is used together with the population sizes
# without the ancestral ones.




###################### POPULATION CONFIGURATIONS, SAMPLES ###############################
#########################################################################################


[n_ext, M_ext] = extendMigMatrix(M)
print ('extenden matrix size', n_ext)  #this is all original cells plus all non-spatially explicit ones ancestral




########################No idea what this is and why we need it###########################




ndip=generateSampleList(celllist, n_ext, n)  #### somehow this returns wrong type, does not work. 


#exit()

#########################################################################################

########################No idea what this is and why we need it###########################

print ('this is cellist', celllist)  #this is coordinates (cell number) of each sample



#exit()
######### DEFINE DEMOGRAPHIC EVENTS THROUGH THE POPULATION HISTORY ######################
#########################################################################################
checkParameters(celllist, tree_num)


population_configurations=[]
demographic_events = []

DefinePopParamInFirstStep(migList.copy())       
bMigList = DefinePopParamForMiddleSteps2( migList.copy())
populationMerge(migList.copy(), tree_num[0], genTime)
MassMigrationToAncPop()



##################################### Initial pop config ########################################           
   

presentExtendedN=generateExtPopN(tree_num, anc_num, n_ext, ngens)
      
bMigList=generateBackwardMigList(migList.copy(), tree_num[ngens-1,:], tree_num[ngens-2,:])
BMM=generateExtendedBMMatrix(n_ext, bMigList)  

population_configurations = [
    msprime.PopulationConfiguration(sample_size=ndip[k],
        initial_size = max(minSize, presentExtendedN[k]),
        growth_rate=0) for k in range(n_ext) ]
#makes a list of length n_ext and for every cell size it sets sample size and cell size, only time when you do pop configuration This is a starting point for the msprime simulations. 



        
##################################### SIMULATION ########################################
#########################################################################################




print("sim start")
dd = msprime.DemographyDebugger(
    population_configurations=population_configurations,
    migration_matrix=BMM,#
    demographic_events=demographic_events)

log_file = open("DebuggerNew.log","w")
dd.print_history(output=log_file)


replicates = msprime.simulate(
	population_configurations=population_configurations,
	demographic_events=demographic_events,
	migration_matrix=BMM,  #this is weird why is it here if it is not necessary, check
	recombination_rate=args.recomb,
	mutation_rate=args.mut)


########################CALCULATE EXPECTED PAIRWISE COALESCENCE TIMES####################
#########################################################################################
subpops = defineSubpopulations(ndip)  #this is necessary for results analysis

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
endT = time.time()
print ('All is finished, hurray,  in time ', endT - startT)



exit()

################################ THIS IS TO BE THROWN AWAY ######################################
#################################################################################################
#################################################################################################

def DefinePopParamForMiddleSteps(BM, migList):

    t= ngens-2  
    
 
    x_last=tree_num[t]
    BM_last=BM
    x=tree_num[t]
    while t>=1:
        print('this is time', t)
        thisN=tree_num[t]
        prevN=tree_num[t-1]  #this is more in the 
       
        x_next=tree_num[t-1]


        t_ago = dt * (ngens - t)   #dt * (ngens -1- t)
        N=(M.T*x_next)  #why is this needed?
        mig=[]
        for i in range(n):
            mig.append([x[i]])
        #print (mig)  #no idea what this is for

        popT=np.zeros(np.shape(tree_num[ngens-1]))   #######maybe this whole chunk could be a special function?
        if migList[0,2]!=0:  # this checks if migration is zero, only if not, this generates a vector - for each cell it contains sum of all neighbours in the previous step. 
            for i in range(len(migList[:,0])):
                cellOrigin=int(migList[i,0])
                cellTarget=int(migList[i,1])
                popT[cellOrigin]=popT[cellOrigin]+prevN[cellTarget]  #CHECK THIS


##############
        for i in range(n):
            if prevN[i]==0 and thisN[i]!=0:#this means that a cell previously empty is now colonised. 
                cellColonised=i
                for j in range(len(migList[:,0])):
                    if cellColonised==int(migList[j,0]):
                        cellSource=int(migList[j,1])
                        demographic_events.append(msprime.PopulationParametersChange(  
                            time=t_ago/genTime, initial_size=max(minSize, prevN[cellSource]),
                                population_id=cellSource, growth_rate=0))
                         #this will increase source pop size, so it is filled before mass migrating the stuff into it
                        if popT[cellColonised]==0.0:  # this means that sum od neighbours is 0 - sad stuff
                            prop=0  #this is where sad stuff happens. Means that a lineage gets stuck somewhere until the end of time.                
                        else:
                            prop=prevN[cellSource]/popT[cellColonised]

                        demographic_events.append(msprime.MassMigration(
                            time=t_ago/genTime, source=cellColonised, destination=cellSource, proportion=min(1,prop)))
                        demographic_events.append(msprime.PopulationParametersChange(
                            time=t_ago/genTime, initial_size=max(minSize, prevN[cellColonised]),
                                population_id=cellColonised, growth_rate=0))
            else:  #some of these mau be already altered, but it should not matter. 
                demographic_events.append(msprime.PopulationParametersChange(
                            time=t_ago/genTime, initial_size=max(minSize, prevN[i]),
                                population_id=i, growth_rate=0))

             
        BM=N/mig  #why? 
        BM[np.isnan(BM)] = 0.0
        BM[BM == inf] = 0.0
        x=x_next.copy()  #this should be .copy() ???
        BM_last=BM
        print ('this N is', thisN)
        print ('prev N is', prevN)
        ######################## This modifies migrations #################
        bMigList=generateBackwardMigList(migList.copy(), thisN.copy(), prevN.copy())  
        #here add stuff for checking the lists! only updating when necessary. 
 
        for i in range(len(bMigList[:,0])):
            sourceCell=int(bMigList[i,0])  #source cell when thinking forward in time - PROBABLY NOT?!
            targetCell=int(bMigList[i,1])
            bmigRate=(bMigList[i,2])
            demographic_events.append(msprime.MigrationRateChange(
                         time=t_ago/genTime, rate=bmigRate, matrix_index=(sourceCell, targetCell)))
            
            print ('in gen ',t_ago/genTime,'modified for cell', sourceCell, 'to cell',  targetCell,' rate ',bmigRate)
        t-=1
    exit()
    
    return x, BM_last, bMigList
