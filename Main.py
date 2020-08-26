##########################################  FUNCTIONS #################################
#######################################################################################

def parseInputFiles(parser):
    parser.add_argument("--pop_sizes", "-pop", type=str, dest="tree_nums_file", required=True,
                    help="name of file to load tree number history from")
    #this loads the file with the demographic history, all cells population sizes in time

    parser.add_argument("--migr_mat", "-mig", type=str, dest="migr_mat_file", required=True,
                    help="name of file containing migration matrix")  
    #this loads migration matrix
    #make this not obligatory a file - if single number is supplied, make symetric matrix. 
    
    #parser.add_argument("--adj_mat", "-am", type=str, dest="adj_mat_file", required=True,
    #                help="name of file containing adjacency matrix")
    #this loads adjacency matrix for when we need migration = 0

    parser.add_argument("--sample_coords", "-sam", type=str, dest="sample_coords_file", required=True,
                    help="name of file containing locations of sampling coordinates")
    #this loads sample coordinates -  i am not sure what is this for. 
    
    parser.add_argument("--row_number", "-row", type=int, dest="row_num", required=True,
                    help="number of rows in the originical grid")
    #number of rows in the original grid. Just for the first check
    
    #############################################  OPTIONAL stuff ########################################## 
    ######################################################################################################

    parser.add_argument("--generation_time", "-gen", type=int, dest="gen_time", default = 1,
                    help="generation time, default 1")
    #generation time 
    
    parser.add_argument("--ancpop_list", "-apl", type=str, dest="ancpop_list", default = 'nan',
                    help="name of file to load a list with the destination ancestral population for each gridcell")
    #this defines ancestral population for each cell. 
    parser.add_argument("--ancpop_size", "-aps",  type=int, nargs='+', dest="ancpop_size", default = [1], 
                    help="list of the sizes of the ancestral populations, default 1 population of size 1")
    #this should give sizes of ancestral populations, should be of dimension of the number of anc pop  %this is wrong, takes integer

    
    parser.add_argument("--dt", "-dt", type=float, dest="dt", default = 10,
                        help="time interval between steps, default 10")
    #This gives dt- time between two known points

     #############################################  SYS stuff.  check if this is needed ########################################## 
    ######################################################################################################
   
    parser.add_argument("--serial", "-ser", type=int, dest="serial",
                    help="serial number")
    #I think this only assignes simulation serial number so it is saved into a different output file. 

    parser.add_argument("--output_dir", "-od", type=str, dest="basedir",
                    help="name of directory to save output files to.")
    #this is important!!!  NOT USED AT THE MOMENT?

    parser.add_argument("--output_logfile", "-of", type=float, dest="logfile",
                    help="name of log file")
    #define log file
    
    parser.add_argument("--prit_debugger", "-pdeb", type=bool, dest="print_deb", default=False, 
                    help="printing debugger option")
    #define log file
 
    args = parser.parse_args()
    m = args.migr_mat_file
    #am = args.adj_mat_file
    ps = args.tree_nums_file
    serial = args.serial
    sample_list=args.sample_coords_file
    dt = args.dt  #timestep
    genTime=args.gen_time
    ancPopSizes=args.ancpop_size
    ancpop_list=args.ancpop_list
    print_deb=args.print_deb
 
    return args, m, ps, serial,  dt, genTime, sample_list, ancPopSizes, ancpop_list, print_deb




def read_input_data(args, m, ps, serial,  dt, genTime,  sample_list, ancPopSizes, ancpop_list):

    tree_num = np.loadtxt("{}.tsv".format(ps))
    sampled_demes = np.loadtxt("{}.tsv".format(sample_list), dtype='int')
    migList = np.loadtxt("{}.tsv".format(m), dtype='float')
    #ADJ = read_migration_matrix("{}.tsv".format(am))

    anc_num=len(ancPopSizes)   #what is this for?
       
    ngens = tree_num.shape[0]
    n = tree_num.shape[1]
    n_ext=n+anc_num
    
    if ancpop_list == 'nan':
        ancpop_list=np.zeros(n)+1
        
        ancpop_list=ancpop_list.astype(int)
    else:
        ancpop_list=np.loadtxt("{}.tsv".format(ancpop_list), dtype='int')   #this is a list assigning ancestral pop to each cell of the grid.
    #values, counts = np.unique(words, return_counts=True)
    sampled_demes.sort()
    samp=list(set(tree_num[-1][sampled_demes]))  #all this fun is just to ceck if we don't sample an empty cell. Why 'set'?
    
    minSize = 1e-10 #
    
    assert tree_num.shape[1]%args.row_num == 0, "Not a rectangular map. Number of cells is not divisible by the number of rows!" 
    assert samp[0]>0, "You cannot sample from cells that are empty!"
    #print ('All is fine')
    ndip = np.zeros(n_ext)
    
    ndip[sampled_demes]=2
    ndip=ndip.astype(int)
    
    return tree_num, migList, anc_num, ancpop_list, ngens, n, n_ext, minSize, sampled_demes, ndip



def generateExtendedBMMatrix(n_ext, bMigList):
    BMM=np.zeros([n_ext, n_ext])
    for i in range(len(bMigList[:,0])):
        sourceCell=int(bMigList[i,0])  #source cell when thinking forward in time
        targetCell=int(bMigList[i,1])
        bmigRate=(bMigList[i,2])
        BMM[sourceCell,targetCell]= bmigRate
    return BMM

def generateExtPopN(tree_num, n, n_ext, ngens):
    presentPop=np.zeros(n_ext)
    for i in range(n):
        presentPop[i]=tree_num[ngens-1,i]
    return presentPop


def generateBackwardMigList(migList, RecentNe, PastNe):  
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




def MassMigrationToAncPop(demographic_events, dt, anc_num, ngens, genTime, ancPopSizes,ancpop_list, n, minSize):

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
                        time=dt * (ngens )/genTime, rate=1e-8, matrix_index=(n+i, n+j)))
    




##################################### FIRST TIMESTEP ####################################
# if a deme has individuals now but in the previous timestep was empty, we need a mass migration (backward in time):
# mass migration from this cell to the neighbors, proportions are weighted by the sizes of the population sizes in the previous generation
# pT, popT is to decide about these ratios when mass migration is needed
# when mass migration happens, a population goes empty, so the order of events (may?) matter
# if there is mass migration from i to j, we first set the new size of j, mass migrate the lineages, set the size of i to minSize
# migration matrix (i,j): how many individuals come from cell j to i, divided by the size of deme i

def DefinePopParamInFirstStep(demographic_events, migList, tree_num, ngens, dt, genTime, minSize):  #this is time step between presence and the recent past
    n=len(tree_num[0])
    popT=np.zeros(np.shape(tree_num[ngens-1]))
    if migList[0,2]!=0:  #this generate a vector - for each cell it contains sum of all neighbours in the penultimate step. 
        for i in range(len(migList[:,0])):
            cellOrigin=int(migList[i,0])
            cellTarget=int(migList[i,1])
            popT[cellOrigin]=popT[cellOrigin]+tree_num[ngens-2,cellTarget]
        #print (popT)
    
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
                        print ('Some lineages may be stuck, as cell no ', cellColonised, 'has no neighbours in the past that could colonize it. ')
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







def DefinePopParamForMiddleSteps2(demographic_events, migList, tree_num, ngens, dt, genTime, minSize):

    t= ngens-2  
    n=len(tree_num[0])
    for t in range(ngens-2, 0, -1):
        #print('this is time', t)
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
            
            #print ('in gen ',dt*(ngens - t-1)/genTime,'modified for cell', sourceCell, 'to cell',  targetCell,' rate ',bmigRate)
        
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

             


    return bMigList



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
            
            
def populationMerge(demographic_events, migList, initialN, genTime, dt,  ngens):  #this is recalculating BM so it is ready for the last step. 
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
            
            #print ('in gen ',t_ago/genTime,'modified for cell', sourceCell, 'to cell',  targetCell,' rate ',bmigRate)
 
    
####################################### STARING TO READ IN DATA ################################################
#######################################################################################

def extract_results(sim_results, subpops, file1):
    tree = sim_results.first()
    aa=sim_results.simplify().first()
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
    

def run_gridcoal():

    startT = time.time()
    print("sim start")
    seed_no = np.random.get_state()
   #!/usr/bin/env python3
    description = '''Simulate a coalescent within the history of Abies expansion.'''

    parser = argparse.ArgumentParser(description = description)
    [args, m, ps, serial,  dt, genTime, sample_list, ancPopSizes, ancpop_list, print_deb] = parseInputFiles(parser)

    [tree_num,  migList,  anc_num, ancpop_list, ngens,  n, n_ext, minSize,sampled_demes, ndip]=read_input_data(args, m, ps, serial, dt, genTime,  sample_list, ancPopSizes, ancpop_list)
    
    


    file1 = "CoalTime_M_{}.AncPop_{}.Dem_{}.dt{}.serial{}.tsv".format(m,anc_num,ps,dt,serial)  
    print (serial)
    if serial ==1: 
        orig_stdout = sys.stdout
        f = open('OutputFile.txt','w')
        sys.stdout = f
        print ('INPUT FILES')
        print('using input file for population sized:')
        print (ps)
        print ('input sizes are:', tree_num)
        print ('total number of cells is:', n)
        print ('number of known steps:',ngens)
        print ('time between steps:', dt)
        print('min population size is', minSize)
        print ('total known time is:', dt*ngens)
        print ('')

        print ('SAMPLING')
        print('using input file for sample sizes: ', ps)
        print ('sample sizes taken:', ndip)
        print ('')

        print ('MIGRATION')
        print('using input file for migration', m)
        print('migration is ["from" "to" "rate"]:')
        print (migList)
        print ('')
        print ('ANCESTRAL POPULATION')

        print ('ancestral population IDs:')
        print (ancpop_list)
        print ('ancestral population sizes:', ancPopSizes)
        print ('extenden population number:', n_ext) 
        print ('')


        print ('output file name is:', file1)
        
        
        print ('')
        
    
        if print_deb == False:
            
            sys.stdout = orig_stdout
            f.close()

    ######### DEFINE DEMOGRAPHIC EVENTS THROUGH THE POPULATION HISTORY ######################
    #########################################################################################

    population_configurations=[]
    demographic_events = []

    DefinePopParamInFirstStep(demographic_events, migList.copy(), tree_num.copy(), ngens, dt,  genTime, minSize)       
    bMigList = DefinePopParamForMiddleSteps2(demographic_events, migList.copy(), tree_num.copy(), ngens, dt, genTime, minSize)
    populationMerge(demographic_events, migList.copy(), tree_num[0], genTime, dt,  ngens)

    MassMigrationToAncPop(demographic_events, dt, anc_num, ngens, genTime, ancPopSizes, ancpop_list, n, minSize)


    ################################# Initial pop config ########################################           
   

    presentExtendedN=generateExtPopN(tree_num, n, n_ext, ngens)
      
    bMigList=generateBackwardMigList(migList.copy(), tree_num[ngens-1,:], tree_num[ngens-2,:])
    BMM=generateExtendedBMMatrix(n_ext, bMigList)  

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=ndip[k],
            initial_size = max(minSize, presentExtendedN[k]),
            growth_rate=0) for k in range(n_ext) ]

        
    ##################################### SIMULATION ########################################
    #########################################################################################

    
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=BMM,#
        demographic_events=demographic_events)
    if serial ==1: 
    #    orig_stdout = sys.stdout
    #    f = open('OutputFile.txt','w')
    #    sys.stdout = f
        if print_deb == True:
            print ('DEMOGRAPHY DEBUGGER')
            dd.print_history(output=f)
            print ('')
            sys.stdout = orig_stdout
            f.close()


    sim_results = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=BMM)  



    ########################CALCULATE EXPECTED PAIRWISE COALESCENCE TIMES####################
    #########################################################################################
    subpops = defineSubpopulations(ndip)  #this is necessary for results analysis
    extract_results(sim_results, subpops, file1)
    
    endT = time.time()
    print ('Simulation ended in time', endT - startT, 'seconds' )
 

    

###################################################################################################################
############################################# MAIN FUNCTION STARTING HERE #########################################
###################################################################################################################

import msprime
import numpy as np
import os
import argparse
import sys

import time

run_gridcoal()



