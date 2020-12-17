"""Main programm running gridCoal."""

import argparse
import sys
import os
import time
import numpy as np
import msprime

##########################################  FUNCTIONS ####################
##########################################################################

def parse_input_files(parser):
    """Parse input arguments and return parameters or input files."""
    parser.add_argument(
        "--pop_sizes",
        "-pop",
        type=str,
        dest="tree_nums_file",
        required=True,
        help="name of file to load tree number history from")
    # this loads the file with the demographic history, all cells population
    # sizes in time

    parser.add_argument(
        "--row_number",
        "-row",
        type=int,
        dest="row_num",
        required=True,
        help="number of rows in the originical grid")
    # number of rows in the original grid. Just for the first check

####################################  OPTIONAL parameters ###############
#########################################################################
    parser.add_argument(
        "--migration_matrix",
        "-mig",
        type=str,
        dest="migr_mat_file",
        required=True,
        help="name of file containing migration matrix")
    # this loads migration matrix


    parser.add_argument(
        "--output_dir",
        "-odir",
        type=str,
        dest="basedir",
        default='OUTPUT',
        help="name of directory to save output files to.")

    parser.add_argument(
        "--sample_coords",
        "-sam",
        type=str,
        dest="sample_coords_file",
        default='nan',
        help="name of file containing locations of sampling coordinates")
    # this loads sample coordinates -  i am not sure what is this for.


    parser.add_argument(
        "--ancpop_list",
        "-apl",
        type=str,
        dest="ancpop_list",
        default='nan',
        help="name of file with a list with the ancestral population for each gridcell")
    # this defines ancestral population for each cell.

    parser.add_argument(
        "--ancpop_size",
        "-aps",
        type=str,
        #nargs='+',
        dest="ancpop_size",
        default='nan',
        help="list of the sizes of the ancestral populations, default 1 population of size 1")
    # this should give sizes of ancestral populations, should be of dimension
    # of the number of anc pop  %this is wrong, takes integer

    parser.add_argument(
        "--generation_time",
        "-gen",
        type=int,
        dest="gen_time",
        default=1,
        help="generation time, default 1")
    # generation time

    parser.add_argument("--delta_t", "-dt", type=float, dest="delta_t", default=10,
                        help="time interval between steps, default 10")
    # This gives delta_t- time between two known points

    parser.add_argument("--replicate", "-rep", type=int, dest="serial", default=1,
                        help="serial number")
    # Assignes simulation serial number so it is saved into a different output file.

    parser.add_argument(
        "--print_debugger",
        "-pdeb",
        type=bool,
        dest="print_deb",
        default=False,
        help="printing debugger option")
    # define log file

    parser.add_argument(
        "--set_seed",
        "-seed",
        type=int,
        dest="set_seed",
        default=0,
        help="setting seed number")
    # define log file

    args = parser.parse_args()
    mig_file = args.migr_mat_file
    demography_file = args.tree_nums_file
    serial = args.serial
    sample_list = args.sample_coords_file
    delta_t = args.delta_t  # timestep
    gen_time = args.gen_time
    anc_pop_sizes = args.ancpop_size
    ancpop_list = args.ancpop_list
    print_deb = args.print_deb
    output_dir=args.basedir

    seed_no=args.set_seed
    if seed_no == 0:
        seed_no = np.random.randint(1000000)

    return args, mig_file, demography_file, serial, delta_t, gen_time, sample_list, anc_pop_sizes, ancpop_list, print_deb, output_dir, seed_no

def sample_all_cells(tree_num):
    today=tree_num[-1, :]
    samplelist=[]
    for i in range(len(today)):
        if today[i]>0:
            samplelist.append(i)
    samplelist=np.array(samplelist)
    return samplelist

def generate_migration_list(rows, cols, m):
    n = rows*cols
    add_mat = np.zeros([n,n])
    for r in range(rows):
        for c in range(cols):
            i = r*cols + c
            # Two inner diagonals
            if c > 0: add_mat[i-1,i] = add_mat[i,i-1] = 1
            # Two outer diagonals
            if r > 0: add_mat[i-cols,i] = add_mat[i,i-cols] = 1
    m_list=[]
    for i in range(n):
        for j in range (n):
            if add_mat[i,j]==1:
                m_list.append([int(i),int(j),m])
    m_list=np.array(m_list)
    return m_list

def read_input_data(args, mig_file, demography_file, sample_list, anc_pop_sizes, ancpop_list):
    """Reads input data from txt files and returns them in variables."""

    tree_num = np.loadtxt("{}".format(demography_file))

    try:
        migration_rate=float(mig_file)
        mig_list=generate_migration_list(args.row_num, int(tree_num.shape[1]/args.row_num), migration_rate)
    except:
        mig_list = np.loadtxt("{}".format(mig_file), dtype='float')

    cell_num = tree_num.shape[1]
    ngens = tree_num.shape[0]
    if ancpop_list == 'nan':
        ancpop_list = np.zeros(cell_num) + 1

        ancpop_list = ancpop_list.astype(int)
    else:
        # this is a list assigning ancestral pop to each cell of the grid.
        ancpop_list = np.loadtxt("{}".format(ancpop_list), dtype='int')
    #values, counts = np.unique(words, return_counts=True)

    if anc_pop_sizes == 'nan':
        anc_pop_sizes= np.zeros(len(np.unique(ancpop_list))) + 1

        anc_pop_sizes = anc_pop_sizes.astype(int)
    else:
        # this is a list assigning ancestral pop to each cell of the grid.
        anc_pop_sizes = np.loadtxt("{}".format(anc_pop_sizes), dtype='int')
    #values, counts = np.unique(words, return_counts=True)
    print (anc_pop_sizes)


    if sample_list == 'nan':
        sampled_demes = sample_all_cells(tree_num)
    else:
        sampled_demes = np.loadtxt("{}".format(sample_list), dtype='int')

    anc_num = len(anc_pop_sizes)  # what is this for?
    n_ext = cell_num + anc_num
    sampled_demes.sort()
    # all this fun is just to ceck if we don't sample an empty cell. Why 'set'?
    samp = list(set(tree_num[-1][sampled_demes]))

    min_size = 1e-10

    assert tree_num.shape[1] % args.row_num == 0, "Number of cells is not divisible by the number of rows!"
    assert samp[0] > 0, "You cannot sample from cells that are empty!"
    #print ('All is fine')
    ndip = np.zeros(n_ext)

    ndip[sampled_demes] = 2
    ndip = ndip.astype(int)

    return tree_num, mig_list, anc_num, ancpop_list, ngens, cell_num, n_ext, min_size, ndip, anc_pop_sizes

def generate_extended_bm_matrix(n_ext, b_mig_list):
    """Generates extended migration matrix."""
    backw_mig_mat = np.zeros([n_ext, n_ext])
    for i in range(len(b_mig_list[:, 0])):
        # source cell when thinking forward in time
        source_cell = int(b_mig_list[i, 0])
        target_cell = int(b_mig_list[i, 1])
        b_mig_rate = (b_mig_list[i, 2])
        backw_mig_mat[source_cell, target_cell] = b_mig_rate
    return backw_mig_mat

def generate_ext_pop_n(tree_num, cell_num, n_ext, ngens):
    """Returns populations with added non spatial ancestral populations."""
    present_pop = np.zeros(n_ext)
    for i in range(cell_num):
        present_pop[i] = tree_num[ngens - 1, i]
    return present_pop

def generate_backward_mig_list(mig_list, recent_ne, past_ne):
    """Generate backward migration list"""
    b_mig_list = []
    for i in range(len(mig_list[:, 0])):
        # source cell when thinking forward in time
        source_cell = int(mig_list[i, 0])
        # target cell when thinking forward in time
        target_cell = int(mig_list[i, 1])
        if recent_ne[target_cell] > 0:
            bw_mig_rate = mig_list[i, 2] * past_ne[source_cell] / \
                recent_ne[target_cell]  # it take -1 because
        else:
            bw_mig_rate = 0
        # note the opposite direction
        b_mig = [target_cell, source_cell, bw_mig_rate]
        b_mig_list.append(b_mig)
    b_mig_list = np.array(b_mig_list)
    b_mig_list = b_mig_list[np.argsort(b_mig_list[:, 0])]
    return b_mig_list

def mass_migration_to_anc_pop(
        demographic_events,
        delta_t,
        anc_num,
        ngens,
        gen_time,
        anc_pop_sizes,
        ancpop_list,
        cell_num,
        min_size):
    """Mass migration of all remaining populations to ancestral non spatial population."""

    # Creates ancestral populations that were zeros until now
    for i in range(anc_num):
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=delta_t * (ngens) / gen_time,
                initial_size=anc_pop_sizes[i],
                population_id=cell_num + i,
                growth_rate=0))
        # I think this puts initial pop sizes (scaled) to all new ancestral populations

    # And this mass migrates all the other lineages into these ancestral
    # populations
    for i in range(cell_num):
        demographic_events.append(msprime.MassMigration(time=delta_t *
                                                        (ngens) /
                                                        gen_time, source=i, destination=cell_num +
                                                        int(ancpop_list[[i]]) -
                                                        1, proportion=1.0))
        demographic_events.append(msprime.PopulationParametersChange(
            time=delta_t * (ngens) / gen_time, initial_size=min_size,
            population_id=i, growth_rate=0))

    # Migration rates are 0s for the grid, and for the ancestral population
    # it's something small
    for i in range(cell_num):
        for j in range(cell_num):
            if i != j:
                demographic_events.append(
                    msprime.MigrationRateChange(
                        time=delta_t *
                        (ngens) /
                        gen_time,
                        rate=0.0,
                        matrix_index=(
                            i,
                            j)))

    for i in range(anc_num):
        for j in range(anc_num):
            if i != j:
                demographic_events.append(
                    msprime.MigrationRateChange(
                        time=delta_t * (ngens) / gen_time,
                        rate=1e-8,
                        matrix_index=(
                            cell_num + i,
                            cell_num + j)))


##################################### FIRST TIMESTEP #####################
# if a deme has individuals now but in the previous timestep was empty, we need a mass migration (backward in time):
# mass migration from this cell to the neighbors, proportions are weighted by the sizes of the population sizes in the previous generation
# pT, pop_t is to decide about these ratios when mass migration is needed
# when mass migration happens, a population goes empty, so the order of events (may?) matter
# if there is mass migration from i to j, we first set the new size of j, mass migrate the lineages, set the size of i to min_size
# migration matrix (i,j): how many individuals come from cell j to i,
# divided by the size of deme i

# this is time step between presence and the recent past
def define_pop_param_in_first_step(
        demographic_events, mig_list, tree_num, ngens, delta_t, gen_time, min_size):
    """defines demopgraphy and its changes for msprime in the most recent step."""
    cell_num = len(tree_num[0])
    pop_t = np.zeros(np.shape(tree_num[ngens - 1]))
    # this generate a vector - for each cell it contains sum of all neighbours
    # in the penultimate step.
    if mig_list[0, 2] != 0:
        for i in range(len(mig_list[:, 0])):
            cell_origin = int(mig_list[i, 0])
            cell_target = int(mig_list[i, 1])
            pop_t[cell_origin] = pop_t[cell_origin] + \
                tree_num[ngens - 2, cell_target]
        #print (pop_t)

    for i in range(cell_num):
        # this means that a cell previously empty is now colonised.
        if tree_num[ngens - 2, i] == 0 and tree_num[ngens - 1, i] != 0:
            cell_colonised = i
            for j in range(len(mig_list[:, 0])):
                if cell_colonised == int(mig_list[j, 0]):
                    cell_source = int(mig_list[j, 1])
                    demographic_events.append(msprime.PopulationParametersChange(
                        time=delta_t / gen_time, initial_size=max(min_size, tree_num[ngens - 2, cell_source]),
                        population_id=cell_source, growth_rate=0))
                    # this will increase source pop size, so it is filled
                    # before mass migrating the stuff into it
                    # this means that sum od neighbours is 0 - sad stuff
                    if pop_t[cell_colonised] == 0.0:
                        # this is where sad stuff happens. Means that a lineage
                        # gets stuck somewhere until the end of time.
                        prop = 0
                        print(
                            'Some lineages may be stuck, as cell no ',
                            cell_colonised,
                            'has no neighbours in the past that could colonize it. ')
                    else:
                        prop = tree_num[ngens - 2,
                                        cell_source] / pop_t[cell_colonised]

                    demographic_events.append(
                        msprime.MassMigration(
                            time=delta_t / gen_time,
                            source=cell_colonised,
                            destination=cell_source,
                            proportion=min(
                                1,
                                prop)))
 # maybe the problem here is that it is hapenning to early. They should be
 # allowed to move freely first then to mass migrate. Maybe this could be
 # corrected by changing the time for mass migration further back? But this
 # seems fine, they should be able to move around for time delta_t/gen_time...
 # Why do they get stuck?

                    demographic_events.append(msprime.PopulationParametersChange(
                        time=delta_t / gen_time, initial_size=max(min_size, tree_num[ngens - 2, cell_colonised]),
                        population_id=cell_colonised, growth_rate=0))
        else:  # some of these mau be already altered, but it should not matter.
            demographic_events.append(msprime.PopulationParametersChange(
                time=delta_t / gen_time,
                initial_size=max(min_size, tree_num[ngens - 2, i]),
                population_id=i,
                growth_rate=0))

def define_pop_param_for_middle_steps(
        demographic_events, mig_list, tree_num, ngens, delta_t, gen_time, min_size):
    """Defines demography and its changes for msprime in all midlle steps."""

    #t = ngens - 2
    cell_num = len(tree_num[0])
    for t in range(ngens - 2, 0, -1):
        #print('this is time', t)
        this_n = tree_num[t]
        prev_n = tree_num[t - 1]  # this is more in the

       ######################## This modifies migrations #################
        b_mig_list = generate_backward_mig_list(
            mig_list.copy(), this_n.copy(), prev_n.copy())
        # here add stuff for checking the lists! only updating when necessary.

        for i in range(len(b_mig_list[:, 0])):
            # source cell when thinking forward in time -  PROBABLY NOT?!!!
            source_cell = int(b_mig_list[i, 0])
            target_cell = int(b_mig_list[i, 1])
            b_mig_rate = (b_mig_list[i, 2])
            demographic_events.append(msprime.MigrationRateChange(
                time=delta_t * (ngens - t - 1) / gen_time, rate=b_mig_rate, matrix_index=(source_cell, target_cell)))

        t_ago = delta_t * (ngens - t)  # delta_t * (ngens -1- t)
        # maybe this whole chunk could be a special function?
        pop_t = np.zeros(np.shape(tree_num[ngens - 1]))
        # this checks if migration is zero, only if not, this generates a
        # vector - for each cell it contains sum of all neighbours in the
        # previous step.
        if mig_list[0, 2] != 0:
            for i in range(len(mig_list[:, 0])):
                cell_origin = int(mig_list[i, 0])
                cell_target = int(mig_list[i, 1])
                pop_t[cell_origin] = pop_t[cell_origin] + \
                    prev_n[cell_target]  # CHECK THIS


##############
        for i in range(cell_num):
            # this means that a cell previously empty is now colonised.
            if prev_n[i] == 0 and this_n[i] != 0:
                cell_colonised = i
                for j in range(len(mig_list[:, 0])):
                    if cell_colonised == int(mig_list[j, 0]):
                        cell_source = int(mig_list[j, 1])
                        demographic_events.append(
                            msprime.PopulationParametersChange(
                                time=t_ago / gen_time,
                                initial_size=max(
                                    min_size,
                                    prev_n[cell_source]),
                                population_id=cell_source,
                                growth_rate=0))
                        # this will increase source pop size, so it is filled
                        # before mass migrating the stuff into it
                        # this means that sum od neighbours is 0 - sad stuff
                        if pop_t[cell_colonised] == 0.0:
                            # this is where sad stuff happens. Means that a
                            # lineage gets stuck somewhere until the end of
                            # time.
                            prop = 0
                        else:
                            prop = prev_n[cell_source] / pop_t[cell_colonised]

                        demographic_events.append(
                            msprime.MassMigration(
                                time=t_ago / gen_time,
                                source=cell_colonised,
                                destination=cell_source,
                                proportion=min(
                                    1,
                                    prop)))
                        demographic_events.append(
                            msprime.PopulationParametersChange(
                                time=t_ago / gen_time,
                                initial_size=max(
                                    min_size,
                                    prev_n[cell_colonised]),
                                population_id=cell_colonised,
                                growth_rate=0))
            else:  # some of these mau be already altered, but it should not matter.
                demographic_events.append(msprime.PopulationParametersChange(
                    time=t_ago / gen_time, initial_size=max(min_size, prev_n[i]),
                    population_id=i, growth_rate=0))

    return b_mig_list


# these subpopulations are needed for final results, getting T between and
# T within.
def define_subpopulations(ndip):
    """Defines all subpopulations."""
    subpops = []
    counter = 0
    for i in range(
            len(ndip)):  # this generates list of subpopulations, with indexed samples (0-59)
        if ndip[i] != 0:
            elem = []
            for i in range(counter, counter + ndip[i]):
                elem.append(i)
            subpops.append(elem)
            counter = elem[-1] + 1
    return subpops


# this is recalculating BM so it is ready for the last step.
def population_merge(demographic_events, mig_list,
                    initial_n, gen_time, delta_t, ngens):
    """Merge populations at the final point in the past."""
    # this cheats with the backward migration - assumes that the population is
    # fixed, takes the same time point (the last one) and calculates bwmig
    # between those.
    b_mig_list = generate_backward_mig_list(
        mig_list.copy(),
        initial_n,
        initial_n)  # twice the same thing
    t_ago = delta_t * (ngens - 1)
    #print ('new way in generation', t_ago/gen_time)

    for i in range(len(b_mig_list[:, 0])):
        # source cell when thinking forward in time
        source_cell = int(b_mig_list[i, 0])
        target_cell = int(b_mig_list[i, 1])
        b_mig_rate = (b_mig_list[i, 2])
        demographic_events.append(
            msprime.MigrationRateChange(
                time=t_ago / gen_time,
                rate=b_mig_rate,
                matrix_index=(
                    source_cell,
                    target_cell)))

####################################### STARING TO READ IN DATA ##########
##########################################################################

def extract_results(sim_results, subpops, file1):
    """Extracting results of the simulations."""
    tree = sim_results.first()
    sime_res_simple = sim_results.simplify().first()
    pop_length = len(subpops)

    result = []
    result = np.zeros((pop_length, pop_length))

    for i in range(pop_length):
        for j in range(pop_length):
            if i != j:
                result[i,
                       j] = (sime_res_simple.tmrca(subpops[i][1],
                                      subpops[j][1]) + sime_res_simple.tmrca(subpops[i][0],
                                                                subpops[j][0]) + sime_res_simple.tmrca(subpops[i][0],
                                                                                          subpops[j][1]) + sime_res_simple.tmrca(subpops[i][1],
                                                                                                                    subpops[j][0])) / 4
            else:
                result[i, j] = sime_res_simple.tmrca(subpops[i][0], subpops[j][1])

    np.savetxt(file1, result)


def run_gridcoal():
    """Main function to run gridcoal simulations."""

    start_t = time.time()
    print("sim start")

   #!/usr/bin/env python3
    description = '''Simulate a coalescent within the history of Abies expansion.'''

    parser = argparse.ArgumentParser(description=description)
    [args, mig_file, demography_file, serial, delta_t, gen_time, sample_list, anc_pop_sizes,
        ancpop_list, print_deb, dir_name, seed_no] = parse_input_files(parser)

    np.random.seed(seed_no)

    [tree_num,
     mig_list,
     anc_num,
     ancpop_list,
     ngens,
     cell_num,
     n_ext,
     min_size,
     ndip, anc_pop_sizes] = read_input_data(args,
                             mig_file,
                             demography_file,
                             sample_list,
                             anc_pop_sizes,
                             ancpop_list)

    file1 = str(dir_name+'/'+dir_name+'CoalTimes'+str(serial)+'.txt')

    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        print("Directory " , dir_name ,  " Created ")


    if serial ==1:
        out_file = open(dir_name+'/'+dir_name+'Output.txt','w')
        out_file.write(
"""
INPUT FILES
using input file for population sizes:
%s
input sizes: %s
total number of cells: %s
number of known steps: %s
time between steps: %s
min population size: %s
total known time: %s

SAMPLING
using input file for sample sizes: %s
sample sizes taken: %s


MIGRATION
using input file for migration: %s
migration is ["from" "to" "rate"]:
%s

ANCESTRAL POPULATION

ancestral population IDs:
%s
ancestral population sizes:%s
extenden population number: %s


output file name: %s
random seed number: %s
"""% (demography_file, tree_num, cell_num, ngens, delta_t, min_size, (delta_t * ngens), sample_list, ndip, mig_file, mig_list, ancpop_list, anc_pop_sizes, n_ext, file1, seed_no))
        out_file.close()


    ######### DEFINE DEMOGRAPHIC EVENTS THROUGH THE POPULATION HISTORY #######
    ##########################################################################

    population_configurations = []
    demographic_events = []

    define_pop_param_in_first_step(
        demographic_events,
        mig_list.copy(),
        tree_num.copy(),
        ngens,
        delta_t,
        gen_time,
        min_size)
    b_mig_list = define_pop_param_for_middle_steps(
        demographic_events,
        mig_list.copy(),
        tree_num.copy(),
        ngens,
        delta_t,
        gen_time,
        min_size)
    population_merge(
        demographic_events,
        mig_list.copy(),
        tree_num[0],
        gen_time,
        delta_t,
        ngens)

    mass_migration_to_anc_pop(
        demographic_events,
        delta_t,
        anc_num,
        ngens,
        gen_time,
        anc_pop_sizes,
        ancpop_list,
        cell_num,
        min_size)

    ################################# Initial pop config #####################

    present_extended_n = generate_ext_pop_n(tree_num, cell_num, n_ext, ngens)

    b_mig_list = generate_backward_mig_list(
        mig_list.copy(), tree_num[ngens - 1, :], tree_num[ngens - 2, :])
    backw_mig_mat = generate_extended_bm_matrix(n_ext, b_mig_list)

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=ndip[k],
                                        initial_size=max(
                                            min_size, present_extended_n[k]),
                                        growth_rate=0) for k in range(n_ext)]

    ##################################### SIMULATION #########################
    ##########################################################################

    if ((serial == 1) and (print_deb)):
        dem_deb = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=backw_mig_mat,
            demographic_events=demographic_events)
        out_file = open(dir_name+'/'+dir_name+'DemographyDebugger.txt','w')
        dem_deb.print_history(output=out_file)
            #sys.stdout = orig_stdout
        out_file.close()

    sim_results = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=backw_mig_mat)

    ########################CALCULATE EXPECTED PAIRWISE COALESCENCE TIMES#####
    ##########################################################################
    # this is necessary for results analysis
    subpops = define_subpopulations(ndip)
    extract_results(sim_results, subpops, file1)

    end_t = time.time()
    print('Simulation ended in time', end_t - start_t, 'seconds')


##########################################################################
# MAIN FUNCTION STARTING HERE
##########################################################################

run_gridcoal()
