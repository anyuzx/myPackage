import numpy as np
import pandas as pd
import create_atoms as ca
import sys

# ------------------------------------------------------------
# give initial configuration for loop/epigenetic chain with
# or without binding proteins
# ------------------------------------------------------------

# ------------------------------------------------------------
# Usage:
#    1.Generate a loop/epigenetic chain using the given loop 
#      anchors and epigenetic state data.
#       myPackage.create.lem_chain(loop_base = loop_base_data, epigenetic_state = epigenetic_state_data)
#    2.Generate a loop/epigenetic chain using the experiment data
#       myPackage.create.lem_chain(filename = file_name, chrom=chrom_name, start=start_loci, end=end_loci, bin_size=binsize)
#
#    loop_base: two dimensional array for loop bases. Example: [[1,3],[5,9],[10,20]]
#    epigenetic_state: a dictionary for monomer states. Example: {'state1': [1,2,5,10],'state2':[3,4,6,7,8],'state3':[9]}
#    filename: directory for the CTCF motif data file.(2014 Cell paper)
#    chrom: chromosome name. Example: 'Chr1', 'Chr2' , 'ChrX', etc
#    start: the starting point of sequence.
#    end: the ending point of sequence.
#    bin_size: the number of base pairs each bead in the chain represents

# ------------------------------------------------------------
# define several function used in the code

# this function determine where a pair is in the group or not
# Example: suppose we have loop pairs [1,5],[1,10],[2,5],[5,10],[12,18],[18,30]
# In the above configurations, the loop [1,5],[1,10],[2,5],[5,10] are in the same "group".
# And [12,18],[18,30] are in the same group.
# Each group form a connected graph. Different groups have no links between them.
def same_group(group, pair):
    for item in pair:
        if item in np.asarray(group).flatten():
            return True
    return False


# this function convert the same group of loops into pairs
# Example: suppose we have two loops [1,5] and [1,10]. We generate
# the initial configuration is such way: bead '1' is 'looped' to bead '5'
# and bead '6' is looped to bead '10'. 
# This function is used to generate pairs [1,5] and [6,10] given the input [1,5,10]
def convert_to_pair(group):
    temp = np.asarray(group)
    temp = np.sort(np.unique(temp))
    number = len(temp)
    pair = []
    pair.append([temp[0], temp[1]])
    for i in range(1,number-1):
        pair.append([temp[i]+1, temp[i+1]])
    return np.asarray(pair)


# Given the input of original loop base data, this function
# group the loops which belong to the same group.
# Example: loop_base = [[1,5],[1,10],[5,10],[12,18],[18,30]]
#          group(loop_base) returns [[1,5,10],[12,18,30]]
def group(loop_base):
    group = [[list(loop_base[0])]]
    loop_base_temp = np.copy(loop_base[1:])
    check = False

    while len(loop_base_temp) > 0:
        check = False
        for pair in loop_base_temp:
            index = np.where((loop_base_temp[:,0] == pair[0]) & (loop_base_temp[:,1] == pair[1]))[0]
            if same_group(group[-1],pair):
                group[-1].append(list(pair))
                loop_base_temp = np.delete(loop_base_temp, index,axis = 0)
                check = True
        if not check:
            group.append([])
            group[-1].append(list(loop_base_temp[0]))
            loop_base_temp = np.delete(loop_base_temp, 0,axis = 0)
    
    return group


# This function generate the modified loop data given the original loop data
# Example: loop_base = [[1,5],[1,10],[5,10],[12,18],[18,30]]
#          give_pair(loop_base) returns [[1,5],[6,10],[12,18],[19,30]]
# The reason for doing this is to generate the initial configuration easily
# But the actual bond information in the data file is still the original loop data.
def give_pair(loop_base):
    f_group = group(loop_base)
    
    pair = []
    i = 0
    for item1 in f_group:
        for item2 in convert_to_pair(item1):
            pair.append(item2)
        i += 1
    
    return np.asarray(pair)



# This function give the loop data from the experimental data
# Argument:
#           filename: directory of data file
#           chrom: name of chromosome
#           start: starting point of the sequence
#           end: ending point of the sequence
#           bin_size: the number of base pairs each bead represent
def get_loop(filename, chrom, start, end, bin_size):
    # Read experimental data from Rao, Lieberman 2014 Cell
    # GM12878 Cell Line
    motif_data = np.loadtxt(filename,dtype='S',skiprows=1,usecols=(0,3,20,21,23,24,25,26,28,29))
    # For our purpose, only keep the loops has two unique CTCF BDS at its both anchors
    motif_data = motif_data[np.where((motif_data[:,5]=='u') & (motif_data[:,9]=='u'))]
    # Keep the information of Chromosome #, position of both motifs, and their orientations
    motif_data = motif_data[:,np.array([0,1,2,3,4,6,7,8])]

    # Read experimental data of subcompartment from Rao, Lieberman 2014 Cell
    # GM12878 Cell Line
    #subcmptmnt_data = np.loadtxt('/Users/Shi/Downloads/GSE63525_GM12878_subcompartments.bed',dtype='S',usecols=(0,1,2,3,4))
    #subcmptmnt_data = subcmptmnt_data[np.where(subcmptmnt_data[:,3]!='NA')]

    # first delete the loop which are not in the same chromosome
    motif_data = motif_data[np.where(motif_data[:,0] == motif_data[:,1])]

    new_motif_data = np.copy(motif_data)
    for row in np.arange(motif_data.shape[0]):
        if motif_data[row,4] == 'p':
            new_motif_data[row,4] = 1
        elif motif_data[row,4] == 'n':
            new_motif_data[row,4] = -1
        
        if motif_data[row,7] == 'p':
            new_motif_data[row,7] = 1
        elif motif_data[row,7] == 'n':
            new_motif_data[row,7] = -1

    # Store each chromosome in a dictionary
    Chr_dic = {}
    for row in np.arange(new_motif_data.shape[0]):
        Chr_name = new_motif_data[row,0]
        if 'Chr'+Chr_name not in Chr_dic.keys():
            Chr_dic['Chr'+Chr_name] = []
            Chr_dic['Chr'+Chr_name].append(new_motif_data[row,2:])
        else:
            Chr_dic['Chr'+Chr_name].append(new_motif_data[row,2:])

    for key in Chr_dic.keys():
        Chr_dic[key] = np.array(Chr_dic[key],dtype=np.int)

    data_exp = Chr_dic[chrom][np.where((Chr_dic[chrom][:,0]>=start) & (Chr_dic[chrom][:,3]<=end))]
    loop_exp = np.asarray(zip((data_exp[:,0]+data_exp[:,1])/2.0,(data_exp[:,3]+data_exp[:,4])/2.0))
    loop_sim = (np.int_(loop_exp)-start)/bin_size

    n = (int(end)-int(start))/bin_size
    return loop_sim, n

# define function to convert the subcompartment or chromatin state to number
# if providing the Encode Chromatin state data, state number <=1 --> open chromatin
#                                                                    state 1
# if providing the subcompartment state data, A subcompartment --> open chromatin
#                                                                  state 1
def convert_state_2_number(state, data_source):
    if data_source == 'Encode_Chrom_State':
        for i in range(1,16):
            if str(i) == state.split('_')[0]:
                if i <= 11:
                    return 1  # 1 represents state 1. open chromatin
                else:
                    return 2  # 2 represents state 2. close chromatin
    elif data_source == 'Cell_Subcompartment':
        if state == 'NA':
            return np.random.choice([1,2],1)[0]
        elif 'A' in state:
            return 1
        elif 'B' in state:
            return 2


# This function gives epigenetic state from the experimental data
# Argument:
#           filename: directory for experimental data file.
#           chrom: name of chromosome
#           start: starting point of the sequence
#           end: ending point of the sequence
#           bin_size: the number of base pairs each bead represent
#           data_source: Now two kinds of files are accepted:
#                        1. Cell Rao, etc paper
#                        2. Encode Chromatin state
def get_state(filename, chrom, start, end, bin_size, data_source):
    # pandas dataframe for data
    data = pd.read_csv(filename, names = ['chrom','chromStart','chromEnd','name','score',\
                       'strand','thickStart','thickEnd','itemRgb'], header = None, \
                        delim_whitespace=True)

    if data_source != 'Encode_Chrom_State' and data_source != 'Cell_Subcompartment':
        sys.stdout.write("Please specify correct data source type:\n 1. 'Encode_Chrom_State'\n 2. 'Cell_Subcompartment'\n")
        raise

    chrom = chrom.lower()

    # convert the experimental subcompartment or chromatin state to simulation \
    # monomer state (monomer atom type)
    data['name'] = data['name'].fillna('NA')
    data['name'] = data['name'].apply(convert_state_2_number, args=(data_source,))

    simulation_data = np.linspace(start, end, int((end-start)/bin_size)+1, dtype=np.int)
    # 10000000 is just a number to make sure that selected experiment data starts from
    # before the actual start point of sequence specified in the function
    experiment_data = data[(data.chrom == chrom) & (data.chromStart >= start - 10000000) & \
                           (data.chromEnd <= end + 10000000)]
    experiment_start_end_state_array = experiment_data[['chromStart', 'chromEnd', \
                                                        'name']].values

    state_of_monomer_dic = {'state 1':[], 'state 2':[]} # initialize the dictionary
    # loop the monomer to check the state of the monomer
    for i in range(len(simulation_data)-1):
        monomer_start = simulation_data[i]
        monomer_end = simulation_data[i+1]
        #np.searchsorted(data[data.chrom == chrom][[]].values, monomer_start)
        state_number_array = np.zeros(2)  # 2 means two distinct epeigenetic states.\
                                          # more states may be implemented in the future.
        for state in experiment_start_end_state_array:
            # case where the subcompartment is inside a monomer
            if monomer_start < state[0] and monomer_end > state[1]:
                state_number_array[state[2]-1] += state[1] - state[0]
                break
            # case where the monomer is inside the subcompartment
            elif monomer_start >= state[0] and monomer_end <= state[1]:
                state_number_array[state[2]-1] += monomer_end - monomer_start
                break
            # case where the left part of monomer is in the subcompartment
            elif monomer_start >= state[0] and monomer_start <= state[1]:
                state_number_array[state[2]-1] += state[1] - monomer_start
                break
            # case where the right part of monomer is in the subcompartment
            elif monomer_end >= state[0] and monomer_end <= state[1]:
                state_number_array[state[2]-1] += monomer_end - state[0]
                break
        
        # get the state of the monomer
        # criterion: whichever the state has more base pairs in one monomer
        state_of_monomer = np.argmax(state_number_array) + 1

        # create the dictionary of state of monomers
        if state_of_monomer == 1:
            state_of_monomer_dic['state 1'].append(i+1)
        elif state_of_monomer == 2:
            state_of_monomer_dic['state 2'].append(i+1)
    
    return state_of_monomer_dic


def get_loop_epigenetic(loop_filename = None, epigenetic_filename = None, chrom = None, start = None, end = None, \
                        bin_size = None, data_source = None):
    if loop_filename == None and epigenetic_filename == None:
        sys.stdout.write('At least one of loop or epigenetic data should be provided\n')
        raise
    elif loop_filename != None and epigenetic_filename == None:
        loop_sim, n = get_loop(filename = loop_filename, chrom = chrom, start = start, end = \
                               end, bin_size = bin_size)
        n = (int(end)-int(start))/bin_size
        epigenetic_state = {'state1': np.arange(1,n+1)}
        return loop_sim, n, epigenetic_state
    elif loop_filename != None and epigenetic_filename != None:
        if data_source == None:
            sys.stdout.write('please specify the type of subcompartment/chromatin state data.\n')
            raise
        else:
            loop_sim, n = get_loop(filename = loop_filename, chrom = chrom, start = start, end = \
                               end, bin_size = bin_size)
            epigenetic_state = get_state(filename = epigenetic_filename, chrom = chrom, start = start,\
                                         end = end, bin_size = bin_size, data_source = data_source)
            return loop_sim, n, epigenetic_state
    elif loop_filename == None and epigenetic_filename != None:
        sys.stdout.write("haven't implemented only epigenetic state no loop feature yet.")
        raise


# ------------------------------------------------------------

__all__ = ['lem_chain']

# ------------------------------------------------------------
# loop epigenetic model
# take the loop bases and epigenetic state as input
# output the initial configuration to run simulation
# arguments:
#   loop_base: 2D array of loop base [[1,2],[3,4],[5,6]]
#   epigenetic_state: dictionary of epigenetic_state of each monomers:
#       {'state 1':[1,2,3], 'state 2':[3,4,5]}
class lem_chain:

    # ---------------------------------------
    def __init__(self, loop_base = None, epigenetic_state = None, loop_filename = None, epigenetic_filename = None, \
                 chrom = None, start = None, end = None, bin_size = None, data_source = None):
        if loop_base != None and epigenetic_state != None and loop_filename == None and epigenetic_filename == None and \
                        chrom == None and start == None and end == None and bin_size == None and data_source ==None:
            self.loop_base = np.asarray(loop_base)
            self.loop_base_simulation = np.copy(self.loop_base)
            self.epigenetic_state = epigenetic_state.copy()
            self.chain_length = 0
            for key in self.epigenetic_state.keys():
                self.epigenetic_state[key] = np.asarray(self.epigenetic_state[key])
                self.chain_length += len(self.epigenetic_state[key])
        elif loop_base == None and epigenetic_state == None and chrom != None and start != None and end != None and \
                                         bin_size != None and (loop_filename != None or epigenetic_filename != None):
            self.loop_base, self.chain_length, self.epigenetic_state = get_loop_epigenetic(loop_filename, epigenetic_filename,\
                                                                       chrom, start, end, bin_size, data_source)
            self.loop_base_simulation = give_pair(self.loop_base)
            for key in self.epigenetic_state.keys():
                self.epigenetic_state[key] = np.asarray(self.epigenetic_state[key])
        else:
            raise 'The arguments are not corrected. Please check the arguments\n'

        self.loop_base.sort()
        self.loop_base_simulation.sort()      
        self.bind_atoms_config = None

    def create_config(self, bond_length=1.0, straight_line = False):
        self.bond_length = bond_length
        if not straight_line:
            rest = self.chain_length - np.sum(self.loop_base_simulation[:,1] - self.loop_base_simulation[:,0] -1)
            loop_index = []
            for i in range(len(self.loop_base_simulation)):
                loop_index = np.concatenate((loop_index, np.arange(self.loop_base_simulation[i,0], self.loop_base_simulation[i,1]-1)))
            no_loop_index = np.delete(np.arange(1,self.chain_length+1), loop_index)
            self.config = np.asarray([[i*bond_length, 0, 0] for i in range(rest)])
            self.config = np.hstack((no_loop_index.reshape((len(no_loop_index),1)), self.config))

            for i in range(len(self.loop_base_simulation)):
                temp1 = int(self.loop_base_simulation[i,1] - self.loop_base_simulation[i,0] - 1)/2
                temp2 = int(self.loop_base_simulation[i,1] - self.loop_base_simulation[i,0] - 1) - temp1

                temp_config1 = np.asarray([[self.loop_base_simulation[i,0]+1+j, self.config[np.where(self.config[:,0] == self.loop_base_simulation[i,0])][0,1],(j+1)*bond_length,0] for j in range(temp1)])
                temp_config2 = np.asarray([[self.loop_base_simulation[i,0]+temp1+j+1,self.config[np.where(self.config[:,0] == self.loop_base_simulation[i,1])][0,1],(temp1+1-j)*bond_length,0] for j in range(temp2)])
                self.config = np.vstack((self.config, temp_config1))
                self.config = np.vstack((self.config, temp_config2))

            self.config = self.config[self.config[:,0].argsort()]
        else:
            self.config = np.asarray([[i+1, i*self.bond_length, 0, 0] for i in range(self.chain_length)])
        
        self.state = range(self.chain_length)
        for key in self.epigenetic_state.keys():
            for i in self.epigenetic_state[key]:
                self.state[i-1] = key

        self.config = pd.DataFrame({'index':self.config[:,0],
                           'state':self.state,
                           'x':self.config[:,1],
                           'y':self.config[:,2],
                           'z':self.config[:,3]})

        lx_min = np.min(self.config['x'].values)
        lx_max = np.max(self.config['x'].values)
        ly_min = np.min(self.config['y'].values)
        ly_max = np.max(self.config['y'].values)
        lz_min = np.min(self.config['z'].values)
        lz_max = np.max(self.config['z'].values)
        self.box_dimension = np.array([[lx_min-20.0*bond_length, lx_max+20.0*bond_length],[ly_min-20.0*bond_length, ly_max+20.0*bond_length],[-(ly_max-ly_min)/2.0-20.0*bond_length, (ly_max-ly_min)/2.0+20.0*bond_length]])
    
    # add binding proteins
    # make the simulation box cubic and periodic boundary
    def add_bind_protein(self, atom_array):
        self.bind_atoms_config = ca.create_atoms(atom_array, self.box_dimension, self.bond_length, initial = self.config.values[:,2:])
        self.number_bind_atoms = np.sum(np.array(atom_array))
        self.number_type_bind_atoms = len(atom_array)

    # plot
    def plot(self):
        try:
            self.config
        except NameError:
            sys.stdout.write('ERROR: No existing chain found. Please create the chain first.\n')

        import matplotlib.pyploy as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplots(111,projection='3d')
        if self.bind_atoms_config == None:
            ax.scatter(chain.config['x'].values, chain.config['y'].values, chain.config['z'].values, c='r')
        else:
            ax.scatter(chain.config['x'].values, chain.config['y'].values, chain.config['z'].values, c='r')
            ax.scatter(chain.bind_atoms_config['x'].values, chain.bind_atoms_config['y'].values, chain.bind_atoms_config['z'].values, c='b')
