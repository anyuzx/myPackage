import myPackage as mp
import numpy as np


chain = mp.lammps.data()
loop_datafile = '/Users/Shi/Downloads/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt'
epigenetic_state_datafile = '/Users/Shi/Desktop/LEM/GM12878/wgEncodeBroadHmmGm12878HMM.bed'
chain.writetofile('LEM_chain','lem_chain_chr5.dat',bond_length=1.13,straight_line=False,\
	              loop_filename=loop_datafile, epigenetic_filename = epigenetic_state_datafile, chrom='Chr5',\
	              start=145870001,end=157870001,bin_size=1200, data_source = 'Encode_Chrom_State')


'''
chain2 = mp.lammps.data()
loop_base = [[201,400],[601,800],[1001,1200],[1401,1600]]
state1 = np.concatenate((np.arange(1,201), np.arange(601,1001)))
state1 = np.concatenate((state1, np.arange(1401,1601)))
state2 = np.concatenate((np.arange(201,601), np.arange(1001,1401)))
state2 = np.concatenate((state2, np.arange(1601,1801)))
epigenetic_state = {'state1':state1,'state2':state2}

chain2.writetofile('LEM_chain','lem_chain_manual.dat',bond_length=1.13,straight_line=False,binding_protein=[200,200],loop_base=loop_base,epigenetic_state=epigenetic_state)
'''