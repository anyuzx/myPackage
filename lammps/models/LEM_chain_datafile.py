from ... import create
import numpy as np
import sys

# write the lammps data file for loop/epigenetic(binding proteins) chain
# argument list:

# This function convert the epigenetic state into the atom type in the data file
def convert_state(state):
    state = state.lower()
    if state == 'state1' or state == 'state 1':
        return 1
    elif state == 'state2' or state == 'state 2':
        return 2
    else:
        sys.stdout.write('Undefined monomer epigenetic state\n')
        raise ValueError

def writefile(foutname,bond_length = 1.0, straight_line = False, binding_protein = None, *argv,**kwargs):
    chain = create.lem_chain(*argv, **kwargs)
    chain.create_config(bond_length = bond_length, straight_line = straight_line)
    if binding_protein is not None:
        chain.add_bind_protein(binding_protein)

    mass = 1.0
    with open(foutname, 'w') as f:
        f.write("Data file for LEM polymer chain   Volume Exclusion:" + str('True')+"   Chain persistence:"+str('False')+"\n\n")
        if binding_protein is None:
            f.write(str(chain.chain_length)+'  atoms   # number of monomers\n')
        else:
            f.write(str(chain.chain_length+chain.number_bind_atoms)+'  atoms   # number of monomers\n')
        f.write(str(chain.chain_length-1+len(chain.loop_base))+'  bonds   # number of bonds\n')
        f.write('\n')
        if binding_protein is None:
            f.write(str(len(chain.epigenetic_state.keys()))+'  atom types   # number of epigenetic type\n')
        else:
            f.write(str(len(chain.epigenetic_state.keys())+chain.number_type_bind_atoms)+'  atom types   # number of epigenetic type\n')
        f.write('2  bond types    # number of bond types\n')
        f.write('\n')
        f.write(str("%.4f" % chain.box_dimension[0,0]) + " " + str("%.4f" % chain.box_dimension[0,1]) + " xlo" + " xhi\n")
        f.write(str("%.4f" % chain.box_dimension[1,0]) + " " + str("%.4f" % chain.box_dimension[1,1]) + " ylo" + " yhi\n")
        f.write(str("%.4f" % chain.box_dimension[2,0]) + " " + str("%.4f" % chain.box_dimension[2,1]) + " zlo" + " zhi\n")
        f.write('\n')
        f.write('Masses\n\n')
        if binding_protein is None:
            for i in range(len(chain.epigenetic_state.keys())):
                f.write(str(i+1).ljust(10)+"%.2f" % mass + '\n')
        else:
            for i in range(len(chain.epigenetic_state.keys())+chain.number_type_bind_atoms):
                f.write(str(i+1).ljust(10)+"%.2f" % mass + '\n')
        f.write('\n')
        f.write('Atoms\n\n')
        for index, atom in chain.config.iterrows():
            f.write(str(int(atom['index'])).ljust(10)+str(1).ljust(10)+str(convert_state(atom['state'])).ljust(10)+("%.4f" % atom['x']).ljust(15)+("%.4f" % atom['y']).ljust(15)+("%.4f" % atom['z']).ljust(15)+"\n")
        if binding_protein is not None:
            for index, atom in chain.bind_atoms_config.iterrows():
                f.write(str(int(atom['atom_index']+chain.chain_length)).ljust(10)+str(1+int(atom['atom_index'])).ljust(10)+str(int(atom['atom_type'])+len(chain.epigenetic_state.keys())).ljust(10)+("%.4f" % atom['x']).ljust(15)+("%.4f" % atom['y']).ljust(15)+("%.4f" % atom['z']).ljust(15)+"\n")
        f.write('\n')
        f.write('Bonds\n\n')
        for i in range(1, chain.chain_length):
            f.write(str(i).ljust(10)+str(1).ljust(10)+str(i).ljust(10)+str(i+1).ljust(10)+"\n")
        for item in chain.loop_base:
            i += 1
            f.write(str(i).ljust(10)+str(2).ljust(10)+str(item[0]).ljust(10)+str(item[1])+'\n')

