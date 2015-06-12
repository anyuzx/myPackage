from math import *
from ... import create
import numpy as np

# write the lammps data file for lattice chain
# argument list:
#       1st: number of monomers
#       2nd: bond length
#       3rd: Volume Exclusion
#       4th: Number of Iterations
def writefile(foutname,*argv,**kwargs):
    chain = create.chain(*argv,**kwargs)
    chain_config, box_dimension = chain.create_chain()
    VolumeExclusion = kwargs.get('VolumeExclusion',True)
    N = len(chain_config)
    mass = 1.0

    with open(foutname,'w') as f:
        f.write('Data file of lattice chain. Volume Exclusion: '+str(VolumeExclusion)+'\n')
        f.write(str(N)+'  atoms   # number of monomers\n')
        f.write(str(N-1)+'  bonds   # number of bonds between monomers\n')
        f.write('\n')

        f.write("1  atom types     # number of atoms types\n")
        f.write("1  bond types    # number of bond types\n")
        f.write("\n")

        f.write(str("%.4f" % box_dimension[0,0]) + " " + str("%.4f" % box_dimension[0,1]) + " xlo " + "xhi\n")
        f.write(str("%.4f" % box_dimension[1,0]) + " " + str("%.4f" % box_dimension[1,1]) + " ylo " + "yhi\n")
        f.write(str("%.4f" % box_dimension[2,0]) + " " + str("%.4f" % box_dimension[2,1]) + " zlo " + "zhi\n")

        f.write("\n")
        f.write("Masses\n\n")
        f.write(str(1).ljust(10)+"%.4f" % mass + "\n")
        f.write("\n")

        f.write("Atoms\n\n")
        index_monomer = 1
        for monomer in chain_config:
            f.write(str(index_monomer).ljust(10) + str(1).ljust(10) + str(1).ljust(10) + \
                    ("%.4f" % monomer[0]).ljust(15) + ("%.4f" % monomer[1]).ljust(15) + ("%.4f" % monomer[2]).ljust(15) + '\n')
            index_monomer += 1

        f.write("\n")
        f.write("Bonds\n\n")
        for i in range(1,N):
            f.write(str(i).ljust(10) + str(1).ljust(10) + str(i).ljust(10) + str(i+1).ljust(10) + '\n')
