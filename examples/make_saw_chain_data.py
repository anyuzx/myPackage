import myPackage as mp
import numpy as np

n = 10000

chain = mp.lammps.data()
chain.writetofile('lattice_chain','IDEAL_harmonic_n'+str(n)+'.dat',n,l=1.0,VolumeExclusion=False,NumberIteration=1000000)
