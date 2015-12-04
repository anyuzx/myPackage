from . import dump
from .. import *
import numpy as np
import sys
import time

__all__ = ['ps']

def init_dump(finname,select):
    f = dump(finname)
    f.gettime()
    tsteps = f.timesteps

    sys.stdout.write("First timestep: {}\n".format(str(tsteps[0])))
    sys.stdout.write("Last timestep: {}\n".format(str(tsteps[-1])))
    sys.stdout.flush()

    f.tselect(tsteps[0])
    snap = f.nextSnap()
    natoms = snap.natoms
    index = []

    for item in ['id','x','y','z']:
        index.append(snap.descriptor.index(item))

    index = np.array(index)

    if select == 'all':
        timeselect = tsteps[:]
    else:
        timeselect = tsteps[select[0]:select[1]]

    f.tselect(timeselect)
    sys.stdout.write('time selection is done\n')
    sys.stdout.flush()

    return f, timeselect, index, natoms

def ps(finname,foutname,select = 'all'):
    sys.stdout.write('Analysing file:{}\n'.format(finname))
    sys.stdout.flush()
    f, timeselect, index, natoms= init_dump(finname,select)

    deltat = 0.0
    with open(foutname,'w') as fout:
        for i in range(len(timeselect)):
            t1 = time.time()
            snap = f.nextSnap()
            timestep = snap.time
            sys.stdout.write('Timestep analyzed {}\n'.format(timestep))
            sys.stdout.flush()
            snap = snap.atoms[:,index]
            snap = snap[snap[:,0].argsort()]
            dist = snap[:,1:]
            contact_prob = analysis.polymertool.contactprob(dist,2.0)
            fout.write('Timestep: {}\n'.format(timestep))
            for j in range(natoms-1):
                fout.write(str(j+1).ljust(10)+'{:.10e}'.format(contact_prob[j])+'\n')
            fout.flush()
            t2 = time.time()
            deltat += t2-t1
            sys.stdout.write('Estimated finished time left: {}\n'.format((deltat/(i+1))*(len(timeselect)-i)))
            sys.stdout.flush()



    sys.stdout.write('Total number of snapshots analyzed: {}\n'.format(len(timeselect)))
    sys.stdout.flush()
