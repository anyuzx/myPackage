# created by Guang Shi 09.13.2015
from . import dump
from .. import *
import numpy as np
import sys
import time

__all__ = ['cmap']

def init_dump(finname,select):
	f = dump(finname)
	f.gettime()
	tsteps = f.timesteps

	sys.stdout.write('First timestep: {}\n'.format(str(tsteps[0])))
	sys.stdout.write('Last timestep: {}\n'.format(str(tsteps[-1])))
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

def cmap(finname,foutname,rc,select='all',split='none'):
	sys.stdout.write('Analysing file:{}\n'.format(finname))
	sys.stdout.flush()
	f, timeselect, index, natoms = init_dump(finname,select)

	if split == 'none':
		pass
	else:
		timesplit = np.array_split(range(len(timeselect)),split)

	deltat = 0.0
	contact_map = np.zeros((natoms,natoms),dtype=np.int32)
	if split == 'none':
		for i in range(len(timeselect)):
			t1 = time.time()
			snap = f.nextSnap()
			timestep = snap.time
			sys.stdout.write('Timestep analyzed {}\n'.format(timestep))
			sys.stdout.flush()
			snap = snap.atoms[:,index]
			snap = snap[snap[:,0].argsort()]
			dist = pdist(snap[:,1:])
			contact_map += analysis.polymertool.contactmap(dist,rc,mode='half')
			t2 = time.time()
			deltat += t2-t1
			sys.stdout.write('Estimated finished time left: {}\n'.format((deltat/(i+1))*(len(timeselect)-i)))
			sys.stdout.flush()
		with open(foutname,'w') as fout:
			for item in contact_map:
				for c in item:
					fout.write(str(c)+'  ')
				fout.write('\n')
	else:
		k = 0
		for i in range(len(timeselect)):
			t1 = time.time()
			snap = f.nextSnap()
			timestep = snap.time
			sys.stdout.write('Timestep analyzed {}\n'.format(timestep))
			sys.stdout.flush()
			snap = snap.atoms[:,index]
			snap = snap[snap[:,0].argsort()]
			dist = pdist(snap[:,1:])
			contact_map += analysis.polymertool.contactmap(dist,2.0,mode='half')
			if i == timesplit[k][-1]:
				with open(foutname+'_'+str(i),'w') as fout:
					for item in contact_map:
						for c in item:
							fout.write(str(c)+'  ')
						fout.write('\n')
				k += 1
			t2 = time.time()
			deltat += t2-t1
			sys.stdout.write('Estimated finished time left: {}\n'.format((deltat/(i+1))*(len(timeselect)-i)))
			sys.stdout.flush()

	sys.stdout.write('Total number of snapshots analyzed: {}\n'.format(len(timeselect)))
	sys.stdout.flush()
