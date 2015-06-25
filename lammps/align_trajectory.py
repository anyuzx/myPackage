from . import dump
from .. import *
import numpy as np
import sys
import time

__all__ = ['Align']

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

	return f, timeselect, index

def Align(finname,foutname,select = 'all'):
	f, timeselect, index = init_dump(finname,select)

	snap0 = f.nextSnap()
	timestep = snap0.time
	natoms = snap0.natoms
	box = snap0.box
	snap0 = snap0.atoms[:,index]
	snap0 = snap0[snap0[:,0].argsort()]
	snap0 = np.array(snap0[:,1:],dtype = np.float32)
	snap0 = snap0 - np.mean(snap0,axis=0)

	deltat = 0.0
	with open(foutname,'w') as fout:
		fout.write('ITEM: TIMESTEP\n')
		fout.write('{}\n'.format(timestep))
		fout.write('ITEM: NUMBER OF ATOMS\n')
		fout.write('{}\n'.format(natoms))
		fout.write('ITEM: BOX BOUNDS pp pp pp\n')
		fout.write('{} {}\n'.format(box[0,0],box[0,1]))
		fout.write('{} {}\n'.format(box[1,0],box[1,1]))
		fout.write('{} {}\n'.format(box[2,0],box[2,1]))
		fout.write('ITEM: ATOMS id x y z\n')

		j = 0
		for atom in snap0:
			fout.write('{} {:.4f} {:.4f} {:.4f}\n'.format(j+1,atom[0],atom[1],atom[2]))
			j += 1

		for i in range(len(timeselect)-1):
			t1 = time.time()
			snap = f.nextSnap()
			box = snap.box
			timestep = snap.time
			sys.stdout.write('Timestep analyzed {}\n'.format(timestep))
			sys.stdout.flush()
			snap = snap.atoms[:,index]
			snap = snap[snap[:,0].argsort()]
			snap = np.array(snap[:,1:],dtype=np.float32)

			# calculate the rotate matrix
			newsnap = analysis.generaltool.optimal_rotate(snap,snap0)
			newsnap = np.array(newsnap)

			fout.write('ITEM: TIMESTEP\n')
			fout.write('{}\n'.format(timestep))
			fout.write('ITEM: NUMBER OF ATOMS\n')
			fout.write('{}\n'.format(natoms))
			fout.write('ITEM: BOX BOUNDS pp pp pp\n')
			fout.write('{} {}\n'.format(box[0,0],box[0,1]))
			fout.write('{} {}\n'.format(box[1,0],box[1,1]))
			fout.write('{} {}\n'.format(box[2,0],box[2,1]))
			fout.write('ITEM: ATOMS id x y z\n')

			j = 0
			for atom in newsnap:
				fout.write('{} {:.4f} {:.4f} {:.4f}\n'.format(j+1, atom[0], atom[1], atom[2]))
				j += 1

			t2 = time.time()
			deltat += t2 - t1
			sys.stdout.write('Estimated finished time left: {}\n'.format((deltat/(i+1))*(len(timeselect)-i)))
			sys.stdout.flush()

		sys.stdout.write('Total number of snapshots analyzed: {}\n'.format(len(timeselect)))
		sys.stdout.flush()
		sys.stdout.write('Finished\n')
		sys.stdout.flush()







