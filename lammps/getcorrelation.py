from . import dump
from .. import *
import numpy as np
import sys
import time

__all__ = ['ProcessDumpCorrelation']

def init_dump(finname,mode,select):
	f = dump(finname)
	f.gettime()
	tsteps = f.timesteps

	sys.stdout.write("First timestep: {}\n".format(str(tsteps[0])))
	sys.stdout.write("Last timestep: {}\n".format(str(tsteps[-1])))
	sys.stdout.flush()

	f.tselect(tsteps[0])
	snap = f.nextSnap()
	box = snap.box
	box_size = np.array([np.abs(box[0,0] - box[0,1]),np.abs(box[1,0] - box[1,1]),np.abs(box[2,0] - box[2,1])])
	natoms = snap.natoms
	index = []

	if 'VACF' in mode:
		for item in ['id','x','y','z','vx','vy','vz']:
			index.append(snap.descriptor.index(item))
	else:
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

	return f,timeselect,index

def ProcessDumpCorrelation(finname,foutname,mode,ave_freq,window,select='all'):
	sys.stdout.write('Analysing file:{}\n'.format(finname))
	sys.stdout.flush()
	func_lst = []
	for item in mode:
		if item == 'MSD':
			func_lst.append(analysis.generaltool.msd)
		elif item == 'VACF':
			func_lst.append(analysis.generaltool.vacf)
		else:
			raise NameError('No available correlation mode\n')

	f,timeselect,index = init_dump(finname,mode,select)
	#if timeselect[-1]%time_lapse == 0:
	#	time_init_lst = np.linspace(0,timeselect[-1]-time_lapse,(timeselect[-1]/time_lapse),dtype=int)
	#else:
	#	raise ValueError('Number of Snapshots analyzed must be a multiple of time_lapse')
	time_init_lst = np.linspace(0,ave_freq*(len(timeselect)/ave_freq),len(timeselect)/ave_freq+1,dtype=int)
	deltat = 0
	cf_data = []
	cf_data_std = []
	for item in mode:
		cf_data.append([])
		cf_data_std.append([])
		cf_data[-1] = [np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]) for i in range(window+1)]
		cf_data_std[-1] = [[] for i in range(window+1)]
	snap_init_lst = []

	for i in range(len(timeselect)):
		t1 = time.time()
		snap = f.nextSnap()
		timestep = snap.time
		sys.stdout.write('Timestep analyzed {}\n'.format(timestep))
		sys.stdout.flush()
		snap = snap.atoms[:,index]
		snap = np.array(snap)
		snap = snap[snap[:,0].argsort()]
		snap = np.float32(snap)

		# store the snap as a initial time if the timestep 
		# of the snap is in time_init_lst:
		if i in time_init_lst:
			snap_init_lst.append(snap)

		for j in range(len(snap_init_lst)):
			for k in range(len(mode)):
				if mode[k] == 'MSD':
					data_temp = func_lst[k](snap[:,1:4],snap_init_lst[j][:,1:4])
				elif mode[k] == 'VACF':
					data_temp = func_lst[k](snap[:,4:7],snap_init_lst[j][:,4:7])
				if i - time_init_lst[j] <= window:
					cf_data[k][i - time_init_lst[j]] += np.array([timeselect[i]-timeselect[time_init_lst[j]],data_temp[0],data_temp[1],data_temp[2],(data_temp**2)[0],(data_temp**2)[1],(data_temp**2)[2],1])

		t2 = time.time()
		deltat += t2 - t1
		sys.stdout.write('Estimated finished time left: {}\n'.format((deltat/(i+1))*(len(timeselect)-i)))
		sys.stdout.flush()

	sys.stdout.write('Total number of snapshots analyzed: {}\n'.format(len(timeselect)))
	sys.stdout.flush()
	sys.stdout.write('Start to write data into file\n')
	sys.stdout.flush()

	output  = []
	for k in range(len(mode)):
		output.append([])
		for i in range(len(cf_data[k])):
			cf_data[k][i] = np.array(cf_data[k][i])
                        temp_std1 = np.sqrt(cf_data[k][i][4]/cf_data[k][i][-1] - (cf_data[k][i][1]/cf_data[k][i][-1])**2)
                        temp_std2 = np.sqrt(cf_data[k][i][5]/cf_data[k][i][-1] - (cf_data[k][i][2]/cf_data[k][i][-1])**2)
                        temp_std3 = np.sqrt(cf_data[k][i][6]/cf_data[k][i][-1] - (cf_data[k][i][3]/cf_data[k][i][-1])**2)
			cf_data_std[k][i] = np.array([temp_std1,temp_std2,temp_std3])
			cf_data[k][i] = np.array([cf_data[k][i][0]/cf_data[k][i][-1],cf_data[k][i][1]/cf_data[k][i][-1],cf_data[k][i][2]/cf_data[k][i][-1],cf_data[k][i][3]/cf_data[k][i][-1]])
			
		cf_data[k] = np.array(cf_data[k])
		cf_data_std[k] = np.array(cf_data_std[k])
		#cf_data_std[k] = cf_data_std[k].reshape(len(cf_data_std[k]),1)
		output[k] = np.hstack((cf_data[k],cf_data_std[k]))


	with open(foutname,'w') as f:
		f.write('Time'.ljust(10)+'          std           '.join([item for item in mode]) + '          std\n')
		for i in range(len(output[0])):
			for k in range(len(mode)):
				item = output[k][i]
				f.write('{:<20e}{:<20e}{:<20e}{:<20e}{:<20e}{:<20e}{:<20e}'.format(item[0],item[1],item[2],item[3],item[4],item[5],item[6]))
			f.write('\n')

	sys.stdout.write('Finished.\n')
	sys.stdout.flush()

