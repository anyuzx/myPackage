import numpy as np
import cython
cimport numpy as np
from libc.math cimport sqrt,fabs

DTYPE1 = np.float32
DTYPE2 = np.int32
ctypedef np.float32_t DTYPE1_t
ctypedef np.int32_t DTYPE2_t


@cython.boundscheck(False)
@cython.wraparound(False)
def hbond(np.ndarray[DTYPE1_t,ndim=3] X,np.ndarray[DTYPE1_t,ndim=1] BOX, double r_c, double cos_angle_c):
	# the function takes a numpy array
	# stores coordinates of all atoms
	# the array is 3 dimensional
	# number of molecules * number of atom in one molecules * 3
	# the order of atoms need to be oxygen,hydrogen,hydrogen
	#
	# r_c: the distance criterion for counting neareast neighbors and hydrogen bonds
	# cos_angle_c: the cosine angle criterion for counting hydrogen bonds

	cdef int N = X.shape[0]

	# hydrogen bond list for each water molecules
	cdef np.ndarray[DTYPE2_t,ndim=1] hbond_lst = np.zeros(N,dtype=DTYPE2)
	cdef np.ndarray[DTYPE2_t,ndim=1] neighbors = np.zeros(N,dtype=DTYPE2)
	cdef np.ndarray[DTYPE1_t,ndim=1] roo_vector = np.zeros(3,dtype=DTYPE1)
	cdef np.ndarray[DTYPE1_t,ndim=1] roh_vector = np.zeros(3,dtype=DTYPE1)

	cdef DTYPE1_t tmp,d,cos_angle,roo,dx,dy,dz
	cdef int i,j,k,p
	for i in xrange(N-1):
		for j in xrange(i+1,N):
			d = 0.0
			for k in range(3):
				tmp = X[i,0,k] - X[j,0,k]
				if fabs(tmp) >= BOX[k]*0.5: tmp = -tmp*(BOX[k]/fabs(tmp)-1)
				roo_vector[k] = tmp
				d += tmp * tmp
			roo = sqrt(d)
			if roo < rc:
				# count the nearest neighbors for hydrogen bond efficiency.
				neighbors[i] += 1
				neighbors[j] += 1

				# evaluate the hydrogens for one oxygen
				for p in xrange(1,3):
					dx = X[j,p,0] - X[j,0,0]
					dy = X[j,p,1] - X[j,0,1]
					dz = X[j,p,2] - X[j,0,2]
					if fabs(dx) >= BOX[0]*0.5: dx = -dx*(BOX[0]/fabs(dx)-1)
					if fabs(dy) >= BOX[1]*0.5: dy = -dy*(BOX[1]/fabs(dy)-1)
					if fabs(dz) >= BOX[2]*0.5: dz = -dz*(BOX[2]/fabs(dz)-1)
					
					roh_vector[0] = dx
					roh_vector[1] = dy
					roh_vector[2] = dz
					#roh_vector = vector(X[j,0,:],X[j,p,:],BOX)
					#roh_vector = np.array([X[j,p,0] - X[j,0,0], X[j,p,1] - X[j,0,1], X[j,p,2] - X[j,0,2]],dtype=DTYPE1_t)
					cos_angle = (roh_vector[0]*roo_vector[0] + roh_vector[1]*roo_vector[1] + roh_vector[2]*roo_vector[2])/roo
					if cos_angle >= cos_angle_c:
						hbond_lst[i] += 1
						hbond_lst[j] += 1
						break
				# evaluate the hydrogens for the other oxygen
				for p in xrange(1,3):
					dx = X[i,p,0] - X[i,0,0]
					dy = X[i,p,1] - X[i,0,1]
					dz = X[i,p,2] - X[i,0,2]
					if fabs(dx) >= BOX[0]*0.5: dx = -dx*(BOX[0]/fabs(dx)-1)
					if fabs(dy) >= BOX[1]*0.5: dy = -dy*(BOX[1]/fabs(dy)-1)
					if fabs(dz) >= BOX[2]*0.5: dz = -dz*(BOX[2]/fabs(dz)-1)
					
					roh_vector[0] = dx
					roh_vector[1] = dy
					roh_vector[2] = dz
					#roh_vector = np.array([dx,dy,dz],dtype=DTYPE1_t)
					#roh_vector = vector(X[i,0,:],X[i,p,:],BOX)
					#roh_vector = np.array([X[i,p,0] - X[i,0,0], X[i,p,1] - X[i,0,1], X[i,p,2] - X[i,0,2]],dtype=DTYPE1_t)
					cos_angle = (-roh_vector[0]*roo_vector[0] - roh_vector[1]*roo_vector[1] - roh_vector[2]*roo_vector[2])/roo
					if cos_angle >= cos_angle_c:
						hbond_lst[i] += 1
						hbond_lst[j] += 1
						break

	return hbond_lst,neighbors



