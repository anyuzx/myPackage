import numpy as np
import cython
cimport numpy as np

DTYPE = np.int32
DTYPE2 = np.float32
ctypedef np.int32_t DTYPE_t
ctypedef np.float32_t DTYPE2_t

def cmap(np.ndarray[DTYPE2_t,ndim=2] D,double rc):
	cdef int N1 = D.shape[0]
	cdef int N2 = D.shape[1]
	cdef int i,j
	cdef np.ndarray[DTYPE_t,ndim=2] cmap = np.zeros((N1,N2),dtype=DTYPE)
	for i in xrange(N1):
		for j in xrange(N2):
			if D[j,i] <= rc:
				cmap[j,i] += 1

	return cmap

def cmaphalf(np.ndarray[DTYPE2_t,ndim=2] D,double rc):
	cdef int N = D.shape[0]
	cdef int i,j
	cdef np.ndarray[DTYPE_t,ndim=2] cmap = np.zeros((N,N),dtype=DTYPE)
	for i in xrange(N):
		for j in xrange(i,N):
			if D[j,i] <= rc:
				cmap[j,i] += 1

	return cmap