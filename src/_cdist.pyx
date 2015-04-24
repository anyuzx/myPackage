import numpy as np
import cython
cimport numpy as np
from libc.math cimport sqrt,fabs

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def cdist(np.ndarray[DTYPE_t,ndim=2] A,np.ndarray[DTYPE_t,ndim=2] B):
	cdef int N1 = A.shape[0]
	cdef int N2 = B.shape[0]
	cdef int dim = A.shape[1]

	cdef np.ndarray[DTYPE_t,ndim=2] D = np.zeros((N1,N2),dtype=DTYPE)
	cdef DTYPE_t tmp,d
	cdef int i,j,k
	for i in xrange(N1):
		for j in xrange(N2):
			d = 0.0
			for k in range(dim):
				tmp = A[i,k] - B[j,k]
				d += tmp * tmp
			D[j,i] = sqrt(d)

	return D

@cython.boundscheck(False)
@cython.wraparound(False)
def cdistwithpc(np.ndarray[DTYPE_t,ndim=2] A,np.ndarray[DTYPE_t,ndim=2] B,np.ndarray[DTYPE_t,ndim=1] BOX):
	cdef int N1 = A.shape[0]
	cdef int N2 = B.shape[0]
	cdef int dim = A.shape[1]

	cdef np.ndarray[DTYPE_t,ndim=2] D = np.zeros((N1,N2),dtype=DTYPE)
	cdef DTYPE_t tmp,d
	cdef int i,j,k
	for i in xrange(N1):
		for j in xrange(N2):
			d = 0.0
			for k in range(dim):
				tmp = A[i,k] - B[j,k]
				if fabs(tmp) >= BOX[k]*0.5: tmp = -tmp*(BOX[k]/fabs(tmp)-1)
				d += tmp * tmp
			D[j,i] = sqrt(d)

	return D