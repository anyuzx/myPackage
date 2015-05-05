import numpy as np
import cython
cimport numpy as np
from libc.math cimport sqrt,fabs


DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def pdist(np.ndarray[DTYPE_t,ndim=2] X):
	cdef int N = X.shape[0]
	cdef int dim = X.shape[1]

	cdef np.ndarray[DTYPE_t,ndim=2] D = np.zeros((N,N),dtype=DTYPE)
	cdef DTYPE_t tmp,d
	cdef int i,j,k
	for i in xrange(N-1):
		for j in xrange(i+1,N):
			d = 0.0
			for k in range(dim):
				tmp = X[i,k] - X[j,k]
				d += tmp * tmp
			D[j,i] = sqrt(d)

	return D

@cython.boundscheck(False)
@cython.wraparound(False)
def pdistwithpc(np.ndarray[DTYPE_t,ndim=2] X,np.ndarray[DTYPE_t,ndim=1] BOX):
	cdef int N = X.shape[0]
	cdef int dim = X.shape[1]

	cdef np.ndarray[DTYPE_t,ndim=2] D = np.zeros((N,N),dtype=DTYPE)
	cdef DTYPE_t tmp,d
	cdef int i,j,k
	for i in xrange(N-1):
		for j in xrange(i+1,N):
			d = 0.0
			for k in range(dim):
				tmp = X[i,k] - X[j,k]
				if fabs(tmp) >= BOX[k]*0.5: tmp = -tmp*(BOX[k]/fabs(tmp)-1)
				d += tmp * tmp
			D[j,i] = sqrt(d)

	return D