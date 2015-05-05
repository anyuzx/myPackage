import numpy as np
import cython
cimport numpy as np
from libc.math cimport sqrt

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def cp(np.ndarray[DTYPE_t,ndim=2] D,double rc):
	cdef int N = D.shape[0]

	cdef np.ndarray[DTYPE_t,ndim=1] contact_lst = np.zeros(N-1,dtype=DTYPE)
	cdef int i,j,k
	for i in xrange(N-1):
		for j in xrange(i+1,N):
			if D[j,i] <= rc:
				contact_lst[j-i-1] += 1.0/(N-(j-i))

	return contact_lst
