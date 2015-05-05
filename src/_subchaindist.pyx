import numpy as np
import cython
cimport numpy as np
from libc.math cimport sqrt

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def sbd(np.ndarray[DTYPE_t,ndim=2] D):
	cdef int N = D.shape[0]

	cdef np.ndarray[DTYPE_t,ndim=1] sbd_lst = np.zeros(N-1,dtype=DTYPE)
	cdef int i,j
	for i in xrange(N-1):
		for j in xrange(i+1,N):
			sbd_lst[j-i-1] += D[j,i]/(N-(j-i))

	return sbd_lst