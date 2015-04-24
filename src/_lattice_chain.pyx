import cython
import numpy as np
cimport numpy as np

cdef extern from "c_lattice_chain.h":
	void c_lattice_chain(double* chain, int N, double l0, int ve, int t)

@cython.boundscheck(False)
@cython.wraparound(False)
def lattice_chain(int N, double l0, int ve, int t):
	cdef np.ndarray[double,ndim = 1,mode="c"] chain = np.zeros(N*3)
	c_lattice_chain(&chain[0],N,l0,ve,t)
	return chain