import numpy as np
import cython
cimport numpy as np

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

#######################################################################################
# @cython.boundscheck(False)                                                          #
# @cython.wraparound(False)                                                           #
# def normmap(np.ndarray[DTYPE_t,ndim=2] cmap,np.ndarray[DTYPE_t,ndim=1] p):          #
#       cdef int i                                                                    #
#       cdef int j                                                                    #
#       cdef int N = cmap.shape[0]                                                    #
#                                                                                     #
#       cdef np.ndarray[DTYPE_t,ndim=2] output = np.zeros((N,N),dtype=DTYPE)          #
#       for i in xrange(N):                                                           #
#               for j in xrange(i+1):                                                 #
#                       output[i,j] = cmap[i,j]/p[(i-j)]                              #
#                                                                                     #
#       for i in xrange(N-1):                                                         #
#               for j in xrange(i+1,N):                                               #
#                       output[i,j] = output[j,i]                                     #
#                                                                                     #
#       return output                                                                 #
#######################################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
def normmatrix(np.ndarray[DTYPE_t,ndim=2] cmap, int a):
        cdef int i,j
        cdef int N = cmap.shape[0]
        # a is normalization factor
        cdef int b = N/a

        cdef np.ndarray[DTYPE_t,ndim=2] output = np.zeros((b,b),dtype=DTYPE)
        cdef DTYPE_t tmp
        for i in xrange(b):
                for j in xrange(i):
                        output[i,j] = np.sum(cmap[i*a:(i+1)*a,j*a:(j+1)*a])

        return output

@cython.boundscheck(False)
@cython.wraparound(False)
def normmatrix_OE(np.ndarray[DTYPE_t,ndim=2] cmap, np.ndarray[DTYPE_t,ndim=1] p,int a):
        cdef int i,j
        cdef int N = cmap.shape[0]
        # a is normalization factor
        cdef int b = N/a

        cdef np.ndarray[DTYPE_t,ndim=2] output = np.zeros((b,b),dtype=DTYPE)
        cdef DTYPE_t exptected, observed
        for i in xrange(b):
                for j in xrange(i):
                        expected = np.sum(np.linspace(1,a,a)*p[(np.abs(j-i)*a - a):(np.abs(j-i)*a)] * np.linspace(N - np.abs(j-i)*a + a - 1,N - np.abs(j-i)*a - 1, a+1))
                        observed = np.sum(cmap[i*a:(i+1)*a,j*a:(j+1)*a])
                        output[i,j] = observed/expected

        ####################################################
        # for i in xrange(b):                              #
        #       for j in xrange(i+1,b):                    #
        #               output[i,j] = output[j,i]          #
        ####################################################

        return output






