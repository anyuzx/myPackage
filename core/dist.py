import _pdist # _pdist module is the c module wrapped by cython.
import _cdist # _cdist module is the c module wrapped by cython. 
import numpy as np
import sys

# ------------------------------------------------------------
# pdist: calculate the pair-wise distance given one list of coordinates 
# cdist: calculate the pair-wise distance given two list of coordinates
# ------------------------------------------------------------

# ------------------------------------------------------------
# Usage:
# 		mypackage.pdist(object,PC=False)
# 		object: a numpy array or python list
#       PC: periodic boundary. default value False(No periodic boundary condition)
#
# Example: 
#		import numpy as np
#		import mypackage as mp
#		x = np.random.rand(1000,3,dtype=np.float32)
# 		mp.pdist(x)
#		mp.pdist(x,PC=[1.0,1.0,1.0])
# ============================================================
# Usage:
#		mypackage.cdist(object,PC=False)
#		object: a numpy array or python list
#		PC: periodic boundary. default value False(No periodic boundary condition)
#
# Example:
#		import numpy as np
#		import mypackage as mp
#		a = np.random.rand(1000,3,dtype=np.float32)
#		b = np.random.rand(1000,3.dtype=np.float32)
#		mp.cdist(a,b)
#		mp.cdist(a,b,PC=[2.0,1.0,3.0])

# _pdist and _cdist modules are doing the actual computation work
# the reason I write these modules myself is that on some computers 
# or clusters, the pdist and cdist function of scipy package is 
# VERY slow for some reasons I don't know so I decided to write 
# my own pdist function

__all__ = ['pdist','cdist']

def _convert_float32(x):
	# convert argument to numpy.float32 type
	if type(x) == np.ndarray:
		if x.dtype != np.float32:
			try:
				x = x.astype(np.float32,copy=False)
			except:
				sys.stdout.write('ERROR:Cannot convert to np.float32 type')
	else:
		try:
			x = np.array(x,dtype=np.float32)
		except:
			sys.stdout.write('ERROR:Cannot convert to np.float32 type')
	return x

def pdist(X,pc=False):
	x = _convert_float32(X)

	if type(pc) == bool:
		if pc:
			raise ValueError('ERROR: pc can either be False or a box size list')
		cond = False
	else:
		PC = _convert_float32(pc)
		cond = True

	return _pdist.pdist(x) if not cond else _pdist.pdistwithpc(x,PC)

def cdist(A,B,pc=False):
	a = _convert_float32(A)
	b = _convert_float32(B)
	if a.shape[1] != b.shape[1]:
		raise ValueError('ERROR: Dimensions of two lists must be the same')

	if type(pc) == bool:
		if pc:
			raise ValueError('ERROR: pc can either be False or a box size list')
		cond = False
	else:
		PC = _convert_float32(pc)
		cond = True

	return _cdist.cdist(a,b) if not cond else _cdist.cdistwithpc(a,b,PC)


