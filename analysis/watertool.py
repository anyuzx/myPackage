import _hbond
import numpy as np
import sys
from ..core import dist

__all__ = ['calhbond','tetrahedral_order']

def _convert_float32(x):
	# convert argument to numpy.float32 type
	if type(x) == np.ndarray:
		if x.dtype != np.float32:
			try:
				x = x.astype(np.float32,copy=False)
			except:
				raise
	else:
		try:
			x = np.array(x,dtype=np.float32)
		except:
			raise 
	return x

# ============================================================
# ------------------------------------------------------------
# calculate the number of hydrogen bonds based on the following criterion
# distance between two oxygen is less than 3.5 Angstrom. And the angle
# between oxygen-hydrogen bond and oxygen-oxygen bond is less than 30 degrees
# ------------------------------------------------------------

# ------------------------------------------------------------
# Usage:
#		mypackage.analysis.calhbond(object,box,rc,angle)
#		object: water molecules array.(for now it's only for SPC/E water)
#				this array needs to be N*3*3, where N is the number of water
#				molecules and 3 means 3 atoms in one molecules and coordinate 
#				for x,y,z dimension
#		box: box size array [box_size_x,box_size_y,box_size_z]
#		rc: the distance criterion for counting nearest neighbors and hydrogen bonds
#		angle: the angle criterion for counting hydrogen bonds
#       
#		return 1) a array of number of hydrogen bonds for each molecule 
#			   2) a array of number of nearest neighbors for each molecule
# ------------------------------------------------------------

def calhbond(X,BOX,RC,ANGLE):
	x = _convert_float32(X)
	box = _convert_float32(BOX)

	try:
		rc = np.float_(RC)
	except:
		raise TypeError('distance criterion need to be a real number')

	try:
		cos_angle_c = np.cos(ANGLE)
	except:
		raise TypeError('Angle criterion need to be a real number')
	cos_angle_c = np.float_(cos_angle_c)

	if x.shape[1] != 3 or x.shape[2] != 3:
		sys.stdout.write('Wrong molecules array. Should be N*3*3 array')
	if box.shape != (3,):
		sys.stdout.write('Wrong box size array dimension')
		
	return _hbond.hbond(x,box,rc,cos_angle_c)


# --------------------------------------------------------------
def my_cosine(a,b,c,box_size):
	# calculate the cosine of angle between a-b-c
	# a,b,c are coordiantes

	vector1 = np.zeros(3)
	vector2 = np.zeros(3)
	for i in range(3):
		tmp1 = a[i] - b[i]
		tmp2 = c[i] - b[i]
		if np.fabs(tmp1) >= box_size[i]*0.5: tmp1 = -tmp1*(box_size[i]/np.fabs(tmp1)-1)
		if np.fabs(tmp2) >= box_size[i]*0.5: tmp2 = -tmp2*(box_size[i]/np.fabs(tmp2)-1)
		vector1[i] = tmp1
		vector2[i] = tmp2

	return np.sum((vector1*vector2))/(np.linalg.norm(vector1)*np.linalg.norm(vector2))

def cal_order_param(x,box_size):
	# x is the array of coordinates of center atom
	# y is the array of the rest four entries coordinates of the four nearest neighbors

	tmp = 0.0
	for j in range(1,4):
		for k in range(j+1,5):
			tmp += (my_cosine(x[j],x[0],x[k],box_size)+(1.0/3.0))**2

	q = 1.0 - (3.0/8.0)*tmp
	return q

def tetrahedral_order(atoms,box_size):
	D = dist.pdist(atoms,pc=box_size)
	D = D + D.T
	N = atoms.shape[0]

	q = []
	for i in range(N):
		q.append(cal_order_param(atoms[D[:,i].argsort()[:5]],box_size))

	q = np.array(q)
	return q