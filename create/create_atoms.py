import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

# ------------------------------------------------------------
# create random particles on lattice
# ------------------------------------------------------------

# ------------------------------------------------------------
# Usage:
#       myPackage.create.create_atoms.create_atoms(atom_array, box_dimension, \
#                                                  lattice_constant, \
#                                                  initial array)
#       atom_array: int array
#          array of all types of particles
#          Ex. [1000,500] means create 1000 atoms of type 1 and 500 type 2 atoms
#       box_dimension: float array
#          the dimension of box you want put the particles in
#       lattice_constant: float
#          lattice constant. the size of lattice
#       initial array: float array
#          if there is already particles in the system, specify it
#
# Example:
#       import numpy as np
#       import myPackage as mp
#       mp.create.create_atoms.create_atoms([10,10],\
#                                            box_dimension = [[0,100],[0,100],\
#                                                            [0,100]],\
#                                            lattice_constant = 1.13)
#       mp.create.create_atoms.create_atoms([10,10],\
#                                            box_dimension = [[0,100],[0,100],\
#                                                            [0,100]],\
#                                            lattice_constant = 1.13,\
#                                            initial = [[1,1,1],[2,2,2],[3,3,3]])
#
# second command means that there are already three particles whose position are
# [1,1,1],[2,2,2] and [3,3,3] in the box. so the code will make sure that the
# particles generated will not overlap with these existing particles.
# ------------------------------------------------------------

# this function gives non-duplicated numpy array from a duplicated array
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# this function generate random position on lattice
def make_point(box_dimension, lattice_constant):
    nx = (box_dimension[0,1] - box_dimension[0,0])/lattice_constant
    ny = (box_dimension[1,1] - box_dimension[1,0])/lattice_constant
    nz = (box_dimension[2,1] - box_dimension[2,0])/lattice_constant
    x = np.random.randint(1, high = nx + 1)
    y = np.random.randint(1, high = ny + 1)
    z = np.random.randint(1, high = nz + 1)
    x = box_dimension[0,0] + (x-1)*lattice_constant
    y = box_dimension[1,0] + (y-1)*lattice_constant
    z = box_dimension[2,0] + (z-1)*lattice_constant
    return np.array([x,y,z])

# this function check whether any two points are at the same location
def check_duplicate(point_coord, points_array):
    if len(points_array) != 0:
        old = np.vstack((points_array, point_coord))
        new = unique_rows(old)
        if len(new) != len(old):
            return 'duplicate'
        else:
            return 'no duplicate'
    else:
        return 'no duplicate'

# create atoms on lattice randomly.
def create_atoms(atom_array, box_dimension, lattice_constant, initial = np.array([])):
    box_dimension = np.array(box_dimension)
    atom_array = np.array(atom_array)
    bind_atom_config = []

    if len(initial) != 0:
        config = np.array(initial)
    else:
        config = np.array([])

    type_index = 1
    atom_type = np.array([])
    for n in atom_array:
        for i in range(n):
            point = make_point(box_dimension, lattice_constant)
            while check_duplicate(point, config) == 'duplicate':
                point = make_point(box_dimension, lattice_constant)
            if len(config) == 0:
                config = np.hstack((config, point))
                bind_atom_config.append(point)
            else:
                config = np.vstack((config, point))
                bind_atom_config.append(point)

        atom_type_temp = np.full(n,type_index,dtype=np.int)
        atom_type = np.concatenate((atom_type, atom_type_temp))
        type_index += 1

    bind_atom_config = np.array(bind_atom_config)
    bind_atom_config = np.hstack((atom_type.reshape(len(atom_type),1),bind_atom_config))
    bind_atom_config = pd.DataFrame({'atom_type':bind_atom_config[:,0],
                           'atom_index':np.arange(1,len(bind_atom_config)+1),
                           'x':bind_atom_config[:,1],
                           'y':bind_atom_config[:,2],
                           'z':bind_atom_config[:,3]})
    return bind_atom_config
