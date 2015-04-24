import numpy as np
from itertools import groupby
import sys
from ..core import dist # use my own module to calculate the pair distance with or without periodic condition

__all__ = ['ice1h']

# ---------------------------------------
# define a unit cell of oxygen. We will use this to generate oxygen-oxygen tedrahedral network.
# a good description of ice1h structure can be found here: http://www.uwgb.edu/dutchs/petrology/Ice%20Structure.HTM
# There are 8 oxygen atoms in one unit cell.
def unit_cell(a):
    a = np.float_(a)
    unit_cell = np.array([[0,0,0],[0,0,-a],\
                         [-np.sqrt(6)*a/3,5*np.sqrt(2)*a/3,(1.0/3.0)*a],[0,2*np.sqrt(2)*a/3,(1.0/3.0)*a],\
                         [0,2*np.sqrt(2)*a/3,(1.0/3.0+1)*a],[-np.sqrt(6)*a/3,np.sqrt(2)*a,0],\
                         [-np.sqrt(6)*a/3,5*np.sqrt(2)*a/3,(1.0/3.0+1)*a],[-np.sqrt(6)*a/3,np.sqrt(2)*a,-a]])
    # unit cell dimension
    base_i = np.array([2*np.sqrt(6)*a/3,0,0]) 
    base_j = np.array([0,2*np.sqrt(2)*a,0])
    base_k = np.array([0,0,(8.0/3.0)*a])
    base = np.array([base_i,base_j,base_k])
    return unit_cell,base


# ---------------------------------------
# give the coordinate of hydrogen based on two neighbor oxygen
# there are twoo alernative position for hydrogen on each oxygen-oxygen bond
# this function give one given the index 
# the function itself is periodic boundary implemented.
# however for vapor-interface situation, the coordinates will be corrected later in other part of the code
def givehydrogen(O1,O2,a,b,dim1,dim2,dim3,index):
    a = np.float_(a)
    b = np.float_(b)
    d = O2-O1 # O1,O2 are the positions of oxygen
    dx = d[0]
    dy = d[1]
    dz = d[2]

    # x_size,y_size,z_size are the periodic boundary size for ice1h.
    # dim1,dim2,dim3 are the replicate number on three directions.
    x_size = dim1*(2*np.sqrt(6)*a/3)
    y_size = (dim2)*(2*a*np.sqrt(2))
    z_size = (8.0/3.0)*a*dim3

    if np.abs(dx) > x_size * 0.5: dx = -dx*(x_size/np.abs(dx)-1)
    if np.abs(dy) > y_size * 0.5: dy = -dy*(y_size/np.abs(dy)-1)
    if np.abs(dz) > z_size * 0.5: dz = -dz*(z_size/np.abs(dz)-1)

    list = np.array([[(b/a)*dx+O1[0], (b/a)*dy+O1[1] , (b/a)*dz+O1[2]],\
                    [((a-b)/a)*dx+O1[0] , ((a-b)/a)*dy+O1[1] , ((a-b)/a)*dz+O1[2]]])
    return list[index]

# ---------------------------------------
# this function create oxygen-oxygen tetrahedral network 
# by duplicating the unit cell
# duplicating number is dim1,dim2,dim3 along x,y,z direction
def makeOxygen(a,dim1,dim2,dim3):
    unit,base = unit_cell(a)

    oxygen = []
    for i in range(dim1):
        for j in range(dim2):
            for k in range(dim3):
                oxygen.append(unit+i*base[0]+j*base[1]+k*base[2])

    oxygen = np.array(oxygen)
    oxygen = np.reshape(oxygen,((dim1*dim2*dim3)*8,3))

    return oxygen

# ---------------------------------------
# this function get all the bonds from position of oxygens
# calculate distances between any oxygens. if the distance equal to a(oxygen-oxygen distance in ideal lattice),
# then store that o-o pair as a bond.
# bond array looks like this: [[o1,o2],[o3,o4],[o5,o6],[o7,o8],...]
# this means o1 and o2 are one pair,o3 and o4 are one pair and so on.
def getbond(oxygen,a,dim1,dim2,dim3):
    # distance is calculated using myPackage.basic_tools.pdist with periodic boundary condition
    dist = dist.pdist(oxygen,pc=[dim1*(2*np.sqrt(6)*a/3),(dim2)*(2*a*np.sqrt(2)),(8.0/3.0)*a*dim3])
    temp = np.where(np.abs(dist - a) < 0.01)
    bonds = np.dstack((temp[0],temp[1]))[0]
    return bonds

# ---------------------------------------
# this function put all hydrogen on each o-o bond
# randomly choose one from two possible position for each hydrogen
# the function return:
# hydrogen: array for position of all hydrogen
# chemical_coord: array for the number of hydrogen for each oxygen
# hydrogen_index: array for indicating whether the hydrogen is belong to which oxygen on that bond
# the order of hydrogen and hydrogen_index has the same order as bonds
# the order of chemical_coord has the same order as oxygen
#
# Ex. hydrogen_index = 0, then this means the hydrogen of this bond is connected to the first oxygen of that o-o pair.
#     if hydrogen_index = 1, then the hydrogen of this bond is connected to the second oxygen of that o-o pair.
def makeHydrogen(oxygen,bonds,a,b,dim1,dim2,dim3):
    hydrogen = []
    chemical_coord = np.zeros(len(oxygen))
    hydrogen_index = []
    for bond in bonds:
        opair = oxygen[bond]
        index = np.random.randint(2)
        h = givehydrogen(opair[0],opair[1],a,b,dim1,dim2,dim3,index)
        hydrogen.append(h)
        chemical_coord[bond[index]] += 1
        hydrogen_index.append(index)
    hydrogen = np.array(hydrogen)
    hydrogen_index = np.array(hydrogen_index)
    return hydrogen,chemical_coord,hydrogen_index

class ice1h:

    # ---------------------------------------
    def __init__(self,*argv,**kwargs):
        self.vapor_interface = kwargs.get('vapor_interface',False)
        if len(argv) < 5:
            raise Exception('ERROR: Please specify the arguments needed.\n')
        self.a = argv[0]
        self.b = argv[1]
        self.dim1 = argv[2]
        self.dim2 = argv[3]
        self.dim3 = argv[4]
        self.phase = '1h'
        self.ice1h_config = []
        self.box = []


    # ---------------------------------------
    def create_ice1h(self):
        pcsize = np.array([self.dim1*(2*np.sqrt(6)*self.a/3),self.dim2*(2*self.a*np.sqrt(2)),self.dim3*(8.0/3.0)*self.a])
        box = np.array([[0,self.dim1*(2*np.sqrt(6)*self.a/3)],[0,self.dim2*(2*self.a*np.sqrt(2))],[0,self.dim3*(8.0/3.0)*self.a]])

        oxygen = makeOxygen(self.a,self.dim1,self.dim2,self.dim3)
        bonds = getbond(oxygen,self.a,self.dim1,self.dim2,self.dim3)

        hydrogen, chemical_coord, hydrogen_index = makeHydrogen(oxygen,bonds,self.a,self.b,self.dim1,self.dim2,self.dim3)

        while not np.array_equal(chemical_coord,np.full(len(oxygen),2.0)): 
            chemical_coord_temp = np.copy(chemical_coord)
            pick = np.random.randint(len(bonds))
            while np.array_equal(chemical_coord[bonds[pick]],np.array([2.0,2.0])):
                pick = np.random.randint(len(bonds))
            coord_diff = chemical_coord[bonds[pick]][0]-chemical_coord[bonds[pick]][1]
            if hydrogen_index[pick] == 0:
                if np.abs(coord_diff-2) < np.abs(coord_diff):
                    opair = oxygen[bonds[pick]]
                    hydrogen[pick] = givehydrogen(opair[0],opair[1],self.a,self.b,self.dim1,self.dim2,self.dim3,1)
                    chemical_coord[bonds[pick][0]] -= 1
                    chemical_coord[bonds[pick][1]] += 1
                    hydrogen_index[pick] = 1
                elif np.abs(coord_diff-2) == np.abs(coord_diff):
                    opair = oxygen[bonds[pick]]
                    if np.random.rand() > 0.5:
                        hydrogen[pick] = givehydrogen(opair[0],opair[1],self.a,self.b,self.dim1,self.dim2,self.dim3,1)
                        chemical_coord[bonds[pick][0]] -= 1
                        chemical_coord[bonds[pick][1]] += 1
                        hydrogen_index[pick] = 1
            elif hydrogen_index[pick] == 1:
                if np.abs(coord_diff+2) < np.abs(coord_diff):
                    opair = oxygen[bonds[pick]]
                    hydrogen[pick] = givehydrogen(opair[0],opair[1],self.a,self.b,self.dim1,self.dim2,self.dim3,0)
                    chemical_coord[bonds[pick][0]] += 1
                    chemical_coord[bonds[pick][1]] -= 1
                    hydrogen_index[pick] = 0
                elif np.abs(coord_diff+2) == np.abs(coord_diff):
                    opair = oxygen[bonds[pick]]
                    if np.random.rand() > 0.5:
                        hydrogen[pick] = givehydrogen(opair[0],opair[1],self.a,self.b,self.dim1,self.dim2,self.dim3,0)
                        chemical_coord[bonds[pick][0]] += 1
                        chemical_coord[bonds[pick][1]] -= 1
                        hydrogen_index[pick] = 0

        ice1h_config = []
        for o in oxygen:
            ice1h_config.append(o)

        for i in range(len(bonds)):
            index = bonds[i][hydrogen_index[i]]
            ice1h_config[index] = np.vstack((ice1h_config[index],hydrogen[i]))

        ice1h_config = np.array(ice1h_config)
        if not self.vapor_interface:
            temp = [0,1,2]
            for i in range(len(ice1h_config)):
                molecule = ice1h_config[i]
                for j in range(3):
                    for dimension in temp:
                        if molecule[j][dimension] > box[dimension,1]:
                            ice1h_config[i,j,dimension] = -pcsize[dimension] + molecule[j][dimension]
                        elif molecule[j][dimension] < box[dimension,0]:
                            ice1h_config[i,j,dimension] = pcsize[dimension] + molecule[j][dimension]

            for dimension in temp:
                temp_min = np.min(ice1h_config[:,:,dimension].flatten())
                temp_max = np.max(ice1h_config[:,:,dimension].flatten())
                temp_lim1 = box[dimension,0] - (pcsize[dimension]-(temp_max-temp_min))/2.0
                temp_lim2 = box[dimension,1] - (pcsize[dimension]-(temp_max-temp_min))/2.0
                box[dimension,0] = temp_lim1
                box[dimension,1] = temp_lim2
        else:
            temp = [0,1,2]
            def tempfunc(a):
                if a == 'x':
                    return 0
                elif a == 'y':
                    return 1
                elif a == 'z':
                    return 2

            # correcting the coordinate of hydrogen for vaper interface direction
            for dimension in self.vapor_interface:
                dim = tempfunc(dimension)
                temp.remove(dim)
                index = np.where(np.dstack((np.abs(ice1h_config[:,0][:,dim] - ice1h_config[:,1][:,dim]),np.abs(ice1h_config[:,0][:,dim] - ice1h_config[:,2][:,dim])))[0]>(self.b+0.1))
                index = zip(index[0],index[1])
                for item in index:
                    dd = ice1h_config[item[0],item[1]+1,dim] - ice1h_config[item[0],0,dim]
                    dd = -dd*(pcsize[dim]/np.abs(dd)-1)
                    ice1h_config[item[0],item[1]+1,dim] = ice1h_config[item[0],0,dim] + dd

                temp_min = np.min(ice1h_config[:,:,dim].flatten())
                temp_max = np.max(ice1h_config[:,:,dim].flatten())
                box[dim,0] = temp_min - 10.0
                box[dim,1] = temp_max + 10.0
            
            # correct the coordinate of hydrogen for periodic connected direction
            for i in range(len(ice1h_config)):
                molecule = ice1h_config[i]
                for j in range(3):
                    for dimension in temp:
                        if molecule[j][dimension] > box[dimension,1]:
                            ice1h_config[i,j,dimension] = -pcsize[dimension] + molecule[j][dimension]
                        elif molecule[j][dimension] < box[dimension,0]:
                            ice1h_config[i,j,dimension] = pcsize[dimension] + molecule[j][dimension]

            # correct the box limits for periodic connected direction
            for dimension in temp:
                temp_min = np.min(ice1h_config[:,:,dimension].flatten())
                temp_max = np.max(ice1h_config[:,:,dimension].flatten())
                temp_lim1 = box[dimension,0] - (pcsize[dimension]-(temp_max-temp_min))/2.0
                temp_lim2 = box[dimension,1] - (pcsize[dimension]-(temp_max-temp_min))/2.0
                box[dimension,0] = temp_lim1
                box[dimension,1] = temp_lim2

        self.ice1h_config = ice1h_config
        self.box = box
        return self.ice1h_config,self.box

    def plot(self,bond=False):
        # import matplotlib package
        # if not installed, ice can still be made
        # this is only for visualizing ice in 3D
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        oxygen = self.ice1h_config[:,0,:]
        hydrogen = self.ice1h_config[:,1:,:].reshape((2*len(oxygen),3))
        if bond:        
            temp = ['x','y','z']
            if self.vapor_interface:
                for item in self.vapor_interface:
                    temp.remove(item)

            if len(temp) != 0:
                sys.stdout.write(','.join(temp)+' direction has no vapor interface.')

            for molecule in self.ice1h_config:
                for i in range(1,3):
                    temp = molecule[np.array([0,i]),:]
                    ax.plot(temp[:,0],temp[:,1],temp[:,2],c='g')

        ax.scatter(oxygen[:,0],oxygen[:,1],oxygen[:,2],c='b',s=30)
        ax.scatter(hydrogen[:,0],hydrogen[:,1],hydrogen[:,2],c='r',s=15)
        ax.set_aspect('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        xlim1 = np.min(self.ice1h_config[:,:,0].flatten())
        xlim2 = np.max(self.ice1h_config[:,:,0].flatten())
        ylim1 = np.min(self.ice1h_config[:,:,1].flatten())
        ylim2 = np.max(self.ice1h_config[:,:,1].flatten())
        zlim1 = np.min(self.ice1h_config[:,:,2].flatten())
        zlim2 = np.max(self.ice1h_config[:,:,2].flatten())
        ax.auto_scale_xyz([xlim1,xlim2],[ylim1,ylim2],[zlim1,zlim2])
        plt.show()
