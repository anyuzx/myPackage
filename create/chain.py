import _lattice_chain # _lattice_SAW is the c module wrapped by cython
import numpy as np
import sys

# ------------------------------------------------------------
# give a random chain on 3D lattice
# chain can be ideal or self-avoiding
# ------------------------------------------------------------

# ------------------------------------------------------------
# Usage:
#       myPackage.create.chain(N,l=1.0,VolumeExclusion=True,NumberIteration=None)
#       N: int
#          number of steps(sites) of chain
#       l: float or int
#          the length of one step or distance between two neighbor sites
#       VolumeExclusion: bool
#          default value: True.
#       NumberIteration: int
#          default value: None. But the function will use the same number of steps
#                         if NumberIteration is not set
#                         This means how many successful pivot moves.
#
# Example: 
#       import numpy as np
#       import myPackage as mp
#       mp.create.chain(1000,l=1.5,VolumeExclusion=True,NumberIteration=5000)
# ------------------------------------------------------------

__all__ = ['chain']

class chain:

    # ---------------------------------------
    def __init__(self,N,**kwargs):
        if not isinstance(N,int):
            raise Exception('ERROR: Number of sites must be integer\n')
        else:
            self.N = N
        
        self.l = kwargs.get('l',1.0)
        try:
            self.l = float(self.l)
        except:
            sys.stdout.write('ERROR: the length of one step must be a real number\n')

        self.VolumeExclusion = kwargs.get('VolumeExclusion',True)
        if not isinstance(self.VolumeExclusion,bool):
            raise Exception('ERROR: VolumeExclusion must be boolean type\n')

        self.NumberIteration = kwargs.get('NumberIteration',None)
        if self.NumberIteration == None:
            self.NumberIteration = self.N
        else:
            if not isinstance(self.NumberIteration,int):
                raise Exception('ERROR: Number of iterations must be integer\n')

        self.chain_config = []

    def create_chain(self):
        self.chain_config = _lattice_chain.lattice_chain(self.N,self.l,int(self.VolumeExclusion==True),self.NumberIteration).reshape(self.N,3)
        box_xmax = np.max(self.chain_config[:,0]) + 10.0
        box_xmin = np.min(self.chain_config[:,0]) - 10.0
        box_ymax = np.max(self.chain_config[:,1]) + 10.0
        box_ymin = np.min(self.chain_config[:,1]) - 10.0
        box_zmax = np.max(self.chain_config[:,2]) + 10.0
        box_zmin = np.min(self.chain_config[:,2]) - 10.0
        self.box = np.array([[box_xmin,box_xmax],[box_ymin,box_ymax],[box_zmin,box_zmax]])
        return self.chain_config, self.box

    def plot(self,bond=True):
        try:
            self.chain_config
        except NameError:
            sys.stdout.write('ERROR: No existing chain found. Please create the chain first.\n')

        if len(self.chain_config) == 0:
            raise ValueError('ERROR: No existing chain found. Please create the chain first.\n')

        # import matplotlib package
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        if bond:
            plt.plot(self.chain_config[:,0],self.chain_config[:,1],self.chain_config[:,2],'o--')
        else:
            plt.scatter(self.chain_config[:,0],self.chain_config[:,1],self.chain_config[:,2])
        
        ax.set_aspect('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        xlim1 = np.min(self.chain_config[:,0].flatten())
        xlim2 = np.max(self.chain_config[:,0].flatten())
        ylim1 = np.min(self.chain_config[:,1].flatten())
        ylim2 = np.max(self.chain_config[:,1].flatten())
        zlim1 = np.min(self.chain_config[:,2].flatten())
        zlim2 = np.max(self.chain_config[:,2].flatten())
        ax.auto_scale_xyz([xlim1,xlim2],[ylim1,ylim2],[zlim1,zlim2])
        plt.show()

