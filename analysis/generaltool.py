import numpy as np

__all__ = ['rdf']

def rdf(dist,n_bins,r_cut,N_a,rou_b):
	hist,bins_edge = np.histogram(dist,bins = np.linspace(0,r_cut,n_bins+1))
	g = hist/(4*np.pi*rou_b*(bins_edge[1:]**2)*(r_cut/n_bins)*N_a)
	return g,bins_edge[1:]

def msd(snap_init,snap_t):
	return np.sum(np.mean(np.power(snap_t-snap_init,2),axis=0))

def vacf(snap_init,snap_t):
	return np.mean(np.sum(snap_t*snap_init,axis=1))
