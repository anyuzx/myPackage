import numpy as np

__all__ = ['rdf','msd','vacf']

def rdf(dist,n_bins,r_cut,N_a,rou_b):
	
	r_cut = np.float_(r_cut)
	rou_b = np.float_(rou_b)
	N_a = np.float_(N_a)
	n_bins = np.float_(n_bins)

	hist,bins_edge = np.histogram(dist,bins = np.linspace(0,r_cut,n_bins+1))
	bins_edge = np.diff(bins_edge)/2.0 + bins_edge[:-1]
	g = hist/(4*np.pi*rou_b*(bins_edge**2)*(r_cut/n_bins)*N_a)
	return g,bins_edge

def msd(snap_init,snap_t):
	return np.sum(np.mean(np.power(snap_t-snap_init,2),axis=0))

def vacf(snap_init,snap_t):
	return np.mean(np.sum(snap_t*snap_init,axis=1))
