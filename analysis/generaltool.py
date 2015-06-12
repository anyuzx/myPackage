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
        center_mass_init = np.mean(snap_init,axis=0)
        center_mass_t = np.mean(snap_t,axis=0)
        snap_t_center = snap_t - center_mass_t
        snap_init_center = snap_init - center_mass_init
	return np.array([np.sum(np.mean(np.power(snap_t-snap_init,2),axis=0)), np.sum(np.mean(np.power(snap_t_center - snap_init_center,2),axis=0)), np.sum(np.power(center_mass_t - center_mass_init,2))])

def vacf(snap_init,snap_t):
	return np.mean(np.sum(snap_t*snap_init,axis=1))
