import myPackage as mp

mp.lammps.ProcessDumpCorrelation('IDEAL_harmonic_n300_msd_precise_1_traj','test.dat',['MSD','VACF'],100,select=[0,1001])