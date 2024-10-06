# Converts the time frames in GSD to dump file format, printing the center of mass coordinates of the probe
# Coordinates are unwrapped and then printed to dump files

import numpy as np
import gsd.hoomd
import sys
import os

N_chains = 300
N_protein = 250 #Beads per chain
N_probe = 736
N_atoms = N_chains*N_protein

os.makedirs('dumpTraj_long_probe', exist_ok = True)

os.chdir('dumpTraj_long_probe')

traj = gsd.hoomd.open("../long_dumps_probe.gsd")
timesteps = np.array([s.configuration.step for s in traj], dtype=np.int64)
start_tstep = traj[0].configuration.step
interval_tstep = traj[1].configuration.step - traj[0].configuration.step

for frames in np.arange(0,len(traj),1):
    molid = 1
    print(int(frames))
    s = traj[int(frames)]
    f= open('dump.traj.'+str(int(start_tstep + frames*interval_tstep)),'w')

    f.write('ITEM: TIMESTEP\n')
    f.write('%d \n' % s.configuration.step)
    f.write('ITEM: NUMBER OF ATOMS \n')
    f.write('%d \n' % (s.particles.N/N_probe))
    f.write('ITEM: BOX BOUNDS pp pp pp \n')
    f.write('%f %f \n' % (-s.configuration.box[0]/2, s.configuration.box[0]/2))
    f.write('%f %f \n' % (-s.configuration.box[1]/2, s.configuration.box[1]/2))
    f.write('%f %f \n' % (-s.configuration.box[2]/2, s.configuration.box[2]/2))
    f.write('ITEM: ATOMS id mol type q mass xu yu zu ix iy iz \n')

    add_xu = 0
    add_yu = 0
    add_zu = 0
    add_m = 0
    for i in np.arange(0,s.particles.N):
        xu = s.particles.position[i,:][0] + (s.configuration.box[0]*s.particles.image[i,:][0])
        yu = s.particles.position[i,:][1] + (s.configuration.box[1]*s.particles.image[i,:][1])
        zu = s.particles.position[i,:][2] + (s.configuration.box[2]*s.particles.image[i,:][2])

        add_xu = add_xu + xu*s.particles.mass[i]
        add_yu = add_yu + yu*s.particles.mass[i]
        add_zu = add_zu + zu*s.particles.mass[i]
        add_m = add_m + s.particles.mass[i]
 
        if ((i+1)%N_probe == 0):
            final_xu = add_xu/add_m
            final_yu = add_yu/add_m
            final_zu = add_zu/add_m

            f.write('%d %d %d %.8f %.8f %.8f %.8f %.8f %d %d %d \n' % (molid, molid, 1, 0, add_m, final_xu, final_yu, final_zu, 0, 0, 0))
            molid = molid+1

            add_xu = 0
            add_yu = 0
            add_zu = 0
            add_m = 0

    f.close()