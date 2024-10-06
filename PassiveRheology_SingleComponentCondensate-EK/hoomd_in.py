import hoomd
hoomd.context.initialize()
from hoomd import md
from hoomd import azplugins
import numpy as np

## Simulation informaiton
T = 300.0 #Kelvin
kB = 0.0019872067 #Kcal/mol/Kelvin
kT = kB*T #Kcal/mol
dt = 0.2045814 #10fs
run_eq = 1500000 #15 nanoseconds
veryshort_period = 10
veryshort_period = 1000
short_period = 10000 
long_period = 400000 
veryshortrun_nve = 100000 
shortrun_nve = 1000000 
longrun_nve = 40000000 #0.4 microseconds

## Initialization
system = hoomd.init.read_gsd('initial_configuration.gsd')
sys_snapshot = system.take_snapshot(all = True)
lbox = sys_snapshot.box.Lx #Assuming the simulation box is cubic

#print(sys_snapshot.particles.typeid)
#num_particles=len(system.particles)
#for i in range(num_particles):
#   if (sys_snapshot.particles.typeid[i] == 2):
#       sys_snapshot.particles.mass[i] = 50
#system.restore_snapshot(sys_snapshot)

## Interactions
### Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', 'body'])

### Bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('1', k=20.000000, r0=3.800000)
harmonic.bond_coeff.set('2', k=250.000000, r0=1.500481) ##########

### Pairwise interactions
nb = azplugins.pair.ashbaugh(r_cut = 0.0, nlist = nl)
nb.pair_coeff.set('1', '1', epsilon = 0.200000, sigma = 5.920000, lam = 0.459459, r_cut = 23.680000)
nb.pair_coeff.set('1', '2', epsilon = 0.200000, sigma = 6.140000, lam = 0.486486, r_cut = 24.560000)
nb.pair_coeff.set('2', '2', epsilon = 0.200000, sigma = 6.360000, lam = 0.513514, r_cut = 25.440000)
nb.pair_coeff.set('3', '3', epsilon = 0.200000, sigma = 1.500000, lam = 0.000000, r_cut = 6.000000)
nb.pair_coeff.set('1', '3', epsilon = 0.200000, sigma = 3.710000, lam = 1.000000, r_cut = 14.840000)
nb.pair_coeff.set('2', '3', epsilon = 0.200000, sigma = 3.930000, lam = 1.000000, r_cut = 15.720000)

### Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut = 0.0, nlist = nl)
yukawa.pair_coeff.set(['1', '2', '3'], ['1', '2', '3'], epsilon=0.000000, kappa=0.000000, r_cut=0.000000) ##########
yukawa.pair_coeff.set('1','1', epsilon = 4.150796, kappa = 0.10, r_cut = 35.00)
yukawa.pair_coeff.set('1','2', epsilon = -4.150796, kappa = 0.10, r_cut = 35.00)
yukawa.pair_coeff.set('2','2', epsilon = 4.150796, kappa = 0.10, r_cut = 35.00)

## Integration
#Particles MD Integrator
solute = hoomd.group.all()
probe = hoomd.group.type('3')

#Zeroing the linear momentum to avoid any possible COM drift
fix_mum=hoomd.md.update.zero_momentum(period=1)

## Running 2 microsecond long LD simulation
hoomd.md.integrate.mode_standard(dt=dt) #time step is 10fs
integrator2 = hoomd.md.integrate.langevin(group = solute, kT = kT, seed = np.random.randint(1,1000000)) #Temp is kT/0.0019872067
integrator2.set_gamma('1',gamma=0.006312)
integrator2.set_gamma('2',gamma=0.006268)
integrator2.set_gamma('3',gamma=0.004890)

## writing outputs and run the simulation
eqthermo = hoomd.analyze.log(filename='afteq_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
hoomd.run(run_eq)
eqthermo.disable()

#Very short, short, and long simulation runs
longthermo = hoomd.analyze.log("long_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = long_period*5, overwrite = True, header_prefix = '#')
longgsd_file = hoomd.dump.gsd("long_dumps.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
longgsd_file_probe = hoomd.dump.gsd("long_dumps_probe.gsd", period = long_period/1000, dynamic = ['property','momentum'], group = probe, overwrite=True)
writegsd = hoomd.dump.gsd("long_restart.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True, truncate=True)
hoomd.run(tsteps=longrun_nve)
longthermo.disable()
longgsd_file.disable()

shortthermo = hoomd.analyze.log("short_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
shortgsd_file = hoomd.dump.gsd("short_dumps.gsd", period = short_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
shortgsd_file_probe = hoomd.dump.gsd("short_dumps_probe.gsd", period = short_period/1000, dynamic = ['property','momentum'], group = probe, overwrite=True)
hoomd.run(tsteps=shortrun_nve)
shortthermo.disable()
shortgsd_file.disable()

veryshortthermo = hoomd.analyze.log("veryshort_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = veryshort_period*10, overwrite = True, header_prefix = '#')
veryshortgsd_file = hoomd.dump.gsd("veryshort_dumps.gsd", period = veryshort_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
veryshortgsd_file_probe = hoomd.dump.gsd("veryshort_dumps_probe.gsd", period = veryshort_period/1000, dynamic = ['property','momentum'], group = probe, overwrite=True)
hoomd.run(tsteps=veryshortrun_nve)
veryshortthermo.disable()
veryshortgsd_file.disable()