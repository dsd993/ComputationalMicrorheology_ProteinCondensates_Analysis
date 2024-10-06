import hoomd
hoomd.context.initialize()
from hoomd import md
from hoomd import azplugins
import numpy as np

#Probe information
probe_tag1 = int(75000 + ((736-6)/2))
probe_tag2 = int(75736 + ((736-6)/2))
probe_tag3 = int(76472 + ((736-6)/2))
probe_zloc1 = 12
probe_zloc2 = -60
probe_zloc3 = -120

## Simulation informaiton
T = 300.0 #Kelvin
kB = 0.0019872067 #Kcal/mol/Kelvin
kT = kB*T #Kcal/mol
dt = 0.2045814 #10fs
run_eq = 20000000 #200 nanoseconds
veryshort_period = 1000
short_period = 10000 
long_period = 400000 
veryshortrun_nve = 100000 
shortrun_nve = 1000000 
longrun_nve = 40000000 #0.4 microseconds

## Initialization
system = hoomd.init.read_gsd('initial_configuration.gsd', frame=-1)
sys_snapshot = system.take_snapshot(all = True)
lbox = sys_snapshot.box.Lx #Assuming the simulation box is cubic

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
probe_center1 = hoomd.group.tags(tag_min = probe_tag1,  update = True)
probe_center2 = hoomd.group.tags(tag_min = probe_tag2,  update = True)
probe_center3 = hoomd.group.tags(tag_min = probe_tag3,  update = True)

#Zeroing the linear momentum to avoid any possible COM drift
fix_mum=hoomd.md.update.zero_momentum(period=1)

## Running 1 microsecond long LD simulation
hoomd.md.integrate.mode_standard(dt=dt) #time step is 10fs
integrator = hoomd.md.integrate.langevin(group = solute, kT = kT, seed = np.random.randint(1,1000000)) #Temp is kT/0.0019872067
integrator.set_gamma('1',gamma=0.006312)
integrator.set_gamma('2',gamma=0.006268)
integrator.set_gamma('3',gamma=0.004890)

## push the probe to desired location
pushhp1 = azplugins.restrain.plane(group = probe_center1, point = (0,0,probe_zloc1), normal = (0,0,1), k = 20)
pushhp2 = azplugins.restrain.plane(group = probe_center2, point = (0,0,probe_zloc2), normal = (0,0,1), k = 20)
pushhp3 = azplugins.restrain.plane(group = probe_center3, point = (0,0,probe_zloc3), normal = (0,0,1), k = 20)
pushthermo = hoomd.analyze.log(filename='pushz_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = run_eq/2000, overwrite = True, header_prefix = '#')
hoomd.run(run_eq/200)
pushhp1.disable()
pushhp2.disable()
pushhp3.disable()
pushthermo.disable()

## writing outputs and run the simulation
hp1 = azplugins.restrain.plane(group = probe_center1, point = (0,0,probe_zloc1), normal = (0,0,1), k = 20)
hp2 = azplugins.restrain.plane(group = probe_center2, point = (0,0,probe_zloc2), normal = (0,0,1), k = 20)
hp3 = azplugins.restrain.plane(group = probe_center3, point = (0,0,probe_zloc3), normal = (0,0,1), k = 20)

longthermo = hoomd.analyze.log("long_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = long_period*5, overwrite = True, header_prefix = '#')
longgsd_file = hoomd.dump.gsd("long_dumps.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
longgsd_file_probe = hoomd.dump.gsd("long_dumps_probe.gsd", period = long_period/100, dynamic = ['property','momentum'], group = probe, overwrite=True)
writegsd = hoomd.dump.gsd("long_restart.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True, truncate=True)
hoomd.run(tsteps=longrun_nve)
longthermo.disable()
longgsd_file.disable()
longgsd_file_probe.disable()

shortthermo = hoomd.analyze.log("short_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
shortgsd_file_probe = hoomd.dump.gsd("short_dumps_probe.gsd", period = short_period/100, dynamic = ['property','momentum'], group = probe, overwrite=True)
hoomd.run(tsteps=shortrun_nve)
shortthermo.disable()
shortgsd_file_probe.disable()

veryshortthermo = hoomd.analyze.log("veryshort_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = veryshort_period*10, overwrite = True, header_prefix = '#')
veryshortgsd_file_probe = hoomd.dump.gsd("veryshort_dumps_probe.gsd", period = veryshort_period/100, dynamic = ['property','momentum'], group = probe, overwrite=True)
hoomd.run(tsteps=veryshortrun_nve)
veryshortthermo.disable()
veryshortgsd_file_probe.disable()