import hoomd
hoomd.context.initialize()
from hoomd import md
from hoomd import azplugins
import numpy as np

## Simulation informaiton
T = 245.0 #Kelvin
kB = 0.0019872067 #Kcal/mol/Kelvin
kT = kB*T #Kcal/mol
dt = 0.2045814 #10fs
run_eq = 2500000 #25 nanoseconds
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

## Interactions
### Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', 'body'])

### Bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('1', k=20.000000, r0=3.800000)
harmonic.bond_coeff.set('2', k=250.000000, r0=1.500481)

### Pairwise interactions
nb = azplugins.pair.ashbaugh(r_cut = 0.0, nlist = nl)
nb.pair_coeff.set('1', '1', epsilon = 0.200000, sigma = 4.500000, lam = 0.493530, r_cut = 18.000000)
nb.pair_coeff.set('1', '2', epsilon = 0.200000, sigma = 4.840000, lam = 0.500883, r_cut = 19.360000)
nb.pair_coeff.set('1', '3', epsilon = 0.200000, sigma = 5.340000, lam = 0.545000, r_cut = 21.360000)
nb.pair_coeff.set('1', '4', epsilon = 0.200000, sigma = 4.770000, lam = 0.508236, r_cut = 19.080000)
nb.pair_coeff.set('1', '5', epsilon = 0.200000, sigma = 5.260000, lam = 0.486177, r_cut = 21.040000)
nb.pair_coeff.set('1', '6', epsilon = 0.200000, sigma = 5.530000, lam = 0.486177, r_cut = 22.120000)
nb.pair_coeff.set('1', '7', epsilon = 0.200000, sigma = 5.090000, lam = 0.500883, r_cut = 20.360000)
nb.pair_coeff.set('1', '8', epsilon = 0.200000, sigma = 5.430000, lam = 0.618530, r_cut = 21.720000)
nb.pair_coeff.set('1', '9', epsilon = 0.200000, sigma = 5.040000, lam = 0.353824, r_cut = 20.160000)
nb.pair_coeff.set('1', '10', epsilon = 0.200000, sigma = 5.480000, lam = 0.655295, r_cut = 21.920000)
nb.pair_coeff.set('1', '11', epsilon = 0.200000, sigma = 5.030000, lam = 0.586177, r_cut = 20.120000)
nb.pair_coeff.set('1', '12', epsilon = 0.200000, sigma = 5.430000, lam = 0.397942, r_cut = 21.720000)
nb.pair_coeff.set('2', '2', epsilon = 0.200000, sigma = 5.180000, lam = 0.508236, r_cut = 20.720000)
nb.pair_coeff.set('2', '3', epsilon = 0.200000, sigma = 5.680000, lam = 0.552354, r_cut = 22.720000)
nb.pair_coeff.set('2', '4', epsilon = 0.200000, sigma = 5.110000, lam = 0.515589, r_cut = 20.440000)
nb.pair_coeff.set('2', '5', epsilon = 0.200000, sigma = 5.600000, lam = 0.493530, r_cut = 22.400000)
nb.pair_coeff.set('2', '6', epsilon = 0.200000, sigma = 5.870000, lam = 0.493530, r_cut = 23.480000)
nb.pair_coeff.set('2', '7', epsilon = 0.200000, sigma = 5.430000, lam = 0.508236, r_cut = 21.720000)
nb.pair_coeff.set('2', '8', epsilon = 0.200000, sigma = 5.770000, lam = 0.625883, r_cut = 23.080000)
nb.pair_coeff.set('2', '9', epsilon = 0.200000, sigma = 5.380000, lam = 0.361178, r_cut = 21.520000)
nb.pair_coeff.set('2', '10', epsilon = 0.200000, sigma = 5.820000, lam = 0.662648, r_cut = 23.280000)
nb.pair_coeff.set('2', '11', epsilon = 0.200000, sigma = 5.370000, lam = 0.593530, r_cut = 21.480000)
nb.pair_coeff.set('2', '12', epsilon = 0.200000, sigma = 5.770000, lam = 0.405295, r_cut = 23.080000)
nb.pair_coeff.set('3', '3', epsilon = 0.200000, sigma = 6.180000, lam = 0.596471, r_cut = 24.720000)
nb.pair_coeff.set('3', '4', epsilon = 0.200000, sigma = 5.610000, lam = 0.559707, r_cut = 22.440000)
nb.pair_coeff.set('3', '5', epsilon = 0.200000, sigma = 6.100000, lam = 0.537648, r_cut = 24.400000)
nb.pair_coeff.set('3', '6', epsilon = 0.200000, sigma = 6.370000, lam = 0.537648, r_cut = 25.480000)
nb.pair_coeff.set('3', '7', epsilon = 0.200000, sigma = 5.930000, lam = 0.552354, r_cut = 23.720000)
nb.pair_coeff.set('3', '8', epsilon = 0.200000, sigma = 6.270000, lam = 0.670000, r_cut = 25.080000)
nb.pair_coeff.set('3', '9', epsilon = 0.200000, sigma = 5.880000, lam = 0.405295, r_cut = 23.520000)
nb.pair_coeff.set('3', '10', epsilon = 0.200000, sigma = 6.320000, lam = 0.706765, r_cut = 25.280000)
nb.pair_coeff.set('3', '11', epsilon = 0.200000, sigma = 5.870000, lam = 0.637648, r_cut = 23.480000)
nb.pair_coeff.set('3', '12', epsilon = 0.200000, sigma = 6.270000, lam = 0.449413, r_cut = 25.080000)
nb.pair_coeff.set('4', '4', epsilon = 0.200000, sigma = 5.040000, lam = 0.522942, r_cut = 20.160000)
nb.pair_coeff.set('4', '5', epsilon = 0.200000, sigma = 5.530000, lam = 0.500883, r_cut = 22.120000)
nb.pair_coeff.set('4', '6', epsilon = 0.200000, sigma = 5.800000, lam = 0.500883, r_cut = 23.200000)
nb.pair_coeff.set('4', '7', epsilon = 0.200000, sigma = 5.360000, lam = 0.515589, r_cut = 21.440000)
nb.pair_coeff.set('4', '8', epsilon = 0.200000, sigma = 5.700000, lam = 0.633236, r_cut = 22.800000)
nb.pair_coeff.set('4', '9', epsilon = 0.200000, sigma = 5.310000, lam = 0.368531, r_cut = 21.240000)
nb.pair_coeff.set('4', '10', epsilon = 0.200000, sigma = 5.750000, lam = 0.670000, r_cut = 23.000000)
nb.pair_coeff.set('4', '11', epsilon = 0.200000, sigma = 5.300000, lam = 0.600883, r_cut = 21.200000)
nb.pair_coeff.set('4', '12', epsilon = 0.200000, sigma = 5.700000, lam = 0.412648, r_cut = 22.800000)
nb.pair_coeff.set('5', '5', epsilon = 0.200000, sigma = 6.020000, lam = 0.478824, r_cut = 24.080000)
nb.pair_coeff.set('5', '6', epsilon = 0.200000, sigma = 6.290000, lam = 0.478824, r_cut = 25.160000)
nb.pair_coeff.set('5', '7', epsilon = 0.200000, sigma = 5.850000, lam = 0.493530, r_cut = 23.400000)
nb.pair_coeff.set('5', '8', epsilon = 0.200000, sigma = 6.190000, lam = 0.611177, r_cut = 24.760000)
nb.pair_coeff.set('5', '9', epsilon = 0.200000, sigma = 5.800000, lam = 0.346471, r_cut = 23.200000)
nb.pair_coeff.set('5', '10', epsilon = 0.200000, sigma = 6.240000, lam = 0.647942, r_cut = 24.960000)
nb.pair_coeff.set('5', '11', epsilon = 0.200000, sigma = 5.790000, lam = 0.578824, r_cut = 23.160000)
nb.pair_coeff.set('5', '12', epsilon = 0.200000, sigma = 6.190000, lam = 0.390589, r_cut = 24.760000)
nb.pair_coeff.set('6', '6', epsilon = 0.200000, sigma = 6.560000, lam = 0.478824, r_cut = 26.240000)
nb.pair_coeff.set('6', '7', epsilon = 0.200000, sigma = 6.120000, lam = 0.493530, r_cut = 24.480000)
nb.pair_coeff.set('6', '8', epsilon = 0.200000, sigma = 6.460000, lam = 0.611177, r_cut = 25.840000)
nb.pair_coeff.set('6', '9', epsilon = 0.200000, sigma = 6.070000, lam = 0.346471, r_cut = 24.280000)
nb.pair_coeff.set('6', '10', epsilon = 0.200000, sigma = 6.510000, lam = 0.647942, r_cut = 26.040000)
nb.pair_coeff.set('6', '11', epsilon = 0.200000, sigma = 6.060000, lam = 0.578824, r_cut = 24.240000)
nb.pair_coeff.set('6', '12', epsilon = 0.200000, sigma = 6.460000, lam = 0.390589, r_cut = 25.840000)
nb.pair_coeff.set('7', '7', epsilon = 0.200000, sigma = 5.680000, lam = 0.508236, r_cut = 22.720000)
nb.pair_coeff.set('7', '8', epsilon = 0.200000, sigma = 6.020000, lam = 0.625883, r_cut = 24.080000)
nb.pair_coeff.set('7', '9', epsilon = 0.200000, sigma = 5.630000, lam = 0.361178, r_cut = 22.520000)
nb.pair_coeff.set('7', '10', epsilon = 0.200000, sigma = 6.070000, lam = 0.662648, r_cut = 24.280000)
nb.pair_coeff.set('7', '11', epsilon = 0.200000, sigma = 5.620000, lam = 0.593530, r_cut = 22.480000)
nb.pair_coeff.set('7', '12', epsilon = 0.200000, sigma = 6.020000, lam = 0.405295, r_cut = 24.080000)
nb.pair_coeff.set('8', '8', epsilon = 0.200000, sigma = 6.360000, lam = 0.743530, r_cut = 25.440000)
nb.pair_coeff.set('8', '9', epsilon = 0.200000, sigma = 5.970000, lam = 0.478825, r_cut = 23.880000)
nb.pair_coeff.set('8', '10', epsilon = 0.200000, sigma = 6.410000, lam = 0.780295, r_cut = 25.640000)
nb.pair_coeff.set('8', '11', epsilon = 0.200000, sigma = 5.960000, lam = 0.711177, r_cut = 23.840000)
nb.pair_coeff.set('8', '12', epsilon = 0.200000, sigma = 6.360000, lam = 0.522942, r_cut = 25.440000)
nb.pair_coeff.set('9', '9', epsilon = 0.200000, sigma = 5.580000, lam = 0.214119, r_cut = 22.320000)
nb.pair_coeff.set('9', '10', epsilon = 0.200000, sigma = 6.020000, lam = 0.515589, r_cut = 24.080000)
nb.pair_coeff.set('9', '11', epsilon = 0.200000, sigma = 5.570000, lam = 0.446472, r_cut = 22.280000)
nb.pair_coeff.set('9', '12', epsilon = 0.200000, sigma = 5.970000, lam = 0.258237, r_cut = 23.880000)
nb.pair_coeff.set('10', '10', epsilon = 0.200000, sigma = 6.460000, lam = 0.817059, r_cut = 25.840000)
nb.pair_coeff.set('10', '11', epsilon = 0.200000, sigma = 6.010000, lam = 0.747942, r_cut = 24.040000)
nb.pair_coeff.set('10', '12', epsilon = 0.200000, sigma = 6.410000, lam = 0.559707, r_cut = 25.640000)
nb.pair_coeff.set('11', '11', epsilon = 0.200000, sigma = 5.560000, lam = 0.678824, r_cut = 22.240000)
nb.pair_coeff.set('11', '12', epsilon = 0.200000, sigma = 5.960000, lam = 0.490589, r_cut = 23.840000)
nb.pair_coeff.set('12', '12', epsilon = 0.200000, sigma = 6.360000, lam = 0.302354, r_cut = 25.440000)
nb.pair_coeff.set('13', '13', epsilon = 0.200000, sigma = 1.500000, lam = 0.000000, r_cut = 6.000000)
nb.pair_coeff.set('1', '13', epsilon = 0.200000, sigma = 3.000000, lam = 1.000000, r_cut = 12.000000)
nb.pair_coeff.set('2', '13', epsilon = 0.200000, sigma = 3.340000, lam = 1.000000, r_cut = 13.360000)
nb.pair_coeff.set('3', '13', epsilon = 0.200000, sigma = 3.840000, lam = 1.000000, r_cut = 15.360000)
nb.pair_coeff.set('4', '13', epsilon = 0.200000, sigma = 3.270000, lam = 1.000000, r_cut = 13.080000)
nb.pair_coeff.set('5', '13', epsilon = 0.200000, sigma = 3.760000, lam = 1.000000, r_cut = 15.040000)
nb.pair_coeff.set('6', '13', epsilon = 0.200000, sigma = 4.030000, lam = 1.000000, r_cut = 16.120000)
nb.pair_coeff.set('7', '13', epsilon = 0.200000, sigma = 3.590000, lam = 1.000000, r_cut = 14.360000)
nb.pair_coeff.set('8', '13', epsilon = 0.200000, sigma = 3.930000, lam = 1.000000, r_cut = 15.720000)
nb.pair_coeff.set('9', '13', epsilon = 0.200000, sigma = 3.540000, lam = 1.000000, r_cut = 14.160000)
nb.pair_coeff.set('10', '13', epsilon = 0.200000, sigma = 3.980000, lam = 1.000000, r_cut = 15.920000)
nb.pair_coeff.set('11', '13', epsilon = 0.200000, sigma = 3.530000, lam = 1.000000, r_cut = 14.120000)
nb.pair_coeff.set('12', '13', epsilon = 0.200000, sigma = 3.930000, lam = 1.000000, r_cut = 15.720000)

### Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut = 0.0, nlist = nl)
yukawa.pair_coeff.set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13'], ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13'], epsilon=0.000000, kappa=0.000000, r_cut=0.000000)
yukawa.pair_coeff.set('6','6', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('6','9', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('6','12', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('9','6', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('9','9', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('9','12', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('12','6', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('12','9', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('12','12', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)

## Integration
#Particles MD Integrator
solute = hoomd.group.all()
probe = hoomd.group.type('13')

#Zeroing the linear momentum to avoid any possible COM drift
fix_mum=hoomd.md.update.zero_momentum(period=1)

## Running 2 microsecond long LD simulation
hoomd.md.integrate.mode_standard(dt=dt) #time step is 10fs
integrator2 = hoomd.md.integrate.langevin(group = solute, kT = kT, seed = np.random.randint(1,1000000)) #Temp is kT/0.0019872067
integrator2.set_gamma('1',gamma=0.002790)
integrator2.set_gamma('2',gamma=0.004258)
integrator2.set_gamma('3',gamma=0.006415)
integrator2.set_gamma('4',gamma=0.003476)
integrator2.set_gamma('5',gamma=0.006264)
integrator2.set_gamma('6',gamma=0.007638)
integrator2.set_gamma('7',gamma=0.005579)
integrator2.set_gamma('8',gamma=0.007198)
integrator2.set_gamma('9',gamma=0.005628)
integrator2.set_gamma('10',gamma=0.007980)
integrator2.set_gamma('11',gamma=0.004749)
integrator2.set_gamma('12',gamma=0.006268)
integrator2.set_gamma('13',gamma=0.004890)

## writing outputs and run the simulation
eqthermo = hoomd.analyze.log(filename='afteq_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
hoomd.run(run_eq)
eqthermo.disable()

#Very short, short, and long simulation runs
longthermo = hoomd.analyze.log("long_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = long_period*5, overwrite = True, header_prefix = '#')
longgsd_file = hoomd.dump.gsd("long_dumps.gsd", period = long_period*4, dynamic = ['property','momentum'], group = solute, overwrite=True)
longgsd_file_probe = hoomd.dump.gsd("long_dumps_probe.gsd", period = long_period/20, dynamic = ['property','momentum'], group = probe, overwrite=True)
writegsd = hoomd.dump.gsd("long_restart.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True, truncate=True)
hoomd.run(tsteps=longrun_nve)
longthermo.disable()
longgsd_file.disable()

shortthermo = hoomd.analyze.log("short_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
#shortgsd_file = hoomd.dump.gsd("short_dumps.gsd", period = short_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
shortgsd_file_probe = hoomd.dump.gsd("short_dumps_probe.gsd", period = short_period/20, dynamic = ['property','momentum'], group = probe, overwrite=True)
hoomd.run(tsteps=shortrun_nve)
shortthermo.disable()
shortgsd_file.disable()

veryshortthermo = hoomd.analyze.log("veryshort_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = veryshort_period*10, overwrite = True, header_prefix = '#')
#veryshortgsd_file = hoomd.dump.gsd("veryshort_dumps.gsd", period = veryshort_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
veryshortgsd_file_probe = hoomd.dump.gsd("veryshort_dumps_probe.gsd", period = veryshort_period/20, dynamic = ['property','momentum'], group = probe, overwrite=True)
hoomd.run(tsteps=veryshortrun_nve)
veryshortthermo.disable()
veryshortgsd_file.disable()