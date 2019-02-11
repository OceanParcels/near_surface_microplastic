"""
Influence of near-surface currents on the global dispersal of marine microplastic
-------------------------------------------------------------------------
David Wichmann, Philippe Delandmeter, Erik van Sebille

d.wichmann@uu.nl

##########################################################################

Code for the computation of particle trajectories
Advection of particles with nemo 1/12 degree. Possible simulation modes are
1. Simulations at fixed depth
2. Simulation with uniform mixing 
3. Simulation with wind-stress based mixing of Kukulka

The latter two are constrained to the interval [0,120] m.
"""

import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4
from argparse import ArgumentParser
from datetime import timedelta
from datetime import datetime
from glob import glob
import math

datadir = '/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/' #Directory for nemo data
outputdir = '/scratch-shared/wichmann/SubSurfaceOutput/' #Directory for output files
griddir = '/home/wichmann/SubsurfaceTransport/ParticleGrid/Global02grid/' #Directory for initial particle distribution

def kernel_uniform_mixing(particle, fieldset, time, dt):  
    """
    Kernel that randomly allocates a particle in the vertical according to a uniform distribution
    """
    particle.depth=random.uniform(0.,120.)

def kernel_kukulka_mixing(particle, fieldset, time, dt):  
    """
    :Kernel that randomly distributes particles along the vertical according to an expovariate distribution.
    :Parameterization according to Kukulka et al. (2012): The effect of wind mixing on the vertical distribution of buoyant plastic debris
    :Comment on dimensions: tau needs to be in Pa
    """
    stress =  fieldset.TAU[particle.time,particle.lon,particle.lat, 0.]
    A0=0.31 * math.pow(stress,1.5)
    l=fieldset.wrise/A0
    d=random.expovariate(l) #Kukulka formula. Used depths of [0 ... 120.] m
    if d>120.:
        particle.depth=120.
    else:
        particle.depth=d

def DeleteParticle(particle, fieldset, time, dt):
    """Kernel for deleting particles if they are out of bounds."""
    particle.delete()

def periodicBC(particle, fieldset, time, dt):
    """
    Kernel for periodic values in longitude
    """
    if particle.lon < 0.:
        particle.lon += 360.
    elif particle.lon > 360.:
        particle.lon -= 360.


def p_advect(ptype=JITParticle, outname='noname', pos=0, y=2001, m=1, d=1, simdays=90, particledepth=0, uniform_mixing=False, kukulka_mixing=False, wrise=None):
    """
    Main function for execution
        - outname: name of the output file. Note that all important parameters are also in the file name.
        - pos: Execution is manually parallelized over different initial position grids. These are indexed.
        - y, m, d: year, month an day of the simulation start
        - simdays: number of days to simulate
        - particledepth: for fixed-depth simulations. Index of nemo depth grid
        - uniform_mixing: True for uniform mixing simulation
        - kukulka_mixing: True for kukulka mixing simulation
        - wrise: Rise velocity, needed for kukulka mixing simulation
    """
    
    print '-------------------------'
    print 'Start run... Parameters: '
    print '-------------------------'
    print 'Initial time (y, m, d): ', (y, m, d)
    print 'Simulation days', simdays
    print '-------------------------'
    
    #Load initial particle positions (grids) from external files
    lons = np.load(griddir + 'Lons' + str(pos) + '.npy')
    lats = np.load(griddir + 'Lats' + str(pos) + '.npy') 
    times = [datetime(y, m, d)]*len(lons)
    print 'Number of particles: ', len(lons)

    ufiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05U.nc'))
    vfiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05V.nc'))
    wfile =  datadir+'means/ORCA0083-N06_20000105d05W.nc' #file for the definition of depths
    
    mesh_mask = datadir + 'domain/coordinates.nc'

    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfile, 'data': ufiles},
                 'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfile, 'data': vfiles}}

    variables = {'U': 'uo', 'V': 'vo'}
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                  'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}

    outfile = outputdir + outname + '_y'+ str(y) + '_m' + str(m) + '_d' + str(d)  + '_simdays' + str(simdays) + '_layer' + str(particledepth) + '_pos' + str(pos)

    #Select indices w.r.t. depth to load only the needed subsets of data.
    if not (uniform_mixing or kukulka_mixing):
        indices = {'depth': [particledepth]}
    else:
        indices = {'depth': [i for i in range(27)]}

    if  kukulka_mixing:
         taufiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05T.nc'))
         filenames['TAU'] = {'lon': mesh_mask, 'lat': mesh_mask, 'data': taufiles}
         variables['TAU']='taum'
         dimensions['TAU'] = {'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_centered'}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions,  allow_time_extrapolation=False, indices=indices)

    if not (uniform_mixing or kukulka_mixing):
        d=fieldset.U.depth #We loaded data from only one depth, which is set to the particle depth.
        depths = [d]*len(lons)
    else:
        depths = [0.]*len(lons) #Place particles initially at the surface for mixing simulations

    fieldset.U.vmax = 10
    fieldset.V.vmax = 10

    pset = ParticleSet(fieldset=fieldset, pclass=ptype, lon=lons, lat=lats, time=times, depth=depths)

    kernels= pset.Kernel(AdvectionRK4) + pset.Kernel(periodicBC)

    if uniform_mixing:
        print 'UNIFORM MIXING'
        kernels+=pset.Kernel(kernel_uniform_mixing)
        outfile+='_uniform_mix'
    elif kukulka_mixing:
        print 'KUKULKA MIXING'
        print 'wrise: ', wrise
        fieldset.wrise = wrise
        fieldset.TAU.vmax = 10
        kernels+=pset.Kernel(kernel_kukulka_mixing)
        outfile+='_kukulka_mix'
    else:
        print 'FIXED DEPTH'

    #Trajectory computation
    pset.execute(kernels, runtime=timedelta(days=simdays), dt=timedelta(minutes=10), 
                 output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(days=15)),
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}, verbose_progress=False)

if __name__=="__main__":
    ptype = {'scipy': ScipyParticle, 'jit': JITParticle}
    p = ArgumentParser(description="""Global advection of different particles""")
    p.add_argument('-ptype', '--ptype',choices=('scipy', 'jit'), nargs='?', default='jit',help='execution mode')
    p.add_argument('-name', '--name', default='noname',help='name of output file')
    p.add_argument('-posidx', '--posidx', type=int,default=0,help='label of lon/lat initial array')
    p.add_argument('-y', '--y', type=int,default=2000,help='year of simulation start')
    p.add_argument('-m', '--m', type=int,default=1,help='month of simulation start')
    p.add_argument('-d', '--d', type=int,default=5,help='day of simulation start')
    p.add_argument('-simdays', '--simdays', type=int,default=100,help='simulation days')
    p.add_argument('-depth', '--depth', type=int,default=0,help='index of depth layer in nemo')
    p.add_argument('-uniformmixing', '--uniformmixing', type=bool,default=False,help='True if we simulate uniform mixing')
    p.add_argument('-kukulkamixing', '--kukulkamixing', type=bool,default=False,help='True if we simulate kukulka mixing')
    p.add_argument('-wrise', '--wrise', type=float,default=None,help='R=rise velocity for Kukulka mixing')
    args = p.parse_args()
    p_advect(ptype=ptype[args.ptype],outname=args.name, pos=args.posidx, y=args.y, m=args.m, d=args.d, simdays=args.simdays, particledepth=args.depth, uniform_mixing=args.uniformmixing, kukulka_mixing = args.kukulkamixing, wrise=args.wrise)