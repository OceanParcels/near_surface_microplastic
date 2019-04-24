"""
Influence of near-surface currents on the global dispersal of marine microplastic
-------------------------------------------------------------------------
David Wichmann, Philippe Delandmeter, Erik van Sebille

d.wichmann@uu.nl

##########################################################################

Code for the computation of 3D passive particle trajectories

"""

import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, ErrorCode, AdvectionRK4_3D
from argparse import ArgumentParser
from datetime import timedelta
from datetime import datetime
from glob import glob

datadir = '/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/' #Directory for nemo data
outputdir = '/scratch-shared/wichmann/SubSurfaceOutput/' #Directory for output files
griddir = '/home/wichmann/SubsurfaceTransport/ParticleGrid/Global02grid/' #Directory for initial particle distribution

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

def p_advect(ptype=JITParticle, outname='noname', pos=0, y=2001, m=1, d=1, simdays=90):
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
    depths = [1.5] * len(lons)
    times = [datetime(y, m, d)]*len(lons)
    print 'Number of particles: ', len(lons)

    ufiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05U.nc'))
    vfiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05V.nc'))
    wfiles =  sorted(glob(datadir+'means/ORCA0083-N06_200?????d05W.nc'))
    mesh_mask = datadir + 'domain/coordinates.nc'

    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
             'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}}

    variables = {'U': 'uo',
             'V': 'vo',
             'W': 'wo'}
    
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)

    outfile = outputdir + outname + '3D_y'+ str(y) + '_m' + str(m) + '_d' + str(d)  + '_simdays' + str(simdays) + '_pos' + str(pos)

    fieldset.U.vmax = 10
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10

    pset = ParticleSet(fieldset=fieldset, pclass=ptype, lon=lons, lat=lats, time=times, depth=depths)

    kernels= pset.Kernel(AdvectionRK4_3D) + pset.Kernel(periodicBC)
    
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
    args = p.parse_args()
    p_advect(ptype=ptype[args.ptype],outname=args.name, pos=args.posidx, y=args.y, m=args.m, d=args.d, simdays=args.simdays)