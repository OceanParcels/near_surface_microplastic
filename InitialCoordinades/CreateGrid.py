"""
Create a uniformly spaced (lon,lat) grid of initial particle locations based on nemo bathymetry
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap

griddir = './'

spacing = 0.2 #spacing between particles
name = '02'

plotspacing=2 #For binning of final plot
outdir=griddir + 'Global' + name + 'grid/'

def create_particles():
    #Create uniform grid of particles
    filename=griddir + 'bathy_level.nc'
    data = Dataset(filename,'r')
    bathy=np.array(data['Bathy_level'][0])
    lon=np.array([data['nav_lon']][0])
    lat=np.array([data['nav_lat']][0])
        
    grid=np.mgrid[-180:180:spacing,-90:90:spacing]
    n=grid[0].size;
    lons=np.reshape(grid[0],n)
    lats=np.reshape(grid[1],n)
       
    bathy_points = griddata(np.array([lon.flatten(), lat.flatten()]).T, bathy.flatten(), (lons, lats), method='nearest')
    
    lons_new=np.array([lons[i] for i in range(len(lons)) if bathy_points[i]!=0])
    lats_new=np.array([lats[i] for i in range(len(lats)) if bathy_points[i]!=0])
    
    lons_new[lons_new<0]+=360
        
    np.save(outdir+'Lons_full' + str(name),lons_new)
    np.save(outdir+'Lats_full' + str(name),lats_new)

#create_particles()

def Plot_particles():
    #Plot to check if everything went well
    lons=np.load(outdir + 'Lons_full' + str(name) + '.npy')
    lats=np.load(outdir + 'Lats_full' + str(name) + '.npy')
    
    assert (len(lons)==len(lats))
    
    print 'Number of particles: ', len(lons)
    fig = plt.figure(figsize=(25, 30))
    ax = fig.add_subplot(211)
    ax.set_title("Particles")
    
    m = Basemap(projection='robin',lon_0=180,resolution='c')
    m.drawcoastlines()
    xs, ys = m(lons, lats)
    m.scatter(xs,ys)
    
    ax = fig.add_subplot(212)
    ax.set_title("Particles per bin. Should be constant everywhere but on land.")

    m = Basemap(projection='robin',lon_0=180,resolution='c')
    m.drawcoastlines()

    lon_bin_edges = np.arange(0, 360+spacing, plotspacing)
    lat_bins_edges = np.arange(-90, 90+spacing, plotspacing)
    
    density, _, _ = np.histogram2d(lats, lons, [lat_bins_edges, lon_bin_edges])
    
    lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bin_edges, lat_bins_edges)
    xs, ys = m(lon_bins_2d, lat_bins_2d)

    plt.pcolormesh(xs, ys, density,cmap=plt.cm.RdBu_r)
    cbar = plt.colorbar(orientation='vertical', shrink=0.625, aspect=20, fraction=0.2,pad=0.02)
    cbar.set_label('Particles per bin',size=8)

#Plot_particles()

def split_grid():
    #split up the big grid into 40 smaller grids for parcels to be faster    
    for name in ['02']:
        Lons = np.load('Global' + name + 'grid//Lons_full' + name + '.npy')
        Lats = np.load('Global' + name + 'grid//Lats_full' + name + '.npy')
       
        print 'Total number of particles: ', len(Lons)
        
        outdir='./Global' + name + 'grid/'
        N=40 #Number of sub-grids
        
        k = len(Lons)//N+1 #Number of particles per file
        print k
        
        for i in range(0,len(Lons)//k+1):
            lo = Lons[i*k:(i+1)*k]
            la = Lats[i*k:(i+1)*k]
            np.save(outdir + 'Lons' + str(i), lo)
            np.save(outdir + 'Lats' + str(i), la)
            print 'lo shape: ', lo.shape
            print 'la shape: ', la.shape

#split_grid()
