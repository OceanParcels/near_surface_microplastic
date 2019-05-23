"""
Potential influence of near-surface currents on the global dispersal of marine microplastic
-------------------------------------------------------------------------
David Wichmann, Philippe Delandmeter, Erik van Sebille

d.wichmann@uu.nl

##########################################################################

Objects for data analysis, used for creating figures
"""

import numpy as np
from netCDF4 import Dataset

#Domain of interest
minlon=0.
maxlon=360.
minlat=-90.
maxlat=90


def regions(ddeg, so_lat=-60, na_lat=65):
    """
    Function to return the different region definitions as array
    - ddeg: binning (square, in degree)
    - so_lat: latitutde separating the Southern Ocean from the southern basins
    - na_lat: latitude separating North Atlantic and Arctic regions
    """
    
    Lons = np.arange(0.,360.,ddeg)
    Lats = np.arange(-90.,90.,ddeg)
    
    basinregion={'NA1': [260,360,0,na_lat],
                 'NA2': [0,30,0,na_lat],
                 'NA3': [30,55,30,na_lat],
                 'NP': [115,260,0,65],
                 'SA1': [290,360,so_lat,0],
                 'SA2': [0,20,so_lat,10],
                 'SP': [130,290,so_lat,0],
                 'IO': [20,130,so_lat,30],
                 'SO': [0,360,-90,so_lat],
                 'PNA1': [260,360,na_lat,90],
                 'PNA2': [0,100,na_lat,90],
                 'PNP': [100,260,65,90]}

    region = np.zeros((len(Lons),len(Lats)))

    for i in range(len(Lats)):
        for j in range(len(Lons)):

            la=Lats[i]
            lo=Lons[j]
            
            (NPminlon, NPmaxlon, NPminlat, NPmaxlat)=basinregion['NP']
            if (lo<NPmaxlon and lo>=NPminlon and la<NPmaxlat and la>=NPminlat or (lo>=260 and lo<285 and la<-4./5.*lo+230) or (lo>=250 and lo<290 and la>=0 and la<8.) or (lo>=280 and lo<281 and la>=0 and la<9.)):
                region[j,i]=1
            (NAminlon, NAmaxlon, NAminlat, NAmaxlat)=basinregion['NA1']
            if (lo<NAmaxlon and lo>=NAminlon and la<NAmaxlat and la>=NAminlat and region[j,i]==0):
                region[j,i]=2
            (NAminlon, NAmaxlon, NAminlat, NAmaxlat)=basinregion['NA2']
            if (lo<NAmaxlon and lo>=NAminlon and la<NAmaxlat and la>=NAminlat and region[j,i]==0):
                region[j,i]=2
            (NAminlon, NAmaxlon, NAminlat, NAmaxlat)=basinregion['NA3']
            if (lo<NAmaxlon and lo>=NAminlon and la<NAmaxlat and la>=NAminlat and region[j,i]==0):
                region[j,i]=2
            (SPminlon, SPmaxlon, SPminlat, SPmaxlat)=basinregion['SP']
            if (lo<SPmaxlon and lo>=SPminlon and la<SPmaxlat and la>=SPminlat):
                region[j,i]=3
            (SAminlon, SAmaxlon, SAminlat, SAmaxlat)=basinregion['SA1']
            if (lo<SAmaxlon and lo>=SAminlon and la<SAmaxlat and la>=SAminlat):
                region[j,i]=4
            (SAminlon, SAmaxlon, SAminlat, SAmaxlat)=basinregion['SA2']
            if (lo<SAmaxlon and lo>=SAminlon and la<SAmaxlat and la>=SAminlat and region[j,i]==0):
                region[j,i]=4
            (IOminlon, IOmaxlon, IOminlat, IOmaxlat)=basinregion['IO']
            if (lo<IOmaxlon and lo>=IOminlon and la<IOmaxlat and la>=IOminlat and region[j,i]==0):
                region[j,i]=5                
            (PNPminlon, PNPmaxlon, PNPminlat, PNPmaxlat)=basinregion['PNP']
            if (lo<PNPmaxlon and lo>=PNPminlon and la<PNPmaxlat and la>=PNPminlat):
                region[j,i]=6
            (PNA1minlon, PNA1maxlon, PNA1minlat, PNA1maxlat)=basinregion['PNA1']
            if (lo<PNA1maxlon and lo>=PNA1minlon and la<PNA1maxlat and la>=PNA1minlat):
                region[j,i]=6
            (PNA2minlon, PNA2maxlon, PNA2minlat, PNA2maxlat)=basinregion['PNA2']
            if (lo<PNA2maxlon and lo>=PNA2minlon and la<PNA2maxlat and la>=PNA2minlat):
                region[j,i]=6
            (SOminlon, SOmaxlon, SOminlat, SOmaxlat)=basinregion['SO']
            if (lo<SOmaxlon and lo>=SOminlon and la<SOmaxlat and la>=SOminlat):
                region[j,i]=7


    region_names = ['North Pacific','North Atlantic','South Pacific','South Atlantic','Indian Ocean','Artic','Southern Ocean']
    region=region.T
    region=region.ravel()
    return region, region_names


class ParticleData(object):
    """
    Class that containing 2D particle data and functions to analyse it
    """

    def __init__(self, lons, lats, times, depths=None):
        """
        -params lons, lats, times, depths: arrays containing the data
        """
        
        print '---------------------'
        print 'Particle data created'
        print '---------------------'
        print 'Particles: ', len(lons)
        print 'Snapshots: ', len(lons[0])
        print '---------------------'
        
        self.lons=lons
        self.lats=lats
        self.times=times
        self.depths=depths
        
    def __del__(self):
        print "Particle Data deleted"
        
    def remove_nans(self): #For removing Nans in problematic regions
        print 'Removing NaNs...'
        nan_entries = np.argwhere(np.isnan(self.lons))[:,0]
        indices = [i for i in range(len(self.lons)) if i not in nan_entries]
        print 'Removed number of NaN values: ', len(self.lons)-len(indices)
        self.lons = self.lons[indices]
        self.lats = self.lats[indices]
        self.times = self.times[indices]
        
        if self.depths is not None:
            self.depths = self.depths[indices]
        
        print 'NaNs are removed'
        print '---------------------'
        
    def get_distribution(self, t, ddeg):
        """
        Calculate the particle distribution at time t. 
        - param t: integer time from loaded particles
        - param ddeg: binning
        """
        
        lon_edges=np.linspace(minlon,maxlon,int((maxlon-minlon)/ddeg)+1)        
        lat_edges=np.linspace(minlat,maxlat,int((maxlat-minlat)/ddeg)+1)  
        d , _, _ = np.histogram2d(self.lats[:,t], self.lons[:,t], [lat_edges, lon_edges])
        return d

    @classmethod
    def from_nc(cls, pdir, fname, tload=None, Ngrids=40):
        """
        Load 2D data from netcdf particle output files. We assume that there are 
        several output files, each for different initial distributions that have to be merged
        :param pdir: directory of files
        :param fname: file name in pdir
        :param tload: array of times for which we load the data (indices, not actual times)
        :Ngrids: number of different output files to be merged (40 in our case)
        """

        print 'Loading data from files: ', pdir + fname
        
        #Load data from first grid-array
        i = 0
        print 'Load grid no: ', i
        pfile = pdir + fname + str(i)+'.nc'     
        data = Dataset(pfile,'r')
        
        times=data.variables['time'][:,tload]
        lons=data.variables['lon'][:,tload]
        lats=data.variables['lat'][:,tload]
        
        #Load data from other grid-arrays        
        for i in range(1,Ngrids):
            print 'Load grid no: ', i
            pfile = pdir + fname + str(i)+'.nc'  
            data = Dataset(pfile,'r')
            times=np.vstack((times, data.variables['time'][:,tload]))
            lons=np.vstack((lons, data.variables['lon'][:,tload]))
            lats=np.vstack((lats, data.variables['lat'][:,tload]))
        times/=86400. #Convert to days

        return cls(lons=lons, lats=lats, times=times)
    


    @classmethod
    def from_nc_3d(cls, pdir, fname, tload=None, Ngrids=40):
        """
        Load 2D data from netcdf particle output files. We assume that there are 
        several output files, each for different initial distributions that have to be merged
        :param pdir: directory of files
        :param fname: file name in pdir
        :param tload: array of times for which we load the data (indices, not actual times)
        :Ngrids: number of different output files to be merged (40 in our case)
        """

        print 'Loading data from files: ', pdir + fname
        
        #Load data from first grid-array
        i = 0
        print 'Load grid no: ', i
        pfile = pdir + fname + str(i)+'.nc'     
        data = Dataset(pfile,'r')
        
        times=data.variables['time'][:,tload]
        lons=data.variables['lon'][:,tload]
        lats=data.variables['lat'][:,tload]
        depths=data.variables['z'][:,tload]
        
        #Load data from other grid-arrays        
        for i in range(1,Ngrids):
            print 'Load grid no: ', i
            pfile = pdir + fname + str(i)+'.nc'  
            data = Dataset(pfile,'r')
            times=np.vstack((times, data.variables['time'][:,tload]))
            lons=np.vstack((lons, data.variables['lon'][:,tload]))
            lats=np.vstack((lats, data.variables['lat'][:,tload]))
            depths=np.vstack((depths, data.variables['z'][:,tload]))
        times/=86400. #Convert to days

        return cls(lons=lons, lats=lats, times=times, depths=depths)
    
    def set_region_labels(self, ddeg, t, so_lat=-60, na_lat=65):
        """
        Function to give particles labels according to the region they started in
        -param ddeg: binning
        -param t: time for labelling
        """
        r, region_names = regions(ddeg, so_lat=so_lat, na_lat=na_lat)
        
        N=360//ddeg
        label=np.array([int(((la-minlat)//ddeg)*N+(lo-minlon)//ddeg) for la,lo in zip(self.lats[:,t],self.lons[:,t])])
        
        region_label = [r[label[i]] for i in range(len(label))]
        return region_label
        
    

class oceanvector(object):
    """
    Class for ocean vectors. Can be 1d or 2d (lon x lat)
    """
    
    def __init__(self,vec,minlon=minlon,maxlon=maxlon, minlat=minlat,maxlat=maxlat,ddeg=2.):
        """
        -val: A value (e.g. an eigenvalue) which is written as a plot title for figures
        """
        
        self.minlon=minlon
        self.maxlon=maxlon
        self.minlat=minlat
        self.maxlat=maxlat
        self.ddeg = ddeg
        
        #Bin Edges
        self.Lons_edges=np.linspace(minlon,maxlon,int((maxlon-minlon)/ddeg)+1)        
        self.Lats_edges=np.linspace(minlat,maxlat,int((maxlat-minlat)/ddeg)+1)
        
        #Bin centers. This is the format of a field as well.
        self.Lons_centered=np.array([(self.Lons_edges[i]+self.Lons_edges[i+1])/2. for i in range(len(self.Lons_edges)-1)])
        self.Lats_centered=np.array([(self.Lats_edges[i]+self.Lats_edges[i+1])/2. for i in range(len(self.Lats_edges)-1)])        
        
        if vec.ndim==1:
            v1d = vec
            v2d = vec.reshape((len(self.Lats_centered),len(self.Lons_centered)))
        else:
            v1d = vec.ravel()
            v2d = vec
        
        self.V1d = v1d
        self.V2d = v2d

