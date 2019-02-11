"""
Influence of near-surface currents on the global dispersal of marine microplastic
-------------------------------------------------------------------------
David Wichmann, Philippe Delandmeter, Erik van Sebille

d.wichmann@uu.nl

##########################################################################

Code to create figures and tables in paper and supplements
"""

import numpy as np
from AnaObjects import ParticleData, oceanvector, regions
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

datadir = '/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' #Data directory #directory of the data.
outdir_paper = './paper_figures/' #directory for saving figures for the paper
outdir_supplementary = './supplementary_material/' #directory for saving figures for the paper

nemo_depth =[0, 1.023907, 2.10319, 3.251309, 4.485053, 5.825238, 7.297443, 
             8.932686, 10.7679, 12.84599, 15.21527, 17.92792, 21.03757, 24.59599, 
             28.64965, 33.23697, 38.3871, 44.12101, 50.45447, 57.40257, 64.9846, 
             73.2287, 82.17556, 91.88141, 102.4202, 113.8852, 126.3909]

#Choose middle of the cells. Note that with C-grid interpolation, the horizontal velocities are the same within each cell
nemo_depth = [int((nemo_depth[i+1]+nemo_depth[i])/2) for i in range(len(nemo_depth)-1)]

models_all = ['Layer0', 'Layer2', 'Layer4', 'Layer7', 'Layer10','Layer13',
          'Layer16', 'Layer19', 'Layer22', 'Layer23', 'Layer25', 'Uniform', 
          'Kukulka']
filenames_all = ['SubSurf_y2000_m1_d5_simdays3650_layer0_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer2_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer4_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer7_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer10_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer13_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer16_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer19_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer22_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer23_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer25_pos', 'Uniform_y2000_m1_d5_simdays3650_layer0_pos',
              'Kukulka_y2000_m1_d5_simdays3650_layer0_pos']
folders_all = ['Layer0', 'Layer2', 'Layer4', 'Layer7', 'Layer10','Layer13',
           'Layer16', 'Layer19', 'Layer22', 'Layer23', 'Layer25', 'UniMix',
           'KukMix']
folders_all = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders_all]


"""
FIGURES PAPER
"""

def FIG1_surface_distribution():
    
    pdir = '/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/Layer0/' #Data directory
    filename = 'SubSurf_y2000_m1_d5_simdays3650_layer0_pos'
    
    #load particle data and get distribution
    pdata=ParticleData.from_nc(pdir, filename, tload=[0,-1], Ngrids=40)
    d_full=pdata.get_distribution(t=-1, ddeg=2.).flatten()
    d=oceanvector(d_full, ddeg = 2.)        
    lon_bins_2d,lat_bins_2d = np.meshgrid(d.Lons_edges, d.Lats_edges)
    
    #plot figure
    fig = plt.figure(figsize = (12,8))
    m = Basemap(projection='robin', lon_0=180, resolution='l')
    m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=1.2, size=10)
    m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=1.2, size=10)
    m.drawcoastlines(linewidth=0.4)
    m.fillcontinents(color='dimgrey')
    xs, ys = m(lon_bins_2d, lat_bins_2d) 
    plt.pcolormesh(xs, ys, d.V2d,cmap='plasma', norm=colors.LogNorm(), rasterized=True)
    cbar=plt.colorbar(orientation='vertical',shrink=0.5)
    cbar.ax.tick_params(labelsize=10, width=0.05)
    cbar.set_label('Particles per bin', size=10)
    fig.savefig(outdir_paper + 'surface_distribution.pdf', dpi=900, bbox_inches='tight')

FIG1_surface_distribution()

def FIG2_distributions_different_models():
    
    folders = ['Layer10', 'Layer16', 'Layer19', 'UniMix']
    folders = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders]
    labels = ['a) z = ' + str(nemo_depth[10]) + ' m', 
    'b) z = ' + str(nemo_depth[16]) + ' m', 
    'c) z = ' + str(nemo_depth[19]) + ' m', 
    'd) Uniform'] #figure titles    
    filenames = ['SubSurf_y2000_m1_d5_simdays3650_layer10_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer16_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer19_pos', 'Uniform_y2000_m1_d5_simdays3650_layer0_pos']
    
    fig = plt.figure(figsize = (14,8)) 
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(wspace=0.01, hspace=0.0)

    for i in range(len(folders)):

        #load particle data and get distribution
        filename = filenames[i]
        pdata=ParticleData.from_nc(folders[i] + '/', filename, tload=[0,-1], Ngrids=40)
        d_full=pdata.get_distribution(t=-1, ddeg=2.).flatten()
        d=oceanvector(d_full, ddeg = 2.)        
        lon_bins_2d,lat_bins_2d = np.meshgrid(d.Lons_edges,d.Lats_edges)

        #create subplot
        plt.subplot(gs1[i])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.drawcoastlines()
        m.fillcontinents(color='dimgrey')
        xs, ys = m(lon_bins_2d, lat_bins_2d)         
        plt.pcolormesh(xs, ys, d.V2d,cmap='plasma', norm=colors.LogNorm(), rasterized=True)
        cbar=plt.colorbar(orientation='vertical',shrink=0.5)
        cbar.ax.tick_params(labelsize=7, width=0.03)
        plt.title(labels[i], size=10, y=1.)
    
    fig.savefig(outdir_paper + 'distributions_different_models.pdf', dpi=900, bbox_inches='tight')

FIG2_distributions_different_models()

def FIG3_scatter_final_basin():

    folders = ['Layer0', 'Layer16', 'Layer25', 'UniMix']
    folders = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders]
    labels = ['a) z = ' + str(nemo_depth[0]) + ' m', 
    'b) z = ' + str(nemo_depth[16]) + ' m', 
    'c) z = ' + str(nemo_depth[25]) + ' m', 
    'd) Uniform'] #figure titles    
    filenames = ['SubSurf_y2000_m1_d5_simdays3650_layer0_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer16_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer25_pos', 'Uniform_y2000_m1_d5_simdays3650_layer0_pos']

    #for the figure
    fig = plt.figure(figsize = (14,8)) 
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(wspace=0.1, hspace=0.0)
    
    #get region names and constrain to the subtropical basins
    _, region_names = regions(2.)
    region_names=region_names[0:5]
    basin_regions = [1,2,3,4,5]

    for j in range(len(folders)):
        
        #load particle data        
        filename = filenames[j]
        pdata=ParticleData.from_nc(folders[j] + '/', filename, tload=[0,-1], Ngrids=40)
        pdata.remove_nans()
        region_label=pdata.set_region_labels(2., -1)
        
        #select the particles in the basins, and label the particles according to their final basin
        selection = [k for k in range(len(region_label)) if region_label[k] in basin_regions]
        lons=pdata.lons[selection]
        lats=pdata.lats[selection]
        region_label = [region_label[i] for i in range(len(region_label)) if region_label[i] in basin_regions]
        
        #plot
        plt.subplot(gs1[j])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.fillcontinents(color='dimgray')
        xs, ys = m(lons[:,0], lats[:,0])
        cmap = plt.get_cmap('jet', 5)
        scat=m.scatter(xs, ys, marker='.', s=1, c=region_label, cmap=cmap, rasterized=True)
        plt.title(labels[j], size=10, y=1.)
    
    #color bar on the right
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.822, 0.25, 0.015, 0.5])
    cbar=fig.colorbar(scat, cax=cbar_ax, ticks=[1,2,3,4,5])
    cbar.ax.tick_params(labelsize=10)
    tick_locs = 1+0.4 * np.array([1,3,5,7,9])
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(region_names)
    fig.savefig(outdir_paper + 'transport_basins_scatter.pdf', dpi=900, bbox_inches='tight')

FIG3_scatter_final_basin()

def FIG4_scatter_final_poles():

    folders = ['Layer0', 'Layer16', 'Layer25', 'UniMix']
    folders = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders]
    labels = ['a) z = ' + str(nemo_depth[0]) + ' m', 
    'b) z = ' + str(nemo_depth[16]) + ' m', 
    'c) z = ' + str(nemo_depth[25]) + ' m', 
    'd) Uniform'] #figure titles    
    filenames = ['SubSurf_y2000_m1_d5_simdays3650_layer0_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer16_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer25_pos', 'Uniform_y2000_m1_d5_simdays3650_layer0_pos']

    #for the figure
    fig = plt.figure(figsize = (14,8)) 
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(wspace=0.1, hspace=0.0)

    #binning for the regions defining functions. this does not affect the result, as the region boundaries are a multiple of 2
    ddeg=2.
    
    #get region names and constrain to the subtropical basins
    _, region_names = regions(ddeg)
    region_names=region_names[5:]
    basin_regions = [6,7]

    for j in range(len(folders)):
        
        #load particle data        
        filename = filenames[j]
        pdata=ParticleData.from_nc(folders[j] + '/', filename, tload=[0,-1], Ngrids=40)
        pdata.remove_nans()
        region_label=pdata.set_region_labels(ddeg, -1)
        
        #select the particles in the basins, and label the particles according to their final basin
        selection = [k for k in range(len(region_label)) if region_label[k] in basin_regions]
        lons=pdata.lons[selection]
        lats=pdata.lats[selection]
        region_label = [region_label[i] for i in range(len(region_label)) if region_label[i] in basin_regions]
        
        #plot
        plt.subplot(gs1[j])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.fillcontinents(color='dimgray')
        xs, ys = m(lons[:,0], lats[:,0])
        cmap = plt.get_cmap('jet', 2)
        scat=m.scatter(xs, ys, marker='.', s=1, c=region_label, cmap=cmap, rasterized=True)
        plt.title(labels[j], size=10, y=1.)
    
    #color bar on the right
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.822, 0.4, 0.015, 0.2])
    cbar=fig.colorbar(scat, cax=cbar_ax, ticks=[6,7])
    cbar.ax.tick_params(labelsize=10)
    tick_locs = 6+0.25 * np.array([1,3])
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(region_names)

    #save figure
    fig.savefig(outdir_paper + 'transport_poles_scatter.pdf', dpi=900, bbox_inches='tight')

FIG4_scatter_final_poles()


def FIG5_RegionalTransport():
    """
    We first create the transport matrix n_{ij} for each layer and save it. Then create the figures from it.
    Function create_transport_matrix needs to be executed first.
    """

    def create_transport_matrix(so_lat=-60, na_lat=65):
        n_matrix={} #Matrix with entries n_{ij} in paper as a list for each layer
        
        ddeg = 2. #Choice of binning, required for the region labels. ddeg=2 has no influence here on the result as we are perfectly on the region boundaries
        
        for j in range(len(folders_all)):            
            pdata=ParticleData.from_nc(folders_all[j] + '/', filenames_all[j], tload=[0,-1], Ngrids=40)
            pdata.remove_nans() #Remove some Nans
            
            #Give the particles labels according to the basin they are in at different times
            initial_region = pdata.set_region_labels(ddeg, 0, so_lat=so_lat, na_lat=na_lat)
            final_region = pdata.set_region_labels(ddeg, -1, so_lat=so_lat, na_lat=na_lat)
            
            n=np.zeros((8,8)) #n-matrix for the specific layer
            
            for i in range(len(initial_region)):
                n[int(initial_region[i]),int(final_region[i])]+=1
    
            n_matrix[models_all[j]]=n
        
        np.save(outdir_paper + 'n_matrix_solat_'+str(so_lat) + '_na_lat_'+str(na_lat), n_matrix)
    
    for so_lat in [-56, -60, -62]:  
        create_transport_matrix(so_lat=so_lat)

    for na_lat in [60, 70]:  
        create_transport_matrix(na_lat=na_lat)
    
    
    def polward_transport():
        n_matrix = np.load(outdir_paper + 'n_matrix_solat_-60_na_lat_65.npy').tolist()

        n_matrix_so56 = np.load(outdir_paper + 'n_matrix_solat_-56_na_lat_65.npy').tolist()
        n_matrix_so62 = np.load(outdir_paper + 'n_matrix_solat_-62_na_lat_65.npy').tolist()

        n_matrix_na60 = np.load(outdir_paper + 'n_matrix_solat_-60_na_lat_60.npy').tolist()
        n_matrix_na70 = np.load(outdir_paper + 'n_matrix_solat_-60_na_lat_70.npy').tolist()

        layers = [0, 2, 4, 7, 10, 13, 16, 19, 22, 23, 25]
        depths = [nemo_depth[d] for d in layers]
        depths.append(140)
        depths.append(160)
        
        f_to_southern_ocean=[]
        f_to_southern_ocean_56=[]
        f_to_southern_ocean_62=[]
        
        f_to_arctic=[]
        f_to_arctic_60=[]
        f_to_arctic_70=[]
        
        _, region_names = regions(2.)
        
        for i in range(len(folders_all)):
            print models_all[i]
            n=n_matrix[models_all[i]]
            n56=n_matrix_so56[models_all[i]]
            n62=n_matrix_so62[models_all[i]]
            
            n60=n_matrix_na60[models_all[i]]
            n70=n_matrix_na70[models_all[i]]
                        
            #row-normalize the matrices x to get F
            
            #(so_lat,sa_lat) = (-60, 65)
            s=np.sum(n,axis=1)
            f=np.zeros(n.shape)            
            for j in range(len(n)):
                if s[j]!=0:
                    f[j:,]=n[j,:]/s[j]
            f_to_southern_ocean.append(f[:,7])
            f_to_arctic.append(f[:,6])
            
            #(so_lat,sa_lat) = (-56, 65)
            s=np.sum(n56,axis=1)
            f=np.zeros(n56.shape)            
            for j in range(len(n56)):
                if s[j]!=0:
                    f[j:,]=n56[j,:]/s[j]
            f_to_southern_ocean_56.append(f[:,7])

            #(so_lat,sa_lat) = (-62, 65)
            s=np.sum(n62,axis=1)
            f=np.zeros(n62.shape)            
            for j in range(len(n62)):
                if s[j]!=0:
                    f[j:,]=n62[j,:]/s[j]

            f_to_southern_ocean_62.append(f[:,7])

            #(so_lat,sa_lat) = (-60, 60)
            s=np.sum(n60,axis=1)
            f=np.zeros(n60.shape)            
            for j in range(len(n60)):
                if s[j]!=0:
                    f[j:,]=n60[j,:]/s[j]

            f_to_arctic_60.append(f[:,6])

            #(so_lat,sa_lat) = (-60, 70)
            s=np.sum(n70,axis=1)
            f=np.zeros(n70.shape)            
            for j in range(len(n70)):
                if s[j]!=0:
                    f[j:,]=n70[j,:]/s[j]

            f_to_arctic_70.append(f[:,6])
        
        #To percent
        f_to_southern_ocean=100*np.array(f_to_southern_ocean)
        f_to_southern_ocean_56=100*np.array(f_to_southern_ocean_56)
        f_to_southern_ocean_62=100*np.array(f_to_southern_ocean_62)
        f_to_arctic=100*np.array(f_to_arctic)
        f_to_arctic_60=100*np.array(f_to_arctic_60)
        f_to_arctic_70=100*np.array(f_to_arctic_70)
        
        
        #Create figures
        plt.figure(figsize = (12,10)) 
        gs1 = gridspec.GridSpec(2, 2)
        gs1.update(wspace=.25, hspace=.35)
        
        #SO F
        plt.subplot(gs1[0])
        plt.grid(linestyle='--', linewidth=1)
        plt.title('a) Transport to Southern Ocean', size=12, y=1.01)
        plt.ylabel('$F_{basin,southern}$ [%]', size=12)
                
        plt.plot(depths[0:-2], f_to_southern_ocean[:,5][0:-2], label = region_names[4], marker='o', linestyle=':', c='b', markersize=5)
        plt.fill_between(depths[0:-2], f_to_southern_ocean_62[:,5][0:-2], f_to_southern_ocean_56[:,5][0:-2] , color='b', alpha=0.2)
        yerr1 = f_to_southern_ocean[:,5][-2:]-f_to_southern_ocean_62[:,5][-2:]
        yerr2 = f_to_southern_ocean_56[:,5][-2:]-f_to_southern_ocean[:,5][-2:]        
        plt.errorbar(depths[-2:], f_to_southern_ocean[:,5][-2:], yerr=[yerr1,yerr2], c='b', fmt='o', capsize=5, markersize=5)
        
        plt.plot(depths[0:-2], f_to_southern_ocean[:,4][0:-2], label = region_names[3], marker='o', linestyle=':', c='g', markersize=5)
        plt.fill_between(depths[0:-2], f_to_southern_ocean_62[:,4][0:-2], f_to_southern_ocean_56[:,4][0:-2] , color='g', alpha=0.2)                
        yerr1 = f_to_southern_ocean[:,4][-2:]-f_to_southern_ocean_62[:,4][-2:]
        yerr2 = f_to_southern_ocean_56[:,4][-2:]-f_to_southern_ocean[:,4][-2:]        
        plt.errorbar(depths[-2:], f_to_southern_ocean[:,4][-2:], yerr=[yerr1,yerr2], c='g', fmt='o', capsize=5, markersize=5)

        plt.plot(depths[0:-2], f_to_southern_ocean[:,3][0:-2], label = region_names[2], marker='o', linestyle=':', c='r', markersize=5)
        plt.fill_between(depths[0:-2], f_to_southern_ocean_62[:,3][0:-2], f_to_southern_ocean_56[:,3][0:-2] , color='r', alpha=0.2)               
        yerr1 = f_to_southern_ocean[:,3][-2:]-f_to_southern_ocean_62[:,3][-2:]
        yerr2 = f_to_southern_ocean_56[:,3][-2:]-f_to_southern_ocean[:,3][-2:]        
        plt.errorbar(depths[-2:], f_to_southern_ocean[:,3][-2:], yerr=[yerr1,yerr2], c='r', fmt='o', capsize=5, markersize=5)

        plt.grid(linestyle='--', linewidth=1)
        leg = plt.legend(prop={'size': 10})
        leg.set_title('Basin of origin', prop = {'size':12})
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.yticks(np.arange(0,10,2))        
        xlabels = np.arange(0,140,20)
        xlabels = np.append(xlabels, ['U','K'])
        plt.xticks(np.arange(0,180,20), xlabels)
        
        #F_SOSO
        plt.subplot(gs1[1])
        plt.title('b) Transport within Southern Ocean', size=12, y=1.01)
        plt.ylabel(r'$F_{southern,southern}$ [%]', size=12)
        plt.plot(depths[0:-2], f_to_southern_ocean[:,7][0:-2], marker='o', linestyle=':', c='k', markersize=5)
        plt.fill_between(depths[0:-2], f_to_southern_ocean_62[:,7][0:-2], f_to_southern_ocean_56[:,7][0:-2] , color='k', alpha=0.2)
        yerr1 = f_to_southern_ocean[:,7][-2:]-f_to_southern_ocean_62[:,7][-2:]
        yerr2 = f_to_southern_ocean_56[:,7][-2:]-f_to_southern_ocean[:,7][-2:]        
        plt.errorbar(depths[-2:], f_to_southern_ocean[:,7][-2:], yerr=[yerr1,yerr2], c='k', fmt='o', capsize=5, markersize=5)
        plt.grid(linestyle='--', linewidth=1)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.yticks(np.arange(0,90,10))
        plt.xticks(np.arange(0,180,20), xlabels)

        #Arctic F
        plt.subplot(gs1[2])
        plt.title(r'c) Transport to Arctic', size=12, y=1.01)        
        plt.ylabel(r'$F_{basin,arctic}$ [%]', size=12)
        plt.plot(depths[0:-2], f_to_arctic[:,2][0:-2], label = region_names[1], marker='o', linestyle=':', c='g', markersize=5)
        plt.fill_between(depths[0:-2], f_to_arctic_60[:,2][0:-2], f_to_arctic_70[:,2][0:-2] , color='g', alpha=0.2)        
        yerr1 = f_to_arctic_60[:,2][-2:]-f_to_arctic[:,2][-2:]
        yerr2 = f_to_arctic[:,2][-2:]-f_to_arctic_70[:,2][-2:]        
        plt.errorbar(depths[-2:], f_to_arctic[:,2][-2:], yerr=[yerr1,yerr2], c='g', fmt='o', capsize=5, markersize=5)
        plt.plot(depths[0:-2], f_to_arctic[:,1][0:-2], label = region_names[0], marker='o', linestyle=':', c='r', markersize=5)
        plt.fill_between(depths[0:-2], f_to_arctic_60[:,1][0:-2], f_to_arctic_70[:,1][0:-2] , color='r', alpha=0.2)
        yerr1 = f_to_arctic_60[:,1][-2:]-f_to_arctic[:,1][-2:]
        yerr2 = f_to_arctic[:,1][-2:]-f_to_arctic_70[:,1][-2:]
        plt.errorbar(depths[-2:], f_to_arctic[:,1][-2:], yerr=[yerr1,yerr2], c='r', fmt='o', capsize=5, markersize=5)
        plt.grid(linestyle='--', linewidth=1)
        leg = plt.legend(prop={'size': 10})
        leg.set_title('Basin of origin', prop = {'size':12})
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.yticks(np.arange(0,16,2))
        plt.xticks(np.arange(0,180,20), xlabels)
        plt.xlabel('depth [m]', size=10)

        #F_AA
        plt.subplot(gs1[3])
        plt.title('d) Transport within Arctic', size=12, y=1.01)
        plt.ylabel(r' $F_{arctic,arctic}$ [%]', size=12)
        plt.plot(depths[0:-2], f_to_arctic[:,6][0:-2], label = region_names[5], marker='o', linestyle=':', c='k', markersize=5)
        plt.fill_between(depths[0:-2], f_to_arctic_60[:,6][0:-2], f_to_arctic_70[:,6][0:-2] , color='k', alpha=0.2)
        yerr1 = f_to_arctic_60[:,6][-2:]-f_to_arctic[:,6][-2:]
        yerr2 = f_to_arctic[:,6][-2:]-f_to_arctic_70[:,6][-2:]
        plt.errorbar(depths[-2:], f_to_arctic[:,6][-2:], yerr=[yerr1,yerr2], c='k', fmt='o', capsize=5, markersize=5)        
        plt.grid(linestyle='--', linewidth=1)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.yticks(np.arange(40,90,10))
        plt.xticks(np.arange(0,180,20), xlabels)
        plt.xlabel('depth [m]', size=10)

        plt.savefig(outdir_paper + 'regional_transport_to_poles.pdf', dpi=900, bbox_inches='tight')
        
    polward_transport()

FIG5_RegionalTransport()


"""
FIGURES SUPPLEMENTARY
"""
        
        
def S1_surface_distribution_different_times():
    """
    Surface distributions at different times to see that distribution comes to an equilibrium
    """
    
    #load data at several times
    fname = 'SubSurf_y2000_m1_d5_simdays3650_layer0_pos' #File names, apart from the label of initial grid
    tload = range(0,730,73)+[729] #Times when data is loaded (indices, not actual times)
    SurfaceData = ParticleData.from_nc(datadir + 'Layer0/', fname, tload=tload, Ngrids=40) #Load data

    #for the figure
    fig = plt.figure(figsize = (14,21)) 
    gs1 = gridspec.GridSpec(6, 2)
    gs1.update(wspace=0.01, hspace=0.0)
    
    for i in range(len(tload)):

        #get distribution at specific time
        d=SurfaceData.get_distribution(i, 2.) #Compute distribution
        d=oceanvector(d,ddeg=2.) #Create oceanvector object from it for plotting
        lon_bins_2d,lat_bins_2d = np.meshgrid(d.Lons_edges,d.Lats_edges)

        #for the plot        
        plt.subplot(gs1[i])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.drawcoastlines()
        m.fillcontinents(color='dimgrey')
        xs, ys = m(lon_bins_2d, lat_bins_2d) 
        
        #plot distribution
        plt.pcolormesh(xs, ys, d.V2d,cmap='plasma', norm=colors.LogNorm(), rasterized=True)
        cbar=plt.colorbar(orientation='vertical',shrink=0.5)
        cbar.ax.tick_params(labelsize=7, width=0.03)
        plt.title('Year ' + str(int(round(tload[i]*5/365.,0))), size=10)
    
    fig.savefig(outdir_supplementary + 'distribution_surface_in_time.pdf', dpi=900, bbox_inches='tight')
    
S1_surface_distribution_different_times()
    
def S2_plot_regions():

    r, region_names = regions(2., so_lat=-60, na_lat=65)
    r = np.array([r[i] if r[i] != 0. else 1 for i in range(len(r))])
    r=oceanvector(r, ddeg=2.)
    
    fig = plt.figure(figsize = (12,8))
    m = Basemap(projection='robin',lon_0=180,resolution='l')
    m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=1.2, size=10)
    m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=1.2, size=10)
    m.drawcoastlines()
    m.fillcontinents(color='dimgrey')
    lon_bins_2d,lat_bins_2d = np.meshgrid(r.Lons_edges,r.Lats_edges)
    xs, ys = m(lon_bins_2d, lat_bins_2d) 
    cmap = plt.get_cmap('jet', 7)
    plt.pcolormesh(xs, ys, r.V2d, cmap=cmap, rasterized=True)
    cbar=plt.colorbar(orientation='vertical', shrink=0.6)
    cbar.ax.tick_params(labelsize=10)
    tick_locs = 1+3./7 * np.array(range(1,16,2))
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(region_names)
    
    fig.savefig(outdir_supplementary + 'regions.pdf', dpi=900, bbox_inches='tight')

S2_plot_regions()


def S3_distributions_different_models():   

    filenames = ['SubSurf_y2000_m1_d5_simdays3650_layer2_pos','SubSurf_y2000_m1_d5_simdays3650_layer4_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer7_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer13_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer22_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer23_pos', 
              'SubSurf_y2000_m1_d5_simdays3650_layer25_pos','Kukulka_y2000_m1_d5_simdays3650_layer0_pos']
    folders = ['Layer2', 'Layer4', 'Layer7','Layer13','Layer22', 'Layer23', 'Layer25','KukMix']    
    folders = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders]           
    labels = ['a) z = ' + str(nemo_depth[2]) + ' m', 
    'b) z = ' + str(nemo_depth[4]) + ' m', 
    'c) z = ' + str(nemo_depth[7]) + ' m', 
    'd) z = ' + str(nemo_depth[13]) + ' m', 
    'e) z = ' + str(nemo_depth[22]) + ' m', 
    'f) z = ' + str(nemo_depth[23]) + ' m', 
    'g) z = ' + str(nemo_depth[25]) + ' m', 
    'h) Kukulka mixing'] #figure titles    
    
    #for the figure
    fig = plt.figure(figsize = (14,14)) 
    gs1 = gridspec.GridSpec(4, 2)
    gs1.update(wspace=0.01, hspace=0.0)

    for i in range(len(folders)):
        
        #load particle data
        filename = filenames[i]
        pdata=ParticleData.from_nc(folders[i] + '/', filename, tload=[0,-1], Ngrids=40)
               
        #create subplot
        plt.subplot(gs1[i])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.drawcoastlines()
        m.fillcontinents(color='dimgrey')
        
        #get distribution
        d_full=pdata.get_distribution(t=-1, ddeg=2.).flatten()
        d=oceanvector(d_full, ddeg = 2.)        
        lon_bins_2d,lat_bins_2d = np.meshgrid(d.Lons_edges,d.Lats_edges)
        xs, ys = m(lon_bins_2d, lat_bins_2d) 
        
        #plot distribution
        plt.pcolormesh(xs, ys, d.V2d,cmap='plasma', norm=colors.LogNorm(), rasterized=True)
        cbar=plt.colorbar(orientation='vertical',shrink=0.5)
        cbar.ax.tick_params(labelsize=7, width=0.03)
        plt.title(labels[i], size=10, y=1.)
    
        fig.savefig(outdir_supplementary + 'distributions_other_scenarios.pdf', dpi=900, bbox_inches='tight')

S3_distributions_different_models()

        
def S4_basintransport_other():   
    
    filenames = ['SubSurf_y2000_m1_d5_simdays3650_layer2_pos','SubSurf_y2000_m1_d5_simdays3650_layer4_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer7_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer10_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer13_pos','SubSurf_y2000_m1_d5_simdays3650_layer19_pos', 
              'SubSurf_y2000_m1_d5_simdays3650_layer22_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer23_pos',
              'Kukulka_y2000_m1_d5_simdays3650_layer0_pos']
    folders = ['Layer2', 'Layer4' , 'Layer7', 'Layer10','Layer13','Layer19', 'Layer22', 'Layer23','KukMix']
    folders = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders]
    labels = ['a) z = ' + str(nemo_depth[2]) + ' m', 
    'b) z = ' + str(nemo_depth[4]) + ' m', 
    'c) z = ' + str(nemo_depth[7]) + ' m', 
    'd) z = ' + str(nemo_depth[10]) + ' m', 
    'e) z = ' + str(nemo_depth[13]) + ' m', 
    'f) z = ' + str(nemo_depth[19]) + ' m', 
    'g) z = ' + str(nemo_depth[22]) + ' m', 
    'h) z = ' + str(nemo_depth[23]) + ' m', 
    'i) Kukulka mixing'] #figure titles    
    
    #for the figure
    fig = plt.figure(figsize = (14,18)) 
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(wspace=0.1   , hspace=0.)

    #get region names and constrain to the subtropical basins
    _, region_names = regions(2.)
    region_names=region_names[0:5]
    basin_regions = [1,2,3,4,5]

    for j in range(len(folders)):
        
        #load particle data        
        filename = filenames[j]
        pdata=ParticleData.from_nc(folders[j] + '/', filename, tload=[0,-1], Ngrids=40)
        pdata.remove_nans()
        region_label=pdata.set_region_labels(2., -1)
        
        #select the particles in the basins, and label the particles according to their final basin
        selection = [k for k in range(len(region_label)) if region_label[k] in basin_regions]
        lons=pdata.lons[selection]
        lats=pdata.lats[selection]
        region_label = [region_label[i] for i in range(len(region_label)) if region_label[i] in basin_regions]
        
        #plot
        plt.subplot(gs1[j])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.fillcontinents(color='dimgray')
        xs, ys = m(lons[:,0], lats[:,0])
        cmap = plt.get_cmap('jet', 5)
        scat=m.scatter(xs, ys, marker='.', s=1, c=region_label, cmap=cmap, rasterized=True)
        plt.title(labels[j], size=10, y=1.0)
    
    #color bar on the right
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.822, 0.43, 0.015, 0.3])
    cbar=fig.colorbar(scat, cax=cbar_ax, ticks=[1,2,3,4,5])
    cbar.ax.tick_params(labelsize=10)
    tick_locs = 1+0.4 * np.array([1,3,5,7,9])
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(region_names)

    #save figure
    fig.savefig(outdir_supplementary + 'transport_basins_scatter_other.pdf', dpi=900, bbox_inches='tight')
    plt.close()

S4_basintransport_other()


def S5_polattransport_other():   
    
    filenames = ['SubSurf_y2000_m1_d5_simdays3650_layer2_pos','SubSurf_y2000_m1_d5_simdays3650_layer4_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer7_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer10_pos',
              'SubSurf_y2000_m1_d5_simdays3650_layer13_pos','SubSurf_y2000_m1_d5_simdays3650_layer19_pos', 
              'SubSurf_y2000_m1_d5_simdays3650_layer22_pos', 'SubSurf_y2000_m1_d5_simdays3650_layer23_pos',
              'Kukulka_y2000_m1_d5_simdays3650_layer0_pos']
    folders = ['Layer2', 'Layer4' , 'Layer7', 'Layer10','Layer13','Layer19', 'Layer22', 'Layer23','KukMix']
    folders = ['/Users/wichmann/Simulations/Proj1_SubSurface_Mixing/' + f for f in folders]        
    labels = ['a) z = ' + str(nemo_depth[2]) + ' m', 
    'b) z = ' + str(nemo_depth[4]) + ' m', 
    'c) z = ' + str(nemo_depth[7]) + ' m', 
    'd) z = ' + str(nemo_depth[10]) + ' m', 
    'e) z = ' + str(nemo_depth[13]) + ' m', 
    'f) z = ' + str(nemo_depth[19]) + ' m', 
    'g) z = ' + str(nemo_depth[22]) + ' m', 
    'h) z = ' + str(nemo_depth[23]) + ' m', 
    'i) Kukulka mixing'] #figure titles        
    
    #for the figure
    fig = plt.figure(figsize = (14,18)) 
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(wspace=0.1   , hspace=0.)
    
    #get region names and constrain to the poles
    _, region_names = regions(2.)
    region_names=region_names[5:]
    basin_regions = [6,7]

    for j in range(len(folders)):
        
        #load particle data        
        filename = filenames[j]
        pdata=ParticleData.from_nc(folders[j] + '/', filename, tload=[0,-1], Ngrids=40)
        pdata.remove_nans()
        region_label=pdata.set_region_labels(2., -1)
        
        #select the particles in the basins, and label the particles according to their final basin
        selection = [k for k in range(len(region_label)) if region_label[k] in basin_regions]
        lons=pdata.lons[selection]
        lats=pdata.lats[selection]
        region_label = [region_label[i] for i in range(len(region_label)) if region_label[i] in basin_regions]
        
        #plot
        plt.subplot(gs1[j])
        m = Basemap(projection='robin',lon_0=180,resolution='l')
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.drawparallels([-60,-30,0,30,60], labels=[True, False, False, True], linewidth=.8, size=7)
        m.drawmeridians([50,150,250,350], labels=[False, False, False, True], linewidth=.8, size=7)
        m.fillcontinents(color='dimgray')
        xs, ys = m(lons[:,0], lats[:,0])
        cmap = plt.get_cmap('jet', 2)
        scat=m.scatter(xs, ys, marker='.', s=1, c=region_label, cmap=cmap, rasterized=True)
        plt.title(labels[j], size=10, y=1.0)
    
    #color bar on the right
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.822, 0.5, 0.015, 0.15])
    cbar=fig.colorbar(scat, cax=cbar_ax, ticks=[6,7])
    cbar.ax.tick_params(labelsize=10)
    tick_locs = 6+0.25 * np.array([1,3])
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(region_names)

    #save figure
    fig.savefig(outdir_supplementary + 'transport_poles_scatter_other.pdf', dpi=900, bbox_inches='tight')
    plt.close()

S5_polattransport_other()


"""
TABLES SUPPLEMENTARY
"""

def T1_13_create_transport_tables():
    
    labels = ['z = ' + str(nemo_depth[0]) + ' m',
              'z = ' + str(nemo_depth[2]) + ' m',
              'z = ' + str(nemo_depth[4]) + ' m',
              'z = ' + str(nemo_depth[7]) + ' m',
              'z = ' + str(nemo_depth[10]) + ' m',
              'z = ' + str(nemo_depth[13]) + ' m',
              'z = ' + str(nemo_depth[16]) + ' m',
              'z = ' + str(nemo_depth[19]) + ' m',
              'z = ' + str(nemo_depth[22]) + ' m',
              'z = ' + str(nemo_depth[23]) + ' m',
              'z = ' + str(nemo_depth[25]) + ' m',
              'Random Uniform',
              'Kukulka mixing']
              
    #Load n_{ij}. If it does not exist, first run create_transport_matrix in FIG5_RegionalTransport
    n_matrix = np.load(outdir_paper + 'n_matrix_solat_-60_na_lat_65.npy').tolist()    
    region_names = ['', 'NP', 'NA', 'SP', 'SA', 'IO', 'A', 'SO']
        
    for i in range(len(models_all)):
        n=n_matrix[models_all[i]]

        #row-normalize the matrix to get F
        s=np.sum(n,axis=1)
        f=np.zeros(n.shape)            
        for j in range(len(n)):
            f[j:,]=n[j,:]/s[j]
        f*=100 
        
        #write in latex format
        with open('./supplementary_material/tables/F_' + models_all[i] + '.txt', "w") as output:
            output.write('\\caption{F-matrix [\%]: ' + labels[i] + '}' + '\n')
            output.write('\\centering' + '\n')
            output.write('\\begin{tabular}{| l | c | c | c | c | c | c | c |}' + '\n')
            output.write('\\hline' + '\n')
            
            for k in range(len(f)):
                if k<len(f)-1:
                    output.write(str(region_names[k]) + ' & ')
                else:
                    output.write(str(region_names[k]) + ' \\\\' + '\n')
            
            output.write('\\hline' + '\n')
            for k in range(1,len(f)):
                output.write(str(region_names[k]) + ' & ')
                for l in range(1,len(f)):
                    if l<len(f)-1:                    
                        output.write(str(round(f[k][l],1)) + ' & ')
                    else:
                        output.write(str(round(f[k][l],1)) + ' \\\\' + '\n')

T1_13_create_transport_tables()