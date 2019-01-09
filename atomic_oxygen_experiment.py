import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import colors

deg = unichr(176)
delta = unichr(916)
k_B = 1.38064852e-23
species_list = ['atomic_oxygen', 'ozone', 'atomic_hydrogen']
symbol_list = ['O', 'O3', 'H']

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/zonal_means/f.e20.FXSD.f19_f19.001.cam.h0.2000-01.nc', 'r', format='NETCDF4')
lats = fname_uni.variables['lat'][:]
fname_uni.close()

def calc_z3_zon_av(month, symbol, levs):  
    if levs == 70:
        fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/Tao_Na_Paper_FW5_JCMIFd5.cam.h0.0001-0%s.nc' %month, 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2001-0%s.nc' %month, 'r', format='NETCDF4')
    z3 = np.zeros([1,levs,96,144])
    z3[:] = fname.variables['Z3'][:]*(1.e-3)
    z3_zon_av = np.mean(z3[:], axis=3)
    z3_zon_mer_av = np.mean(z3_zon_av[:], axis=2)
    z3_zon_mer_t_av = np.mean(z3_zon_mer_av[:], axis=0) 
    fname.close()
    return z3_zon_mer_t_av

def calc_species_zon_av(month, symbol, levs):  
    if levs == 70:
        fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/Tao_Na_Paper_FW5_JCMIFd5.cam.h0.0001-0%s.nc' %month, 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2001-0%s.nc' %month, 'r', format='NETCDF4') 
    species_dat = np.zeros([1,levs,96,144]) 
    species_dat = fname.variables[symbol][:]*(1.e6)
    species_tm = np.mean(species_dat[:], axis=0)
    species_zon_av = np.mean(species_tm[:], axis=2)
    fname.close()
    return species_zon_av

def interp_waccmx_species(z3_1, z3_2, species_2):
    species_2_int = np.zeros([70,96])
    species_2_int_rev = np.zeros([70,96])
    z3_1_rev = z3_1[::-1]
    z3_2_rev = z3_2[::-1]
    species_2_rev = species_2[::-1]
    for i in range(0,70):  
        for j in range(0,96):
            species_2_int[i,j] = np.interp(z3_1_rev[i], z3_2_rev[:], species_2_rev[:,j]) 
            species_2_int_rev = species_2_int[::-1]
    return species_2_int_rev

def calc_diff(param1, param2): 
    diff = np.zeros([70,96])
    for i in range(0, 70):
        for j in range(0, 96): 
            diff[i,j] = ( (param2[i,j] - param1[i,j]) / param1[i,j] ) * 100
            #if diff[i,j] > 500:
                #diff[i,j] = 500
    return diff

def plot_2d(name, z3, species, plot_no):
    plt.subplot(gs1[plot_no])
    x, y = np.meshgrid(lats, z3)
    plt.xlabel('Latitude [%s]' %deg, fontsize=12)
    plt.xticks(np.arange(-90,120,30), fontsize=12) 
    plt.yticks(np.arange(0,220,20), fontsize=12)   
    plt.ylim(90,200)
    plt.axhline(y=waccm_z3[0], color='w', linewidth=1, linestyle=':')
    if name == 'atomic_oxygen':
        diffs = [1.e+3, 175.e+1, 25.e+2, 325.e+1, 4.e+3, 475.e+1, 55.e+2, 625.e+1, 7.e+3, 775.e+1, 85.e+2, 925.e+1, 1.e+4, 175.e+2, 25.e+3, 325.e+2, 4.e+4, 475.e+2, 55.e+3, 625.e+2, 7.e+4, 775.e+2, 85.e+3, 925.e+2, 1.e+5, 175.e+3, 25.e+4, 325.e+3, 4.e+5, 475.e+3, 55.e+4, 625.e+3, 7.e+5, 775.e+3, 85.e+4, 925.e+3, 1.e+6]
        cbar_ticks = [1.e+3, 1.e+4, 1.e+5, 1.e+6]
        plot = 'log'
    elif name == 'ozone':
        diffs = [1.e-8, 25.e-9, 4.e-8, 55.e-9, 7.e-8, 85.e-9, 1.e-7, 25.e-8, 4.e-7, 55.e-8, 7.e-7, 85.e-8, 1.e-6, 25.e-7, 4.e-6, 55.e-7, 7.e-6, 85.e-7, 1.e-5, 25.e-6, 4.e-5, 55.e-6, 7.e-5, 85.e-6, 1.e-4, 25.e-5, 4.e-4, 55.e-5, 7.e-4, 85.e-5, 1.e-3, 25.e-4, 4.e-3, 55.e-4, 7.e-3, 85.e-4, 1.e-2, 25.e-3, 4.e-2, 55.e-3, 7.e-2, 85.e-3, 1.e-1, 25.e-2, 4.e-1, 55.e-2, 7.e-1, 85.e-2, 1.e+0, 25.e-1, 4.e+0, 55.e-1, 7.e+0, 85.e-1, 1.e+1]
        cbar_ticks = [1.e-8, 1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.e+0, 1.e+1]
        plot = 'log'
    elif name == 'atomic_hydrogen':
        diffs = np.arange(1,61,1)
        cbar_ticks = np.arange(5,65,5)
        plot = 'linear'
    diffs_per = np.arange(-400,401,1)
    if plot_no == 0:
        if plot == 'linear':
            ax = plt.contourf(x[:,:], y[:,:], species[:,:], diffs)
        elif plot == 'log':
            ax = plt.contourf(x[:,:], y[:,:], species[:,:], diffs, norm=colors.LogNorm())
        plt.title('WACCM')
        plt.ylabel('Altitude [km]', fontsize=12)
    elif plot_no == 1:
        if plot == 'linear':
            ax = plt.contourf(x[:,:], y[:,:], species[:,:], diffs)
        elif plot == 'log':
            ax = plt.contourf(x[:,:], y[:,:], species[:,:], diffs, norm=colors.LogNorm())
        plt.title('WACCM-X')
        plt.tick_params(labelleft='off')
        cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(ax, cax=cbar_ax, ticks=cbar_ticks, orientation='vertical')
        cbar.set_label('%s [ppmv]' %name, fontsize=12)
        cbar.ax.tick_params(labelsize=12)
    elif plot_no == 2:
        ax2 = plt.contourf(x[:,:], y[:,:], species[:,:], diffs_per, cmap=plt.get_cmap('seismic'))
        plt.title('WACCM to WACCM-X Difference')
        plt.tick_params(labelleft='off')
        cbar_ax2 = fig.add_axes([1.05, 0.15, 0.02, 0.7])
        #cbar2 = fig.colorbar(ax2, cax=cbar_ax2, ticks=np.arange(-200,250,50), orientation='vertical')
        cbar2 = fig.colorbar(ax2, cax=cbar_ax2, ticks=np.arange(-400,500,100), orientation='vertical')
        cbar2.set_label('[%]', fontsize=12)
        cbar2.ax.tick_params(labelsize=12)
    return

month = 7
name = species_list[2]
symbol = symbol_list[2]

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(11,5))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(wspace=0.1, hspace=0.1)

if month == 1:
    fig.suptitle('January', fontsize=16)
elif month == 7:
    fig.suptitle('July', fontsize=16)

waccm_z3 = calc_z3_zon_av(month, symbol, 70)
waccmx_z3 = calc_z3_zon_av(month, symbol, 145)
waccm_species = calc_species_zon_av(month, symbol, 70)
waccmx_species = calc_species_zon_av(month, symbol, 145)
waccmx_species_int = interp_waccmx_species(waccm_z3, waccmx_z3, waccmx_species)
diff = calc_diff(waccm_species, waccmx_species_int)
plot_2d(name, waccm_z3, waccm_species, 0)
plot_2d(name, waccmx_z3, waccmx_species, 1)
plot_2d(name, waccm_z3, diff, 2)

#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/%s_month%s.jpg' %(name, month), bbox_inches='tight', dpi=300)
plt.show()