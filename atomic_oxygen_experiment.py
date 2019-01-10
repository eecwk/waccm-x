import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec
import math
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
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    z3 = np.zeros([1,levs,96,144])
    z3[:] = fname.variables['Z3'][:]*(1.e-3)
    z3_zon_av = np.mean(z3[:], axis=3)
    z3_zon_mer_av = np.mean(z3_zon_av[:], axis=2)
    z3_zon_mer_t_av = np.mean(z3_zon_mer_av[:], axis=0) 
    fname.close()
    return z3_zon_mer_t_av

def calc_species_zon_av(month, symbol, levs):  
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4') 
    species_dat = np.zeros([1,levs,96,144]) 
    species_dat = fname.variables[symbol][:]*(1.e6)
    species_tm = np.mean(species_dat[:], axis=0)
    species_zon_av = np.mean(species_tm[:], axis=2)
    fname.close()
    return species_zon_av

def interp_waccmx_species(z3_1, z3_2, species_2):
    species_2_int = np.zeros([88,96])
    species_2_int_rev = np.zeros([88,96])
    z3_1_rev = z3_1[::-1]
    z3_2_rev = z3_2[::-1]
    species_2_rev = species_2[::-1]
    for i in range(0,88):  
        for j in range(0,96):
            species_2_int[i,j] = np.interp(z3_1_rev[i], z3_2_rev[:], species_2_rev[:,j]) 
            species_2_int_rev = species_2_int[::-1]
    return species_2_int_rev

def calc_diff(param1, param2): 
    diff = np.zeros([88,96])
    for i in range(0,88):
        for j in range(0,96): 
            diff[i,j] = ( (param2[i,j] - param1[i,j]) / param1[i,j] ) * 100
    return diff

def calc_cos_factor(param, levs, lowlat, highlat):
    param_weighted = np.zeros(levs)    
    for j in range (0, levs):    
        sig_cos_x = 0
        sig_cos = 0
        for k in range (lowlat, highlat):
            sig_cos_x = sig_cos_x + (math.cos(math.radians(lats[k])) * param[j][k])
            sig_cos = sig_cos + math.cos(math.radians(lats[k]))         
            if  k == (highlat - 1):
                param_weighted[j] = sig_cos_x / sig_cos
    return param_weighted

def calc_profiles(param, levs, lowlat, highlat):
    param_weighted = calc_cos_factor(param, levs, lowlat, highlat)
    return param_weighted
    
def plot_1d(name, config, z3, species, lowlat, highlat, color, plot_no):
    if plot_no > 2:
        plot_no = plot_no - 3
    plt.subplot(gs1[plot_no])
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = species[::-1]
    y = z3[::-1]
    plt.plot(x, y, color=color, label=config)
    plt.xlabel('%s [ppmv]' %name, fontsize=12)
    plt.ylim(90,150)
    if plot_no == 0:
        plt.ylabel('Altitude [km]', fontsize=12)
    if plot_no == 1:
        plt.tick_params(labelleft='off')
    if plot_no == 2:
        plt.tick_params(labelleft='off')
    if name == 'atomic_oxygen':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlim(0,500000)
    if name == 'ozone':
        plt.xscale('log')
    if name == 'atomic_hydrogen':
        plt.xlim(0,30)
    if config == 'waccm-x' and plot_no == 2:
        plt.legend(loc=4)
    return
    

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
        diffs = np.arange(1,91,1)
        cbar_ticks = np.arange(10,100,10)
        plot = 'linear'
    diffs_per = np.arange(-200,201,1)
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
        ax2 = plt.contourf(x[:,:], y[:,:], species[:,:], diffs_per, extend='both', cmap=plt.get_cmap('seismic'))
        plt.title('WACCM to WACCM-X Difference')
        plt.tick_params(labelleft='off')
        cbar_ax2 = fig.add_axes([1.05, 0.15, 0.02, 0.7])
        cbar2 = fig.colorbar(ax2, cax=cbar_ax2, ticks=np.arange(-200,250,50), orientation='vertical')
        cbar2.cmap.set_under('#001648')
        cbar2.set_label('[%]', fontsize=12)
        cbar2.ax.tick_params(labelsize=12)
    return

year = 2009
month = 7
name = species_list[2]
symbol = symbol_list[2]

waccm_z3 = calc_z3_zon_av(month, symbol, 88)
waccmx_z3 = calc_z3_zon_av(month, symbol, 145)
waccm_species = calc_species_zon_av(month, symbol, 88)
waccmx_species = calc_species_zon_av(month, symbol, 145)
waccmx_species_int = interp_waccmx_species(waccm_z3, waccmx_z3, waccmx_species)
diff = calc_diff(waccm_species, waccmx_species_int)

# 1D Plot Code
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(11,5))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(wspace=0.1, hspace=0.1)

if month == 1:
    plt.suptitle('January %s' %year, fontsize=16)
elif month == 7:
    plt.suptitle('July %s' %year, fontsize=16)

step = 16

for i in range(3,6):
    
    lowlat = i * step
    highlat = (i * step) + step
    lowlat_no = int((lowlat * 1.875) - 90)
    highlat_no = int((highlat * 1.875) - 90)

    waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
    waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat)
    plot_1d(name, 'waccm', waccm_z3, waccm_species_profile, lowlat, highlat, 'k', i)
    plot_1d(name, 'waccm-x', waccmx_z3, waccmx_species_profile, lowlat, highlat, 'b', i)
#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/solarmin_2009/%s_month%s_profile_SH_bands.jpg' %(name, month), bbox_inches='tight', dpi=300)
plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/solarmin_2009/%s_month%s_profile_NH_bands.jpg' %(name, month), bbox_inches='tight', dpi=300)

# 2D Plot Code
#fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(11,5))
#gs1 = gridspec.GridSpec(1, 3)
#gs1.update(wspace=0.1, hspace=0.1)

#if month == 1:
#    fig.suptitle('January 2004', fontsize=16)
#elif month == 7:
#    fig.suptitle('July 2004', fontsize=16)

#plot_2d(name, waccm_z3, waccm_species, 0)
#plot_2d(name, waccmx_z3, waccmx_species, 1)
#plot_2d(name, waccm_z3, diff, 2)
#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s_month%s.jpg' %(name, month), bbox_inches='tight', dpi=300)

plt.show()