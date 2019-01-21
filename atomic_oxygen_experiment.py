import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec
import math
from matplotlib import colors

deg = unichr(176)
delta = unichr(916)
k_B = 1.38064852e-23
N_A = 6.02214086e+23
R = 8.3144598
species_list = ['atomic_oxygen', 'ozone', 'atomic_hydrogen', 'temperature', 'density']
symbol_list = ['O', 'O3', 'H', 'T', 'n']
units_list = ['ppmv', '$\mathregular{cm^{-3}}$', 'K']

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/zonal_means/f.e20.FXSD.f19_f19.001.cam.h0.2000-01.nc', 'r', format='NETCDF4')
lats = fname_uni.variables['lat'][:]
fname_uni.close()

def calc_z3_zon_mer_t_av(levs):  
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

def calc_z3_zon_t_av(levs):  
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    z3 = np.zeros([1,levs,96,144])
    z3[:] = fname.variables['Z3'][:]*(1.e-3)
    z3_zon_av = np.mean(z3[:], axis=3)
    z3_zon_t_av = np.mean(z3_zon_av[:], axis=0) 
    fname.close()
    return z3_zon_t_av

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

def calc_species_zon_av(symbol, levs):  
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4') 
    species_dat = np.zeros([1,levs,96,144])
    if symbol == 'T':
        species_dat = fname.variables[symbol][:]
    elif symbol == 'n':
        species_dat = fname.variables['T'][:]
    else:
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

def calc_z3_zon_t_av_weighted(levs, lowlat, highlat):
    z3_zon_t_av = calc_z3_zon_t_av(levs)
    z3_zon_t_av_weighted = calc_cos_factor(z3_zon_t_av, levs, lowlat, highlat)
    return z3_zon_t_av_weighted

def get_lev(levs):
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4') 
    lev = np.zeros([levs])
    lev = fname.variables['lev'][:]
    fname.close()
    return lev

def calc_density(T, levs, lowlat, highlat):
    lev = get_lev(levs)
    T_weighted = calc_cos_factor(T, levs, lowlat, highlat)
    n_weighted = np.zeros(levs)
    for i in range(0,levs):
        n_weighted[i] = (N_A * 100 * lev[i]) / (R * T_weighted[i]) * (1.e-6)
    return n_weighted

def calc_profiles(param, levs, lowlat, highlat):
    if symbol == 'n':
        param_weighted = calc_density(param, levs, lowlat, highlat)
    else:
        param_weighted = calc_cos_factor(param, levs, lowlat, highlat)
    return param_weighted

def calc_conc_profiles(param, levs, lowlat, highlat):
    lev = get_lev(levs)
    T_zon_t_av = calc_species_zon_av('T', levs)
    T_zon_t_av_weighted = calc_cos_factor(T_zon_t_av, levs, lowlat, highlat)    
    param_weighted = calc_cos_factor(param, levs, lowlat, highlat)
    param_weighted_conc = np.zeros(levs)  
    for i in range(0,levs):
        param_weighted_conc[i] = (param_weighted[i] * 1.e-6 * N_A * 100 * lev[i]) / (R * T_zon_t_av_weighted[i]) * (1.e-6)
    return param_weighted_conc

def plot_1d_global(name, config, units, z3, species, color, plot_no):
    x = species[::-1]
    y = z3[::-1]
    plt.plot(x, y, color=color, label=config)
    plt.xlabel('%s [%s]' %(name, units), fontsize=12)
    plt.ylabel('Altitude [km]', fontsize=12)
    plt.ylim(60,160)
    if name == 'atomic_oxygen':
        if units == 'ppmv':
            plt.xlim(0,500000)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,8.e+11)
    if name == 'ozone':
        plt.xscale('log')
        if units == 'ppmv':
            plt.xlim(1.e-8,1.e+1)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(1.e-4,1.e+12)
    if name == 'atomic_hydrogen':
        if units == 'ppmv':
            plt.xlim(0,20)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,5.e+8)
    if name == 'temperature':
        1==1
    if name == 'density':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ylim(100,160)
        plt.xlim(0,3.e+12)
    if plot_no == 1:
        plt.legend()
    return

def plot_1d_multi(name, config, units, z3, species, lowlat, highlat, color, plot_no):
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = species[::-1]
    y = z3[::-1]
    plt.plot(x, y, color=color, label=config)
    plt.ylim(60,160)
    if plot_no == 0:
        plt.ylabel('Altitude [km]', fontsize=12)
        plt.tick_params(labelbottom='off')
    if plot_no == 1:
        plt.tick_params(labelleft='off')
        plt.tick_params(labelbottom='off')
    if plot_no == 2:
        plt.tick_params(labelleft='off')
        plt.tick_params(labelbottom='off')
    if plot_no == 3:
        plt.xlabel('%s [%s]' %(name, units), fontsize=12)
        plt.ylabel('Altitude [km]', fontsize=12)
    if plot_no == 4:
        plt.xlabel('%s [%s]' %(name, units), fontsize=12)
        plt.tick_params(labelleft='off')
    if plot_no == 5:
        plt.xlabel('%s [%s]' %(name, units), fontsize=12)
        plt.tick_params(labelleft='off')
    if name == 'atomic_oxygen':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        if units == 'ppmv':
            plt.xlim(0,500000)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,8.e+11)
    if name == 'ozone':
        plt.xscale('log')
        if units == 'ppmv':
            plt.xlim(1.e-8,1.e+1)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(1.e-4,1.e+12)
    if name == 'atomic_hydrogen':
        if units == 'ppmv':
            plt.xlim(0,20)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,5.e+8)
    if name == 'temperature':
        1==1
    if name == 'density':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ylim(100,160)
        plt.xlim(0,3.e+12)
    if config == 'waccm-x' and plot_no == 2:
        plt.legend(loc=1)
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
        cbar_ticks = np.arange(0,100,10)
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

year = 2014
month = 1
name = species_list[4]
symbol = symbol_list[4]
units = units_list[1]
chemistry = False
global_only = True
save = False

if units == 'ppmv':
    units_print = 'ppmv'
elif units == '$\mathregular{cm^{-3}}$':
    units_print = 'cm-3'

if global_only == True:
    if month == 1:
        plt.title('January %s Global' %year , fontsize=16)
    elif month == 7:
        plt.title('July %s Global' %year, fontsize=16)
    step = 96
    a = 0
    b = 1
else:
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
    gs1 = gridspec.GridSpec(2, 3)
    gs1.update(wspace=0.1, hspace=0.1)
    if month == 1:
        plt.suptitle('January %s' %year, fontsize=16)
    elif month == 7:
        plt.suptitle('July %s' %year, fontsize=16)
    step = 16
    a = 0
    b = 6

#waccm_z3 = calc_z3_zon_mer_t_av(88)
#waccmx_z3 = calc_z3_zon_mer_t_av(145)

waccm_species = calc_species_zon_av(symbol, 88)
waccmx_species = calc_species_zon_av(symbol, 145)

#waccmx_species_int = interp_waccmx_species(waccm_z3, waccmx_z3, waccmx_species)
#diff = calc_diff(waccm_species, waccmx_species_int)

# 1D Plot Code
for i in range(a,b):    
    lowlat = i * step
    highlat = (i * step) + step
    lowlat_no = int((lowlat * 1.875) - 90)
    highlat_no = int((highlat * 1.875) - 90)
    waccm_z3_weighted = calc_z3_zon_t_av_weighted(88, lowlat, highlat)
    waccmx_z3_weighted = calc_z3_zon_t_av_weighted(145, lowlat, highlat)    
    if chemistry == True:
        if units == 'ppmv':
            waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat)
        elif units == '$\mathregular{cm^{-3}}$':
            waccm_species_profile = calc_conc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_conc_profiles(waccmx_species, 145, lowlat, highlat)
    else:
        if symbol == 'T':
            waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat) 
        elif symbol == 'n':
            waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat)        
    if global_only == True:
        plot_1d_global(name, 'waccm', units, waccm_z3_weighted, waccm_species_profile, 'k', 0)
        plot_1d_global(name, 'waccm-x', units, waccmx_z3_weighted, waccmx_species_profile, 'b', 1)
    else:
        plot_1d_multi(name, 'waccm', units, waccm_z3_weighted, waccm_species_profile, lowlat, highlat, 'k', i)
        plot_1d_multi(name, 'waccm-x', units, waccmx_z3_weighted, waccmx_species_profile, lowlat, highlat, 'b', i)
if save == True:
    if global_only == True:
        plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s/%s_month%s_profile_global_%s.jpg' %(year, name, month, units_print), bbox_inches='tight', dpi=300)
    else:
        plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s/%s_month%s_profile_lat_bands_%s.jpg' %(year, name, month, units_print), bbox_inches='tight', dpi=300)
'''
# 2D Plot Code
# Currently for chemistry only
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(11,5))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(wspace=0.1, hspace=0.1)

if month == 1:
    fig.suptitle('January %s' %year, fontsize=16)
elif month == 7:
    fig.suptitle('July %s' %year, fontsize=16)

plot_2d(name, waccm_z3, waccm_species, 0)
plot_2d(name, waccmx_z3, waccmx_species, 1)
plot_2d(name, waccm_z3, diff, 2)
plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s/%s_month%s.jpg' %(year, name, month), bbox_inches='tight', dpi=300)
'''
plt.show()