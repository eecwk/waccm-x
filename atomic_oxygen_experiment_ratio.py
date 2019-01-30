import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec
import math

deg = unichr(176)
delta = unichr(916)
k_B = 1.38064852e-23
N_A = 6.02214086e+23
R = 8.3144598
species_list = ['atomic_oxygen', 'ozone', 'atomic_hydrogen', 'carbon_dioxide', 'carbon_monoxide', 'temperature', 'density']
symbol_list = ['O', 'O3', 'H', 'CO2', 'CO', 'T', 'n']
units_list = ['ppmv', '$\mathregular{cm^{-3}}$', 'K']

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/zonal_means/f.e20.FXSD.f19_f19.001.cam.h0.2000-01.nc', 'r', format='NETCDF4')
lats = fname_uni.variables['lat'][:]
lons = fname_uni.variables['lon'][:]
fname_uni.close()

def calc_z3_zon_t_av(symbol, levs):  
    if levs == 88:
        if symbol == 'CO':
            fname = netCDF4.Dataset('/nfs/a249/earfw/NCAS/NOHO_mee_2013_2015_sofie/NOHO_mee_2013_2015_sofie.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
        if symbol == 'CO2':
            fname = netCDF4.Dataset('/nfs/ncas/earfw/NOHO_PAPER/CO2/NOHO_mee_2013_2015_CO2_4DAN/NOHO_mee_2013_2015_CO2_4DAN.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    z3 = np.zeros([1,levs,96,144])
    z3[:] = fname.variables['Z3'][:]*(1.e-3)
    z3_zon_av = np.mean(z3[:], axis=3)
    z3_zon_t_av = np.mean(z3_zon_av[:], axis=0) 
    fname.close()
    return z3_zon_t_av

def calc_species_zon_av(symbol, levs):  
    if levs == 88:
        if symbol == 'CO':
            fname = netCDF4.Dataset('/nfs/a249/earfw/NCAS/NOHO_mee_2013_2015_sofie/NOHO_mee_2013_2015_sofie.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
        if symbol == 'CO2':
            fname = netCDF4.Dataset('/nfs/ncas/earfw/NOHO_PAPER/CO2/NOHO_mee_2013_2015_CO2_4DAN/NOHO_mee_2013_2015_CO2_4DAN.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4') 
    species_dat = np.zeros([1,levs,96,144])
    species_dat = fname.variables[symbol][:]*(1.e6) 
    species_tm = np.mean(species_dat[:], axis=0)
    species_zon_av = np.mean(species_tm[:], axis=2)
    fname.close()
    return species_zon_av

def calc_ratio(param1, param2, levs):
    ratio = np.zeros([levs,96])
    for i in range(0,levs):
        for j in range(0,96): 
            ratio[i,j] = param1[i,j] / param2[i,j]
    return ratio

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
    
def calc_z3_zon_t_av_weighted(levs, lowlat, highlat):
    z3_zon_t_av = calc_z3_zon_t_av(symbol, levs)
    z3_zon_t_av_weighted = calc_cos_factor(z3_zon_t_av, levs, lowlat, highlat)
    return z3_zon_t_av_weighted

def calc_profiles(param, levs, lowlat, highlat):
    param_weighted = calc_cos_factor(param, levs, lowlat, highlat)
    return param_weighted

def plot_1d_ratio(name, config, units, z3, species, lowlat, highlat, color, plot_no):
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = species[::-1]
    y = z3[::-1]
    plt.plot(x, y, color=color, label=config)
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
        plt.xlabel('%s' %name, fontsize=12)
        plt.ylabel('Altitude [km]', fontsize=12)
    if plot_no == 4:
        plt.xlabel('%s' %name, fontsize=12)
        plt.tick_params(labelleft='off')
    if plot_no == 5:
        plt.xlabel('%s' %name, fontsize=12)
        plt.tick_params(labelleft='off')
    plt.ylim(60,160)
    plt.xlim(0,10)
    if config == 'waccm-x' and plot_no == 2:
        plt.legend(loc=1)
    return

year = 2013
month = 7
name = species_list[4]
symbol = symbol_list[4]
units = units_list[0]
chemistry = True
global_only = False
save = True
# For ratio:
name2 = species_list[3]
symbol2 = symbol_list[3]

if units == 'ppmv':
    units_print = 'ppmv'
elif units == '$\mathregular{cm^{-3}}$':
    units_print = 'cm-3'
elif units == 'K':
    units_print = 'K'

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

# 1D Ratio Plot Code:
waccm_species = calc_species_zon_av(symbol, 88)
waccm_species2 = calc_species_zon_av(symbol2, 88)
waccm_ratio = calc_ratio(waccm_species, waccm_species2, 88)

waccmx_species = calc_species_zon_av(symbol, 145)
waccmx_species2 = calc_species_zon_av(symbol2, 145)
waccmx_ratio = calc_ratio(waccmx_species, waccmx_species2, 145)

for i in range(a,b):    
    lowlat = i * step
    highlat = (i * step) + step
    lowlat_no = int((lowlat * 1.875) - 90)
    highlat_no = int((highlat * 1.875) - 90)
    waccm_z3_weighted = calc_z3_zon_t_av_weighted(88, lowlat, highlat)
    waccmx_z3_weighted = calc_z3_zon_t_av_weighted(145, lowlat, highlat)
    waccm_species_profile = calc_profiles(waccm_ratio, 88, lowlat, highlat) 
    waccmx_species_profile = calc_profiles(waccmx_ratio, 145, lowlat, highlat)
    plot_1d_ratio('%s / %s ratio' %(symbol, symbol2), 'waccm', units, waccm_z3_weighted, waccm_species_profile, lowlat, highlat, 'k', i)      
    plot_1d_ratio('%s / %s ratio' %(symbol, symbol2), 'waccm-x', units, waccmx_z3_weighted, waccmx_species_profile, lowlat, highlat, 'b', i)
if save == True:
    plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/CO_CO2_ratio/%s/%s_%s_ratio_month%s_profile_lat_bands.jpg' %(year, name, name2, month), bbox_inches='tight', dpi=300)

plt.show()