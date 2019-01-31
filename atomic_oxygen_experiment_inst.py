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

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.2014-01-04-00000.nc', 'r', format='NETCDF4')
lats = fname_uni.variables['lat'][:]
lons = fname_uni.variables['lon'][:]
fname_uni.close()

def calc_z3_zon_av(levs):  
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-0%s-00000.nc' %(year, month, start_day), 'r', format='NETCDF4')
    z3 = np.zeros([7,levs,96,144])
    z3[:] = fname.variables['Z3'][:]*(1.e-3)
    z3_t = z3[2,:,:,:]
    z3_zon_av = np.mean(z3_t, axis=2)
    fname.close()    
    return z3_zon_av

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
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-0%s-00000.nc' %(year, month, start_day), 'r', format='NETCDF4') 
    species_dat = np.zeros([7,levs,96,144])
    if symbol == 'T':
        species_dat = fname.variables[symbol][:]
    elif symbol == 'n':
        species_dat = fname.variables['T'][:]
    else:
       species_dat = fname.variables[symbol][:]*(1.e6) 
    species_t = species_dat[2,:,:,:]
    species_zon_av = np.mean(species_t[:], axis=2)
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
    
def calc_z3_zon_t_av_weighted(levs, lowlat, highlat):
    z3_zon_t_av = calc_z3_zon_av(levs)
    z3_zon_t_av_weighted = calc_cos_factor(z3_zon_t_av, levs, lowlat, highlat)
    return z3_zon_t_av_weighted

def get_lev(levs):
    if levs == 88:
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-0%s.nc' %(year, month), 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-0%s-00000.nc' %(year, month, start_day), 'r', format='NETCDF4') 
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

def plot_1d_multi(name, config, units, z3, species, lowlat, highlat, color, plot_no):
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = species[::-1]
    y = z3[::-1]
    plt.plot(x, y, color=color, label=config)
    plt.ylim(77,100)
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
            plt.xlim(0,50000)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,8.e+11)
    if name == 'ozone':
        if units == 'ppmv':
            plt.xlim(1.e-8,1.e+1)
        if units == '$\mathregular{cm^{-3}}$':
            plt.xlim(0,7.e+8)
    if name == 'atomic_hydrogen':
        if units == 'ppmv':
            plt.xlim(0,20)
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,5.e+8)
    if name == 'carbon_dioxide':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        if units == 'ppmv':
            plt.xlim(0,500)
            1==1
        if units == '$\mathregular{cm^{-3}}$':    
            plt.xlim(0,1.e+13)    
    if name == 'temperature':
        plt.xlim(0,800)
    if name == 'density':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlim(0,3.e+12)
    if config == 'waccm-x' and plot_no == 2:
        plt.legend(loc=1)
    return

year = 2014
month = 1
start_day = 4
day = 6
name = species_list[0]
symbol = symbol_list[0]
units = units_list[0]
chemistry = True
save = False

if units == 'ppmv':
    units_print = 'ppmv'
elif units == '$\mathregular{cm^{-3}}$':
    units_print = 'cm-3'
elif units == 'K':
    units_print = 'K'

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.1, hspace=0.1)
plt.suptitle('%s, DOY=%s, night' %(year, day), fontsize=16)
step = 16
a = 0
b = 6

#waccm_species = calc_species_zon_av(symbol, 88)
waccmx_species = calc_species_zon_av(symbol, 145)

# 1D Plot Code
for i in range(a,b):    
    lowlat = i * step
    highlat = (i * step) + step
    lowlat_no = int((lowlat * 1.875) - 90)
    highlat_no = int((highlat * 1.875) - 90)
    #waccm_z3_weighted = calc_z3_zon_t_av_weighted(88, lowlat, highlat)
    waccmx_z3_weighted = calc_z3_zon_t_av_weighted(145, lowlat, highlat)    
    if chemistry == True:
        if units == 'ppmv':
            #waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat)
        elif units == '$\mathregular{cm^{-3}}$':
            #waccm_species_profile = calc_conc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_conc_profiles(waccmx_species, 145, lowlat, highlat)
    else:
        if symbol == 'T':
            #waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat) 
        elif symbol == 'n':
            #waccm_species_profile = calc_profiles(waccm_species, 88, lowlat, highlat)
            waccmx_species_profile = calc_profiles(waccmx_species, 145, lowlat, highlat)        
    #plot_1d_multi(name, 'waccm', units, waccm_z3_weighted, waccm_species_profile, lowlat, highlat, 'k', i)
    plot_1d_multi(name, 'waccm-x', units, waccmx_z3_weighted, waccmx_species_profile, lowlat, highlat, 'b', i)

if save == True:
    plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s/%s_month%s_profile_lat_bands_%s.jpg' %(year, name, month, units_print), bbox_inches='tight', dpi=300)
        
plt.show()