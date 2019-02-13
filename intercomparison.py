import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import math

N_A = 6.02214086e+23
R = 8.3144598

deg = unichr(176)
set_year = 2014
set_month = 1
set_day = 6
waccmx_start_day = 4

# NRLMSISE
fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/Meteostuff-master/NRLMSIS_F90/output/20140106-00000_1/nrlmsise_atomicoxygen_20140106-00000.nc', 'r', format='NETCDF4')
msis_dat = np.zeros([1,145,96,144])
msis_dat = fname.variables['O'][:]
msis_o = msis_dat[0,:,:,72]
fname.close()

# SABER
fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/SABER/NIGHTLY/atox_athy_night_YY%s_V5.3_fixedfnight_SV2.nc' %set_year, 'r', format='NETCDF4')
o = fname.variables['qatox'][:]
T = fname.variables['ktemp'][:]
p = fname.variables['pressure'][:]
lats_saber = fname.variables['lat'][:]
lons_saber = fname.variables['lon'][:]
alt = fname.variables['alt'][:]
year = fname.variables['year'][:]
day = fname.variables['day'][:]
fname.close
lats_bin_mid_saber = np.arange(-87.5, 92.5, 5)

def make_saber_array(tracer, set_year, set_day):
    tracer_bin = np.zeros([16,36])
    for j in range(0,16):
        for k in range(0,36):            
            k_min = (k*5) - 90
            k_max = (k*5) - 85      
            #time_lat_bin = np.append(np.where((lats_saber > k_min) & (lats_saber < k_max) & (year == set_year) & (day == set_day)), True)
            time_lat_bin = np.append(np.where((lons_saber > 160) & (lons_saber < 200) & (lats_saber > k_min) & (lats_saber < k_max) & (year == set_year) & (day == set_day)), True) 
            tracer_scans = tracer[time_lat_bin, j]
            tracer_scans_good = np.array([])
            for i in range(0, len(tracer_scans)):
                if tracer_scans[i] > 0:
                    tracer_scans_good = np.append(tracer_scans_good, tracer_scans[i])
            if len(tracer_scans_good) < 2:
                tracer_scans_good = 0
            tracer_scans_mean = np.mean(tracer_scans_good)
            tracer_bin[j,k] = tracer_scans_mean                               
    return tracer_bin      

def calc_saber_cos_factor(tracer_bin, lowlat, highlat):
    tracer_weighted = np.zeros(16)    
    for j in range(0,16):    
        sig_cos_x = 0
        sig_cos = 0
        for k in range(lowlat, highlat):
            sig_cos_x = sig_cos_x + (math.cos(math.radians(lats_bin_mid_saber[k])) * tracer_bin[j][k])
            sig_cos = sig_cos + math.cos(math.radians(lats_bin_mid_saber[k]))         
            if  k == (highlat - 1):
                tracer_weighted[j] = sig_cos_x / sig_cos
    return tracer_weighted

def calc_saber_conc_profile(tracer, lowlat, highlat):
    T_saber = make_saber_array(T, set_year, set_day)
    T_saber_weighted = calc_saber_cos_factor(T_saber, lowlat, highlat)
    #p_saber = make_saber_array(p, set_year, set_day)
    #p_saber_weighted = calc_saber_cos_factor(p_saber, lowlat, highlat)
    tracer_weighted = calc_saber_cos_factor(tracer, lowlat, highlat)
    tracer_weighted_conc = np.zeros(16)  
    for i in range(0,16):
        if T_saber_weighted[i] > 0:
            tracer_weighted_conc[i] = (tracer_weighted[i] * 1.e-6 * N_A * 100 * p[i]) / (R * T_saber_weighted[i])
        else:
            tracer_weighted_conc[i] = np.nan
    return tracer_weighted_conc

# WACCM-X
fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.2014-01-04-00000.nc', 'r', format='NETCDF4')
levs_waccmx = fname.variables['lev'][:]
lats_waccmx = fname.variables['lat'][:]
fname.close()

def make_waccmx_array(symbol, factor):  
    fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-0%s-00000.nc' %(set_year, set_month, waccmx_start_day), 'r', format='NETCDF4') 
    tracer_dat = np.zeros([7,145,96,144])
    tracer_dat = fname.variables[symbol][:]*factor
    tracer = tracer_dat[2,:,:,72]
    fname.close()
    return tracer

def calc_waccmx_cos_factor(tracer, lowlat, highlat):
    tracer_weighted = np.zeros(145)    
    for j in range(0,145):    
        sig_cos_x = 0
        sig_cos = 0
        for k in range(lowlat, highlat):
            sig_cos_x = sig_cos_x + (math.cos(math.radians(lats_waccmx[k])) * tracer[j][k])
            sig_cos = sig_cos + math.cos(math.radians(lats_waccmx[k]))         
            if  k == (highlat - 1):
                tracer_weighted[j] = sig_cos_x / sig_cos
    return tracer_weighted

def calc_waccmx_conc_profile(tracer, lowlat, highlat):
    T_waccmx = np.zeros([1,145,96,144])
    T_waccmx = make_waccmx_array('T', 1)
    T_waccmx_weighted = calc_waccmx_cos_factor(T_waccmx, lowlat, highlat)    
    tracer_weighted = calc_waccmx_cos_factor(tracer, lowlat, highlat)
    tracer_weighted_conc = np.zeros(145)  
    for i in range(0,145):
        tracer_weighted_conc[i] = (tracer_weighted[i] * 1.e-6 * N_A * 100 * levs_waccmx[i]) / (R * T_waccmx_weighted[i])
    return tracer_weighted_conc

# Plot
def plot_1d(tracer_weighted, alt_weighted, factor, name, lowlat, highlat, config, units, color, plot_no):
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    lowlat_no = int((lowlat * factor) - 90)
    highlat_no = int((highlat * factor) - 90) 
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = tracer_weighted[::-1]
    y = alt_weighted[::-1]
    plt.plot(x, y, color=color, label=config)
    plt.xlim(0,8.e+11)
    plt.ylim(50,200)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
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
    if config == 'msis' and plot_no == 2:
        plt.legend(loc=1)
    return

def setup_plot_1d(tracer, alt, step, factor, name, config, units, color):
    for i in range(0,6):    
        lowlat = i * step
        highlat = (i * step) + step
        if config == 'saber':
            #saber_tracer_weighted = calc_saber_cos_factor(tracer, lowlat, highlat)
            saber_tracer_weighted_conc = calc_saber_conc_profile(tracer, lowlat, highlat)
            saber_alt_weighted = calc_saber_cos_factor(alt, lowlat, highlat)
            #plot_1d(saber_tracer_weighted, saber_alt_weighted, factor, name, lowlat, highlat, config, units, color, i)
            plot_1d(saber_tracer_weighted_conc, saber_alt_weighted, factor, name, lowlat, highlat, config, units, color, i)
        if config == 'waccm-x':
            #waccmx_tracer_weighted = calc_waccmx_cos_factor(tracer, lowlat, highlat)
            waccmx_tracer_weighted_conc = calc_waccmx_conc_profile(tracer, lowlat, highlat)
            waccmx_alt_weighted = calc_waccmx_cos_factor(alt, lowlat, highlat)       
            #plot_1d(waccmx_tracer_weighted, waccmx_alt_weighted, factor, name, lowlat, highlat, config, units, color, i)
            plot_1d(waccmx_tracer_weighted_conc, waccmx_alt_weighted, factor, name, lowlat, highlat, config, units, color, i)
        if config == 'msis':
            msis_tracer_weighted = calc_waccmx_cos_factor(tracer, lowlat, highlat)
            waccmx_alt_weighted = calc_waccmx_cos_factor(alt, lowlat, highlat)
            plot_1d(msis_tracer_weighted, waccmx_alt_weighted, factor, name, lowlat, highlat, config, units, color, i)  
    return

saber_o = make_saber_array(o, set_year, set_day)
saber_alt = make_saber_array(alt, set_year, set_day)
waccmx_o =  make_waccmx_array('O',1)
waccmx_alt =  make_waccmx_array('Z3',1.e-3)

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.1, hspace=0.1)
plt.suptitle('%s, DOY=%s, night, 180%sE' %(set_year, set_day, deg), fontsize=16)
setup_plot_1d(saber_o, saber_alt, 6, 5, 'atomic_oxygen', 'saber', 'ppmv', 'm')
setup_plot_1d(waccmx_o, waccmx_alt, 16, 1.875, 'atomic_oxygen', 'waccm-x', 'ppmv', 'k')
setup_plot_1d(msis_o, waccmx_alt, 16, 1.875, 'atomic_oxygen', 'msis', 'cm-3', 'b')

#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/waccmx_v_saber_comparisons/atomic_oxygen_%s-%s-%s_night.jpg' %(set_year, set_month, set_day), bbox_inches='tight', dpi=300)     
plt.show()