import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import math

deg = unichr(176)
set_year = 2014
set_month = 1
set_day = 6
waccmx_start_day = 4

# SABER
fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/SABER/NIGHTLY/atox_athy_night_YY%s_V5.3_fixedfnight_SV2.nc' %set_year, 'r', format='NETCDF4')
o = fname.variables['qatox'][:]*(1.e+6)
lats_saber = fname.variables['lat'][:]
alt = fname.variables['alt'][:]
year = fname.variables['year'][:]
day = fname.variables['day'][:]
fname.close
lats_bin_mid_saber = np.arange(-87.5, 92.5, 5)

def make_retrieval_arrays(tracer, set_year, set_day):
    tracer_bin = np.zeros([16,36])
    for j in range(0,16):
        for k in range(0,36):            
            k_min = (k*5) - 90
            k_max = (k*5) - 85      
            time_lat_bin = np.append(np.where((lats_saber > k_min) & (lats_saber < k_max) & (year == set_year) & (day == set_day)), True)        
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

def calc_cos_factor_saber(tracer_bin, lowlat, highlat):
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

# WACCM-X
fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.2014-01-04-00000.nc', 'r', format='NETCDF4')
lats_waccmx = fname.variables['lat'][:]
fname.close()

def make_waccmx_array(symbol, factor):  
    fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-0%s-00000.nc' %(set_year, set_month, waccmx_start_day), 'r', format='NETCDF4') 
    tracer_dat = np.zeros([7,145,96,144])
    tracer_dat = fname.variables[symbol][:]*factor
    tracer_t = tracer_dat[2,:,:,:]
    tracer_zon_av = np.mean(tracer_t[:], axis=2)
    fname.close()
    return tracer_zon_av

def calc_cos_factor_waccmx(tracer, lowlat, highlat):
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

# Plot
def plot_1d(tracer_weighted, alt_weighted, factor, name, lowlat, highlat, config, units, color, plot_no):
    lowlat_no = int((lowlat * factor) - 90)
    highlat_no = int((highlat * factor) - 90) 
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = tracer_weighted[::-1]
    y = alt_weighted[::-1]
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
    if config == 'waccm-x' and plot_no == 2:
        plt.legend(loc=1)
    return

def setup_plot_1d(tracer_bin, alt_bin, step, factor, name, config, units, color):
    for i in range(0,6):    
        lowlat = i * step
        highlat = (i * step) + step
        if config == 'saber':
            tracer_weighted_saber = calc_cos_factor_saber(tracer_bin, lowlat, highlat)
            alt_weighted_saber = calc_cos_factor_saber(alt_bin, lowlat, highlat)
            plot_1d(tracer_weighted_saber, alt_weighted_saber, factor, name, lowlat, highlat, config, units, color, i)
        if config == 'waccm-x':
            tracer_weighted_waccmx = calc_cos_factor_waccmx(tracer_bin, lowlat, highlat)
            alt_weighted_waccmx = calc_cos_factor_waccmx(alt_bin, lowlat, highlat)       
            plot_1d(tracer_weighted_waccmx, alt_weighted_waccmx, factor, name, lowlat, highlat, config, units, color, i)
    return

saber_o_bin = make_retrieval_arrays(o, set_year, set_day)
saber_alt_bin = make_retrieval_arrays(alt, set_year, set_day)
waccmx_o =  make_waccmx_array('O',1.e+6)
waccmx_alt =  make_waccmx_array('Z3',1.e-3)

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.1, hspace=0.1)
plt.suptitle('%s, DOY=%s, night' %(set_year, set_day), fontsize=16)
setup_plot_1d(saber_o_bin, saber_alt_bin, 6, 5, 'atomic_oxygen', 'saber', 'ppmv', 'k')
setup_plot_1d(waccmx_o, waccmx_alt, 16, 1.875, 'atomic_oxygen', 'waccm-x', 'ppmv', 'b')

#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/.jpg', bbox_inches='tight', dpi=300)     
plt.show()