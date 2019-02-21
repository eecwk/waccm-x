import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import math

N_A = 6.02214086e+23
R = 8.3144598
deg = unichr(176)

years = [2009, 2009, 2009, 2009, 2014, 2014, 2014, 2014]
months = [3, 6, 9, 12, 3, 6, 9, 12]
days = [20, 21, 22, 21, 20, 21, 23, 21]
waccmx_start_days = [14, 20, 19, 19, 15, 21, 20, 20]
doy = [79, 172, 265, 355, 79, 172, 266, 355]

event = 0
set_year = years[event]
set_month = months[event]
set_day = days[event]
waccmx_start_day = waccmx_start_days[event]
saber_doy = doy[event]
waccmx_select_day = set_day - waccmx_start_day

# SABER
fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/SABER/NIGHTLY/atox_athy_night_YY%s_V5.3_fixedfnight_SV2.nc' %set_year, 'r', format='NETCDF4')
o = fname.variables['qatox'][:]
h = fname.variables['qathy'][:]
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
            time_lat_bin = np.append(np.where((lons_saber > 160) & (lons_saber < 200) & (lats_saber > k_min) & (lats_saber < k_max) & (year == set_year) & (day == saber_doy)), True) 
            tracer_scans = tracer[time_lat_bin, j]
            tracer_scans_good = np.array([])
            for i in range(0, len(tracer_scans)):
                if tracer_scans[i] > 0:
                    tracer_scans_good = np.append(tracer_scans_good, tracer_scans[i])
            if len(tracer_scans_good) < 2:
                tracer_scans_good = np.nan
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

def calc_saber_lat_means(tracer_bin, lowlat, highlat):
    tracer_select = np.zeros(16)
    for j in range(0,16):
        tracer_select[j] = np.mean(tracer_bin[j,lowlat:highlat])
    return tracer_select

def calc_saber_conc_profile(tracer, lowlat, highlat):
    saber_T_lat_band = calc_saber_lat_means(saber_T, lowlat, highlat)
    tracer_weighted = calc_saber_cos_factor(tracer, lowlat, highlat)
    tracer_weighted_conc = np.zeros(16)  
    for i in range(0,16):
        if saber_T_lat_band[i] > 0:
            tracer_weighted_conc[i] = (tracer_weighted[i] * 1.e-6 * N_A * 100 * p[i]) / (R * saber_T_lat_band[i])
        else:
            tracer_weighted_conc[i] = np.nan
    return tracer_weighted_conc

# WACCM-X
fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.2014-01-04-00000.nc', 'r', format='NETCDF4')
levs_waccmx = fname.variables['lev'][:]
lats_waccmx = fname.variables['lat'][:]
fname.close()

def make_waccmx_array(symbol, factor):  
    if set_month < 10:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-%s-00000.nc' %(set_year, set_month, waccmx_start_day), 'r', format='NETCDF4') 
    else:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-%s-%s-00000.nc' %(set_year, set_month, waccmx_start_day), 'r', format='NETCDF4') 
    tracer_dat = np.zeros([7,145,96,144])
    tracer_dat = fname.variables[symbol][:]*factor
    tracer = np.mean(tracer_dat[waccmx_select_day,:,:,64:81], axis=2)
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

def calc_waccmx_lat_means(tracer_bin, lowlat, highlat):
    tracer_select = np.zeros(145)
    for j in range(0,145):
        tracer_select[j] = np.mean(tracer_bin[j,lowlat:highlat])
    return tracer_select

def calc_waccmx_conc_profile(tracer, lowlat, highlat):
    waccmx_T_lat_band = calc_waccmx_lat_means(waccmx_T, lowlat, highlat)
    tracer_weighted = calc_waccmx_cos_factor(tracer, lowlat, highlat)
    tracer_weighted_conc = np.zeros(145)
    for i in range(0,145):
        tracer_weighted_conc[i] = (tracer_weighted[i] * 1.e-6 * N_A * 100 * levs_waccmx[i]) / (R * waccmx_T_lat_band[i])
    return tracer_weighted_conc

# NRLMSISE
fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/output_data/nrlmsise_ghp7_20090320-00000.nc', 'r', format='NETCDF4')
lats_msis = fname.variables['latitude'][:]
fname.close()

def make_msis_array(symbol):  
    if set_month < 10:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/output_data/nrlmsise_ghp7_%s0%s%s-00000.nc' %(set_year, set_month, set_day), 'r', format='NETCDF4') 
    else:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/output_data/nrlmsise_ghp7_%s%s%s-00000.nc' %(set_year, set_month, set_day), 'r', format='NETCDF4') 
    tracer_dat = np.zeros([1,145,96,144])
    tracer_dat = fname.variables[symbol][:]
    tracer = np.mean(tracer_dat[0,:,:,64:81], axis=2)
    fname.close()
    return tracer

def calc_msis_cos_factor(tracer, lowlat, highlat):
    tracer_weighted = np.zeros(145)    
    for j in range(0,145):    
        sig_cos_x = 0
        sig_cos = 0
        for k in range(lowlat, highlat):
            sig_cos_x = sig_cos_x + (math.cos(math.radians(lats_msis[k])) * tracer[j][k])
            sig_cos = sig_cos + math.cos(math.radians(lats_msis[k]))         
            if  k == (highlat - 1):
                tracer_weighted[j] = sig_cos_x / sig_cos
    return tracer_weighted

def calc_msis_lat_means(tracer_bin, lowlat, highlat):
    tracer_select = np.zeros(145)
    for j in range(0,145):
        tracer_select[j] = np.mean(tracer_bin[j,lowlat:highlat])
    return tracer_select

# Plot
def plot_1d(tracer_weighted, alt_weighted, factor, xlim, name, lowlat, highlat, config, units, color, plot_no):
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    lowlat_no = abs(int((lowlat * factor) - 90))
    highlat_no = abs(int((highlat * factor) - 90)) 
    if plot_no < 3:
        plt.title('%s%sS to %s%sS' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    if plot_no > 2:
        plt.title('%s%sN to %s%sN' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = tracer_weighted[::-1]
    y = alt_weighted[::-1]
    plt.plot(x, y, color=color, label=config)
    #
    plt.xlim(0,xlim)
    plt.ylim(50,200)
    #
    #plt.xlim(1.e+6,1.e+10)
    #plt.ylim(150,500)
    #plt.xscale('log')
    #
    #plt.xlim(1.e+4,1.e+7)
    #plt.ylim(150,500)
    #plt.xscale('log')
    #
    #plt.xlim(500,1200)
    #plt.ylim(150,500)
    #
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
    if config == 'saber' and plot_no == 2:
        plt.legend(loc=1)
    return

def setup_plot_1d_chem(tracer, alt, step, factor, xlim, name, config, units, color):
    for i in range(0,6):    
        lowlat = i * step
        highlat = (i * step) + step
        if config == 'waccm-x':
            waccmx_tracer_weighted_conc = calc_waccmx_conc_profile(tracer, lowlat, highlat)
            waccmx_alt_lat_band = calc_waccmx_lat_means(alt, lowlat, highlat)
            plot_1d(waccmx_tracer_weighted_conc, waccmx_alt_lat_band, factor, xlim, name, lowlat, highlat, config, units, color, i)
        if config == 'msis':
            msis_tracer_weighted = calc_msis_cos_factor(tracer, lowlat, highlat)
            msis_alt_lat_band = calc_msis_lat_means(alt, lowlat, highlat)
            plot_1d(msis_tracer_weighted, msis_alt_lat_band, factor, xlim, name, lowlat, highlat, config, units, color, i)
        if config == 'saber':
            saber_tracer_weighted_conc = calc_saber_conc_profile(tracer, lowlat, highlat)
            saber_alt_lat_band = calc_saber_lat_means(alt, lowlat, highlat)
            plot_1d(saber_tracer_weighted_conc, saber_alt_lat_band, factor, xlim, name, lowlat, highlat, config, units, color, i)    
            if set_month < 10:
                plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/%s_%s-0%s-%s.jpg' %(name, set_year, set_month, set_day), bbox_inches='tight', dpi=300)
            else:
                plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/%s_%s-%s-%s.jpg' %(name, set_year, set_month, set_day), bbox_inches='tight', dpi=300)                
    return

def setup_plot_1d_phys(tracer, alt, step, factor, xlim, name, config, units, color):
    for i in range(0,6):    
        lowlat = i * step
        highlat = (i * step) + step
        if config == 'waccm-x':
            waccmx_tracer_lat_band = calc_waccmx_lat_means(tracer, lowlat, highlat)
            waccmx_alt_lat_band = calc_waccmx_lat_means(alt, lowlat, highlat)
            plot_1d(waccmx_tracer_lat_band, waccmx_alt_lat_band, factor, xlim, name, lowlat, highlat, config, units, color, i)
        if config == 'msis':
            msis_tracer_lat_band = calc_msis_lat_means(tracer, lowlat, highlat)
            msis_alt_lat_band = calc_msis_lat_means(alt, lowlat, highlat)
            plot_1d(msis_tracer_lat_band, msis_alt_lat_band, factor, xlim, name, lowlat, highlat, config, units, color, i)
        if config == 'saber':
            saber_tracer_lat_band = calc_saber_lat_means(tracer, lowlat, highlat)
            saber_alt_lat_band = calc_saber_lat_means(alt, lowlat, highlat)
            plot_1d(saber_tracer_lat_band, saber_alt_lat_band, factor, xlim, name, lowlat, highlat, config, units, color, i)            
            if set_month < 10:
                plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/%s_%s-0%s-%s.jpg' %(name, set_year, set_month, set_day), bbox_inches='tight', dpi=300)
            else:
                plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/%s_%s-%s-%s.jpg' %(name, set_year, set_month, set_day), bbox_inches='tight', dpi=300)
    return

saber_alt = make_saber_array(alt, set_year, set_day)
saber_o = make_saber_array(o, set_year, set_day)
saber_h = make_saber_array(h, set_year, set_day)
saber_T = make_saber_array(T, set_year, set_day)

waccmx_alt =  make_waccmx_array('Z3',1.e-3)
waccmx_o =  make_waccmx_array('O',1)
waccmx_h = make_waccmx_array('H',1)
waccmx_T = make_waccmx_array('T', 1)

msis_alt = make_msis_array('Z3')
msis_o = make_msis_array('O')
msis_h = make_msis_array('H')
msis_T = make_msis_array('T')

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.1, hspace=0.1)
if set_month < 10:
    plt.suptitle('%s/0%s/%s, night, 20%sW to 20%sE' %(set_year, set_month, set_day, deg, deg), fontsize=16)
else:
    plt.suptitle('%s/%s/%s, night, 20%sW to 20%sE' %(set_year, set_month, set_day, deg, deg), fontsize=16)

setup_plot_1d_chem(waccmx_o, waccmx_alt, 16, 1.875, 8.e+11, 'atomic_oxygen', 'waccm-x', 'cm-3', 'k')
setup_plot_1d_chem(msis_o, msis_alt, 16, 1.875, 8.e+11, 'atomic_oxygen', 'msis', 'cm-3', 'b')
setup_plot_1d_chem(saber_o, saber_alt, 6, 5, 8.e+11, 'atomic_oxygen', 'saber', 'cm-3', 'm')

#setup_plot_1d_chem(waccmx_h, waccmx_alt, 16, 1.875, 4.e+8, 'atomic_hydrogen', 'waccm-x', 'cm-3', 'k')
#setup_plot_1d_chem(msis_h, msis_alt, 16, 1.875, 4.e+8, 'atomic_hydrogen', 'msis', 'cm-3', 'b')
#setup_plot_1d_chem(saber_h, saber_alt, 6, 5, 4.e+8, 'atomic_hydrogen', 'saber', 'cm-3', 'm')

#setup_plot_1d_phys(waccmx_T, waccmx_alt, 16, 1.875, 1000, 'temperature', 'waccm-x', 'K', 'k')
#setup_plot_1d_phys(msis_T, msis_alt, 16, 1.875, 1000, 'temperature', 'msis', 'K', 'b')
#setup_plot_1d_phys(saber_T, saber_alt, 6, 5, 1000, 'temperature', 'saber', 'K', 'm')

plt.show()