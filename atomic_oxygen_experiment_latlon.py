import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

N_A = 6.02214086e+23
R = 8.3144598
deg = unichr(176)

years = [2009, 2009, 2009, 2009, 2014, 2014, 2014, 2014]
months = [3, 6, 9, 12, 3, 6, 9, 12]
days = [20, 21, 22, 21, 20, 21, 23, 21]
waccmx_start_days = [14, 20, 19, 19, 15, 21, 20, 20]
doy = [79, 172, 265, 355, 79, 172, 266, 355]

altitude = 85
event = 6
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

def make_saber_array(tracer):
    tracer_bin = np.zeros([16,36,36])
    for j in range(0,16):
        for k in range(0,36):
            k_min = (k*5) - 90
            k_max = (k*5) - 85
            for l in range(0,36):            
                l_min = (l*5)
                l_max = (l*10) + 10
                time_lat_lon_bin = np.append(np.where((lats_saber > k_min) & (lats_saber < k_max) & (lons_saber > l_min) & (lons_saber < l_max) & (year == set_year) & (day == saber_doy)), True) 
                tracer_scans = tracer[time_lat_lon_bin, j]                
                tracer_scans_good = np.array([])
                for i in range(0, len(tracer_scans)):
                    if tracer_scans[i] > 0:
                        tracer_scans_good = np.append(tracer_scans_good, tracer_scans[i])
                if len(tracer_scans_good) < 2:
                    tracer_scans_good = np.nan
                tracer_scans_mean = np.mean(tracer_scans_good)
                tracer_bin[j,k,l] = tracer_scans_mean                               
    return tracer_bin

def interp_saber_p(altitude, gpheight, pressure):        
    p_int = np.zeros([36,36])
    for i in range(0,36):
        for j in range(0,36):
            p_int[i,j] = np.interp(altitude, gpheight[:,i,j][::-1], pressure[:][::-1])     
    return p_int

def interp_saber_tracer(tracer, pressure, p_int):
    tracer_int = np.zeros([36,36])
    for i in range(0,36):
        for j in range(0,36):    
            tracer_int[i,j] = np.interp(p_int[i,j], pressure[:], tracer[:,i,j])
    return tracer_int

def calc_saber_conc_profile(tracer, p_int, T_int):
    tracer_conc = np.zeros([36,36])  
    for i in range(0,36):
        for j in range(0,36):
            if T_int[i,j] > 0:
                tracer_conc[i,j] = (tracer[i,j] * 1.e-6 * N_A * 100 * p_int[i,j]) / (R * T_int[i,j])
            else:
                tracer_conc[i,j] = np.nan
    return tracer_conc

# WACCM-X
fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.2014-01-04-00000.nc', 'r', format='NETCDF4') 
levs_waccmx = fname.variables['lev'][:]
lats_waccmx = fname.variables['lat'][:]
lons_waccmx = fname.variables['lon'][:]
fname.close()

def make_waccmx_array(symbol, factor):
    if set_month < 10:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-0%s-%s-00000.nc' %(set_year, set_month, waccmx_start_day), 'r', format='NETCDF4') 
    else:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-%s-%s-00000.nc' %(set_year, set_month, waccmx_start_day), 'r', format='NETCDF4')    
    tracer_dat = np.zeros([7,145,96,144])
    tracer_dat = fname.variables[symbol][:]*factor    
    tracer = tracer_dat[waccmx_select_day,:,:,:]    
    fname.close()
    return tracer

def interp_waccmx_p(altitude, gpheight, levs):        
    p_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):
            p_int[i,j] = np.interp(altitude, gpheight[:,i,j][::-1], levs[:][::-1])     
    return p_int

def interp_waccmx_tracer(tracer, pressure, p_int):
    tracer_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):    
            tracer_int[i,j] = np.interp(p_int[i,j], pressure[:], tracer[:,i,j])
    return tracer_int

def calc_waccmx_conc_profile(tracer, p_int, T_int):
    tracer_conc = np.zeros([96,144])  
    for i in range(0,96):
        for j in range(0,144):
            tracer_conc[i,j] = (tracer[i,j] * 1.e-6 * N_A * 100 * p_int[i,j]) / (R * T_int[i,j])
    return tracer_conc

# NRLMSISE
fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/output_data/nrlmsise_ghp7_20090320-00000.nc', 'r', format='NETCDF4')
levs_msis = fname.variables['level'][:]
lats_msis = fname.variables['latitude'][:]
lons_msis = fname.variables['longitude'][:]
fname.close()

def make_msis_array(symbol):  
    if set_month < 10:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/output_data/nrlmsise_ghp7_%s0%s%s-00000.nc' %(set_year, set_month, set_day), 'r', format='NETCDF4') 
    else:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/nrlmsise/output_data/nrlmsise_ghp7_%s%s%s-00000.nc' %(set_year, set_month, set_day), 'r', format='NETCDF4') 
    tracer_dat = np.zeros([1,145,96,144])
    tracer_dat = fname.variables[symbol][:]
    tracer = tracer_dat[0,:,:,:]
    fname.close()
    return tracer

def interp_msis_p(altitude, gpheight, levs):        
    p_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):
            p_int[i,j] = np.interp(altitude, gpheight[:,i,j][::-1], levs[:][::-1])     
    return p_int

def interp_msis_tracer(tracer, pressure, p_int):
    tracer_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):    
            tracer_int[i,j] = np.interp(p_int[i,j], pressure[:], tracer[:,i,j])
    return tracer_int

def calc_msis_conc_profile(tracer, p_int, T_int):
    tracer_conc = np.zeros([96,144])  
    for i in range(0,96):
        for j in range(0,144):
            tracer_conc[i,j] = (tracer[i,j] * 1.e-6 * N_A * 100 * p_int[i,j]) / (R * T_int[i,j])
    return tracer_conc

# Plot
def plot_2d_latlon_sub(tracer, lats, lons, diffs, config, name, units, plot_no):
    plt.subplot(gs1[plot_no])
    x, y = np.meshgrid(lons, lats)
    ax = plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    #ax = plt.contourf(x[:,:], y[:,:], tracer[:,:])
    plt.xlim(0,360)
    plt.ylim(-90,90)
    plt.xlabel('Longitude [%s]' %deg, fontsize=12)
    plt.title('%s' %config, fontsize=14)
    if plot_no == 0:
        plt.ylabel('Latitude [%s]' %deg, fontsize=12)
    if plot_no == 1:
        plt.tick_params(labelleft='off')
    if plot_no == 2:
        plt.tick_params(labelleft='off')
        cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(ax, cax=cbar_ax, orientation='vertical')
        cbar.set_label('%s [%s]' %(name, units), fontsize=12)
        cbar.ax.tick_params(labelsize=12)
    return

saber_alt = make_saber_array(alt)
saber_o = make_saber_array(o)
saber_h = make_saber_array(h)
saber_T = make_saber_array(T)
saber_p_int = interp_saber_p(altitude, saber_alt, p)
saber_o_int = interp_saber_tracer(saber_o, p, saber_p_int)
saber_h_int = interp_saber_tracer(saber_h, p, saber_p_int)
saber_T_int = interp_saber_tracer(saber_T, p, saber_p_int)
saber_o_int_conc = calc_saber_conc_profile(saber_o_int, saber_p_int, saber_T_int)
saber_h_int_conc = calc_saber_conc_profile(saber_h_int, saber_p_int, saber_T_int)

waccmx_alt = make_waccmx_array('Z3', 1.e-3)
waccmx_o = make_waccmx_array('O', 1)
waccmx_h = make_waccmx_array('H', 1)
waccmx_T = make_waccmx_array('T', 1)
waccmx_p_int = interp_waccmx_p(altitude, waccmx_alt, levs_waccmx)
waccmx_o_int = interp_waccmx_tracer(waccmx_o, levs_waccmx, waccmx_p_int)
waccmx_h_int = interp_waccmx_tracer(waccmx_h, levs_waccmx, waccmx_p_int)
waccmx_T_int = interp_waccmx_tracer(waccmx_T, levs_waccmx, waccmx_p_int)
waccmx_o_int_conc = calc_waccmx_conc_profile(waccmx_o_int, waccmx_p_int, waccmx_T_int)
waccmx_h_int_conc = calc_waccmx_conc_profile(waccmx_h_int, waccmx_p_int, waccmx_T_int)

msis_alt = make_msis_array('Z3')
msis_o = make_msis_array('O')
msis_h = make_msis_array('H')
msis_T = make_msis_array('T')
msis_p_int = interp_msis_p(altitude, msis_alt, levs_msis)
msis_o_int = interp_msis_tracer(msis_o, levs_msis, msis_p_int)
msis_h_int = interp_msis_tracer(msis_h, levs_msis, msis_p_int)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(wspace=0.1, hspace=0.5)
if set_month < 10:
    plt.suptitle('%s/0%s/%s, night, (%s km)' %(set_year, set_month, set_day, altitude), fontsize=16)
else:
    plt.suptitle('%s/%s/%s, night, (%s km)' %(set_year, set_month, set_day, altitude), fontsize=16)

'''
diffs = np.arange(1.e+11,7.e+11,2.e+10)
plot_2d_latlon_sub(waccmx_o_int_conc, lats_waccmx, lons_waccmx, diffs, 'WACCM-X', 'atomic_oxygen', 'cm-3', 0)
plot_2d_latlon_sub(msis_o_int, lats_msis, lons_msis, diffs, 'MSIS', 'atomic_oxygen', 'cm-3', 1)
plot_2d_latlon_sub(saber_o_int_conc, np.arange(-87.5,92.5,5), np.arange(5,365,10), diffs, 'SABER', 'atomic_oxygen', 'cm-3', 2)
if set_month < 10:
    plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/atomic_oxygen_%s-0%s-%s.jpg' %(set_year, set_month, set_day), bbox_inches='tight', dpi=300)
else:
    plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/atomic_oxygen_%s-%s-%s.jpg' %(set_year, set_month, set_day), bbox_inches='tight', dpi=300)
'''

diffs = np.arange(5.e+7,40.e+7,1.e+7)
plot_2d_latlon_sub(waccmx_h_int_conc, lats_waccmx, lons_waccmx, diffs, 'WACCM-X', 'atomic_hydrogen', 'cm-3', 0)
plot_2d_latlon_sub(msis_h_int, lats_msis, lons_msis, diffs, 'MSIS', 'atomic_hydrogen', 'cm-3', 1)
plot_2d_latlon_sub(saber_h_int_conc, np.arange(-87.5,92.5,5), np.arange(5,365,10), diffs, 'SABER', 'atomic_hydrogen', 'cm-3', 2)
if set_month < 10:
    plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/atomic_hydrogen_%s-0%s-%s_%skm.jpg' %(set_year, set_month, set_day, altitude), bbox_inches='tight', dpi=300)
else:
    plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/intercomparison/atomic_hydrogen_%s-%s-%s_%skm.jpg' %(set_year, set_month, set_day, altitude), bbox_inches='tight', dpi=300)

plt.show()