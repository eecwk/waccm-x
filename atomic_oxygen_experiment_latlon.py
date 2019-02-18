import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec

N_A = 6.02214086e+23
R = 8.3144598
deg = unichr(176)
p_int_w = np.zeros([96,144])
p_int_wx = np.zeros([96,144])

def get_fixed_variables(config, year, month, symbol):
    if config == 'waccmx':
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-%s-0%s-00000.nc' %(year, month, start_day), 'r', format='NETCDF4') 
    var = fname.variables[symbol][:]
    fname.close()
    return var

def get_variables(config, levs, year, month, symbol):
    tracer_dat = np.zeros([1,levs,96,144])
    tracer = np.zeros([levs,96,144])
    if config == 'waccmx':
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_daily_inst/f.e20.FXSD.f19_f19.001.cam.h2.%s-%s-0%s-00000.nc' %(year, month, start_day), 'r', format='NETCDF4')
    tracer_dat = fname.variables[symbol][:]
    tracer = tracer_dat[2,:,:,:]
    fname.close()
    return tracer

def interp_p(altitude, gpheight, levs):        
    p_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):
            p_int[i,j] = np.interp(altitude, gpheight[:,i,j][::-1], levs[:][::-1])     
    return p_int

def interp_tracer(tracer, levs, p_int):
    tracer_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):    
            tracer_int[i,j] = np.interp(p_int[i,j], levs[:], tracer[:,i,j])
    return tracer_int

def calc_conc_profile(tracer):
    T_wx = np.zeros([1,145,96,144])
    T_wx = get_variables('waccmx', 145, year, month, 'T')
    T_wx = interp_tracer(T_wx, levs_wx, p_int_wx)
    tracer_conc = np.zeros([96,144])  
    for i in range(0,96):
        for j in range(0,144):
            tracer_conc[i,j] = (tracer[i,j] * 1.e-6 * N_A * 100 * p_int_wx[i,j]) / (R * T_wx[i,j])
    return tracer_conc

def plot_2d_latlon(tracer, diffs, name, units):
    x, y = np.meshgrid(lons, lats)
    #plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    plt.contourf(x[:,:], y[:,:], tracer[:,:])
    plt.xlabel('Longitude [%s]' %deg, fontsize=12)
    plt.ylabel('Latitude [%s]' %deg, fontsize=12)
    altitude_km = altitude/1000
    plt.title('%s/%s/0%s, night, (%s km)' %(year, month, set_day, altitude_km), fontsize=14)
    cbar = plt.colorbar(orientation='vertical', cmap=plt.get_cmap('viridis'))
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.show()

def plot_2d_latlon_sub(config, tracer, diffs, name, units, plot_no):
    x, y = np.meshgrid(lons, lats)
    #ax = plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    ax = plt.contourf(x[:,:], y[:,:], tracer[:,:])
    plt.xlabel('Longitude [%s]' %deg, fontsize=12)
    altitude_km = altitude/1000
    plt.title('%s %s-%s (%s km)' %(config, year, month, altitude_km), fontsize=14)
    if plot_no == 0:
        plt.ylabel('Latitude [%s]' %deg, fontsize=12)
    if plot_no == 1:
        plt.tick_params(labelleft='off')
        cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(ax, cax=cbar_ax, orientation='vertical')
        cbar.set_label('%s [%s]' %(name, units), fontsize=12)
        cbar.ax.tick_params(labelsize=12)

year = '2014'
month = '01'
start_day = 4
set_day = 6
altitude = 100000

levs_wx = get_fixed_variables('waccmx', year, month, 'lev')
lats = get_fixed_variables('waccmx', year, month, 'lat')
lons = get_fixed_variables('waccmx', year, month, 'lon')

z3_wx = get_variables('waccmx', 145, year, month, 'Z3')
p_int_wx[:,:] = interp_p(altitude, z3_wx, levs_wx)
o_wx = get_variables('waccmx', 145, year, month, 'O')
#o_wx_int = interp_tracer(o_wx*1.e+6, levs_wx, p_int_wx)
o_wx_int = interp_tracer(o_wx, levs_wx, p_int_wx)
o_wx_int_conc = calc_conc_profile(o_wx_int)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9,6))
gs1 = gridspec.GridSpec(1, 1)
gs1.update(wspace=0.1, hspace=0.1)
diffs = np.arange(160000,360000,10000)
#plot_2d_latlon(o_wx_int, diffs, 'atomic_oxygen', 'ppmv')
plot_2d_latlon(o_wx_int_conc, diffs, 'atomic_oxygen', 'cm-3')
plt.show()