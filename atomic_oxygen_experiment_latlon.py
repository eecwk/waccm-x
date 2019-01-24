import matplotlib.pyplot as plt
import netCDF4
import numpy as np

deg = unichr(176)

levs = np.zeros([145])
lats = np.zeros([96])
lons = np.zeros([144])
p_int = np.zeros([96,144])

def get_fixed_variables(year, month, symbol):
    fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-%s.nc' %(year, month), 'r', format='NETCDF4')
    var = fname.variables[symbol][:]
    fname.close()
    return var

def get_variables(config, year, month, symbol):
    if config == 'waccmx':
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-%s.nc' %(year, month), 'r', format='NETCDF4')
        tracer_dat = np.zeros([1,145,96,144])
        tracer = np.zeros([145,96,144])
    if config == 'waccm':
        fname == 1
        tracer_dat = np.zeros([1,88,96,144])
        tracer = np.zeros([88,96,144])
    tracer_dat = fname.variables[symbol][:]
    tracer = np.mean(tracer_dat, axis=0)
    fname.close()
    return tracer

def interp_p(altitude, gpheight):        
    for i in range(0,96):
        for j in range(0,144):
            p_int[i,j] = np.interp(altitude, gpheight[:,i,j][::-1], levs[:][::-1])     
    return p_int

def interp_tracer(tracer):
    tracer_int = np.zeros([96,144])
    for i in range(0,96):
        for j in range(0,144):    
            tracer_int[i,j] = np.interp(p_int[i,j], levs[:], tracer[:,i,j])
    return tracer_int

def plot_2d_latlon(tracer, diffs, name, units):
    x, y = np.meshgrid(lons, lats)
    plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    plt.xlabel('Longitude [%s]' %deg, fontsize=12)
    plt.ylabel('Latitude [%s]' %deg, fontsize=12)
    altitude_km = altitude/1000
    plt.title('%s km' %altitude_km, fontsize=14)
    cbar = plt.colorbar(orientation='vertical', cmap=plt.get_cmap('viridis'))
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.show()

year = '2014'
month = '01'
altitude = 120000

levs = get_fixed_variables(year, month, 'lev')
lats = get_fixed_variables(year, month, 'lat')
lons = get_fixed_variables(year, month, 'lon')

z3 = get_variables('waccmx', year, month, 'Z3')
p_int[:,:] = interp_p(altitude, z3)

o = get_variables('waccmx', year, month, 'O')
o_int = interp_tracer(o*1.e+6)
diffs = np.arange(180000,340000,5000)
plot_2d_latlon(o_int, diffs, 'atomic_oxygen', 'ppmv')

#o3 = get_variables('waccmx', year, month, 'O')
#o3_int = interp_tracer(o3*1.e+6)
#diffs = np.arange(0,3.0,0.1)
#plot_2d_latlon(o3_int, diffs, 'ozone', 'ppmv')