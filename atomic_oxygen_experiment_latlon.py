import matplotlib.pyplot as plt
import netCDF4
import numpy as np

deg = unichr(176)

levs = np.zeros([145])
lats = np.zeros([96])
lons = np.zeros([144])
p_int = np.zeros([96,144])

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2000-01.nc', 'r', format='NETCDF4')
levs = fname_uni.variables['lev'][:]
lats = fname_uni.variables['lat'][:]
lons = fname_uni.variables['lon'][:]
fname_uni.close()

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
    
def calc_tracer_tmean(symbol):
    tracer_dat = np.zeros([1,145,96,144])
    tracer = np.zeros([145,96,144])
    fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2014-01.nc', 'r', format='NETCDF4')
    tracer_dat = fname.variables[symbol][:]*(1.e+6)
    tracer = np.mean(tracer_dat, axis=0)
    return tracer
    
def make_arrays(fixed_altitude):
    z3_dat = np.zeros([1,145,96,144])
    z3 = np.zeros([145,96,144])
    fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2014-01.nc', 'r', format='NETCDF4')
    z3_dat = fname.variables['Z3'][:]
    z3 = np.mean(z3_dat, axis=0)
    fname.close()            
    p_int[:,:] = interp_p(fixed_altitude, z3)
    return p_int

def plot_2d_latlon(tracer, name, units):
    x, y = np.meshgrid(lons, lats)
    diffs = np.arange(0,3.0,0.1)
    plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    plt.xlabel('Longitude [%s]' %deg, fontsize=12)
    plt.ylabel('Latitude [%s]' %deg, fontsize=12)
    altitude_km = altitude/1000
    plt.title('%s km' %altitude_km, fontsize=14)
    cbar = plt.colorbar(orientation='vertical', cmap=plt.get_cmap('viridis'))
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.show()

altitude = 95000
p_int = make_arrays(altitude)
ozone = calc_tracer_tmean('O3')
ozone_int = interp_tracer(ozone)
plot_2d_latlon(ozone_int, 'ozone', 'ppmv')