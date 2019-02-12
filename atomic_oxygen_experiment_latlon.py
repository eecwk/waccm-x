import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec

deg = unichr(176)
p_int_w = np.zeros([96,144])
p_int_wx = np.zeros([96,144])

def get_fixed_variables(config, year, month, symbol):
    if config == 'waccmx':
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-%s.nc' %(year, month), 'r', format='NETCDF4')
    if config == 'waccm':
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-%s.nc' %(year, month), 'r', format='NETCDF4')  
    var = fname.variables[symbol][:]
    fname.close()
    return var

def get_variables(config, levs, year, month, symbol):
    tracer_dat = np.zeros([1,levs,96,144])
    tracer = np.zeros([levs,96,144])
    if config == 'waccmx':
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s-%s.nc' %(year, month), 'r', format='NETCDF4')
    if config == 'waccm':
        fname = netCDF4.Dataset('/nfs/a265/earfw/SD_WACCM4/john_ca_paper_JDmif_nad4cad7.cam2.h0.%s-%s.nc' %(year, month), 'r', format='NETCDF4')
    tracer_dat = fname.variables[symbol][:]
    tracer = np.mean(tracer_dat, axis=0)
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

def plot_2d_latlon(tracer, diffs, name, units):
    x, y = np.meshgrid(lons, lats)
    plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    plt.xlabel('Longitude [%s]' %deg, fontsize=12)
    plt.ylabel('Latitude [%s]' %deg, fontsize=12)
    altitude_km = altitude/1000
    plt.title('%s-%s (%s km)' %(year, month, altitude_km), fontsize=14)
    cbar = plt.colorbar(orientation='vertical', cmap=plt.get_cmap('viridis'))
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.show()

def plot_2d_latlon_sub(config, tracer, diffs, name, units, plot_no):
    plt.subplot(gs1[plot_no])
    x, y = np.meshgrid(lons, lats)
    ax = plt.contourf(x[:,:], y[:,:], tracer[:,:], diffs)
    #ax = plt.contourf(x[:,:], y[:,:], tracer[:,:])
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
#altitude = 75000
#altitude = 85000
altitude = 90000
#altitude = 95000
#altitude = 120000

levs_w = get_fixed_variables('waccm', year, month, 'lev')
levs_wx = get_fixed_variables('waccmx', year, month, 'lev')
lats = get_fixed_variables('waccmx', year, month, 'lat')
lons = get_fixed_variables('waccmx', year, month, 'lon')

z3_w = get_variables('waccm', 88, year, month, 'Z3')
p_int_w[:,:] = interp_p(altitude, z3_w, levs_w)
z3_wx = get_variables('waccmx', 145, year, month, 'Z3')
p_int_wx[:,:] = interp_p(altitude, z3_wx, levs_wx)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9,3))
gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.1, hspace=0.1)
'''
#diffs = np.arange(10000,90000,2500)
diffs = np.arange(160000,360000,10000)
o_w = get_variables('waccm', 88, year, month, 'O')
o_w_int = interp_tracer(o_w*1.e+6, levs_w, p_int_w)
plot_2d_latlon_sub('waccm', o_w_int, diffs, 'atomic_oxygen', 'ppmv', 0)
o_wx = get_variables('waccmx', 145, year, month, 'O')
o_wx_int = interp_tracer(o_wx*1.e+6, levs_wx, p_int_wx)
plot_2d_latlon_sub('waccmx', o_wx_int, diffs, 'atomic_oxygen', 'ppmv', 1)
plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s/atomic_oxygen_month%s_120km_latlon_ppmv.jpg' %(year, month), bbox_inches='tight', dpi=300)
'''

#diffs = np.arange(0,3.2,0.1)
#diffs = np.arange(0,1.2,0.05)
diffs = np.arange(0,3.6,0.1)
o3_w = get_variables('waccm', 88, year, month, 'O3')
o3_w_int = interp_tracer(o3_w*1.e+6, levs_w, p_int_w)
plot_2d_latlon_sub('waccm', o3_w_int, diffs, 'ozone', 'ppmv', 0)
o3_wx = get_variables('waccmx', 145, year, month, 'O3')
o3_wx_int = interp_tracer(o3_wx*1.e+6, levs_wx, p_int_wx)
plot_2d_latlon_sub('waccmx', o3_wx_int, diffs, 'ozone', 'ppmv', 1)
#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/atomic_oxygen_experiment/john_ca_paper_JDmif_nad4cad7/%s/ozone_month%s_90km_latlon_ppmv.jpg' %(year, month), bbox_inches='tight', dpi=300)

plt.show()