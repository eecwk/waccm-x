import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import colors

deg = unichr(176)
delta = unichr(916)
k_B = 1.38064852e-23

def waccm_meridional_slice(month, species, symbol, levs, plot_no):
    
    if levs == 70:
        fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/Tao_Na_Paper_FW5_JCMIFd5.cam.h0.0001-0%s.nc' %month, 'r', format='NETCDF4')
    if levs == 145:
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2001-0%s.nc' %month, 'r', format='NETCDF4')
   
    name = species
    z3 = np.zeros([1,levs ,96,144])  
    species_dat = np.zeros([1,levs,96,144]) 
    lats = fname.variables['lat'][:]
    z3[:] = fname.variables['Z3'][:]*(1.e-3)
    z3_zon_av = np.mean(z3[:], axis=3)
    z3_zon_mer_av = np.mean(z3_zon_av[:], axis=2)
    z3_zon_mer_t_av = np.mean(z3_zon_mer_av[:], axis=0) 
    species_dat = fname.variables[symbol][:]*(1.e6)
    species_tm = np.mean(species_dat[:], axis=0)
    species_zon_av = np.mean(species_tm[:], axis=2)
    fname.close()
  
    plt.subplot(gs1[plot_no])
    x, y = np.meshgrid(lats, z3_zon_mer_t_av)
    plt.xlabel('Latitude [%s]' %deg, fontsize=12)
    plt.xticks(np.arange(-90,120,30), fontsize=12) 
    plt.yticks(np.arange(0,220,20), fontsize=12)   
    plt.ylim(80,200)
    plt.axhline(y=146, color='w', linewidth=1, linestyle=':')
    
    if name == 'atomic_oxygen':
        diffs = [1.e+2, 175.e+0, 25.e+1, 325.e+0, 4.e+2, 475.e+0, 55.e+1, 625.e+0, 7.e+2, 775.e+0, 85.e+1, 925.e+0, 1.e+3, 175.e+1, 25.e+2, 325.e+1, 4.e+3, 475.e+1, 55.e+2, 625.e+1, 7.e+3, 775.e+1, 85.e+2, 925.e+1, 1.e+4, 175.e+2, 25.e+3, 325.e+2, 4.e+4, 475.e+2, 55.e+3, 625.e+2, 7.e+4, 775.e+2, 85.e+3, 925.e+2, 1.e+5, 175.e+3, 25.e+4, 325.e+3, 4.e+5, 475.e+3, 55.e+4, 625.e+3, 7.e+5, 775.e+3, 85.e+4, 925.e+3, 1.e+6]
        ax = plt.contourf(x[:,:], y[:,:], species_zon_av[:,:], diffs, norm=colors.LogNorm())
    elif name == 'ozone':
        diffs = np.arange(0,10.6,0.2)        
        ax = plt.contourf(x[:,:], y[:,:], species_zon_av[:,:], diffs)

    if plot_no == 0:
        plt.title('WACCM')
        plt.ylabel('Altitude [km]', fontsize=12)
    elif plot_no == 1:
        plt.title('WACCM-X')
        plt.tick_params(labelleft='off')
        cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(ax, cax=cbar_ax, ticks=[1.e+0, 1.e+1, 1.e+2, 1.e+3, 1.e+4, 1.e+5, 1.e+6], orientation='vertical')
        cbar.set_label('%s vmr [ppmv]' %species, fontsize=12)
        cbar.ax.tick_params(labelsize=12)
    return

month = 1

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(11,5))
gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.1, hspace=0.1)

if month == 1:
    fig.suptitle('January', fontsize=16)
elif month == 7:
    fig.suptitle('July', fontsize=16)

waccm_meridional_slice(month, 'atomic_oxygen', 'O', 70, 0)
waccm_meridional_slice(month, 'atomic_oxygen', 'O', 145, 1)

#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/.jpg' %name, dpi=300)
plt.show()