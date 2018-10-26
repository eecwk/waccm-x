import matplotlib
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib import colors, dates
from datetime import datetime, timedelta, date
import calendar
import matplotlib.gridspec as gridspec

deg = unichr(176)

def add_months(sourcedate,months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month // 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return date(year,month,day)

def waccm_time_series(species, symbol, plot_no):
    
    name = species
    solar_dat = np.zeros([177])
    lats_waccm = np.zeros([96])
    z3 = np.zeros([1,145,96,144])
    z3_zon_av = np.zeros([1,145,96])
    z3_zon_mer_av = np.zeros([1,145])
    z3_zon_mer_t_av = np.zeros([145])
    species_series = np.zeros([177,96])
    species_dat = np.zeros([1,145,96,144])
    species_tm = np.zeros([145,96,144])

    for i in range(0,177):    
        fname_waccm = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s.nc' %fdates[i], 'r', format='NETCDF4')
        solar_dat[i] = fname_waccm.variables['f107'][:]
        lats_waccm[:] = fname_waccm.variables['lat'][:]
        z3[:] = fname_waccm.variables['Z3'][:]
        z3_zon_av[:] = np.mean(z3[:], axis=3)
        z3_zon_mer_av[:] = np.mean(z3_zon_av[:], axis=2)
        z3_zon_mer_t_av[:] = np.mean(z3_zon_mer_av[:], axis=0)
        species_dat = fname_waccm.variables[symbol][:]*(1.e6)
        species_tm = np.mean(species_dat[:], axis=0)
        species_zon_av = np.mean(species_tm[:], axis=2)
        species_zon_av_120 = species_zon_av[44]
        species_series[i] = species_zon_av_120

    fname_waccm.close()
    n = np.arange(0,180,12)
    markers = ndates[n]
    labels = ['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
  
    y, x = np.meshgrid(lats_waccm, ndates)

    if name == 'f107':
        plt.subplot(gs1[plot_no])
        plt.plot(ndates[:], solar_dat[:], color='#4B0082')
        plt.ylabel('f107 [$\mathregular{10^{-22}}$ $\mathregular{W}$ $\mathregular{m^{-2}}$ $\mathregular{Hz^{-1}}$]')
        plt.xticks(markers, labels, rotation='45')
        plt.margins(x=0)
        plt.tick_params(labelbottom='off')
    
    else:    
        plt.subplot(gs1[plot_no])
        
        if name == 'atomic_oxygen':
            diffs = np.arange(80000,365000,5000)
            ax = plt.contourf(x[:,:], y[:,:], species_series[:,:], diffs, cmap=plt.get_cmap('inferno'))
            plt.xticks(markers, labels, rotation='45')
            plt.tick_params(labelbottom='off')
            plt.ylabel('Latitude [%s]' %deg)
            cbaxes = fig.add_axes([0.83, 0.552, 0.015, 0.16]) 
            cbar = plt.colorbar(ax, cax=cbaxes, ticks=[np.arange(5.e+4,4.e+5,5.e+4)], orientation='vertical')

        elif name == 'ozone':
            diffs = np.arange(0,0.00162,0.00002)
            ax = plt.contourf(x[:,:], y[:,:], species_series[:,:], diffs, cmap=plt.get_cmap('inferno'))
            plt.xticks(markers, labels, rotation='45')
            plt.tick_params(labelbottom='off')     
            plt.ylabel('Latitude [%s]' %deg)
            cbaxes = fig.add_axes([0.83, 0.334, 0.015, 0.16]) 
            cbar = plt.colorbar(ax, cax=cbaxes, ticks=[np.arange(0,0.0018,0.0003)], orientation='vertical')
     
        elif name == 'atomic_hydrogen':
            diffs = np.arange(6,16.6,0.1)
            ax = plt.contourf(x[:,:], y[:,:], species_series[:,:], diffs, cmap=plt.get_cmap('inferno'))
            plt.xticks(markers, labels, rotation='45')
            plt.xlabel('Year')
            plt.ylabel('Latitude [%s]' %deg)
            cbaxes = fig.add_axes([0.83, 0.122, 0.015, 0.16]) 
            cbar = plt.colorbar(ax, cax=cbaxes, ticks=[np.arange(6,18,2)], orientation='vertical')
        
        cbar.set_label('%s vmr [ppmv]' %species)

start = date(2000,1,1)
date_list = []
fdates = []
ndates = np.zeros(177)

for i in range(0,177):
    date_i = add_months(start,i)
    date_list.append(date_i)
    ndate = matplotlib.dates.date2num(date_i)
    ndates[i] = ndate
    fdates.append(date_list[i].strftime('%Y-%m'))

fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(9,12))
fig.subplots_adjust(right=0.8)
fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.1)
gs1 = gridspec.GridSpec(4, 1)
gs1.update(wspace=0.1, hspace=0.1)

waccm_time_series('f107', 'O', 0)
waccm_time_series('atomic_oxygen', 'O', 1)
waccm_time_series('ozone', 'O3', 2)
waccm_time_series('atomic_hydrogen', 'H', 3)

plt.savefig('/nfs/a328/eecwk/waccm-x/figures/time_series.jpg', dpi=300)
plt.show()