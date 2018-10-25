import matplotlib
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib import colors, dates
from datetime import datetime, timedelta, date
import calendar

deg = unichr(176)

def add_months(sourcedate,months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month // 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return date(year,month,day)

def waccm_time_series(species, symbol, plot_no):
    
    name = species
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
    plt.figure(figsize=(9,4))
    y, x = np.meshgrid(lats_waccm, ndates)
  
    if name == 'atomic_oxygen':
        diffs = np.arange(80000,365000,5000)
        plt.contourf(x[:,:], y[:,:], species_series[:,:], diffs)
        cbar = plt.colorbar(ticks=[np.arange(5.e+4,4.e+5,5.e+4)], orientation='vertical')

    elif name == 'ozone':
        diffs = np.arange(0,0.00162,0.00002)
        plt.contourf(x[:,:], y[:,:], species_series[:,:], diffs)
        cbar = plt.colorbar(orientation='vertical')

    cbar.set_label('%s vmr [ppmv]' %species)
    n = np.arange(0,180,12)
    markers = ndates[n]
    labels = ['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
    plt.xticks(markers, labels, rotation='45')
    plt.xlabel('Year')
    plt.ylabel('Latitude [%s]' %deg)
    plt.title('120 km')
    plt.show()

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

#waccm_time_series('atomic_oxygen', 'O', 0)
waccm_time_series('ozone', 'O3', 1)