import matplotlib
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from datetime import date
import calendar

deg = unichr(176)
delta = unichr(916)
k_B = 1.38064852e-23

lats = np.zeros([96])
levs = np.zeros([145])
p_int = np.zeros([96])
T_int = np.zeros([96])
n = np.zeros([96])
p_time_series = np.zeros([177,96])
T_time_series = np.zeros([177,96])
n_time_series = np.zeros([177,96])

p_seasonal_av = np.zeros([12,96])
T_seasonal_av = np.zeros([12,96])
n_seasonal_av = np.zeros([12,96])
p_seasonal_av_diff = np.zeros([14,12])
T_seasonal_av_diff = np.zeros([14,12])
n_seasonal_av_diff = np.zeros([14,12])

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/zonal_means/f.e20.FXSD.f19_f19.001.cam.h0.2000-01.nc', 'r', format='NETCDF4')
levs = fname_uni.variables['lev'][:]
lats = fname_uni.variables['lat'][:]
fname_uni.close()

def add_months(sourcedate, months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month // 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year, month)[1])
    return date(year, month, day)

def add_years(sourcedate, years):
    year = sourcedate.year - 1 + years
    return date(year)

def interp_p(altitude, gpheight):
    for i in range(0,96):
        p_int[i] = np.interp(altitude, gpheight[:,i][::-1], levs[:][::-1])     
    return p_int

def interp_T(pressure_int, temp):
    for i in range(0,96):
        T_int[i] = np.interp(pressure_int[i], levs[:], temp[:,i])
    return T_int

def calc_n(pressure_int, temp_int):
    for i in range(0,96):
        n[i] = (pressure_int[i] * 100) / (temp_int[i] * k_B)
    return n

def calc_mthly_means(param, month):
    start = month
    end = (13 * 12) + month
    step = np.arange(start, end, 12)
    param_mthly = np.zeros([len(step),96])
    for i in range(0,len(step)):
        year = step[i]
        param_mthly[i,:] = param[year,:]
    param_mthly_mean = np.mean(param_mthly, axis=0)
    return param_mthly_mean

def make_time_series_arrays(param, altitude):
    param_time_series = np.zeros([177,96])
    z3_dat = np.zeros([1,145,96])
    z3 = np.zeros([145,96])
    T = np.zeros([96])
    start = date(2000,1,1)
    date_list = []
    fdates = []
    ndates = np.zeros(177)
    for i in range(0,177):
        date_i = add_months(start, i)
        date_list.append(date_i)
        ndate = matplotlib.dates.date2num(date_i)
        ndates[i] = ndate
        fdates.append(date_list[i].strftime('%Y-%m'))
    for i in range(0,177):
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/zonal_means/f.e20.FXSD.f19_f19.001.cam.h0.%s.nc' %fdates[i], 'r', format='NETCDF4')
        z3_dat = fname.variables['Z3'][:]
        z3 = np.mean(z3_dat, axis=0)
        T_dat = fname.variables['T'][:]
        T = np.mean(T_dat, axis=0)
        fname.close()
        p_time_series[i,:] = interp_p(altitude, z3)
        T_time_series[i,:] = interp_T(p_time_series[i,:], T)
        n_time_series[i,:] = calc_n(p_time_series[i,:], T_time_series[i,:])
    if param == 'p':          
        param_time_series[:,:] = p_time_series[:,:]
    if param == 'T':
        param_time_series[:,:] = T_time_series[:,:]
    if param == 'n':
        param_time_series[:,:] = n_time_series[:,:]
    return param_time_series
    
def make_seasonal_av_arrays(param):
    param_seasonal_av = np.zeros([12,96])
    for i in range(0,12):
        param_seasonal_av[i,:] = calc_mthly_means(param, i)
    return param_seasonal_av

def make_seasonal_av_diff_arrays(param0, param1):
    baseline = np.zeros([12])
    param_time_series_set_year = np.zeros([14,12,96])
    param_time_series_set_year_lat_av = np.zeros([14,12])
    param_seasonal_av_diff = np.zeros([14,12]) 
    baseline = np.mean(param1, axis=1)
    step = np.arange(0, 180, 12)
    for i in range(0,14): 
        param_time_series_set_year[i,:,:] = param0[step[i]:step[i+1],:]
    param_time_series_set_year_lat_av = np.mean(param_time_series_set_year, axis=2)
    for i in range(0,14):
        param_seasonal_av_diff[i,:] = param_time_series_set_year_lat_av[i,:] - baseline
    return param_seasonal_av_diff

def plot_1d(param0, param1, label, altitude):
    plt.figure(figsize=(8,4))
    #baseline = np.zeros([12])
    #baseline = np.mean(param0, axis=1)
    months = np.arange(0,12,1)
    x = months
    years = ['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013']
    colors_light = ['#f44336', '#e91e63', '#9c27b0', '#3f51b5', '#2196f3', '#00bcd4', '#009688', '#4caf50', '#cddc39', '#ffeb3b', '#ff9800', '#795548', '#9e9e9e', '#607d8b']
    #plt.plot(x, baseline)
    for i in range(0,14):
        y = param1[i,:]
        plt.plot(x, y, color=colors_light[i], label=years[i])
    labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    plt.xticks(months, labels, rotation='45')
    plt.yticks()
    plt.xlabel('Month')
    plt.ylabel('%s' %label)
    plt.legend(bbox_to_anchor=(1.05, 1.01))
    plt.title('Difference from 2000-14 seasonal mean (%s m)' %altitude)
    #plt.savefig('/nfs/a328/eecwk/waccm-x/figures/diff_from_seasonal_mean', bbox_inches='tight', dpi=300)
    plt.show()
    return

def plot_2d(param, label):
    months = np.arange(0,12,1)
    x, y = np.meshgrid(months, lats)
    plt.figure(figsize=(6,4))
    #diffs = np.arange(700,1100,10)
    labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    plt.xticks(months, labels, rotation='45')
    plt.yticks(np.arange(-90,120,30))
    z = param[:,:]
    zT = z.transpose()
    #plt.contourf(x[:,:], y[:,:], zT[:,:], diffs, cmap=plt.get_cmap('seismic'))
    plt.contourf(x[:,:], y[:,:], zT[:,:], cmap=plt.get_cmap('seismic'))
    cbar = plt.colorbar()
    cbar.set_label('%s' %label)
    plt.xlabel('Month')
    plt.ylabel('Latitude [%s]' %deg)
    plt.title('Seasonal mean 2000-14 (%s m)' %altitude)
    #plt.savefig('/nfs/a328/eecwk/waccm-x/figures/%s_200km_%s' %label %tscale, dpi=300)
    plt.show()
    return

altitude = 200000

p_time_series = make_time_series_arrays('p', altitude)
T_time_series = make_time_series_arrays('T', altitude)
n_time_series = make_time_series_arrays('n', altitude)

p_seasonal_av[:,:] = make_seasonal_av_arrays(p_time_series)
T_seasonal_av[:,:] = make_seasonal_av_arrays(T_time_series)
n_seasonal_av[:,:] = make_seasonal_av_arrays(n_time_series)

p_seasonal_av_diff[:,:] = make_seasonal_av_diff_arrays(p_time_series, p_seasonal_av)
T_seasonal_av_diff[:,:] = make_seasonal_av_diff_arrays(T_time_series, T_seasonal_av)
n_seasonal_av_diff[:,:] = make_seasonal_av_diff_arrays(n_time_series, n_seasonal_av)

plot_1d(p_seasonal_av, p_seasonal_av_diff, 'p [HPa]', altitude)
plot_1d(T_seasonal_av, T_seasonal_av_diff, '%sT [K]' %delta, altitude)
plot_1d(n_seasonal_av, n_seasonal_av_diff, 'n [$\mathregular{m^{-3}}$]', altitude)

#plot_2d(p_seasonal_av, 'p', altitude)
#plot_2d(T_seasonal_av, 'T', altitude)
#plot_2d(n_seasonal_av, 'n', altitude)