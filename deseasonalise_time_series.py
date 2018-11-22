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
p_yrly = np.zeros([177,96])
T_yrly = np.zeros([177,96])
n_yrly = np.zeros([177,96])
p_all_mthly = np.zeros([12,96])
T_all_mthly = np.zeros([12,96])
n_all_mthly = np.zeros([12,96])

p_indv_yrs = np.zeros([14,12,96])
T_indv_yrs = np.zeros([14,12,96])
n_indv_yrs = np.zeros([14,12,96])

p_indv_yrs_diff = np.zeros([14,12])
T_indv_yrs_diff = np.zeros([14,12])
n_indv_yrs_diff = np.zeros([14,12])

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
        T_int[i] = np.interp(pressure_int[i], levs[:][::-1], temp[:,i][::-1])
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

def slice_yrly_array(param):
    yrly_slice = np.zeros([14,12,96])
    step = np.arange(0, 180, 12)
    for i in range(0,14):
        yrly_slice[i,:,:] = param[step[i]:step[i+1],:]
    return yrly_slice

def make_yrly_diff_arrays(param0, param1):
    baseline = np.zeros([12])
    param_indv_yrs_global = np.zeros([14,12])
    param_indv_yrs_diff = np.zeros([14,12])
    baseline = np.mean(param0, axis=1)
    param_indv_yrs_global = np.mean(param1, axis=2)
    for i in range(0,14):
        param_indv_yrs_diff[i,:] = param_indv_yrs_global[i,:] - baseline
    return param_indv_yrs_diff

def make_desolar_arrays(output):
    param_all_mthly = np.zeros([12,96])
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
        p_yrly[i,:] = interp_p(200000, z3)
        T_yrly[i,:] = interp_T(p_int, T)
        n_yrly[i,:] = calc_n(p_int, T_int)      
    for i in range(0,12):
        if output == 'p':
            param_all_mthly[i,:] = calc_mthly_means(p_yrly, i)
        if output == 'T':
            param_all_mthly[i,:] = calc_mthly_means(T_yrly, i)
        if output == 'n':
            param_all_mthly[i,:] = calc_mthly_means(n_yrly, i)
    return param_all_mthly

def plot_1d(param0, param1, label):
    plt.figure(figsize=(6,4))
    months = np.arange(0,12,1)
    x = months
    years = ['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013']
    #colors_dark = ['#b71c1c', '#880e4f', '#4a148c', '#1a237e', '#0d47a1', '#006064', '#004d40', '#1b5e20', '#827717', '#f57f17', '#e65100', '#3e2723', '#212121', '#263238']
    colors_light = ['#f44336', '#e91e63', '#9c27b0', '#3f51b5', '#2196f3', '#00bcd4', '#009688', '#4caf50', '#cddc39', '#ffeb3b', '#ff9800', '#795548', '#9e9e9e', '#607d8b']
    for i in range(0,14):
        y1 = param1[i,:]
        plt.plot(x, y1, color=colors_light[i], label=years[i])
    labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    plt.xticks(months, labels, rotation='45')
    plt.yticks()
    plt.xlabel('Month')
    plt.ylabel('%s' %label)
    plt.legend(bbox_to_anchor=(1.05, 1.01))
    plt.title('Difference from 2000-14 seasonal mean (200 km)')
    #plt.savefig('/nfs/a328/eecwk/waccm-x/figures/%s_200km_%s' %label %tscale, dpi=300)
    plt.show()

def plot_2d(param, label):
    months = np.arange(0,12,1)
    x, y = np.meshgrid(months, lats)
    plt.figure(figsize=(6,4))
    diffs = np.arange(700,1100,10)
    labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    plt.xticks(months, labels, rotation='45')
    plt.yticks(np.arange(-90,120,30))
    z = param[:,:]
    zT = z.transpose()
    plt.contourf(x[:,:], y[:,:], zT[:,:], diffs, cmap=plt.get_cmap('seismic'))
    cbar = plt.colorbar()
    cbar.set_label('%s' %label)
    plt.xlabel('Month')
    plt.ylabel('Latitude [%s]' %deg)
    plt.title('De-seasonalised monthly 2000-2014 means 200 km')
    #plt.savefig('/nfs/a328/eecwk/waccm-x/figures/%s_200km_%s' %label %tscale, dpi=300)
    plt.show()


p_all_mthly[:,:] = make_desolar_arrays('p')
T_all_mthly[:,:] = make_desolar_arrays('T')
n_all_mthly[:,:] = make_desolar_arrays('n')

p_indv_yrs[:,:,:] = slice_yrly_array(p_yrly)
T_indv_yrs[:,:,:] = slice_yrly_array(T_yrly)
n_indv_yrs[:,:,:] = slice_yrly_array(n_yrly)

p_indv_yrs_diff[:,:] = make_yrly_diff_arrays(p_all_mthly, p_indv_yrs)
T_indv_yrs_diff[:,:] = make_yrly_diff_arrays(T_all_mthly, T_indv_yrs)
n_indv_yrs_diff[:,:] = make_yrly_diff_arrays(n_all_mthly, n_indv_yrs)

plot_1d(T_all_mthly, T_indv_yrs_diff, '%sT' %delta)
#plot_2d(T_all_mthly, 'Temperature')