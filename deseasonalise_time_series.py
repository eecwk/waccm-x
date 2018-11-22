import matplotlib
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from datetime import date
import calendar

deg = unichr(176)
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

def calc_mthly_means(param_yrly, month):
    start = month
    end = (13 * 12) + month
    step = np.arange(start, end, 12)
    param_mthly = np.zeros([len(step),96])
    for i in range(0,len(step)):
        year = step[i]
        param_mthly[i,:] = param_yrly[year,:]
    param_mthly_mean = np.mean(param_mthly, axis=0)
    return param_mthly_mean

def make_deseason_arrays(output):
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
        p_all_mthly[i] = calc_mthly_means(p_yrly, i)
        T_all_mthly[i] = calc_mthly_means(T_yrly, i)
        n_all_mthly[i] = calc_mthly_means(n_yrly, i)    
    for i in range(0,12):
        if output == 'p':
            param_all_mthly[i] = calc_mthly_means(p_yrly, i)
        if output == 'T':
            param_all_mthly[i] = calc_mthly_means(T_yrly, i)
        if output == 'n':
            param_all_mthly[i] = calc_mthly_means(n_yrly, i)
    return param_all_mthly

def plot(param, label):
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

p_all_mthly = make_deseason_arrays('p')
T_all_mthly = make_deseason_arrays('T')
n_all_mthly = make_deseason_arrays('n')

plot(T_all_mthly, 'Temperature')

'''
days = np.arange(0,365,1)
mth_mids = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349]
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
z = np.zeros([365,96])
for i in range(0,365):
    if i < 31:
        z[i,:] = T_all_mthly[0,:]
    if 31 <= i < 59:
        z[i,:] = T_all_mthly[1,:]
    if 59 <= i < 90:
        z[i,:] = T_all_mthly[2,:]    
    if 90 <= i < 120:
        z[i,:] = T_all_mthly[3,:]  
    if 120 <= i <151:
        z[i,:] = T_all_mthly[4,:]
    if 151 <= i < 181:
        z[i,:] = T_all_mthly[5,:]
    if 181 <= i < 212:
        z[i,:] = T_all_mthly[6,:]    
    if 212 <= i < 243:
        z[i,:] = T_all_mthly[7,:]  
    if 243 <= i < 273:
        z[i,:] = T_all_mthly[6,:]
    if 273 <= i < 304:
        z[i,:] = T_all_mthly[9,:]
    if 304 <= i < 334:
        z[i,:] = T_all_mthly[10,:]    
    if 334 <= i < 365:
        z[i,:] = T_all_mthly[11,:]

x, y = np.meshgrid(days, lats)
diffs = np.arange(700,1100,10)
plt.figure(figsize=(6,4))
plt.xticks(mth_mids, labels, rotation='45')
plt.yticks(np.arange(-90,120,30))
zT = z.transpose()
plt.contourf(x[:,:], y[:,:], zT[:,:], diffs, cmap=plt.get_cmap('inferno'))
#cbar = plt.colorbar(ticks=np.arange(-100,125,25))
cbar = plt.colorbar()
cbar.set_label('Temperature [K]')
plt.xlabel('Month')
plt.ylabel('Latitude [%s]' %deg)
plt.title('De-seasonalised monthly 2000-2014 means 200 km')
#plt.savefig('/nfs/a328/eecwk/waccm-x/figures/%s_200km_%s' %label %tscale, dpi=300)
plt.show()
'''