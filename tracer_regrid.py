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
lons = np.zeros([144])
p_int = np.zeros([96,144])
T_int = np.zeros([96,144])
n = np.zeros([96,144])
p_yrly = np.zeros([115,96,144])
T_yrly = np.zeros([115,96,144])
n_yrly = np.zeros([115,96,144])

fname_uni = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.2000-01.nc', 'r', format='NETCDF4')
levs = fname_uni.variables['lev'][:]
lats = fname_uni.variables['lat'][:]
lons = fname_uni.variables['lon'][:]
fname_uni.close()

def add_months(sourcedate, months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month // 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year, month)[1])
    return date(year, month, day)

def interp_p(altitude, gpheight):
    for i in range(0,96):
        for j in range(0,144):
            p_int[i,j] = np.interp(altitude, gpheight[:,i,j][::-1], levs[:][::-1])     
    return p_int

def interp_T(pressure_int, temp):
    for i in range(0,96):
        for j in range(0,144):
            T_int[i,j] = np.interp(pressure_int[i,j], levs[:][::-1], temp[:,i,j][::-1])
    return T_int

def calc_n(pressure_int, temp_int):
    for i in range(0,96):
        for j in range(0,144):
            n[i,j] = (pressure_int[i,j] * 100) / (temp_int[i,j] * k_B)
    return n

def calc_diff(param, start, end):
    p1 = param[start]
    p2 = param[end]       
    diff = ( (p2 - p1) / p1 ) * 100
    return diff

def make_diff_arrays(tscale, output):
    z3_dat = np.zeros([1,145,96,144])
    z3 = np.zeros([145,96,144])
    T = np.zeros([96,144])
    start = date(2000,1,1)
    date_list = []
    fdates = []
    ndates = np.zeros(115)

    for i in range(0,115):
        date_i = add_months(start, i)
        date_list.append(date_i)
        ndate = matplotlib.dates.date2num(date_i)
        ndates[i] = ndate
        fdates.append(date_list[i].strftime('%Y-%m'))

    for i in range(0,115):
        fname = netCDF4.Dataset('/nfs/a328/eecwk/earth_system_grid/ccsm4_monthly_ave/f.e20.FXSD.f19_f19.001.cam.h0.%s.nc' %fdates[i], 'r', format='NETCDF4')
        z3_dat = fname.variables['Z3'][:]
        z3 = np.mean(z3_dat, axis=0)
        T_dat = fname.variables['T'][:]
        T = np.mean(T_dat, axis=0)
        fname.close()            
        p_yrly[i,:,:] = interp_p(200000, z3)
        T_yrly[i,:,:] = interp_T(p_int, T)
        n_yrly[i,:,:] = calc_n(p_int, T_int)

    if tscale == 'seasonal':
        start = 0
        end = 6
    if tscale == 'solar':
        start = 0
        end = 108
    if output == 'T':
        diff = calc_diff(T_yrly, start, end)
    if output == 'n':
        diff = calc_diff(n_yrly, start, end)
    return diff

def plot(param, label, tscale):
    x, y = np.meshgrid(lons, lats)
    plt.figure(figsize=(6,4))
    diffs = np.arange(-100,101,1)
    plt.xticks(np.arange(0,420,60))
    plt.yticks(np.arange(-90,120,30))
    z = param[:,:]
    plt.contourf(x[:,:], y[:,:], z[:,:], diffs, cmap=plt.get_cmap('seismic'))
    cbar = plt.colorbar(ticks=np.arange(-100,125,25))
    cbar.set_label('%s difference' %label)
    plt.xlabel('Longitude [%s]' %deg)
    plt.ylabel('Latitude [%s]' %deg)
    plt.title('200 km %s' %tscale)
    #plt.savefig('/nfs/a328/eecwk/waccm-x/figures/%s_200km_%s' %label %tscale, dpi=300)
    plt.show()

tscale = ['seasonal', 'solar']

for i in range(0, len(tscale)):
    T_diff = make_diff_arrays(tscale[i], 'T')
    n_diff = make_diff_arrays(tscale[i], 'n')
    plot(T_diff, 'Temperature', tscale[i])
    plot(n_diff, 'Density', tscale[i]);