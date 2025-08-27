"""
Simplified methane plume setup.
Profiles are slightly idealised, from 04:00 UTC ERA5 data.
Wind is rotated to be perfectly westerly.
"""

import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc
import os

pl.close('all'); pl.ion()

float_type = 'f8'

Mair = 28.9664
Mch4 = 16.04
T0 = 273.15
Rd = 287.04
Rv = 461.5
ep = Rd/Rv

def calc_esat(T):
    a_ab = 611.21; b_ab = 18.678; c_ab = 234.5; d_ab = 257.14
    return a_ab * np.exp((b_ab - ((T-T0) / c_ab)) * ((T-T0) / (d_ab + (T-T0))))

def calc_qsat(T, p):
    esat = calc_esat(T)
    return ep*esat / (p-(1.-ep)*esat)


# Get number of vertical levels and size from .ini file
with open('ch4_plumes.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])

# Vertical grid LES
dz = zsize / kmax
z = np.arange(0.5*dz, zsize, dz)

# idealised profiles of potential temperature
v_thl = np.array([285.7, 291.9, 293, 297.4, 307])
z_thl = np.array([0, 400, 2000, 2500, 5000])
thl = np.interp(z, z_thl, v_thl)

# idealised profiles of specific humidity
z_qt = np.array([0, 400, 2000, 2500, 5000])
v_qt = np.array([6.2, 4.93, 3.61, 1, 0.3])/1000
qt = np.interp(z, z_qt, v_qt)

# idealised profiles of zonal wind component
z_u = np.array([0, 270, 3000, 5000])
v_u = np.array([2.3, 8.5, 0.6, 5.7])
u = np.interp(z, z_u, v_u)

# no meridional wind and no ch4 background concentration
v = np.zeros(kmax)
ch4 = np.zeros(kmax)

# Surface fluxes, again idealised from ERA5.
t0 = 4*3600
t1 = 16*3600
td1 = 12*3600
td2 = 14*3600

time = np.linspace(t0, t1, 32)
wthl = 0.17   * np.sin(np.pi * (time-t0) / td1)
wqt  = 8.3e-5 * np.sin(np.pi * (time-t0) / td2)

# Write input NetCDF file
nc_file = nc.Dataset('ch4_plumes_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', kmax)
nc_z = nc_file.createVariable('z' , float_type, ('z'))

nc_group_init = nc_file.createGroup('init');
nc_u = nc_group_init.createVariable('u' , float_type, ('z'))
nc_v = nc_group_init.createVariable('v' , float_type, ('z'))
nc_th = nc_group_init.createVariable('thl', float_type, ('z'))
nc_qt = nc_group_init.createVariable('qt', float_type, ('z'))
nc_ch4 = nc_group_init.createVariable('ch4', float_type, ('z'))
nc_ch4_inflow = nc_group_init.createVariable('ch4_inflow', float_type, ('z'))

nc_z [:] = z[:]
nc_u [:] = u[:]
nc_v [:] = v[:]
nc_th[:] = thl[:]
nc_qt[:] = qt[:]
nc_ch4[:] = ch4[:]
nc_ch4_inflow[:] = ch4[:]

nc_group_tdep = nc_file.createGroup('timedep')
nc_group_tdep.createDimension("time_surface", time.size)
nc_time_surface = nc_group_tdep.createVariable("time_surface", float_type, ("time_surface"))
nc_thl_sbot = nc_group_tdep.createVariable("thl_sbot", float_type, ("time_surface"))
nc_qt_sbot = nc_group_tdep.createVariable("qt_sbot" , float_type, ("time_surface"))

nc_time_surface[:] = time
nc_thl_sbot[:] = wthl
nc_qt_sbot[:] = wqt

nc_file.close()

# Print .ini settings emissions:
# Coordinates of central cooling tower (m):
x0 = 1000
y0 = ysize/2.

# Std-dev of plume widths:
sigma_x = 0.5
sigma_y = 0.5


z0 = 10
sigma_z = 0.5

# Strength of plumes
# lower expected detection limit of AVIRIS-4: 10 kg/h / 3600 s/h / Mch4 kg/kmol
strength_ch4 = 10 / 3600 / Mch4   # kmol(CH4) s-1


print('sourcelist={}'.format('ch4'))

print('source_x0={}'.format(x0))
print('source_y0={}'.format(y0))
print('source_z0={}'.format(z0))

print('sigma_x={}'.format(sigma_x))
print('sigma_y={}'.format(sigma_y))
print('sigma_z={}'.format(sigma_z))

print('strength={}'.format(strength_ch4))
print('swvmr={}'.format('true'))

print('line_x={}'.format(0))
print('line_y={}'.format(0))
print('line_z={}'.format(0))
