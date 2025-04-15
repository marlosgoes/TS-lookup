#load_woa13_pad2.py

import numpy as np
import netCDF4
from netCDF4 import Dataset
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def load_woa13_pad(latitude, longitude, month=None, modo='pad'):
    if month is None:
        month = 13
    if isinstance(month, int):
        month = np.full(len(latitude), month)
    print(month)    
    offset = 12
    mm0 = 12 + np.ceil(month / 3).astype(int)
    if modo == 'shallow':
        offset = 0
        mm0 = month
    
    num = f'{offset + 1:02d}'
    direc = 'WOA18/' #'/phodnet/data/WOA13/'
    preffix = '/decav/0.25/woa18_decav_'
    sfile = f'{direc}salt{preffix}s{num}_04.nc'
    print(sfile) 
    # Read latitude, longitude and depth from the netCDF file
    with Dataset(sfile, 'r') as ncid:
        ncid.set_auto_mask(False)
        lat = ncid.variables['lat'][:]
        lon = ncid.variables['lon'][:]
        lon[lon > 180] -= 360
        depth = ncid.variables['depth'][:]
    
    LONG_MIN = min(longitude) - 0.5
    LONG_MAX = max(longitude) + 0.5
    LAT_MIN = min(latitude) - 0.5
    LAT_MAX = max(latitude) + 0.5
    
    ilat = np.where((lat > LAT_MIN) & (lat < LAT_MAX))[0]
    ilon = np.where((lon > LONG_MIN) & (lon < LONG_MAX))[0]
    lon = lon[ilon]
    lat = lat[ilat]
    
    nz = len(depth)
    nlon = len(ilon)
    nlat = len(ilat)
    
    T_lev = np.full((nz, len(latitude)), np.nan)
    S_lev = np.full((nz, len(latitude)), np.nan)
    D_lev = depth
    Y1, X1 = np.meshgrid(lat, lon)
    
    Nm = np.unique(mm0)
    for mm in Nm:
        mon = f'{mm+1:02d}'
        ai = np.where(mm0 == mm)[0]
        latitude2 = latitude[ai]
        longitude2 = longitude[ai]
        
        sfile = f'{direc}salt{preffix}s{mon}_04.nc'
        tfile = f'{direc}temp{preffix}t{mon}_04.nc'
        print(sfile)
        # Load salinity
        with Dataset(sfile, 'r') as ncid:
            s_an = ncid.variables['s_an'][0,:,ilat,ilon].data
            miss = ncid.variables['s_an']._FillValue
            s_an[s_an == miss] = np.nan

        # Load temperature
        with Dataset(tfile, 'r') as ncid:
            t_an = ncid.variables['t_an'][0,:,ilat,ilon].data
            t_an[t_an == miss] = np.nan

        t_an = np.transpose(t_an, (0, 2, 1))  #(2, 0, 1))
        s_an = np.transpose(s_an, (0, 2, 1))  #2, 0, 1))
        mask = np.where(~np.isnan(t_an[0, :, :]), 1, np.nan)
        X = (X1 * mask).flatten()
        Y = (Y1 * mask).flatten()

        t_an = t_an.reshape(nz, nlon * nlat)
        s_an = s_an.reshape(nz, nlon * nlat)

        for ii in range(len(latitude2)):
            dist = np.abs(X - longitude2[ii]) + np.abs(Y - latitude2[ii])
            aa = np.nanargmin(dist)
            T_lev[:, ai[ii]] = t_an[:, aa]
            S_lev[:, ai[ii]] = s_an[:, aa]
 
    if modo == 'pad':
        Dbase = np.tile(D_lev[:, None], (1, len(latitude)))
        Dbase[np.isnan(T_lev)] = 0 #np.nan
        dmax, dvec = np.nanmax(Dbase, axis=0), np.nanargmax(Dbase, axis=0)
        
        # Dummy bathymetry data
      #  LONG, LAT = np.meshgrid(np.linspace(LONG_MIN, LONG_MAX, 100), np.linspace(LAT_MIN, LAT_MAX, 100))
      #  ELEV = np.random.random(LONG.shape) * -5000  # To be replaced with actual bathymetry data
      #  F = griddata((LONG.flatten(), LAT.flatten()), ELEV.flatten(), (longitude, latitude), method='cubic')
        if 0:
         F = dmax[:,np.newaxis].T
         D_EXT = np.arange(5600, 6501, 50)
         D_lev = np.concatenate((D_lev, D_EXT))
         nz2 = len(D_lev)
         bath = np.tile(F, (nz2, 1))

         Tdep = ([T_lev[dvec[kk],kk] for kk in range(len(latitude))])
         Sdep = ([S_lev[dvec[kk],kk] for kk in range(len(latitude))])
         Tdep2 = np.nan * np.ones((1,len(latitude))) # 
         Sdep2 = Tdep2
         Tdep2[:] = Tdep[:]
         Sdep2[:] = Sdep[:]
         for jj in range(nz2-nz):
                #tap = T_lev[dvec[kk],kk]
                T_lev = np.concatenate((T_lev, Tdep2), axis=0)
                S_lev = np.concatenate((S_lev, Sdep2), axis=0)

         mask = np.tile(D_lev[:, None],(1,len(latitude))) > bath
         T_lev[mask] = np.nan
         S_lev[mask] = np.nan

#        print(S_lev.shape)
#        print(T_lev.shape)
#        print(D_lev.shape) 
#        print(len(latitude))
   
#        X, Y = np.meshgrid(np.arange(0,len(latitude),1),D_lev) 
#        plt.figure(5)
#        plt.clf()
#        plt.contourf(X, Y, S_lev, 20, cmap='viridis')
#        plt.show() 
    return T_lev, S_lev, D_lev
