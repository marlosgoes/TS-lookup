#April 14, 2025
#By: Marlos Goes
#from MATLAB function [y2,y3,TT,PP]=Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(T,P,lat,lon,time,pad,Pout,method)
#TS lookup for the Atlantic basin
#INPUT:  T: Temperature (degC) 
#        P: Depth       (m, positive)
#        Lat:   (Degree north)
#        Lon:   (Degree East [-100:20])
#        time:   E.g. 200410 (yyymm)
#        pad:    0 [Default] or 1
#        Pout: Depth output (m, positive)
#        method: 'Thacker'[Default] / 'Goes' / 'annual' / 'stommel' / 'svd'
#OUTPUT: y2: Salinity from the TS method chosen
#        y3: smoothed version of y2 
#        TT: Temperature output using Pout as vertical axis if defined
#        PP: Depth defined for the profiles


import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.interpolate import interp1d
import load_woa13_pad2
import array

def make_sal(TT, PP, lat, lon, timet, pad=None, Pout=None, method='thacker'):
    if pad is None:
        pad = False
    if Pout is None:
        IntP = False
    else:
        IntP = True
    if method is None:
        method = 'thacker'
    else:
        method = method.lower()
    if method not in ['goes', 'thacker', 'annual', 'stommel', 'svd']:
        raise ValueError('MethodNotFound')
    Dref = 800
    if pad > 1:
        Dref = pad
        pad = 1

    lon2 = lon.copy()
    lon2[lon2 < 0] += 360

    Z = np.concatenate((np.arange(7.5, 301, 10), np.arange(320, 1001, 20), np.arange(1050, 2001, 50)))

    if len(timet) == 1:
        timet = np.full(lat.shape, timet)
    timet = timet.flatten()
    timet[np.isnan(timet)] = 0
    
    basin = np.empty(lat.shape, dtype='<U3')

    sai = (lat >= -45) & (lat < -15) & (lon > -70) & (lon < 20)
    basin[sai] = 'SA'

    tai = (lat >= -15) & (lat < 10) & (lon > -80) & (lon < 20)
    basin[tai] = 'TA'

    nai = (lat >= 10) & (lat < 45) & (lon > -98) & (lon < 0)
    basin[nai] = 'NA'

    fci = (lat >= 26.5) & (lat < 27.2) & (lon > -80.2) & (lon < -78.9)
    basin[fci] = 'FC'
    nai[fci] = False

    mdi = (lat >= 30) & (lat <= 45) & (lon > -6) & (lon <= 30)
    basin[mdi] = 'MD'

    noi = (lat >= 45) & (lat < 65) & (lon > -100) & (lon < 30)
    basin[noi] = 'NO'

    soi = (lat >= -70) & (lat < -45) & (lon > -70) & (lon <= 20)
    basin[soi] = 'SO'

    tii = (lat >= -45) & (lat < 30) & (lon > 20) & (lon < 120)
    basin[tii] = 'TI'

    soii = (lat >= -75) & (lat < -45) & (lon > 20) & (lon < 120)
    basin[soii] = 'SOI'

    spi = (lat >= -45) & (lat < -15) & (lon2 > 120) & (lon2 < 290)
    basin[spi] = 'SP'

    tpi = (lat >= -15) & (lat < 15) & (lon2 > 120) & (lon2 < 285)
    basin[tpi] = 'TP'

    npi = (lat >= 15) & (lat < 45) & (lon2 > 120) & (lon2 < 260)
    basin[npi] = 'NP'

    nopi = (lat >= 45) & (lat < 65) & (lon2 > 120) & (lon2 < 260)
    basin[nopi] = 'NOP'

    sopi = (lat >= -75) & (lat < -45) & (lon2 > 120) & (lon2 <= 290)
    basin[sopi] = 'SOP'

    n0, nla1 = TT.shape   
    nmo1 = 1
    if len(PP) == PP.size:
        PP = np.tile(PP[:, np.newaxis], (1, nla1)) #, nmo1))
       # P = np.tile(P.reshape(-1,1,1),(1, nla1, nmo1))

    P2 = PP.copy()
    T2 = TT.copy()
    T_interp = np.nan* np.ones((85, nla1))
    for ii in range(nla1):
            t = TT[:, ii]
            p = PP[:, ii]
            acha = ~np.isnan(t + p)
            if not np.any(acha):
                continue
            if np.count_nonzero(acha) > 1:
                T_interp[:, ii] = interp1d(p[acha], t[acha], bounds_error=False, fill_value=np.nan)(Z)
            elif np.count_nonzero(acha) == 1 or np.max(p[acha]) < Z[0]:
                imin = np.argmin(np.abs(Z - p[acha]))
                aux = np.nan * np.ones(Z.shape)
                aux[imin] = t[acha]
                T_interp[:, ii] = aux
    TT0 = TT
    TT = T_interp
    nz, nla1 = TT.shape
    if len(PP) == PP.size:
        PP = Z
    else:
        PP = np.tile(Z[:, np.newaxis], (1, nla1))

    nt = nla1  # * nmo1
    y2 = np.nan * np.ones((nz, nla1))  #, nmo1))
    stnai = np.nonzero([tai.any(), sai.any(), nai.any(), soi.any(), noi.any(), tii.any(), soii.any(), nopi.any(), npi.any(), tpi.any(), spi.any(), sopi.any(), fci.any(), mdi.any()])[0]
    if method == 'svd':
        preffix = 'fit_sigma_cora_argo_extended_svd_'
        suffix = '.mat'
    else:
        preffix = 'fit_sigma_cora_argo_extended_fitlm_'
        suffix = '_thackergoes_noyear.mat'
    
    for bb in stnai:

        if bb == 0:
            filename = preffix+'TA'+suffix
            ai = tai
        elif bb == 1:
            filename = preffix+'SA'+suffix
            ai = sai
        elif bb == 2:
            filename = preffix+'NA'+suffix
            ai = nai
        elif bb == 3:
            filename = preffix+'SO'+suffix
            ai = soi
        elif bb == 4:
            filename = preffix+'NO'+suffix
            ai = noi
        elif bb == 5:
            filename = preffix+'TI'+suffix
            ai = tii
        elif bb == 6:
            filename = preffix+'SOI'+suffix
            ai = soii
        elif bb == 7:
            filename = preffix+'NOP'+suffix
            ai = nopi
        elif bb == 8:
            filename = preffix+'NP'+suffix
            ai = npi
        elif bb == 9:
            filename = preffix+'TP'+suffix
            ai = tpi
        elif bb == 10:
            filename = preffix+'SP'+suffix
            ai = spi
        elif bb == 11:
            filename = preffix+'SOP'+suffix
            ai = sopi
        elif bb == 12:
            filename = preffix+'FC'+suffix[:-4]+'_ryanonly.mat'
            ai = fci
        elif bb == 13:
            filename = preffix+"MD2"+suffix
            ai = mdi

        print(filename)
        data = loadmat(filename)
        xmean_ii = data['xmean_ii']
        ymean_ii = data['ymean_ii']
        X = data['X']
        Y = data['Y']
        T = TT[:, ai]
        P = PP[:, ai]

        if (bb < 8 or bb>11):
            longitude = lon[ai]
        else:
            longitude = lon2[ai]

        latitude = lat[ai]
        time = timet[ai]
        alnan = time == 0

        nt = len(time)
        time = time.astype(str).T
        year = np.nan * np.ones((nt)) #, nmo1))
        month = year.copy()          #month = np.full(nt, np.nan)
        auxy = time[~alnan]
        auxy2=([int(s[:4]) for s in auxy])
        year[~alnan] = np.clip(auxy2, 1993, 2015) 
        auxy2=([int(s[5:6]) for s in auxy])
#        year[~alnan] = np.clip(time[~alnan].astype(int)[:, 0:4], 1993, 2015)
#        month[~alnan] = time[~alnan].astype(int)[:, 4:6]
        month[~alnan] = auxy2 

        nz, nla = T.shape

        yy = np.tile(year, (nz, 1))
        mo = np.tile(month, (nz, 1))
        pp = P.copy()

      #  T = T.reshape(nz, nla * nmo)

        pred3 = np.nan * np.ones(T.shape)
        aan = np.nonzero(~alnan)[0]
        for aa in aan:
            la  = latitude[aa]
            lo = longitude[aa]
            dflat = np.abs(latitude[aa] - Y) + np.abs(longitude[aa] - X)
            ll = np.argmin(dflat)
            xmean = xmean_ii[ll][0]
            
            if dflat[ll] > 2 or xmean.size == 0:
                 continue

            ymean = ymean_ii[ll][0]
            t0 = T[:, aa]
            t2 = t0[:,np.newaxis] - xmean #.astype(np.float64) #[np.newaxis,:])
            if method == 'svd':
               beta = data['CS_ii'][ll,:] #[0]
               U = data['U_ii'][ll][0]
               V = data['V_ii'][ll][0]
               t2[np.isnan(t2)] = 0
               yn = 0
               for ee in range(15):
                   Cs = beta[ee]
                   aux1 = U[:, ee] #.reshape(1,1)
                   aux2 = V[:, ee] #.reshape(1,1)
                   UV = np.matmul(aux1[:,np.newaxis],aux2[:,np.newaxis].T)
                   #yn += Cs * (U[:, ee] @ V[:, ee].T @ t2)
                   #yn += Cs * (UV @ t2) + ymean
                   yn += Cs * (UV @ t2) 
               yn += ymean    
            else:
                beta_ii = data['beta_ii']
                beta1 = beta_ii[ll, 0]
                beta2 = beta_ii[ll, 1]
                lat_ii = np.tile(latitude[aa,np.newaxis],(1,85))
                lon_ii = np.tile(longitude[aa,np.newaxis],(1,85))
                
                x0 = X[ll]
                y0 = Y[ll]
                x12 = np.column_stack([np.ones(85), t2, t2 ** 2, np.cos(2 * np.pi * mo[:, aa] / 12),
                    np.sin(2 * np.pi * mo[:, aa] / 12),
                    np.cos(4 * np.pi * mo[:, aa] / 12),
                    np.sin(4 * np.pi * mo[:, aa] / 12),
                    lat_ii.T - y0, lon_ii.T - x0])

                betam = 0*beta1    #annual
                if method == 'goes':
                    betam = beta2
                elif method == 'thacker':
                    betam = beta1

                if method == 'stommel':
                    nnan = ~np.isnan(t0)
                    tnan = ~np.isnan(ymean)
                    t3 = np.clip(t0[nnan], np.min(xmean), np.max(xmean))
                    try:
                       pred3[:np.count_nonzero(nnan), aa] = interp1d(xmean[tnan], ymean[tnan], bounds_error=False, fill_value=np.nan)(t3)
                    except:
                          print('except')
                          xmean, xunic = np.unique(xmean, return_index=True)
                          ymean = ymean[xunic]
                          tnan = ~np.isnan(ymean) & ~np.isnan(xmean)
                          pred3[:np.count_nonzero(nnan), aa] = interp1d(xmean[tnan], ymean[tnan], bounds_error=False, fill_value=np.nan)(t3)
                else:
                   for kk in range(nz):
                       aw = np.nonzero(pp[:, aa] == Z[kk])[0]
                       pred3[aw, aa] = np.matmul(x12[aw,:],betam[kk,:].T) + ymean[kk]
                        
                yn = pred3[:, aa]

       #END method 
            aza = np.nonzero((yn > 42) | (yn < 5))[0]
            if aza.size > 0:
               print('change to clim')

            yn[aza] = np.nan
            pred3[:, aa] = yn.flatten()
       #END method #AAN

        y2[:, ai] = pred3
       #END BB

    nz = P2.shape[0]
    PP = P2[:nz, :]
    nlon = y2.shape[1]

#  INTERPOLATE BACK TO ORIGINAL
    Y = np.nan * np.ones((nz, nlon))
    for ii in range(nlon):
        y = y2[:, ii]
        p = PP[:, ii]
        acha = ~np.isnan(y)
        acha2 = ~np.isnan(p)
        if np.count_nonzero(acha) > 1:
            f = interp1d(Z[acha], y[acha], bounds_error=False, fill_value=np.nan)
            Y[acha2, ii] = f(p[acha2])
        elif np.count_nonzero(acha) == 1:
            imin = np.argmin(np.abs(Z[acha] - p[acha2]))
            aux = np.full(p.shape, np.nan)
            aux[imin] = y[acha]
            Y[acha2, ii] = aux[acha2]

    y2 = Y
    TT = T2.copy()

#SLAB APPROXIMATION
    nlon = y2.shape[1]
    for jj in range(nlon):
        s2 = y2[:,jj]
        p = PP[:,jj]
        if any(~np.isnan(s2)):
           facha1 = np.where(~np.isnan(s2))[0][0]
           if facha1 > 0:
              s2[:facha1] = s2[facha1]
           acha2 = ~np.isnan(s2+p)
           if np.count_nonzero(acha2) > 2 and np.count_nonzero(acha2) < nz:
              f = interp1d(p[acha2], s2[acha2], bounds_error=False, fill_value=np.nan)
              s2[~acha2] = f(p[~acha2])
              s2[np.isnan(p)] = np.nan

        y2[:,jj] = s2


# REDO TIME ARRAY
    alnan = timet == 0
    nt = len(timet)
    time = timet.astype(str).T
    month = np.nan * np.ones((nt)) #, nmo1))
    auxy = time[~alnan]
    auxy2=([int(s[5:6]) for s in auxy])
    month[~alnan] = auxy2


#TESTE NEW PADDING
    if pad and 1:
        T_lev, S_lev, D_lev = load_woa13_pad2.load_woa13_pad(lat, lon, month)
        Pmax = np.nanmax(PP.flatten())
        P = np.unique(np.round(PP.flatten()))
        P = P[~np.isnan(P)]
        Dmax = np.argmin(np.abs(P - Dref))
        P = P[:Dmax]
        yi = np.nan * np.ones((Dmax, nt))
        Ti = yi.copy()
        
        for jj in range(nt):
            nonan = ~np.isnan(TT[:, jj])
            if np.count_nonzero(nonan) > 1:
                Ti[:Dmax, jj] = interp1d(PP[nonan,jj], TT[nonan,jj], bounds_error=False, fill_value=np.nan)(P)
                yi[:Dmax, jj] = interp1d(PP[nonan,jj], y2[nonan,jj], bounds_error=False, fill_value=np.nan)(P)
        y2 = yi
        T = Ti
        Dmin = D_lev > Dref
        Z_I = np.concatenate((P, D_lev[Dmin])) 
        nnanz = ~np.isnan(Z_I)
        Z_I = Z_I[nnanz]
        P = Z_I
        PP = np.tile(Z_I[:, np.newaxis], (1, nt))
        T_lev2 = np.nan * np.ones(PP.shape)
        S_lev2 = T_lev2.copy()
        for jj in range(nt):
             nonan = ~np.isnan(T_lev[:, jj])
             f = interp1d(D_lev[nonan], T_lev[nonan,jj], bounds_error=False,kind='nearest', fill_value='extrapolate')
             T_lev2[:,jj] = f(Z_I) #[nnanz])   
             f = interp1d(D_lev[nonan], S_lev[nonan,jj], bounds_error=False,kind='nearest', fill_value='extrapolate')
             S_lev2[:,jj] = f(Z_I)

        zrange = np.arange(Dmax, len(Z_I))
        for jj in zrange:
            ari = np.nan * np.ones((1,nt))
            y2 = np.concatenate((y2, ari),axis=0)
            T = np.concatenate((T, ari),axis=0)

        ddnan = np.isnan(y2) | (PP > Dref)
        y2[ddnan] = S_lev2[ddnan]
        ddnan = np.isnan(T) | (PP > Dref)
        T[ddnan] = T_lev2[ddnan]

    if not IntP:
        Pout = P
#Smoothing Filter
    nz = Pout.shape[0]
    nlon = y2.shape[1]
    y3 = np.nan * np.ones((nz,nlon)) #_like(y2)
    y4 = np.nan * np.ones((nz,nlon)) #_like(y2)
    T2 = np.nan * np.ones((nz,nlon)) #_like(y3)
    print('Interpolation and Smoothing')
    for jj in range(nlon):
        s3 = y2[:, jj]
        t3 = TT[:, jj]
        P = PP[:, jj]
        facha2 = np.nonzero((~np.isnan(s3)) & (s3 > 0))[0]
        if facha2.size < 2:
            continue

        if IntP:
 
            s3 = interp1d(P[facha2], s3[facha2], bounds_error=False, fill_value=np.nan)(Pout)
            t3 = interp1d(P[facha2], t3[facha2], bounds_error=False, fill_value=np.nan)(Pout)
            facha2 = np.nonzero(~np.isnan(s3))[0]
            if facha2.size == 0:
                continue

        tfilt2 = s3.copy()
        tfilt2[tfilt2 == 0] = np.nan

        if 1:   #nargout > 1:
            eeo = tfilt2.copy()
            dif1 = eeo[1:] - eeo[:-1]
            oi = dif1.copy()
            oi[np.isnan(oi)] = 0
            oi=oi[:-1]
             
            epsl = 0.1
            nei = np.abs(s3[2:] - s3[1:-1] - (s3[1:-1] - s3[:-2])) > epsl
            s3_diff = np.abs(np.diff(s3, n=2))
            nei = s3_diff > epsl
            f = interp1d(np.nonzero(~nei)[0], oi[~nei], bounds_error=False, fill_value=np.nan)
            oi[nei] = f(np.nonzero(nei)[0])
            oi[np.isnan(oi)] = 0

            e0 = s3[facha2[-1]]
            oii = np.append(-oi, e0, axis=None) 

            eeoi = np.flip(np.cumsum(np.flip(oii)))
            tfilt2[:len(facha2)] = eeoi[:len(facha2)]

        y3[:tfilt2.size, jj] = tfilt2
        T2[:tfilt2.size, jj] = t3[:tfilt2.size]
        y4[:tfilt2.size, jj] = s3

    TT = T2
    PP = Pout
    y3 = y3[:nz, :nla1]
    y2 = y4[:nz, :nla1]


    return y2, y3, TT, PP
