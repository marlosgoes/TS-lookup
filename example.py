#example.py

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import subprocess
import Calc_sal_Thacker_Goes_EmDr_Stom_svd2
import Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe


# Load data
data = scipy.io.loadmat('data_example.mat')
S = data['S']
T = data['T']
P = data['P'].flatten()
time = data['time'].flatten()
latitude = data['latitude'].flatten()
longitude = data['longitude'].flatten()
acha = np.where(~np.isnan(T).all(axis=0))[0].tolist() # ~np.isnan(T)
S = S[:,acha]
T = T[:,acha]

latitude = latitude[acha]
longitude = longitude[acha]
time = time[acha]

# Resample data
nl = 5000 #100 
aa = np.random.randint(len(time), size=nl)
S = S[:, aa]
T = T[:, aa]
time = time[aa]
latitude = latitude[aa]
longitude = longitude[aa]
Pout = np.arange(0, 6001, 2)
Pad = 0
method = 'thacker'  #'annual' 'Thacker' #'svd' #stommel' #'Goes'

# Calculate salinity
S2, S3, TT, PP = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe.make_sal(T, P, latitude, longitude, time, Pad, Pout, method)  

# TS Plot
plt.figure(1)
plt.clf()
aa, = plt.plot(S.flatten(), T.flatten(), 'r.', markersize=5)
bb, = plt.plot(S2.flatten(), TT.flatten(), 'b.', markersize=5)
plt.legend([aa, bb], ['Original TS', 'Reconstructed TS'], loc='upper left')
plt.savefig('TS_orig_reconstruct_SA.jpg', dpi=300)
print(S.shape)

# Scatter Plot
dep = 5
dep2 = np.argmin(np.abs(Pout - P[dep]))
s1 = S[dep, :].squeeze()
s2 = S2[dep2, :].squeeze()
s1[np.isnan(s2)] = np.nan

plt.figure(2)
plt.clf()
plt.plot(s1, s2, '.')
plt.xlabel('original')
plt.ylabel('prediction')
x_lim = plt.xlim()
y_lim = plt.ylim()
plt.xlim([min(x_lim[0], y_lim[0]), max(x_lim[1], y_lim[1])])
plt.ylim(plt.xlim())
plt.savefig('TS_orig_reconstruct_scatter_SA.jpg', dpi=300)

# Map Plot
s1 = S[dep,  :].squeeze()
s2 = S2[dep2,  :].squeeze()
t1 = T[dep,  :].squeeze()
coast = scipy.io.loadmat('coast.mat')

plt.figure(3)
plt.clf()
plt.subplot(2, 1, 2)
h1 = plt.scatter(longitude, latitude, c = s2, cmap='plasma')
plt.title('Prediction: Depth =' + str(P[dep]) + 'm')
plt.colorbar(label='S[psu]');
plt.axis('image')
c_ax = h1.get_clim()
plt.plot(coast['long'], coast['lat'], 'k.')
plt.box(True)

plt.subplot(2, 1, 1)
plt.scatter(longitude, latitude, c = s1, cmap='plasma')
plt.title('Original: Depth =' + str(Pout[dep2]) + ' m')
plt.colorbar(label='S[psu]') 
plt.axis('image')
plt.clim(c_ax)
plt.plot(coast['long'], coast['lat'], 'k.')
plt.box(True)
plt.savefig(f'TS_orig_reconstruct_map_{P[dep]}m_scatter_SA.jpg', dpi=300)
#plt.savefig('TS_orig_reconstruct_map_m_scatter_SA.jpg', dpi=300)

plt.show()


