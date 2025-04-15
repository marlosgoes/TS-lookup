%EXAMPLE

clear
close all

load data_example.mat
nl   = 5000;      %RESAMPLE DATA WITH 5000 profiles
aa   = randi(length(time),nl,1);
S    = S(:,aa);
T    = T(:,aa);
time = time(aa);
latitude  = latitude(aa);
longitude = longitude(aa);

Pout = 0:2:6000;
Pad = 1;
method = 'Goes'
%CALC SALINITY
[S2,S3,TT,PP]=Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(T,P,latitude,longitude,time,Pad,Pout,method);

%SOME PLOTS

%TS PLOT
figure(1),clf 
 aa=plot(S(:,:),T(:,:),'r.','markersize',5);
 hold on
 bb=plot(S2(:,:),TT(:,:),'.b','markersize',5);
 aa = aa(1); bb = bb(1);
 legend([aa bb],'Original TS','Reconstructed TS','location','northwest')
figname = 'TS_orig_reconstruct_SA.jpg';
print('-djpeg','-r300',figname)

%SCATTER PLOT
%Choose a depth
dep = 5;% 10;
[~,dep2] = min(abs(P(dep)-Pout));

s1 = squeeze(S(dep,:,:));
s2 = squeeze(S2(dep2,:,:));
s1(isnan(s2))=nan;

figure(2),clf
plot(s1,s2(:),'.'),xlabel('original'),ylabel('prediction')
x_lim=xlim;y_lim=ylim;
xlim([min([x_lim y_lim]),max([x_lim y_lim])]);ylim(xlim);
figname = 'TS_orig_reconstruct_scatter_SA.jpg';
print('-djpeg','-r300',figname)

%MAP PLOT
%dep = 5;%5;% 100;
s1 = squeeze(S(dep,:,:));
s2 = squeeze(S2(dep2,:,:));
t1 = squeeze(T(dep,:,:));
cc = load('coast.mat');

figure(3),clf
subplot(2,1,2)
h2=color_line(longitude,latitude,s2,'.');
title(sprintf('Prediction: Depth = %g m',Pout(dep2))),h1=colorbar;
c_ax = caxis;
axis image;axis normal
ax1 = axis;
hold on
plot(cc.long,cc.lat,'k.')
box on
axis(ax1);

subplot(2,1,1)
h1=color_line(longitude,latitude,s1,'.');
title(sprintf('Original: Depth = %g m',P(dep))),h2=colorbar;
axis image;axis normal
caxis(c_ax)
hold on
plot(cc.long,cc.lat,'k.')
box on
axis(ax1);
figname = sprintf('TS_orig_reconstruct_map_%gm_scatter_SA.jpg',P(dep));
print('-djpeg','-r300',figname)
