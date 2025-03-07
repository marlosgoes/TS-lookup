%EXAMPLE
addpath /data/mgoes/matlab/m_map/

clear
close all

filename = '/data/mgoes/ARGO/IPRC/argo_2005-2015_grd.nc';
TIME      = ncread(filename,'TIME'); 
longitude = ncread(filename,'LONGITUDE');
latitude  = ncread(filename,'LATITUDE');
P     = ncread(filename,'LEVEL');
longitude(longitude>180)=longitude(longitude>180)-360;

   T      = ncread(filename,'TEMP');
   T     = permute(T,[3 2 1 4]);
  mask = squeeze(mean(T(1,:,:,:),4));mask = ~isnan(mask);
   T2 = T(:,:,:);%,100);
   T2 = T2(:,mask);
[X,Y] = meshgrid(longitude,latitude);


clear time
years = 2005:2015;
count = 0;
for yy=1:length(years)
 for ii=1:12;
   count = count+1;
   mon = sprintf('0%i',ii);mon = mon(end-1:end);
   time(count,:) = sprintf('%i%s%i',years(yy),mon,15);
 end
end
time = str2num(time);
time = time(1:length(TIME));
nl   = length(T2);      %RESAMPLE DATA WITH 5000 profiles
time = time(1)*ones(1,nl);
X = X(mask);Y=Y(mask);

keyboard
%load data_example.mat
%aa   = 1:nl; %randi(length(time),nl,1);
%S    = S(:,aa);
%T    = T(:,aa);
%time = time(aa);
%latitude  = latitude(aa);
%longitude = longitude(aa);

Pout = [0:2:10 100:500:5000]; %:2:6000;
Pad = 0;
method = 'annual'%Goes'
%CALC SALINITY

[S2,S3,TT,PP]=Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(T2,P,Y,X,time,Pad,Pout,method);

S4 = T*nan;S4 = S4(1:length(Pout),:,:,1);
S4(:,mask) = S2;

%SOME PLOTS
longitude(longitude<0) = longitude(longitude<0)+360;
figure;clf
set(gcf,'PaperPosition',[1 1 7 5]);
set(gcf,'DefaultLineLineWidth',1.5)
wysiwyg;
%subplot('position',[.15 .45 .8 .5])
m_proj('miller','long',[0 360],'lat',[-60 70])
m_pcolor(longitude,latitude,squeeze(S4(1,:,:)));
shading flat
%caxis([-.5 .5])
colormap(jet(17))
m_grid('box','fancy')
m_coast('patch',[.5 .5 .5]);



stop

%TS PLOT
figure(1),clf 
 aa=plot(S(:,:),T(:,:),'r.','markersize',5);
 hold on
 bb=plot(S2(:,:),TT(:,:),'.b','markersize',5);
 aa = aa(1); bb = bb(1);
 legend([aa bb],'Original TS','Reconstructed TS','location','northwest')
figname = 'TS_orig_reconstruct_SA.jpg';
%print('-djpeg','-r300',figname)

%SCATTER PLOT
%Choose a depth
dep = 5;% 10;
s1 = squeeze(S(dep,:,:));
s2 = squeeze(S2(dep,:,:));
s1(isnan(s2))=nan;

figure(2),clf
plot(s1,s2(:),'.'),xlabel('original'),ylabel('prediction')
x_lim=xlim;y_lim=ylim;
xlim([min([x_lim y_lim]),max([x_lim y_lim])]);ylim(xlim);
figname = 'TS_orig_reconstruct_scatter_SA.jpg';
%print('-djpeg','-r300',figname)

%MAP PLOT
dep = 5;%5;% 100;
s1 = squeeze(S(dep,:,:));
s2 = squeeze(S2(dep,:,:));
t1 = squeeze(T(dep,:,:));
cc = load('coast.mat');

figure(3),clf
subplot(2,1,2)
h2=color_line(longitude,latitude,s2,'.');
title(sprintf('Prediction: Depth = %g m',P(dep))),h1=colorbar;
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
%print('-djpeg','-r300',figname)
