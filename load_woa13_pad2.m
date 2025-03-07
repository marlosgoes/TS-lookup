%USES M_MAP TOOLBOX
%FIXED SALINITY ncid name Apr 2017
%THIS FUNCTION INTERPOLATES LEVITUS WOA13 to a section
%USE:

%[T_lev,S_lev,D_lev]=load_woa13_pad(latitude,longitude,month,modo)
%latitude, longitude: vectors of lat,lon (1xm)
%month: month number (1 or 1xm) :ex: 2
%modo: 'pad'[Default] or 'shallow'; Uses deep seasonal (0:5000) or 
%        shallow monthly(0:1500)

function [T_lev,S_lev,D_lev]=load_woa13_pad(latitude,longitude,month,modo)

%addpath /data/mgoes/matlab/netcdf/
%addpath /data/mgoes/matlab/netcdf/mexnc
%addpath /data/mgoes/matlab/netcdf/ncsource
%addpath /data/mgoes/matlab/netcdf/nctype
%addpath /data/mgoes/matlab/netcdf/ncutility

if nargin < 4
    modo = 'pad';
end
disp(modo)

if nargin < 3
    month = 13;
end

if length(month)==1;
    month = month*ones(1,length(latitude));
end

offset = 12;
mm0 = 12 + ceil(month/3);
if strcmp(modo,'shallow')
    offset=0;
    mm0 = month;
end
num = sprintf('0%i',offset + 1);num = num(end-1:end);

%
%mm0 =  offset + ceil(month/3);

direc = '/phodnet/data/WOA13/';
preffix = '/decav/0.25/woa13_decav_';
sfile = sprintf('%s%s%s%s%s%s',direc,'salt',preffix,'s',num,'_04.nc');

%sfile = sprintf('%s%s%s%s%i%s',direc,'salt',preffix,'s',13,'_04.nc');
%tfile = sprintf('%s%s%s%s%i%s',direc,'temp',preffix,'t',13,'_04.nc');

ncid = netcdf.open(sfile,'NC_NOWRITE');

varid = netcdf.inqVarID(ncid,'lat');
  lat = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'lon');
  lon = netcdf.getVar(ncid,varid);    lon(lon>180)=lon(lon>180)-360;
varid = netcdf.inqVarID(ncid,'depth');
  depth = netcdf.getVar(ncid,varid);
netcdf.close(ncid);

%FROM THE DATA
    LONG_MIN = min(longitude)  -0.5;
    LONG_MAX = max(longitude)  +0.5;
    LAT_MIN  = min(latitude) -0.5;
    LAT_MAX  = max(latitude) +0.5;
    
ilat = find(lat>LAT_MIN & lat < LAT_MAX);
ilon = find(lon>LONG_MIN & lon < LONG_MAX);

%ilat = find(lat>min(latitude)-0.5 & lat < max(latitude)+0.5);
%ilon = find(lon>min(longitude)-0.5 & lon < max(longitude)+0.5);

%South Atlantic region
%ilat = find(lat>-50 & lat < -10);
%ilon = find(lon>-70 & lon < 25);

lon = lon(ilon);
lat = lat(ilat);

%Define size
nz = length(depth);
nlon = length(ilon);
nlat = length(ilat);

%OUTPUT VARIABLES
T_lev = nan*zeros(nz,length(latitude));
S_lev=T_lev;
D_lev = depth;
[Y1,X1] = meshgrid(lat,lon);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nm = unique(mm0); 
for jj = 1:length(Nm)
    mm = Nm(jj);
    mon = sprintf('0%i',mm);mon = mon(end-1:end);
    ai = find(mm0==mm);
    latitude2 = latitude(ai);
    longitude2 = longitude(ai);
    
%sfile = sprintf('%s%s%s%s%i%s',direc,'salt',preffix,'s',mm,'_04.nc');
%tfile = sprintf('%s%s%s%s%i%s',direc,'temp',preffix,'t',mm,'_04.nc');
sfile = sprintf('%s%s%s%s%s%s',direc,'salt',preffix,'s',mon,'_04.nc');
tfile = sprintf('%s%s%s%s%s%s',direc,'temp',preffix,'t',mon,'_04.nc');

disp('load salinity...')
ncid = netcdf.open(sfile,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'s_an');  
 %s_an = netcdf.getVar(ncid,varid,[ilat(1)-1 ilon(1)-1 1 1],[nlat nlon nz 1]); 
s_an = ncread(sfile,'s_an',[ilon(1) ilat(1) 1 1],[nlon nlat nz 1]);
miss = netcdf.getAtt(ncid,varid,'_FillValue');
mask = s_an==miss;
s_an(mask)=nan;
netcdf.close(ncid);

disp('load temperature...')
ncid = netcdf.open(tfile,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'t_an');  
t_an = ncread(tfile,'t_an',[ilon(1) ilat(1) 1 1],[nlon nlat nz 1]);
miss = netcdf.getAtt(ncid,varid,'_FillValue');
t_an(mask)=nan;
netcdf.close(ncid);


t_an = shiftdim(t_an,2);
s_an = shiftdim(s_an,2);


mask = single(~isnan(squeeze(t_an(1,:,:))));
mask(mask==0)=nan;

X=X1.*mask;
Y=Y1.*mask;
X=X(:);Y=Y(:);
t_an = reshape(t_an,nz,nlon*nlat);
s_an = reshape(s_an,nz,nlon*nlat);

for ii=1:length(latitude2)
    dist = abs(X-longitude2(ii))+abs(Y-latitude2(ii));
    [~,aa]= min(dist);
    T_lev(:,ai(ii)) = squeeze(t_an(:,aa));
    S_lev(:,ai(ii)) = squeeze(s_an(:,aa));
end  %ii
end %jj

if strcmp(modo,'pad')
%EXTRAP TO 6500
Dbase = D_lev(:,ones(1,length(latitude)));Dbase(isnan(T_lev))=nan;
    [dmax,dvec] = max(Dbase);
  [ELEV,LONG,LAT]=m_tbase([LONG_MIN LONG_MAX LAT_MIN LAT_MAX]);
    F = scatteredInterpolant(LONG(:),LAT(:),ELEV(:));
D_EXT = [5600:50:6500]';
D_lev = cat(1,D_lev,D_EXT);
  nz2 = length(D_lev);
bath  = ones(nz2,1)*F(longitude(:),latitude(:))'*-1;

%Uses the maximum depth for extrapolation
for kk = 1:length(latitude)
  T_lev(nz+1:nz2,kk) = ones(nz2-nz,1)*T_lev(dvec(kk),kk);
  S_lev(nz+1:nz2,kk) = ones(nz2-nz,1)*S_lev(dvec(kk),kk);
   
end
%Uses the last depth for extrapolation
%T_lev(nz+1:nz2,:) = ones(nz2-nz,1)*T_lev(nz,:);
%S_lev(nz+1:nz2,:) = ones(nz2-nz,1)*S_lev(nz,:);

  mask = D_lev(:,ones(1,length(latitude)))>bath;
  T_lev(mask)=nan;
  S_lev(mask)=nan;
end   %pad


return
