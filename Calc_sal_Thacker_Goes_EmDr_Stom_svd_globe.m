%function [y2,y3,TT,PP]=Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(T,P,lat,lon,time,pad,Pout,method)
%TS lookup for the Atlantic basin
%INPUT:  T: Temperature (degC) 
%        P: Depth       (m, positive)
%        Lat:   (Degree north)
%        Lon:   (Degree East [-100:20])
%        time:   E.g. 200410 (yyymm)
%        pad:    0 [Default] or 1
%        Pout: Depth output (m, positive)
%        method: 'Thacker'[Default] / 'Goes' / 'annual' / 'stommel' / 'svd'
%OUTPUT: y2: Salinity from the TS method chosen
%        y3: smoothed version of y2 
%        TT: Temperature output using Pout as vertical axis if defined
%        PP: Depth defined for the profiles

%EXAMPLE

%load data_example.mat
%aa = randi(length(time),10000,1); %assign 10000 random samples 
%%P = P'*ones(1,10000);  %not needed, but works either way
%S = S(:,aa);
%T = T(:,aa);
%latitude  = latitude(aa);
%longitude = longitude(aa);
%time = time(aa);
%[S1,S2]=Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(T,P,latitude,longitude,time);

%UPDATES:
%June 1 - added Emery&Dewar (annual) to the reconstruction
%June 19 - corrected interpolated depth
%          declare pred5 variable
%          delete values below 5 psu, if happens
%June 26 - fixed varargin
%          exclude pred4 and pred5        
%Aug 17  - include nans where P is nan    
%        - change padding to 800m instead of Pmax  
%        - replace outliers with mean instead of nans
%Oct 4   - added acha2 to l.255 interpolate back
%        - added t0 for mean temp
%        -INCLUDE STOMMEL METHOD
%Oct 10  - added p to insnan in line 93
%Oct 25  - Padding fixed to allow multiple depths
%Nov27   - Relaxed aza from S-Smean>1 to S-Smean>2
%June6,2018 - fixed tnan in stommel   (CHECK!!!!)
%FEB 2019  - added nnz(nonan)>1 for padding interpolation
%           - added Dref as a substitute for pad . 
%             Just choose pad > 1 and this will be the reference depth
%           - in smoothing changed isempty(facha2) to facha2<2
%           line 344: acha2 = (~isnan(s2+P));
%MARCH 6 -MAJOR: corrected padding grid. Works great now.
%MARCH 26 - MAJOR: FIXED padding Dmax definition
%          - ddnan changed to & instead of |
%NOV 2020 - NPI,TPI changed to 15N instead of 10N
% July 2022 - Adjusted interpolation for one dimensional data (bin1d-l.142)
% July 2022 - Adjusted interpolation for one dimensional data (interp1q-l.376)
% OCT  2022 - Added FC - 93,163,244,259
% DEC  2022 - Correct adjusted interpolation to avoid NaNs (elseif l.389)
% DEC  2022 - Include LAT>26.5  (26.8 before)for FC (l. 97).
% APR  2023 - FC from Ryan's data. make*_ryan.m . %_FC files updated
% DEC  2023 - Changed FC to _ryanonly  make*_ryanonly.m  %_FC files updated - no cora/argo
% Mar  2024 - Included aux(acha2) in l. 403. Was aux;
% Feb  2025 - Added MEDITERRANEAN
% Jul  2025 - Corrected call for Pacific (l.278 ) and added Zcorr for interpolation to maximum depth (l.155-156)
function [y2,y3,TT,PP]=Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(varargin) %TT,PP,lat,lon,timet,pad,Pout,method)
InpVars = cell(1,8);
InpVars(1:nargin) = varargin;
[TT,PP,lat,lon,timet,pad,Pout,method] = InpVars{:};
IntP = true;
if isempty(pad), pad    = false;    end
if isempty(Pout), IntP = false;     end
if isempty(method), method = 'thacker';else method = lower(method);end
if ~any(strcmp(method,{'goes','thacker','annual','stommel','svd'}))
    error('MethodNotFound');
end

Dref = 800; if pad > 1, Dref=pad; pad = 1;end             %ADDED PAD DEPTH

%DEFINE LON:0-360
lon2 = lon; lon2(lon2<0)=lon2(lon2<0)+360;

disp(method)

Z = [7.5:10:300 320:20:1000 1050:50:2000]';
disp('Contruct year, month arrays')

%Contruct year, month arrays
if length(timet)==1; timet=timet*ones(size(lat));end
timet = timet(:);
timet(isnan(timet))=0;

%DEFINE REGIONAL INDICES
disp('Load Fitted parameters')
sai = find(lat >=-45 & lat <-15 & lon >-70 & lon < 20);
    basin(sai,1:2) = repmat('SA',length(sai),1);

tai = find(lat >=-15 & lat <10 & lon>-80 & lon < 20);
    basin(tai,1:2) = repmat('TA',length(tai),1);

nai = find(lat >=10 & lat <45 & lon >-98 & lon < 0);
    basin(nai,1:2)  = repmat('NA',length(nai),1);

fci =  find(lat >=26.5 & lat <27.2 & lon >-80.2 & lon < -78.9);  %WAS 26.8 
    basin(fci,1:2)  = repmat('FC',length(fci),1);nai(fci)=[];

mdi = find(lat >=30 & lat<=45 & lon > -6 & lon <= 30);           %ADDED MD Jan2025
    basin(mdi,1:2)  = repmat('MD',length(mdi),1);

noi = find(lat >=45 & lat <65 & lon >-100 & lon < 30);
    basin(noi,1:2)  = repmat('NO',length(noi),1);

soi = find(lat >=-70 & lat <-45 & lon >-70 & lon <= 20);
    basin(soi,1:2)  = repmat('SO',length(soi),1);

tii = find(lat >=-45 & lat < 30 & lon >20 & lon < 120);
    basin(tii,1:2)  = repmat('TI',length(tii),1);

soii = find(lat >= -75 & lat < -45 & lon > 20 & lon < 120);
    basin(soii,1:3)  = repmat('SOI',length(soii),1);

%NEED TO DEFINE LON2: 0-360
spi = find(lat >=-45 & lat <-15 & lon2 >120 & lon2 < 290);
    basin(spi,1:2) = repmat('SP',length(spi),1);

tpi = find(lat >=-15 & lat <15 & lon2>120 & lon2 < 285);
    basin(tpi,1:2) = repmat('TP',length(tpi),1);

npi = find(lat >=15 & lat <45 & lon2 >120 & lon2 < 260);
    basin(npi,1:2)  = repmat('NP',length(npi),1);

nopi = find(lat >=45 & lat <65 & lon2 >120 & lon2 < 260);
    basin(nopi,1:3)  = repmat('NOP',length(nopi),1);

sopi = find(lat >=-75 & lat <-45 & lon2 >120 & lon2 <= 290);
    basin(sopi,1:3)  = repmat('SOP',length(sopi),1);


    [~,nla1,nmo1]=size(TT);
    
    if length(PP)==numel(PP)
   PP = repmat(PP(:),[1,nla1,nmo1]);
    end
  %ADDED TO AVOID NAN's at the surface
  %PP(PP<Z(1)) = Z(1);
    
PP2 = PP;
TT2 = TT;
T = nan*ones(85,nla1);
for ii=1:nla1
    for jj = 1:nmo1
        t = TT(:,ii,jj);
        p = PP(:,ii,jj);
        acha = ~isnan(t+p);
              if ~any(acha);continue;end
        if nnz(acha)>1
            Zcorr = Z;
            Zcorr(find(Z>max(p),1)) = max(p); %Correct for maximum depth %Jul 2025
            T(:,ii,jj) = interp1(p(acha),t(acha),Zcorr);
        elseif nnz(acha)==1 || max(p(acha))<Z(1)
            [~,imin]=min(abs(Z-p(acha)));                    %DEC 2023
            aux = Z*nan;aux(imin)=t(acha);
            T(:,ii,jj)= aux;
%          aux = bin1d(p(acha),t(acha),Z);  %July 2022
%          T(:,ii,jj)= aux.mean;
        end
    end
end
TT=T;clear T;

    [nz,nla1,nmo1]=size(TT);
%if length(PP)==numel(PP) %ORIGINAL  %prod(size(PP))
if size(PP,1)==numel(PP) %July 2022
PP = Z;
else
PP = Z(:,ones(1,nla1),ones(1,nmo1));
end 

 nt = nla1 * nmo1;

   y2 = nan*ones(nz,nla1,nmo1);
 stnai = find([~isempty(tai) ~isempty(sai) ~isempty(nai) ~isempty(soi) ...
    ~isempty(noi) ~isempty(tii) ~isempty(soii) ~isempty(nopi) ...
    ~isempty(npi) ~isempty(tpi) ~isempty(spi) ~isempty(sopi) ...
    ~isempty(fci) ~isempty(mdi)] );

for bb = stnai 

%LOAD FITTED PARAMETERS   
if strcmp(method,'svd')
   preffix = 'fit_sigma_cora_argo_extended_svd_';
   suffix  = '.mat';
else
   preffix = 'fit_sigma_cora_argo_extended_fitlm_';
   suffix  = '_thackergoes_noyear.mat';
end

switch bb
    case 1
     disp('Load TA')   
filename = sprintf('%sTA%s',preffix,suffix);
load(filename)
ai = tai;
    case 2
     disp('Load SA')   
filename = sprintf('%sSA%s',preffix,suffix);
load(filename)
ai = sai;
    case 3
     disp('Load NA')   
filename = sprintf('%sNA%s',preffix,suffix);
load(filename)
ai = nai;
    case 4
     disp('Load SO')   
filename = sprintf('%sSO%s',preffix,suffix);
load(filename)
ai = soi;
    case 5
     disp('Load NO')   
filename = sprintf('%sNO%s',preffix,suffix);
load(filename)
ai = noi;
   case 6
%INDIAN
     disp('Load TI')
filename = sprintf('%sTI%s',preffix,suffix);
load(filename)
ai = tii;
    case 7
     disp('Load SOI')
filename = sprintf('%sSOI%s',preffix,suffix);
load(filename)
ai = soii;
    case 8
%PACIFIC
     disp('Load NOP')
filename = sprintf('%sNOP%s',preffix,suffix);
load(filename)
ai = nopi;
    case 9
     disp('Load NP')
filename = sprintf('%sNP%s',preffix,suffix);
load(filename)
ai = npi;
    case 10
     disp('Load TP')
filename = sprintf('%sTP%s',preffix,suffix);
load(filename)
ai = tpi;
    case 11
     disp('Load SP')
filename = sprintf('%sSP%s',preffix,suffix);
load(filename)
ai = spi;
    case 12
     disp('Load SOP')
filename = sprintf('%sSOP%s',preffix,suffix);
load(filename)
ai = sopi;
   case 13
     disp('Load FC')
filename = sprintf('%sFC%s_ryanonly.mat',preffix,suffix(1:end-4));
load(filename)
ai = fci;
    case 14
        disp('Load MED')
filename = sprintf('%sMD2%s',preffix,suffix)
load(filename)
ai = mdi;

end
clear rmse_ii std_er_ii AICBIC_ii LAMBDA_ii RMS_ii

T = TT(:,ai);
P = PP(:,ai);

%longitude = lon(ai);
%DIFFER FOR PACIFIC
if (bb < 8 || b > 12)
   longitude = lon(ai);
else
   disp('Pacific')
   longitude = lon2(ai);
end

latitude  = lat(ai);
     time = timet(ai,:);
    alnan = time==0;

  nt = length(time);
time = num2str(time);

 year  =  nan*ones(nt,1);
 month =  nan*ones(nt,1);
% 
 year(~alnan)   = str2num(time(~alnan,1:4)); 
 year(~alnan)   = min(max(year(~alnan),1993),2015);        %CORRECT YEAR RANGE
 month(~alnan)  = str2num(time(~alnan,5:6));

[nz,nla,nmo]=size(T);

 yy = repmat(year',nz,1);            size(yy);    
 mo = repmat(month',nz,1);           size(mo);    
 pp = P;                             size(pp);    
 latitude  = latitude(:);        
 longitude = longitude(:); 
 T    = reshape(T,nz,nla*nmo);        size(T);
 
 %%%NEW%%%%%%%%%%%%%%%%%%%%%%%%%
 pred3 = nan*ones(size(T));

 aan = find(~alnan);
 for aa = aan'%1:length(s)
    dflat = abs(latitude(aa)-Y(:))+abs(longitude(aa)-X(:));
    [~,ll] = min(dflat);
  xmean = squeeze(xmean_ii{ll});
   if dflat(ll) > 2 || isempty(xmean),continue,end
  ymean = ymean_ii{ll};    %Annual climatology
  t0 = T(:,aa);
  t2 = t0 - xmean;

 if strcmp(method,'svd')      
    beta  = CS_ii(ll,:);    %Coefficients Salinity
    U     = U_ii{ll,:};    %Coefficients Salinity
    V     = V_ii{ll,:};    %Coefficients Salinity
   t2(isnan(t2))=0;
   yn = 0;
  for ee=1:15
     Cs(ee) = beta(ee); %regress(a(ee,:)',b(ee,:)');
     %Ct(ee) = regress(b(ee,:)',a(ee,:)');
     %for jj=1:np
     yn=yn+Cs(ee)*U(:,ee)*V(:,ee)'*t2(:);
     %yt(:,jj)=yt(:,jj)+Ct(ee)*V(:,ee)*U(:,ee)'*PS(:,jj);
  end

yn = yn + ymean;


else

  beta1  = beta_ii{ll,1};   %Coefficients Thacker
  beta2  = beta_ii{ll,2};   %Coefficients Goes
 
      lat_ii = latitude(aa*ones(1,85));
      lon_ii = longitude(aa*ones(1,85));
      x0 = X(ll);y0 = Y(ll);
  x12 = [ones(85,1) t2(:) t2(:).^2 cos(2*pi*mo(:,aa)/12),sin(2*pi*mo(:,aa)/12), ...
    cos(4*pi*mo(:,aa)/12), sin(4*pi*mo(:,aa)/12) lat_ii(:)-y0,lon_ii(:)-x0];
     
     betam = 0*beta1;
 if strcmp(method,'goes')
     betam = beta2;
 elseif strcmp(method,'thacker')
     betam = beta1;
 end
    if strcmp(method,'stommel')
        %truncate temp
        nnan = ~isnan(t0);
        tnan = ~isnan(ymean);
        t3 = max(t0(nnan),min(xmean));t3 = min(t3,max(xmean));
        try
        pred3(1:nnz(nnan),aa) = interp1(xmean(tnan),ymean(tnan),t3);
        catch me
            [xmean,xunic] = unique(xmean);
            ymean=ymean(xunic);
            tnan = ~isnan(ymean)&~isnan(xmean);%&[0; diff(xmean)]<0;
         pred3(1:nnz(nnan),aa) = interp1(xmean(tnan),ymean(tnan),t3);
        end
    else
        
        for kk = 1:nz %length(P)
            aw = find(pp(:,aa) == Z(kk));
          pred3(aw,aa) = x12(aw,:)*betam(kk,:)' + ymean(kk);

         
        end
    end
          yn = pred3(:,aa);

end %if method
         aza = find(yn > 42 | yn < 5);% | (abs(yn-ymean)>2 & pp(:,aa)>600));
         if any(aza),disp('change to clim');end
         yn(aza) = nan;
         %yn(aza) = ymean(aza);
         pred3(:,aa) = yn;
  
 end
  
% clear yy mo latitude longitude pp t

y2(:,ai)   = pred3;%reshape(y,nz,nla);
end
    
 nz   = size(PP2,1);
 PP   = PP2(1:nz,:);
 nlon = size(y2,2);

 %INTERPOLATE BACK TO ORIGINAL
Y=nan*ones(nz,nlon);%nla1*nmo1);%nla);
for ii=1:nlon%nla   
        y = y2(:,ii);
        p = PP(:,ii);
        acha = ~isnan(y);
        acha2 = ~isnan(p);                       %added improves interpolation
        if nnz(acha)>1
   Y(acha2,ii) = interp1(Z(acha),y(acha),p(acha2));
        elseif nnz(acha)==1
            [~,imin]=min(abs(Z(acha)-p(acha2)));       %DEC 2023
            aux = p*nan;aux(imin)=y(acha);
            Y(acha2,ii) = aux(acha2);% MAR 2024        Y(:,ii) = nan; %y(acha);
        end
end
y2 = Y; 
TT=TT2; %PP=PP2;
clear TT2 PP2 Y z

%SLAB approximation to surface values%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nlon = size(y2,2);
 for jj =1:nlon
     s2 = y2(:,jj);
     P = PP(:,jj);                     %new
     facha1 = find(~isnan(s2));
         if ~isempty(facha1)
          if facha1(1)>1
         s2(1:facha1(1)) =   s2(facha1(1));
          end
         end
        %INTERPOLATE NAN"S WITHIN
        % acha2 = (~isnan(s2));
         acha2 = (~isnan(s2+P));
        if (nnz(acha2)>2 && nnz(acha2)<nz)
           s2(~acha2) = interp1(P(acha2),s2(acha2),P(~acha2));
           s2(isnan(P))=nan;           %new Aug17,17
        end 
         y2(:,jj) = s2;

 end
pai = find(PP>2000);
y2(pai)=nan;
y2(isnan(TT))=nan;                  %ADD NAN's where T's are


%ADD PADDING AND BATHYMETRY
%nz2 = size(Pout,1);
time = num2str(timet);
month  = str2num(time(:,5:6));
nt = length(month);



%TESTE NEW
if pad & 1
    disp('Load Climatology...')
    [T_lev,S_lev,D_lev]=load_woa13_pad2(lat,lon,month);%,'pad'); %NOW WITH BATHYMETRY INCLUDED
    Pmax = max(PP(:));
        P=unique(round(PP(:)));P=P(~isnan(P));
    [~,Dmax] = min(abs(P-Dref)); %NOW USES DREF input from pad is available
    P    = P(1:Dmax);
    yi = nan*ones(Dmax,nt);Ti=yi;
    for jj=1:nt
        nonan = ~isnan(TT(:,jj));
        if nnz(nonan)>1                %ADDED FEB  2019
           yi(1:Dmax,jj)  = interp1(PP(nonan,jj),y2(nonan,jj),P);
           Ti(1:Dmax,jj)  = interp1(PP(nonan,jj),TT(nonan,jj),P);
        end
    end
    y2 = yi;TT=Ti;
    Dmin  = (D_lev>Dref);% Pmax);%     %THIS IS NOW 850m; change Pmax if necessary
    Z_I   = cat(1,P,D_lev(Dmin));
    PP    = Z_I*ones(1,nt);   
    T_lev = interp1(D_lev,T_lev,Z_I);
    S_lev = interp1(D_lev,S_lev,Z_I);
    y2(Dmax+1:length(Z_I),:) = nan;
    TT(Dmax+1:length(Z_I),:) = nan;
    ddnan = isnan(y2) | PP>Dref;%Pmax;%800;  %CHANGED TO 800 FOR LOWER PADDING
    y2(ddnan) = S_lev(ddnan);
    ddnan = isnan(TT) | PP>Dref;%Pmax;%800;  %CHANGED TO 800 FOR LOWER PADDING
    TT(ddnan) = T_lev(ddnan);
end
clear Ti yi


 if ~IntP, Pout=PP;end


%Smoothing Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz   = size(Pout,1);
nlon = size(y2,2);
y3   = nan*ones(nz,nlon);
%y3 = y2;
TT2  = y3*nan;
disp('INTERPOLATION/SMOOTHING...')

for jj =nlon:-1:1  %:nlon
     s3 = y2(:,jj);
     t3 = TT(:,jj);
     P  = PP(:,jj);                     %new
     facha2 = find(~isnan(s3) & s3>0);
     %if isempty(facha2),continue,end
     if length(facha2)<2,continue,end    %CHANGED FEB2019

     %INTERPOLATE TO DEFINE DEPTHS

    if IntP
      s3 = interp1(P(facha2),s3(facha2),Pout(:));
      t3 = interp1(P(facha2),t3(facha2),Pout(:));
      facha2 = find(~isnan(s3));
      if isempty(facha2),continue,end
    end


  %  tfilt2 = medfilt1(s3,7,'omitnan');
    tfilt2 = s3;%ORIGINAL
       nz2 = size(tfilt2,1);

    tfilt2(tfilt2==0)=nan;

 if nargout > 1

    %This is to filter spikes (NEW()
    eeo = tfilt2;%s3;%eeo = ee;   


    dif1 = eeo(2:end)-eeo(1:end-1);
      oi = dif1;%./dz;

     oi(isnan(oi))=0;

      epsl = 0.1;  %Threshold for gradient
      nei = (abs(s3(3:end)-s3(2:end-1) - (s3(2:end-1)-s3(1:end-2)))>epsl);       %NEW METHOD

      oi(nei) = interp1(find(~nei),oi(~nei),find(nei));
      oi(isnan(oi))=0;

%       %inverse cummulative sum
         e0 =  s3(facha2(end));
        oii = [-1*oi;e0]; %if use dif1
       eeoi = flipud(cumsum(flipud(oii)));

      %UPDATE SALINITY 
      tfilt2(1:facha2(end)) = eeoi(1:facha2(end));
    %  tfilt2 = medfilt1(tfilt2,7);%,'omitnan');  %ORIGINAL
 end
       y3(1:nz2,jj)  = tfilt2;
       TT2(1:nz2,jj)  = t3(1:nz2);
        y2(1:nz2,jj)  = s3;               %UPDATE Y2 to Interpolated 
end
[nz,nla] = size(TT2);
TT = TT2;
PP = Pout;
y3 = y3(1:nz,1:nla);y3(y3==0)=nan;
y2 = y2(1:nz,1:nla);y2(y2==0)=nan;

%y3(pai,:)=nan;
%y3(isnan(TT))=nan;                      %ADD NAN's where T is NAN
%y2(isnan(TT))=nan;                  %ADDed here NAN's where T's are
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear beta_ii xmean_ii ymean_ii

return
