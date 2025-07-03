#!\bin\csh

#This code is to install the TS lookup package from Goes et al., 2018.
#Run 'csh install.sh' to install this package
#Then run example.m to test the package

tar -xvzf TS_PACKAGE_ftp_Feb25.tar.gz

ftp https://www.aoml.noaa.gov/ftp/pub/phod/mgoes/TS/MATFILES/matfile.tar.gz
ftp https://www.aoml.noaa.gov/ftp/pub/phod/mgoes/TS/MATFILES/WOA18.tar.gz

cd TS_PACKAGE_THACKER_noyear_globe_ftp/

tar -xvzf matfile.tar.gz
tar -xzvf WOA18.tar.gz


 
