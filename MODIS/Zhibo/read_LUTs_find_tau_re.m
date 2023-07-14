%Have theses for SZA=20 (from IDL .sav file) and SZA=79 (NetCDF file)

%% IDL .sav file - read in using Matlab function - think is for nadir and
%% rza=0
sav_data=restore_idl('/home/disk/eos1/d.grosvenor/idl/LUTs/MODIS_LUT.sav');

%Have to double() all the data, though

R086=double(sav_data.REF_086);
R21=double(sav_data.REF_213);
R37=double(sav_data.REF_370);

re=double(sav_data.REG);
tau=double(sav_data.TAUG);
[re2,tau2]=meshgrid(re,tau);

tau_ret_sav = griddata(R086,R21,tau2,0.43,0.23)
re_ret_sav = griddata(R086,R21,re2,0.43,0.23)
%Answer agrees well with Fig. 2 of Zhang (2016)

% sav_data = 
% 
%            NTAUG: 35
%             TAUG: [35x1 single]
%             NREG: 18
%              REG: [18x1 single]
%             NMUG: 28
%              MUG: [28x1 single]
%     EMS_CLOUD_11: [28x35x18 single]
%      EMS_SURF_11: [28x35x18 single]
%     EMS_CLOUD_37: [28x35x18 single]
%      EMS_SURF_37: [28x35x18 single]
%          REF_086: [35x18 single]
%          REF_164: [35x18 single]
%          REF_213: [35x18 single]
%          REF_370: [35x18 single]

%%
%SZA=79 file (netCDF)

ncfile='/home/disk/eos1/d.grosvenor/idl/LUTs/MODIS_LUT_verylowsun2_netCDF3.nc'
nc=netcdf(ncfile);


bands=nc{'Band Center'}(:);   % =   0.8650 1.6400 2.1300 3.7500
rza=nc{'Relative Azimuth'}(:); % = 30 150
vza=nc{'Viewing Zenith Angle'}(:); % = 0 10 20 30 40 50

R=nc{'Reflectance'}(:); %size(R)=     [4    18    34     2     6]  (band, re, tau, phi, vza)

R086 = squeeze(R(1,:,:,1,1));
R21 = squeeze(R(3,:,:,1,1));
R37 = squeeze(R(4,:,:,1,1));

re=nc{'Cloud Droplet Effective Radius'}(:); %1D vector of length 18
tau=nc{'Cloud Optical Thickness'}(:); %length 34

[re2,tau2]=meshgrid(re,tau); %mesh these into 2d array, but need to transpose to line up with reflectance arrays
re2=re2'; tau2=tau2';


%Can now interpolate to find tau and re for a given reflectance pair, e.g. 
tau_ret = griddata(R086,R37,tau2,0.4,0.08)
re_ret = griddata(R086,R37,re2,0.4,0.08)

% ncdump -h /home/disk/eos1/d.grosvenor/idl/LUTs/MODIS_LUT_verylowsun2_netCDF3.nc
% netcdf MODIS_LUT_verylowsun2_netCDF3 {
% dimensions:
% 	vza = 6 ;
% 	sza = 1 ;
% 	phi = 2 ;
% 	tau - Anchored to 0.650 band = 34 ;
% 	re - micrometers = 18 ;
% 	band - micron = 4 ;
% 	tau = 34 ;
% 	re = 18 ;
% 	band = 4 ;
% variables:
% 	double Viewing Zenith Angle(vza) ;
% 	double Solar Zenith Angle(sza) ;
% 	double Relative Azimuth(phi) ;
% 	double Cloud Optical Thickness(tau - Anchored to 0.650 band) ;
% 	double Cloud Droplet Effective Radius(re - micrometers) ;
% 	double Band Center(band - micron) ;
% 	double Reflectance(band, re, tau, phi, vza) ;
% }
