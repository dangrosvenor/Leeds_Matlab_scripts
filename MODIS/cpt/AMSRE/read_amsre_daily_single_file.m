%single read of AMSRE daily file

%filedir = '/home/disk/eos10/robwood/AMSR/bmaps_v05/y2007/m01/';
filedir = '/home/disk/eos5/d.grosvenor/AMSRE/y2002/m06/';
%file_name2 = 'amsre_200701v5.gz';
file_name2 = 'amsre_200701v7'; %these are monthly values - should be 6 maps - however they've averaged ascending and descending
file_name2 = 'amsre_20070131v5'; %thess are the daily maps (should be 7 maps - the extra one is time UTC,
%which is the first one)
%file_name2 = 'amsre_20070131v5_d3d'; %d3d are the 3 day averages
%file_name2 = 'amsre_20111001v7'; %weekly files - same as monthly (i.e. have also averaged day and night)
file_name2 = 'amsre_20020612v7';

filedir = '/home/disk/eos5/d.grosvenor/AMSRE/y2008/m10/'; file_name2 = 'amsre_20081026v7';

file_name =[filedir file_name2];


[time,sst,wspdLF,wspdMF,vapor,amsre_lwp,rain]=read_amsr_day_v7_Dan(file_name);
%   [mingmt,sst,wspdLF,wspdMF,vapor,cloud,rain,wspdAW,wdir]
%   mingmt is gmt time in hours
%   sst  is surface water temperature at depth of about 1 mm in deg C
%   wspdLF is 10 meter surface wind in m/s made using 10.7 GHz channel and above, low frequency channels
%   wspdMF is 10 meter surface wind in m/s made using 18.7 GHz channel and above, medium frequency channels
%   vapor is atmospheric water vapor in millimeters
%   cloud is liquid cloud water in millimeters
%   rain  is rain rate in millimeters/hour
%   wspdAW is 10 meter surface wind for all weather conditions made using 3 algorithms
%   wdir is wind direction oceanographic convention, blowing North = 0 in degrees
%
%  The center of the first cell of the 1440 column and 720 row map is at 0.125 E longitude and -89.875 latitude.
%  The center of the second cell is 0.375 E longitude, -89.875 latitude.
% 		XLAT=0.25*ILAT-90.125
%		XLON=0.25*ILON-0.125

%convert LWP from mm to g/m2. Mass_water/area = depth*rhoW
%amsre_lwp = amsre_lwp * 1e-3 * 1e3 * 1e3; %i.e. convert mm to m, multiply by rhoW and convert kg into g
amsre_lwp = amsre_lwp * 1e3;

