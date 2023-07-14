%read Rob's LWP netCDF file


ncdir='/home/disk/eos10/robwood/VOCALS/REMSS/';
ncfile='vocals_microwave.nc';

nc = netcdf([ncdir ncfile]);

lat_lwp = nc{'Latitude'}(:);
lon_lwp = nc{'Longitude'}(:);
time_lwp = nc{'Time'}(:);

LWP_lwp = nc{'LWP'}(:);
LWPanom_lwp = nc{'LWPanom'}(:);


% netcdf vocals_microwave {
% dimensions:
%         Dimension_1 = 200 ;
%         Dimension_2 = 160 ;
%         Dimension_3 = 61 ;
%         Dimension_4 = 4 ;
% variables:
%         float LWP(Dimension_4, Dimension_3, Dimension_2, Dimension_1) ;
%         float LWPanom(Dimension_4, Dimension_3, Dimension_2, Dimension_1) ;
%         float Latitude(Dimension_2) ;
%         float Longitude(Dimension_1) ;
%         float Time,days(Dimension_3) ;
%         float Localtime(Dimension_4) ;
% }
