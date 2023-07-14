%MODIS L2 variable list

clear modis_var

istring=1;
modis_var{istring}='scantime';         istring=istring+1;
modis_var{istring}='scantime_matlab';  istring=istring+1;
modis_var{istring}='solar_zenith';     istring=istring+1;
modis_var{istring}='sensor_zenith';    istring=istring+1;
modis_var{istring}='tau_5km';    istring=istring+1;
modis_var{istring}='re_5km';    istring=istring+1;
modis_var{istring}='cf';    istring=istring+1;
modis_var{istring}='t_top';    istring=istring+1; %Cloud Top Temperature
modis_var{istring}='p_top';    istring=istring+1;
%modis_var{istring}='tau_un_5km';    istring=istring+1; %uncertainties
%modis_var{istring}='re_un_5km';    istring=istring+1;