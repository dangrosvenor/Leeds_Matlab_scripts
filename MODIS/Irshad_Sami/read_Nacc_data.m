%read the aerosol data from Irshad
% see >>> read_weather_data_Irshad <<< for the weather data script and
% >>> Puijo_station_cloud_event_info <<< for the CDP data. 

load('/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/DMPS_80_100_150nm.mat');
% Name                       Size              Bytes  Class     Attributes
% 
%   DMPS_80_100_150nm      43184x10            3454720  double  


date_DMPS = DMPS_80_100_150nm(:,1); %is based on mid-points of the hour (10:30, 11:30, etc)
Nacc_DMPS = DMPS_80_100_150nm(:,9); %accumulation mode aerosol concentration (cm3).


Nacc = NaN*ones(size(Date_Time_Swath.timeseries3));

dtol = 0.5/24; %match within one hour - needs to be in days. Will be +/- this value

for i=1:length(Date_Time_Swath.timeseries3)
   [minval,imin] = min( abs(date_DMPS -  Date_Time_Swath.timeseries3(i)) );
   if minval<dtol
       Nacc(i) = Nacc_DMPS(imin);
   end
end

NaccDMPS.timeseries3 = Nacc; %for saving purposes