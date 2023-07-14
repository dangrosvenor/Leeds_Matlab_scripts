%define the variables that will be made, so that they can be cleared easily at the
%start
clear modis_var; istring=1;


        
% modis_var{istring}='Cloud_Water_Path_Liquid_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1;        
% 

 modis_var{istring}='Cloud_Fraction_Liquid'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
 modis_var{istring}='Cloud_Fraction_Liquid_Pixel_Counts'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Fraction_Ice'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Fraction_Undetermined'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
 modis_var{istring}='Cloud_Fraction_Combined'; field_string{istring}='.timeseries3'; istring=istring+1;%the normal CF
% modis_var{istring}='Cloud_Fraction_NoOpt'; field_string{istring}='.timeseries3'; istring=istring+1;%fraction of pixels that were designated cloud, but had no optical retrievals
% 
% 
 modis_var{istring}='Cloud_Effective_Radius_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Effective_Radius_Liquid_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Effective_Radius_Liquid_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Effective_Radius_Liquid_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Effective_Radius_Liquid_Mean_Uncertainty'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Effective_Radius_Liquid_Log_Mean_Uncertainty'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% 

 modis_var{istring}='Cloud_Optical_Thickness_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=337)
 modis_var{istring}='Cloud_Optical_Thickness_Liquid_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Optical_Thickness_Liquid_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Optical_Thickness_Liquid_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Optical_Thickness_Liquid_Mean_Uncertainty'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% modis_var{istring}='Cloud_Optical_Thickness_Liquid_Log_Mean_Uncertainty'; field_string{istring}='.timeseries3'; istring=istring+1;%(nvar=397)
% 
% 
 modis_var{istring}='Cloud_Top_Temperature_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
 modis_var{istring}='Cloud_Top_Temperature_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1;
% 
%         
%         
% %clear modis_var; istring=1;
 modis_var{istring}='Solar_Zenith_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
 modis_var{istring}='Solar_Azimuth_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
% modis_var{istring}='Sensor_Zenith_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
 modis_var{istring}='Sensor_Azimuth_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;

%        modis_var{istring}='Solar_Zenith_Maximum2'; field_string{istring}='.timeseries3'; istring=istring+1;
% 
  modis_var{istring}='Cloud_Top_Temperature_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1;
%  modis_var{istring}='Cloud_Top_Temperature_Day_ice_liq_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
%  modis_var{istring}='Cloud_Top_Temperature_Day_ice_liq_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1;
% 
% 
% 
% %modis_var{istring}='MLAT'; field_string{istring}=''; istring=istring+1; 
% %modis_var{istring}='MLON'; field_string{istring}=''; istring=istring+1; 
% 
% %modis_var{istring}='Cloud_Water_Path_Liquid'; field_string{istring}='.timeseries3'; istring=istring+1;
% modis_var{istring}='Cloud_Water_Path_Liquid_Log_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
% 
%  modis_var{istring}='Cloud_Optical_Thickness_Liquid_Log_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
%    


% modis_var{istring}='CTH'; field_string{istring}='.timeseries3'; istring=istring+1;
% modis_var{istring}='N_time3'; field_string{istring}='.timeseries3'; istring=istring+1; 
% modis_var{istring}='N_time3_16'; field_string{istring}='.timeseries3'; istring=istring+1;  
% modis_var{istring}='N_time3_37'; field_string{istring}='.timeseries3'; istring=istring+1;  
% modis_var{istring}='W_time3'; field_string{istring}='.timeseries3'; istring=istring+1;   
 
 