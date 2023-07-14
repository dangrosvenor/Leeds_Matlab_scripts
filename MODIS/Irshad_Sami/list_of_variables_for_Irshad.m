% %savevarname = '/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/11km_matches_Dan.mat';
% 
%      16 means 1.6 micron and 37 3.7 micron
%      MODIS_swath_filename = the MODIS L2 file used for each datapoint
%      (cell data, so e.g. MODIS_swath_filename{1})
%      Nd_CDP = Nd from the CDP
%      Droplet_Number_Concentration  =  the mean droplet concentration
%           calculated using the Nd of the individual pixels in the 11x11km
%           square
%      Standard_Deviation_Droplet_Number_Concentration  =  as above but the
%           std. dev


%   Absolute_Combined_Error_Mean_Droplet_Number_Concentration      70x1               560  double              
%   Cloud_Effective_Radius_16_Liquid_Log_Mean                      70x1               560  double              
%   Cloud_Effective_Radius_16_Liquid_Maximum                       70x1               560  double              
%   Cloud_Effective_Radius_16_Liquid_Mean                          70x1               560  double              
%   Cloud_Effective_Radius_16_Liquid_Minimum                       70x1               560  double              
%   Cloud_Effective_Radius_16_Liquid_Standard_Deviation            70x1               560  double              
%   Cloud_Effective_Radius_37_Liquid_Log_Mean                      70x1               560  double              
%   Cloud_Effective_Radius_37_Liquid_Maximum                       70x1               560  double              
%   Cloud_Effective_Radius_37_Liquid_Mean                          70x1               560  double              
%   Cloud_Effective_Radius_37_Liquid_Minimum                       70x1               560  double              
%   Cloud_Effective_Radius_37_Liquid_Standard_Deviation            70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Log_Mean                         70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Log_Mean_Uncertainty             70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Maximum                          70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Mean                             70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Mean_Uncertainty                 70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Minimum                          70x1               560  double              
%   Cloud_Effective_Radius_Liquid_Standard_Deviation               70x1               560  double              
%   Cloud_Fraction_Combined                                        70x1               560  double              
%   Cloud_Fraction_Ice                                             70x1               560  double              
%   Cloud_Fraction_Liquid                                          70x1               560  double              
%   Cloud_Fraction_Liquid2                                         70x1               560  double              
%   Cloud_Fraction_Liquid_Pixel_Counts                             70x1               560  double              
%   Cloud_Fraction_Liquid_Pixel_Counts2                            70x1               560  double              
%   Cloud_Fraction_NoOpt                                           70x1               560  double              
%   Cloud_Fraction_Undetermined                                    70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Log_Mean                        70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Log_Mean_Uncertainty            70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Maximum                         70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Mean                            70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Mean_Uncertainty                70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Minimum                         70x1               560  double              
%   Cloud_Optical_Thickness_Liquid_Standard_Deviation              70x1               560  double              
%   Cloud_Top_Temperature_Day_CF_Mean                              70x1               560  double              
%   Cloud_Top_Temperature_Day_CF_Standard_Deviation                70x1               560  double              
%   Cloud_Top_Temperature_Day_Maximum                              70x1               560  double              
%   Cloud_Top_Temperature_Day_Mean                                 70x1               560  double              
%   Cloud_Top_Temperature_Day_Minimum                              70x1               560  double              
%   Cloud_Top_Temperature_Day_Standard_Deviation                   70x1               560  double              
%   Cloud_Top_Temperature_Day_ice_liq_Mean                         70x1               560  double              
%   Cloud_Top_Temperature_Day_ice_liq_Standard_Deviation           70x1               560  double              
%   Cloud_Top_Temperature_Ice_Day_Maximum                          70x1               560  double              
%   Cloud_Top_Temperature_Ice_Day_Mean                             70x1               560  double              
%   Cloud_Top_Temperature_Ice_Day_Minimum                          70x1               560  double              
%   Cloud_Water_Path_Liquid                                        70x1               560  double              
%   Cloud_Water_Path_Liquid_Log_Mean                               70x1               560  double              
%   Cloud_Water_Path_Liquid_Standard_Deviation                     70x1               560  double              
%   Date_Time_Swath                                                70x1               560  double              
%   Droplet_Number_Concentration                                   70x1               560  double              
%   Droplet_Number_Concentration_16                                70x1               560  double              
%   Droplet_Number_Concentration_37                                70x1               560  double              
%   LWC_CDP                                                        70x1               560  double              
%   MLAT                                                            1x1                 8  double              
%   MLON                                                            1x1                 8  double              
%   MODIS_swath_filename                                           70x1             15462  cell                
%   Nd_CDP                                                         70x1               560  double              
%   Percent_Error_Mean_Droplet_Number_Concentration                70x1               560  double              
%   Pixel_fraction_CF                                              70x1               560  double              
%   Pixel_fraction_Nd                                              70x1               560  double              
%   Sensor_Azimuth_Mean                                            70x1               560  double              
%   Sensor_Zenith_Mean                                             70x1               560  double              
%   Solar_Azimuth_Mean                                             70x1               560  double              
%   Solar_Zenith_Mean                                              70x1               560  double              
%   Standard_Deviation_Droplet_Number_Concentration                70x1               560  double              
%   Standard_Deviation_Droplet_Number_Concentration_16             70x1               560  double              
%   Standard_Deviation_Droplet_Number_Concentration_37             70x1               560  double              
%   Total_pixels                                                   70x1               560  double              
%   stdNd_CDP                                                      70x1               560  double              
%   weather_index                                                  70x1               560  double  