%Irshad weather data

%weather_data_type = 'new with all Tout data';

%switch weather_data_type
%    case 'new with all Tout data'
%this contains data only from 19th June, 2008 to 10th Feb, 2013        
        
        
%        Before 7.11.2007 there are 32 columns:
%   1 = Matlab date
%   2 = Year
%   3 = Month
%   4 = Day
%   5 = Hour
%   6 = Minute
%   7 = Second
%   8 = Temperature out [C]
%   9 = Temperature in [C]
%   10 = RH out [%]
%   11 = RH in [%]
%   12 = Temperature of anemometer [C]
%   13 = 1-min average wind speed [m/s]
%   14 = 1-min minimum wind speed [m/s]
%   15 = 1-min maximum wind speed [m/s]
%   16 = 1-min deviation in wind speed [m/s]
%   17 = 1-min maximum wind speed in gusts [m/s]
%   18 = 1-min minimum wind speed in gusts [m/s]
%   19 = 1-min average wind direction [degrees]
%   20 = 1-min deviation in wind direction [degrees]
%   21 = VisibilityAlarm [code]
%   22 = 1-min average visibility [m]
%   23 = 10-min average visibility [m]
%   24 = InstantPresentWeather [code]
%   25 = 15minPresentWeather [code]
%   26 = 1hourPresentWeather [code]
%   27 = Precipitation intensity [mm/h]
%   28 = Cumulative Water Sum [mm]
%   29 = Cumulative Snow Sum [mm]
%   30 = Pressure [hPa]
%   31 = TT[C] (NOT TO BE USED)
%   32 = Error[code] (NOT TO BE USED)
% 
% Between 7.11.2007 and 17.6.2008 there are 35 columns:
%   1 = Matlab date
%   2 = Year
%   3 = Month
%   4 = Day
%   5 = Hour
%   6 = Minute
%   7 = Second
%   8 = Temperature out [C]
%   9 = Temperature in [C]
%   10 = RH out [%]
%   11 = RH in [%]
%   12 = Temperature of anemometer [C]
%   13 = 1-min average wind speed [m/s]
%   14 = 1-min minimum wind speed [m/s]
%   15 = 1-min maximum wind speed [m/s]
%   16 = 1-min deviation in wind speed [m/s]
%   17 = 1-min maximum wind speed in gusts [m/s]
%   18 = 1-min minimum wind speed in gusts [m/s]
%   19 = 1-min average wind direction [degrees]
%   20 = 1-min deviation in wind direction [degrees]
%   21 = VisibilityAlarm [code]
%   22 = 1-min average visibility [m]
%   23 = 10-min average visibility [m]
%   24 = InstantPresentWeather [code]
%   25 = 15minPresentWeather [code]
%   26 = 1hourPresentWeather [code]
%   27 = Precipitation intensity [mm/h]
%   28 = Cumulative Water Sum [mm]
%   29 = Cumulative Snow Sum [mm]
%   30 = Pressure [hPa]
%   31 = TT[C] (NOT TO BE USED)
%   32 = Error[code] (NOT TO BE USED)
%   33 = M1[g] (AN ICING PARAMETER)
%   34 = M2[g] (AN ICING PARAMETER)
%   35 = M3[C] (AN ICING PARAMETER)
% 
% 18.6.2008 was a bad day and the data are missing.
% 
% From 19.6.2008 on there are 46 columns:
%   1 = Matlab date
%   2 = Year
%   3 = Month
%   4 = Day
%   5 = Hour
%   6 = Minute
%   7 = Second
%   8 = (OBSOLETE, USE COLUMN 37) [Temperature out [C]]
%   9 = Temperature in [C]
%   10 = (OBSOLETE, USE COLUMN 36) [RH out [%]]
%   11 = RH in [%]
%   12 = Temperature of anemometer [C]
%   13 = 1-min average wind speed [m/s]
%   14 = 1-min minimum wind speed [m/s]
%   15 = 1-min maximum wind speed [m/s]
%   16 = 1-min deviation in wind speed [m/s]
%   17 = 1-min maximum wind speed in gusts [m/s]
%   18 = 1-min minimum wind speed in gusts [m/s]
%   19 = 1-min average wind direction [degrees]
%   20 = 1-min deviation in wind direction [degrees]
%   21 = VisibilityAlarm [code]
%   22 = 1-min average visibility [m]
%   23 = 10-min average visibility [m]
%   24 = InstantPresentWeather [code]
%   25 = 15minPresentWeather [code]
%   26 = 1hourPresentWeather [code]
%   27 = Precipitation intensity [mm/h]
%   28 = Cumulative Water Sum [mm]
%   29 = Cumulative Snow Sum [mm]
%   30 = Pressure [hPa]
%   31 = TT[C] (NOT TO BE USED)
%   32 = Error[code] (NOT TO BE USED)
%   33 = M1[g] (AN ICING PARAMETER)
%   34 = M2[g] (AN ICING PARAMETER)
%   35 = M3[C] (AN ICING PARAMETER)
%   36 = RH out [%]
%   37 = Temperature out [C]
%   38 = Tdf[C] (dew point or frost point temperature)
%   39 = Td[C] (dew point temperature)
%   40 = a[g/m3] (absolute humidity)
%   41 = x[g/kg] (mixing ratio)
%   42 = Tw[C] (wet temperature)
%   43 = H2O[ppmV] (ratio of humid air volume to dry air volume)
%   44 = pw[hPa] (water vapour pressure)
%   45 = pws[hPa] (saturated water vapour pressure)
%   46 = h[kJ/kg] (entalphy)
%   47 = dT[C] (difference between T and Tdf)

file_dir_Irshad = '/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/';
load([file_dir_Irshad 'weather_1hAVG_20080619_with_Tout.mat']);

        
weather_new = weather;

%the data from 19th June, 2008 - will add it to main array later
RHout_recent = weather_new(:,36);
Tout_recent = weather_new(:,37);
        
%    case 'old'

file_dir_Irshad = '/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/';
load([file_dir_Irshad 'weather_1hAVG_2.mat']);
weather_old = weather;

%loads a variable called "weather"
% � 1 = matlab date/time
% � 2 = Year
% � 3 = Month
% � 4 = Day
% � 5 = Hour
% � 6 = Minute
% � 7 = Second
% � 8 = Temperature out [C]
% � 9 = Temperature in [C]
% � 10 = RH out [%]
% � 11 = RH in [%]
% � 12 = Temperature of anemometer [C]
% � 13 = 1-min average wind speed [m/s]
% � 14 = 1-min minimum wind speed [m/s]
% � 15 = 1-min maximum wind speed [m/s]
% � 16 = 1-min deviation in wind speed [m/s]
% � 17 = 1-min maximum wind speed in gusts [m/s]
% � 18 = 1-min minimum wind speed in gusts [m/s]
% � 19 = 1-min average wind direction [degrees]
% � 20 = 1-min deviation in wind direction [degrees]
% � 21 = VisibilityAlarm [code]
% � 22 = 1-min average visibility [m]
% � 23 = 10-min average visibility [m]
% � 24 = InstantPresentWeather [code]
% � 25 = 15minPresentWeather [code]
% � 26 = 1hourPresentWeather [code]
% � 27 = Precipitation intensity [mm/h]


%end


%recent data goes from 19th June, 2008 to 10th Feb, 2013
%old data runs from to 31st Oct, 2012, 23:00
istart_old=find(weather_old(:,1)==datenum('19-Jun-2008 00:00'));
iend_new=find(weather_new(:,1)==datenum('31-Oct-2012 23:00'));

weather(istart_old:end,8) = Tout_recent(1:iend_new);
weather(istart_old:end,10) = RHout_recent(1:iend_new);

disp('Done read');
%contains an array 'weather' consisting of [n 26]
%Wind speed = column 13
%Wind dir = column 19
%Visibility = column 22






if ~exist('ioverride_read_weather_Irshad') | ioverride_read_weather_Irshad==0;
    draw_square=1;
    iproduce_CTH=1;
end

%18 km north, south east and west of the station 62.909,27.656
%corresponds to 63.07 and 62.747 N and 27.3 and 28.01 E


switch iproduce_CTH
    case 1
        clear Tin Tout RHin RHout distCTT T_CB
        
        dtol = 0.5/24; %weather data is every hour, so make this +/- 30 mins
        f = 1.6094e+06;
        for i=1:length(Cloud_Top_Temperature_Day_Mean.timeseries3)
            [minval,imin] = min( abs( weather(:,1) -  Date_Time_Swath.timeseries3(i)) );
            if minval<dtol
                Tin(i) = weather(imin,9) + 273.15; %K
                Tout(i) = weather(imin,8) + 273.15; %K          
                RHin(i) = weather(imin,11); % RH in %
                RHout(i) = weather(imin,10); % RH in %
                
                T = Tout(i);
                P = 980e2; 
                RH = RHout(i);
                
%                 qsat = SatVapPress(T,'goff','liq',P,1)/f; %kg/kg
% %                  th_parcel = potemp(T,P);
%                   q_parcel = RH/100 * qsat;
%         
% %         pgrid_parcel = [1080:-1:100];
%     %Tad_parcel=th_parcel./(1000./pgrid_parcel).^0.286; %calculate the temperatures during a constant dry adiabatic ascent    
%     t_dew_parcel=Tdew(q_parcel,P); %the dew point temperature for our q
%     %[minT_parcel imin_parcel]=min(abs(Tad_parcel-273.15-t_dew_parcel));
%     %cb_pressure_parcel = pgrid_parcel(imin_parcel)
%     
%     LR_dry = 9.81/1005;  % K/m

               if T>173    %seems like there are some points with T = 0 K...
                   [distCTT(i),T_CB(i)] = dist_to_CTT(T,P,RH,Cloud_Top_Temperature_Day_Mean.timeseries3(i));
               else
                   distCTT(i) = NaN;
                   T_CB(i) = NaN;
               end
    
    
                
                
            end            
        end
        
      
        
        %now estimate CTH from CTT and ground temperature
        
        
        
end

switch draw_square
    case 1
        %m_line draws one line per column
%         box_lines_lon = [27.3 28.01; 28.01 28.01; 28.01 27.3; 27.3 27.3]';
%         box_lines_lat = [63.07 63.07; 63.07 62.747; 62.747 62.747; 62.747 63.07]';        
%         m_line(box_lines_lon,box_lines_lat,'color','w');


swath_index=2; %used to match MOD06_L2.A2006260.1055.051.2010310011457.hdf (when had 
% load('/home/disk/eos8/d.grosvenor/saved_data_L2/Puijo_Sami/Puijo/terra/Puijo_terra_2006_L2_Puijo_matches_20130315T082055.mat')

swath_index=2;
%load('/home/disk/eos8/d.grosvenor/saved_data_L2/Puijo_Sami/Puijo/aqua/Puijo_aqua_2006_L2_Puijo_matches_20130315T121530.mat');


dt3=(Date_Time_Swath.timeseries3(swath_index) - (weather(weather_index.timeseries3(swath_index)-5:weather_index.timeseries3(swath_index),1)+0.5/24))*24*3600;
dtt=-diff([dt3; 0]);
wind_speed2=weather(weather_index.timeseries3(swath_index)-5:weather_index.timeseries3(swath_index),13);
wind_dir2=weather(weather_index.timeseries3(swath_index)-5:weather_index.timeseries3(swath_index),19);
[u2,v2]=uv_from_winddir(wind_speed2,wind_dir2);
dx3 = cumsum(flipdim(u2.*dtt,1));
dy3 = cumsum(flipdim(v2.*dtt,1));
[lon_new3,lat_new3]=calc_lat_lon_change_for_dx_dy(27.656,62.909,dx3,dy3);

lat_Puj = 62.909;
lon_Puj = 27.656;

m_plot(lon_new3,lat_new3,'ws');
m_plot(lon_Puj,lat_Puj,'kx'); m_plot(lon_Puj,lat_Puj,'wo');
        
cen_lat = lat_new3(1);
cen_lon = lon_new3(1);
dlat = 0.5;
dlon = 1;
         %m_line draws one line per column
        box_lines_lon = [cen_lon-dlon/2 cen_lon+dlon/2; cen_lon+dlon/2 cen_lon+dlon/2; cen_lon+dlon/2 cen_lon-dlon/2; cen_lon-dlon/2 cen_lon-dlon/2]';
        box_lines_lat = [cen_lat+dlat/2 cen_lat+dlat/2; cen_lat+dlat/2 cen_lat-dlat/2; cen_lat-dlat/2 cen_lat-dlat/2; cen_lat-dlat/2 cen_lat+dlat/2]';        
        m_line(box_lines_lon,box_lines_lat,'color','w');
        
   iav = find( Plat2D >= cen_lat-dlat/2 & Plat2D <= cen_lat+dlat/2 & Plon2D >= cen_lon-dlon/2 & Plon2D <= cen_lon+dlon/2 );
   
   [Pav_box,N_box,std_box] = meanNoNan(P(iav),1)
        
        
end
