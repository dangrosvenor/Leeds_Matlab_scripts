%filename='/home/mbexddg5/work/manchester_flt19.txt';
filename='c:/documents and settings/dan/my documents/Antarctica_general/Rothera_Feb09_flights/flight90_part1';

fid=fopen(filename,'rt');

dat = fscanf(fid,'%f');

dat_flt90=( reshape(dat,[16 size(dat,1)/16]) )';


%time_flt19=dat_flt19(:,1)/3.6e6; %convert from ms to hours

%td_flt19=dat_flt19(inds,7)+273.15;  %frost point instrument ***wasn't working properly***
%td2_flt19=dat_flt19(inds,8)+273.15; %so use this humicap data

f=1e6*28.97/18; 

fclose(fid);

%qv_flt19 = SatVapPress(td_flt19,'goff','liq',100*dat_flt19(inds,4),1) / f;
%qv2_flt19 = SatVapPress(td2_flt19,'goff','liq',100*dat_flt19(inds,4),1) / f;

disp('Done Feb09 reading aircraft data');


% 1) time
% 2) lat
% 3) lon
% 4) static pressure
% 5) infra red thermometer
% 6) normal temp
% 7) frost point hygrometer
% 8) dew point from humicap
% 9) wind speed
% 10) wind direction
% 11) radar altimeter
% 12) gps height
% 13) sw radiometer out
% 14) sw in
% 15) ir out
% 16) ir in

filename='c:/documents and settings/dan/my documents/Antarctica_general/Rothera_Feb09_flights/flight90_part2';

fid=fopen(filename,'rt');

dat = fscanf(fid,'%f');

dat_flt90_2=( reshape(dat,[6 size(dat,1)/6]) )';

% 1) time
% 2) u
% 3) v
% 4) w
% 5) fast temp1
% 6) fast temp2

fclose(fid);




