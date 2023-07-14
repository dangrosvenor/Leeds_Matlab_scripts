rothera_file = '/home/disk/eos1/d.grosvenor/Antarctica/Radiosonde_hires_data_from_JK_20060105.txt';
%rothera_file = '/home/disk/eos1/d.grosvenor/Antarctica/Radiosonde_hires_data_from_JK_20060106.txt';

fid=fopen(rothera_file,'rt'); %is quicker to use fscanf that dlmread
data_rs = fgetl(fid);
dat_rs_IN = fscanf(fid,'%f',[9 inf]);
fclose(fid);
dat_rs_IN = dat_rs_IN';

% Data is arranged as [variable time]
% Var 1  = Time (mins?)
% Var 2  = Time (secs?)
% Var 3  = Pressure (hPa)
% Var 4  = Height (m)
% Var 5  = Temperature?
% Var 6  = RH?
% Var 7  = Dew point Temperature?
% Var 8  = wind direction
% Var 9  = wind speed (knots?)

dat_rs = dat_rs_IN;

dat_rs(7,:) = dat_rs_IN(7,:)*0.514444444; %convert from knots to m/s
%time_rs_matlab = datenum(dat_rothera(1,:),dat_rothera(2,:),dat_rothera(3,:),dat_rothera(4,:),0,0);

