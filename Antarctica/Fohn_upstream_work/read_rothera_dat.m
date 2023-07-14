rothera_file = '/home/disk/eos1/d.grosvenor/Antarctica/Rothera_data_from_John_for_ref_comments - flight19_r.txt'

fid=fopen(rothera_file,'rt'); %is quicker to use fscanf that dlmread
dat_rothera_IN = fscanf(fid,'%f',[8 inf]);
fclose(fid);

% Data is arranged as [variable time]
% Var 1  = year
% Var 2  = month
% Var 3  = day
% Var 4  = hour
% Var 5  = mslp
% Var 6  = temperature
% Var 7  = wind speed (knots)
% Var 8  = wind direction

dat_rothera = dat_rothera_IN;

dat_rothera(7,:) = dat_rothera_IN(7,:)*0.514444444; %convert from knots to m/s
time_roth_matlab = datenum(dat_rothera(1,:),dat_rothera(2,:),dat_rothera(3,:),dat_rothera(4,:),0,0);

