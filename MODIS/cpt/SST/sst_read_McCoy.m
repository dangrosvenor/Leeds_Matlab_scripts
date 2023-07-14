%Load data - note, there is also an SST2.mat file - not sure what this is?
%load('/home/disk/eos7/dtmccoy/research/SST/SST.mat');
fn_dir = '/home/disk/eos8/d.grosvenor/'; %contains up to Oct 2015
%fn_dir = '/home/disk/eos15/d.grosvenor/SSTs/' %updated to run to 30 Apr, 2017
% from this file fn=[fn_dir 'sst.wkmean.1990-present.nc'];
save_name = [fn_dir 'SST2.mat'];
sst_load = load(save_name);

%The Matlab file SST.mat contains the read in version with an SST array that is 360x180x1253. The 1253 is 1253/52 = 24.0962 years, which makes sense since is labelled as 1990 to present (i.e. 24 years if present is 2014).
%This is the same size as for the .nc file (sst.wkmean.1990-present.nc). Time convention is 
%  days since 1800-1-1 00:00:00
% lat and lon are 1x1 deg.

sst_time = sst_load.time_matlab; %
lat_sst = sst_load.lat;
lon_sst = sst_load.lon;

SST = 0.01*sst_load.SST; %Convert to degC

