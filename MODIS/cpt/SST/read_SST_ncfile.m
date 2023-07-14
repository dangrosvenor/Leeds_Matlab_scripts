function[]=read_SST_ncfile()
% Converted from Daniel McCoy's read_process_SST file in
% /home/disk/eos7/dtmccoy/research/SST/
fn_dir = '/home/disk/eos8/d.grosvenor/';
fn=[fn_dir 'sst.wkmean.1990-present.nc'];
save_name = [fn_dir 'SST2.mat'];
nc = netcdf(fn);

%fn='sst.mnmean.nc';
[SST]=nc{'sst'}(:);
lat=nc{'lat'}(:);
lon=nc{'lon'}(:);
time=nc{'time'}(:);

time_matlab = time + datenum('01-Jan-1800'); %This gives the first day as 31st Dec 1989 - is this right?
    % Yes - see http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html

nc2 = netcdf([fn_dir 'lsmask.nc']); lsmask=nc2{'mask'};
nc3 = netcdf([fn_dir 'icec.mnmean.nc']); icec=nc3{'icec'}; timeice=nc3{'time'};

% [SST]=ncread(fn,'sst');
% lat=ncread(fn,'lat');
% lon=ncread(fn,'lon');
% time=ncread(fn,'time');
% lsmask=ncread('lsmask.nc','mask');
% icec=ncread('icec.mnmean.nc','icec');
% timeice=ncread('icec.mnmean.nc','time');
%save('SST')
lon(lon>180)=lon(lon>180)-360;
%% first, we go through and get the times

addpath /home/disk/eos7/dtmccoy/research/scripts/
doy=time*NaN;
YMD=[time time time]*NaN;

for i=1:length(time)
[doy(i),YMD(i,:),decy(i)]=dayssince2doy(time(i),[1800 1 0]);
end
for i=1:length(timeice)
[doyice(i),YMDice(i,:),decyice(i)]=dayssince2doy(timeice(i),[1800 1 0]);
end
%save('SSTMON')
save(save_name)
end
