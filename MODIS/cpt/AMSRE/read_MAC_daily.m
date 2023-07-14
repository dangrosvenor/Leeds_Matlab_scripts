function [lwp_out,tlwp_out,time_out,lat,lon,lat2d,lon2d] = read_MAC_daily(mac_path)

[lwp_out,time_out,lat,lon] = mac_read_var(mac_path,'LWP');
[tlwp_out,time_out,lat,lon] = mac_read_var(mac_path,'TLWP');

lon = [-179.5:179.5];
lat = [89.5:-1:-89.5];

[lon2d, lat2d] = meshgrid(lon,lat);



% fnames=dir([mac_path '/tlwp*.nc']);
% tlwp_out=[];
% time=[];
% for i=1:length(fnames)
%     file_path = [mac_path '/' fnames(i).name];
%     nc=netcdf(file_path);     
%     
%     tlwp_out_month = nc{'totallwp'}(:);    
%     tlwp_out = cat(1,tlwp_out,tlwp_out_month/1e3);    
% end
% 

function [lwp_out,time_out,lat,lon] = mac_read_var(mac_path,varname)

switch varname
    case 'LWP'
        nc_var='cloudlwp';
        fname_str = 'clwp';
    case 'TLWP'  
        nc_var='totallwp';
        fname_str = 'tlwp';
end


fnames=dir([mac_path '/' fname_str '*.nc']);
lwp_out=[];
time_out=[];
for i=1:length(fnames)
    file_path = [mac_path '/' fnames(i).name];
    nc=netcdf(file_path);
    
    istr=findstr('.nc',file_path);
    month_str = file_path(istr-2:istr-1);
    month = str2num(month_str);
    year_str = file_path(istr-6:istr-3);
    year = str2num(year_str);

    time_read=nc{'time'}(:);
    time = datenum(year,month,1+time_read); %Matlab time
    lat=nc{'lat'}(:);
    lon=nc{'lon'}(:);
    
    lwp_read = nc{nc_var}(:);
    inan = find(lwp_read<-99.8);
    lwp_read(inan)=NaN;
    
    %Make it from -180 to 180 as for MODIS
    dat_west = lwp_read(:,:,1:180);
    dat_east = lwp_read(:,:,181:end);        
    lwp_out_month = cat(3,dat_east,dat_west);   
    lwp_out_month = flipdim(lwp_out_month,2); %flip lat lon to be
    %consistent with AMSR-E
    
    lwp_out = cat(1,lwp_out,lwp_out_month/1e3); %convert to kg/m2 for consistency with AMSRE
    time_out = cat(1,time_out,time);
end







