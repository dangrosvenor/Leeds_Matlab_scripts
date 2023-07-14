%output the output from read_Matts_precip_cloudsat in a nicer form to give
%to Matt W

save_file_MattW = '/home/disk/eos1/d.grosvenor/cloudsat_data_for_MattW/CloudSat_precip_global_2007-2011.mat';

% the code below is from plot_global_maps to see how data works

%   units_str_plot='mm hr^{-1}'; 
% 
%             
%             dat_modis = rain_warm_ocean;
% 
% switch asc_desc
%                 case 'ascending'                   
%                     dat_modis = dat_modis(1:2:end,:,:);
% %                    time_inds_average = time_inds_average(1:2:end);
%                     iasc_desc=1; %set this in case want to run case 115 of waterVap (to get the right
%                     %time inds in time_inds_modisL3_timeseries3
%                 case 'descending'
%                     dat_modis = dat_modis(2:2:end,:,:);    
% %                    time_inds_average = time_inds_average(2:2:end);
%                     iasc_desc=2;
%             end

% [Plon2D_matt_edges,Plat2D_matt_edges] = meshgrid(lon_matt_edges,lat_matt_edges);
% 
% lat_matt_centres = [-89:2:89]; %90 faces
% lon_matt_centres = [-178:4:178]; %91 edges
% %i180=find(lon_matt_centres>180);
% %lon_matt_centres(i180)=lon_matt_centres(i180)-360;
% 
% [Plon2D_matt_centres,Plat2D_matt_centres] = meshgrid(lon_matt_centres,lat_matt_centres);


%want the data for Oct/Nov all years
cloudsat_month = month_calipso_matt(1:2:end);
cloudsat_year = year_calipso_matt(1:2:end);
cloudsat_precip_daytime = rain_warm_ocean(1:2:end,:,:);
cloudsat_precip_nighttime = rain_warm_ocean(2:2:end,:,:);
cloudsat_lat2D_centres = Plat2D_matt_centres;
cloudsat_lon2D_centres = Plon2D_matt_centres;
cloudsat_lat2D_edges = Plat2D_matt_edges;
cloudsat_lon2D_edges = Plon2D_matt_edges;


save(save_file_MattW,'cloudsat_month','cloudsat_year','cloudsat_precip_daytime','cloudsat_precip_nighttime',...
    'cloudsat_lat2D_centres','cloudsat_lon2D_centres','cloudsat_lat2D_edges','cloudsat_lon2D_edges','-V7.3');






