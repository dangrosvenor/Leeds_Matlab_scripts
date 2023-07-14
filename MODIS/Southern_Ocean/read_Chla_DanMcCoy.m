% N.B. -this doesn't seem to work on Matlab 2007b (7.5.0)
% So, on Challenger, load using matlab-r2013a
%Read daily Chla files from Dan McCoy

%Are two different file types that use different averaging alogirithms, e.g. :-
%  L3m_20141219__GLOB_100_GSM-MODVIR_CHL1_DAY_00.nc
%  L3m_20141219__GLOB_100_AVW-MODVIR_CHL1_DAY_00.nc

% Looks like both combine MODIS and VIRS obs
% But the AVW uses weighted averaging and the GSM one some other method
% (Garver-Siegel-Maritorena method - ?)
% Actually, reading the documentation found here - http://www.globcolour.info/CDR_Docs/GlobCOLOUR_PUG.pdf
% it seems that it uses MODIS (Aqua), MERIS, SeaWifs and VIIRS.

file_type_chla = 'AVW'; %or GSM

chla_dir = '/home/disk/eos7/dtmccoy/research/GLOBCOLOUR/dataset_1_day/';

save_file_Chla = ['/home/disk/eos8/d.grosvenor/Chla/Daily_saved_Chla_' file_type_chla '_DJF_2006-2007_' datestr(now,30) '.mat'];

%Will plot Dec 2006 to end Feb 2007 for now.
%year_range = [2006 2007];
%day_range = {[335:365],[1:60]};

%All of 2006
year_range = [2006];
day_range = {[1:365]};

isave_plot=1;
iplot_polar=1; %Whether to plot in polar stereographic projection

seaice_fileload_1deg_global ='/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2007_global_20140508T094510.mat';
load(seaice_fileload_1deg_global);
seaice_array_datenum = seaice_array_1deg_datenum;

% N.B. will need to run 
%   seaice_match_times 
% afterwards to match the times to the MODIS data


daynum_timeseries3=[]; modisyear_timeseries3=[];
inds=0;
for iyear_chla=1:length(year_range)
    daynum_timeseries3 = [daynum_timeseries3 day_range{iyear_chla}];
    inds = [inds(end)+1:inds(end)+length(day_range{iyear_chla})];
    modisyear_timeseries3(inds) = year_range(iyear_chla);
end
daynum_timeseries3_Chla = daynum_timeseries3;
modisyear_timeseries3_Chla = modisyear_timeseries3;

iday_linear = 0;
for iyear_chla=1:length(year_range)
    
    day_range2 = day_range{iyear_chla};
    
    for iday_clha=1:length(day_range2)
        iday_linear = iday_linear + 1;
        
        day_chla = day_range2(iday_clha);
        year_chla = year_range(iyear_chla);
        
        date_str = date_from_day_of_year_func(day_chla,year_chla);
        [Y,M,D] = datevec(date_str);
        
        day_str = num2str(D,'%02g');
        mon_str = num2str(M,'%02g');
        year_str = num2str(Y);
        
        file_chla_dir = dir([chla_dir 'L3m_' year_str mon_str day_str '__GLOB_100_' file_type_chla '*.nc']);
        
        if length(file_chla_dir)>0
            
        file_chla = [chla_dir file_chla_dir(1).name];
        
        %       nc_chla = netcdf(file_chla);
        
        % Using the Matlab NetCDF reader
        lat_chla = ncread(file_chla,'lat'); MLAT = double(lat_chla');
        lon_chla = ncread(file_chla,'lon'); MLON = double(lon_chla');
        chla_mean = ncread(file_chla,'CHL1_mean');
        chla_error = ncread(file_chla,'CHL1_error');
        
        dat_modis = double(chla_mean');
        
        if iday_linear == 1;
           Cloud_Fraction_Liquid.timeseries3 = NaN*ones([size(dat_modis) length(daynum_timeseries3)]);
           chla_3D = NaN*ones([size(dat_modis) length(daynum_timeseries3)]);
           chla_error_3D = NaN*ones([size(dat_modis) length(daynum_timeseries3)]);
           seaice_match_times
        end
        
        %save the data in a big array
        chla_3D(:,:,iday_linear) = dat_modis;
        chla_error_3D(:,:,iday_linear) = double(chla_error');
        
        %Making a "new" gcm_str to avoid using timeseries3 type data in
        %pdf2D
        %         gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
        %         gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
        %         gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
        %         gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;
        %Actually can't really define this it varies over the globe for
        %ASMRE
        %         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre);
        gcm_time_matlab_AMSRE2 = 0;
        gcm_time_UTC_AMSRE2 = 0;
        daynum_timeseries3_AMSRE2 = 1;
        modisyear_timeseries3_AMSRE2 = 1;
        
        gcm_time_matlab_MODIS = 0;
        gcm_time_UTC_MODIS = 0;
%        daynum_timeseries3_MODIS = 1;
        %modisyear_timeseries3 = 1;
        
        %--- run the file to set up the defaults
        plot_global_maps_defaults
        
        %--- set some options for these particular plot loops

        
        cont_dat=flipdim(seaice_time3(:,:,iday_linear),1);
        cont_dat(isnan(cont_dat))=0; %Doing this has the effect of putting
        %in the coastline too - but this makes it look like the ice-shelf
        %regions might be ok for retreivals (there is no sea-ice, but
        %retrievals are still bad)
        %Actually decided to leave it since otherwise get strange bits of
        %contour.
        
        %Sea-ice data is in opposite direction to CHLA
        cont_ints = [0.01 0.01];
        icontour=1;
        contour_case='Specified data';
        cont_col_str='c';  %N.B. - 'w' colour doesn't show up when saving for some reason?
        %Use 'c' (cyan) instead
        
        if iplot_polar==1
            proj_type='polar';
            stereo_str01='lat_polar=-90;'; %Antarctic
            stereo_str02='m_proj(''stereographic'',''lat'',lat_polar,''lon'',0,''rad'',50);';
        else
            irestrict_domain=1;
            thresh_LAT=[-90 -40];
            thresh_LON=[-180 180];
        end
        
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';
        
        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=1;
        
        
        
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';
        
        
        titlenam_driver = ['Chla for ' date_str file_type_chla ' '];
        units_str_plot = 'mg m^{-3}';
        
        %mod_data_type='AMSRE';
        mod_data_type='timeseries3 lambert';
        gcm_str_select='AMSRE2';
        
        
        
        %        i_dpcolor=1;
        ifull_swath=0;
        igcm_screen=0;
        
        
        
        
        %--- Apply override flags
        ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
        % iocean_only=1;
        ioverride_time_selection=0; %Override the times to include
        ioverride_plotglobal_loc=1; %Override the location of the plot window
        ioverride_years_time_screen=0; %Override years for screening?
        
        %---  Run plot script and save
        plot_global_maps
        
        %Add a box for the UM simulation region
        %plot_box_on_map
        
        
        
        %Save the figure
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
        
        end  %if length        
        
        
    end
    
    
    
end



seaice_time3_Chla=flipdim(seaice_time3,1);
save(save_file_Chla,'-V7.3','chla_3D','chla_error_3D','daynum_timeseries3_Chla', ...
    'modisyear_timeseries3_Chla','seaice_time3_Chla','MLAT','MLON' ...
    )


% AVW ncdump :-

% netcdf L3m_20141219__GLOB_100_AVW-MODVIR_CHL1_DAY_00 {
% dimensions:
%         lat = 180 ;
%         lon = 360 ;
% variables:
%         float lat(lat) ;
%                 lat:long_name = "latitude" ;
%                 lat:units = "degrees_north" ;
%                 lat:axis = "Y" ;
%         float lon(lon) ;
%                 lon:long_name = "longitude" ;
%                 lon:units = "degrees_east" ;
%                 lon:axis = "X" ;
%         float CHL1_mean(lat, lon) ;
%                 CHL1_mean:standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water" ;
%                 CHL1_mean:long_name = "Chlorophyll concentration - Mean of the binned pixels" ;
%                 CHL1_mean:_FillValue = -999.f ;
%                 CHL1_mean:units = "mg/m3" ;
%                 CHL1_mean:pct_characterised_error = 43.31f ;
%         short CHL1_flags(lat, lon) ;
%                 CHL1_flags:long_name = "Chlorophyll concentration - Flags" ;
%                 CHL1_flags:_FillValue = 0s ;
%         short CHL1_error(lat, lon) ;
%                 CHL1_error:long_name = "Chlorophyll concentration - Error estimation" ;
%                 CHL1_error:_FillValue = -32768s ;
%                 CHL1_error:units = "0.01%" ;
%
% // global attributes:
%                 :Conventions = "CF-1.4" ;
%                 :title = "GlobColour daily merged MODIS/VIIRSN product" ;
%                 :product_name = "L3m_20141219__GLOB_100_AVW-MODVIR_CHL1_DAY_00.nc" ;
%                 :product_type = "day" ;
%                 :product_version = "2014.0" ;
%                 :product_level = 3s ;
%                 :parameter_code = "CHL1" ;
%                 :parameter = "Chlorophyll concentration" ;
%                 :parameter_algo_list = "OC3v5,OC3v5" ;
%                 :site_name = "GLOB" ;
%                 :sensor_name = "WEIGHTED_AVERAGING" ;
%                 :sensor = "Merged data - weighted mean" ;
%                 :sensor_name_list = "MOD,VIR" ;
%                 :start_time = "20141219T010509Z" ;
%                 :end_time = "20141220T062650Z" ;
%                 :duration_time = "PT105702S" ;
%                 :period_start_day = "20141219" ;
%                 :period_end_day = "20141219" ;
%                 :period_duration_day = "P1D" ;
%                 :grid_type = "Equirectangular" ;
%                 :spatial_resolution = 111.3195f ;
%                 :nb_equ_bins = 360 ;
%                 :registration = 5 ;
%                 :lat_step = 1.f ;
%                 :lon_step = 1.f ;
%                 :earth_radius = 6378.137 ;
%                 :max_north_grid = 90.f ;
%                 :max_south_grid = -90.f ;
%                 :max_west_grid = -180.f ;
%                 :max_east_grid = 180.f ;
%                 :northernmost_latitude = 47.f ;
%                 :southernmost_latitude = -79.f ;
%                 :westernmost_longitude = -180.f ;
%                 :easternmost_longitude = 180.f ;
%                 :nb_grid_bins = 64800 ;
%                 :nb_bins = 64800 ;
%                 :pct_bins = 100.f ;
%                 :nb_valid_bins = 13384 ;
%                 :pct_valid_bins = 20.65432f ;
%                 :software_name = "globcolour_l3_reproject" ;
%                 :software_version = "2014.0" ;
%                 :institution = "ACRI" ;
%                 :processing_time = "20150217T012430Z" ;
%                 :netcdf_version = "4.3.2 of Jan 20 2015 16:19:57 $" ;
%                 :DPM_reference = "GC-UD-ACRI-PUG" ;
%                 :IODD_reference = "GC-UD-ACRI-PUG" ;
%                 :references = "http://www.globcolour.info" ;
%                 :contact = "service@globcolour.info" ;
%                 :copyright = "Copyright ACRI-ST - GlobColour. GlobColour has been originally funded by ESA with data from ESA, NASA, NOAA and GeoEye. This reprocessing version has received funding from the European Community’s Seventh Framework Programme ([FP7/2007-2013]) under grant agreement n° 282723 [OSS2015 project]." ;





% GSM ncdump :-
%
% netcdf L3m_20141219__GLOB_100_GSM-MODVIR_CHL1_DAY_00 {
% dimensions:
%         lat = 180 ;
%         lon = 360 ;
% variables:
%         float lat(lat) ;
%                 lat:long_name = "latitude" ;
%                 lat:units = "degrees_north" ;
%                 lat:axis = "Y" ;
%         float lon(lon) ;
%                 lon:long_name = "longitude" ;
%                 lon:units = "degrees_east" ;
%                 lon:axis = "X" ;
%         float CHL1_mean(lat, lon) ;
%                 CHL1_mean:standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water" ;
%                 CHL1_mean:long_name = "Chlorophyll concentration - Mean of the binned pixels" ;
%                 CHL1_mean:_FillValue = -999.f ;
%                 CHL1_mean:units = "mg/m3" ;
%                 CHL1_mean:pct_characterised_error = 43.31f ;
%         short CHL1_flags(lat, lon) ;
%                 CHL1_flags:long_name = "Chlorophyll concentration - Flags" ;
%                 CHL1_flags:_FillValue = 0s ;
%         short CHL1_error(lat, lon) ;
%                 CHL1_error:long_name = "Chlorophyll concentration - Error estimation" ;
%                 CHL1_error:_FillValue = -32768s ;
%                 CHL1_error:units = "0.01%" ;
%
% // global attributes:
%                 :Conventions = "CF-1.4" ;
%                 :title = "GlobColour daily merged MODIS/VIIRSN product" ;
%                 :product_name = "L3m_20141219__GLOB_100_GSM-MODVIR_CHL1_DAY_00.nc" ;
%                 :product_type = "day" ;
%                 :product_version = "2014.0" ;
%                 :product_level = 3s ;
%                 :parameter_code = "CHL1" ;
%                 :parameter = "Chlorophyll concentration" ;
%                 :parameter_algo_list = "GSM method" ;
%                 :site_name = "GLOB" ;
%                 :sensor_name = "GSM" ;
%                 :sensor = "Garver-Siegel-Maritorena" ;
%                 :sensor_name_list = "MOD,VIR" ;
%                 :start_time = "20141219T010509Z" ;
%                 :end_time = "20141220T062650Z" ;
%                 :duration_time = "PT105702S" ;
%                 :period_start_day = "20141219" ;
%                 :period_end_day = "20141219" ;
%                 :period_duration_day = "P1D" ;
%                 :grid_type = "Equirectangular" ;
%                 :spatial_resolution = 111.3195f ;
%                 :nb_equ_bins = 360 ;
%                 :registration = 5 ;
%                 :lat_step = 1.f ;
%                 :lon_step = 1.f ;
%                 :earth_radius = 6378.137 ;
%                 :max_north_grid = 90.f ;
%                 :max_south_grid = -90.f ;
%                 :max_west_grid = -180.f ;
%                 :max_east_grid = 180.f ;
%                 :northernmost_latitude = 47.f ;
%                 :southernmost_latitude = -79.f ;
%                 :westernmost_longitude = -180.f ;
%                 :easternmost_longitude = 180.f ;
%                 :nb_grid_bins = 64800 ;
%                 :nb_bins = 64800 ;
%                 :pct_bins = 100.f ;
%                 :nb_valid_bins = 13155 ;
%                 :pct_valid_bins = 20.30093f ;
%                 :software_name = "globcolour_l3_reproject" ;
%                 :software_version = "2014.0" ;
%                 :institution = "ACRI" ;
%                 :processing_time = "20150217T012429Z" ;
%                 :netcdf_version = "4.3.2 of Jan 20 2015 16:19:57 $" ;
%                 :DPM_reference = "GC-UD-ACRI-PUG" ;
%                 :IODD_reference = "GC-UD-ACRI-PUG" ;
%                 :references = "http://www.globcolour.info" ;
%                 :contact = "service@globcolour.info" ;
%                 :copyright = "Copyright ACRI-ST - GlobColour. GlobColour has been originally funded by ESA with data from ESA, NASA, NOAA and GeoEye. This reprocessing version has received funding from the European Community’s Seventh Framework Programme ([FP7/2007-2013]) under grant agreement n° 282723 [OSS2015 project]." ;
