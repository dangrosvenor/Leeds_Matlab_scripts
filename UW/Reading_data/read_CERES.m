% This can be used to read in both the Aqua and Terra data (concatenates
% them.)

%Read CERES data
% First data was for VOCALS case & stored in  /home/disk/eos8/d.grosvenor/CERES/12thNov_VOCALS/

grid_data=0; %Whether to grid the irregular data onto a regular lat/lon grid. Takes a little time, but might be useful for plotting
%this data as a map. However, there is some distortion near the edge of
%swaths. Could perhaps remove data with the higher VZA values to get rid of
%this?
vza_thresh=65.5;
%vza_thresh=9e9;

plotting_stuff_to_move=0; %Now have the Boutle*SWLW*Map script

data_dir = '/home/disk/eos8/d.grosvenor/CERES/12thNov_VOCALS/';

clear filename  %N.B. Reads in ALL of the filenames left uncommented
ifile=1;
filename{ifile} = 'CERES_SSF_XTRK-MODIS_Edition4A_Subset_2008111200-2008111422.nc'; ifile=ifile+1;
filename{ifile} = 'TERRA/CERES_SSF_XTRK-MODIS_Edition4A_Subset_2008111202-2008111421.nc'; ifile=ifile+1;

for itime_read_CERES=1:length(filename)

    filepath = [data_dir filename{itime_read_CERES}];
    nc = netcdf(filepath);


    time_CERES = nc{'time'}(:); %"days since 1970-01-01 00:00:00"
    %Data is stored in one long vector with different times, lons, lats, etc.

    %conver to Matlab time
    time_CERES = time_CERES + datenum('01-Jan-1970');

    lat_CERES = nc{'lat'}(:);
    lon_CERES = nc{'lon'}(:);

    %% Separate all the swaths (south-to-north segments)
    %  Will be a big jump in latitude when move to the next one
    %dlat_max = maxALL(lat_CERES) - minALL(lat_CERES);
    %dlat_new_swath = 2; %If jumps more than 2 degrees then assume is a new swath - might not work well at teh poles, though

    %Use the jump in time
    dt=diff(time_CERES);
    dmins=1; %jump in time in hours required for new swath - set to 1 min since the gap between scanlines should be much smaller.
    iswaths=find(dt>dmins/60/24); %location of the start of the new swaths
    iswaths=cat(1,iswaths,length(time_CERES)); %adds the last position to the vector

    %itime=find(  abs(time_CERES - datenum('13-Nov-2008 18:45')) < 0.5/24 );
    %iloc = find(  abs(lat_CERES - -20 ) < 0.5 & abs(lon_CERES - -86 ) < 10);

    SW_TOA_up_CERES = nc{'CERES_SW_TOA_flux___upwards'}(:);
    SW_TOA_down_CERES = nc{'TOA_Incoming_Solar_Radiation'}(:);
    LW_TOA_up_CERES = nc{'CERES_LW_TOA_flux___upwards'}(:);

    VZA_CERES = nc{'CERES_viewing_zenith_at_surface'}(:);



    vars_CERES = {'lat_CERES','lon_CERES','VZA_CERES','SW_TOA_up_CERES','SW_TOA_down_CERES','LW_TOA_up_CERES','time_CERES'};
    vars_CERES_reg = vars_CERES(3:end);

    dlat=0.25; dlon=0.25;
    %Make a regular grid
    lat_reg = [minALL(lat_CERES):dlat:maxALL(lat_CERES)];
    lon_reg = [minALL(lon_CERES):dlon:maxALL(lon_CERES)];
    [lats,lons]=meshgrid(lat_reg,lon_reg);

    gcm_Plat2D_CERES=lats;
    gcm_Plon2D_CERES=lons;


    lat_reg_edges = [minALL(lat_CERES)-dlat/2:dlat:maxALL(lat_CERES)+dlat/2];
    lon_reg_edges = [minALL(lon_CERES)-dlon/2:dlon:maxALL(lon_CERES)+dlon/2];
    [gcm_Plat2D_edges_CERES,gcm_Plon2D_edges_CERES] = meshgrid(lat_reg_edges,lon_reg_edges);
    gcm_str='CERES';

    daynum_timeseries3_CERES = 1;
    gcm_time_UTC_CERES = 1;
    month_amsre=1;
    year_amsre=1;

    for i=1:length(vars_CERES)
        eval_str = ['ii=find(' vars_CERES{i} '>3.4e38);']; eval(eval_str);
        if length(ii)>0
            eval_str = [vars_CERES{i} '(' vars_CERES{i} '>3.4e38) = NaN;'];eval(eval_str);
        end



        if itime_read_CERES>1  %If not the first iteration of the loop
            eval_str = [vars_CERES{i} '= cat(1,' vars_CERES{i} '_old,' vars_CERES{i} ');']; eval(eval_str);
        end
        eval_str = [vars_CERES{i} '_old = ' vars_CERES{i} ';']; eval(eval_str);

    end


end


if grid_data==1



    for i=1:length(vars_CERES_reg)


        dat=eval(vars_CERES_reg{i});

        eval(['clear ' vars_CERES_reg{i} '_reg']);
        clear time_swath

        i0=1;

        %       istart=10;
        %       i0=iswaths(istart-1)+1;
        %       for is=istart:length(iswaths)-1
        for is=1:length(iswaths)-1
            inds=i0:iswaths(is);
            time_swath(is)=time_CERES(i0);
            i0=iswaths(is)+1; %reset i0

            inan=find(isnan(dat)==0);
            inan=find(dat>0);
            if length(inan)>0

                dat2=griddata(lat_CERES(inds),lon_CERES(inds),dat(inds),lats,lons);

                %Trim off the edge data due to distortion, except for the VZA itself
                switch(vars_CERES_reg{i})
                    case 'VZA_CERES'
                    otherwise
                        dat2(VZA_CERES_reg{is}>vza_thresh)=NaN;
                end

                eval([vars_CERES_reg{i} '_reg{is} = dat2;']);

            else
                eval([vars_CERES_reg{i} '_reg{is}= NaN*ones(size(lats));']);
            end





        end


    end

end


%% Combining swaths, etc.

if plotting_stuff_to_move==1

    swaths_combine=[1:5]; %10pm LST nightime Terra daytime overpasses for 11th Nov
    swaths_combine=[6:11]; %10am Terra daytime overpasses for 12th Nov
    %   swaths_combine=[12:17]; %10pm Terra nightime overpasses for 12th Nov
    %   swaths_combine=[18:22]; %10am Terra daytime overpasses for 13th Nov

    %Aqua
    %   swaths_combine=[4:8]; %2am on 12th
    %   swaths_combine=[9:13]; %2pm on 12th
    %   swaths_combine=[14:18]; %2am on 13th
    %N.B. - 2pm on 13th is not available
    %   swaths_combine=[21:23]; %2am on 14th - NO data for the near-coastal region east of around 100W
    %   swaths_combine=[24:26]; %2pm on 14th - NO data further west than around 120W

    [SW_CERES_combined]=combine_overlapping_data_with_Nans(SW_TOA_up_CERES_reg,swaths_combine);
    [LW_CERES_combined]=combine_overlapping_data_with_Nans(LW_TOA_up_CERES_reg,swaths_combine);
    [time_CERES_combined]=combine_overlapping_data_with_Nans(time_CERES_reg,swaths_combine);

    times_swaths_CERES = time_swath(swaths_combine); %all the times

end



% netcdf CERES_SSF_XTRK-MODIS_Edition4A_Subset_2008111200-2008111422 {
% dimensions:
%         time = UNLIMITED ; // (416372 currently)
%         The_8_most_prevalent_surface_types = 8 ;
%         Conditions_clear__lower__upper__upper_over_lower = 4 ;
%         Lower_and_Upper_Cloud_Layers = 2 ;
% variables:
%         double time(time) ;
%                 time:long_name = "time" ;
%                 time:units = "days since 1970-01-01 00:00:00" ;
%                 time:_FillValue = 1.79769313486232e+308 ;
%                 time:valid_range = 0., 39412.5 ;
%         float lon(time) ;
%                 lon:long_name = "longitude" ;
%                 lon:units = "degrees_east" ;
%                 lon:_FillValue = 3.402823e+38f ;
%                 lon:valid_range = -180.f, 180.f ;
%         float lat(time) ;
%                 lat:long_name = "latitude" ;
%                 lat:units = "degrees_north" ;
%                 lat:_FillValue = 3.402823e+38f ;
%                 lat:valid_range = -90.f, 90.f ;
%         double Time_of_observation(time) ;
%                 Time_of_observation:orig_name = "Time of observation" ;
%                 Time_of_observation:units = "day" ;
%                 Time_of_observation:format = "F18.9" ;
%                 Time_of_observation:_FillValue = 1.79769313486232e+308 ;
%                 Time_of_observation:valid_range = 2440000., 2480000. ;
%         float Longitude_of_CERES_FOV_at_surface(time) ;
%                 Longitude_of_CERES_FOV_at_surface:orig_name = "Longitude of CERES FOV at surface" ;
%                 Longitude_of_CERES_FOV_at_surface:units = "degrees" ;
%                 Longitude_of_CERES_FOV_at_surface:format = "F18.9" ;
%                 Longitude_of_CERES_FOV_at_surface:_FillValue = 3.402823e+38f ;
%                 Longitude_of_CERES_FOV_at_surface:valid_range = 0.f, 360.f ;
%         float Colatitude_of_CERES_FOV_at_surface(time) ;
%                 Colatitude_of_CERES_FOV_at_surface:orig_name = "Colatitude of CERES FOV at surface" ;
%                 Colatitude_of_CERES_FOV_at_surface:units = "degrees" ;
%                 Colatitude_of_CERES_FOV_at_surface:format = "F18.9" ;
%                 Colatitude_of_CERES_FOV_at_surface:_FillValue = 3.402823e+38f ;
%                 Colatitude_of_CERES_FOV_at_surface:valid_range = 0.f, 180.f ;
%         float Longitude_of_subsatellite_point_at_surface_at_observation(time) ;
%                 Longitude_of_subsatellite_point_at_surface_at_observation:orig_name = "Longitude of subsatellite point at surface at observation" ;
%                Longitude_of_subsatellite_point_at_surface_at_observation:_FillValue = 3.402823e+38f ;
%                 Longitude_of_subsatellite_point_at_surface_at_observation:valid_range = 0.f, 360.f ;
%         float Colatitude_of_subsatellite_point_at_surface_at_observation(time) ;
%                 Colatitude_of_subsatellite_point_at_surface_at_observation:orig_name = "Colatitude of subsatellite point at surface at observation" ;
%                 Colatitude_of_subsatellite_point_at_surface_at_observation:units = "degrees" ;
%                 Colatitude_of_subsatellite_point_at_surface_at_observation:format = "F18.9" ;
%                 Colatitude_of_subsatellite_point_at_surface_at_observation:_FillValue = 3.402823e+38f ;
%                 Colatitude_of_subsatellite_point_at_surface_at_observation:valid_range = 0.f, 180.f ;
%         float CERES_viewing_zenith_at_surface(time) ;
%                 CERES_viewing_zenith_at_surface:orig_name = "CERES viewing zenith at surface" ;
%                 CERES_viewing_zenith_at_surface:units = "degrees" ;
%                 CERES_viewing_zenith_at_surface:format = "F18.9" ;
%                 CERES_viewing_zenith_at_surface:_FillValue = 3.402823e+38f ;
%                 CERES_viewing_zenith_at_surface:valid_range = 0.f, 90.f ;
%         float CERES_relative_azimuth_at_surface(time) ;
%                 CERES_relative_azimuth_at_surface:orig_name = "CERES relative azimuth at surface" ;
%                 CERES_relative_azimuth_at_surface:units = "degrees" ;
%                 CERES_relative_azimuth_at_surface:format = "F18.9" ;
%                 CERES_relative_azimuth_at_surface:_FillValue = 3.402823e+38f ;
%                 CERES_relative_azimuth_at_surface:valid_range = 0.f, 360.f ;
%   float CERES_solar_zenith_at_surface(time) ;
% 		CERES_solar_zenith_at_surface:orig_name = "CERES solar zenith at surface" ;
% 		CERES_solar_zenith_at_surface:units = "degrees" ;
% 		CERES_solar_zenith_at_surface:format = "F18.9" ;
% 		CERES_solar_zenith_at_surface:_FillValue = 3.402823e+38f ;
% 		CERES_solar_zenith_at_surface:valid_range = 0.f, 180.f ;
% 	float CERES_SW_TOA_flux___upwards(time) ;
% 		CERES_SW_TOA_flux___upwards:orig_name = "CERES SW TOA flux - upwards" ;
% 		CERES_SW_TOA_flux___upwards:units = "Watts per square meter" ;
% 		CERES_SW_TOA_flux___upwards:format = "F18.9" ;
% 		CERES_SW_TOA_flux___upwards:_FillValue = 3.402823e+38f ;
% 		CERES_SW_TOA_flux___upwards:valid_range = 0.f, 1400.f ;
% 	float CERES_SW_radiance___upwards(time) ;
% 		CERES_SW_radiance___upwards:orig_name = "CERES SW radiance - upwards" ;
% 		CERES_SW_radiance___upwards:units = "Watts per square meter per steradian" ;
% 		CERES_SW_radiance___upwards:format = "F18.9" ;
% 		CERES_SW_radiance___upwards:_FillValue = 3.402823e+38f ;
% 		CERES_SW_radiance___upwards:valid_range = -10.f, 510.f ;
% 	float CERES_LW_TOA_flux___upwards(time) ;
% 		CERES_LW_TOA_flux___upwards:orig_name = "CERES LW TOA flux - upwards" ;
% 		CERES_LW_TOA_flux___upwards:units = "Watts per square meter" ;
% 		CERES_LW_TOA_flux___upwards:format = "F18.9" ;
% 		CERES_LW_TOA_flux___upwards:_FillValue = 3.402823e+38f ;
% 		CERES_LW_TOA_flux___upwards:valid_range = 0.f, 500.f ;
% 	float CERES_LW_radiance___upwards(time) ;
% 		CERES_LW_radiance___upwards:orig_name = "CERES LW radiance - upwards" ;
% 		CERES_LW_radiance___upwards:units = "Watts per square meter per steradian" ;
% 		CERES_LW_radiance___upwards:format = "F18.9" ;
% 		CERES_LW_radiance___upwards:_FillValue = 3.402823e+38f ;
% 		CERES_LW_radiance___upwards:valid_range = 0.f, 200.f ;
% 	float CERES_WN_TOA_flux___upwards(time) ;
% 		CERES_WN_TOA_flux___upwards:orig_name = "CERES WN TOA flux - upwards" ;
% 		CERES_WN_TOA_flux___upwards:units = "Watts per square meter" ;
% 		CERES_WN_TOA_flux___upwards:format = "F18.9" ;
% 		CERES_WN_TOA_flux___upwards:_FillValue = 3.402823e+38f ;
% 		CERES_WN_TOA_flux___upwards:valid_range = 0.f, 200.f ;
% 	float CERES_WN_radiance___upwards(time) ;
% 		CERES_WN_radiance___upwards:orig_name = "CERES WN radiance - upwards" ;
% 		CERES_WN_radiance___upwards:units = "Watts per square meter per steradian" ;
% 		CERES_WN_radiance___upwards:format = "F18.9" ;
% 		CERES_WN_radiance___upwards:_FillValue = 3.402823e+38f ;
% 		CERES_WN_radiance___upwards:valid_range = 0.f, 60.f ;
% 	float TOA_Incoming_Solar_Radiation(time) ;
% 		TOA_Incoming_Solar_Radiation:orig_name = "TOA Incoming Solar Radiation" ;
% 		TOA_Incoming_Solar_Radiation:units = "Watts per square meter" ;
% 		TOA_Incoming_Solar_Radiation:format = "F18.9" ;
% 		TOA_Incoming_Solar_Radiation:_FillValue = 3.402823e+38f ;
% 		TOA_Incoming_Solar_Radiation:valid_range = 0.f, 1400.f ;
% 	float CERES_downward_SW_surface_flux___Model_B(time) ;
% 		CERES_downward_SW_surface_flux___Model_B:orig_name = "CERES downward SW surface flux - Model B" ;
% 		CERES_downward_SW_surface_flux___Model_B:units = "Watts per square meter" ;
% 		CERES_downward_SW_surface_flux___Model_B:format = "F18.9" ;
% 		CERES_downward_SW_surface_flux___Model_B:_FillValue = 3.402823e+38f ;
% 		CERES_downward_SW_surface_flux___Model_B:valid_range = 0.f, 1400.f ;
% 	float CERES_downward_SW_surface_flux___Model_B__clearsky(time) ;
% 		CERES_downward_SW_surface_flux___Model_B__clearsky:orig_name = "CERES downward SW surface flux - Model B, clearsky" ;
% 		CERES_downward_SW_surface_flux___Model_B__clearsky:units = "Watts per square meter" ;
% 		CERES_downward_SW_surface_flux___Model_B__clearsky:format = "F18.9" ;
% 		CERES_downward_SW_surface_flux___Model_B__clearsky:_FillValue = 3.402823e+38f ;
% 		CERES_downward_SW_surface_flux___Model_B__clearsky:valid_range = 0.f, 1400.f ;
% 	float CERES_net_SW_surface_flux___Model_B(time) ;
% 		CERES_net_SW_surface_flux___Model_B:orig_name = "CERES net SW surface flux - Model B" ;
% 		CERES_net_SW_surface_flux___Model_B:units = "Watts per square meter" ;
% 		CERES_net_SW_surface_flux___Model_B:format = "F18.9" ;
% 		CERES_net_SW_surface_flux___Model_B:_FillValue = 3.402823e+38f ;
% 		CERES_net_SW_surface_flux___Model_B:valid_range = 0.f, 1400.f ;
% 	float CERES_downward_LW_surface_flux___Model_B(time) ;
% 		CERES_downward_LW_surface_flux___Model_B:orig_name = "CERES downward LW surface flux - Model B" ;
% 		CERES_downward_LW_surface_flux___Model_B:units = "Watts per square meter" ;
% 		CERES_downward_LW_surface_flux___Model_B:format = "F18.9" ;
% 		CERES_downward_LW_surface_flux___Model_B:_FillValue = 3.402823e+38f ;
% 		CERES_downward_LW_surface_flux___Model_B:valid_range = 0.f, 700.f ;
% 	float CERES_downward_LW_surface_flux___Model_B__clearsky(time) ;
% 		CERES_downward_LW_surface_flux___Model_B__clearsky:orig_name = "CERES downward LW surface flux - Model B, clearsky" ;
% 		CERES_downward_LW_surface_flux___Model_B__clearsky:units = "Watts per square meter" ;
% 		CERES_downward_LW_surface_flux___Model_B__clearsky:format = "F18.9" ;
% 		CERES_downward_LW_surface_flux___Model_B__clearsky:_FillValue = 3.402823e+38f ;
% 		CERES_downward_LW_surface_flux___Model_B__clearsky:valid_range = 0.f, 700.f ;
% 	float CERES_net_LW_surface_flux___Model_B(time) ;
% 		CERES_net_LW_surface_flux___Model_B:orig_name = "CERES net LW surface flux - Model B" ;
% 		CERES_net_LW_surface_flux___Model_B:units = "Watts per square meter" ;
% 		CERES_net_LW_surface_flux___Model_B:format = "F18.9" ;
% 		CERES_net_LW_surface_flux___Model_B:_FillValue = 3.402823e+38f ;
% 		CERES_net_LW_surface_flux___Model_B:valid_range = -250.f, 50.f ;
% 	float Altitude_of_surface_above_sea_level(time) ;
% 		Altitude_of_surface_above_sea_level:orig_name = "Altitude of surface above sea level" ;
% 		Altitude_of_surface_above_sea_level:units = "meters" ;
% 		Altitude_of_surface_above_sea_level:format = "F18.9" ;
% 		Altitude_of_surface_above_sea_level:_FillValue = 3.402823e+38f ;
% 		Altitude_of_surface_above_sea_level:valid_range = -1000.f, 10000.f ;
% 	short Surface_type_index(time, The_8_most_prevalent_surface_types) ;
% 		Surface_type_index:orig_name = "Surface type index" ;
% 		Surface_type_index:units = "N/A" ;
% 		Surface_type_index:format = "I10" ;
% 		Surface_type_index:_FillValue = 32767s ;
% 		Surface_type_index:valid_range = 1s, 20s ;
% 	short Surface_type_percent_coverage(time, The_8_most_prevalent_surface_types) ;
% 		Surface_type_percent_coverage:orig_name = "Surface type percent coverage" ;
% 		Surface_type_percent_coverage:units = "N/A" ;
% 		Surface_type_percent_coverage:format = "I10" ;
% 		Surface_type_percent_coverage:_FillValue = 32767s ;
% 		Surface_type_percent_coverage:valid_range = 0s, 100s ;
% 	float CERES_LW_surface_emissivity(time) ;
% 		CERES_LW_surface_emissivity:orig_name = "CERES LW surface emissivity" ;
% 		CERES_LW_surface_emissivity:units = "N/A" ;
% 		CERES_LW_surface_emissivity:format = "F18.9" ;
% 		CERES_LW_surface_emissivity:_FillValue = 3.402823e+38f ;
% 		CERES_LW_surface_emissivity:valid_range = 0.f, 1.f ;
% 	float CERES_WN_surface_emissivity(time) ;
% 		CERES_WN_surface_emissivity:orig_name = "CERES WN surface emissivity" ;
% 		CERES_WN_surface_emissivity:units = "N/A" ;
% 		CERES_WN_surface_emissivity:format = "F18.9" ;
% 		CERES_WN_surface_emissivity:_FillValue = 3.402823e+38f ;
% 		CERES_WN_surface_emissivity:valid_range = 0.f, 1.f ;
% 	float CERES_broadband_surface_albedo(time) ;
% 		CERES_broadband_surface_albedo:orig_name = "CERES broadband surface albedo" ;
% 		CERES_broadband_surface_albedo:units = "N/A" ;
% 		CERES_broadband_surface_albedo:format = "F18.9" ;
% 		CERES_broadband_surface_albedo:_FillValue = 3.402823e+38f ;
% 		CERES_broadband_surface_albedo:valid_range = 0.f, 1.f ;
% 	float Surface_wind___U_vector(time) ;
% 		Surface_wind___U_vector:orig_name = "Surface wind - U-vector" ;
% 		Surface_wind___U_vector:units = "meters per second" ;
% 		Surface_wind___U_vector:format = "F18.9" ;
% 		Surface_wind___U_vector:_FillValue = 3.402823e+38f ;
% 		Surface_wind___U_vector:valid_range = -100.f, 100.f ;
% 	float Surface_wind___V_vector(time) ;
% 		Surface_wind___V_vector:orig_name = "Surface wind - V-vector" ;
% 		Surface_wind___V_vector:units = "meters per second" ;
% 		Surface_wind___V_vector:format = "F18.9" ;
% 		Surface_wind___V_vector:_FillValue = 3.402823e+38f ;
% 		Surface_wind___V_vector:valid_range = -100.f, 100.f ;
% 	float Surface_skin_temperature(time) ;
% 		Surface_skin_temperature:orig_name = "Surface skin temperature" ;
% 		Surface_skin_temperature:units = "Kelvin" ;
% 		Surface_skin_temperature:format = "F18.9" ;
% 		Surface_skin_temperature:_FillValue = 3.402823e+38f ;
% 		Surface_skin_temperature:valid_range = 175.f, 375.f ;
% 	float Surface_pressure(time) ;
% 		Surface_pressure:orig_name = "Surface pressure" ;
% 		Surface_pressure:units = "hectoPascal" ;
% 		Surface_pressure:format = "F18.9" ;
% 		Surface_pressure:_FillValue = 3.402823e+38f ;
% 		Surface_pressure:valid_range = 0.f, 1100.f ;
% 	float Column_averaged_relative_humidity(time) ;
% 		Column_averaged_relative_humidity:orig_name = "Column averaged relative humidity" ;
% 		Column_averaged_relative_humidity:units = "N/A" ;
% 		Column_averaged_relative_humidity:format = "F18.9" ;
% 		Column_averaged_relative_humidity:_FillValue = 3.402823e+38f ;
% 		Column_averaged_relative_humidity:valid_range = 0.f, 100.f ;
% 	float Surface_minus_750_mb_air_temperature_difference(time) ;
% 		Surface_minus_750_mb_air_temperature_difference:orig_name = "Surface minus 750 mb air temperature difference" ;
% 		Surface_minus_750_mb_air_temperature_difference:units = "Kelvin" ;
% 		Surface_minus_750_mb_air_temperature_difference:format = "F18.9" ;
% 		Surface_minus_750_mb_air_temperature_difference:_FillValue = 3.402823e+38f ;
% 		Surface_minus_750_mb_air_temperature_difference:valid_range = -200.f, 200.f ;
% 	float Estimated_Inversion_Strength(time) ;
% 		Estimated_Inversion_Strength:orig_name = "Estimated Inversion Strength" ;
% 		Estimated_Inversion_Strength:units = "Kelvin" ;
% 		Estimated_Inversion_Strength:format = "F18.9" ;
% 		Estimated_Inversion_Strength:_FillValue = 3.402823e+38f ;
% 		Estimated_Inversion_Strength:valid_range = -200.f, 200.f ;
% 	float \750_mb_minus_surface_air_potential_temperature_difference(time) ;
% 		\750_mb_minus_surface_air_potential_temperature_difference:orig_name = "750 mb minus surface air potential temperature difference" ;
% 		\750_mb_minus_surface_air_potential_temperature_difference:units = "Kelvin" ;
% 		\750_mb_minus_surface_air_potential_temperature_difference:format = "F18.9" ;
% 		\750_mb_minus_surface_air_potential_temperature_difference:_FillValue = 3.402823e+38f ;
% 		\750_mb_minus_surface_air_potential_temperature_difference:valid_range = -200.f, 200.f ;
% 	float Precipitable_water(time) ;
% 		Precipitable_water:orig_name = "Precipitable water" ;
% 		Precipitable_water:units = "centimeters" ;
% 		Precipitable_water:format = "F18.9" ;
% 		Precipitable_water:_FillValue = 3.402823e+38f ;
% 		Precipitable_water:valid_range = 0.001f, 10.f ;
% 	short Cloud_mask_clear_strong_percent_coverage(time) ;
% 		Cloud_mask_clear_strong_percent_coverage:orig_name = "Cloud-mask clear-strong percent coverage" ;
% 		Cloud_mask_clear_strong_percent_coverage:units = "N/A" ;
% 		Cloud_mask_clear_strong_percent_coverage:format = "I10" ;
% 		Cloud_mask_clear_strong_percent_coverage:_FillValue = 32767s ;
% 		Cloud_mask_clear_strong_percent_coverage:valid_range = 0s, 100s ;
% 	short Cloud_mask_clear_weak_percent_coverage(time) ;
% 		Cloud_mask_clear_weak_percent_coverage:orig_name = "Cloud-mask clear-weak percent coverage" ;
% 		Cloud_mask_clear_weak_percent_coverage:units = "N/A" ;
% 		Cloud_mask_clear_weak_percent_coverage:format = "I10" ;
% 		Cloud_mask_clear_weak_percent_coverage:_FillValue = 32767s ;
% 		Cloud_mask_clear_weak_percent_coverage:valid_range = 0s, 100s ;
% 	float CWG_surface_skin_temperature(time) ;
% 		CWG_surface_skin_temperature:orig_name = "CWG surface skin temperature" ;
% 		CWG_surface_skin_temperature:units = "Kelvin" ;
% 		CWG_surface_skin_temperature:format = "F18.9" ;
% 		CWG_surface_skin_temperature:_FillValue = 3.402823e+38f ;
% 		CWG_surface_skin_temperature:valid_range = 175.f, 375.f ;
% 	float Clear_layer_overlap_percent_coverages(time, Conditions_clear__lower__upper__upper_over_lower) ;
% 		Clear_layer_overlap_percent_coverages:orig_name = "Clear/layer/overlap percent coverages" ;
% 		Clear_layer_overlap_percent_coverages:units = "N/A" ;
% 		Clear_layer_overlap_percent_coverages:format = "F18.9" ;
% 		Clear_layer_overlap_percent_coverages:_FillValue = 3.402823e+38f ;
% 		Clear_layer_overlap_percent_coverages:valid_range = 0.f, 100.f ;
% 	float Mean_visible_optical_depth_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_visible_optical_depth_for_cloud_layer:orig_name = "Mean visible optical depth for cloud layer" ;
% 		Mean_visible_optical_depth_for_cloud_layer:units = "N/A" ;
% 		Mean_visible_optical_depth_for_cloud_layer:format = "F18.9" ;
% 		Mean_visible_optical_depth_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_visible_optical_depth_for_cloud_layer:valid_range = 0.f, 512.f ;
% 	float Stddev_of_visible_optical_depth_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_visible_optical_depth_for_cloud_layer:orig_name = "Stddev of visible optical depth for cloud layer" ;
% 		Stddev_of_visible_optical_depth_for_cloud_layer:units = "N/A" ;
% 		Stddev_of_visible_optical_depth_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_visible_optical_depth_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_visible_optical_depth_for_cloud_layer:valid_range = 0.f, 300.f ;
% 	float Mean_logarithm_of_visible_optical_depth_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer:orig_name = "Mean logarithm of visible optical depth for cloud layer" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer:units = "N/A" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer:format = "F18.9" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer:valid_range = -6.f, 6.3f ;
% 	float Stddev_of_logarithm_of_visible_optical_depth_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_logarithm_of_visible_optical_depth_for_cloud_layer:orig_name = "Stddev of logarithm of visible optical depth for cloud layer" ;
% 		Stddev_of_logarithm_of_visible_optical_depth_for_cloud_layer:units = "N/A" ;
% 		Stddev_of_logarithm_of_visible_optical_depth_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_logarithm_of_visible_optical_depth_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_logarithm_of_visible_optical_depth_for_cloud_layer:valid_range = 0.f, 6.f ;
% 	float Mean_cloud_infrared_emissivity_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_infrared_emissivity_for_cloud_layer:orig_name = "Mean cloud infrared emissivity for cloud layer" ;
% 		Mean_cloud_infrared_emissivity_for_cloud_layer:units = "N/A" ;
% 		Mean_cloud_infrared_emissivity_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_infrared_emissivity_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_infrared_emissivity_for_cloud_layer:valid_range = 0.f, 2.f ;
% 	float Stddev_of_cloud_infrared_emissivity_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_cloud_infrared_emissivity_for_cloud_layer:orig_name = "Stddev of cloud infrared emissivity for cloud layer" ;
% 		Stddev_of_cloud_infrared_emissivity_for_cloud_layer:units = "N/A" ;
% 		Stddev_of_cloud_infrared_emissivity_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_cloud_infrared_emissivity_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_cloud_infrared_emissivity_for_cloud_layer:valid_range = 0.f, 2.f ;
% 	float Mean_cloud_top_pressure_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_top_pressure_for_cloud_layer:orig_name = "Mean cloud top pressure for cloud layer" ;
% 		Mean_cloud_top_pressure_for_cloud_layer:units = "hectoPascals" ;
% 		Mean_cloud_top_pressure_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_top_pressure_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_top_pressure_for_cloud_layer:valid_range = 0.f, 1100.f ;
% 	float Stddev_of_cloud_top_pressure_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_cloud_top_pressure_for_cloud_layer:orig_name = "Stddev of cloud top pressure for cloud layer" ;
% 		Stddev_of_cloud_top_pressure_for_cloud_layer:units = "hectoPascals" ;
% 		Stddev_of_cloud_top_pressure_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_cloud_top_pressure_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_cloud_top_pressure_for_cloud_layer:valid_range = 0.f, 600.f ;
% 	float Mean_cloud_top_temperature_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_top_temperature_for_cloud_layer:orig_name = "Mean cloud top temperature for cloud layer" ;
% 		Mean_cloud_top_temperature_for_cloud_layer:units = "Kelvin" ;
% 		Mean_cloud_top_temperature_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_top_temperature_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_top_temperature_for_cloud_layer:valid_range = 100.f, 350.f ;
% 	float Mean_cloud_top_height_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_top_height_for_cloud_layer:orig_name = "Mean cloud top height for cloud layer" ;
% 		Mean_cloud_top_height_for_cloud_layer:units = "kilometers" ;
% 		Mean_cloud_top_height_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_top_height_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_top_height_for_cloud_layer:valid_range = 0.f, 20.f ;
% 	float Mean_cloud_effective_pressure_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_effective_pressure_for_cloud_layer:orig_name = "Mean cloud effective pressure for cloud layer" ;
% 		Mean_cloud_effective_pressure_for_cloud_layer:units = "hectoPascals" ;
% 		Mean_cloud_effective_pressure_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_effective_pressure_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_effective_pressure_for_cloud_layer:valid_range = 0.f, 1100.f ;
% 	float Stddev_of_cloud_effective_pressure_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_cloud_effective_pressure_for_cloud_layer:orig_name = "Stddev of cloud effective pressure for cloud layer" ;
% 		Stddev_of_cloud_effective_pressure_for_cloud_layer:units = "hectoPascals" ;
% 		Stddev_of_cloud_effective_pressure_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_cloud_effective_pressure_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_cloud_effective_pressure_for_cloud_layer:valid_range = 0.f, 500.f ;
% 	float Mean_cloud_effective_temperature_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_effective_temperature_for_cloud_layer:orig_name = "Mean cloud effective temperature for cloud layer" ;
% 		Mean_cloud_effective_temperature_for_cloud_layer:units = "Kelvin" ;
% 		Mean_cloud_effective_temperature_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_effective_temperature_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_effective_temperature_for_cloud_layer:valid_range = 100.f, 350.f ;
% 	float Stddev_of_cloud_effective_temperature_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_cloud_effective_temperature_for_cloud_layer:orig_name = "Stddev of cloud effective temperature for cloud layer" ;
% 		Stddev_of_cloud_effective_temperature_for_cloud_layer:units = "Kelvin" ;
% 		Stddev_of_cloud_effective_temperature_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_cloud_effective_temperature_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_cloud_effective_temperature_for_cloud_layer:valid_range = 0.f, 150.f ;
% 	float Mean_cloud_effective_height_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_effective_height_for_cloud_layer:orig_name = "Mean cloud effective height for cloud layer" ;
% 		Mean_cloud_effective_height_for_cloud_layer:units = "kilometers" ;
% 		Mean_cloud_effective_height_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_effective_height_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_effective_height_for_cloud_layer:valid_range = 0.f, 20.f ;
% 	float Stddev_of_cloud_effective_height_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_cloud_effective_height_for_cloud_layer:orig_name = "Stddev of cloud effective height for cloud layer" ;
% 		Stddev_of_cloud_effective_height_for_cloud_layer:units = "kilometers" ;
% 		Stddev_of_cloud_effective_height_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_cloud_effective_height_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_cloud_effective_height_for_cloud_layer:valid_range = 0.f, 12.f ;
% 	float Mean_cloud_base_pressure_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_base_pressure_for_cloud_layer:orig_name = "Mean cloud base pressure for cloud layer" ;
% 		Mean_cloud_base_pressure_for_cloud_layer:units = "hectoPascals" ;
% 		Mean_cloud_base_pressure_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_base_pressure_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_base_pressure_for_cloud_layer:valid_range = 0.f, 1100.f ;
% 	float Stddev_of_cloud_base_pressure_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_cloud_base_pressure_for_cloud_layer:orig_name = "Stddev of cloud base pressure for cloud layer" ;
% 		Stddev_of_cloud_base_pressure_for_cloud_layer:units = "hectoPascals" ;
% 		Stddev_of_cloud_base_pressure_for_cloud_layer:format = "F18.9" ;
% 		Stddev_of_cloud_base_pressure_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Stddev_of_cloud_base_pressure_for_cloud_layer:valid_range = 0.f, 600.f ;
% 	float Mean_cloud_base_temperature_for_cloud_layer(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_base_temperature_for_cloud_layer:orig_name = "Mean cloud base temperature for cloud layer" ;
% 		Mean_cloud_base_temperature_for_cloud_layer:units = "Kelvin" ;
% 		Mean_cloud_base_temperature_for_cloud_layer:format = "F18.9" ;
% 		Mean_cloud_base_temperature_for_cloud_layer:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_base_temperature_for_cloud_layer:valid_range = 100.f, 350.f ;
% 	float Mean_liquid_water_path_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_liquid_water_path_for_cloud_layer__3_7_:orig_name = "Mean liquid water path for cloud layer (3.7)" ;
% 		Mean_liquid_water_path_for_cloud_layer__3_7_:units = "grams per square meter" ;
% 		Mean_liquid_water_path_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Mean_liquid_water_path_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Mean_liquid_water_path_for_cloud_layer__3_7_:valid_range = 0.f, 10000.f ;
% 	float Stddev_of_liquid_water_path_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_liquid_water_path_for_cloud_layer__3_7_:orig_name = "Stddev of liquid water path for cloud layer (3.7)" ;
% 		Stddev_of_liquid_water_path_for_cloud_layer__3_7_:units = "grams per square meter" ;
% 		Stddev_of_liquid_water_path_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Stddev_of_liquid_water_path_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Stddev_of_liquid_water_path_for_cloud_layer__3_7_:valid_range = 0.f, 8000.f ;
% 	float Mean_ice_water_path_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_ice_water_path_for_cloud_layer__3_7_:orig_name = "Mean ice water path for cloud layer (3.7)" ;
% 		Mean_ice_water_path_for_cloud_layer__3_7_:units = "grams per square meter" ;
% 		Mean_ice_water_path_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Mean_ice_water_path_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Mean_ice_water_path_for_cloud_layer__3_7_:valid_range = 0.f, 10000.f ;
% 	float Stddev_of_ice_water_path_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_ice_water_path_for_cloud_layer__3_7_:orig_name = "Stddev of ice water path for cloud layer (3.7)" ;
% 		Stddev_of_ice_water_path_for_cloud_layer__3_7_:units = "grams per square meter" ;
% 		Stddev_of_ice_water_path_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Stddev_of_ice_water_path_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Stddev_of_ice_water_path_for_cloud_layer__3_7_:valid_range = 0.f, 8000.f ;
% 	float Mean_water_particle_radius_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_water_particle_radius_for_cloud_layer__3_7_:orig_name = "Mean water particle radius for cloud layer (3.7)" ;
% 		Mean_water_particle_radius_for_cloud_layer__3_7_:units = "microns" ;
% 		Mean_water_particle_radius_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Mean_water_particle_radius_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Mean_water_particle_radius_for_cloud_layer__3_7_:valid_range = 0.f, 40.f ;
% 	float Stddev_of_water_particle_radius_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_water_particle_radius_for_cloud_layer__3_7_:orig_name = "Stddev of water particle radius for cloud layer (3.7)" ;
% 		Stddev_of_water_particle_radius_for_cloud_layer__3_7_:units = "microns" ;
% 		Stddev_of_water_particle_radius_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Stddev_of_water_particle_radius_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Stddev_of_water_particle_radius_for_cloud_layer__3_7_:valid_range = 0.f, 20.f ;
% 	float Mean_ice_particle_effective_radius_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__3_7_:orig_name = "Mean ice particle effective radius for cloud layer (3.7)" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__3_7_:units = "microns" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__3_7_:valid_range = 0.f, 300.f ;
% 	float Stddev_of_ice_particle_effective_radius_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Stddev_of_ice_particle_effective_radius_for_cloud_layer__3_7_:orig_name = "Stddev of ice particle effective radius for cloud layer (3.7)" ;
% 		Stddev_of_ice_particle_effective_radius_for_cloud_layer__3_7_:units = "microns" ;
% 		Stddev_of_ice_particle_effective_radius_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Stddev_of_ice_particle_effective_radius_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Stddev_of_ice_particle_effective_radius_for_cloud_layer__3_7_:valid_range = 0.f, 200.f ;
% 	float Mean_cloud_particle_phase_for_cloud_layer__3_7_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_cloud_particle_phase_for_cloud_layer__3_7_:orig_name = "Mean cloud particle phase for cloud layer (3.7)" ;
% 		Mean_cloud_particle_phase_for_cloud_layer__3_7_:units = "N/A" ;
% 		Mean_cloud_particle_phase_for_cloud_layer__3_7_:format = "F18.9" ;
% 		Mean_cloud_particle_phase_for_cloud_layer__3_7_:_FillValue = 3.402823e+38f ;
% 		Mean_cloud_particle_phase_for_cloud_layer__3_7_:valid_range = 1.f, 2.f ;
% 	float Mean_water_particle_radius_for_cloud_layer__1_2_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_water_particle_radius_for_cloud_layer__1_2_:orig_name = "Mean water particle radius for cloud layer (1.2)" ;
% 		Mean_water_particle_radius_for_cloud_layer__1_2_:units = "microns" ;
% 		Mean_water_particle_radius_for_cloud_layer__1_2_:format = "F18.9" ;
% 		Mean_water_particle_radius_for_cloud_layer__1_2_:_FillValue = 3.402823e+38f ;
% 		Mean_water_particle_radius_for_cloud_layer__1_2_:valid_range = 0.f, 40.f ;
% 	float Mean_ice_particle_effective_radius_for_cloud_layer__1_2_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__1_2_:orig_name = "Mean ice particle effective radius for cloud layer (1.2)" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__1_2_:units = "microns" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__1_2_:format = "F18.9" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__1_2_:_FillValue = 3.402823e+38f ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__1_2_:valid_range = 0.f, 300.f ;
% 	float Mean_logarithm_of_visible_optical_depth_for_cloud_layer__1_2_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__1_2_:orig_name = "Mean logarithm of visible optical depth for cloud layer (1.2)" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__1_2_:units = "N/A" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__1_2_:format = "F18.9" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__1_2_:_FillValue = 3.402823e+38f ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__1_2_:valid_range = -6.f, 6.3f ;
% 	float Mean_water_particle_radius_for_cloud_layer__2_1_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_water_particle_radius_for_cloud_layer__2_1_:orig_name = "Mean water particle radius for cloud layer (2.1)" ;
% 		Mean_water_particle_radius_for_cloud_layer__2_1_:units = "microns" ;
% 		Mean_water_particle_radius_for_cloud_layer__2_1_:format = "F18.9" ;
% 		Mean_water_particle_radius_for_cloud_layer__2_1_:_FillValue = 3.402823e+38f ;
% 		Mean_water_particle_radius_for_cloud_layer__2_1_:valid_range = 0.f, 40.f ;
% 	float Mean_ice_particle_effective_radius_for_cloud_layer__2_1_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__2_1_:orig_name = "Mean ice particle effective radius for cloud layer (2.1)" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__2_1_:units = "microns" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__2_1_:format = "F18.9" ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__2_1_:_FillValue = 3.402823e+38f ;
% 		Mean_ice_particle_effective_radius_for_cloud_layer__2_1_:valid_range = 0.f, 225.f ;
% 	float Mean_logarithm_of_visible_optical_depth_for_cloud_layer__2_1_(time, Lower_and_Upper_Cloud_Layers) ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__2_1_:orig_name = "Mean logarithm of visible optical depth for cloud layer (2.1)" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__2_1_:units = "N/A" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__2_1_:format = "F18.9" ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__2_1_:_FillValue = 3.402823e+38f ;
% 		Mean_logarithm_of_visible_optical_depth_for_cloud_layer__2_1_:valid_range = -6.f, 6.3f ;
% 	short Percentage_of_CERES_FOV_with_MODIS_land_aerosol(time) ;
% 		Percentage_of_CERES_FOV_with_MODIS_land_aerosol:orig_name = "Percentage of CERES FOV with MODIS land aerosol" ;
% 		Percentage_of_CERES_FOV_with_MODIS_land_aerosol:units = "N/A" ;
% 		Percentage_of_CERES_FOV_with_MODIS_land_aerosol:format = "I10" ;
% 		Percentage_of_CERES_FOV_with_MODIS_land_aerosol:_FillValue = 32767s ;
% 		Percentage_of_CERES_FOV_with_MODIS_land_aerosol:valid_range = 0s, 100s ;
% 	short PSF_wtd_MOD04_cloud_fraction_land(time) ;
% 		PSF_wtd_MOD04_cloud_fraction_land:orig_name = "PSF-wtd MOD04 cloud fraction land" ;
% 		PSF_wtd_MOD04_cloud_fraction_land:units = "N/A" ;
% 		PSF_wtd_MOD04_cloud_fraction_land:format = "I10" ;
% 		PSF_wtd_MOD04_cloud_fraction_land:_FillValue = 32767s ;
% 		PSF_wtd_MOD04_cloud_fraction_land:valid_range = 0s, 100s ;
% 	int PSF_wtd_MOD04_aerosol_types_land(time) ;
% 		PSF_wtd_MOD04_aerosol_types_land:orig_name = "PSF-wtd MOD04 aerosol types land" ;
% 		PSF_wtd_MOD04_aerosol_types_land:units = "N/A" ;
% 		PSF_wtd_MOD04_aerosol_types_land:format = "I10" ;
% 		PSF_wtd_MOD04_aerosol_types_land:_FillValue = 2147483647 ;
% 		PSF_wtd_MOD04_aerosol_types_land:valid_range = 0, 9999 ;
% 	float PSF_wtd_MOD04_corrected_optical_depth_land__0_550_(time) ;
% 		PSF_wtd_MOD04_corrected_optical_depth_land__0_550_:orig_name = "PSF-wtd MOD04 corrected optical depth land (0.550)" ;
% 		PSF_wtd_MOD04_corrected_optical_depth_land__0_550_:units = "N/A" ;
% 		PSF_wtd_MOD04_corrected_optical_depth_land__0_550_:format = "F18.9" ;
% 		PSF_wtd_MOD04_corrected_optical_depth_land__0_550_:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_corrected_optical_depth_land__0_550_:valid_range = 0.f, 5.f ;
% 	short Percentage_of_CERES_FOV_with_MODIS_deep_blue_aerosol(time) ;
% 		Percentage_of_CERES_FOV_with_MODIS_deep_blue_aerosol:orig_name = "Percentage of CERES FOV with MODIS deep blue aerosol" ;
% 		Percentage_of_CERES_FOV_with_MODIS_deep_blue_aerosol:units = "N/A" ;
% 		Percentage_of_CERES_FOV_with_MODIS_deep_blue_aerosol:format = "I10" ;
% 		Percentage_of_CERES_FOV_with_MODIS_deep_blue_aerosol:_FillValue = 32767s ;
% 		Percentage_of_CERES_FOV_with_MODIS_deep_blue_aerosol:valid_range = 0s, 100s ;
% 	float PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_(time) ;
% 		PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_:orig_name = "PSF-wtd MOD04 deep blue aerosol optical depth land (0.550)" ;
% 		PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_:units = "N/A" ;
% 		PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_:format = "F18.9" ;
% 		PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_:valid_range = 0.f, 5.f ;
% 	float PSF_wtd_MOD04_deep_blue_angstrom_exponent_land(time) ;
% 		PSF_wtd_MOD04_deep_blue_angstrom_exponent_land:orig_name = "PSF-wtd MOD04 deep blue angstrom exponent land" ;
% 		PSF_wtd_MOD04_deep_blue_angstrom_exponent_land:units = "N/A" ;
% 		PSF_wtd_MOD04_deep_blue_angstrom_exponent_land:format = "F18.9" ;
% 		PSF_wtd_MOD04_deep_blue_angstrom_exponent_land:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_deep_blue_angstrom_exponent_land:valid_range = 0.5f, 5.f ;
% 	float PSF_wtd_MOD04_deep_blue_single_scattering_albedo_land__0_470_(time) ;
% 		PSF_wtd_MOD04_deep_blue_single_scattering_albedo_land__0_470_:orig_name = "PSF-wtd MOD04 deep blue single scattering albedo land (0.470)" ;
% 		PSF_wtd_MOD04_deep_blue_single_scattering_albedo_land__0_470_:units = "N/A" ;
% 		PSF_wtd_MOD04_deep_blue_single_scattering_albedo_land__0_470_:format = "F18.9" ;
% 		PSF_wtd_MOD04_deep_blue_single_scattering_albedo_land__0_470_:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_deep_blue_single_scattering_albedo_land__0_470_:valid_range = 0.f, 1.f ;
% 	short Percentage_of_CERES_FOV_with_MODIS_ocean_aerosol(time) ;
% 		Percentage_of_CERES_FOV_with_MODIS_ocean_aerosol:orig_name = "Percentage of CERES FOV with MODIS ocean aerosol" ;
% 		Percentage_of_CERES_FOV_with_MODIS_ocean_aerosol:units = "N/A" ;
% 		Percentage_of_CERES_FOV_with_MODIS_ocean_aerosol:format = "I10" ;
% 		Percentage_of_CERES_FOV_with_MODIS_ocean_aerosol:_FillValue = 32767s ;
% 		Percentage_of_CERES_FOV_with_MODIS_ocean_aerosol:valid_range = 0s, 100s ;
% 	short PSF_wtd_MOD04_cloud_fraction_ocean(time) ;
% 		PSF_wtd_MOD04_cloud_fraction_ocean:orig_name = "PSF-wtd MOD04 cloud fraction ocean" ;
% 		PSF_wtd_MOD04_cloud_fraction_ocean:units = "N/A" ;
% 		PSF_wtd_MOD04_cloud_fraction_ocean:format = "I10" ;
% 		PSF_wtd_MOD04_cloud_fraction_ocean:_FillValue = 32767s ;
% 		PSF_wtd_MOD04_cloud_fraction_ocean:valid_range = 0s, 100s ;
% 	float PSF_wtd_MOD04_cloud_condensation_nuclei_ocean__average(time) ;
% 		PSF_wtd_MOD04_cloud_condensation_nuclei_ocean__average:orig_name = "PSF-wtd MOD04 cloud condensation nuclei ocean, average" ;
% 		PSF_wtd_MOD04_cloud_condensation_nuclei_ocean__average:units = "CCN per square centimeter" ;
% 		PSF_wtd_MOD04_cloud_condensation_nuclei_ocean__average:format = "F18.9" ;
% 		PSF_wtd_MOD04_cloud_condensation_nuclei_ocean__average:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_cloud_condensation_nuclei_ocean__average:valid_range = 0.f, 1.e+11f ;
% 	float PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_(time) ;
% 		PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_:orig_name = "PSF-wtd MOD04 effective optical depth average ocean (0.550)" ;
% 		PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_:units = "N/A" ;
% 		PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_:format = "F18.9" ;
% 		PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_:valid_range = 0.f, 5.f ;
% 	float PSF_wtd_MOD04_optical_depth_small_average_ocean__0_550_(time) ;
% 		PSF_wtd_MOD04_optical_depth_small_average_ocean__0_550_:orig_name = "PSF-wtd MOD04 optical depth small average ocean (0.550)" ;
% 		PSF_wtd_MOD04_optical_depth_small_average_ocean__0_550_:units = "N/A" ;
% 		PSF_wtd_MOD04_optical_depth_small_average_ocean__0_550_:format = "F18.9" ;
% 		PSF_wtd_MOD04_optical_depth_small_average_ocean__0_550_:_FillValue = 3.402823e+38f ;
% 		PSF_wtd_MOD04_optical_depth_small_average_ocean__0_550_:valid_range = 0.f, 5.f ;
%
% // global attributes:
% 		:Conventions = "CF-1.0" ;
% 		:Subsetter_title = "ASDC CERES Subset" ;
% 		:Subsetter_version = "2.8.b1" ;
% 		:Subsetter_institution = "Atmospheric Science Data Center (ASDC) http://eosweb.larc.nasa.gov" ;
% 		:Subsetter_history = "2015-12-18T08:48:15 -0500 SubsetCeresSsf" ;
% 		:Subsetter_temporalFilter = "2008-11-12T00:00:00.000000Z to 2008-11-14T23:59:59.000000Z" ;
% 		:Subsetter_spatialFilter = "POLYGON ((-160 -40, -160 10, -60 10, -60 -40, -160 -40))" ;
% 		:Subsetter_parameterFilter = "none" ;
% 		:history = "Fri Dec 18 08:59:20 2015: ncrcat -o CERES_SSF_XTRK-MODIS_Edition4A_Subset_2008111200-2008111422.nc" ;
% 		:nco_input_file_number = 33 ;
