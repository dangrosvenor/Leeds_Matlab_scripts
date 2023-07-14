try
    %If just plotting then change what to plot and loat/lon in
    %plot_global_maps
    
    %read and average multiple modis files from a directory
    %see MODIS_multi_dir_process for the list of the directories to process,
    %etc.
    %Is also used to make mock L3 data from L2 swaths - search for 'make mock L3 data'
    % ***   make_mockL3_variables is important - this sets the variables that are
    % saved!! ***
    % make_MODIS_variables_01_mockL3 -  this is where the data is added to the large array
    % Call to open_L2_MODIS_file_01, etc. made here. filtering_data_L2
    % called from open_L2_MODIS_file_01.


    % if ~exist('make_mockL3_variables_script_name')
    %     make_mockL3_variables_script_name='make_mockL3_variables';
    % end

    suppress_output=1;
    Nswath_dim=20; %size of the swath dimension (max number of swaths).
    %N.B. needs to be at least 18!! (even though should only be 16?)

    action2 = action; %for make_MODIS_variables_01_mockL3

%% Go through and check all of the files that we will actually use (that are actually hdf files)    
    skip_filecheck=0;
    if skip_filecheck==0

        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/'; %for the plots



        switch action
            case 'store data etc'
                %these are generic sensor angles and scantimes that should apply for all
                %swaths, I think
                load('/home/disk/eos1/d.grosvenor/matlab/work/MODIS.mat','scantime_L2','sensor_zenith_L2');

            case 'draw and save plots'
                %don't do anything

        end


        %choose type of MODIS file to read - moved this to MODIS_multi_dir_process
        %MODIS_filetype = 'L2'; %standard L2 file
        %MODIS_filetype = 'L2 Joint'; %Joint L2 file - subsampled dataset (5 km)




        switch time_restrict_case
            %set the times required for the time restriction (if required)
            case 'choose here'
                %the default case - just make the savename
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_data_' action '_' datestr(now,30)]; %now gives the data and time

                times_restrict = datenum(2010,2,6,14,15,0); %flight 99 14:15 overpass
                times_restrict = datenum(2010,2,9,14,45,0); %flight 101 14:45 overpass


            case 'VOCALS cloud segments (Chris T)'
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_VOCALS_matches_' datestr(now,30)]; %now gives the data and time

                %        isave_MODIS_data=0;
                VOCALS_cloud_leg_info  %(external script)
                %N = DATENUM(Y,MO,D,H,MI,S)  (in days since the reference time of 01-Jan 00:00)
                %list of the vocals cloud leg times in Matlab datenum time
                times_restrict = datenum(year_voc,month_voc,day_voc,hours_voc,mins_voc,secs_voc);
                timLON = lon_voc;

            case 'Puijo station'
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_15mins' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_2mins' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_15mins_5km' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_2mins_11km_assumeUTC' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_2mins_11km_UTC_minus_10mins' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_2mins_5km_single_liquid_layers_only' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_2mins_5km' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_30mins_11km_DMPS' datestr(now,30)]; %now gives the data and time
                savevarname = [savedir_var remove_character(modis_dir,'/','_') 'L2_Puijo_matches_30mins_0.5x1deg_DMPS' datestr(now,30)]; %now gives the data and time
                
                %        isave_MODIS_data=0;
                datatype_Puijo = 'CDP';
                datatype_Puijo = 'DMPS';
                switch datatype_Puijo
                    case 'CDP'
                        % -------------------------------------------------------------
                        assume_UTC = 'no';
                        %        assume_UTC = 'yes';  %assume that the CDP was in UTC (not in UTC+2 as suggested by Irshad) for a test
                        %        assume_UTC = 'UTC_minus_10mins';
                        Puijo_station_cloud_event_info  %(external script)
                        % -------------------------------------------------------------
                        %N = DATENUM(Y,MO,D,H,MI,S)  (in days since the reference time of 01-Jan 00:00)
                        %list of the vocals cloud leg times in Matlab datenum time
                        %        times_restrict = datenum(year_voc,month_voc,day_voc,hours_voc,mins_voc,secs_voc);
                        times_restrict = MatlabTime_Puijo;
                        %these are now in UTC, but have selected the mid-point of each hour
                        %(e.g. 01:30, 02:30, etc.), which would be averages from 01:00 to
                        %02:00, etc.
                        %        timLON = lon_voc;

                    case 'DMPS'
                        %read the aerosol data from Irshad

                        load('/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/DMPS_80_100_150nm.mat');
                        times_restrict = DMPS_80_100_150nm(:,1);
                        Nacc_DMPS = DMPS_80_100_150nm(:,9); %accumulation mode aerosol concentration (cm3).
                end

                ioverride_read_weather_Irshad = 1;
                draw_square=0;
                % -------------------------------------------------------------
                read_weather_data_Irshad %external script
                % -------------------------------------------------------------

        end


        files_mod = dir([filedir modis_dir]);

        %loop through files and store decide which ones we want to process
        imod_store=0;
        ireject=0;
        clear files_to_read time_match_cloud_legs imod_success rejected_swaths
        for imod_read=1:length(files_mod)
            if length( strfind(files_mod(imod_read).name,'.hdf') )>0
                success=1;
                switch MODIS_filetype
                    case 'L2 Joint'
                        %test to see if the file contains solar_zenith angle
                        [success]=L2_joint_test_file([filedir modis_dir files_mod(imod_read).name]);
                    case 'L2'
                        %test to see if the file contains any usable solar_zenith
                        %angles
                        [success]=L2_test_file([filedir modis_dir files_mod(imod_read).name]);
                    otherwise %add AMSRE_L2
                        success=1;
                end


                if success==1 %if the file contains solar_zenith - the nightime-only ones don't bother with SZA

                    if time_restrict==1

                        file_name_h5=files_mod(imod_read).name;


                        iday=findstr(file_name_h5,'.A');
                        modis_year_str=file_name_h5(iday+2:iday+5);
                        modis_day_str=file_name_h5(iday+6:iday+8);
                        modis_time_str=file_name_h5(iday+10:iday+13);
                        aq_terr_str = file_name_h5(iday-8:iday-1);
                        date_str=datestr(datenum(['01-Jan-' modis_year_str])+str2num(modis_day_str)-1,1);

                        times_swath = datenum([date_str ' ' modis_time_str(1:2) ':' modis_time_str(3:4)]);

                        %absolute difference between all the vocals cloud leg times in hours
                        dtime = 24*abs(times_swath - times_restrict);
                        ires = find(dtime<=dtime_restrict);

                        if length(ires)>1
                            [minval,iclosest] = min(dtime(ires));
                            %May have a few datapoints that match within tolerance
                            %(although a max of 2 that will be equidistant if we are doing +/- 30 mins)
                            %it's possible that there will be 2 CDP datapoints that are
                            %equidistant - i.e. if the MODIS swath is exactly halfway.
                            %Should really average them, but will just take one for now
                            %                 ires = find(abs(dtime(ires)-minval)<1e-10);
                            ires = ires(iclosest);

                        end

                        %if we have a time match then will load the file (and check for
                        %a location match)
                        if length(ires)>0

                            imod_store=imod_store+1;
                            files_to_read(imod_store)=imod_read;
                            time_match_cloud_legs{imod_store}=ires;
                            time_of_swath{imod_store}=times_swath;

                        else
                            ireject = ireject + 1;
                            rejected_swaths{ireject} = times_swath;

                        end
                    else
                        imod_store=imod_store+1;
                        files_to_read(imod_store)=imod_read;
                    end %time_restrict

                    imod_success(imod_read)=1;
                else
                    imod_success(imod_read)=0;
                end %success
            end
        end

    end

    % files_mod(1:2)=[]; %these are the . and .. listings
    %
    % files_to_read=[1:length(files_mod)-2];  %

    if exist('files_to_read')
        nMOD_av = length(files_to_read);
    else
        nMOD_av = 0;
    end


    if exist('SD_id')
        status = hdfsd('end',SD_id);
    end

    %put the required variable names into modis_var{i}
    set_MODIS_L2_variable_names

%% ______________________________________
    % **** now process the files ****
    % ______________________________________

    nday_old=0;
    first_plot=1;
    
    N_save =[]; lat_save=[]; lon_save=[];



    for imod_read=1:nMOD_av
        %    imod_read
        if suppress_output==0
            fprintf(1,'\n%d of %d ',imod_read,nMOD_av);
        end
        imr=files_to_read(imod_read);

        file_name_h5 = [modis_dir files_mod(imr).name];

        imodis_file_override=1;
        night_time=0;

%% read the L2 file
        switch MODIS_filetype
            case 'L2'
                success_tau=1;
                open_L2_MODIS_file_01
            case 'L2_C6'
                success_tau=1;
                open_L2_C6_MODIS_file_01  
                N = N37;
            case 'L2 Joint'
                success_tau=1;
                open_L2_Joint_MODIS_file_02
            case 'AMSRE_L2'
                success_tau=1;
                open_L2_AMSRE_file_01
        end

        nday=num2str(modis_day_str);

        switch action
            case 'make mock L3 data daily' %put into Nswath_dim swaths for each day instead of one
                %map for each swath
                switch multiL2L3_project
                    case {'Puijo_Sami','Puijo_Sami_jointL2','Puijo_Sami_DMPS_L2'}
                        %For the Puijo case where we are tyring to match a station value to a particular area
                        %of the swath will just do mock_L3_from_L2 for the region
                        %(based on wind drift between the swath time and cloud event
                        %time)
                        %Can set LAT_min and LON_min, etc. and LAT_step, LON_step to
                        %make a box (just one square) around the desired location.
                        %also store Nd from the CDP here -  Nd_Puijo(time_match_cloud_legs{imod_read})
                        %time difference between the swath and the cloud event
                        time_diff = time_of_swath{imod_read} - times_restrict(time_match_cloud_legs{imod_read});
                        %if positive, the swath happened after the cloud event.
                        %So then the area of cloud will have drifted WITH (not against)
                        %the wind direction

                        CDP_av_time = 1; %hours
                        CDP_av_time = 5/60; %hours

                        %find the corresponding weather entry
                        iweather = find(abs(times_restrict(time_match_cloud_legs{imod_read})  - ( weather(:,1)+ 0.5/24) ) < 0.01/24  );
                        % Should really average together the CDP values if we are equidistant, but
                        % leave for now
                        %                     for imulti_CDP = 1:length(time_match_cloud_legs{imod_read})
                        %                         iweather2(imulti_CDP) = find(abs(times_restrict(time_match_cloud_legs{imod_read}(imulti_CDP)) - 0.5/24 - weather(:,1)) < 0.01/24 );
                        %                     end

                        if length(iweather)==0
                            iweather=NaN;
                        end

                        if use_wind_Irshad==1
                            if length(iweather==1)
                                wind_speed = weather(iweather,13);
                                wind_dir = weather(iweather,19);
                            else %if no weather data
                                continue
                            end
                            if isnan(wind_speed) | isnan(wind_dir)
                                continue
                            end


                            %total distance moved in metres
                            dt = time_diff*24*3600; %
                            [u,v]=uv_from_winddir(wind_speed,wind_dir); %wind_speed in m/s, wind_dir in degrees
                            %clockwise from north
                            dx = u.*dt;
                            dy = v.*dt;


                            [lon_new,lat_new]=calc_lat_lon_change_for_dx_dy(27.656,62.909,dx,dy);

                        else
                            lon_new = 27.656;
                            lat_new = 62.909;


                        end

                        switch Puijo_box_size
                            case '0.5 x 1 deg (lat x lon)'

                                %now just set a box centred around this location that
                                %is the right size for our averaging - we want an
                                %average for an approx an hour - set to be around
                                %10*3600/1e3 = 36x36km (assuming approx wind speed of
                                %10 m/s). At this latitude 1 deg east is 50.7 km and
                                %north is 111.3 km. So could do one degree east-west
                                %and 0.5 degrees north-south

                                override_mockL3_options=1;
                                box_type = 'regular lat lon grid';

                                LAT_min = lat_new-0.25; LAT_max = lat_new+0.25; LAT_step = 0.5;
                                LON_min = lon_new-0.5; LON_max = lon_new+0.5; LON_step = 1;

                            case 'NxNkm'
                                override_mockL3_options=1;
                                box_type = 'Selected NxN pixel region';
                                Npix=Npix_Puijo; %makes squares of Npix*Npix


                                %                    LAT_min = lat_new-0.25; LAT_max = lat_new+0.25; LAT_step = 0.5;
                                %                    LON_min = lon_new-0.5; LON_max = lon_new+0.5; LON_step = 1;

                        end


                        action2 = 'make mock L3 all days';

                        Nswath_dim = nMOD_av;
                        %here are using the iswath dimension instead for the
                        %days - are storing all datapoints
                        
                    case 'Iceland_30W30E_40N90N'    
                        override_mockL3_options=1;
%                        box_type = 'Selected NxN pixel region';
%                        box_type = 'NxN pixel square';
%                        Npix=Npix_Iceland; %makes squares of Npix*Npix
% Couldn't get the above to work - but perhaps better on a regular lat long
% grid anyway?
                        
                         box_type = 'regular lat lon grid';
                         LAT_step=0.25;
                         LON_step=0.25;
                         
                         LAT_step=1.0; %Will stick to this for PDFs etc. first for speed
                         LON_step=1.0;                         

                    otherwise
                        %set these flags to override the defaults in mock_L3_from_L2
                        override_mockL3_options=1;
                        box_type = 'regular lat lon grid';
                        LAT_step=1;
                        LON_step=1;

                end



                %on first pass define memory for all of the timeseries3 variables
                if imod_read==1
                    save_or_load = 'save';

                    %--------------------------------------------------------------------------
                    make_mockL3_variables; %just puts the variables we want to extract
                    %into modis_var array
                    %--------------------------------------------------------------------------

                    %now call mock_L3_from_L2 - but this is just a call to get the size of the arrays that
                    %will be produced
                    ijust_size_check=1; %flag to tell it this
                    if exist('irun_from_modis_multi_DAY_processL3L2') & irun_from_modis_multi_DAY_processL3L2==1
                        override_mockL3_options2=1;
                        override_mockL3_options=1;
                    end
                    mock_L3_from_L2
                    for istring=1:length(modis_var)
                        %if we process one day at a time then we know we only
                        %have one day of data!
                        if strcmp(modis_var{istring},'MODIS_swath_filename')==1
                            %make an array of empty cells
                            eval_str = [modis_var{istring} '.timeseries3 = cell(length(LATS)-1,length(LONS)-1,' num2str(Nswath_dim) ');'];
                        else
                            eval_str = [modis_var{istring} '.timeseries3 = NaN*ones(length(LATS)-1,length(LONS)-1,' num2str(Nswath_dim) ');'];
                        end
                        eval(eval_str);
                    end
                    %number of swaths written to a particular lat lon - start
                    %at zero
                    nswath = zeros(length(LATS)-1,length(LONS)-1);
                end

                if success_tau==1 %some files had no tau??

                    %average the data into 1x1 degree regions
                    if exist('irun_from_modis_multi_DAY_processL3L2') & irun_from_modis_multi_DAY_processL3L2==1
                        override_mockL3_options=1;
                        override_mockL3_options2=1;
                    end
                    %---------------------------------------
                    mock_L3_from_L2
                    %---------------------------------------
                    %rename the variables to be consistent with timeseries3 data
                    iswath=imod_read;
                    %this is where the data is added to the large array
                    make_MODIS_variables_01_mockL3

                end

            case 'make mock L3 data'
                %set these flags to override the defaults in mock_L3_from_L2
                override_mockL3_options=1;
                box_type = 'regular lat lon grid';
                LAT_step=1;
                LON_step=1;

                %on first pass define memory for all of the timeseries3 variables
                if imod_read==1
                    save_or_load = 'save';

                    %--------------------------------------------------------------------------
                    make_mockL3_variables; %just puts the variables we want to extract
                    %in modis_var array
                    %--------------------------------------------------------------------------

                    %now call mock_L3_from_L2 - but this is just a call to get the size of the arrays that
                    %will be produced
                    ijust_size_check=1; %flag to tell it this
                    mock_L3_from_L2
                    for istring=1:length(modis_var)
                        eval_str = [modis_var{istring} '.timeseries3 = NaN*ones(length(LATS)-1,length(LONS)-1,nMOD_av);'];
                        eval(eval_str);
                    end
                end


                %average the data into 1x1 degree regions
                mock_L3_from_L2

                %rename the variables to be consistent with timeseries3 data
                iswath=imod_read;
                %this is where the data is added to the large array
                make_MODIS_variables_01_mockL3

            case 'store data etc'
                MODIS_average_L2_data

            case 'draw and save plots'
                limit_time=0;
                if limit_time==1
                   med_lon = prctile(lon(:),50); %Find the approx centre lon of the track
                   [Y,M,D,H,MIN]=datevec(modis_date_time); %UTC time of overpass
                   local_time = H+MIN/60 + med_lon/15;
                   
                   scantime_matlab; %
                    
                end
                
                reject_swath=0;
                
                reject_if_over_plotting=0;
                if reject_if_over_plotting==1
                    %Do not plot a swath if it will overplot onto stuff
                    %already plotted anywhere in the box of interest
                    
                    if first_plot==1
                        thresh_LAT = [48 80];  thresh_LON = [-70 40]; %The box being considered
                        lats = thresh_LAT(1):0.1:thresh_LAT(2);
                        lons = thresh_LON(1):0.1:thresh_LON(2);
                        [lon_box,lat_box]=meshgrid(lons,lats);
                        %Create the array saying whether points have been
                        %plotted on or not.
                        points_written = zeros(size(lon_box));
                        %Sub-set the points actually within our box
                        %ilat_box = find(lat>=thresh_LAT(1) & lat<=thresh_LAT(2) & lon>=thresh_LON(1) & lon<=thresh_LON(2));
                        
                        scantime_matlab_old = scantime_matlab;                        
                    end
                    
                    if abs(scantime_matlab - scantime_matlab_old) < 6/60/24 %less than 6 mins apart
                        concurrent_swaths=1; %N.B. - this relies on the filenames from dir command beign in time order - but seems to be 
                        %the case, so hopefully will be ok.
                    else
                        concurrent_swaths=0;
                    end

                    %Interpolate some data onto the regular grid of our
                    %domain - should make it NaN outside of the swath
                    dat_grid = griddata(lon,lat,scantime,lon_box,lat_box);
                    inan=isnan(dat_grid);
                    dat_grid(inan)=0;
                    dat_grid(~inan)=1;
                    sum_points = points_written + dat_grid;
                    

                    if maxALL(sum_points)==2 & concurrent_swaths==0  %If have overplotted onto previous plots
                                              %Do not update points_written
                                              %since we are rejecting this
                                              %swath. But plot if they are
                       %concurrent plots since otherwise they don't get
                       %plotted due to the relative coarseness of the
                       %lon_box, lat_box grid.
                        reject_swath=1;
                    else
                        %Otherwise accepth the swath and make the combined
                        %no. points the new total
                        points_written = sum_points; 
                        %Reduce back to one since may get some overlap in
                        %cases of concurrent swaths
                        points_written(points_written>1)=1;
                    end
                    
                    scantime_matlab_old = scantime_matlab;
                    
                end
                
                if night_time==0
                    supress_colorbar=0; %Setting to one stops the colorbar from appearing
                    
                    if composite_single_days==1 & (nday_old==nday) %so, not a new day - keep plotting on the old figure
                        inew_figure=0; %add to the old figure
                        supress_colorbar=1; %Setting to one stops the colorbar from appearing
                        
                    elseif first_plot==1 %except on first pass we want a new figure
                        inew_figure=1;
                    else %Started a new day, so save and close the old figure, and open a new one.
                        inew_figure=1; %save the figure and open a new one
                        %                    caxis(climits);
                        
                        increase_font_size_map_figures
                        
                        %Save without the PDF/eps versin since this is slow for
                        %MODIS L2 plots.
                        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'','/home/disk/eos8/d.grosvenor/misc_files/default_fig_notes',0);
                        if iclose_figures==1
                            close(gcf);
                        end
                      

                    end



                    nday_old=nday;

                    %the plot gets made here - set up the desired plot in plot_global_maps
                    if night_time==0 & reject_swath==0
                        
                        isave_all_data_then_plot=0;
                        if isave_all_data_then_plot==1
                            N_save = cat(1,N_save,N);
                            lat_save = cat(1,lat_save,Plat2_L2);
                            lon_save = cat(1,lon_save,Plon2_L2);
%                            lat2_edges_save = cat(1,lat2_save,Plat2D);
%                            lon2_edges_save = cat(1,lon2_save,Plon2D);                            
                            if imod_read==nMOD_av
                                N = N_save;
                                Plat2_L2 = lat_save;
                                Plon2_L2 = lon_save;
                                iover_ride_plot_global=1;
                                i_increase_font_size_map_figures_OFF=1; %don't do i_increase_font_size_map_figures_OFF
                                inew_figure=1;
                                plot_global_maps;
                            end
                        else
                            iover_ride_plot_global=1;
                            i_increase_font_size_map_figures_OFF=1; %don't do i_increase_font_size_map_figures_OFF
                            plot_global_maps;
                        end
                        first_plot=0;
                    end

                end






        end       %switch action
        %% end of switch action

        status = hdfsd('end',SD_id); %end access to the file - think otherwise things start to go wrong
        
       


    end  %for imod_read=1:nMOD_av
    %%

    if nMOD_av>0 %only save if had files to process
        switch action
            case 'draw and save plots'                
%                saveas_ps_fig_emf(gcf,[savename]);
                %Save without the PDF/eps versin since this is slow for
                %MODIS L2 plots.
                saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'','/home/disk/eos8/d.grosvenor/misc_files/default_fig_notes',0);
                if iclose_figures==1
                    close(gcf);
                end
                
            case 'store data etc'
                if isave_MODIS_data==1
                    switch time_restrict_case
                        case 'choose here'
                            save_MODIS_L2_data
                        case 'VOCALS cloud segments (Chris T)'
                            save(savevarname,'VOCALS_leg_number_matches');
                    end
                end

            case {'make mock L3 data','make mock L3 data daily'}
                timLAT = LATS;
                timLON = LONS;         %the regular grid used in mock_L3_from_L2
                save_MODIS_L2_data
        end
    end

    clear suppress_output
catch multiL2_error
    clear supresss_output
    rethrow(multiL2_error)
end