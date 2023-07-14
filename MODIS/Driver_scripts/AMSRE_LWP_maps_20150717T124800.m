% LWP map plots for AMSRE data
%The level of zoom is controled by thresh_LAT and LON and restrict_domain
%in plot_global
%Define LAT_val and LON_val elsewhere to plot a box highlighting a region
%Will also calculate the overpass time for this box

%Runs multi_read_amsre_daily

%Will setup up to loop over several days and satellites (AMSRE, SSMI,
%Windsat, TMI) and loop over ascending and descending to make plots of all

i_mask_low_LWP=1; %Make any values below thresh_LWP equal to NaN for the plotting
thresh_LWP_mask = 20;

isave_plot=0;
%savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/GOES/';
                            
sats={'AMSRE','SSMI-f13','SSMI-f15','SSMI-f16','SSMI-f17'};
sats={'SSMI-f13','SSMI-f15'};
sats={'AMSRE','SSMI-f13','SSMI-f15','SSMI-f16','SSMI-f17','TMI','Windsat'};
sats={'AMSRE'};
%sats={'TMI'};
%sats={'Windsat'};

isave_vars_driver = 1;
save_mat_file_remss = ['/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_' datestr(now,30) '.mat'];

days_lwp = [13]; %[12:15];
months_lwp = 11; iyear_lwp=1;
years_lwp = 2008; imonth_lwp=1;

idriver=0;

time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)

var_names={'lwp','time'};  %Leaving out sst here

%thresh_frac_covered = 0.99;
thresh_frac_covered = 0.85;


        for iday_lwp=1:length(days_lwp)

            date_str = [num2str(years_lwp(iyear_lwp)) num2str(months_lwp(imonth_lwp)) num2str(days_lwp(iday_lwp))];
            year_single = ['y' num2str(years_lwp(iyear_lwp))];
            month_single = ['m' num2str(months_lwp(imonth_lwp))];




%% Set a box to draw onto the map and to use to determine whether there is
%% enough satellite data
            % 12th Nov case
            %Need to come up with a new corner for teh 12th Nov case since this is
            %centred at 76W, 20S (not at 82W).
            % Domain runs from -78.93 to -73.08 W and -22.70 to -17.28 (based on the edges)
            LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER = [-78.93 -73.08];
            % Smaller region to account for boundary inflow of LWP and
            % spin-up during advection.
            %Looks like this mainly affects the south of the domain and to
            %the east (for 26th Oct POC case the east was also affected).
            %Also remove a bit for the boundary itself (around 0.25 deg
            %should be enough).
            LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];
            
            LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov            
%            LAT_val_DRIVER = [-25.70 -14.28]; LON_val_DRIVER =[-81.93 -70.08]; %Wider regions around UM domain         
%            LAT_val_DRIVER = [-30 -10]; LON_val_DRIVER = [-90 -70]; %GOES map Wider view2
            
            

%% Be sure to add a case here for the type of sat data - also need ot unzip
%% the files when reading daily ones.
            for isat=1:length(sats)

                sat = sats{isat};
                switch sat
                    case 'AMSRE'
                        filedir = '/home/disk/eos5/d.grosvenor/AMSRE/';
                        sat_str = 'amsre';
                        vers_str='v7';                        
                    case 'SSMI-f13'
                        filedir = '/home/disk/eos5/d.grosvenor/SSMI/f13/';
                        sat_str = 'f13';
                        vers_str='v7';
                    case 'SSMI-f15'
                        filedir = '/home/disk/eos5/d.grosvenor/SSMI/f15/';
                        sat_str = 'f15';
                        vers_str='v7';
                    case 'SSMI-f16'
                        filedir = '/home/disk/eos5/d.grosvenor/SSMI/f16/';
                        sat_str = 'f16';
                        vers_str='v7';
                    case 'SSMI-f17'
                        filedir = '/home/disk/eos5/d.grosvenor/SSMI/f17/';
                        sat_str = 'f17';
                        vers_str='v7';
                    case 'TMI'
                        filedir = '/home/disk/eos5/d.grosvenor/TMI/';
                        sat_str = 'tmi'; %The string that is actually in the filename
                        vers_str='v7.1';  
                    case 'Windsat'
                        filedir = '/home/disk/eos5/d.grosvenor/SSMI/windsat/';
                        sat_str = 'wsat'; %The string that is actually in the filename
                        vers_str='v7.0.1';                         
                    otherwise
                        fprintf(1,'\n\n*** This type of data has not been added yet!! ***\n\n');
                        return
                end
                
                sat_remss{isat} = sat_str;





                    %% Will first run :-
                    %   *** multi_read_amsre_daily ***, selecting right day

                    ioverride_read_amsre=1;
                    amsre_action='read_daily';
                    case_study = 'Choose days';

                    iday=1;
                    chosen_files(iday).name = [sat_str '_' date_str vers_str];  iday=iday+1;%12th Nov, 2008 for UM Boutle case - select 'Choose days' above

% --- Run script to get AMSRE/REMSS data :-                    
                    multi_read_amsre_daily

                    %Making a "new" gcm_str to avoid using timeseries3 type data in
                    %pdf2D
                    gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
                    gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
                    gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
                    gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;
                    %Actually can't really define this it varies over the globe for
                    %ASMRE
                    %         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre);
                    gcm_time_matlab_AMSRE2 = 0;
                    gcm_time_UTC_AMSRE2 = 0;
                    daynum_timeseries3_AMSRE2 = 1;
                    modisyear_timeseries3_AMSRE2 = 1;
                    
                    
                    for iasc_desc_DRIVER=1:2
                        
                        idriver = idriver + 1;

                    %asc_str = 'Ascending'; %Daytime
                    %asc_str = 'Descending'; %Nighttime

                    % switch asc_str
                    %     case 'Ascending'
                    %         iasc_desc_DRIVER = 1;
                    %     case 'Descending'
                    %         iasc_desc_DRIVER = 2;
                    % end

                    switch iasc_desc_DRIVER
                        case 1
                            asc_str = 'Ascending';  %Daytime
                        case 2
                            asc_str = 'Descending'; %Nighttime
                    end




                    % -- For option setting see inside the loops



                    %--- Load and process the data
                    %dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
                    %dirUM='/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/';
                    clear fileUM xdat_import ydat_import
                    idat=1;

                    fileUM{idat} = 'GOES'; labs_import(idat).l = sat; pole_lat=70; pole_lon=284; idat=idat+1;

                    for idat=1:length(fileUM)

                        %For AMSRE teh time of day depends on where the satellite was in its
                        %orbit, so just set the date here
                        time_driver = datenum(year_amsre,month_amsre,day_amsre);
                        nt_driver=1;



                        for it_driver=1:nt_driver

                            %         %------- Calculate the data to plot
                            %          %read in the GOES data for the specific time
                            %          ioverride_goes = 1;
                            %          goes_action = 'load a particular file';
                            %          read_GOES_vocals_netcdf_files

                            %--- run the file to set up the defaults
                            plot_global_maps_defaults

                            %--- set some options for these particular plot loops
                            irestrict_domain=1;
                            
                            LAT_val = LAT_val_DRIVER;
                            LON_val = LON_val_DRIVER;
                            
                            thresh_LAT=LAT_val;
                            thresh_LON=LON_val;

                            set_screening = {'none'};
                            modis_data_plot = 'Map of 2D data from outside driver script';

                            iset_min_clim=1;
                            clim_min=0;
                            iset_max_clim=1;
                            clim_max=300;

%                            isave_plot=0;

                            
                            %% Find the overpass time range
                            iloc = find(gcm_Plat2D_AMSRE2>LAT_val(1) & gcm_Plat2D_AMSRE2<LAT_val(2) & gcm_Plon2D_AMSRE2>LON_val(1) & gcm_Plon2D_AMSRE2<LON_val(2));
                            times = squeeze(time_amsre(:,:,1,iasc_desc_DRIVER));
                            min_time = minALL(times(iloc));
                            max_time = maxALL(times(iloc));
                            frac_covered(idriver)=length(find(isnan(times(iloc))==0)) ./ length(iloc);
                            
                            iloc_lat = find(lat_centres_sub>LAT_val(1) & lat_centres_sub<LAT_val(2)); nlat=length(iloc_lat);
                            iloc_lon = find(lon_centres_sub2>LON_val(1) & lon_centres_sub2<LON_val(2)); nlon=length(iloc_lon);                            
                            
%% Store some data for saving and future use
                            if idriver==1
                                %Make some empty arrays on first loop. One
                                %for each sat of dimensions [nlat nlon ndays nasc_desc nsat] 
                                % where nsat are the number of satellites     
                                iv=1; clear vars_remss
                                lwp_remss = NaN*ones([nlat nlon length(days_lwp) 2 length(sats)]);   vars_remss{iv}='lwp_remss'; iv=iv+1;                          
                                timeUTC_remss = NaN*ones([length(days_lwp) 2 length(sats)]);       vars_remss{iv}='timeUTC_remss'; iv=iv+1;                                                                                          
                                date_remss = NaN*ones([length(days_lwp) 2 length(sats)]);       vars_remss{iv}='date_remss'; iv=iv+1;
                                lwp_mean_remss = NaN*ones([length(days_lwp) 2 length(sats)]);   vars_remss{iv}='lwp_mean_remss'; iv=iv+1;   
                                lwp_std_remss = NaN*ones([length(days_lwp) 2 length(sats)]);   vars_remss{iv}='lwp_std_remss'; iv=iv+1;                               
                                
                                dat_present_remss = zeros([length(days_lwp) 2 length(sats)]);   vars_remss{iv}='dat_present_remss'; iv=iv+1;                                
                                vars_remss{iv}='sats'; iv=iv+1;
                                
                                gcm_Plat2D_REMSS = gcm_Plat2D_AMSRE(iloc_lat,iloc_lon);
                                gcm_Plon2D_REMSS = gcm_Plon2D_AMSRE(iloc_lat,iloc_lon);                                
                                gcm_Plat2D_edges_REMSS = gcm_Plat2D_edges_AMSRE(iloc_lat(1):iloc_lat(end)+1,iloc_lon(1):iloc_lon(end)+1);                                
                                gcm_Plon2D_edges_REMSS = gcm_Plon2D_edges_AMSRE(iloc_lat(1):iloc_lat(end)+1,iloc_lon(1):iloc_lon(end)+1);                                
                                
                                daynum_timeseries3_REMSS=[1:nmonths];
                                gcm_time_UTC_REMSS=zeros([1 nmonths]);
                                
                            end
                            
                            if frac_covered(idriver)>thresh_frac_covered
                                lwp_dat = 1e3*squeeze(lwp_amsre(iloc_lat,iloc_lon,:,iasc_desc_DRIVER));
                                lwp_remss(:,:,iday_lwp,iasc_desc_DRIVER,isat) = lwp_dat;
                                [lwp_mean_remss(iday_lwp,iasc_desc_DRIVER,isat), N_temp, lwp_std_remss(iday_lwp,iasc_desc_DRIVER,isat)] = meanNoNan(lwp_dat(:),1);                                
                                timeUTC_remss(iday_lwp,iasc_desc_DRIVER,isat) = min_time;
                                dat_present_remss(iday_lwp,iasc_desc_DRIVER,isat) = 1;                                
                                date_remss(iday_lwp,iasc_desc_DRIVER,isat) = datenum(years_lwp(iyear_lwp),months_lwp(imonth_lwp), days_lwp(iday_lwp));                                              
                            end

                            %Calculate the data to plot
                            dat_modis = 1e3*lwp_amsre(:,:,:,iasc_desc_DRIVER);  %N.B. - use this and not amsre_lwp!
                            %        dat_modis = goes_Vis;
                            
                            if i_mask_low_LWP
                                dat_modis(dat_modis<thresh_LWP_mask)=NaN;
                            end

                            %Set various things
                            time = time_driver(it_driver);
                            %Round to the nearest minute as sometimes get 18:59:59
                            time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy');
                            time_str = [time_str ' ' num2str(floor(min_time)) ':' num2str( 60*(min_time-floor(min_time)), '%02.0f' ) ' UTC'];
                            titlenam_driver = ['LWP for ' time_str ' ' asc_str ' ' labs_import(idat).l ' '];
                            units_str_plot = 'g m^{-2}';

                            mod_data_type='AMSRE';
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
                            plot_box_on_map
                            
                          %% Calculate the overall cloud fraction (or rather LWP>20 fraction)  
                            %The region being considered, cut out from the
                            %big array
                            dat_region = P(iloc_lat,iloc_lon);
                            inan=find(isnan(dat_region)==0);
                            CF(idriver) = length(inan)./length(dat_region(:));

                            %Save the figure
                            if isave_plot==1
                                savedir = savedir_driver;
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                                close(gcf);
                            end



                        end


                    end
                    %    xdat_import(idat).x =



                end %    for iasc_desc_DRIVER=1:2

            end % for isat


        end %iday_lwp
        
        if isave_vars_driver==1

            save(save_mat_file_remss,'LAT_val','LON_val','-V7.3');
            for ivar=1:length(vars_remss)
                save(save_mat_file_remss,vars_remss{ivar},'-V7.3','-APPEND');
            end

        end
        
        

