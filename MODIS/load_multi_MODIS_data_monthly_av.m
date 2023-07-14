function test()
% Runs multiple instances of MODIS_multi_DAY_processL3L2 for mutliple years
% (one at a time).
% Set MODIS_multi_DAY_processL3L2 with make_mockL3_variables_script_name = 'make_mockL3_variables_just_for_Nd_May2017';
% and multiL2L3_case = 'load L3 and concatenate'; for one year.
clear gca

isave_plots=1;
clims = [50 150];
ibrown_cbar=0;

iload_seaice_screened = 1; %Whether to load the data that has been screened for sea-ice.
iload_seaice_screened = 0; %Not currently available for CF>0 data.

load_method = 'MODIS_multi_DAY_processL3L2'; %loading the mock L3 data from daily files
load_method = 'saved_daily_files'; %loading more processsed versions of the above (e.g. as created for online datasets - may also have sea-ice screening for example)
%load_method = 'saved_monthly_files'; %loading more processsed versions of the above (e.g. as created for online datasets - may also have sea-ice screening for example)

iL3=1; %Whether are using Level-3 data (otherwise mock L3)

CF_str = '0.8'; CF_str2 = '80';
CF_str = '0.0'; CF_str2 = '0'; %N.B. - don't currently have the Cloud_Fraction_Liquid variable for CF>0 dataset for mock L3...

av_method = 'equal weighting to each year'; Ndays_thresh = 5; %For load_method = 'saved_monthly_files' - Min no. days per month to consider a monthly average valid
%av_method = 'weight by days';


seaice_str='';
seaice_max_2week_str='';

%ONly have the sea-ice screened data for up to 2014 not 2015
%direcs={'aqua'};  years_multi=[2003:2015]; days_multi={[1:365],[1:366],[1:365],[1:365],[1:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365],[1:365]}; %,[1:365]};%
%direcs={'aqua'};  years_multi=[2014:2015]; days_multi={[1:365],[1:365]}; %,[1:365]};%
direcs={'aqua'};  years_multi=[2003:2014]; days_multi={[1:365],[1:366],[1:365],[1:365],[1:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365],[1:365]}; %,[1:365]};%
%direcs={'aqua'};  years_multi=[2009]; days_multi={[1:365]}; %,[1:365]};%

direcs={'aqua'};  years_multi=[2016]; days_multi={[336:366]};%
direcs={'aqua'};  years_multi=[2017]; days_multi={[1:365]};%

icount=0;
for iyear_Nd_dataset = 1:length(years_multi)
    icount=icount+1;
    
    years = years_multi(iyear_Nd_dataset);
    days = days_multi(iyear_Nd_dataset);
    switch load_method
        case 'MODIS_multi_DAY_processL3L2'
            %Run the script to load the data
            i_override_MODIS_multi=1;
            file_dir=['/home/disk/eos15/d.grosvenor/mock_L3/CF_' CF_str '_meanCTT_173_meanCTH_3.2km_SZA_65/'];
            MODIS_multi_DAY_processL3L2
            %%
        case 'saved_daily_files'
            if iL3==1
                file_dir01 = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L3/Nd_from_L3/';
            else
                file_dir01 = '/home/disk/eos1/d.grosvenor/mock_L3/';
            end
            
            file_dir_save=[file_dir01 'CF_' CF_str '_meanCTT_173_meanCTH_3.2km_SZA_65/'];            
            file_dir = file_dir_save;
            
            if iload_seaice_screened == 1
                seaice_str = '_screened_for_seaice_';
                seaice_max_2week_str = '_2week_max';
            else
                if iL3~=1                    
                    file_dir = ['/home/disk/eos15/d.grosvenor/mock_L3/CF_' CF_str '_meanCTT_173_meanCTH_3.2km_SZA_65/'];
                end
                seaice_str = '';
                seaice_max_2week_str = '';
            end
            
            channel_str = '21';
            channel_str = '37';
            
            str_1km_1deg = '1km';  %Based on 1km tau and reff data.
            str_1km_1deg = '1deg'; %Based on the 1x1 deg re and tau - but note that originally the files were not labelled as '_1deg' for this (e.g. just Nd_monthly_37_2015_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat)
            
            
            %The Nd_re files are daily files - here we just average them to
            %monthly and re-package.
            file_load = [file_dir 'Nd_re_tau_cf_' channel_str '_' num2str(years) '_SZA_LT_65_CF_GT_' CF_str2 '_CTH_LT_3.2km' seaice_str seaice_max_2week_str '.mat'];
            modisL3_dat = load(file_load);
            
            modis_datenum = datenum(modisL3_dat.time(:,1),modisL3_dat.time(:,2),modisL3_dat.time(:,3));
            modisyear_timeseries3 = modisL3_dat.time(:,1)';
            daynum_timeseries3_MODIS = day_of_year_from_date_func(modis_datenum);
            mod_data_type = 'timeseries3';
            MLAT = modisL3_dat.lat; MLAT = flipdim(MLAT,2); LAT=MLAT; %MLAT seems to be the wrong way around
            MLON = modisL3_dat.lon; LON = MLON;
            [lon2d,lat2d] = meshgrid(LON,LAT);
            make_mockL3_variables_script_name = 'Loaded from .mat files using load_multi_MODIS_data_monthly_av.m with load_method = ''saved_daily_files'' ';
            multiL2L3_case=''; multiL2L3_project=''; daily_averaged_files_loc2=''; modis_data_plot='';
            
            ioverride_monthly_options=1;
            %set the values controlled by ioverride_monthly_options
            multi_case = 'Monthly';
            plot_script = 'plot_global_maps';
            isave=0;
            %--- run the file to set up the defaults
            plot_global_maps_defaults
            irestrict_domain=0; %Whether to restrict to a domain or not
            data_select='specific_modis_data';
            years_required_for_mean = unique(modisyear_timeseries3);
            screen_type ='none';
            %    proj_type = 'global oval';
            ifilter_ndays=0;
            
            Notes = modisL3_dat.Notes;
            
            switch channel_str
                case '21'
                    switch str_1km_1deg
                        case '1km';
                            N_time3 = permute(modisL3_dat.Nd_1km,[2 1 3]); %1km Nd
                        case '1deg'; %Based on the 1x1 deg re and tau - but note that originally the files were not labelled as '_1deg' for this (e.g. just Nd_monthly_37_2015_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat)
                            N_time3 = permute(modisL3_dat.Nd,[2 1 3]); %1deg Nd
                    end
                    modis_data_plot = 'Number of droplets cell values time mean'; %using 1x1 deg re and tau to make Nd
                    
                    
                case '37'
                    switch str_1km_1deg
                        case '1km';
                            N_time3_37 = permute(modisL3_dat.Nd_1km,[2 1 3]); %1km Nd
                        case '1deg'; %Based on the 1x1 deg re and tau - but note that originally the files were not labelled as '_1deg' for this (e.g. just Nd_monthly_37_2015_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat)
                            N_time3_37 = permute(modisL3_dat.Nd,[2 1 3]); %1deg Nd
                    end
                    modis_data_plot = 'Number of droplets cell values 3.7um time mean'; %using 1x1 deg re and tau to make Nd
                    
            end
            
            
        case 'saved_monthly_files'
            seaice_str = '_screened_for_seaice_';
            seaice_max_2week_str = '_2week_max';
            channel_str = '21';
            channel_str = '37';
            
            
            switch channel_str
                case '21'
                    channel_str2 = '2.1 \mum';
                case '37'
                    channel_str2 = '3.7 \mum';
            end
            
            
            
            
            file_dir=['/home/disk/eos1/d.grosvenor/mock_L3/CF_' CF_str '_meanCTT_173_meanCTH_3.2km_SZA_65/'];
            
            file_load = [file_dir 'Nd_monthly_' channel_str '_' num2str(years) '_SZA_LT_65_CF_GT_' CF_str2 '_CTH_LT_3.2km' seaice_str seaice_max_2week_str '.mat'];
            modisL3_dat = load(file_load);
            
            if icount==1
                Nd_multi_annual = zeros(size(modisL3_dat.Nd_1deg_mean));
                weight = zeros(size(modisL3_dat.Nd_1deg_mean));
            end
            
            
            switch av_method
                case 'equal weighting to each year'
                    i = find(modisL3_dat.Nd_1deg_Ndatap < Ndays_thresh);
                    Nd_multi_annual(i) = NaN; %remove points where don't have at least Ndays_thresh days in every month
                    Nd_multi_annual  = Nd_multi_annual + modisL3_dat.Nd_1deg_mean; %running sum
                    weight = weight + 1;
                case 'weight by days'
                    Nd_multi_annual  = Nd_multi_annual + modisL3_dat.Nd_1deg_mean .* modisL3_dat.Nd_1deg_Ndatap;
                    weight = weight + modisL3_dat.Nd_1deg_Ndatap;
            end
            
    end
    
    %%
    switch load_method
        case {'MODIS_multi_DAY_processL3L2','saved_daily_files'}
            %Run the monthly mean generator - set irestrict_domain=0 if want global
            %data
            monthly_means_from_plot_global
            close all
            
            %    [date_str,date_num] = date_from_day_of_year_func(daynum_timeseries3,modisyear_timeseries3);
            %    time = datevec(date_num); time=time(:,1:3);
            lat=LAT;
            lon=LON;
            
            year=years_required_for_mean*ones([1 12]);
            month=[1:12];
            
            clear var_str
            var_str{1} = ['Nd_' str_1km_1deg '_mean']; eval_str = [var_str{1} '= season_mean_ALL;']; eval(eval_str);
            var_str{2} = ['Nd_' str_1km_1deg '_Ndatap'];eval_str = [var_str{2} '= season_Ndatap_ALL;']; eval(eval_str);
            var_str{3} = ['Nd_' str_1km_1deg '_std_dev'];eval_str = [var_str{3} '= season_std_ALL;']; eval(eval_str);
            %Nd_1km = permute(Droplet_Number_Concentration_37.timeseries3,[2 1 3]);
            file_save = [file_dir_save 'Nd_monthly_' channel_str '_' str_1km_1deg '_' num2str(year(1)) '_SZA_LT_65_CF_GT_' CF_str2 '_CTH_LT_3.2km' seaice_str seaice_max_2week_str '.mat'];
            
            save(file_save,'-V7.3',var_str{1},var_str{2},var_str{3},'lon','lat','lat2d','lon2d','year','month','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2','modis_data_plot','Notes');
            mat2nc_Dan(file_save,[file_save '.nc']);
            
    end
    
    
end


%% Plot all 12 months (mulit-annual mean)
switch load_method
    case 'saved_monthly_files'
        
        modis_datenum = datenum(modisL3_dat.year,modisL3_dat.month,1);
        modisyear_timeseries3 = modisL3_dat.year;
        daynum_timeseries3_MODIS = day_of_year_from_date_func(modis_datenum);
        mod_data_type = 'timeseries3';
        MLAT = modisL3_dat.lat; LAT = MLAT;
        MLON = modisL3_dat.lon; LON = MLON;
        make_mockL3_variables_script_name = 'Loaded from .mat files using load_multi_MODIS_data_monthly_av.m with load_method = ''saved_monthly_files'' ';
        multiL2L3_case=''; multiL2L3_project=''; daily_averaged_files_loc2=''; modis_data_plot='';
        
        for imon=1:12
            
            mon_str = datestr(datenum(2008,imon,1),'mmm');
            
            dat_modis = Nd_multi_annual(:,:,imon) ./ weight(:,:,imon); %permute(modisL3_dat.Nd_1km,[2 1 3]);
            
            
            
            
            
            %--- run the file to set up the defaults
            plot_global_maps_defaults
            ioverride_plotglobal_loc=1;
            ioverride_plotglobal_thresh=1;  %comment out if want to use screenings set in plot_global
            titlenam_driver = ['N_d multi-year ' mon_str ' mean, ' channel_str2 ' ' av_method];
            %time override should aready be set (ioverride_time_selection)
            ioverride_years_time_screen=1; %required to specify the different years
            inew_cticks=0;  %colorbar is of the non-linear type
            %modis_data_plot = 'Generic plot specified outside of script';
            modis_data_plot = 'Map of 2D data from outside driver script';
            proj_type='polar'; stereo_str01='lat_polar=-90'; stereo_str02='m_proj(''stereographic'',''lat'',lat_polar,''lon'',-98,''rad'',55)';
            plot_global_maps
            
            if ibrown_cbar==1
                lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
                caxis(clims);
            else
                caxis(clims); caxis(clims);
            end
            
            if isave_plots==1
                savename=['/home/disk/eos1/d.grosvenor/modis_work/Nd_Antarctica/' titlenam_driver];
                clear opts
                %    opts.iplot_png=1;  %Don't set this since this creates a PNG using
                %    Matlab's driver - but better to convert from eps since otherwise font
                %    is too small
                opts.iplot_eps=1;
                saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
                
            end
            
            
        end
        
end