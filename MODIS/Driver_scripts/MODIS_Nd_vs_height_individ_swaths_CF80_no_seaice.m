%Run pdf2D_plot_commands to produce Nd vs CTH plots
% Saves this data using save_vars_mat_run

fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot

%        modis_data_plots={'MODIS LWP minus AMSRE'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';
%        modis_data_plots={'Number of droplets cell values time mean - specific days'}; fraction_plotcase='Fraction of days remaining after screening';
%         modis_data_plots={'LWP AMSRE time3'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';

% -- Specify the plots to call in plot_global etc.
%        modis_data_plots={'Dummy string'}; %just put somethign in here to make the loop work
modis_data_plots ={'Number of droplets 1.6um cell values time mean - specific days', ...
    'Number of droplets cell values time mean - specific days', ...
    'Number of droplets 3.7um cell values time mean - specific days', ...
    };

%        zero_CF_out_of_range=0;
time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';

set_screening = {'none','NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
%Mean CTH:-
set_screening = {'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
%Max CTH:-
%        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};
%Mean CTH with sea-ice screening:-
set_screening = {'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff + seaice'};
%        set_screening = {'Dummy string'}; %Dummy string to enter loop, but will set this later - acutally might not
%be used for the GCM plot

%set_screening = {'none'}; %if using mockL3 would generally use this as the screening has already
%been done. However, if using concatenated L3 from separate overpasses then would want
%to set some screening


%set a default screening - then can change individual ones in the loop
screening_eval_str=[ ...
    'thresh_ndays=1;' ...
    'thresh_SZA=[0 65];' ...
    'thresh_CF=[0.8 1.00001];' ...
    'thresh_NP=50;' ...
    'thresh_sensZA=[0 90];'     ...
    'thresh_CTT = [173 273+100];'  ...
    'thresh_relAZ = [0 180];'   ...
    'thresh_CTH = [-0.01 1e9];' ...
    'thresh_stdW = [0 1e9];' ...
    'thresh_sigCTT = [0 1e9];' ...
    'thresh_reff = [0 30];' ...
    'thresh_seaice=[-1 1e-9];'...
    ];




%just ones that alter it from the above baseline
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];'};
screen_edits = {'thresh_CF=[0.8 1.00001]','thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];','thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];','thresh_CF=[0.8 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_reff=[0 30];'};
screen_edits = {''};

screen_edits = {'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %all CF, all re
%screen_edits = {'thresh_reff=[0 30];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %CF>80 and all CF (all re)

screen_edits = {'thresh_CF=[0.1 0.4]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
    ,'thresh_CF=[0.4 0.6]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
    ,'thresh_CF=[0.6 0.8]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
    ,'thresh_CF=[0.8 1.00001]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
    };  %ranges of CF (applies to model too when choosing re CF range)

screen_edits = {''};

set_years = [2007];

ifilter_ndays=0;




set_MOD = {'MOD06 liquid','MOD35'};
set_MOD = {'MOD35'};  %redundant for Nd
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
set_real_abs = {'absolute','percentage'};
set_real_abs = {'percentage'};   %redundant for Nd
%        set_real_abs = {'absolute'};

%what to plot for the black contours. - %redundant for Nd
cont_dat_choose = 'calipso mid+highCF';

%what to use for filtering - %redundant for Nd
filt_dat_choose = 'calipso mid+highCF';
ifilter_clhcalipso=1; %whether to do the filtering or not
thresh_clh = [-0.01 0.3];

iuse_MODIS_CF_error=0; %if =1 then uses both the AMSRE comparison error and the change due to changing CF
%from MODIS to plot the contours for the more dubious area.
%Otherwise just uses the LWP error. Decided that the latter is
%better as there will be a real component ot the change in Nd with
%CF.

iuse_saved_mask=0; %redundant for Nd
saved_mask_name = 'imask_ndays15_cal_midhigh';

clear RMSE_DRIVER mean_bias_DRIVER season_mean_driver season_Ndatap_driver season_std_driver mean_sig2_driver Nspatial_driver season_timestd_driver season_datatype_driver season_vals_timelabs_driver screen_type_driver

icount_DRIVER=0;

for iloop_plot=1:length(modis_data_plots)
    modis_data_plot = modis_data_plots{iloop_plot};

    for iscreen_set=1:length(set_screening)
        screen_type = set_screening{iscreen_set};

        for iscreen_thresh=1:length(screen_edits)
            eval(screening_eval_str); %eval the basline screening
            eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline

            for iloop_MOD=1:length(set_MOD)
                MODCF_type = set_MOD{iloop_MOD}

                for iloop_abs=1:length(set_real_abs)
                    error_type= set_real_abs{iloop_abs};
                    for iloop_years=1:length(set_years)
                        years_required_for_mean = set_years(iloop_years)
                        icount_DRIVER = icount_DRIVER+1;


                        %             modis_data_plot = modis_data_plots{};  %set above



                        ioverride_plotglobal_thresh=1;
                        %                      iocean_only=1;
                        ioverride_time_selection=1;
                        ioverride_plotglobal_loc=1;
                        mod_data_type='timeseries3 lambert';
                        %                            modis_data_plot = modis_data_plots{iloop_plot};

                        ioverride_years_time_screen=1; %required for the different years
                        ioverride_monthly_options=1;
                        multi_case = 'Monthly';
                        plot_script = 'plot_global_maps';
                        ioverride_location_selection2=1;
                        LAT_vals = [-70 -60;]; LON_vals = [-180 -40;]; %Southern Ocean


                        iset_min_clim=0;
                        clim_min=0;
                        iset_max_clim=0;
                        %        clim_max=225;
                        clim_max=150;

                        %                            times_required = [0:24]; %need to specify all hours (in 1 hour increments)
                        %                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                        % The above should be set in monthly_means_from_plot_global

                        isave=0; %flag for whether monthly_means_from_plot_global should plot or not



                        % plot_global_maps
                        monthly_means_from_plot_global
                        %                            Psave_MODIS = P;
                        %                            Npoints_save_MODIS = Npoints;
                        %                            thresh_str_MODIS = thresh_str;
                        %                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
                        %                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                        thresh_str_driver{icount_DRIVER} = thresh_str; %save the thresh_str used
                        screen_type_driver{icount_DRIVER} = screen_type;
                        years_required_for_mean_driver{icount_DRIVER} = years_required_for_mean;

                        season_mean_driver{icount_DRIVER} = season_mean;
                        season_Ndatap_driver{icount_DRIVER} = season_Ndatap;
                        season_std_driver{icount_DRIVER} = season_std;
                        mean_sig2_driver{icount_DRIVER} = mean_sig2;
                        Nspatial_driver{icount_DRIVER} = Nspatial;
                        season_timestd_driver{icount_DRIVER} = season_timestd;

                        season_datatype_driver{icount_DRIVER} = season_datatype;
                        season_vals_timelabs_driver{icount_DRIVER} = season_vals_timelabs;

                        if exist('thresh_str_multiL2')
                            thresh_str_multiL2_driver{icount_DRIVER} = thresh_str_multiL2;
                        else
                            thresh_str_multiL2_driver{icount_DRIVER} = '';
                        end

                        if exist('LON_str')
                            season_LON_str{iseason} = LON_str;
                            season_LAT_str{iseason} = LAT_str;
                        end



                        %                            if isave_vals==1
                        %                                fprintf(fid_mtab,'%f ',mean_bias_MODIS);
                        %                               fprintf(fid_Rtab,'%f ',RMSE_MODIS);
                        %                            end

                        if isave_plot==1
                            saveas_ps_fig_emf(gcf,[savename],'',0,1);
                        end


                    end

                end

            end

        end

    end
end

            