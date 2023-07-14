% -- Plot of SW vs Nd, along with the cloud properties (LWP, Nd and CF).
%  for the cloud properties will weight the averages by the mean SW_up, so
%  that the times of day that contribute most to the SW count more towards
%  the average.

%% Choose run sets
% Run sets are defined here :- UM_case_select_runs
% Choose the one required further below in this script

% Select the variables, lat/lon, UM run set, time range, etc. below


% --- N.B. - add a case in UM_maps_generic_time_loop_20141126T041216_v1.m
% ---    for the requried variable !!

isave_UM=1;

clear var_to_calc_multi
i=1;
% --- N.B. - add a case in UM_maps_generic_time_loop_20141126T041216_v1.m
var_to_calc_multi{i}='SW_up_TOA'; i=i+1;
%var_to_calc_multi{i}='SW_down_TOA'; i=i+1;
%var_to_calc_multi{i}='SW_down_surf_LWP_LT_0pt1'; i=i+1;
%var_to_calc_multi{i}='Transmission_down_surf_LWP_LT_0pt1'; i=i+1; %SW_down_surf ./ SW_down_TOA (fraction of SW that
%makes it to the surface)
var_to_calc_multi{i}='LWP'; i=i+1;
var_to_calc_multi{i}='RWP'; i=i+1;
var_to_calc_multi{i}='Nd'; i=i+1;
%var_to_calc_multi{i}='CF_0pt25_LWP_4km_20'; i=i+1;
var_to_calc_multi{i}='LWP_incloud_20gmsq'; i=i+1;
var_to_calc_multi{i}='CF_LWP_20'; i=i+1; %CF based on LWP that has not been coarsened and just calculating the CF over the whole domain
%var_to_calc_multi{i}='CF_LWP_150'; i=i+1; %CF based on LWP that has not been coarsened and just calculating the CF over the whole domain
%var_to_calc_multi{i}='LWP_incloud_150gmsq'; i=i+1;

% for CF will need to process a .mat file for the CF as for Lwp, etc.


for ivar_multi=1:length(var_to_calc_multi)

    %    UM_time_in.varname = var_to_calc_multi{ivar_multi};

    var_timser = var_to_calc_multi{ivar_multi};



    %% mat file for saving results to
    savefile_UM = ['/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/SW_values_for_ACI_' datestr(now,30) '.mat'];


    %% Run for SW first - this will give us the weights for averaging

    clear UM_time_in

    %% Set the details on what data to get

    %Choose the runs to plot/calculate for - see UM_case_select_runs
    UM_time_in.UM_cases = '12th Nov case, as of May 2016';
    UM_time_in.UM_cases = '12th Nov case, as of May 2016 adhoc';

    %Choose the variable
    %UM_time_in.varname = 'SW_up_TOA';
    UM_time_in.varname = var_to_calc_multi{ivar_multi}; %from outside script
    UM_time_in.var_units_str = ''; %may want to set this where set nc_varname, etc.

    %Choose the lat lon required
    UM_time_in.LAT_val_DRIVER = [-22.70 -17.28]; UM_time_in.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
    UM_time_in.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

    %Choose times
    %For time specify either a range (will pick out all times within the range)
    %or a specific set of times (can be a single time). BUT DON'T SET NOT BOTH!!
    %UM_time_in.time_range =  [datenum('01-Sep-2014 12:00') datenum('01-Sep-2014 12:00')];
    %UM_time_in.time_specific(1) =  datenum('01-Sep-2014 12:00');
    %UM_time_in.time_specific =  [datenum('12-Nov-2008 22:00') datenum('13-Nov-2008 12:00')];
    UM_time_in.time_range =[ datenum('12-Nov-2008 06:00')  datenum('14-Nov-2008 00:00') ];

    %% Set the flags
    UM_time_in.action='timeseries domain means';  %save_all
    UM_time_in.noplot = 0; %set to 1 for timeseries domain means for faster processing
    UM_time_in.iclose_figs=1; %whether to close the figures after plotting
    UM_time_in.isave_plot=0;
    UM_time_in.iplot_maps=1;

    UM_time_in.i_mask_low_LWP=0; %Make any values below thresh_LWP equal to NaN
    UM_time_in.thresh_LWP_mask = 20;

    %Stuff for coarsening
    UM_time_in.icoarsen = 0;
    %UM_time_in.dlat_target=dlat_target; UM_time_in.dlon_target = dlon_target;

    %Plotting options
    UM_time_in.iset_min_clim=1;
    UM_time_in.clim_min=0;
    UM_time_in.iset_max_clim=1;
    UM_time_in.clim_max=650;

    UM_time_in.iplot_markers=0; %marker points on the map

    %% Run the time looping function
    [UM_out] = UM_maps_generic_time_loop_20141126T041216_v1(UM_time_in);




    if isave_UM==1

        switch UM_time_in.action
            case 'timeseries domain means'


                %% Now save the timeseries to a .mat file in the run directory
                finfo = UM_out{1}.UM_case_out;
                for irun=1:length(finfo.fileUM)
                    %Make the filename
                    save_file_timser = [finfo.dirUM finfo.fileUM{irun}];
                    %        save_file_timser = [ remove_character(save_file_timser,'VAR_NAME',finfo.VAR_NAME_STR) '_timeseries.mat' ];
                    save_file_timser = [ remove_character(save_file_timser,'VAR_NAME',var_timser) '_timeseries.mat' ];

                    %Extract the relevant data
                    timeseries_UM = UM_out{irun}.timeseries;
                    time_UM = UM_out{irun}.time';

                    clear var_list
                    i=1;
                    var_list{i} = 'timeseries_UM'; i=i+1;
                    var_list{i} = 'time_UM'; i=i+1;
                    var_list{i} = 'UM_time_in'; i=i+1;
                    %         var_list{i} = 'gcm_Plat2D_edges_UM'; i=i+1;
                    %         var_list{i} = 'gcm_Plon2D_edges_UM'; i=i+1;
                    %         var_list{i} = 'it'; i=i+1;
                    %         var_list{i} = 'daynum_timeseries3_UM'; i=i+1;
                    %         var_list{i} = 'modisyear_timeseries3_UM'; i=i+1;
                    %         var_list{i} = 'gcm_time_UTC_UM'; i=i+1;
                    %         var_list{i} = 'gcm_time_matlab_UM'; i=i+1;

                    for i=1:length(var_list)
                        if i==1;
                            app_str='';
                        else
                            app_str=[',''-APPEND'''];
                        end
                        eval(['save(''' save_file_timser ''',''' var_list{i} ''',''-V7.3''' app_str ');']);
                    end
                    
                end


                    otherwise

                        %% Save the results and the options used to get it
                        save(savefile_UM,'UM_out','UM_time_in','-V7.3');

             


        end


    end




end