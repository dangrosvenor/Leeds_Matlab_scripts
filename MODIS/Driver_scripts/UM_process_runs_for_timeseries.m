%Top level script for processing UM runs to get out timeseries of the
%domain means.

% For each variable need to add a specific case in :-
% UM_maps_generic_time_loop_20141126T041216_v1.m

%% N.B. - these are NOT the UM files to run - those are set in the scripts
%% specified by these names
%UM_time_processing_case = '12th Nov';
%UM_time_processing_case = 'Iceland';
UM_time_processing_case = 'VOCALS processing runs';

% Run sets are located here :- UM_case_select_runs


% Select the lat/lon, UM run set, time range, etc. here :-  UM_SW_weighted_ACI_plots


% --- N.B. - add a case in UM_maps_generic_time_loop_20141126T041216_v1.m
% ---    for the requried variable !!

clear var_to_calc
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
%var_to_calc_multi{i}='accum_mass_total_column_to_z3000'; i=i+1;
%var_to_calc_multi{i}='accum_number_total_column_to_z3000'; i=i+1;
%var_to_calc_multi{i}='droplet_number_total_column_to_z3000'; i=i+1;
%var_to_calc_multi{i}='accum_mass_column_integrated'; i=i+1;

%var_to_calc_multi{i}='accum_number_column_integrated'; i=i+1;
%var_to_calc_multi{i}='droplet_number_column_integrated'; i=i+1;
%var_to_calc_multi{i}='aitken_number_column_integrated'; i=i+1;
%var_to_calc_multi{i}='coarse_number_column_integrated'; i=i+1;

% var_to_calc_multi{i}='accum_number_total_column_z1500_to_top'; i=i+1; % total aerosol number in the column up to 3km (/m2).
% var_to_calc_multi{i}='droplet_number_total_column_z1500_to_top'; i=i+1; % total droplet number in the column up to 3km (/m2).
% var_to_calc_multi{i}='coarse_number_total_column_z1500_to_top'; i=i+1; % total droplet number in the column up to 3km (/m2).

var_to_calc_multi{i}='air_mass_total_column_to_z1500'; i=i+1;

 var_to_calc_multi{i}='accum_number_total_column_to_z1500'; i=i+1;
 var_to_calc_multi{i}='coarse_number_total_column_to_z1500'; i=i+1;
 var_to_calc_multi{i}='droplet_number_total_column_to_z1500'; i=i+1;
% 
 var_to_calc_multi{i}='act_mass_liq_total_column_to_z1500'; i=i+1; % total aerosol number in the column up to 3km (/m2).
 var_to_calc_multi{i}='act_mass_rain_total_column_to_z1500'; i=i+1; % total droplet number in the column up to 3km (/m2).
% 
 var_to_calc_multi{i}='accum_mass_total_column_to_z1500'; i=i+1; % total aerosol number in the column up to 3km (/m2).
 var_to_calc_multi{i}='coarse_mass_total_column_to_z1500'; i=i+1; % total droplet number in the column up to 3km (/m2).

%var_to_calc_multi{i}='BL height using qL'; i=i+1;
var_to_calc_multi{i}='BL height using RH50';

% for CF will need to process a .mat file for the CF as for Lwp, etc.


for ivar_multi=1:length(var_to_calc_multi)

    %    UM_time_in.varname = var_to_calc_multi{ivar_multi};

    var_timser = var_to_calc_multi{ivar_multi};

    %%  Run the script to get the timeseries for all of the UM runs requested.

    switch UM_time_processing_case
        case '12th Nov'
            UM_SW_weighted_ACI_plots

        case 'Iceland'
            UM_settings_for_multi_var_timeseries_Iceland

        case 'VOCALS processing runs'
            UM_settings_for_multi_var_timeseries_VOCALS_processing_runs

    end

    %% Now save the timeseries to a .mat file in the run directory
    finfo = UM_out{1}.UM_case_out;
    for irun=1:length(finfo.fileUM)
        %Make the filename
        if iscell(finfo.dirUM)==1
            dirUM_i = finfo.dirUM{irun};
        else
            dirUM_i = finfo.dirUM;
        end
        save_file_timser = [dirUM_i finfo.fileUM{irun}];
        %        save_file_timser = [ remove_character(save_file_timser,'VAR_NAME',finfo.VAR_NAME_STR) '_timeseries.mat' ];
        save_file_timser = [ remove_character(save_file_timser,'VAR_NAME',var_timser) '_timeseries.mat' ];

        %Extract the relevant data
        timeseries_UM = UM_out{irun}.timeseries;
        time_UM = UM_out{irun}.time';
        time_mean_2D = UM_out{irun}.datUM_timemean{1};
        gcm_Plat2D_UM = UM_out{irun}.gcm_Plat2D_UM;
        gcm_Plon2D_UM = UM_out{irun}.gcm_Plon2D_UM;
        gcm_Plat2D_edges_UM = UM_out{irun}.gcm_Plat2D_edges_UM;
        gcm_Plon2D_edges_UM = UM_out{irun}.gcm_Plon2D_edges_UM;

        clear var_list
        i=1;
        var_list{i} = 'timeseries_UM'; i=i+1;
        var_list{i} = 'time_UM'; i=i+1;
        var_list{i} = 'UM_time_in'; i=i+1;
        var_list{i} = 'time_mean_2D'; i=i+1;
        var_list{i} = 'gcm_Plat2D_edges_UM'; i=i+1;
        var_list{i} = 'gcm_Plon2D_edges_UM'; i=i+1;
        var_list{i} = 'gcm_Plat2D_UM'; i=i+1;
        var_list{i} = 'gcm_Plon2D_UM'; i=i+1;        
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

end

