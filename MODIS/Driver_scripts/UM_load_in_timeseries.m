% Plot domain and time mean LWP and RWP vs Nd
% Data is loaded from a .mat files for the timeseries, which are generated
% by UM_process_runs_for_timeseries and stored in each directory

isave_plot_overall=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

% --- Add run cases below to UM_case_select_runs
UM_cases = '12th Nov case ACI, as of May 2016';
UM_cases = '12th Nov case ACI, as of May 2016 with sed runs'; %N.B - messes up the bar chart for SW etc.
            %Also changes idat_micro to only plot line for the
            %sedimentation runs
UM_cases =  '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM';

% -- Calls UM_case_select_runs
UM_case_select_RUN; %external script to select the case and put all the variables into current memory space


% save name for the figure generated here :-
savename_ACI = '/home/disk/eos1/d.grosvenor/modis_work/plots/UM/Processing_LWP_vs_Nd';


nsub=0; %counter for which plot we are on


% vars to load timeseries for
clear var_list xdat_save ydat_save
i=1;
%var_list{i} = 'SW_up_TOA'; i=i+1;
var_list{i} = 'LWP'; i=i+1;
%var_list{i} = 'LWP_incloud_20gmsq'; i=i+1;
%var_list{i} = 'RWP'; i=i+1;
%var_list{i} = 'Nd'; i=i+1;
%var_list{i} = 'CF_LWP_20'; i=i+1;
%var_list{i} = 'CF_0pt25_LWP_4km_20'; i=i+1;
%var_list{i} = 'UM_time_in'; i=i+1;
var_list{i} = 'accum_number_total_column_to_z1500'; i=i+1;
var_list{i} = 'droplet_number_total_column_to_z1500'; i=i+1;
var_list{i} = 'coarse_number_total_column_to_z1500'; i=i+1;

time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC) for VOCALS (add this on for LST)


%subplotting settings
subplotting_DRIVER=1;
xsub=length(var_list)-1;
ysub=1;

%Load all the timeseries into e.g. SW_up_TOA(irun).timeseries
for irun=1:length(fileUM)    
    for ivar=1:length(var_list)
        var_str = var_list{ivar};
        if irun==1
            eval(['clear ' var_str]);
        end
        filename_tim = [remove_character( [dirUM{irun} fileUM{irun}], 'VAR_NAME',var_str) '_timeseries.mat'];
        eval_str=[var_str '(irun)=load(filename_tim);'];
        eval(eval_str);        
    end    
end

fprintf(1,'\n Done load of timeseries');






