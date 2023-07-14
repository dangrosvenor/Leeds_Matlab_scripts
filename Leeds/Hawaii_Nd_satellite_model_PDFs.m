% Hawaii_Nd_satellite_model_PDFs.m
% Calls :-
% Hawaii_Nd_satellite_model_PDFs_FUNC.m (mostly a wrapper), which calls
% Hawaii_Nd_satellite_model_pdf (main script)
% Then also calls this for the final plotting of the combined PDF :-
% Hawaii_Nd_satellite_model_PDFs_PLOT.m

%Runs the PDF plotting function for various runs (e.g., different emission strenghts). The plotting function runs for 
%volc ON and OFF vs MODIS both in and out of plume - and plots the in and
%out of plume results in one plot.

% Probably should simplify to give the in and out of plume PDFs for each
% run (rather than having to do a pair each time...)

UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';

%var_str is such to match the filename as follows:-
%match_file = [UM_base_dir um_case_PI '/modis_' var_str '_matches_' coarse_str '.mat'];

%var_str = 'Nd';
var_str = 'Nd_Nd_N37';
var_str = 'Nd_N37_all_times';
var_str = 'Nd_N37_modis_4x4_coarse_grained_keep_NaNs';
var_str = 'Nd_N37_modis_4x4_coarse_grained_keep_NaNs_min_tau_3';
%var_str = 'CF_tau0pt3_cf_all_times';
%var_str = 'CF_tau0pt3_cf_tau0pt3_modis'; 
%var_str = 'CF_tau0pt3_cf_tau0pt3_modis_4x4_coarse_grained';
%var_str = 'CF_tau1_cf_tau1pt0_modis_4x4_coarse_grained';
%var_str = 'CF_tau2_cf_tau2pt0_modis_4x4_coarse_grained';

%var_str= 'CF_LWP5_cf_LWP5_modis_4x4_coarse_grained';
%var_str= 'CF_LWP10_cf_LWP10_modis_4x4_coarse_grained';
%var_str= 'CF_LWP20_cf_LWP20_modis_4x4_coarse_grained';
%var_str= 'CF_LWP30_cf_LWP30_modis_4x4_coarse_grained';
%var_str= 'LWP_LWP_modis_4x4_coarse_grained';
var_str= 'LWP_LWP_modis_4x4_coarse_grained_min_tau_3';
%var_str= 'tau_cloud_tau_modis_4x4_coarse_grained_zero_NaNs';
%var_str= 'tau_cloud_tau_modis_4x4_coarse_grained_keep_NaNs';
%var_str = 'total_max_random_cloud_amount_cf_all_times';
%var_str = 'low_CF_cf_all_times';
%var_str = 'mid_CF_cf_all_times';
%var_str = 'CF_LWP5_cf_all_times';
%var_str = 'SW_TOA_dummy variable';



%% 10x reduced emission rate (relative to u-ch765)
%um_case_PI = 'u-ch765';  %volcano OFF, orog ON
%um_case_PD = 'u-cr138';  %volcano ON, orog ON, 10x reduced emission rate (relative to u-ch765)

um_case_PI = 'u-cj086';  %volcano OFF, orog OFF
um_case_PD = 'u-cr139';  %volcano ON, orog OFF, 10x reduced emission rate (relative to u-ch765)

end_str = '_10x';
eval_str=['[output_in_plume' end_str ',output_out_plume' end_str ']' ...
    '= Hawaii_Nd_satellite_model_PDFs_FUNC(var_str,UM_base_dir,um_case_PI,um_case_PD);']; eval(eval_str);
%Function to do the final plots with in vs out of plume (model and modis)
%altogether
eval_str=['Hawaii_Nd_satellite_model_PDFs_PLOT(output_in_plume' end_str ',output_out_plume' end_str ')']; eval(eval_str);


%% End here
return

%% 5x reduced emission rate (relative to u-ch764)
um_case_PI = 'u-ch765';  %volcano OFF, orog ON
um_case_PD = 'u-co295';  %volcano ON, orog ON, 5x reduced emission rate (relative to u-ch764)

end_str = '_5x';
eval_str=['[xdat_model_PI' end_str ',ydat_model_PI' end_str ',xdat_model_PI_out' end_str ',ydat_model_PI_out' end_str ',xdat_model_PD' end_str ',ydat_model_PD' end_str ',xdat_model_PD_out' end_str ','...
    'ydat_model_PD_out' end_str ', xdat_modis_PD' end_str ',ydat_modis_PD' end_str ', xdat_modis_PD_out' end_str ',ydat_modis_PD_out' end_str ']' ...
    '= Hawaii_Nd_satellite_model_PDFs_FUNC(var_str,UM_base_dir,um_case_PI,um_case_PD);'];
eval(eval_str);


%% Make PDF of 5x reduced emission vs 10x

%% Plot PI, PD and MODIS together for in and out of plume
 figure('color','w'); hold on
 set(gca,'fontsize',16);
 clear leg_str
 i=1;
 plot(xdat_model_PD_5x(1).x,ydat_model_PD_5x(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON 5x less, in plume']; i=i+1;
 plot(xdat_model_PD_out_5x(1).x,ydat_model_PD_out_5x(1).y,'b--','linewidth',3); leg_str{i}=['Model Volc ON 5x less, out of plume']; i=i+1;
 plot(xdat_model_PD_10x(1).x,ydat_model_PD_10x(1).y,'g','linewidth',3); leg_str{i}=['Model Volc ON 10x less, in plume']; i=i+1;
 plot(xdat_model_PD_out_10x(1).x,ydat_model_PD_out_10x(1).y,'g--','linewidth',3); leg_str{i}=['Model Volc ON 10x less, out of plume']; i=i+1;
 
 %title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
 %             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
 
 ylabelstr='AOD 550nm'; 
 legend(leg_str);
 xlabel(ylabelstr);
 if logbin_norm_driver==1
     set(gca,'xscale','log');
     ylabel('df/dlog_{10}(x)');
 else
     ylabel('df/dx');
 end
 %set(gca,'xlim',[-5 500]);
 
 set(gca,'xscale','log');
 %set(gca,'yscale','log');
 set(gca,'xlim',[0 15]);
 title(thresh_str);
 grid on

         
         
 
 
 






