function [output_in_plume,output_out_plume] ...
    = Hawaii_Nd_satellite_model_PDFs_FUNC(var_str,UM_base_dir,um_case_PI,um_case_PD)


%% Load Nd and SO2 data
%var_str = 'Nd';

coarse_str='1x1';
switch coarse_str
    case '1x1'
        coarse_str2 = '';
    otherwise
        coarse_str2=['_' coarse_str];
end
match_file = [UM_base_dir um_case_PI '/modis_' var_str '_matches_' coarse_str '.mat'];
matches_Nd_1x1_PI = load(match_file);
%Load the full so2 data (i.e., all times) from the .mat file.
so2_save_file = [UM_base_dir um_case_PI '/so2.mat'];
model_so2_PI = load(so2_save_file);
eval_str=['model_dat_so2_PI = model_so2_PI.dat' coarse_str2 ';']; eval(eval_str);

match_file = [UM_base_dir um_case_PD '/modis_' var_str '_matches_' coarse_str '.mat'];
matches_Nd_1x1_PD = load(match_file);
%Load the full so2 data (i.e., all times) from the .mat file.
so2_save_file = [UM_base_dir um_case_PD '/so2.mat'];
model_so2_PD = load(so2_save_file);
eval_str=['model_dat_so2_PD = model_so2_PD.dat' coarse_str2 ';']; eval(eval_str);



%In-plume (defined by SO2 column>1e-5)
so2_thresh = [1e-5 1e99];
[output_in_plume] = ...
    Hawaii_Nd_satellite_model_pdf(var_str,matches_Nd_1x1_PD,model_dat_so2_PD,matches_Nd_1x1_PI,so2_thresh);

%Out of plume
so2_thresh = [-1e99 1e-5];
[output_out_plume] = ...
    Hawaii_Nd_satellite_model_pdf(var_str,matches_Nd_1x1_PD,model_dat_so2_PD,matches_Nd_1x1_PI,so2_thresh);




% %% Plot PI, PD and MODIS together for in and out of plume
%  figure('color','w'); hold on
%  set(gca,'fontsize',16);
%  clear leg_str
%  i=1;
%  plot(output_in_plume.xdat_model_PI(1).x,output_in_plume.ydat_model_PI(1).y,'k','linewidth',3); leg_str{i}=['Model Volc OFF, in plume']; i=i+1;
%  plot(output_out_plume.xdat_model_PI(1).x,output_out_plume.ydat_model_PI(1).y,'k--','linewidth',3); leg_str{i}=['Model Volc OFF, out of plume']; i=i+1;
%  plot(output_in_plume.xdat_model_PD(1).x,output_in_plume.ydat_model_PD(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON, in plume']; i=i+1;
%  plot(output_out_plume.xdat_model_PD(1).x,output_out_plume.ydat_model_PD(1).y,'b--','linewidth',3); leg_str{i}=['Model Volc ON, out of plume']; i=i+1;
%  plot(output_in_plume.xdat_modis_PD(1).x,output_in_plume.ydat_modis_PD(1).y,'r','linewidth',3); leg_str{i}=['MODIS, in plume']; i=i+1;
%  plot(output_out_plume.xdat_modis_PD(1).x,output_out_plume.ydat_modis_PD(1).y,'r--','linewidth',3); leg_str{i}=['MODIS, out of plume']; i=i+1;
%  %title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
%  %             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
%  
%  %ylabelstr='AOD 550nm'; 
%  %ylabelstr='N_d (cm^{-3})'; 
%  legend(leg_str);
%  xlabel(output_in_plume.xlabelstr);
%  ylabel(output_in_plume.ylabelstr);
%  
%  %set(gca,'xlim',[-5 500]);
%  
%  set(gca,'xscale','log');
%  %set(gca,'yscale','log');
%  set(gca,'xlim',[0 1000]);
%  title(thresh_str);
%  grid on

         
         
 
 
 

