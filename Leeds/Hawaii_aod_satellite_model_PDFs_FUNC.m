function [xdat_model_PI,ydat_model_PI,xdat_model_PI_out,ydat_model_PI_out,xdat_model_PD,ydat_model_PD,xdat_model_PD_out,...
    ydat_model_PD_out, xdat_modis_PD,ydat_modis_PD, xdat_modis_PD_out,ydat_modis_PD_out] ...
    = Hawaii_aod_satellite_model_PDFs_FUNC(UM_base_dir,um_case_PI,um_case_PD)


%% Plot AOD PDFs
var_str = 'aod';


coarse_str='3x3';
match_file = [UM_base_dir um_case_PI '/modis_' var_str '_matches_' coarse_str '.mat'];
matches_aod_3x3_PI = load(match_file);

coarse_str='3x3';
match_file = [UM_base_dir um_case_PD '/modis_' var_str '_matches_' coarse_str '.mat'];
matches_aod_3x3_PD = load(match_file);

var_str = 'so2';

coarse_str='3x3';
match_file = [UM_base_dir um_case_PI '/modis_' var_str '_matches_' coarse_str '.mat'];
matches_so2_3x3_PI = load(match_file);

coarse_str='3x3';
match_file = [UM_base_dir um_case_PD '/modis_' var_str '_matches_' coarse_str '.mat'];
matches_so2_3x3_PD = load(match_file);

%In-plume (defined by SO2 column>1e-5)
so2_thresh = [1e-5 1e99];
[output_IN,xdat_model_PD,ydat_model_PD,xdat_model_PI,ydat_model_PI,xdat_modis_PD,ydat_modis_PD,ylabelstr,logbin_norm_driver,thresh_str] = Hawaii_aod_satellite_model_pdf(matches_aod_3x3_PD,matches_so2_3x3_PD,matches_aod_3x3_PI,so2_thresh);

%Out of plume
so2_thresh = [-1e99 1e-5];
[output_OUT,xdat_model_PD_out,ydat_model_PD_out,xdat_model_PI_out,ydat_model_PI_out,xdat_modis_PD_out,ydat_modis_PD_out,ylabelstr,logbin_norm_driver,thresh_str] = Hawaii_aod_satellite_model_pdf(matches_aod_3x3_PD,matches_so2_3x3_PD,matches_aod_3x3_PI,so2_thresh);




%% Plot PI, PD and MODIS together for in and out of plume
 figure('color','w'); hold on
 set(gca,'fontsize',16);
 clear leg_str
 i=1;
 plot(xdat_model_PI(1).x,ydat_model_PI(1).y,'k','linewidth',3); leg_str{i}=['Model Volc OFF, in plume \mu=' num2str(output_IN.Y_mean_overall_OFF,'%.2f')]; i=i+1;
 plot(xdat_model_PI_out(1).x,ydat_model_PI_out(1).y,'k--','linewidth',3); leg_str{i}=['Model Volc OFF, out of plume \mu=' num2str(output_OUT.Y_mean_overall_OFF,'%.2f')]; i=i+1;
 plot(xdat_model_PD(1).x,ydat_model_PD(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON, in plume \mu=' num2str(output_IN.Y_mean_overall_ON,'%.2f')]; i=i+1;
 plot(xdat_model_PD_out(1).x,ydat_model_PD_out(1).y,'b--','linewidth',3); leg_str{i}=['Model Volc ON, out of plume \mu=' num2str(output_OUT.Y_mean_overall_ON,'%.2f')]; i=i+1;
 plot(xdat_modis_PD(1).x,ydat_modis_PD(1).y,'r','linewidth',3); leg_str{i}=['MODIS, in plume \mu=' num2str(output_IN.Y_mean_overall_modis,'%.2f')]; i=i+1;
 plot(xdat_modis_PD_out(1).x,ydat_modis_PD_out(1).y,'r--','linewidth',3); leg_str{i}=['MODIS, out of plume \mu=' num2str(output_OUT.Y_mean_overall_modis,'%.2f')]; i=i+1;
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

         
         
 
 
 

