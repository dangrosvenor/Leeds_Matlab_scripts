function Hawaii_Nd_satellite_model_PDFs_PLOT(output_in_plume,output_out_plume,opts)

if exist('opts')
    %Convert all of the variable names in the input structure to actual names
    %for ease of use
    name_struc='opts'; %The name of the structure
    names = eval(['fieldnames(' name_struc ');']);
    for i=1:length(names)
        eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
        eval(eval_str);
    end
    
end

%% Plot PI, PD and MODIS together for in and out of plume
 figure('color','w'); hold on
 set(gca,'fontsize',16);
 clear leg_str
 i=1;
 plot(output_in_plume.xdat_model_OFF(1).x,output_in_plume.ydat_model_OFF(1).y,'k','linewidth',3); leg_str{i}=['Model Volc OFF, in plume \mu=' num2str(output_in_plume.Y_mean_overall_OFF,'%.2f') ]; i=i+1;
 plot(output_out_plume.xdat_model_OFF(1).x,output_out_plume.ydat_model_OFF(1).y,'k--','linewidth',3); leg_str{i}=['Model Volc OFF, out of plume \mu=' num2str(output_out_plume.Y_mean_overall_OFF,'%.2f')]; i=i+1;
 plot(output_in_plume.xdat_model_ON(1).x,output_in_plume.ydat_model_ON(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON, in plume \mu=' num2str(output_in_plume.Y_mean_overall_ON,'%.2f')]; i=i+1;
 plot(output_out_plume.xdat_model_ON(1).x,output_out_plume.ydat_model_ON(1).y,'b--','linewidth',3); leg_str{i}=['Model Volc ON, out of plume \mu=' num2str(output_out_plume.Y_mean_overall_ON,'%.2f')]; i=i+1;
 plot(output_in_plume.xdat_modis(1).x,output_in_plume.ydat_modis(1).y,'r','linewidth',3); leg_str{i}=['MODIS, in plume  \mu=' num2str(output_in_plume.Y_mean_overall_modis,'%.2f')]; i=i+1;
 plot(output_out_plume.xdat_modis(1).x,output_out_plume.ydat_modis(1).y,'r--','linewidth',3); leg_str{i}=['MODIS, out of plume  \mu=' num2str(output_out_plume.Y_mean_overall_modis,'%.2f')]; i=i+1;
 %title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
 %             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
 
 %ylabelstr='AOD 550nm'; 
 %ylabelstr='N_d (cm^{-3})'; 
 legend(leg_str);
 xlabel(output_in_plume.xlabelstr);
 ylabel(output_in_plume.ylabelstr);
 title(output_in_plume.tit_str);
 
 %set(gca,'xlim',[-5 500]);
 
 if exist('i_xscale_log') & i_xscale_log==1
     set(gca,'xscale','log');
 end
 if exist('i_yscale_log') & i_yscale_log==1
     set(gca,'yscale','log');
 end
 if exist('xlims')
     set(gca,'xlim',xlims);
 end
% title(thresh_str);
 grid on


 %% Plot PI, PD and MODIS together for in and out of plume - cumulative
 figure('color','w'); hold on
 set(gca,'fontsize',16);
 clear leg_str
 i=1;
 plot(output_in_plume.xdat_model_OFF_cum(1).x,output_in_plume.ydat_model_OFF_cum(1).y,'k','linewidth',3); leg_str{i}=['Model Volc OFF, in plume \mu=' num2str(output_in_plume.Y_mean_overall_OFF,'%.2f') ]; i=i+1;
 plot(output_out_plume.xdat_model_OFF_cum(1).x,output_out_plume.ydat_model_OFF_cum(1).y,'k--','linewidth',3); leg_str{i}=['Model Volc OFF, out of plume \mu=' num2str(output_out_plume.Y_mean_overall_OFF,'%.2f')]; i=i+1;
 plot(output_in_plume.xdat_model_ON_cum(1).x,output_in_plume.ydat_model_ON_cum(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON, in plume \mu=' num2str(output_in_plume.Y_mean_overall_ON,'%.2f')]; i=i+1;
 plot(output_out_plume.xdat_model_ON_cum(1).x,output_out_plume.ydat_model_ON_cum(1).y,'b--','linewidth',3); leg_str{i}=['Model Volc ON, out of plume \mu=' num2str(output_out_plume.Y_mean_overall_ON,'%.2f')]; i=i+1;
 plot(output_in_plume.xdat_modis_cum(1).x,output_in_plume.ydat_modis_cum(1).y,'r','linewidth',3); leg_str{i}=['MODIS, in plume  \mu=' num2str(output_in_plume.Y_mean_overall_modis,'%.2f')]; i=i+1;
 plot(output_out_plume.xdat_modis_cum(1).x,output_out_plume.ydat_modis_cum(1).y,'r--','linewidth',3); leg_str{i}=['MODIS, out of plume  \mu=' num2str(output_out_plume.Y_mean_overall_modis,'%.2f')]; i=i+1;
 %title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
 %             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
 
 %ylabelstr='AOD 550nm'; 
 %ylabelstr='N_d (cm^{-3})'; 
 legend(leg_str);
 xlabel(output_in_plume.xlabelstr);
 %ylabel(output_in_plume.ylabelstr);
 ylabel('Cumulative frequency');
 
 title(output_in_plume.tit_str);
 
 %set(gca,'xlim',[-5 500]);
 
 if exist('i_xscale_log') & i_xscale_log==1
     set(gca,'xscale','log');
 end
 if exist('i_yscale_log') & i_yscale_log==1
     set(gca,'yscale','log');
 end
 if exist('xlims')
     set(gca,'xlim',xlims);
 end
% title(thresh_str);
 grid on

         
         
 
 
 


         
 
 
 

