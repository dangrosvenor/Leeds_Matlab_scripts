%Runs the PDF plotting function for various runs (e.g., different emission strenghts). The plotting function runs for 
%volc ON and OFF vs MODIS both in and out of plume - and plots the in and
%out of plume results in one plot.

% Probably shoudl simplify to give the in and out of plume PDFs for each
% run (rather than having to do a pair each time...)

%% 10x reduced emission rate (relative to u-ch765)
um_case_PI = 'u-ch765';  %volcano OFF, orog ON
um_case_PD = 'u-cr138';  %volcano ON, orog ON, 10x reduced emission rate (relative to u-ch765)

end_str = '_10x';
eval_str=['[xdat_model_PI' end_str ',ydat_model_PI' end_str ',xdat_model_PI_out' end_str ',ydat_model_PI_out' end_str ',xdat_model_PD' end_str ',ydat_model_PD' end_str ',xdat_model_PD_out' end_str ','...
    'ydat_model_PD_out' end_str ', xdat_modis_PD' end_str ',ydat_modis_PD' end_str ', xdat_modis_PD_out' end_str ',ydat_modis_PD_out' end_str ']' ...
    '= Hawaii_aod_satellite_model_PDFs_FUNC(UM_base_dir,um_case_PI,um_case_PD);'];
eval(eval_str);

%% End here
return

%% 5x reduced emission rate (relative to u-ch764)
um_case_PI = 'u-ch765';  %volcano OFF, orog ON
um_case_PD = 'u-co295';  %volcano ON, orog ON, 5x reduced emission rate (relative to u-ch764)

end_str = '_5x';
eval_str=['[xdat_model_PI' end_str ',ydat_model_PI' end_str ',xdat_model_PI_out' end_str ',ydat_model_PI_out' end_str ',xdat_model_PD' end_str ',ydat_model_PD' end_str ',xdat_model_PD_out' end_str ','...
    'ydat_model_PD_out' end_str ', xdat_modis_PD' end_str ',ydat_modis_PD' end_str ', xdat_modis_PD_out' end_str ',ydat_modis_PD_out' end_str ']' ...
    '= Hawaii_aod_satellite_model_PDFs_FUNC(UM_base_dir,um_case_PI,um_case_PD);'];
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

         
         
 
 
 






