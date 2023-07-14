function [output,xdat_model_ON,ydat_model_ON,xdat_model_OFF,ydat_model_OFF,xdat_modis,ydat_modis,ylabelstr,logbin_norm_driver,thresh_str]...
    = Hawaii_aod_satellite_model_pdf(matches_aod_ON,matches_so2_ON,matches_aod_OFF,so2_thresh)

i_plot_norm_driver=1; %Whether to normalise the PDF by the total bin counts or use raw bin counts

if exist('matches_so2_ON')
   iso2_plume=1;  
else
   iso2_plume=0; 
end

%so2_thresh=1e-5;

model_aod_all_ON=[];
model_aod_all_OFF=[];
model_so2_all_ON=[];
modis_aod_all=[];
%Including all swaths for now - could separate into individual days
%perhaps, or sets of 2 days to look at the progression over time. Or just
%create a timeseries of the mean.
for isw=1:length(matches_aod_ON.model_match)    
    model_aod_all_ON = cat(2,model_aod_all_ON,matches_aod_ON.model_match{isw});    
    model_aod_all_OFF = cat(2,model_aod_all_OFF,matches_aod_OFF.model_match{isw});    
    modis_aod_all = cat(2,modis_aod_all,matches_aod_ON.modis_match{isw});   
    if iso2_plume==1
        model_so2_all_ON = cat(2,model_so2_all_ON,matches_so2_ON.model_match{isw}); 
    end
end

if iso2_plume==1
    inan = find(model_so2_all_ON<so2_thresh(1) | model_so2_all_ON>=so2_thresh(2));
    model_aod_all_ON(inan) = NaN;
    model_aod_all_OFF(inan) = NaN;
    modis_aod_all(inan) = NaN;
    thresh_str_so2 = [num2str(so2_thresh(1)) '<=SO2<' num2str(so2_thresh(2))];
end


 %Make a 1D PDF of the data from the above
 
 logbin_norm_driver=1; %set to one if using log bins and setting the xscale to log so that the area
 %Under the PDFs in log space will be the no. datapoints.
 
 if logbin_norm_driver==0
     dy = 0.02;
     Ybins_DRIVER = [-dy/2:dy:15];
 else
     dy = 0.05
     Ybins_DRIVER = 10.^[-2:dy:1.5];
 end
 ylabelstr='AOD 550nm'; 
 
 Y_driver = model_aod_all_ON;
 %Y_driver = dat_modis(:);
        
 Hawaii_dNd_1D_PDF
 xdat_model_ON = xdat;
 ydat_model_ON = ydat;
 
 output.Y_mean_overall_ON = Y_mean_overall;
 
 %% Make a 1D PDF of the OFF data
 
 logbin_norm_driver=1; %set to one if using log bins and setting the xscale to log so that the area
 %Under the PDFs in log space will be the no. datapoints.
 
 if logbin_norm_driver==0
     dy = 0.02;
     Ybins_DRIVER = [-dy/2:dy:15];
 else
     dy = 0.05
     Ybins_DRIVER = 10.^[-2:dy:1.5];
 end
 ylabelstr='AOD 550nm'; 
 
 Y_driver = model_aod_all_OFF;
 %Y_driver = dat_modis(:);
        
 Hawaii_dNd_1D_PDF
 xdat_model_OFF = xdat;
 ydat_model_OFF = ydat;
 
 output.Y_mean_overall_OFF = Y_mean_overall;
 

 %% PDF of the modis data
 logbin_norm_driver=1; %set to one if using log bins and setting the xscale to log so that the area
 %Under the PDFs in log space will be the no. datapoints.
 
 if logbin_norm_driver==0
     dy = 0.02;
     Ybins_DRIVER = [-dy/2:dy:15];
 else
     dy = 0.05
     Ybins_DRIVER = 10.^[-2:dy:1.5];
 end
 ylabelstr='AOD 550nm'; 
 
 Y_driver = modis_aod_all;
 %Y_driver = dat_modis(:);

         
 Hawaii_dNd_1D_PDF
 xdat_modis = xdat;
 ydat_modis = ydat;
 
 output.Y_mean_overall_modis = Y_mean_overall;
 
 %% Plot all 3 PDFs
        
 figure('color','w'); hold on
 set(gca,'fontsize',16);
 clear leg_str
 i=1;
 plot(xdat_model_OFF(1).x,ydat_model_OFF(1).y,'k','linewidth',3); leg_str{i}=['Model Volc OFF']; i=i+1;
 plot(xdat_model_ON(1).x,ydat_model_ON(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON']; i=i+1;
 plot(xdat_modis(1).x,ydat_modis(1).y,'r','linewidth',3); leg_str{i}=['MODIS']; i=i+1;
 %title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
 %             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
 
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
 title(thresh_str_so2);
 grid on
 

         
         
 
 
 

