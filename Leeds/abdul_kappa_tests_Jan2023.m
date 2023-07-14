function [output] = Hawaii_Nd_satellite_model_pdf(var_str,matches_Nd_ON,model_so2_ON,matches_Nd_OFF,so2_thresh)

 icoarse=0; %whether to coarse grain the data before plotting on PDF
 %May get reset later...!
 

if exist('model_so2_ON')
   iso2_plume=1;  
else
   iso2_plume=0; 
end

f_thresh = 0.1; % min threshold for the fraction of good data values relative to the total number (within the requested SO2 values).
%I.e. =0.1 means that need at least 10% good data.

%so2_thresh=1e-5;

% model_aod_all_ON=[];
% model_aod_all_OFF=[];
% model_so2_all_ON=[];
% modis_aod_all=[];
% %Including all swaths for now - could separate into individual days
% %perhaps, or sets of 2 days to look at the progression over time. Or just
% %create a timeseries of the mean.
% for isw=1:length(matches_Nd_ON.model_match)    
%     model_aod_all_ON = cat(2,model_aod_all_ON,matches_Nd_ON.model_match{isw});    
%     model_aod_all_OFF = cat(2,model_aod_all_OFF,matches_Nd_OFF.model_match{isw});    
%     modis_aod_all = cat(2,modis_aod_all,matches_Nd_ON.modis_match{isw});   
%     if iso2_plume==1
%         model_so2_all_ON = cat(2,model_so2_all_ON,model_so2_ON.model_match{isw}); 
%     end
% end

if iso2_plume==1
    %Select the times where we had matches with MODIS swaths
    so2_model_match_ON = model_so2_ON(:,:,matches_Nd_ON.it_model_match);    
    siz = size(so2_model_match_ON);
    inan = find(so2_model_match_ON<so2_thresh(1) | so2_model_match_ON>=so2_thresh(2));
    matches_Nd_ON.model_match(inan) = NaN;
    matches_Nd_ON.mask_match(inan) = NaN;
    matches_Nd_ON.modis_match(inan) = NaN;
    
    matches_Nd_OFF.model_match(inan) = NaN; 
    matches_Nd_OFF.mask_match(inan) = NaN;
    matches_Nd_OFF.modis_match(inan) = NaN;
    
    thresh_str_so2 = [num2str(so2_thresh(1)) '<=SO2<' num2str(so2_thresh(2))];
    
    %[ix,iy,it_nan] = ind2sub(siz,inan); %convert to ix,iy,ity indices to use for testing the fraction of NaNs for each time
        
    %Check that there is a certain fraction of MODIS data that has been
    %selected from the sub-sample (e.g., in-plume or out-of-plume).
    %If not then don't compare to the model (set model and MODIS to NaNs).
    %mod_dat = matches_Nd_ON.modis_match;
    %mask_dat = matches_Nd_ON.mask_match;
    %mod_dat(inan)=NaN; %delete the values outside of the area to leave only those inside
    %mask_dat(inan)=NaN;
        
    for it=1:siz(3) %Loop over time
%         i = find(it_nan==it);
%         inan2 = sub2ind(siz,ix(i),iy(i),it_nan(i)); %linear indices of area outside of requested region        
%         mod_dat = matches_Nd_ON.modis_match;
%         mask_dat = matches_Nd_ON.mask_match;
%         mod_dat(inan2)=[]; %delete the values outside of the area to leave only those inside
%         mask_dat(inan2)=[];
        
        Ntot_sub = length(find(matches_Nd_ON.mask_match(:,:,it)==1)); %Actually want this to be the number of MODEL data points that were within the swath.
        Nnan_modis = length(find(isnan(matches_Nd_ON.modis_match(:,:,it))==0 & matches_Nd_ON.mask_match(:,:,it)==1 )); %find the number of non-NaN values        
        if Nnan_modis / Ntot_sub < f_thresh %Require the fraction of good data over total to be larger than a min value
            matches_Nd_ON.modis_match(:,:,it) = NaN;
            matches_Nd_OFF.modis_match(:,:,it) = NaN;
            matches_Nd_ON.model_match(:,:,it) = NaN;
            matches_Nd_OFF.model_match(:,:,it) = NaN;
        end
    end
end


%% Set the bins, PDF settings, etc.

 %Make a 1D PDF of the data from the above

 
 logbin_norm_driver=1; %set to one if using log bins and setting the xscale to log so that the area
 %... under the PDFs in log space will be the no. datapoints.
 %Might get overwritten below
 
 fscale=1;
 switch var_str
     case 'Nd'
         if logbin_norm_driver==0
             dy = 0.02;
             Ybins_DRIVER = [-dy/2:dy:15];
         else
             dy = 0.05;
             Ybins_DRIVER = 10.^[0:dy:3];
         end
         ylabelstr='Nd (cm^{-3})';
         
         fscale=1e-6; %convert to per cc from per m3
         
     case 'CF_tau0pt3'
         icoarse=1; M_coarse_grain=10; N_coarse_grain=10;
         logbin_norm_driver=0; %set to one if using log bins and setting the xscale to log so that the area
 %... under the PDFs in log space will be the no. datapoints.
 
         if logbin_norm_driver==0
             dy = 0.02;
             Ybins_DRIVER = [-dy/2:dy:1+dy];
         else
             dy = 0.05;
             Ybins_DRIVER = 10.^[-2:dy:3];
         end
         ylabelstr='CF';
         
 end
 
 %% Make a 1D PDF of the OFF data
 
 Y_driver = matches_Nd_ON.model_match*fscale; 
 %Y_driver = dat_modis(:);
 
 if icoarse==1
     dat = Y_driver;
     clear Y_driver
     for it=1:size(dat,3)
          Y_driver(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);
     end     
 end
 
 Hawaii_dNd_1D_PDF
 
 if logbin_norm_driver==1
     set(gca,'xscale','log');
 end
        
 
 output.xdat_model_ON = xdat;
 output.ydat_model_ON = ydat;
 output.Y_mean_overall_ON = Y_mean_overall;
 
 %% Make a 1D PDF of the OFF data
 
 
 Y_driver = matches_Nd_OFF.model_match*fscale;
 %Y_driver = dat_modis(:);
 
 if icoarse==1
     dat = Y_driver;
     clear Y_driver
     for it=1:size(dat,3)
          Y_driver(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);
     end     
 end
        
 Hawaii_dNd_1D_PDF
 
  if logbin_norm_driver==1
     set(gca,'xscale','log');
  end
 
 output.xdat_model_OFF = xdat;
 output.ydat_model_OFF = ydat;
 output.Y_mean_overall_OFF = Y_mean_overall;
 

 %% PDF of the modis data
 Y_driver = matches_Nd_ON.modis_match;
 %Y_driver = dat_modis(:);

 if icoarse==1
     dat = Y_driver;
     clear Y_driver
     for it=1:size(dat,3)
          Y_driver(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);
     end     
 end
         
 Hawaii_dNd_1D_PDF
 
  if logbin_norm_driver==1
     set(gca,'xscale','log');
  end
 
 output.xdat_modis = xdat;
 output.ydat_modis = ydat;
 output.Y_mean_overall_modis = Y_mean_overall;
 
 %% Plot all 3 PDFs
        
 figure('color','w'); hold on
 set(gca,'fontsize',16);
 clear leg_str
 i=1;
 plot(output.xdat_model_OFF(1).x,output.ydat_model_OFF(1).y,'k','linewidth',3); leg_str{i}=['Model Volc OFF']; i=i+1;
 plot(output.xdat_model_ON(1).x,output.ydat_model_ON(1).y,'b','linewidth',3); leg_str{i}=['Model Volc ON']; i=i+1;
 plot(output.xdat_modis(1).x,output.ydat_modis(1).y,'r','linewidth',3); leg_str{i}=['MODIS']; i=i+1;
 %title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
 %             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
 
 legend(leg_str);
 xlabel(ylabelstr);
 if logbin_norm_driver==1
     set(gca,'xscale','log');
     ylabelstr_actual = 'df/dlog_{10}(x)';
     
 else
     ylabelstr_actual='df/dx';
 end
 
 ylabel(ylabelstr_actual);
 %set(gca,'xlim',[-5 500]);
 
 set(gca,'xscale','log');
 %set(gca,'yscale','log');
 set(gca,'xlim',[0 1000]);
 title(thresh_str_so2);
 grid on

 
output.xlabelstr = ylabelstr; %is actually the x-axis label in this case
output.logbin_norm_driver = logbin_norm_driver;
output.thresh_str = thresh_str;
output.ylabelstr = ylabelstr_actual;
         
         
 
 
 

