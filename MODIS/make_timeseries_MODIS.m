
% timLAT=[-89.5:1:89.5];
% eval([MODIS_varname '.timeseries_LAT = timLAT;']);
% timLON=0.5;
% eval([MODIS_varname '.timeseries_LON = timLON;']);
% 
% if imr==1
%     eval([MODIS_varname '.timeseries=NaN*ones(length(timLAT),length(timLON),nMOD_av); 
%     W_timeseries=NaN*ones(length(timLAT),length(timLON),nMOD_av); 
%     H_timeseries=NaN*ones(length(timLAT),length(timLON),nMOD_av);    
% end
% 
% 
% 
% iLON=findheight_nearest(MLON,timLON);
% %MLON runs from -179.5 to 179.5, so they represent the mid point between
% %whole numbers e.g. -180 to -179 for the first one - so better to specify
% %values starting in between whole numbers
% 
% for itim=1:length(timLAT)
% 
%     iLAT=findheight_nearest(MLAT,timLAT(itim));
% 
%     if ihisto==0
%         eval([MODIS_varname '.timeseries(itim,imr)=modis_data_read(iLAT,iLON);']);
%     else
%         eval([MODIS_varname '.timeseries(:,:,itim,imr)=modis_data_read(:,:,iLAT,iLON);']);
%     end
% 
% end
% 
% timLON=-179.5;
% eval([MODIS_varname '.timeseries2_LON = timLON;']);
% eval([MODIS_varname '.timeseries2_LAT = timLAT;']);
% 
% iLON=findheight_nearest(MLON,timLON);
% %MLON runs from -179.5 to 179.5, so they represent the mid point between
% %whole numbers e.g. -180 to -179 for the first one - so better to specify
% %values starting in between whole numbers
% 
% for itim=1:length(timLAT)
% 
%     iLAT=findheight_nearest(MLAT,timLAT(itim));
%     if ihisto==0
%         eval([MODIS_varname '.timeseries2(itim,imr)=modis_data_read(iLAT,iLON);']);
%     else
%         eval([MODIS_varname '.timeseries(:,:,itim,imr)=modis_data_read(:,:,iLAT,iLON);']);
%     end
% end

% 
% timLAT=[-60.5 -59.5];
% timLON=[-179.5:1:179.5];
% eval([MODIS_varname '.timeseries3_LON = timLON;']);
% eval([MODIS_varname '.timeseries3_LAT = timLAT;']);
% 
% for ilon=1:length(timLON)
% 
%     iLON=findheight_nearest(MLON,timLON(ilon));
%     %MLON runs from -179.5 to 179.5, so they represent the mid point between
%     %whole numbers e.g. -180 to -179 for the first one - so better to specify
%     %values starting in between whole numbers
% 
%     for ilat=1:length(timLAT)
% 
%         iLAT=findheight_nearest(MLAT,timLAT(ilat));
% 
%         if ihisto==0
%             eval([MODIS_varname '.timeseries3(ilat,ilon,imr)=modis_data_read(iLAT,iLON);']);
%         else
%             eval([MODIS_varname '.timeseries3(:,:,ilat,ilon,imr)=modis_data_read(:,:,iLAT,iLON);']);
%         end
% 
%     end            
% 
% end


%% timLAT and timLON allow us to extract a region I think
timLAT=[-60.5 -59.5];

timLAT=[-89.5:1:89.5];
timLON=[-179.5:1:179.5];
%MLON runs from -179.5 to 179.5, so they represent the mid point between
%whole numbers e.g. -180 to -179 for the first one - so better to specify
%values starting in between whole numbers

%preallocate the memory for more efficient processing and memory handling -
%less likely to run out of memory since when the size of an array is continually expanded
%Matlab might have to copy the array in memory doubling the size it takes
%up - it also slows Matlab down

if imr==1  %if on the first time loop
    
    if ihisto==0
        if length(size(modis_data_read))==3
            eval([MODIS_varname '.timeseries3_prof=NaN*ones(length(timLAT),length(timLON),size(modis_data_read,1),nMOD_av);']);
        else
            eval([MODIS_varname '.timeseries3=NaN*ones(length(timLAT),length(timLON),nMOD_av);']);
        end
    else

        clear Nd_timeseries W_timeseries H_timeseries

        Nd_timeseries.mean=NaN*ones(length(timLAT),length(timLON),nMOD_av);
        Nd_timeseries.std_dev=NaN*ones(length(timLAT),length(timLON),nMOD_av);
%        Nd_timeseries.std_dev_norm=NaN*ones(length(timLAT),length(timLON),nMOD_av);

        W_timeseries.mean=NaN*ones(length(timLAT),length(timLON),nMOD_av);
        W_timeseries.std_dev=NaN*ones(length(timLAT),length(timLON),nMOD_av);
%        W_timeseries.std_dev_norm=NaN*ones(length(timLAT),length(timLON),nMOD_av);

        H_timeseries.mean=NaN*ones(length(timLAT),length(timLON),nMOD_av);
        H_timeseries.std_dev=NaN*ones(length(timLAT),length(timLON),nMOD_av);
%        H_timeseries.std_dev_norm=NaN*ones(length(timLAT),length(timLON),nMOD_av);

    end

end
    
eval([MODIS_varname '.timeseries3_LON = timLON;']);
eval([MODIS_varname '.timeseries3_LAT = timLAT;']);

clear iLON iLAT
for i=1:length(timLAT)
    iLAT(i) = find(MLAT==timLAT(i));
end
for i=1:length(timLON)
    iLON(i) = find(MLON==timLON(i));
end


if ihisto==0
    if length(size(modis_data_read))==3
        eval([MODIS_varname '.timeseries3_prof(iLAT,iLON,:,imr)=permute(modis_data_read(:,iLAT,iLON),[2 3 1]);']);
    else
        eval([MODIS_varname '.timeseries3(iLAT,iLON,imr)=modis_data_read(iLAT,iLON);']);
    end
else
    if strmatch(MODIS_varname,'Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius')==1
        
        [histo_output]=...
         histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',modis_data_read(:,:,iLAT,iLON),Cloud_Top_Temperature_Day_Mean.timeseries3(iLAT,iLON,imr));
     
     ilat_inds = [1:length(iLAT)];
     ilon_inds = [1:length(iLON)];
     
         Nd_timeseries.mean(ilat_inds,ilon_inds,imr)=histo_output.N_histo_mean;
         Nd_timeseries.std_dev(ilat_inds,ilon_inds,imr)=histo_output.N_histo_std;        
%         Nd_timeseries.std_dev_norm(ilat_inds,ilon_inds,imr)=histo_output.N_std_norm;
         
         W_timeseries.mean(ilat_inds,ilon_inds,imr)=histo_output.W_histo_mean;
         W_timeseries.std_dev(ilat_inds,ilon_inds,imr)=histo_output.W_histo_std;        
%         W_timeseries.std_dev_norm(ilat_inds,ilon_inds,imr)=histo_output.W_std_norm;
         
         H_timeseries.mean(ilat_inds,ilon_inds,imr)=histo_output.H_histo_mean;
         H_timeseries.std_dev(ilat_inds,ilon_inds,imr)=histo_output.H_histo_std;        
%         H_timeseries.std_dev_norm(ilat_inds,ilon_inds,imr)=histo_output.H_std_norm;

         LWC_timeseries.mean(ilat_inds,ilon_inds,imr)=histo_output.LWC_histo_mean;
         LWC_timeseries.std_dev(ilat_inds,ilon_inds,imr)=histo_output.LWC_histo_std;        
%         H_timeseries.std_dev_norm(ilat_inds,ilon_inds,imr)=histo_output.H
%         _std_norm;

        

    end
%    eval([MODIS_varname '.timeseries3(:,:,ilat_inds,ilon_inds,imr)=modis_data_read(:,:,iLAT,iLON);']);
end







eval([MODIS_varname '.timeseries_days(imr) = str2num(modis_day_str);']);


