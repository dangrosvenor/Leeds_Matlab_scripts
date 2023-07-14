irecalc = 0;
ifind_lats=0;

clear vars_in

mon_thresh=[0.4:0.1:0.9];
Nwindow=[5:5:20];

%Restrict time to these indices
itime_restrict = 121:244; %1st May - 1st Sep




vars_in.ifind_lats=ifind_lats;  %set to one if have requested new lat values
vars_in.gcm_Plat2D_SMOS = gcm_Plat2D_SMOS;
vars_in.gcm_Plon2D_SMOS = gcm_Plon2D_SMOS;
vars_in.ilat2=ilat2;
vars_in.ilon2=ilon2;

lon_region = [66:1:90];
lat_region=[16:1:30];

lon_region = [75:1:88];
lat_region=[16:1:28];

% if irecalc==1
% 
%     clear x timser_meancomp Ncomp
%     for imon=1:length(mon_thresh)
%         for iwin=1:length(Nwindow)
% 
%             mon_thresh_in = mon_thresh(imon);
%             Nwindow_in = Nwindow(iwin);
% 
%             [x{imon,iwin},timser_meancomp{imon,iwin},Ncomp{imon,iwin},ilat2,ilon2]=timeseries_plots_SMOS_func(mon_thresh_in,Nwindow_in,Soil_Moisture.timeseries3(:,:,itime_restrict),vars_in);
% 
%         end
%     end
% 
% end


 if irecalc==1
     clear x timser_meancomp Ncomp
 end



for imon=1:length(mon_thresh)
    mon_thresh_in = mon_thresh(imon);

    titlenam = ['Threshold = ' num2str(mon_thresh_in) ' m^3 m^{-3}'];
    for iwin=1:length(Nwindow)


        Nwindow_in = Nwindow(iwin);

        if irecalc==1
            [x{imon,iwin},timser_meancomp{imon,iwin},Ncomp{imon,iwin},ilat2,ilon2]=timeseries_plots_SMOS_func(mon_thresh_in,Nwindow_in,Soil_Moisture.timeseries3(:,:,itime_restrict),vars_in);
        end
        
        xdat_import(iwin).x = x{imon,iwin};
        ydat_import(iwin).y = timser_meancomp{imon,iwin}; 
        labs_import(iwin).l = ['Nwin=' num2str(Nwindow_in) ', N=' num2str(Ncomp{imon,iwin})];        
        
        
    end
    
    %--- run the file to set up the defaults
    watervap_defaults
    
    xlab = 'Time relative to threshold (days)';
    ylab = 'Composite soil moisture (m^3 m^{-3})';
    
    lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
    
    graph=0; %graph choice in watervap
    % --------- Override flags for watervap --------
    man_choose_water_graph=1;    %for watervap
    waterVapourMay2005
    
    savename = [savedir 'SMOS_composite_timeseries_' titlenam];
    saveas_ps_fig_emf(gcf,savename);
    
    
    
end


%% Do a plot showing all the thresholds, but just one window size
iwin=2; %Set to use just teh N=10 window
Nwindow_in = Nwindow(iwin);

clear xdat_import

for imon=1:4 %length(mon_thresh)
    mon_thresh_in = mon_thresh(imon);
    
        titlenam = ['Window size = ' num2str(Nwindow_in) ' days'];
    
       xdat_import(imon).x = x{imon,iwin};
       ydat_import(imon).y = timser_meancomp{imon,iwin}; 
       labs_import(imon).l = ['Thresh=' num2str(mon_thresh_in) ', N=' num2str(Ncomp{imon,iwin})];          
 
 
end


 man_choose_water_graph=1;    %for watervap
 waterVapourMay2005
 
 savename = [savedir 'SMOS_composite_timeseries_var_thresh_' titlenam];
 saveas_ps_fig_emf(gcf,savename);

