dtau_method = 'constant';
dtau_method = 'function of tauc';

iplot_pdfs=0; %Need to fix nhdhistc for this - or use older Matlab version?
iplot_maps=1;
iplot_seasonal=1;
iplot_std=1;
itable=1;

%Load these two in the first time
ireload_modis=0;
iload_seaice_FILE = 0; %Load in the trimmed data

% Or re-do sea-ice (no need if have saved .mat file already)
iload_seaice_ALL = 0;  %Load all the sea-ice for all years, time match and cut down.


iscreen_seaice=1;%whether to screen sea-ice
iscreen_land=1; %whether to screen land - since it uses CTH I think there is no data over land?

savedir_DRIVER = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/';

clear gca offset_lon_DRIVER offset_lat_DRIVER

i=1;
LATs{i} = [-28 -8]; LONs{i} = [-90 -70]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=2; i2(i)=1; region_name{i2(i)} = 'SE Pacific'; i=i+1; %R04 in Bennartz (2017) (VOCALS)
LATs{i} = [-28 -8]; LONs{i} = [-10 15]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=2; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'SE Atlantic'; i=i+1; % (W. Africa)
LATs{i} = [15 35]; LONs{i} = [-140 -115]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=2; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'California'; i=i+1;% California
LATs{i} = [10 40]; LONs{i} = [105 150]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=0; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'East China Sea'; i=i+1;% China
LATs{i} = [-60 -40]; LONs{i} = [-30 120]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=0; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'Southern Ocean'; i=i+1; % SO
%LATs{i} = [40 60]; LONs{i} = [160 -160]; offset_lon_DRIVER=0; offset_lat_DRIVER(i)=0;i2(i)=i2(i-1)+1; region_name{i2(i)} = 'SE Pacific'; i=i+1; % Alaska
LATs{i} = [40 60]; LONs{i} = [160 181]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=0; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'Bering Sea'; i=i+1; % Alaska (LHS)
LATs{i} = [40 60]; LONs{i} = [-181 -160]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=0; i2(i)=i2(i-1)+0; region_name{i2(i)} = 'Bering Sea'; i=i+1; % Alaska (RHS) - can just combine the PDFs together
LATs{i} = [49 57]; LONs{i} = [-50 -27]; offset_lon_DRIVER(i)=0; offset_lat_DRIVER(i)=0; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'NW Atlantic'; i=i+1; % N. Atlantic
LATs{i} = [72 75]; LONs{i} = [-3 48]; offset_lon_DRIVER(i)=-20; offset_lat_DRIVER(i)=-4; i2(i)=i2(i-1)+1; region_name{i2(i)} = 'Barents Sea'; i=i+1; % Barents Sea

iregions=unique(i2);

%Load the data
if ireload_modis==1
    filename = '/home/disk/eos8/d.grosvenor/mat_files_various/SAVED_ann2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat';
    load(filename,'Cloud_Optical_Thickness_Liquid_Mean','MLAT','MLON','Cloud_Fraction_Liquid',...
        'Cloud_Effective_Radius_Liquid_Mean','Cloud_Effective_Radius_37_Liquid_Mean','Cloud_Top_Temperature_Day_Mean',...
        'modisyear_timeseries3','daynum_timeseries3');

end




%single_hemisphere='single';
single_hemisphere='combine';



if iscreen_land==1  %

    mask_size=2;   tag='_landmask_no_smooth'; %Can just set mask_size to =1 for no smoothing

    landmask_load= load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');

    %landmask = repmat(flipdim(landmask_load.amsre_land_mask,1),[1 1 Nmonths]);
    %landmask = flipdim(landmask_load.amsre_land_mask,1);

    %Smooth out the land mask to make it more aggressive - seemed to be
    %allowing some small bits of land otherwise.
    %Also need to do this for the sea ice
    halo=8;
    fs=fspecial('average',mask_size); %Make the averaging filter
    landmask_tmp = add_halo(landmask_load.amsre_land_mask,halo);
    landmask_smooth = filter2(fs,landmask_tmp,'same');
    landmask_smooth = remove_halo(landmask_smooth,halo);
    %landmask_smooth = repmat(flipdim(landmask_smooth,1),[1 1 Nmonths]);
    landmask_smooth = flipdim(landmask_smooth,1);

    %Land mask is zero where there is land and NaN where not - so can just add to the
    %array to make the non-land regions NaN

    landmask_3d = repmat(landmask_smooth,[1 1 size(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,3)]);

    tau_no_land = Cloud_Optical_Thickness_Liquid_Mean.timeseries3 + landmask_3d;


else
    tau_no_land = Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
end

if iscreen_seaice==1
    if iload_seaice_FILE==1  %Load the pre-processed sea-ice from a file rather than all the daily data
        load('/home/disk/eos8/d.grosvenor/mat_files_various/SAVED_ann2008_sea_ice.mat');
    end


    if iload_seaice_ALL==1
        %% Load and combine NH and SH seaice (is only SH in the file)
        seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2000-2015_all_NHemisphere_20161223T184753.mat';
        seaiceNH = load(seaice_fileload_1deg_NH);
        seaice_fileload_1deg_SH = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2000-2015_all_SHemisphere_20170110T094057.mat';
        seaiceSH = load(seaice_fileload_1deg_SH);

        switch single_hemisphere
            case 'single'
                [seaice_time3,seaice_max_time3] = seaice_match_times_FUNC(seaiceSH,Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3);
            case 'combine'
                %Script to combine the two together - not needed for SO case
                sea_ice_combine_NH_SH
        end

    end



    %find points where there is any sea-ice using the threshold below
    sea_ice_thresh = 1e-5;
    %    sea_ice_thresh = 1e9;
    isea_high = find(seaice_max_time3>sea_ice_thresh); %max over 2 week window
    %isea_high = find(seaice_time3>sea_ice_thresh); %daily

    % Screen the data as required
    %Nd_16(isea_high)=NaN;
    %Nd_21(isea_high)=NaN;
    %Nd_37(isea_high)=NaN;

    tau_no_land(isea_high)=NaN;

end

%% Nd and error calcs

% To get Nd think would have to re-do one of the processing stages;
% or use the later files as given to Michael Diamond - but think
% only have 1km Nd?
[N21,H21,LWP21]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,1e-6*Cloud_Effective_Radius_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);
[N37,H37,LWP37]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,1e-6*Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);



switch dtau_method
    case 'constant'

        %Estimate with dtau removed :-
        dtau21=3.3;%the penetration depth from cloud top that corresponds to the retrieved re, allowing for penetration depth.
        %Fig. 4 of Platnick (2000) sugggests that this is around 3.5 for 2.1um
        dtau37=2.0; %for 3.7um
        [N21_dtau]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau21,1e-6*Cloud_Effective_Radius_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);
        [N37_dtau]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau37,1e-6*Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);

        %limit corrected tau value to a minimum to avoid divide very large
        %errors as tau-dtau approaches zero
        min_tau = 1.0;
        N21_dtau( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau21 < min_tau ) = NaN;
        N37_dtau( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau37 < min_tau ) = NaN;
        N21_2 = N21;
        N21_2( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau21 < min_tau ) = NaN;
        N37_2 = N37;
        N37_2( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau37 < min_tau ) = NaN;

    case 'function of tauc'

        
        tau_extrap_method = 'constant dtau at tau>15';
%        tau_extrap_method = 'zero dtau at tau>15';
        tau_extrap_method = 'Mean curves from Odran';
        tau_extrap_method = 'Mean curves from Odran reff ratio';        
        
        %Functions to interpolate
        tauc = [999 15    10     8     5];
       
        switch tau_extrap_method
            case 'constant dtau at tau>15'
                tau_star21 =[4.6119   4.6119   4.1375   3.7446   2.5126 ];
                tau_star37 = [0.7428 0.7428    1.8654    1.8328    1.8589];

            case 'zero dtau at tau>15'
                tau_star21 =[0   4.6119   4.1375   3.7446   2.5126 ];
                tau_star37 = [0 0.7428    1.8654    1.8328    1.8589];
                
            case 'Mean curves from Odran'
                Y_mean_file = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved.mat';
                load(Y_mean_file);
                tauc = X_mean_all_Nd_21_mum; %2.1 and 3.7um tau values (x-axis values) are the same
                tau_star21 = Y_mean_all_Nd_21_mum;
                tau_star37 = Y_mean_all_Nd_37_mum;      
                tauc(end+1)=999;
                tau_star21(end+1) = tau_star21(end);
                tau_star37(end+1) = tau_star37(end); 
                
            case 'Mean curves from Odran reff ratio'
                Y_mean_file = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved_reff_param.mat';
                load(Y_mean_file);
                tauc = X_mean_all_Nd_21_mum; %2.1 and 3.7um tau values (x-axis values) are the same
                reff_ratio21 = Y_mean_all_Nd_21_mum;
                reff_ratio37 = Y_mean_all_Nd_37_mum;      
                tauc(end+1)=999;
                reff_ratio21(end+1) = reff_ratio21(end);
                reff_ratio37(end+1) = reff_ratio37(end); 
        end
        
        switch tau_extrap_method
            case 'Mean curves from Odran reff ratio'
                
                g_re21 = interp1(tauc,reff_ratio21,Cloud_Optical_Thickness_Liquid_Mean.timeseries3(:));
                g_re21 = reshape(g_re21,size(Cloud_Optical_Thickness_Liquid_Mean.timeseries3));
                
                g_re37 = interp1(tauc,reff_ratio37,Cloud_Optical_Thickness_Liquid_Mean.timeseries3(:));
                g_re37 = reshape(g_re37,size(Cloud_Optical_Thickness_Liquid_Mean.timeseries3));            
                
                [N21_dtau,H21_dtau,LWP21_dtau]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,1e-6*g_re21.*Cloud_Effective_Radius_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);
                [N37_dtau,H37_dtau,LWP37_dtau]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,1e-6*g_re37.*Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);

                itau = find( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 < 5);
                N21_dtau(itau) = NaN; LWP21_dtau(itau) = NaN;
                N37_dtau(itau) = NaN; LWP37_dtau(itau) = NaN;
                

            otherwise


                dtau21 = interp1(tauc,tau_star21,Cloud_Optical_Thickness_Liquid_Mean.timeseries3(:));
                dtau21 = reshape(dtau21,size(Cloud_Optical_Thickness_Liquid_Mean.timeseries3));

                dtau37 = interp1(tauc,tau_star37,Cloud_Optical_Thickness_Liquid_Mean.timeseries3(:));
                dtau37 = reshape(dtau37,size(Cloud_Optical_Thickness_Liquid_Mean.timeseries3));

                [N21_dtau]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau21,1e-6*Cloud_Effective_Radius_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);
                [N37_dtau]=MODIS_N_H_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau37,1e-6*Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,'calc',NaN,Cloud_Top_Temperature_Day_Mean.timeseries3);

        end
        
        %Already filtering for tau>5 in the above since don't have info for
        %tau_star for tau<5
        
%        N21_dtau( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau21 < min_tau ) = NaN;
%        N37_dtau( Cloud_Optical_Thickness_Liquid_Mean.timeseries3 - dtau37 < min_tau ) = NaN;
        N21_2 = N21;
        inan = find(isnan(N21_dtau)==1);
        N21_2( inan ) = NaN;
        
        N37_2 = N37;
        inan = find(isnan(N21_dtau)==1);
        N37_2( inan ) = NaN;

        LWP21_2 = LWP21;
        inan = find(isnan(LWP21_dtau)==1);
        LWP21_2( inan ) = NaN;
        
        LWP37_2 = LWP37;
        inan = find(isnan(LWP37_dtau)==1);
        LWP37_2( inan ) = NaN;       
        
        
end


ifilter_tau_as_Nd = 1; %e.g. to remove tau<5 if also do for Nd
if ifilter_tau_as_Nd==1
    iNd_nan = find(isnan(N21_dtau)==1);
    tau_no_land(iNd_nan) = NaN;
    
    iNd_nan37 = find(isnan(N37_dtau)==1);
    tau_no_land(iNd_nan37) = NaN;
end


%Have just screened tau so far - now screen Nd too
inan=find(isnan(tau_no_land)==1);

N21(inan) = NaN;
N37(inan) = NaN;
N21_dtau(inan) = NaN;
N37_dtau(inan) = NaN;
N21_2(inan) = NaN;
N37_2(inan) = NaN;

LWP21(inan) = NaN;
LWP37(inan) = NaN;
LWP21_dtau(inan) = NaN;
LWP37_dtau(inan) = NaN;
LWP21_2(inan) = NaN;
LWP37_2(inan) = NaN;

if iplot_maps==1

%% Make a map of the number of successful retrievals (indicating Sc regions)
    %Useful for getting an idea of where most of the data (and hence Sc
    %regions are).
    
    caxis_range = [0 150];
    
    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure

    hs01=subplot(xsub,ysub,1);
    titlenam_driver = ['Number of days in dataset'];
    
%    [me,nnums] = meanNoNan(tau_no_land,3);
    [me,nnums] = meanNoNan(N21_2,3);    
    dat_modis = nnums;

%    MODIS_vert_pen_tau_DRIVER_template_global_map
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    region_mean_ndays = mean_reg_val;
    region_max_ndays = max_reg_val;
    
    %Change the colormap to blue to brown :-
    hcbar = find_peer_colorbars_of_an_axes(gca);
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));    
    set(hcbar,'XaxisLocation','bottom');
    
    map_fig = gcf;
            
    %Plot the boxes on the map
    for iLATs=1:length(LATs)
        LAT_val_DRIVER = LATs{iLATs}; LON_val_DRIVER = LONs{iLATs};


        ioverride_box_colour=1;
        itext_in_box = 1;
        irotated_pole_box=0;
        imap=1;
        box_text=num2str(i2(iLATs));
        LAT_val = LAT_val_DRIVER; LON_val = LON_val_DRIVER; %12th Nov UM FULL domain
        col_str='k-';  %
        offset_lon = offset_lon_DRIVER(iLATs);
        offset_lat = offset_lat_DRIVER(iLATs);
        box_lwidth = 1;
        figure(map_fig);
        plot_box_on_map
    end
    
    
    caxis(caxis_range);
%    gca = hs01;
    axes(hs01);
    increase_font_size_map_figures
    
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    hcbar = Hc1s;        
    pos2 = get(Hc1s,'position'); %[left bottom width height]
    pos_new = pos2; %from increase_font_size_map_figures
    pos_new(2)=0.25; %pos_new(3)=0.58;
    set(Hc1s,'position',pos_new);

% ------------------------------------------------------------
% Make the 2nd subplot panel a map of the mean tau
% or could od fraction of days with tau>5?
    %subplot of mean tau
    
    
    caxis_range=[0 25];
    
    hs02=subplot(xsub,ysub,2);
    titlenam_driver = ['Mean optical depth'];

    

    %Since will likely want to restrict to tau>5
    [me,nnums,stdev] = meanNoNan(tau_no_land,3);

    
    dat_modis = me;
    

    %MODIS_vert_pen_tau_DRIVER_template_global_map
    MODIS_vert_pen_NdError_DRIVER_template_global_map
%    region_mean_tau = me_reg_val; %The mean weighted by the number of time values for each column
%    region_std_tau = std_reg_val; %This is the overall std over both time and space for the region (i.e. weighted as above)
    region_mean_tau = mean_reg_val; %The straight mean of the time-averaged tau values for each column (each column has equal weighting)
    region_std_tau = mean_std_equal_weighted_reg_val; %This is a straight mean of the std over time (each column with equal weighting)
    
    caxis(caxis_range);
%    gca = hs02;
    axes(hs02);
    increase_font_size_map_figures
    
    %move colorbar down a bit and expand to full width
    Hc2s = find_peer_colorbars_of_an_axes(gca);
    hcbar = Hc2s;        
    pos2 = get(Hc2s,'position'); %[left bottom width height]
    pos_new = pos2; %from increase_font_size_map_figures
    pos_new(2)=0.25; pos_new(1)=pos2(1)-0.1203;
    set(Hc2s,'position',pos_new);
    
    %adjust the position of the second plot to reduce space between the
    %subplots
    pos=get(hs02,'position');
    pos_new=pos;
    pos_new(1)=0.45;
    set(hs02,'position',pos_new);
    
    %Change the colormap to blue to brown :-
    hcbar = find_peer_colorbars_of_an_axes(gca);
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));
    set(hcbar,'XaxisLocation','bottom');
    
    %    saveas_ps_fig_emf(gcf,[savename],'',0,1);
    
    
    
%% Figure for fraction of days with tau<thresh_tau (not shown in paper, but used in the table)
  %Just used to get frac values for tables.
    thresh_tau = 10.8;
    thresh_tau = 10.0;    
    itau_lt_5 = find(tau_no_land <= thresh_tau);
    tau_gt_5 = tau_no_land; tau_gt_5(itau_lt_5)=NaN;
    [me2,nnums2] = meanNoNan(tau_gt_5,3);
 
    %If doing fraction of days with tau>5
%    dat_modis = nnums2./nnums;
%    dat_modis(nnums<15)=NaN;
    %If doing fraction of days with tau>5
    dat_modis = nnums2./nnums;
    dat_modis(nnums<15)=NaN;
    MODIS_vert_pen_tau_DRIVER_template_global_map
%    MODIS_vert_pen_NdError_DRIVER_template_global_map
    region_mean_fdays = mean_reg_val;
    region_max_fdays = max_reg_val;
    

%% Make a map of the mean error in Nd using certain dtau values
    % For N21_dtau the tau is limited to min_tau(=1.0) - set to NaN if below
    % this - prevents zero taus from being used.
    % N21_2 also has this applied N21 does not.

    %% N.B. Expand to full window size before saving!

    caxis_range = [0 45];

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure

    hs01=subplot(xsub,ysub,1);

    % Time mean of Nd divided by mean of corrected values
%    dat_modis = 100*(meanNoNan(N21_2,3)./meanNoNan(N21_dtau,3) - 1);
    % Or time mean of the individual relative errors
    dat_modis = 100*meanNoNan( N21_2./N21_dtau - 1 , 3 );    

    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product
    titlenam_driver = ['% error for 2.1\mum N_d'];    %So, we overestimate Nd due to this error

    %call script - also calcuulates overall mean for the different regions
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    
    %save mean rel errors for each region    
    rel_err_region_mean_21 = mean_reg_val; %straight mean of the individual location relative errors
    %over the different regions
    
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    hcbar = Hc1s;
    
    
    pos2 = get(Hc1s,'position'); %[left bottom width height]
    pos_new = pos2; %from increase_font_size_map_figures
    pos_new(2)=0.25; pos_new(3)=0.58;
    set(Hc1s,'position',pos_new);

    hs02=subplot(xsub,ysub,2);
    %dat_modis = 100*(meanNoNan(N37_2,3)./meanNoNan(N37_dtau,3) - 1);
    dat_modis = 100*meanNoNan(N37_2 ./ N37_dtau - 1 , 3);
    titlenam_driver = ['% error for 3.7\mum N_d'];
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    
    %save mean rel errors for each region
    rel_err_region_mean_37 = mean_reg_val; %straight mean of the individual location relative errors
    %over the different regions
    
    caxis(caxis_range);
    gca = hs02;
    increase_font_size_map_figures
    %move colorbar down a bit
    delete(Hc1s)  %just need one big colorbar
    % Hc1s = find_peer_colorbars_of_an_axes(gca);
    % pos2 = get(Hc1s,'position');
    % pos_new = pos2;
    % pos_new(2)=0.25;
    % set(Hc1s,'position',pos_new);

    %adjust the position of the second plot to reduce space between the
    %subplots
    pos=get(hs02,'position');
    pos_new=pos;
    pos_new(1)=0.45;
    set(hs02,'position',pos_new);


    %Change the colormap to blue to brown :-
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));    
    set(hcbar,'XaxisLocation','bottom');

%% LWP - make a map of the mean error in LWP using certain dtau values
    % For N21_dtau the tau is limited to min_tau(=1.0) - set to NaN if below
    % this - prevents zero taus from being used.
    % N21_2 also has this applied N21 does not.

    %% N.B. Expand to full window size before saving!

    caxis_range = [-15 0];

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure

    hs01=subplot(xsub,ysub,1);

    % Time mean of Nd divided by mean of corrected values
%    dat_modis = 100*(meanNoNan(N21_2,3)./meanNoNan(N21_dtau,3) - 1);
    % Or time mean of the individual relative errors
    dat_modis = 100*meanNoNan( LWP21_2./LWP21_dtau - 1 , 3 );    

    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product
    titlenam_driver = ['% error for 2.1\mum LWP'];    %So, we overestimate Nd due to this error

    %call script - also calcuulates overall mean for the different regions
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    
    %save mean rel errors for each region    
    rel_err_region_mean_LWP21 = mean_reg_val; %straight mean of the individual location relative errors
    %over the different regions
    
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    hcbar = Hc1s;
    
    
    pos2 = get(Hc1s,'position'); %[left bottom width height]
    pos_new = pos2; %from increase_font_size_map_figures
    pos_new(2)=0.25; pos_new(3)=0.58;
    set(Hc1s,'position',pos_new);

    hs02=subplot(xsub,ysub,2);
    %dat_modis = 100*(meanNoNan(N37_2,3)./meanNoNan(N37_dtau,3) - 1);
    dat_modis = 100*meanNoNan(LWP37_2 ./ LWP37_dtau - 1 , 3);
    titlenam_driver = ['% error for 3.7\mum LWP'];
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    
    %save mean rel errors for each region
    rel_err_region_mean_LWP37 = mean_reg_val; %straight mean of the individual location relative errors
    %over the different regions
    
    caxis(caxis_range);
    gca = hs02;
    increase_font_size_map_figures
    %move colorbar down a bit
    delete(Hc1s)  %just need one big colorbar
    % Hc1s = find_peer_colorbars_of_an_axes(gca);
    % pos2 = get(Hc1s,'position');
    % pos_new = pos2;
    % pos_new(2)=0.25;
    % set(Hc1s,'position',pos_new);

    %adjust the position of the second plot to reduce space between the
    %subplots
    pos=get(hs02,'position');
    pos_new=pos;
    pos_new(1)=0.45;
    set(hs02,'position',pos_new);


    %Change the colormap to blue to brown :-
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));    
    set(hcbar,'XaxisLocation','bottom');
   
    
    
    
    
%% Keep this block statement to allow execution of the above block (due to
% if end statement)

end



%% Make a map of the std dev of the error in Nd using certain dtau values
    % For N21_dtau the tau is limited to min_tau(=1.0) - set to NaN if below
    % this - prevents zero taus from being used.
    % N21_2 also has this applied N21 does not.

    %% N.B. Expand to full window size before saving!
    
if iplot_std==1    

%% Nd std dev
    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure

    hs01=subplot(xsub,ysub,1);

    % Time mean of Nd divided by mean of corrected valeus
    %error_abs = ( N21_2 - N21_dtau );
    %[me,nnums,stdev]=meanNoNan(error_abs,3);
    %dat_modis = 100 * stdev ./ me;  %./ meanNoNan(N21_dtau,3);  %100*(meanNoNan(N21_2,3)./meanNoNan(N21_dtau,3) - 1);
    
    error_rel = 100*(N21_2./N21_dtau - 1);
    [me,nnums,stdev]=meanNoNan(error_rel,3);
    
    %dat_modis = stdev; % normal std dev of percentage errors.
    %titlenam_driver = ['Std dev (%) of percentage errors for 2.1\mum N_d'];    %So, we overestimate Nd due to this error
    %caxis_range = [0 30];
    
    dat_modis = 100 * stdev ./ me; %relative std devs of percentage errors
    titlenam_driver = ['Relative std dev (%) of percentage errors for 2.1\mum N_d'];    %So, we overestimate Nd due to this error        
    caxis_range = [0 80];
    
    mean_error_map_21 = dat_modis;


    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product

    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    hcbar = Hc1s;
    pos2 = get(Hc1s,'position'); %[left bottom width height]
    pos_new = pos2; %from increase_font_size_map_figures
    pos_new(2)=0.25; pos_new(3)=0.58;
    set(Hc1s,'position',pos_new);
    
    %save the std devs.
    std_region_mean_21 = mean_reg_val; %straight mean of the individual location relative std devs
    %over the different regions
    std_region_max_21 = max_reg_val; %straight mean of the individual location relative std devs    
    std_region_min_21 = min_reg_val; %straight mean of the individual location relative std devs            -


    %2nd plot - 3.7um
    hs02=subplot(xsub,ysub,2);

%    error_abs = ( N37_2 - N37_dtau );
%    [me,nnums,stdev]=meanNoNan(error_abs,3);    
%    dat_modis = 100 * stdev ./ me;  %100*(meanNoNan(N21_2,3)./meanNoNan(N21_dtau,3) - 1);

    error_rel = 100*(N37_2./N37_dtau - 1);
    [me,nnums,stdev]=meanNoNan(error_rel,3);
%    dat_modis = stdev; %already in %   % normal std dev of percentage errors.
%    titlenam_driver = ['Std dev (%) of percentage errors for 3.7\mum N_d'];
%    caxis_range = [0 30];

    dat_modis = 100 * stdev ./ me; %relative std devs of percentage errors
    titlenam_driver = ['Relative std dev (%) of percentage errors for 3.7\mum N_d'];
    
    mean_error_map_37 = dat_modis;
    
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs02;
    increase_font_size_map_figures
    %move colorbar down a bit
    delete(Hc1s)  %just need one big colorbar
    % Hc1s = find_peer_colorbars_of_an_axes(gca);
    % pos2 = get(Hc1s,'position');
    % pos_new = pos2;
    % pos_new(2)=0.25;
    % set(Hc1s,'position',pos_new);

    %save std devs
    std_region_mean_37 = mean_reg_val; %straight mean of the individual location relative std devs
    std_region_max_37 = max_reg_val; %straight mean of the individual location relative std devs    
    std_region_min_37 = min_reg_val; %straight mean of the individual location relative std devs        
    %over the different regions


    %adjust the position of the second plot to reduce space between the
    %subplots
    pos=get(hs02,'position');
    pos_new=pos;
    pos_new(1)=0.45;
    set(hs02,'position',pos_new);

    %Change the colormap to blue to brown :-
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));
    
    set(hcbar,'XaxisLocation','bottom');

%% LWP std dev
    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure

    hs01=subplot(xsub,ysub,1);

    % Time mean of Nd divided by mean of corrected valeus
    %error_abs = ( N21_2 - N21_dtau );
    %[me,nnums,stdev]=meanNoNan(error_abs,3);
    %dat_modis = 100 * stdev ./ me;  %./ meanNoNan(N21_dtau,3);  %100*(meanNoNan(N21_2,3)./meanNoNan(N21_dtau,3) - 1);
    
    error_rel = 100*(LWP21_2./LWP21_dtau - 1);
    [me,nnums,stdev]=meanNoNan(error_rel,3);
    
    %dat_modis = stdev; % normal std dev of percentage errors.
    %titlenam_driver = ['Std dev (%) of percentage errors for 2.1\mum N_d'];    %So, we overestimate Nd due to this error
    %caxis_range = [0 30];
    
    dat_modis = 100 * stdev ./ abs(me); %relative std devs of percentage errors
    titlenam_driver = ['Relative std dev (%) of percentage errors for 2.1\mum N_d'];    %So, we overestimate Nd due to this error        
    caxis_range = [0 65];
    
    mean_error_map_LWP21 = dat_modis;


    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product

    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    hcbar = Hc1s;
    pos2 = get(Hc1s,'position'); %[left bottom width height]
    pos_new = pos2; %from increase_font_size_map_figures
    pos_new(2)=0.25; pos_new(3)=0.58;
    set(Hc1s,'position',pos_new);
    
    %save the std devs.
    std_region_mean_LWP21 = mean_reg_val; %straight mean of the individual location relative std devs
    %over the different regions
    std_region_max_LWP21 = max_reg_val; %straight mean of the individual location relative std devs    
    std_region_min_LWP21 = min_reg_val; %straight mean of the individual location relative std devs            -


    %2nd plot - 3.7um
    hs02=subplot(xsub,ysub,2);

%    error_abs = ( N37_2 - N37_dtau );
%    [me,nnums,stdev]=meanNoNan(error_abs,3);    
%    dat_modis = 100 * stdev ./ me;  %100*(meanNoNan(N21_2,3)./meanNoNan(N21_dtau,3) - 1);

    error_rel = 100*(LWP37_2./LWP37_dtau - 1);
    [me,nnums,stdev]=meanNoNan(error_rel,3);
%    dat_modis = stdev; %already in %   % normal std dev of percentage errors.
%    titlenam_driver = ['Std dev (%) of percentage errors for 3.7\mum N_d'];
%    caxis_range = [0 30];

    dat_modis = 100 * stdev ./ abs(me); %relative std devs of percentage errors
    titlenam_driver = ['Relative std dev (%) of percentage errors for 3.7\mum N_d'];
    
    mean_error_map_LWP37 = dat_modis;
    
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs02;
    increase_font_size_map_figures
    %move colorbar down a bit
    delete(Hc1s)  %just need one big colorbar
    % Hc1s = find_peer_colorbars_of_an_axes(gca);
    % pos2 = get(Hc1s,'position');
    % pos_new = pos2;
    % pos_new(2)=0.25;
    % set(Hc1s,'position',pos_new);

    %save std devs
    std_region_mean_LWP37 = mean_reg_val; %straight mean of the individual location relative std devs
    std_region_max_LWP37 = max_reg_val; %straight mean of the individual location relative std devs    
    std_region_min_LWP37 = min_reg_val; %straight mean of the individual location relative std devs        
    %over the different regions


    %adjust the position of the second plot to reduce space between the
    %subplots
    pos=get(hs02,'position');
    pos_new=pos;
    pos_new(1)=0.45;
    set(hs02,'position',pos_new);

    %Change the colormap to blue to brown :-
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));
    
    set(hcbar,'XaxisLocation','bottom');    
    
%% Keep to avoid having to use the if statement if executing block
    
    
end



%% Make a table of the mean and std of relative errors for the different
%% regions. Original table with tau stats and Nd errors
% Make the prc biases table
% 6 columns wide - # & Region name & Mean 2.1 & Std dev 2.1 & Mean 3.7 & Std dev 3.7 

if itable==1

for i=1:length(rel_err_region_mean_21)
    fprintf(1,'%d',i);
    fprintf(1,' & %s ',region_name{i});
    fprintf(1,' & %2.1f & %2.0f',region_mean_ndays(i),region_max_ndays(i)); 
    fprintf(1,' & %.1f & %.2f & %.2f &',region_mean_tau(i),region_std_tau(i),1-region_mean_fdays(i));     
    fprintf(1,' & %2.1f', rel_err_region_mean_21(i));
    fprintf(1,' & %2.1f &', std_region_mean_21(i));
    fprintf(1,' & %2.1f', rel_err_region_mean_37(i));
    fprintf(1,' & %2.1f', std_region_mean_37(i));       
    fprintf(1,'%s\n',' \\');
end

rel_std_21 = std_region_mean_21./rel_err_region_mean_21;
rel_std_37 = std_region_mean_37./rel_err_region_mean_37;


end


%% Make a table of the mean and std of relative errors for the different regions. First new table - just tau stats 
% Make the prc biases table
% 6 columns wide - # & Region name & Mean 2.1 & Std dev 2.1 & Mean 3.7 & Std dev 3.7 

if itable==1

for i=1:length(rel_err_region_mean_21)        
    fprintf(1,'%d',i);
    fprintf(1,' & %s ',region_name{i});
    
    fprintf(1,' & %2.1f & %2.0f',region_mean_ndays(i),region_max_ndays(i)); 
    fprintf(1,' & %.1f & %.2f & %.2f',region_mean_tau(i),region_std_tau(i),1-region_mean_fdays(i));     
    fprintf(1,'%s\n',' \\');
end

rel_std_21 = std_region_mean_21./rel_err_region_mean_21;
rel_std_37 = std_region_mean_37./rel_err_region_mean_37;


end

%% Make a table of the mean and std of relative errors for the different
%  regions. Second new table - just Nd and LWP errors and std devs

% Make the prc biases table
% 10 columns wide - # & Region name & Mean 2.1 & Std dev 2.1 & Mean 3.7 & Std dev 3.7 & Mean 2.1 & Std dev 2.1 & Mean 3.7 & Std dev 3.7 

if itable==1

for i=1:length(rel_err_region_mean_21)
    fprintf(1,'%d',i);
    fprintf(1,' & %s &',region_name{i});

    fprintf(1,' & %2.1f', rel_err_region_mean_21(i));
    fprintf(1,' & %2.1f &', std_region_mean_21(i));
    fprintf(1,' & %2.1f', rel_err_region_mean_37(i));
    fprintf(1,' & %2.1f &', std_region_mean_37(i));       
    
    fprintf(1,' & %2.1f', rel_err_region_mean_LWP21(i));
    fprintf(1,' & %2.1f &', std_region_mean_LWP21(i));
    fprintf(1,' & %2.1f', rel_err_region_mean_LWP37(i));
    fprintf(1,' & %2.1f', std_region_mean_LWP37(i));           
    
    fprintf(1,'%s\n',' \\');
end

rel_std_21 = std_region_mean_21./rel_err_region_mean_21;
rel_std_37 = std_region_mean_37./rel_err_region_mean_37;


end




%% Seasonal map of the Nd error for 2.1um only (can probably infer the
    %% 3.7um map from this?)
    %% Make a map of the mean error in Nd using certain dtau values
    % For N21_dtau the tau is limited to min_tau(=1.0) - set to NaN if below
    % this - prevents zero taus from being used.
    % N21_2 also has this applied N21 does not.

    %% N.B. Expand to full window size before saving!
    
if iplot_seasonal==1

    days=load(filename,'daynum_timeseries3');

    caxis_range = [0 60];
    caxis_range = [0 50];    

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=2; ysub=2; %no. rows, columns
    figure

    %% ----  DJF ---------------------------
    hs01=subplot(xsub,ysub,1);
    gca = hs01;
    season_str='DJF';
    iseason = time_match_season_days(days.daynum_timeseries3,season_str);
    N21_2_season = meanNoNan(N21_2(:,:,iseason),3);
    N21_dtau_season = meanNoNan(N21_dtau(:,:,iseason),3);

    % Time mean of Nd divided by mean of corrected valeus
    dat_modis = 100*( N21_2_season./N21_dtau_season - 1 );

    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product
    titlenam_driver = ['% error for 2.1\mum ' season_str];    %So, we overestimate Nd due to this error
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    %pos2 = get(Hc1s,'position'); %[left bottom width height]
    %pos_new = pos2; %from increase_font_size_map_figures
    %pos_new(2)=0.25; pos_new(3)=0.58;
    %set(Hc1s,'position',pos_new);
    delete(Hc1s)  %just need one big colorbar


    %% ----  MAM ---------------------------
    hs02=subplot(xsub,ysub,2);
    gca = hs02;
    season_str='MAM';
    iseason = time_match_season_days(days.daynum_timeseries3,season_str);
    N21_2_season = meanNoNan(N21_2(:,:,iseason),3);
    N21_dtau_season = meanNoNan(N21_dtau(:,:,iseason),3);

    % Time mean of Nd divided by mean of corrected valeus
    dat_modis = 100*( N21_2_season./N21_dtau_season - 1 );

    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product
    titlenam_driver = ['% error for 2.1\mum ' season_str];    %So, we overestimate Nd due to this error
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs02;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    %pos2 = get(Hc1s,'position'); %[left bottom width height]
    %pos_new = pos2; %from increase_font_size_map_figures
    %pos_new(2)=0.25; pos_new(3)=0.58;
    %set(Hc1s,'position',pos_new);
    delete(Hc1s)  %just need one big colorbar



    %% ----  JJA ---------------------------
    hs03=subplot(xsub,ysub,3);
    gca = hs03;
    season_str='JJA';
    iseason = time_match_season_days(days.daynum_timeseries3,season_str);
    N21_2_season = meanNoNan(N21_2(:,:,iseason),3);
    N21_dtau_season = meanNoNan(N21_dtau(:,:,iseason),3);

    % Time mean of Nd divided by mean of corrected valeus
    dat_modis = 100*( N21_2_season./N21_dtau_season - 1 );

    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product
    titlenam_driver = ['% error for 2.1\mum ' season_str];    %So, we overestimate Nd due to this error
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs03;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    %pos2 = get(Hc1s,'position'); %[left bottom width height]
    %pos_new = pos2; %from increase_font_size_map_figures
    %pos_new(2)=0.25; pos_new(3)=0.58;
    %set(Hc1s,'position',pos_new);
    delete(Hc1s)  %just need one big colorbar


    %% ----  SON ---------------------------
    hs04=subplot(xsub,ysub,4);
    gca = hs04;
    season_str='SON';
    iseason = time_match_season_days(days.daynum_timeseries3,season_str);
    N21_2_season = meanNoNan(N21_2(:,:,iseason),3);
    N21_dtau_season = meanNoNan(N21_dtau(:,:,iseason),3);

    % Time mean of Nd divided by mean of corrected valeus
    dat_modis = 100*( N21_2_season./N21_dtau_season - 1 );

    % Time mean of relative values (probably better to use time mean of
    % absolute values?):-
    %dat_modis = 100*meanNoNan(N21_2./N21_dtau - 1 , 3); %N.B. N21_dtau is the corrected (i.e. "true" value). N21 is the standard biased product
    titlenam_driver = ['% error for 2.1\mum ' season_str];    %So, we overestimate Nd due to this error
    MODIS_vert_pen_NdError_DRIVER_template_global_map
    caxis(caxis_range);
    gca = hs04;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    Hc1s = find_peer_colorbars_of_an_axes(gca);
    %pos2 = get(Hc1s,'position'); %[left bottom width height]
    %pos_new = pos2; %from increase_font_size_map_figures
    %pos_new(2)=0.25; pos_new(3)=0.58;
    %set(Hc1s,'position',pos_new);
    %delete(Hc1s)  %just need one big colorbar - keep the last one





    %make them all bigger - seem to have to increase the height rather than the
    %width for this
    new_height=0.35;
    gca=hs01; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
    gca=hs02; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
    gca=hs03; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
    gca=hs04; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);

    % Move them all down
    dZ=-0.12;
    gca=hs01; pos=get(gca,'position'); pos(2)=pos(2)+dZ;set(gca,'position',pos);
    gca=hs02; pos=get(gca,'position'); pos(2)=pos(2)+dZ;set(gca,'position',pos);
    gca=hs03; pos=get(gca,'position'); pos(2)=pos(2)+dZ;set(gca,'position',pos);
    gca=hs04; pos=get(gca,'position'); pos(2)=pos(2)+dZ;set(gca,'position',pos);

    %adjust the position of the second and 4th plots to reduce space between the
    %subplots
    new_xpos=0.42;

    pos=get(hs02,'position');
    pos_orig=pos;
    pos(1)=new_xpos;
    set(hs02,'position',pos);

    pos=get(hs04,'position');
    pos_orig=pos;
    pos(1)=new_xpos;
    set(hs04,'position',pos);

    %move 3 and 4 upwards
    new_ypos=0.15;
    gca=hs03;
    pos=get(gca,'position');
    pos_new=pos;
    pos_new(2)=new_ypos;
    set(gca,'position',pos_new);

    gca=hs04;
    pos=get(gca,'position');
    pos_new=pos;
    pos_new(2)=new_ypos;
    set(gca,'position',pos_new);


    %Make the colorbar wider and move to below the plot
    pos=get(Hc1s,'position');
    pos(1)=0.16;
    pos(2)=0.06;
    pos(3)=0.5;
    set(Hc1s,'position',pos);

    %Change the colormap to blue to brown :-
    hcbar = find_peer_colorbars_of_an_axes(gca);
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));    
    set(hcbar,'XaxisLocation','bottom');


end

%% Plot the PDFs for the different regions
if iplot_pdfs==1

    for iLATs=1:length(LATs)
        LAT_val_DRIVER = LATs{iLATs}; LON_val_DRIVER = LONs{iLATs};
        %run the script
        MODIS_vert_pen_VOCALS_tau_PDF_DRIVER_template_1D_PDF

        %plot the box on the map at the end
        ioverride_box_colour=1;
        LAT_val = LAT_val_DRIVER; LON_val = LON_val_DRIVER; %12th Nov UM FULL domain
        col_str='k-';  %
        figure(map_fig);
        plot_box_on_map
    end


end



