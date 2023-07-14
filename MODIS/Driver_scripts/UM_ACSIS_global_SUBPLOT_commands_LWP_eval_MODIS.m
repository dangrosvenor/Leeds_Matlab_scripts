clear gca xlab Hc1s hs

isub=0;

iuse_CALIPSO_COSP=1; %Whether to use COSP CALISPO values for the model, or standard model CF

%    isave_plot_LWP_modis = 1;

fsize_latlon = 14;
caxis_range = [0 300];

subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
xsub=1; ysub=3; %no. rows, columns
figure('position',scrsz);

bias_type = 'Absolute';
bias_type = 'Percentage';





%% ----  Model ---------------------------
hs01=subplot(xsub,ysub,1);
isub=isub+1; hs{isub} = hs01;
gca = hs01;

% Do time matching of model to satellite.
%Get the 3d array of local times for each location  in [lat,lon,time]
%order.
[LT] = local_times_from_UTC(HH_UM,gcm_Plon2D_UM);
tol_hrs = 3; %tolerance for the matching (+/- N hours)
switch lwp_day_night
    case 'average'
        inot_tol_LT=find(abs(LT-13.5)>tol_hrs & abs(LT-1.5)>tol_hrs);
    case 'day'
        inot_tol_LT=find(abs(LT-13.5)>tol_hrs);
    case 'night'
        inot_tol_LT=find(abs(LT-1.5)>tol_hrs);
end


switch model_lwp_var
    case 'Model COSP MODIS LWP'      
        var_Latex = 'LwpMODIScosp';        
        um_data_temp = modis_cosp_lwp_PD;
        
    case 'Model LWP 2-391'
        var_Latex = 'LwpUM';        
        um_data_temp = LWP_2_391_PD;
end

subtitle_str=[model_lwp_var ' ' lwp_day_night ' ' modis_CTH_screening ' lwp cf thresh=' num2str(lwp_cf_thresh) ' idiv lwp CF=' num2str(idiv_lwp_CF)];

um_data_temp(inot_tol_LT) = NaN;

um_data = 1e3*meanNoNan(um_data_temp,3);

%Degrade UM data to that of CALIPSO (2x2 deg)
%um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);
%dat_modis = um_cal_res;
dat_modis = um_data;


%Plotting Script
var_UM = subtitle_str;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %If plotting at UM native res
%    CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global    %CALIPSO res
titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error
%title([titlenam_driver ' mean=' num2str(Pmean,'%.2g')]);
caxis(caxis_range);
gca = hs01;
increase_font_size_map_figures
%move colorbar down a bit and expand to full width

%N.B. need to get a new handle each time want to change a different
%colobar  - doesn't seem to allow multiple handles
gca = hs01; Hc1s = find_peer_colorbars_of_an_axes(gca);
delete(Hc1s)  %just need one big colorbar for plots 1 and 2


%% ----  MODIS sat data - not sure if doing correct time matching here ---------------------------
hs02=subplot(xsub,ysub,2);
isub=isub+1; hs{isub} = hs02;
gca = hs02;
subtitle_str='MODIS LWP';

%     switch cloud_alt
%         case 'Low altitude cloud'
%             dat_modis = meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
%         case 'Mid altitude cloud'
%             dat_modis = meanNoNan(clmcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
%         case 'High altitude cloud'
%             dat_modis = meanNoNan(clhcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
%         case 'Total cloud'
%             dat_modis = meanNoNan(cltcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
%     end
%
%     cal_data = dat_modis;


switch modis_CTH_screening
    case '3.2km, CF>80'
        %Degrade data to that of UM (2x2 deg for N96)
        modis_LWP_time_period_mean = 1e3*meanNoNan( modis_loaded.W_time3 ,3);
        dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
    case '3.2km';
        %Degrade data to that of UM (2x2 deg for N96)
        modis_LWP_time_period_mean = 1e3*meanNoNan( modis_loaded_CF0.W_time3 ,3);
        dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
    case 'none'
        %Degrade data to that of UM (2x2 deg for N96)
        modis_LWP_time_period_mean = 1e3*meanNoNan( modis_loaded_CF0_2.W_time3 ,3);
        dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
    case 'AMSR-E'
        dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,1e3*amsre_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
        dat_modis = flipdim(dat_modis,1); %Flip the lat axis since AMSRE is oriented differently to MODIS
        subtitle_str='AMSR-E all-sky LWP';
        
end



%cf_modis_timemean = meanNoNan(modis_loaded_CF0_2.Cloud_Fraction_Liquid.timeseries3,3);

%Regrid to grid of UM
%sat_data = griddata(modis_loaded_CF0_2.gcm_Plat2D_AMSRE,modis_loaded_CF0_2.gcm_Plon2D_AMSRE ,cf_modis_timemean,gcm_Plat2D_UM,gcm_Plon2D_UM);
sat_data = dat_modis;


%Run plotting script
%CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
var_UM = subtitle_str;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%titlenam_driver = [subtitle_str];
%title([titlenam_driver ' mean=' num2str(Pmean,'%.2g')]);
caxis(caxis_range);
gca = hs02;
increase_font_size_map_figures
%move colorbar down a bit and expand to full width

%N.B. need to get a new handle each time want to change a different
%colobar  - doesn't seem to allow multiple handles
Hc1s = find_peer_colorbars_of_an_axes(hs02);







%     %make them all bigger - seem to have to increase the height rather than the
%     %width for this
%     new_height=0.25;
%     gca=hs01; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%     gca=hs02; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%     gca=hs03; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%
%












%% ----  Difference plot ---------------------------

iplot_box_on_map = 0;

hs03=subplot(xsub,ysub,3);
isub=isub+1; hs{isub} = hs03;
gca = hs03;
subtitle_str='Model bias';

%    sat_data = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,cal_data,gcm_Plat2D_UM,gcm_Plon2D_UM);

%abs_bias = um_cal_res - cal_data;
%prc_bias = 100 * abs_bias ./ cal_data;
abs_bias = um_data - sat_data;
prc_bias = 100 * abs_bias ./ sat_data;
bias_str='LWP bias';

switch bias_type
    case 'Absolute'
        dat_modis = abs_bias;
        caxis_range_diff = [-0.3 0.3];
    case 'Percentage'
        dat_modis = prc_bias;
        caxis_range_diff = [-300 300];
        bias_str = [bias_str ' (%)'];
end







%Run plotting script
var_UM = subtitle_str;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error
%title([titlenam_driver ' mean=' num2str(Pmean,'%.2g')]);
caxis(caxis_range_diff);
gca = hs03;
increase_font_size_map_figures
%move colorbar down a bit and expand to full width

%N.B. need to get a new handle each time want to change a different
%colobar  - doesn't seem to allow multiple handles
Hc1s = find_peer_colorbars_of_an_axes(hs03);


%adjust the position of the 2nd and 3rd plots to reduce space between the
%subplots
new_xpos=0.35;
pos=get(hs02,'position');
pos_orig=pos;
pos(1)=new_xpos;
set(hs02,'position',pos);

new_xpos = 0.57;
pos=get(hs03,'position');
pos_orig=pos;
pos(1)=new_xpos;
set(hs03,'position',pos);


%Need to sort move the colourbars separately



%Make the 1st colorbar wider and move to below the plot
Hc1s = find_peer_colorbars_of_an_axes(hs02);
pos=get(Hc1s,'position');
pos(1)=0.12;
pos(2)=0.24; cb_vert_pos = pos(2);
pos(3)=0.4;
set(Hc1s,'position',pos);


%Adjust the 2nd colorbar
Hc1s = find_peer_colorbars_of_an_axes(hs03);
pos=get(Hc1s,'position');
pos(1)=0.5698;
pos(2)=cb_vert_pos;
pos(3)=0.17;
set(Hc1s,'position',pos);

%Change the colormap to blue to brown :-
lb_map = lbmap(16,'brownblue');
colormap(flipdim(lb_map,1));
%    dz_tit = 0.3;
dz_tit = -0.7;

set(Hc1s,'XaxisLocation','bottom');
xlab = xlabel(Hc1s,bias_str,'fontsize',16);
set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
pos = get(xlab,'position');
set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);

Hc1s = find_peer_colorbars_of_an_axes(hs02);

set(Hc1s,'XaxisLocation','bottom');
xlab = xlabel(Hc1s,['Liquid water path (g m^{-2})'],'fontsize',16);
set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
pos = get(xlab,'position');
set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);


if isave_plot_LWP_modis==1
    %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_CF_bias_3panel_' cloud_alt];
    titlenam_driver=['global_CF_bias_3panel_LWP_MODIS'];
    savename=[savedir_date titlenam_driver];
    
    % Run script to calculate biases, etc. for the different regions.
    col_str = 'k';
    axes(hs02); %switch to the satellite plot
    DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region
    
    %do the a,b,c, etc. labels
    label_subfigs
    
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
    
    
end

