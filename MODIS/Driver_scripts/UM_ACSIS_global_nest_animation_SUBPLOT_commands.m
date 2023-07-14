clear gca xlab Hc1s

isave_plot_anim = 1;
fsize_latlon = 14;  
caxis_range = [0 200];

%No. colour grades for the colorbar
ncols=20; %Good ot change this to match the number of tick marks on the colorbar

var_name_units = 'Liquid Water Path (g m^{-2})';


    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure('position',scrsz);

   
    
%% ----  Global model ---------------------------
    hs01=subplot(xsub,ysub,1);
    gca = hs01;


    %um_data = 
    dat_modis = LWP_global_ALL(:,:,it_anim);
    gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
    gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
    gcm_Plat2D_edges_UM = dat_global.gcm_Plat2D_edges_UM;
    gcm_Plon2D_edges_UM = dat_global.gcm_Plon2D_edges_UM; 
      
    subtitle_str=['Global N96 UM'];

    
    %Degrade UM data to that of CALIPSO (2x2 deg)
%    um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);    
    

    um_data = dat_modis; %save for later

    

    %Plotting Script
%    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %If plotting at UM native res
    %MODIS_ACSIS_LWP_global_vs_nest_quick_plot_commands_global    %
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error
    title([titlenam_driver ' mean=' num2str(Pmean,'%.2f')]);    
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    gca = hs01; Hc1s = find_peer_colorbars_of_an_axes(gca);    
    delete(Hc1s)  %just need one big colorbar for plots 1 and 2


%% ----  Nest data ---------------------------

    
    hs02=subplot(xsub,ysub,2);
    gca = hs02;
    subtitle_str='4km Nested UM';   
    
    dat_modis = LWP_nest_ALL(:,:,it_anim);
    gcm_Plat2D_UM = dat_nest.gcm_Plat2D_UM;
    gcm_Plon2D_UM = dat_nest.gcm_Plon2D_UM;
    gcm_Plat2D_edges_UM = dat_nest.gcm_Plat2D_edges_UM;
    gcm_Plon2D_edges_UM = dat_nest.gcm_Plon2D_edges_UM; 
    


    subtitle_str = [subtitle_str];


    
    sat_data = dat_modis; %back up data for future use perhaps


    %Run plotting script
%    MODIS_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)    
    titlenam_driver = [subtitle_str];        
    title([titlenam_driver ' mean=' num2str(Pmean,'%.2f')]);
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



    %adjust the position of the 2nd and 3rd plots to reduce space between the
    %subplots
    new_xpos=0.45;;
    pos=get(hs02,'position');
    pos_orig=pos;
    pos(1)=new_xpos;
    set(hs02,'position',pos);
       

%Need to sort move the colourbars separately    
  
    %Make the 1st colorbar wider and move to below the plot
    Hc1s = find_peer_colorbars_of_an_axes(hs02);    
    pos=get(Hc1s,'position'); % l b w h
    pos(1)=0.12;
    pos(2)=0.24; cb_vert_pos = pos(2);
    pos(3)=0.605;
    set(Hc1s,'position',pos);
    
    


    %Change the colormap to blue to brown :-
    lb_map = lbmap(ncols,'brownblue');
    colormap(flipdim(lb_map,1));    
%    dz_tit = 0.3;
    dz_tit = -0.7;
    
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,var_name_units,'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  

    
    
     if isave_plot_anim==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/' num2str(it_anim,'%02i') '_anim_LWP_' titlenam_driver '_and_global'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%        close(gcf);
    end