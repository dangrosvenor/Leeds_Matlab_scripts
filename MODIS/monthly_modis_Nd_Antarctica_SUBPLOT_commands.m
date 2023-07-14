clear gca xlab Hc1s

plot_global_maps_defaults

%isave_plot_CERES=0;
fsize_latlon = 14;  
caxis_range = [0 150];
    
var_name_units = 'N_d (cm^{-3})';


    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    inew_figure = 0;
    xsub=3; ysub=2; %no. rows, columns
    figure('position',scrsz);

    
%% Loop over months
      for imon=[1:3 10:12]
            
            mon_str = datestr(datenum(2008,imon,1),'mmm');

            switch channel_str
                case '21'

                case '37'
                    dat_modis = Nd_multi_annual(:,:,imon) / icount; %permute(modisL3_dat.Nd_1km,[2 1 3]);

            end

       
    hs01=subplot(xsub,ysub,imon);
    gca = hs01;


      
    subtitle_str=[mon_str];

  %--- run the file to set up the defaults

            ioverride_plotglobal_loc=1;
            ioverride_plotglobal_thresh=1;  %comment out if want to use screenings set in plot_global
            titlenam_driver = ['N_d multi-year ' mon_str ' mean'];
            %time override should aready be set (ioverride_time_selection)
            ioverride_years_time_screen=1; %required to specify the different years
            iover_ride_plot_global=1; %overrides inew_figure=1; supress_colorbar=0; i_increase_font_size_map_figures_OFF = 0;
            inew_cticks=0;  %colorbar is of the non-linear type
            %modis_data_plot = 'Generic plot specified outside of script';
            modis_data_plot = 'Map of 2D data from outside driver script';
            proj_type='polar'; stereo_str01='lat_polar=-90'; stereo_str02='m_proj(''stereographic'',''lat'',lat_polar,''lon'',-98,''rad'',55)';
            plot_global_maps

 
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
      end
    
    %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    gca = hs01; Hc1s = find_peer_colorbars_of_an_axes(gca);    
    delete(Hc1s)  %just need one big colorbar for plots 1 and 2


    
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
    xlab = xlabel(Hc1s,['Bias in ' var_name_units],'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);    

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,var_name_units,'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  

    
    
     
    
     if isave_plot_CERES==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_SW_' titlenam_driver '_vs_CERES'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        opts.isavefig=1; %For Olympus
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%        close(gcf);
    end