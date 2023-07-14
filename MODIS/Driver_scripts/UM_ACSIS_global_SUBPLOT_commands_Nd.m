    clear gca xlab Hc1s

 iload_Nd_time_match=1; %whether to load the data for the time match from disk or recalculate (since it is slow).
    

    fsize_latlon = 14;
    caxis_range = [0 300];    

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=3; %no. rows, columns
    figure('position',scrsz);    

   


    
    
%% ----  Model ---------------------------
    hs01=subplot(xsub,ysub,1);
    gca = hs01;
    subtitle_str=['UM (' Nd_type ')'];

   
    
    %Match the nearest model LT time to 13 LT (set time_choice.find_nearest=1)
    %for each day of th eMODIS data.
    %The can remove the days then MODIS is NaN.
    
    if  i_Iceland==1
        time_match_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS/UM_Nd_time_match_data.mat';
        
        if iload_Nd_time_match==1            
            load(time_match_file);            
        else
            
            [LT] = local_times_from_UTC(HH_UM,gcm_Plon2D_UM);
            %Y_UM_3D=repmat(Y_UM,[1 size(LT,1) size(LT,2)]); Y_UM_3D=repmat([2,3,1]);
            %[Y_UM_3D] = replicate_array_3D(Y_UM,LT);
            LT_matlab_time = datenum(replicate_array_3D(Y_UM,LT),replicate_array_3D(M_UM,LT),replicate_array_3D(D_UM,LT),LT,replicate_array_3D(zeros(size(Y_UM)),LT),replicate_array_3D(zeros(size(Y_UM)),LT));
            LT_matlab_time2 = reshape(LT_matlab_time,size(LT)); %Have checked that this reshaping gives the right ordering
            clear time_choice
            time_choice.find_nearest=0;
            %inot_tol_LT=find(abs(LT-13)>tol_hrs);
            time_choice.tol=3/24; %3 horus (convert to days)
            
            [date_str date_num] = date_from_day_of_year_func(modis_loaded.daynum_timeseries3,modis_loaded.modisyear_timeseries3);
            [Y,M,D]=datevec(date_num);
            time_choice.time_specific = datenum(Y,M,D,13*ones(size(Y)),zeros(size(Y)),zeros(size(Y)));
            
            %The following is quite slow, so better to save the data and load
            %from the file if possible
            dim=3;
            [out_Nd, time_out_Nd, time_inds_Nd, dtime_match_Nd] = get_time_range_of_array(Nd_PD_ALL,LT_matlab_time2,time_choice,dim);
            
            model_res=4; %km
            modis_res=111;%km
            N_coarse=floor(modis_res/model_res); M_coarse=N_coarse;
            
            
            %Re-grid the UM data to that of MODIS (1x1 deg)
            nT_modis = size(out_Nd,3);
            um_data2=NaN*ones(size(gcm_Plat2D_AMSRE,1),size(gcm_Plat2D_AMSRE,2),nT_modis);
            
            gcm_Plat2D_GENERIC = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse,N_coarse);
            gcm_Plon2D_GENERIC = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse,N_coarse);
            
            for it=1:nT_modis
                %Coarse grain from 4km model res to approx 1degree MODIS res first
                temp = reduce_matrix_subsample_mean(out_Nd(:,:,it),N_coarse,M_coarse);
                %Then regrid to the actual MODIS grid
                um_data2(:,:,it) = griddata(gcm_Plat2D_GENERIC,gcm_Plon2D_GENERIC,temp,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE);                
            end
            
            save(time_match_file,'out_Nd','time_out_Nd','dtime_match_Nd','um_data2','-V7.3');
        end
        

        %Model is 4km and MODIS is 1x1deg, so probably shoudl coarse-grain
        %rather than re-grid here?
        
        inan=find(isnan(Droplet_Number_Concentration_37.timeseries3)==1);
        um_data2(inan)=NaN;
        
        
        
        %    um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);
        
        %um_data =
        %
        dat_modis = meanNoNan(um_data2,3) / 1e6; %convert to per cm3
        
        dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        
    else
        dat_modis = meanNoNan(Nd_PD_ALL,3) / 1e6; %convert to per cm3
    end

    um_data = dat_modis; %save for later

    

    %Plotting Script
%    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %If plotting at UM native res

    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global    %
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


    

    %% ----  MODIS sat data ---------------------------
    hs02=subplot(xsub,ysub,2);
    gca = hs02;
    subtitle_str='MODIS';

    Nd_MODIS_dat=modis_loaded.Droplet_Number_Concentration_37.timeseries3;
    
    if i_Iceland==1
        inan = find(isnan(um_data2)==1);
        Nd_MODIS_dat(inan)=NaN;
    end
        
    
    [Nd37_1deg_time_mean, Ndatap_MODIS] = meanNoNan(Nd_MODIS_dat,3);  %N.B. - no time selection here
    
    if i_Iceland==1
        %inan=find(Ndatap_MODIS<2);
        %Nd37_1deg_time_mean(inan)=NaN;                        
    end
    
    
    %Regrid MODIS data to that of UM (2x2 deg for N96, 4km for Iceland)
    dat_modis = griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,Nd37_1deg_time_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);  
    %dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,Ndatap_MODIS,gcm_Plat2D_UM,gcm_Plon2D_UM); 
    
    sat_data = dat_modis; %back up data for future use perhaps


    %Run plotting script
    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
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


    
    
    
%% ----  Difference plot ---------------------------
    caxis_range_diff = [-50 50];  
    
    hs03=subplot(xsub,ysub,3);
    gca = hs03;
    subtitle_str='Bias UM minus satellite';


    dat_modis = um_data - sat_data;            

    %Run plotting script
    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    title([titlenam_driver ' mean=' num2str(Pmean,'%.2f')]);
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
    if  i_Iceland==1
        pos(2)=0.34; cb_vert_pos = pos(2);
    else
        pos(2)=0.24; cb_vert_pos = pos(2);
    end
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
    xlab = xlabel(Hc1s,'Droplet concentration bias (cm^{-3})','fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);    

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,['Droplet concentration (cm^{-3})'],'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  

 if isave_plot_Nd_bias==1
%        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_Nd_bias_3panel_vs_MODIS'];
        titlenam_driver=['global_Nd_bias_3panel'];
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %        close(gcf);
%        
 end