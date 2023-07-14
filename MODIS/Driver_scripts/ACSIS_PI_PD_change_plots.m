icoarse_grain=1; M_coarse_grain=3; N_coarse_grain=3;
isave_plot_global_LWP=1; %for draft paper
icontour_DRIVER=0;
ioverride_ticks_DRIVER=0;


 
    %dat_modis = 100 * (meanNoNan(f1,3)./meanNoNan(f0,3) - 1); var_UM = 'Percentage change in cloud fraction (PD minus PI)';    
    %dat_modis = meanNoNan(f1,3) - meanNoNan(f0,3); var_UM = 'Absolute change in cloud fraction (PD minus PI)';  
    PD_mean = meanNoNan(PD_data,3);
    PI_mean = meanNoNan(PI_data,3);
    dat_modis = 100*(PD_mean - PI_mean) ./ PI_mean; 
    var_UM = ['Percentage change in ' var_str_tit ' (PD minus PI)'];
    %tit_str_clean='% f_c change';
    subtitle_str=['% ' var_str_tit ' change'];
    dat_modis(dat_modis>999)=999; %to deal with gridpoints with Inf values from where f0_mean=0 
    
    um_data = PD_mean;
    sat_data = PI_mean; %for the stats - special case where need regional calcs of two things first before doing percentage diff/change
    dtimemean=dat_modis;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-0.1 0.1]);   
    caxis(crange);      
    
      % Save
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date titlenam_driver];          
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        save_map_plot_data
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);

           
        DRIVER_calc_biases_for_regions
    end 
    
    