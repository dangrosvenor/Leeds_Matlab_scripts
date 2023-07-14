for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    var_UM_clean = remove_character(var_UM ,'.','pt');
    var_UM_clean2 = remove_character(var_UM ,'_',' ');
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM_clean '_PI_ALL'];
    var_name_out = var_name_out_PI;    
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM_clean '_PD_ALL'];
    var_name_out = var_name_out_PD;    
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    eval_str=[var_UM_clean '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
    eval_str=[var_UM_clean '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);
    
    eval_str=['d' var_UM_clean ' = meanNoNan(' var_name_out_PD ' - ' var_name_out_PI ',3);']; eval(eval_str);
    
    
    iplot_SW=1;
    if iplot_SW==1
        figure
        var_UM = var_UM_clean;    
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10; 
        dat_modis = eval([var_UM '_PI']); var_UM = ['Time mean ' var_UM_clean2 ' (Volcano OFF)'];
        %run plotting script        
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 4000]);
        icoarse_grain=0;
        
        figure
        var_UM = var_UM_clean;    
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10; 
        dat_modis = eval([var_UM '_PD']); var_UM = ['Time mean ' var_UM_clean2 ' (Volcano ON)'];
        %run plotting script        
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 4000]);
        icoarse_grain=0;        
        
        
        figure
        %ilabel_colorbar = 1; col_bar_lab_str = '(g m^{-2})';
        var_UM = var_UM_clean;
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;        
        %dat_modis = ( SO2_PD - SO2_PI); var_UM = 'Change in SO_2 (PD minus PI; kg m^{-2})';
        dat_modis = eval(['d' var_UM]); var_UM = ['Change in ' var_UM_clean2 ' (PD minus PI)'];
        %run plotting script        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-400 400]);
        icoarse_grain=0;       
        
        
        var_UM = var_UM_clean;
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = ['% Change in ' var_UM_clean2] % \DeltaLWP_{ic subgrid} (g m^{-2})';
        %var_UM = '';
        dat_modis = eval(['100*(' var_UM '_PD - '  var_UM '_PI )./ ' var_UM '_PI;']);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '(%)';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-50 50]);
        icoarse_grain=0;
        
        var_UM = var_UM_clean;
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = ['Change in timemean ' var_UM_clean2] % \DeltaLWP_{ic subgrid} (g m^{-2})';
        %var_UM = '';
        dat_modis = eval([var_UM '_PD - '  var_UM '_PI;']);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-50 50]);
        icoarse_grain=0;
        
    end
    
end