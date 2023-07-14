thresh=[10 20 30 40 50 60 70 80];


%some default flags - some may be overwritten below
thresh_SZA=60;
thresh_CF=0.8;
thresh_NP=50;
thresh_sensZA=40;
thresh_sensZA=30;
thresh_maxSZA=60;


screen_type='NP + CF + MEAN sensZA';
screen_type='NP + CF + MAX sensZA';

irestrict_domain=0;

                                        
for ithresh=1:length(thresh)
    
    %set the thresholds required and set the flag to override the
    %thresholds in plot_global_maps.m
    thresh_sensZA = thresh(ithresh);
    proj_type='global oval';
    data_select='specific_modis_data';
    mod_data_type='timeseries3';
    modis_data_plot='Number of droplets cell values time mean'; %from timeseries3
    
    ioverride_plotglobal_thresh=1; %using this one flag for all settings - clear on last use within plot_global
    
    %%%%  Have tried to automate to make sure the settings are set to do a plot of droplet number
    %%%%  plot! Quite a few flags to set, though - any been missed?
    plot_global_maps
    
    %now run the bootstrap analysis with 1000 samples with the screened
    %result of plot_global_maps
    boot_out = bootstrap_array(dat_modis,1000,[2.5 5 7.5 10 20 30 50 70 80 90 92.5 95 97.5]);
        
    %save the boot_out array to disk
    ioverride_saveload=1;
    saveload='save'; % - make sure that the flag is set to save!
    saveload_various_MODIS 
    
end
    