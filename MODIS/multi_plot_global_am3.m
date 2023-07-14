low_cf_thresh2 = [0.01 0.2 0.4 0.6 0.8]; %
hi_cf_thresh2 = [0.2 0.4 0.6 0.8 1.01]; %

for imulti_am3=1:length(low_cf_thresh2)
    low_cf_thresh(1)=low_cf_thresh2(imulti_am3);
    low_cf_thresh(2)=hi_cf_thresh2(imulti_am3);
    
    plot_global_maps
    saveas_ps_fig_emf(gcf,savename,'0_150_cscale',0);
    
end