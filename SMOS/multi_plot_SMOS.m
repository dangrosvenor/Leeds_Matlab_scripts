savedir_SMOS='/home/disk/eos8/d.grosvenor/SMOS/SMOS_L3_CPDC/plots/Africa/daily/desc/';
%savedir_SMOS='/home/disk/eos8/d.grosvenor/SMOS/SMOS_L3_CPDC/plots/every_10_days/desc/';

for itime_SMOS=1:1:363
    plot_global_maps;
    savename = [savedir_SMOS 'SMOS_daily_' day_str];
    
saveas_ps_fig_emf(gcf,savename,'');
close(gcf);
    
end