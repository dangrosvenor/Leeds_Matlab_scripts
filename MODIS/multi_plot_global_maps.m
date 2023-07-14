year_multi = 2014;

for days_required_for_mean=245:245   %244:304 
    time_mean_str = ['t=' num2str(days_required_for_mean) ' ' datestr( datenum(year_multi,1,1) + days_required_for_mean - 1) ' '];
    ioverride_time_selection=1;
    plot_global_maps
    saveas_ps_fig_emf(gcf,[savename]);
    close(gcf);
end

