lat_range=-70:-50;

for ilat_range=1:length(lat_range)
    ioverride_location_selection=1;
    LAT_val = [lat_range(ilat_range):lat_range(ilat_range)+1]; LON_val = [-120:-90]; %Antarctica, random region
    plotTimeHeightVap3
    saveas_ps_fig_emf(gcf,[savename]);
    close(gcf);
end

