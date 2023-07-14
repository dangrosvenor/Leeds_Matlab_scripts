
% thresh_LAT = [0 85]; thresh_LON = [-180 180]; %N. hemisphere - whole globe
% thresh_LAT = [30 80];  thresh_LON = [-90 80]; %UK and Iceland
% 
% 
% thresh_LAT = [0 60]; thresh_LON = [90 180]; %China


lats_multi = {[0 85],[30 80],[0 60]};
lons_multi = {[-180 180],[-90 80],[90 180]};

%just China
lats_multi = {[0 60]};
lons_multi = {[90 180]};

%just N. Hemisphere and UK/Europe
lats_multi = {[0 85],[30 80]};
lons_multi = {[-180 180],[-90 80]};
%VOCALS 
lats_multi = {[-22.74 -18],[-22.74 -18],[-22.74 -18],[-22.74 -18]}; 
lons_multi = {[-76.25 -71.25], [-81.25 -76.25], [-87.25 -81.25], [-91.25 -81.25]};
lats_multi = {[-21.5 -18.5],[-21.5 -18.5],[-21.5 -18.5]};
lons_multi = {[-76.25 -71.25], [-81.25 -76.25], [-87.25 -81.25]};

 thresh_ndays=1;
 
 clear location_vals_save

for imulti_loc=1:length(lats_multi)

    thresh_LAT = lats_multi{imulti_loc};
    thresh_LON = lons_multi{imulti_loc};
    
    LAT_val = lats_multi{imulti_loc};
    LON_val = lons_multi{imulti_loc};

    ioverride_location_selection2=1;

    monthly_means_from_plot_global
    
%     this is what monthly_means_from_plot_global saves for each season
%     season_vals_save(iseason).NdPDF = ydat(1).y;
%     season_vals_save(iseason).NdPDF_xbins = xdat(1).x;
%     season_vals_save(iseason).NdPDF_xlab = xlabelstr;
%     season_vals_save(iseason).pdflab = pdflab;

%save these in a cell structure
location_vals_save{imulti_loc} = season_vals_save;
location_vals_save{imulti_loc}(1).lat_lab = LAT_str;
location_vals_save{imulti_loc}(1).lon_lab = LON_str;

%clear season_vals_save %this is done inside monthly_means_from_plot_global
    

end