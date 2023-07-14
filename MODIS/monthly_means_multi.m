

LAT_vals = [-22 -18; -22 -18; -22 -18; -22 -18]; LON_vals = [-76.5 -72; -81 -76.5; -85.5 -81; -90 -85.5];

LAT_vals = [-22 -18; -22 -18; -22 -18; -20 -10]; LON_vals = [-75 -70; -80 -75; -86 -80; -90 -80];

LAT_val2 = [-22.74 -18];
LON_val2 = [-120:1-71.25];
icount=1;
clear LON_vals LAT_vals
for irep=1:length(LON_val2)-1
    LON_vals(icount,:)=LON_val2(irep:irep+1);
    icount=icount+1;
end
LAT_vals = repmat(LAT_val2,[size(LON_vals,1) 1]);

clear time_means_multi LAT_str_multi LON_str_multi val_save lon_save
for ilats_multi=1:size(LAT_vals,1)
    LAT_val=LAT_vals(ilats_multi,:);
    LON_val=LON_vals(ilats_multi,:);
    
    ioverride_location_selection2=1;
    monthly_means_from_plot_global
    
    
   time_means_multi(ilats_multi).dat = seasonal_meanX; 
   LAT_str_multi(ilats_multi).dat = LAT_str; 
   LON_str_multi(ilats_multi).dat = LON_str;    
%   time_labs_multi(ilats).dat = season_vals_save(ilats).pdflab;
   Xmean_multi(ilats_multi).dat = sum(seasonal_NX.*seasonal_meanX)/sum(seasonal_NX);  
   
   %special case
   val_save.modis(ilats_multi)=time_means_multi(ilats_multi).dat(1);
   lon_save.modis(ilats_multi)=mean(LON_vals(ilats_multi,:));
   
end

    
    