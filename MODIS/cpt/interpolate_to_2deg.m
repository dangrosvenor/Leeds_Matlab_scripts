gcm_str_select2 = 'CAM5_CLUBBv2_COSP';
gcm_str_select2 = 'AM3';
gcm_str_select2 = 'CAM5';


time_strs_2deg = {'DJF','MAM','JJA','SON'};

for itime_2deg=1:length(time_strs_2deg)
    
    
    time_mean_str_2deg = time_strs_2deg{itime_2deg};

time_mean_str = time_mean_str_2deg;
 
switch time_mean_str
                case 'DJF'
                    days_required_for_mean = [336:366 1:60];
                case  'MAM'    
                    days_required_for_mean = [61:152];
                case 'JJA'
                    days_required_for_mean = [153:244]; 
                case 'SON'
                    days_required_for_mean = [245:335];     
                case 'ALL'
                    days_required_for_mean = [1:365];                         
end    

 ioverride_time_selection=1;
 ioverride_month_select_multi_years=1;
 ioverride_plotglobal_thresh=1;
 modis_data_plot='LTS ECMWF'; mod_data_type='GCM';   gcm_time_of_day_select=2; 
 gcm_str_select ='ERAInt';
 plot_global_maps
 LTS_ec = P;
 
 ioverride_time_selection=1;
 ioverride_month_select_multi_years=1;
 ioverride_plotglobal_thresh=1;
 modis_data_plot='LTS GCM'; mod_data_type='GCM';   gcm_time_of_day_select=2;
 gcm_str_select = gcm_str_select2;
 plot_global_maps
 LTS_am3 = P;
 


%time averages
%LTS_ec = meanNoNan(ecLTS_mon_Man_ALL,1);
%LTS_am3 = meanNoNan(gcm_LTS_AM3,1);


%ZI = GRIDDATA(X,Y,Z,XI,YI);

MLON = eval([' [floor(minALL(gcm_Plon2D_' gcm_str_select ')) : 2 : ceil(maxALL(gcm_Plon2D_' gcm_str_select ')) ];']);
MLAT = eval(['[floor(minALL(gcm_Plat2D_' gcm_str_select ')) : 2 : ceil(maxALL(gcm_Plat2D_' gcm_str_select ')) ];']);
[X,Y]= meshgrid(MLON,MLAT);

daynum_timeseries3_MODIS = 1;
modisyear_timeseries3 = 1;


LTS_am3_2 = eval(['GRIDDATA(gcm_Plon2D_' gcm_str_select ',gcm_Plat2D_' gcm_str_select ',LTS_am3,X,Y);']);

%ECMWF

LTS_ec_2 = GRIDDATA(gcm_Plon2D_ERAInt,gcm_Plat2D_ERAInt,LTS_ec,X,Y);



 ioverride_time_selection=1;
 time_mean_str=time_mean_str_2deg;
 ioverride_month_select_multi_years=1;
 ioverride_plotglobal_thresh=1;
 days_required_for_mean = 1;
 modis_data_plot='LTS bias from ECMWF'; mod_data_type='Monthly_data';   gcm_time_of_day_select=2;
 gcm_str_select ='ERAInt';
 plot_global_maps
 

 titlenamw=textwrap({short_plot_name},70);
 title(titlenamw,'fontsize',fsize+2);
 
 saveas_ps_fig_emf(gcf,[savename],'',0,1)
 
end
 
 








