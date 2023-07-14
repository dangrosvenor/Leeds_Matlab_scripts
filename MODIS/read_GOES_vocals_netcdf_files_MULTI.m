
save_file_amsre = ['/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_' datestr(now,30) '.mat'];

noplot=1;

iplot_goes_map=1; %whether to plot and save a map of the GOES field

thresh_LWP_DRIVER = 5; %LWP threshold for Nd calculation from GOES (g/m2)

%26th Oct POC case
%LAT_val = [-21 -15.44]; LON_val = [-86.93 -82]; %Same as for revised PDFs post-AGU. Revised smaller region in NW corner  
%date_range = [datenum('25-Oct-2008 12:00') datenum('26-Oct-2008 21:00')];

% 12th Nov case
%Need to come up with a new corner for the 12th Nov case since this is
%centred at 76W, 20S (not at 82W).
% Domain runs from -78.93 to -73.08 W and -22.70 to -17.28 (based on the edges)
LAT_val = [-22.70 -17.28]; LON_val = [-78.93 -73.08]; %FULL UM domain for 12th Nov - Latest file:-
%/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151116T084819.mat  %FULL UM domain for 12-15th Nov
% For the shorter period of [datenum('13-Nov-2008 12:10')
% datenum('13-Nov-2008 21:20')]; use /home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151119T042920.mat

% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
%LAT_val= [-20.5 -17.5]; LON_val = [-78.75 -73.25];
%Latest file for this :- load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151112T064655.mat';

%The location of the ship = 20S, 75W
%One degree region around it
%LAT_val = [-20.5 -19.5]; LON_val = [-75.5 -74.5];

LAT_val = [-1e9 1e9]; LON_val = [-1e9 1e9]; %the whole GOES region
%Latest file = /home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat

date_range = [datenum('12-Nov-2008 00:00') datenum('15-Nov-2008 00:00')];
date_range = [datenum('13-Nov-2008 12:10') datenum('13-Nov-2008 21:20')];
date_range = [datenum('12-Nov-2008 12:10') datenum('12-Nov-2008 21:20')];

goes_dir = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/';

files=dir([goes_dir '/*.nc']);

clear files2 goes_LWP_mean goes_LWP_multi goes_Tau_multi goes_Reff_multi times_GOES_save goes_Nd_mean goes_Nd_multi

i2=0;
for i=1:length(files)
   year = str2num(files(i).name(23:26));
   month = str2num(files(i).name(27:28));
   day = str2num(files(i).name(29:30));   
   hour = str2num(files(i).name(31:32));   
   mins = str2num(files(i).name(33:34)); 
   
   time = datenum(year,month,day,hour,mins,0);   
   [day_of_year,day_str] = day_of_year_from_date_func(time);
   

   
   if time>=date_range(1) & time<=date_range(2)
       i2=i2+1;
       files2(i2) = files(i);    
       times_GOES_save(i2)=time;
       daynum_timeseries3_GOES_save(i2)=day_of_year;
       modisyear_timeseries3_GOES_save(i2)=year;
       gcm_time_UTC_GOES_save(i2)=hour+mins/60;
   end
end


for itime_amsre=1:length(files2)

%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file

goes_file_to_load = files2(itime_amsre).name;




        
%------- Calculate the data to plot
         %read in the GOES data for the specific time
          %read in the GOES data for the specific time
         ioverride_goes = 1;
         goes_action = 'load a particular file';        
         read_GOES_vocals_netcdf_files 
         
if iplot_goes_map==1         
%-------- Plot a map of latest file based on settings in the script below
         ioverride_goes_multi=1;
         isave_plot_driver=1;
         %goes_file_to_load is then used in this script
         POC_26Oct2008_GOES_general_maps_20141125T032943
end
       
% ----- Use pdf2d to pick out the required region of interest and to
% calculate LWP and Nd (better to just do this directly?)
        
% ----- Set various things

%        mod_data_type='AMSRE';
        gcm_str_select='GOES';
        gcm_str='GOES';
       
        month_amsre = goes_month;
        year_amsre = goes_year;

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
%         watervap_defaults
         pdf2D_defaults
         
        
        %--- set some options for these particular plot loops
%        set_screening = {'none'};
%        modis_data_plot = 'Map of 2D data from outside driver script';
        i577 = 'MODIS_plot_UW';

        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=200;
        
        logflag=0;
        dlogflag=0;
        
        isave_plot=0;
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
        
        
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
        y_axis_vals = 'GOES LWP';
        datatype = 'gcm_data';          
        

                                
          axis1D = 'y';                                
                                
       

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;

        % --------- Override flags for watervap --------
        man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
        plotTimeHeightVap3
        close(gcf);
        
        goes_LWP_multi{itime_amsre} = goes_LWP_save;
        goes_Tau_multi{itime_amsre} = goes_Tau;
        goes_Reff_multi{itime_amsre} = goes_Reff;
        goes_LWP_mean(itime_amsre) = Y_mean_overall;
        
% --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;
        
        
        
%--- Nd ---        
        x_axis_vals = 'Dummy data'; %dummy data
        y_axis_vals = 'GOES Nd';        

        % --------- Override flags for watervap --------
        man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
        plotTimeHeightVap3
        close(gcf);
        
        goes_Nd_multi{itime_amsre} = goes_Nd_save;
        goes_Nd_mean(itime_amsre) = Y_mean_overall;    
        
%--- Teff
        Y = goes_Teff;
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        goes_Teff_multi{itime_amsre}=Y;
        
%--- CTH
        Y = goes_CTH;
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        goes_CTH_multi{itime_amsre}=Y;        

%--- Vis
        Y = goes_Vis;
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        goes_Vis_multi{itime_amsre}=Y; 
        
%--- IRemit
        Y = goes_IRemit;
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        goes_IRemit_multi{itime_amsre}=Y;
        
        
        
end


daynum_timeseries3_GOES = daynum_timeseries3_GOES_save;
modisyear_timeseries3_GOES = modisyear_timeseries3_GOES_save;
gcm_time_UTC_GOES = gcm_time_UTC_GOES_save;

gcm_Plat2D_GOES_nonan = inpaint_nans(gcm_Plat2D_GOES);
gcm_Plon2D_GOES_nonan = inpaint_nans(gcm_Plon2D_GOES);


save(save_file_amsre,'-V7.3','times_GOES_save','goes_LWP_mean','thresh_LWP_DRIVER',...
    'goes_Nd_mean','gcm_Plon2D_GOES','gcm_Plon2D_edges_GOES','gcm_Plat2D_edges_GOES','gcm_Plat2D_GOES',...
    'goes_CTH_multi','goes_Vis_multi','goes_IRemit_multi','goes_LWP_multi','goes_Tau_multi','goes_Reff_multi',...
    'goes_Nd_multi','goes_Teff_multi','daynum_timeseries3_GOES','gcm_time_UTC_GOES','gcm_time_UTC_GOES',...
    'gcm_Plat2D_GOES_nonan','gcm_Plon2D_GOES_nonan');


