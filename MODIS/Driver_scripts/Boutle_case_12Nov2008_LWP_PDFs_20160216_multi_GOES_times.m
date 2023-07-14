% PDFs of LWP

% N.B - need to load the AMSRE data in all cases for coarsening

% For AMSRE at the moment will use the loaded file - fix later.
    %    Using 
    %         multi_read_amsre_daily 
    % to load the amsre file
    % Sep 2016 - still need to run multi_read_amsre_daily
    % Also set plot_remss to zero
    
% If want to calculate LWP from qL and rho then use the routine
% calc_LWP_multi.m and save as .mat
% rather than calculating here. Then specify the .mat file here.
% For GOES it loads in the requested file using
    % read_GOES_vocals_netcdf_files
            
%% Note that pole_lon = 284 for this case (centred at 76W)    

% Use data stored in the remss array rather than the AMSR-E data
% Although only have to choose the right file for the UM domain
plot_remss=0;
tol_remss = 1.5/24; %in days
iload_GOES=0;
isave_plot_DRIVER=0;
icoarsen=1;

coarsen_target='AMSRE';
coarsen_target='GOES';


%% Use iasc_desc to switch between day and night, or specific times
iasc_desc=1; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass
%iasc_desc=2; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass
iasc_desc=0; %Choose specific times

% 12th nov time difference to LST :-
% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)



%For 13th Nov the AMSRE overpass was at 06:30 UTC for the descending
%and 18:42 for the ascending, so pick 07:00 and 19:00 UM output times.

%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;
Ybins_DRIVER = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1; %bins used for AGU
Ybins_DRIVER = [-0.01 10.^[log10(30):0.15:log10(2500)]]; ichoose_Ybins=1;   %Revised wider bins at lower LWP
Ybins_DRIVER = [-0.01 10.^[log10(10):0.15:log10(2500)]]; ichoose_Ybins=1;   %Went to lower LWP bin for cumulative plots
Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;   %Went to lower LWP bin for cumulative plots
Ybins_DRIVER = [-0.01 10.^[log10(20):0.05:log10(2500)]]; ichoose_Ybins=1;   %

%% Times get selected here
sat_str_remss='AMSR-E'; %Default
isat_remss=1; %Default
if iasc_desc==1
    i_plot_goes=1; %Whether to plot GOES data or not

    if plot_remss==0
        time_select = datenum('13-Nov-2008 19:00'); %for UM - daytime
    else
        time_select = datenum('13-Nov-2008 13:00'); i_plot_goes=0; %f16 approx 8am LST on 13th (12:36 UTC)
          %Closest match to the 10:35am Terra CERES overpass
        isat_remss = 4; % 1=AMSRE, 2=f13, 3=f15, 4=f16, 5=f17, 6=TMI, 7=Windsat
        sat_str_remss=sats{isat_remss};  
        
        time_select = datenum('13-Nov-2008 19:00'); i_plot_goes=0; %AMSR-E at 13:54 LST (18:42 UTC)
          %Closest match to the 10:35am Terra CERES overpass
        isat_remss = 1; % 1=AMSRE, 2=f13, 3=f15, 4=f16, 5=f17, 6=TMI, 7=Windsat
        sat_str_remss=sats{isat_remss};         
        
%         time_select = datenum('12-Nov-2008 21:00'); i_plot_goes=0; %TMI 21:18 UTC overpass (16:30 LST for 76W) on12th
%             % Closest match to 14:00 LST Aqua Ceres overpass
%         isat_remss = 6; % 1=AMSRE, 2=f13, 3=f15, 4=f16, 5=f17, 6=TMI, 7=Windsat
%         sat_str_remss=sats{isat_remss};  
%         
%         time_select = datenum('12-Nov-2008 13:00'); i_plot_goes=0; %f16 approx 8am LST on 12th (12:48 UTC)
           % Closest match to 10am Terra on 12th for CERES
%         isat_remss = 4; % 1=AMSRE, 2=f13, 3=f15, 4=f16, 5=f17, 6=TMI, 7=Windsat
%         sat_str_remss=sats{isat_remss}; 
%         
%         i_plot_goes=0; %Combination of all the LWP plots that match the 3 Ceres overpasses
%         time_select = [datenum('12-Nov-2008 13:00') datenum('12-Nov-2008 21:00') datenum('13-Nov-2008 13:00')]; 
%         isat_remss = [4 6 4];        
%         sat_str_remss = 'REMSS';
        
    end
elseif iasc_desc==2
    i_plot_goes=0; %Whether to plot GOES data or not
    time_select = datenum('13-Nov-2008 07:00'); %for UM
else
    i_plot_goes=1;
    %Actual CERES times
    time_select = [datenum('12-Nov-2008 15:12') datenum('12-Nov-2008 19:06') datenum('13-Nov-2008 15:55')];
    %Except UM is only every 30 mins...
    time_select = [datenum('12-Nov-2008 15:00') datenum('12-Nov-2008 19:00') datenum('13-Nov-2008 16:00')];
%    time_select = [datenum('12-Nov-2008 15:00')];    
%    time_select = [datenum('12-Nov-2008 19:00')];    
%    time_select = [datenum('13-Nov-2008 16:00')];        
end

%Only using the time function for GOES at the moment
%time_choice.tol = 5/60/24;
time_choice.time_specific = time_select;
time_choice.find_nearest=1; %Means that it will just find the nearest matching time
    
    

load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150727T071106.mat';
%Using full domain
%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T083601.mat';
%Using partial domain to avoid boundary issues
%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T101619.mat';





goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc'; %GOES file
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131845.nc'; %GOES file



idat_driver=0;
clear fileUM xdat_import ydat_import line_pattern_DRIVER  line_colour_DRIVER marker_style_DRIVER


% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];

%FULL UM domain for 12th Nov
  LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; 




LAT_val_UM = LAT_val_DRIVER; LON_val_UM = LON_val_DRIVER;
LAT_val_GOES = LAT_val_DRIVER; LON_val_GOES = LON_val_DRIVER;



% -- For other option setting see inside the loops
pdf_type_driver='normal';
pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths (also normalises)

%--- Load and process the data

% For use later to coarsen the UM and GOES data
         d=abs(diff(gcm_Plat2D_AMSRE,[],1));
         dlat_AMSRE = meanNoNan(meanNoNan(d,1),1);
         d=abs(diff(gcm_Plon2D_AMSRE,[],2));
         dlon_AMSRE = meanNoNan(meanNoNan(d,1),1);
         

if plot_remss==1

%% ------------------------------
% ------ AMSRE & REMSS data --------
% ------------------------------
idat_driver=idat_driver+1;

line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0 0 1]; marker_style_DRIVER(idat_driver).m='*';

LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;
        
%------- Calculate the data to plot         
         
         
         
         %Making a "new" gcm_str to avoid using timeseries3 type data in
         %pdf2D
         gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
         gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
         gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
         gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;
         %Actually can't really define this it varies over the globe for
         %ASMRE
%         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre); 
            gcm_time_matlab_AMSRE2 = 0;
            gcm_time_UTC_AMSRE2 = 0;
            daynum_timeseries3_AMSRE2 = 1;
            modisyear_timeseries3_AMSRE2 = 1;
            
%            ilat_amsre = find(gcm_Plat2D_AMSRE>LAT_val(1) & gcm_Plat2D_AMSRE<LAT_val(2) & gcm_Plon2D_AMSRE>LON_val(1) & gcm_Plon2D_AMSRE<LON_val(2) )

Y_driver=[]; %set to null since will cat together the data
for itime_pdf=1:length(time_select)
    
    isat_remss = isat_remss(itime_pdf);

    if plot_remss==1
        load(load_file_remss)
        datetime_remss = date_remss + timeUTC_remss/24 ;
        it_remss = find( abs(datetime_remss(:,iasc_desc,isat_remss) - time_select(itime_pdf)) < tol_remss );

        if length(it_remss)>0
            Y_temp = squeeze(lwp_remss(:,:,it_remss,iasc_desc,isat_remss));
        else
            error('No REMSS time data for the time selected');
        end

        Y_driver = cat(1,Y_driver,Y_temp);

        %Making a "new" gcm_str to avoid using timeseries3 type data in
        %pdf2D
        gcm_Plat2D_AMSRE2 = gcm_Plat2D_REMSS;
        gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_REMSS;
        gcm_Plon2D_AMSRE2 = gcm_Plon2D_REMSS;
        gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_REMSS;
        %Actually can't really define this it varies over the globe for
        %ASMRE
        %         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre);
        gcm_time_matlab_AMSRE2 = 0;
        gcm_time_UTC_AMSRE2 = 0;
        daynum_timeseries3_AMSRE2 = 1;
        modisyear_timeseries3_AMSRE2 = 1;


        y_axis_vals = 'General y-axis no ilat simple'; %avoid the lat lon designation - not neede since the region is picked out
          %already in the .mat file

    else
        Y_driver = 1e3*squeeze(lwp_amsre(:,:,:,iasc_desc));
        y_axis_vals = 'General GCM-style';
    end


end

% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='AMSRE2';
        gcm_str='AMSRE2';
       
%        month_amsre = goes_month;
%        year_amsre = goes_year;

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
         Ybins = Ybins_DRIVER; ichoose_Ybins=1;
        
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
        

        
        ylabelstr='LWP (g m^{-2})';
        


        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
          axis1D = 'y';                                
                                
          logbin_norm = logbin_norm_driver;
          i_plot_norm=i_plot_norm_driver;
          i_div_bin_widths=i_div_bin_widths_driver;
          pdf_type = pdf_type_driver;
                                
%        gcm_str = gcm_str_last_loaded;        

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;
        datatype = 'gcm_data';        

        % --------- Override flags for watervap --------
        man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
        plotTimeHeightVap3
        close(gcf);
        waterVapourMay2005
        close(gcf);
        
        %store the PDF data
         switch pdf_type_driver
            case 'normal'
                ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
                xdat_import(idat_driver) = xdat_norm;
            case 'cumulative'
                ydat_import(idat_driver) = ydat_cum;
                xdat_import(idat_driver) = xdat_cum;
        end
       
        labs_import(idat_driver).l = sat_str_remss;
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;

        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;     

end


if i_plot_goes==1

%% ------------------------------
% ------ GOES data --------
% ------------------------------
idat_driver=idat_driver+1;


        
% %------- Calculate the data to plot
%          %read in the GOES data for the specific time
%           %read in the GOES data for the specific time
%          ioverride_goes = 1;
%          goes_action = 'load a particular file';        
%          read_GOES_vocals_netcdf_files
         
         
file_goes_multi_whole_dom = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat';
if iload_GOES==1
    load(file_goes_multi_whole_dom);
else
    %load the lat and lon files, though, as these get chagned
    load(file_goes_multi_whole_dom,'gcm_Plat2D_GOES','gcm_Plon2D_GOES');
end

LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;

%ilat_GOES = find(gcm_Plat2D_GOES >= LAT_val(1) & gcm_Plat2D_GOES <= LAT_val(2) ...
%& gcm_Plon2D_GOES>= LON_val(1) & gcm_Plon2D_GOES <= LON_val(2) );

[ilatB,ilonL] = getind_latlon_quick(gcm_Plat2D_GOES,gcm_Plon2D_GOES,LAT_val(1),LON_val(1),0.1);
[ilatT,ilonL] = getind_latlon_quick(gcm_Plat2D_GOES,gcm_Plon2D_GOES,LAT_val(2),LON_val(1),0.1);
[ilatB,ilonL] = getind_latlon_quick(gcm_Plat2D_GOES,gcm_Plon2D_GOES,LAT_val(1),LON_val(1),0.1);
[ilatB,ilonR] = getind_latlon_quick(gcm_Plat2D_GOES,gcm_Plon2D_GOES,LAT_val(1),LON_val(2),0.1);

ilat_GOES = min([ilatB ilatT]):max([ilatB ilatT]);
ilon_GOES = min([ilonL ilonR]):max([ilonL ilonR]);

% Calculate dlon and dloat for GOES data
d=diff(gcm_Plat2D_GOES,[],1);
dlat_GOES = meanNoNan(meanNoNan(d,1),1);
d=diff(gcm_Plon2D_GOES,[],2);
dlon_GOES = meanNoNan(meanNoNan(d,1),1);



%          if length(time_select_GOES)==1
%              itimes_GOES = find(abs(times_GOES_save-time_select_GOES)<5/60/24); %find nebarest one within 5 mins
%                        %Round to the nearest minute as sometimes get 18:59:59
% %                       time_str = datestr(round(time_out(itimes_GOES)*24*60)/24/60,'dd-mmm-yyyy HH:MM');
%          else
%              itimes(1) = find(abs(times_GOES_save-time_select_GOES(1))<5/60/24);
%              itimes(2) = find(abs(times_GOES_save-time_select_GOES(2))<5/60/24); 
%              times_GOES_freq = times_GOES_save(itimes(1)):freq_GOES/24:times_GOES_save(itimes(2));
%              itimes_GOES2 = NaN*ones(size(times_GOES_freq));
%              for i=1:length(times_GOES_freq)
%                  it = find(abs(times_GOES_save-times_GOES_freq(i))<5/60/24);
%                  if length(it)>0                    
%                      itimes_GOES2(i) = find(abs(times_GOES_save-times_GOES_freq(i))<5/60/24);
%                  end
%              end
% %             time_str = [datestr(round(time_out(itimes_GOES(1))*24*60)/24/60,'dd-mmm-yyyy HH:MM') ' to ' datestr(round(time_out(itimes_UM(end))*24*60)/24/60,'dd-mmm-yyyy HH:MM')];
%          end


        
        [out, gcm_time_matlab_GOES, itimes_GOES2] = get_time_range_of_array([],times_GOES_save,time_choice,0);
         
         clear Y_driver
         itime_DR2=0;
         for itime_DR=1:length(itimes_GOES2)

             itimes_GOES = itimes_GOES2(itime_DR);

             lwp_GOES = goes_LWP_multi{itimes_GOES}(ilat_GOES,ilon_GOES);
             
              
             if icoarsen==1  %
                 
                 switch coarsen_target
                     case 'AMSRE'
                         dlat_target = dlat_AMSRE;
                         dlon_target = dlon_AMSRE;
                     case 'GOES'
                         dlat_target = dlat_GOES;
                         dlon_target = dlon_GOES;

                 end

                 %average to the coarser resolution of AMSRE
                 N = ceil(abs(dlat_target/dlat_GOES));
                 M = ceil(abs(dlon_target/dlon_GOES));

                 lwp_GOES = reduce_matrix_subsample_mean(lwp_GOES,N,M);

                 %Am averaging Reff and Tau here rather than the LWP itself
                 %- hopefully is ok.
%                 goes_Reff = reduce_matrix_subsample_mean(goes_Reff,N,M);
%                 goes_Tau = reduce_matrix_subsample_mean(goes_Tau,N,M);


                    if itime_DR==1
                        gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
                        gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
                        %Work out the cell edges (as halfway between the centres)
                        [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);
                    end

             end
             
             if itime_DR==1
                 Y_driver = lwp_GOES;
             else
                 Y_driver = cat(3,lwp_GOES,Y_driver);
             end

         

         end %time loop

         
         

         
        
         
         
        
% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='GOES';
        gcm_str='GOES';
       
%        month_amsre = goes_month;
%        year_amsre = goes_year;

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
         Ybins = Ybins_DRIVER; ichoose_Ybins=1;
                 
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
        
               
        datatype = 'makeshift';
        ylabelstr='LWP (g m^{-2})';
        y_axis_vals = 'General y-axis no ilat simple';
        
        %Ybins = [-0.01 10.^[log10(10):0.15:log10(2500)]]; ichoose_Ybins=1;   %Went to lower LWP bin for cumulative plots
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
          axis1D = 'y';                                
                                
          logbin_norm = logbin_norm_driver;
          i_plot_norm=i_plot_norm_driver;
          i_div_bin_widths=i_div_bin_widths_driver;
          pdf_type = pdf_type_driver;
          

                                
%        gcm_str = gcm_str_last_loaded;        

        
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
        waterVapourMay2005
        close(gcf);
        
        %store the PDF data
         switch pdf_type_driver
            case 'normal'
                ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
                xdat_import(idat_driver) = xdat_norm;
            case 'cumulative'
                ydat_import(idat_driver) = ydat_cum;
                xdat_import(idat_driver) = xdat_cum;
        end
       
        labs_import(idat_driver).l = 'GOES';
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;    

else
        %Quick fix to get the linestyles to be the same as without GOES
    idat_driver=idat_driver+1;
    labs_import(idat_driver).l = ' ';
    xdat_import(idat_driver).x=NaN;
    ydat_import(idat_driver).y=NaN;
    
    Y_mean_import(idat_driver) = NaN;
    X_mean_import(idat_driver) = NaN;    

end

line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0.6 0.6 0.8]; marker_style_DRIVER(idat_driver).m='o';

%% ------------------------------
% ------ UM data --------
% ------------------------------
LAT_val = LAT_val_UM;
LON_val = LON_val_UM;
%N.B. - the regional screening for this is done in pdf2D_plot_commands
%(case 'UM LWP')

dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
idat=1;

for i=1:99
    flag{i}='';
end



% - 12th Nov case
%fileUM{idat} = '/xlhg-u/xlhgu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlhg-w/xlhgw_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-m/xlydm_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-aproc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;



%fileUM{idat} = '/xmmz-u/xmmzu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar(increased)';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xmmz-v/xmmzv_LWP_.pp.nc.mat';  labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xmmz-w/xmmzw_LWP_.pp.nc.mat';  labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  flag{idat} = 'load_mat';  pole_lat=70; pole_lon=284;
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;



%% 

fileUM{idat} = '/xmmz-u/xmmzu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
% fileUM{idat} = '/xmmz-v/xmmzv_LWP_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%     line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
% fileUM{idat} = '/xmmz-w/xmmzw_LWP_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%     line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;    
% %fileUM{idat} = '/xlhg-v/xlhgv_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
% fileUM{idat} = '/xmmz-x/xmmzx_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
% fileUM{idat} = '/xmmz-n/xmmzn_LWP_.pp.nc.mat'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
% %fileUM{idat} = '/xmmz-l/xmmzl_LWP_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
% %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0];
% %    marker_styleUM(idat).m='v'; idat=idat+1;
% 


for idat_UM=1:length(fileUM)
    idat_driver=idat_driver+1;
    
    line_pattern_DRIVER(idat_driver)=line_patternUM(idat_UM);  line_colour_DRIVER(idat_driver)=line_colourUM(idat_UM); marker_style_DRIVER(idat_driver)=marker_styleUM(idat_UM);    
    

    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    

    Y_driver=[];
    gcm_time_save=[];
    for itime_pdf=1:length(time_select)
        
        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        time = time_select(itime_pdf);
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
        
         vars_in.var = 'LWP';
         vars_in.flag = flag{idat_UM};
         vars_in.file_lwp =  [dirUM fileUM{idat_UM}]; 
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time;
    
     [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


     
%          vars_in.var = 'RWP';
%          vars_in.flag = flag{idat_UM};
%          vars_in.file_lwp = filename;
%          vars_in.file_rho = ''; %filename_rho;
%          vars_in.pole_lat = pole_lat;
%          vars_in.pole_lon = pole_lon;
%          vars_in.time_in = time_select;
%     
%      [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
%      
     
%        lwp = (lwp+rwp)*1e3;
        

        if icoarsen==1
            
            
            switch coarsen_target
                case 'AMSRE'
                    dlat_target = dlat_AMSRE;
                    dlon_target = dlon_AMSRE;
                case 'GOES'
                    dlat_target = dlat_GOES;
                    dlon_target = dlon_GOES;

            end
            
        %average to the coarser resolution of goes
       
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_target/dlat_UM));
        
       
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_target/dlon_UM));
        
        
        

        lwp_UM_n5 = reduce_matrix_subsample_mean(lwp,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        else
            lwp_UM_n5 = lwp;
        end
        
       %  Concatenate teh data from multiple times
         Y_driver = cat(3,Y_driver,lwp_UM_n5);
%% Important - also concatenate the times (or just make it an arrya of size
%% ntimes) since this determines the size of Plat3D in pdf2D_plot_commands
%% that sets iregion_lin - otherwise will miss the data from the times
%% except for t=1
        gcm_time_save = cat(1,gcm_time_save,gcm_time_matlab_UM);         
         
    end
    
    lwp_UM_n5 = Y_driver; %variable used by the pdf routine
    gcm_time_matlab_UM = gcm_time_save;   %Needed to set iregion_lin properly
        
% ----- Set various things

          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
%        titlenam_driver = ['LWP for ' time_str ' ' labs_import(idat).l];
%        units_str_plot = 'g m^{-2}';
         
%        mod_data_type='AMSRE';
        gcm_str_select='UM';
        gcm_str='UM';

       
        month_amsre = [1:length(time_matlab)];
        year_amsre = [1:length(time_matlab)];

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;
         
        %--- set some options for these particular plot loops
%        set_screening = {'none'};
%        modis_data_plot = 'Map of 2D data from outside driver script';
        i577 = 'MODIS_plot_UW';
        datatype = 'gcm_data';

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
        y_axis_vals = 'UM LWP';
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
        axis1D = 'y';                                
        
        logbin_norm = logbin_norm_driver;
        i_plot_norm=i_plot_norm_driver;
        i_div_bin_widths=i_div_bin_widths_driver;
        pdf_type = pdf_type_driver;
                                
%        gcm_str = gcm_str_last_loaded;        

        
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
        waterVapourMay2005
        close(gcf);
        
        %store the PDF data
         switch pdf_type_driver
            case 'normal'
                ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
                xdat_import(idat_driver) = xdat_norm;
            case 'cumulative'
                ydat_import(idat_driver) = ydat_cum;
                xdat_import(idat_driver) = xdat_cum;
         end
        
       
        labs_import(idat_driver).l = labs_UM(idat_UM).l;
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;    

end


        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------


%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'LWP PDFs';
xlab='Liquid Water Path (g m^{-2})';
ylab = ylab_import;
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=isave_plot_DRIVER;

%idate_ticks_fix=1;
%iaxis_square=0; %switch to make axis square
ichoose_styles=1;

line_pattern = line_pattern_DRIVER;  line_colour=line_colour_DRIVER; marker_style=marker_style_DRIVER;


%---  Main script to do plots and save
DRIVER_lineplot_watervap

switch pdf_type_driver
    case 'cumulative'
        set(gca,'xscale','log');
        set(gca,'xlim',[10 800]);
    case 'normal'
        %        set(gca,'xscale','log');
        %        set(gca,'xlim',[9 800]);

        switch iasc_desc
            case 1
                set(gca,'xscale','linear');
                set(gca,'xlim',[0 300]);
            case 2
                set(gca,'xscale','linear');
                set(gca,'xlim',[0 550]);
        end

end
        
        
        
        
        if isave_plot==1
            [datestr_now]=saveas_ps_fig_emf(gcf,[savename],'',0,1);
%            close(gcf);
        end
   
     

%    xdat_import(idat).x =







                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
