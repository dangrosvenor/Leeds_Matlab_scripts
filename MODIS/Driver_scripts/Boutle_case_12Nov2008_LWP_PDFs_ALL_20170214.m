% PDFs of LWP for all daytime times and several nightime ones in response
% to ref comments

% UM data for LWP is every 30 mins, so prob frequent enough for comparison to Remss overpassess    
    
%% Note that pole_lon = 284 for this case (centred at 76W)    

iload_GOES=1;

day_or_night='all';
day_or_night='night';
%day_or_night='day';

%For 13th Nov the AMSRE overpass was at 06:30 UTC for the descending
%and 18:42 for the ascending, so pick 07:00 and 19:00 UM output times.

iplot_remss=1; %set to zero for Nd - will keep code in for other parameters.
iplot_goes=1;
iplot_modis=0; %whether to include goes/modis Nd PDF


dlat_GOES = 1.0; %Default values - will be overridden if appropriate sat data provided
dlon_GOES = 1.0;


% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)


Ybins_DRIVER = [-0.01 10.^[log10(30):0.15:log10(2500)]]; ichoose_Ybins=1; %Used this for night plot I think (following the new figure after ref request)
Ybins_DRIVER = [-0.01 10.^[log10(30):0.25:log10(2500)]]; ichoose_Ybins=1; %Trying larger bins to get Poisson error down (partiuculary for lower bins for night)



switch day_or_night
    case 'day 13th Nov'
        time_select = datenum('13-Nov-2008 19:00'); %for UM - set one value or a range
        time_select = [datenum('13-Nov-2008 13:00') datenum('13-Nov-2008 21:00')]; %for UM - set one value or a range
    case 'all'
        time_select = [ datenum('31-Aug-2014 00:00')  datenum('05-Sep-2014 21:02') ];
    case 'night'
        iplot_goes=0;
        %Put in a few ranges that will contain the snapshots near the LWP
        %minima for the microwave instruments (see the LWP timseries)
        clear time_select
        it=1;      
        time_select{it} = [datenum('12-Nov-2008 03:00') datenum('12-Nov-2008 09:00')] - time_shift; it=it+1; %convert to UTC
        time_select{it} = [datenum('12-Nov-2008 20:00') datenum('13-Nov-2008 10:00')] - time_shift; it=it+1; %convert to UTC
        time_select{it} = [datenum('13-Nov-2008 18:30') datenum('13-Nov-2008 20:00')] - time_shift; it=it+1; %convert to UTC     
        
        Ybins_DRIVER = [-0.01 10.^[log10(30):0.15:log10(2500)]]; ichoose_Ybins=1; %Used this for night plot I think (following the new figure after ref request)
    case 'day'
        %Put in a few ranges that will contain the snapshots near the LWP
        %minima for the microwave instruments (see the LWP timseries)
        clear time_select
        it=1;      
        time_select{it} = [datenum('12-Nov-2008 10:00') datenum('12-Nov-2008 18:00')] - time_shift; it=it+1; %convert to UTC
        time_select{it} = [datenum('13-Nov-2008 10:00') datenum('13-Nov-2008 18:00')] - time_shift; it=it+1; %convert to UTC
%        time_select{it} = [datenum('13-Nov-2008 18:30') datenum('13-Nov-2008 20:00')] - time_shift; it=it+1; %convert to UTC  

        Ybins_DRIVER = [-0.01 10.^[log10(30):0.25:log10(2500)]]; ichoose_Ybins=1; %Trying larger bins to get Poisson error down (partiuculary for lower bins for night)

end



icoarsen=1;

idat_driver=0;
clear fileUM xdat_import ydat_import line_pattern_DRIVER  line_colour_DRIVER marker_style_DRIVER


% Smaller region to account for boundary inflow of LWP and
% spin-up during adatestr(time_sel+time_shift)dvection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
LAT_val_DRIVER = [-1e9 1e9]; LON_val_DRIVER = [-1e9 1e9];
LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];
LAT_val_DRIVER = [50.6 79.4]; LON_val_DRIVER = [-40.0 10.0];
LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov



thresh_LWP_DRIVER = 30; %g/m2 threshold to report GOES Nd
thresh_LWP_DRIVER = 5; %g/m2 threshold to report GOES Nd


LAT_val_UM = LAT_val_DRIVER; LON_val_UM = LON_val_DRIVER;
LAT_val_GOES = LAT_val_DRIVER; LON_val_GOES = LON_val_DRIVER;



% -- For other option setting see inside the loops
pdf_type_driver='normal';
%pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths (also normalises)
iarea_normalize_driver=1; %This is for the 2D PDFs I think.

%--- Load and process the data

if iplot_remss==1
%% ------------------------------
% ------ REMSS data --------
% ------------------------------
idat_driver=idat_driver+1;
line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0 0 1]; marker_style_DRIVER(idat_driver).m='*';


% Using full UM domain (the data in this file ies already cut down to the
% UM domain). It is the 0.25x0.25 deg gridded product. LWP array is
% size(lwp_remss) = 22    24     4     2     7 
% [nlat nlon ndays asc/desc nsats]
% size(timeUTC_remss) = 4     2     7
   load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T083601.mat';
   load(load_file_remss);    
   time_remss = date_remss + timeUTC_remss/24;
   
%    % File containing the 0.25 deg AMSRE grid - just used to get dlat and
%    % dlon, so could just specify this!

dlat_AMSRE = 0.25;
dlon_AMSRE = 0.25;

%    gcm_Plat2D_amsre_file = '/home/disk/eos1/d.grosvenor/AMSRE_gcm_Plat2D_etc.mat';
%    load(gcm_Plat2D_amsre_file);
% 
% LAT_val = LAT_val_GOES;
% LON_val = LON_val_GOES;
%         
% %------- Calculate the data to plot         
%          
%          % For use later to coarsen the UM data
%          d=abs(diff(gcm_Plat2D_AMSRE,[],1));
%          dlat_AMSRE = meanNoNan(meanNoNan(d,1),1);
%          d=abs(diff(gcm_Plon2D_AMSRE,[],2));
%          dlon_AMSRE = meanNoNan(meanNoNan(d,1),1);
%          
%          %Making a "new" gcm_str to avoid using timeseries3 type data in
%          %pdf2D
%          gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
%          gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
%          gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
%          gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;
%          %Actually can't really define this it varies over the globe for
%          %ASMRE


%         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre); 
            gcm_time_matlab_AMSRE2 = 0;
            gcm_time_UTC_AMSRE2 = 0;
            daynum_timeseries3_AMSRE2 = 1;
            modisyear_timeseries3_AMSRE2 = 1;
            
%            ilat_amsre = find(gcm_Plat2D_AMSRE>LAT_val(1) & gcm_Plat2D_AMSRE<LAT_val(2) & gcm_Plon2D_AMSRE>LON_val(1) & gcm_Plon2D_AMSRE<LON_val(2) )

Y_driver=[];
%Loop over all the requested time periods
for itime_periods=1:length(time_select)

        itf = find(time_remss >=time_select{itime_periods}(1) & time_remss <=time_select{itime_periods}(2) );
        [i1,i2,i3]=ind2sub(size(time_remss),itf);
        
        ii=1;
        iremove=[];
        tol_same = 30/60/24;
        for i=1:length(itf)
            for j=i+1:length(itf)
                tdiff = abs( time_remss(itf(i)) - time_remss(itf(j)) );
                if tdiff < tol_same 
                    iremove(ii)=i;
                    ii=ii+1;
                end
            end
        end
        if length(iremove)>0
            itf(iremove)=[];  %remove f15 if f13 and f15 both present
        end
        
%         imatch=0;
%         for i=1:length(i3)
%             switch sats{i3(i)}
%                 case {'SSMI-f13','SSMI-f15'}
%                     imatch=imatch+1;
%                     iremove = i;
%             end
%         end
%         if imatch==2
%             itf(iremove)=[];  %remove f15 if f13 and f15 both present
%         end
        dat = lwp_remss(:,:,itf); %in g/m2 already
%        Y_driver = 1e3*squeeze(lwp_amsre(:,:,:,iasc_desc));
        Y_driver = cat(1,Y_driver,dat(:));
         
end

% ----- Set various things

          
         
%        mod_data_type='AMSRE';

%         gcm_str_select='AMSRE2';
%         gcm_str='AMSRE2';
       
%        month_amsre = goes_month;
%        year_amsre = goes_year;

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
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
%        y_axis_vals = 'General GCM-style';       
        y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data
        
        datatype = 'makeshift';
        gcm_str='';
        
        ylabelstr='LWP (g m^{-2})';
        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;
        %Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1; %bins used for AGU
        %Ybins = [-0.01 10.^[log10(30):0.15:log10(2500)]]; ichoose_Ybins=1;   %Revised wider bins at lower LWP     
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;
         
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
          axis1D = 'y';                                
                                
          logbin_norm = logbin_norm_driver;
          i_plot_norm = i_plot_norm_driver;
          i_div_bin_widths=i_div_bin_widths_driver;
          pdf_type = pdf_type_driver;
          iarea_normalize = iarea_normalize_driver;
                                
%        gcm_str = gcm_str_last_loaded;        

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;
%        datatype = 'gcm_data';        

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
       
        labs_import(idat_driver).l = ['REMSS, \mu=' num2str(Y_mean_overall,'%.1f') ];
        xlab_import = xlab;
        ylab_import = ylab;
        
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;   
%        ioverride_savePDF=1;
%        save_1D_pdfs;

end


                            

                            


%% ------------------------------
% ------ GOES data --------
% ------------------------------
if iplot_goes==1
    
    
idat_driver=idat_driver+1;


LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;



        
%------- Calculate the data to plot
         %read in the GOES data for the specific time
          %read in the GOES data for the specific time
%          ioverride_goes = 1;
%          goes_action = 'load a particular file';        
%          read_GOES_vocals_netcdf_files
         
        file_goes_multi_whole_dom = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat';
        if iload_GOES==1
            load(file_goes_multi_whole_dom);
        else
            %Just load the lat and lon due to coarsening
            load(file_goes_multi_whole_dom,'gcm_Plat2D_GOES','gcm_Plon2D_GOES');
        end
        
        iregion_goes = find( gcm_Plat2D_GOES_nonan >= LAT_val_DRIVER(1) & gcm_Plat2D_GOES_nonan <= LAT_val_DRIVER(2) ...
            & gcm_Plon2D_GOES_nonan >= LON_val_DRIVER(1) & gcm_Plon2D_GOES_nonan <= LON_val_DRIVER(2) );
        %Actually, need this to be a regular 2d square for coarse graining.
        %So just just the limits of the region and assume is a rectangle in
        %incdices
        [ilat,ilon]=ind2sub(size(gcm_Plat2D_GOES_nonan),iregion_goes);
        ilat_region = min(ilat):max(ilat);
        ilon_region = min(ilon):max(ilon);
        
        %The other option would be to griddata the linear indices onto a
        %regular grid just using the restricted data, but will leave as a
        %square for now
        
        itimes_UM=[];
        for itime_periods=1:length(time_select)
            itf = find(times_GOES_save' >=time_select{itime_periods}(1) & times_GOES_save' <=time_select{itime_periods}(2) );
            itimes_UM = cat(1,itimes_UM,itf);
        end
         
         
         
         % For use later to coarsen the UM data
         d=diff(gcm_Plat2D_GOES,[],1);
         dlat_GOES = meanNoNan(meanNoNan(d,1),1);
         d=diff(gcm_Plon2D_GOES,[],2);
         dlon_GOES = meanNoNan(meanNoNan(d,1),1);
         
         for itime_DR=1:length(itimes_UM)
             
             isingle = itimes_UM(itime_DR);
             
             %Don't use if SZA>65 deg anywhere in the doimain
             if maxALL(sza_ALL{isingle}(ilat_region,ilon_region)) > 65 
                 continue
             end
             
             goes_LWP = goes_LWP_multi{isingle}(ilat_region,ilon_region);             
             
             
         
         if icoarsen==1

             %            dlat_target = dlat_GOES;
             %            dlon_target = dlon_GOES;

             dlat_target = dlat_AMSRE;
             dlon_target = dlon_AMSRE;

             %average to the coarser resolution of AMSRE

%              d=diff(gcm_Plat2D_UM,[],1);
%              dlat_UM = meanNoNan(meanNoNan(d,1),1);
             N = ceil(abs(dlat_target/dlat_GOES));


%             d=diff(gcm_Plon2D_UM,[],2);
%             dlon_UM = meanNoNan(meanNoNan(d,1),1);
             M = ceil(abs(dlon_target/dlon_GOES));



                %Am averaging Reff and Tau here rather than the LWP itself
                %- hopefully is ok. 
%             goes_Reff = reduce_matrix_subsample_mean(goes_Reff,N,M);
             goes_LWP_coarse = reduce_matrix_subsample_mean(goes_LWP,N,M);    
%              gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
%              gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
%              %Work out the cell edges (as halfway between the centres)
%              [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);
% 
% %          else
% %              lwp_UM_n5 = lwp;
         end
         
         
             if itime_DR==1
                 Y_driver = goes_LWP_coarse;
             else
                 Y_driver = cat(3,goes_LWP_coarse,Y_driver);
             end
         
         
         end
         
         if icoarsen==1

             N = ceil(abs(dlat_target/dlat_GOES));
             M = ceil(abs(dlon_target/dlon_GOES));

             gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
             gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
             %Work out the cell edges (as halfway between the centres)
             [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);
         end
        
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
%         x_axis_vals = 'Dummy data'; %dummy data
%         y_axis_vals = 'GOES LWP';
        
        x_axis_vals = 'Dummy data'; %dummy data
%        y_axis_vals = 'General GCM-style';       
        y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data
        
        datatype = 'makeshift';
        gcm_str='';
        
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
       
        labs_import(idat_driver).l = ['GOES, \mu=' num2str(Y_mean_overall,'%.1f') ];;
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;    

% else
%         %Quick fix to get the linestyles to be the same as without GOES
%     idat_driver=idat_driver+1;
%     labs_import(idat_driver).l = ' ';
%     xdat_import(idat_driver).x=NaN;
%     ydat_import(idat_driver).y=NaN;
%     
%     Y_mean_import(idat_driver) = NaN;
%     X_mean_import(idat_driver) = NaN;    

    
    line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0.6 0.6 0.8]; marker_style_DRIVER(idat_driver).m='o';
    
end






%% ------------------------------
% ------ MODIS data --------
% ------------------------------
if iplot_modis==1
    
idat_driver=idat_driver+1;


LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;



        
%------- Calculate the data to plot
         %read in the mockL3 MODIS data (1x1 deg)
         %CF>80 mock L3 data :-
         modis_Nd = load('/home/disk/eos8/d.grosvenor/mat_files_various/MODIS_Nd_Iceland_30Aug_10Sep_2014_CF_0.8_meanCTT_173_meanCTH_0.5_3.2km_SZA_65_meanTau10.mat');
         
         % For use later to coarsen the UM data
%          d=diff(gcm_Plat2D_GOES,[],1);
%          dlat_GOES = meanNoNan(meanNoNan(d,1),1);
%          d=diff(gcm_Plon2D_GOES,[],2);
%          dlon_GOES = meanNoNan(meanNoNan(d,1),1);
%          
%          sLAT = length(modis_Nd.MLAT);
%          sLON = length(modis_Nd.MLON);
%          
%          gcm_Plat2D_MODIS = repmat(modis_Nd.MLAT,[sLON 1]);
%          
%          [lat_out,lon_out] = get_edges_lat_lon(lat,lon)


         
         
         if icoarsen==99 %don't need to coarsen GOES

             %            dlat_target = dlat_GOES;
             %            dlon_target = dlon_GOES;

             dlat_target = dlat_AMSRE;
             dlon_target = dlon_AMSRE;

             %average to the coarser resolution of AMSRE

%              d=diff(gcm_Plat2D_UM,[],1);
%              dlat_UM = meanNoNan(meanNoNan(d,1),1);
             N = ceil(abs(dlat_target/dlat_GOES));


%             d=diff(gcm_Plon2D_UM,[],2);
%             dlon_UM = meanNoNan(meanNoNan(d,1),1);
             M = ceil(abs(dlon_target/dlon_GOES));



                %Am averaging Reff and Tau here rather than the LWP itself
                %- hopefully is ok. 
             goes_Reff = reduce_matrix_subsample_mean(goes_Reff,N,M);
             goes_Tau = reduce_matrix_subsample_mean(goes_Tau,N,M);    
             gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
             gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
             %Work out the cell edges (as halfway between the centres)
             [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);

%          else
%              lwp_UM_n5 = lwp;
         end
         
         
 
 
        

    

        
        %---  Run plot script and save
        MODIS_Nd_Iceland_DRIVER_template_1D_PDF_with_lat_lon_selection
        
        %store the PDF data
         switch pdf_type_driver
            case 'normal'
                ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
                xdat_import(idat_driver) = xdat_norm;
            case 'cumulative'
                ydat_import(idat_driver) = ydat_cum;
                xdat_import(idat_driver) = xdat_cum;
         end
        
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;   
       
        labs_import(idat_driver).l = ['MODIS, \mu=' num2str(Y_mean_overall,'%.1f') ];
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;
% else
%         %Quick fix to get the linestyles to be the same as without GOES
%     idat_driver=idat_driver+1;
%     labs_import(idat_driver).l = ' ';
%     xdat_import(idat_driver).x=NaN;
%     ydat_import(idat_driver).y=NaN;    

line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0.6 0.6 0.8]; marker_style_DRIVER(idat_driver).m='o';

end





%% ------------------------------
% ------ UM data --------
% ------------------------------
% Select cases from those in UM_case_select_runs

%UM_cases = '12th Nov case, as of May 2016 processing runs PLOTS multi-dirUM';
UM_cases = 'Iceland_9day_runs_Nov2016';
%UM_cases = 'Iceland_9day_runs_Nov2016_low_background_only';
UM_cases = '12th Nov case, as of May 2016';  %standard LWP has data every 30 mins


%% Script to get the UM run details by providing the run set name
%% Provide the case in UM_case_select_runs
UM_case_select_RUN  %runs UM_case_select_runs

LAT_val = LAT_val_UM;
LON_val = LON_val_UM;
%N.B. - the regional screening for this is done in pdf2D_plot_commands
%(case 'UM LWP')


for idat_UM=1:length(fileUM)
    idat_driver=idat_driver+1;
    
   if iscell(dirUM)==1
        dirUM_i = dirUM{idat_UM};
    else
        dirUM_i = dirUM;
    end
    
    line_pattern_DRIVER(idat_driver)=line_patternUM(idat_UM);  line_colour_DRIVER(idat_driver)=line_colourUM(idat_UM); marker_style_DRIVER(idat_driver)=marker_styleUM(idat_UM);    
    

    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    

        

        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        time = time_select;
                
         
         vars_in.var = 'LWP';
         vars_in.flag = 'load_mat';
         vars_in.file_lwp =  [dirUM_i remove_character(fileUM{idat_UM},'VAR_NAME','LWP') '.mat'];   %Keep as file_lwp for loading .mat files
%         vars_in.file_rho = [dirUM_i fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = []; %set to this for all times
    
     [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

    % i_LWP = find(lwp_ALL<thresh_LWP_DRIVER);
    % var_UM(i_LWP)=NaN;
    
    itimes_UM=[];
    for itime_periods=1:length(time_select)

        itf = find(time_matlab >=time_select{itime_periods}(1) & time_matlab <=time_select{itime_periods}(2) );
%        [i1,i2,i3]=ind2sub(size(time_remss),itf);

        itimes_UM = cat(1,itimes_UM,itf);
        
    end
     
%      if length(time_select)==1
%          itimes_UM = find(abs(time_out-time_select)<5/60/24); %find nebarest one within 5 mins
%          %Round to the nearest minute as sometimes get 18:59:59
%          time_str = datestr(round(time_out(itimes_UM)*24*60)/24/60,'dd-mmm-yyyy HH:MM');
%      else
%          % Sometimes get two zero times appearing - (find out why!) - this deals with this issue   
%          itimesA = find(abs(time_out-time_select(1))<5/60/24);
%          if length(itimesA)>1
%              for iA=1:length(itimesA)
%                  inan=find(isnan(maxALL(var_UM(itimesA(iA),:,:))));
%                  if inan==1 %if no NaNs use this time
%                  else
%                      itimes(1)=itimesA(iA);
%                      break
%                  end
%              end
%          else
%             itimes(1)=itimesA; 
%          end
%          itimesA = find(abs(time_out-time_select(2))<5/60/24);
%          if length(itimesA)>1
%              for iA=1:length(itimesA)
%                  inan=find(isnan(maxALL(var_UM(itimesA(iA),:,:))));
%                  if inan==1
%                  else
%                      itimes(2)=itimesA(iA);
%                      break
%                  end
%              end
%          else
%              itimes(2)=itimesA;
%          end
%          
%          itimes_UM = itimes(1):itimes(2);
%          time_str = [datestr(round(time_out(itimes_UM(1))*24*60)/24/60,'dd-mmm-yyyy HH:MM') ' to ' datestr(round(time_out(itimes_UM(end))*24*60)/24/60,'dd-mmm-yyyy HH:MM')];
%      end
     
     % Set this correctly for pdf* to work properly
         gcm_time_matlab_UM = gcm_time_matlab_UM(itimes_UM); 
         
         
 for itime_DR=1:length(itimes_UM) 
     filename = [dirUM fileUM{idat_UM}];
     itime_UM = itimes_UM(itime_DR);     
%     time = time_out(itime_UM);
     var_single = squeeze(var_UM(itime_UM,:,:));
     

        if icoarsen==1
            
%            dlat_target = dlat_GOES;
%            dlon_target = dlon_GOES;
            
            dlat_target = dlat_AMSRE;            
            dlon_target = dlon_AMSRE;
            
        %average to the coarser resolution of goes
       
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_target/dlat_UM));
        
       
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_target/dlon_UM));
        
        
        

        var_UM_n5 = reduce_matrix_subsample_mean(var_single,N,M);
     
        
        else
            var_UM_n5 = var_single;
        end
        
        
% ----- Set various things

          %Round to the nearest minute as sometimes get 18:59:59
%        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
%        titlenam_driver = ['LWP for ' time_str ' ' labs_import(idat).l];
%        units_str_plot = 'g m^{-2}';
         
%        mod_data_type='AMSRE';
        gcm_str_select='UM';
        gcm_str='UM';
        gcm_str='';

       
        month_amsre = [1:length(time_matlab)];
        year_amsre = [1:length(time_matlab)];

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
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
        y_axis_vals = 'UM LWP';
%         y_axis_vals = 'General GCM-style';
         y_axis_vals = 'General y-axis no ilat simple';
            
         datatype = 'gcm_data';
         datatype = 'makeshift';
         
%        ylabelstr = ['Nd for LWC.GT.0.05 g m^{-3}, cm^{-3}'];
        ylabelstr = 'LWP (g m^{-2})';
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;  
        
        if itime_DR==1
            Y_driver = var_UM_n5;
        else
            Y_driver = cat(3,var_UM_n5,Y_driver);
        end
        
        
 end   %time loop
        
         Y_driver=permute(Y_driver,[3 1 2]);
         lwp_UM_n5 = Y_driver;
         UM_overall_mean(idat_driver) = meanNoNan(Y_driver(:),1);
         
         %Coarsen the UM lat and lon
         if icoarsen==1
             gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
             gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
             %Work out the cell edges (as halfway between the centres)
             [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
         end
        
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
        
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;   
        
        labs_import(idat_driver).l = [labs_UM(idat_UM).l ', \mu=' num2str(Y_mean_overall,'%.1f') ];
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;

end %Loop over UM files


        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------


%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
%titlenam = ['Iceland Nd PDFs for LWP.GT.' num2str(thresh_LWP_DRIVER) ' g m^{-2}'];
titlenam = ['LWP PDFs'];
%xlab='Droplet Number Concentration (cm^{-3})';
xlab='Liquid Water Path (g m^{-2})';
ylab = ylab_import;
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=0;

%idate_ticks_fix=1;
%iaxis_square=0; %switch to make axis square
ichoose_styles=1;

line_pattern = line_pattern_DRIVER;  line_colour=line_colour_DRIVER; marker_style=marker_style_DRIVER;


%---  Main script to do plots and save
DRIVER_lineplot_watervap

    
% 
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% set(gca,'xlim',[2 3500]);
        
        
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
   
     

%    xdat_import(idat).x =







                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
