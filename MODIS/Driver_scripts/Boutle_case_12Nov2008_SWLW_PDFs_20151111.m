% PDFs of SW or LW
% Run read_CERES first

    
icut_for_CERES=1; %Cut out code not needed for the CERES analysis

icoarsen=1;

% Select set of UM runs - chosen in UM_case_select_runs function (via UM_case_select_RUN routine)
UM_cases = '12th Nov case, as of May 2016';
UM_cases = '12th Nov case, simple for Steve Abel';

%% Note that pole_lon = 284 for this case (centred at 76W)    


var_DRIVER='SW_up_TOA';
%var_DRIVER='SW_down_TOA';
%var_DRIVER='LW_up_TOA';

%set to =1 for Nd
iasc_desc=1; %Whether to use ascending (=1, daytime) or descending (=2, nighttime) overpass
%iasc_desc=2; %Whether to use ascending (=1, daytime) or descending (=2, nighttime) overpass

%For 13th Nov the AMSRE overpass was at 06:30 UTC for the descending
%and 18:42 for the ascending, so pick 07:00 and 19:00 UM output times.

iplot_ceres=1; 

tol_CERES_hrs=1;


time_range = 'Daytime CERES matches for SW';
%time_range = 'Nighttime CERES matches for LW';

%% Times get selected here
%if iasc_desc==1
%    i_plot_goes=1; %Whether to plot GOES data or not

switch time_range
    case 'Daytime CERES matches for SW'

% For the CERES overpasses, have set it to be 30 mins either side of these
% times. UM times have to be within one min I think.
    time_select = datenum('13-Nov-2008 19:00'); aqua_terra2={'Aqua'}; %Daytime, Aqua
    time_select = datenum('13-Nov-2008 15:00'); aqua_terra2={'Terra'}; %Daytime, Terra - only have 15:00 or 17:00 for the model
                                                                    %17 UTC does not work for the UM domain region (no CERES data)
                                %prob since this is 17-4.8 = 12:12 LST, whereas 15UTC is 10:12 LST (close
                                %to Terra time)
                                
%Combining several overpasses into one PDF (SW - daytime only)                                
%time_select = [datenum('12-Nov-2008 15:00') datenum('12-Nov-2008 19:00') datenum('13-Nov-2008 15:00')]; aqua_terra2={'Terra','Aqua','Terra'}; %Daytime, Terra - only have 15:00 or 17:00 for the model                                
% For the newer UM runs with output every 30 mins.
time_select = [datenum('12-Nov-2008 15:00') datenum('12-Nov-2008 19:00') datenum('13-Nov-2008 16:00')]; aqua_terra2={'Terra','Aqua','Terra'}; %Daytime, Terra - only have 15:00 or 17:00 for the model                                
%time_select = [datenum('12-Nov-2008 19:00')']; aqua_terra2={'Aqua'}; %Daytime, Terra - only have 15:00 or 17:00 for the model                                                                
%time_select = [datenum('12-Nov-2008 15:00')]; aqua_terra2={'Terra'}; %Daytime, Terra - only have 15:00 or 17:00 for the model                                


    case 'Nighttime CERES matches for LW'
% For LW can use nightime too - put nightime only here :-
% 12th Nov 10pm Terra, 12 Nov 2am Aqua, 13th Nov 2am Aqua
time_select = [datenum('13-Nov-2008 03:30') datenum('12-Nov-2008 07:30') datenum('13-Nov-2008 06:30')]; aqua_terra2={'Terra','Aqua','Aqua'}; %Daytime, Terra - only have 15:00 or 17:00 for the model                                
% 2am on 13th looks distorted - think is near the edge of the swath, so
% exclude this one. Actualy the 13th 03UTC looks a bit dodgy too...  :-
%time_select = [datenum('13-Nov-2008 03:00') datenum('12-Nov-2008 07:00')]; aqua_terra2={'Terra','Aqua'}; 

%time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC) =
  % = -0.2 = -4.8 hrs  (for around 76W)
  
    case 'All daytime (not matched to CERES)'
        %May want to weight to times when have most SW when looking at LWP,
        %etc. effect, though.
        
        time_select = [ datenum('12-Nov-2008 06:00')  datenum('14-Nov-2008 00:00') ];
        
    otherwise
%    i_plot_goes=0; %Whether to plot GOES data or not
    time_select = datenum('13-Nov-2008 07:00'); %for UM
end

 i_plot_goes=0;
 
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc'; %GOES file
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131845.nc'; %GOES file


% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)

idat_driver=0;
clear fileUM xdat_import ydat_import line_pattern_DRIVER  line_colour_DRIVER marker_style_DRIVER X_mean_import Y_mean_import


% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.

LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25]; %Top half of UM domain (12th Nov)
LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov


%Ybins_DRIVER = [0:10:290 300:50:950 1000:500:6000];   %need to be replicated in pdf2d_plot_commands



thresh_LWP_DRIVER = 30; %g/m2 threshold to report GOES Nd
thresh_LWP_DRIVER = 5; %g/m2 threshold to report GOES Nd


LAT_val_UM = LAT_val_DRIVER; LON_val_UM = LON_val_DRIVER;
LAT_val_GOES = LAT_val_DRIVER; LON_val_GOES = LON_val_DRIVER;
mean_lat_UM = mean(LAT_val_DRIVER);
mean_lon_UM = mean(LON_val_DRIVER);


% -- For other option setting see inside the loops
pdf_type_driver='normal';
%pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths (also normalises)

%--- Load and process the data

if iplot_ceres==1

%% ------------------------------
% ------ CERES data --------
% ------------------------------
idat_driver=idat_driver+1;

line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[1 0 0]; marker_style_DRIVER(idat_driver).m='o';

LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;

dlat_CERES = 20/distlatlon(mean_lat_UM,mean_lon_UM,mean_lat_UM+1,mean_lon_UM);
dlon_CERES = 20/distlatlon(mean_lat_UM,mean_lon_UM,mean_lat_UM,mean_lon_UM+1); %CERES footprint data is stated as being nominally at 20km
                     %Should also test with teh girdded data
        
%------- Calculate the data to plot 
if icut_for_CERES==0 %Cut out code not needed for the CERES analysis
    if icoarsen==1

        % For use later to coarsen the UM data
        d=abs(diff(gcm_Plat2D_AMSRE,[],1));
        dlat_AMSRE = meanNoNan(meanNoNan(d,1),1);
        d=abs(diff(gcm_Plon2D_AMSRE,[],2));
        dlon_AMSRE = meanNoNan(meanNoNan(d,1),1);

        %Making a "new" gcm_str to avoid using timeseries3 type data in
        %pdf2D
        gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
        gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
        gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
        gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;

    end
end
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
        
          aqua_terra = aqua_terra2{itime_pdf};
         
            itime_CERES = find(  abs(time_CERES - time_select(itime_pdf) ) < tol_CERES_hrs/24 ); %within 0.5hr either side
            if length(itime_CERES)>0
                lat_driver = lat_CERES(itime_CERES);
                lon_driver = lon_CERES(itime_CERES);
                switch var_DRIVER
                    case 'SW_up_TOA'
                        Y_temp = squeeze(SW_TOA_up_CERES(itime_CERES));
                        ylabelstr='SW upwards TOA (W m^{-2})';
                        Ybins_DRIVER = [0:40:900];
                        xlimits_DRIVER=[0 900];
                    case 'SW_down_TOA'
                        Y_temp = squeeze(SW_TOA_down_CERES(itime_CERES));
                        ylabelstr='SW downwards TOA (W m^{-2})';
                    case 'LW_up_TOA'
                        Y_temp = squeeze(LW_TOA_up_CERES(itime_CERES));
                        ylabelstr='LW upwards TOA (W m^{-2})';
                        Ybins_DRIVER = [200:2:320];
                        xlimits_DRIVER=[260 320];

                end

            else
                error('No time data within range!');
            end
            
            ylabelstr_DRIVER = ylabelstr;

            %Find the required location based on LAT_val_DRIVER and
            %LON_val_DRIVER
            iloc = find( lat_driver>LAT_val_DRIVER(1) & lat_driver<LAT_val_DRIVER(2) & lon_driver>LON_val_DRIVER(1) & lon_driver<LON_val_DRIVER(2) );            
            
            if length(iloc)==0
                error('No data within time and lat/lon range');
            end


            Y_driver = cat(1,Y_driver,Y_temp(iloc));
            time_CERES_overpass(itime_pdf) = mean(time_CERES(itime_CERES(iloc)));
            time_CERES_overpass_LST(itime_pdf) = mean(time_CERES(itime_CERES(iloc)))+time_shift;            
            
            
            
            
            
    end

% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='CERES';
        gcm_str='CERES';
        
       
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
        y_axis_vals = 'General y-axis no ilat simple';
        

        
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;
%        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1; %
    
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
          axis1D = 'y';                                
                                
          logbin_norm = logbin_norm_driver;
          i_plot_norm = i_plot_norm_driver;
          i_div_bin_widths=i_div_bin_widths_driver;
          pdf_type = pdf_type_driver;
                                
%        gcm_str = gcm_str_last_loaded;        

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;
%        datatype = 'gcm_data';     
        datatype = 'makeshift'; %Doing all of the time and lat/lon selection here

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
       
%        labs_import(idat_driver).l = ['CERES-' aqua_terra];
        labs_import(idat_driver).l = ['CERES'];        
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
if i_plot_goes==1

idat_driver=idat_driver+1;


LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;



        
%------- Calculate the data to plot
         %read in the GOES data for the specific time
          %read in the GOES data for the specific time
         ioverride_goes = 1;
         goes_action = 'load a particular file';        
         read_GOES_vocals_netcdf_files
         
         % For use later to coarsen the UM data
         d=diff(gcm_Plat2D_GOES,[],1);
         dlat_GOES = meanNoNan(meanNoNan(d,1),1);
         d=diff(gcm_Plon2D_GOES,[],2);
         dlon_GOES = meanNoNan(meanNoNan(d,1),1);
         
         
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
             goes_Reff = reduce_matrix_subsample_mean(goes_Reff,N,M);
             goes_Tau = reduce_matrix_subsample_mean(goes_Tau,N,M);    
             gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
             gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
             %Work out the cell edges (as halfway between the centres)
             [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);

%          else
%              lwp_UM_n5 = lwp;
         end
         
         
        
% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='GOES';
        gcm_str='GOES';
       
        month_amsre = goes_month;
        year_amsre = goes_year;
        
       

        
        
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
        y_axis_vals = 'GOES Nd';
        datatype = 'gcm_data';        
        
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
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;            
%        ioverride_savePDF=1;
%        save_1D_pdfs;
else
        %Quick fix to get the linestyles to be the same as without GOES
%    idat_driver=idat_driver+1;
%    labs_import(idat_driver).l = ' ';
%    xdat_import(idat_driver).x=NaN;
%    ydat_import(idat_driver).y=NaN;    

end

line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0.6 0.6 0.8]; marker_style_DRIVER(idat_driver).m='o';

%% ------------------------------
% ------ UM data --------
% ------------------------------
LAT_val = LAT_val_UM;
LON_val = LON_val_UM;
%N.B. - the regional screening for this is done in pdf2D_plot_commands
%(case 'UM LWP')

% dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
% dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
% idat=1;
% 
% for i=1:99
%     flag{i}='';
% end


% Selecting UM data
%   set flag{idat} to load_mat if loading in from a .mat file. Otherwise
%   can set to ''. 

% - 12th Nov case
%fileUM{idat} = ['/xlhg-u/xlhgu_' var_DRIVER '_.pp.nc']; fileUM_Nd{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = ''; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = ['/xlhg-v/xlhgv_' var_DRIVER '_.pp.nc']; fileUM_Nd{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = ''; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = ['/xlhg-w/xlhgw_' var_DRIVER '_.pp.nc']; fileUM_Nd{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = ''; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = ['/xlyd-m/xlydm_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlyd-w/xlydm_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-10-aproc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = ['/xlyd-z/xlydz_'  var_DRIVER '_.pp.nc']; labs_UM(idat).l = 'CASIM-Ndvar-10-SWNd';  flag{idat} = ''; fileUM_rho{idat} = '/xlyd-z/xlydz_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;


%% New runs with adjusted Nd for final draft version of paper (May 2016).

% fileUM{idat} = ['/xmmz-u/xmmzu_' var_DRIVER '_.pp.nc']; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = ''; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
%     line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
% fileUM{idat} = ['/xmmz-v/xmmzv_' var_DRIVER '_.pp.nc']; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = ''; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%     line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
% fileUM{idat} = ['/xmmz-w/xmmzw_' var_DRIVER '_.pp.nc']; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = ''; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%     line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;    
% %fileUM{idat} = ['/xlhg-v/xlhgv_' var_DRIVER '_.pp.nc']; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = ''; pole_lat=70; pole_lon=284;
% %        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
% fileUM{idat} = ['/xmmz-x/xmmzx_' var_DRIVER '_.pp.nc']; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = ''; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
% fileUM{idat} = ['/xmmz-n/xmmzn_' var_DRIVER '_.pp.nc']; labs_UM(idat).l = 'Old-mphys';  flag{idat} = ''; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
% %fileUM{idat} = ['/xmmz-l/xmmzl_' var_DRIVER '_.pp.nc']; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = ''; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
% %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0];
% %    marker_styleUM(idat).m='v'; idat=idat+1;



% --- New function to slelect the UM runs --- %
    UM_case_select_RUN



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

% For most cases will want to set the file_lwp variable        

switch var_DRIVER
    case 'LW_up_TOA'
         vars_in.var = 'LW_up_TOA';
         vars_in.flag = '';
%         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
%         vars_in.file_lwp =  [dirUM fileUM{idat_UM}];         
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'VAR_NAME',vars_in.var)];            
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time;
    
         [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
     
    case 'SW_up_TOA'

         vars_in.var = 'SW_up_TOA';
         vars_in.flag = '';
%         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
%         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','LWP')];          
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'VAR_NAME',vars_in.var)];                   
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time;
    
         [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

    case 'SW_down_TOA'

         vars_in.var = 'SW_down_TOA';
         vars_in.flag = '';
%         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
%         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','LWP')];       
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'VAR_NAME',vars_in.var)];            
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time_select;
    
         [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
         
         
end

%     i_LWP = find(lwp<thresh_LWP_DRIVER);
%     var_UM(i_LWP)=NaN;
     
      

        if icoarsen==1
            
%            dlat_target = dlat_GOES;
%            dlon_target = dlon_GOES;
            
%            dlat_target = dlat_AMSRE;            
%            dlon_target = dlon_AMSRE;
            
            dlat_target = dlat_CERES;
            dlon_target = dlon_CERES;
            
        %average to the coarser resolution of goes
       
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_target/dlat_UM));
        
       
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_target/dlon_UM));
        
        
        

        var_UM_n5 = reduce_matrix_subsample_mean(var_UM,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        else
            var_UM_n5 = var_UM;
        end
        
        
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
        
        month_amsre = [1:length(it_driver)];
        year_amsre = [1:length(it_driver)];        

        
        
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
%        y_axis_vals = 'UM LWP';
         y_axis_vals = 'General GCM-style';
         datatype = 'gcm_data';
         
        ylabelstr = ['Nd for LWC.GT.0.05 g m^{-3}, cm^{-3}'];
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;   
        
        
        Y_driver = cat(3,var_UM_n5,Y_driver);
%% Important - also concatenate the times (or just make it an arrya of size
%% ntimes) since this determines the size of Plat3D in pdf2D_plot_commands
%% that sets iregion_lin - otherwise will miss the data from the times
%% except for t=1
        gcm_time_save = cat(1,gcm_time_save,gcm_time_matlab_UM(it_driver));
        
    end
    
    gcm_time_matlab_UM = gcm_time_save;
        
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
        Y_mean_import(idat_driver) = Y_mean_overall;
        X_mean_import(idat_driver) = X_mean_overall;          
%        ioverride_savePDF=1;
%        save_1D_pdfs;

end


        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------


%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
%titlenam = ['SWLW PDFs for GOES LWP.GT.' num2str(thresh_LWP_DRIVER) ' g m^{-2}'];
titlenam = ['SWLW PDFs'];
titlenam = [''];
xlab=ylabelstr_DRIVER;
ylab = ylab_import;
xlims=1;
xlimits=xlimits_DRIVER;


izlim=0;
zmin=1500;
zmax=3000;

lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=0;

%idate_ticks_fix=1;
%iaxis_square=0; %switch to make axis square
ichoose_styles=1;

line_pattern = line_pattern_DRIVER;  line_colour=line_colour_DRIVER; marker_style=marker_style_DRIVER;

%Add the means into the legend label
for i=1:length(xdat_import)
   labs_import(i).l = [ labs_import(i).l ' \mu=' num2str(X_mean_import(i),'%2.1f')];
    
end



%---  Main script to do plots and save
DRIVER_lineplot_watervap

    

%set(gca,'xscale','log');
%set(gca,'yscale','log');
%set(gca,'xlim',[5 2000]);
        
    
savename = [savename 'CERES_' var_DRIVER '_PDF_' time_range];

        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
   
     

%    xdat_import(idat).x =







                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
