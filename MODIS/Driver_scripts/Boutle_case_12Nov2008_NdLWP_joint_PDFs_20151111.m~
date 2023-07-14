% PDFs of Nd for GOES or UM
% For GOES it loads in the requested file using
    % read_GOES_vocals_netcdf_files
    
    
%% Note that pole_lon = 284 for this case (centred at 76W)   



iswap_xy_DRIVER=1; %swaps the x and y axes (plots LWP vs Nd)

joint_plot_case =  'LWP vs Nd';
%joint_plot_case =  'LWC vs Nd (3D)';

%set to =1 for Nd
iasc_desc=1; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass
%iasc_desc=2; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass

%For 13th Nov the AMSRE overpass was at 06:30 UTC for the descending
%and 18:42 for the ascending, so pick 07:00 and 19:00 UM output times.

iplot_amsre=0; %set to zero for Nd - will keep code in for other parameters.
i_plot_goes=1;
i_plot_UM = 0;


%% Times get selected here
time_select = datenum('13-Nov-2008 19:00'); %for UM - daytime
%time_select = datenum('12-Nov-2008 03:00'); %

iplot_precip_contours=0;

% if iasc_desc==1
%     i_plot_goes=1; %Whether to plot GOES data or not
%     time_select = datenum('13-Nov-2008 19:00'); %for UM - daytime
% else
%     i_plot_goes=0; %Whether to plot GOES data or not
%     time_select = datenum('13-Nov-2008 07:00'); %for UM
% end

%26th Oct
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261545.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261715.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file

%13th Nov
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131215.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131445.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131645.nc'; %GOES file
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131845.nc'; %GOES file



%Other days
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811121845.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811141845.nc'; %GOES file
%%only have 12, 13 and 14th for November
%Looks like have most days for October
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810311845.nc'; %GOES file
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810191515.nc'; %   Hook case from Rhea's paper
%goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810191845.nc'; %   At low LWP point


time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)

icoarsen=0;

idat_driver=0;
clear fileUM xdat_import ydat_import line_pattern_DRIVER  line_colour_DRIVER marker_style_DRIVER


% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25]; %Partial UM domain (top half)

LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov

%LAT_val_DRIVER = [-30.5 -6.44]; LON_val_DRIVER = [-95.93 -68.08]; %Wider
%view - whole GOES domain
%LAT_val_DRIVER = [-20 -19.9]; LON_val_DRIVER =[-78.93 -73.08]; %Picking slice, 12th Nov

thresh_LWP_DRIVER = 30; %g/m2 threshold to report GOES Nd
thresh_LWP_DRIVER = 5; %g/m2 threshold to report GOES Nd
thresh_LWP_DRIVER = -5; %g/m2 threshold to report GOES Nd

LAT_val_UM = LAT_val_DRIVER; LON_val_UM = LON_val_DRIVER;
LAT_val_GOES = LAT_val_DRIVER; LON_val_GOES = LON_val_DRIVER;


Xbins_DRIVER = [-0.01 10.^[log10(5):0.05:log10(350)]]; Xbins_driver = Xbins_DRIVER;
%Ybins_DRIVER = [0:10:290 300:50:950 1000:500:6000];   %need to be replicated in pdf2d_plot_commands
%Ybins_DRIVER = [0:10:290 300:50:550];
Ybins_DRIVER = [0:10:750];
Ybins_DRIVER = [0:10:350];

% -- For other option setting see inside the loops
pdf_type_driver='normal';
%pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths (also normalises)

%--- Load and process the data

iswap_xy = iswap_xy_DRIVER;

if iplot_amsre==1

%% ------------------------------
% ------ AMSRE data --------
% ------------------------------
idat_driver=idat_driver+1;

line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0 0 1]; marker_style_DRIVER(idat_driver).m='*';

LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;
        
%------- Calculate the data to plot         
         
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
         %Actually can't really define this it varies over the globe for
         %ASMRE
%         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre); 
            gcm_time_matlab_AMSRE2 = 0;
            gcm_time_UTC_AMSRE2 = 0;
            daynum_timeseries3_AMSRE2 = 1;
            modisyear_timeseries3_AMSRE2 = 1;
            
%            ilat_amsre = find(gcm_Plat2D_AMSRE>LAT_val(1) & gcm_Plat2D_AMSRE<LAT_val(2) & gcm_Plon2D_AMSRE>LON_val(1) & gcm_Plon2D_AMSRE<LON_val(2) )
         
         Y_driver = 1e3*squeeze(lwp_amsre(:,:,:,iasc_desc));

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
        y_axis_vals = 'General GCM-style';
        
        ylabelstr='LWP (g m^{-2})';
        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;
        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1; %bins used for AGU
        Ybins = [-0.01 10.^[log10(30):0.15:log10(2500)]]; ichoose_Ybins=1;   %Revised wider bins at lower LWP     
        
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
       
        labs_import(idat_driver).l = 'AMSRE';
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;

end


                            
if i_plot_goes==1

%% ------------------------------
% ------ GOES data --------
% ------------------------------
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
%        clim_min=0;
        clim_min= 1e-4;       
        iset_max_clim=1;
%        clim_max=200;
        clim_max=1;  
        
         iminovr=1;
        mincovOvr = 1e-4;
%        mincovOvr = 1e-2;
        imaxovr=1;
        maxcovOvr = 1;
        

        
        logflag=1;
        dlogflag=0;
        
        isave_plot=0;
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
        
        

        extra_title_info=[' ' datestr(gcm_time_matlab_GOES) '(' datestr(gcm_time_matlab_GOES+time_shift,'HH:MM') ' LST)'];


%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'GOES LWP x-axis'; %
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
        
        % Plot contours of precip

        [LWP,ND]=meshgrid(Xbins,Ybins);
        Precip=24*precip_rate_Wood_2008(LWP,ND);  %mm/day
        
        if iswap_xy_DRIVER==1
            ycont = LWP'; xcont = ND'; zcont = Precip';
        else
            xcont = LWP; ycont = ND; zcont = Precip;
        end
        
        cont_precip=10.^[-2:1:1];
        [cs_precip,h_precip] = contour(xcont,ycont,zcont,cont_precip,'g','linewidth',2);
        clabel(cs_precip,h_precip,'fontsize',25,'color','g','rotation',0);
        

%        
        
%         close(gcf);
%         waterVapourMay2005
%         close(gcf);
%         
%         %store the PDF data
%          switch pdf_type_driver
%             case 'normal'
%                 ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
%                 xdat_import(idat_driver) = xdat_norm;
%             case 'cumulative'
%                 ydat_import(idat_driver) = ydat_cum;
%                 xdat_import(idat_driver) = xdat_cum;
%         end
%        
%         labs_import(idat_driver).l = 'GOES';
%         xlab_import = xlab;
%         ylab_import = ylab;
% %        ioverride_savePDF=1;
% %        save_1D_pdfs;
else
        %Quick fix to get the linestyles to be the same as without GOES
%     idat_driver=idat_driver+1;
%     labs_import(idat_driver).l = ' ';
%     xdat_import(idat_driver).x=NaN;
%     ydat_import(idat_driver).y=NaN;    

end

%line_pattern_DRIVER(idat_driver).p= '-';  line_colour_DRIVER(idat_driver).c=[0.6 0.6 0.8]; marker_style_DRIVER(idat_driver).m='o';

if i_plot_UM == 1

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
%fileUM{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;


fileUM{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Nd_fixed_act-10';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlyd-y/xlydy_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-new';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-y/xlydy_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlyd-p/xlydp_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-NoRain';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-p/xlydp_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-z/xlydz_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-SWNd';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-z/xlydz_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-q/xlydq_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-fixed-act-NoRain';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-q/xlydq_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0.0 0.0]; marker_styleUM(idat).m='<'; idat=idat+1;
%fileUM{idat} = '/xlyd-r/xlydr_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-fixed-act-NoRain-Less-Mixing';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-r/xlydr_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0.0]; marker_styleUM(idat).m='<'; idat=idat+1;
%fileUM{idat} = '/xlyd-s/xlyds_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-NoRain-Less-Mixing';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-s/xlyds_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 1]; marker_styleUM(idat).m='<'; idat=idat+1;
%fileUM{idat} = '/xmmz-d/xmmzd_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Nd-fixed_act-NoRain-Less-Mixing-CF=1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-d/xmmzd_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 1]; marker_styleUM(idat).m='<'; idat=idat+1;
%fileUM{idat} = '/xmmz-h/xmmzh_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Nd-fixed_act-NoRain-Less-Mixing-CF=1-NoSed';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-d/xmmzh_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 1]; marker_styleUM(idat).m='<'; idat=idat+1;




for idat_UM=1:length(fileUM)
    idat_driver=idat_driver+1;
    
    line_pattern_DRIVER(idat_driver)=line_patternUM(idat_UM);  line_colour_DRIVER(idat_driver)=line_colourUM(idat_UM); marker_style_DRIVER(idat_driver)=marker_styleUM(idat_UM);    
    

    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    


        extra_title_info=[' ' labs_UM(idat_UM).l ', ' datestr(time_select) '(' datestr(time_select+time_shift,'HH:MM') ' LST)'];        

        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        time = time_select;
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]

%%
switch joint_plot_case
    case 'LWP vs Nd'
        
        
        Xbins_DRIVER = [-0.01 10.^[log10(5):0.05:log10(300)]]; Xbins_driver = Xbins_DRIVER;
        %Ybins_DRIVER = [0:10:290 300:50:950 1000:500:6000];   %need to be replicated in pdf2d_plot_commands
        %Ybins_DRIVER = [0:10:290 300:50:550];
        Ybins_DRIVER = [0:10:750];
        
         vars_in.var = 'Nd';
         vars_in.flag = flag{idat_UM};
%         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
         vars_in.file_lwp =  [dirUM fileUM{idat_UM}];   
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','Nd_no_min_qL')];            
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time_select;
    
     [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


         vars_in.var = 'LWP';
         vars_in.flag = flag{idat_UM};
%         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','LWP')];          
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time_select;
    
     [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

     i_LWP = find(lwp<thresh_LWP_DRIVER);
     var_UM(i_LWP)=NaN;
     
     
    case 'LWC vs Nd (3D)'
        
        i_slice_only=0;
        

        Xbins_DRIVER = [-0.01 10.^[log10(0.02):0.05:log10(2)]]; Xbins_driver = Xbins_DRIVER;
%        Xbins_DRIVER = [-0.01:0.02:2]; Xbins_driver = Xbins_DRIVER;        
        %Ybins_DRIVER = [0:10:290 300:50:950 1000:500:6000];   %need to be replicated in pdf2d_plot_commands
        %Ybins_DRIVER = [0:10:290 300:50:550];
        Ybins_DRIVER = [0:10:950];
        
        vars_in.var = 'Nd_3D';
        vars_in.flag = ''; %flag{idat_UM};
        %         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}];
        vars_in.file_lwp =  [dirUM fileUM{idat_UM}];
        vars_in.file_lwp =  remove_character(vars_in.file_lwp,'.mat','');        
        %         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
        vars_in.pole_lat = pole_lat;
        vars_in.pole_lon = pole_lon;
        vars_in.time_in = time_select;

        [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

        
        var_UM = var_UM/1e6; %convert to #/mg (aking to per cc)

        vars_in.var = 'qL';
        vars_in.flag = ''; %flag{idat_UM};
        %         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}];
        vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','qL')]; 
        vars_in.file_lwp =  remove_character(vars_in.file_lwp,'.mat','');        
        %         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
        vars_in.pole_lat = pole_lat;
        vars_in.pole_lon = pole_lon;
        vars_in.time_in = time_select;

        [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

        lwp = lwp*1e3; %Convert to g/kg from kg/kg
        
%        ilwc = find(lwp<0.05);
%        lwp(ilwc)=NaN;
%        var_UM(ilwc)=NaN;
        
        
%        i_LWP = find(lwp<thresh_LWP_DRIVER);
%        var_UM(i_LWP)=NaN;

%         if i_slice_only==1
%             %Need to convert this to use a 2D case in PDF2D RE the lat lon
%             %Plat3D etc.
%             %Can just use a narrow lat range at the top of this script
%             %anwyay.
%             LAT_slice = -18;
%             LON_slice = -75;
%             [ilat,ilon] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,LAT_slice_DRIVER,LON_slice,0.1);
% 
%             nz=20; %max height index to plot up to
% 
%             lwp = lwp(1:nz,ilat,ilon);
%             var_UM= var_UM(1:nz,ilat,ilon);
%         end
        
        
        
end
     
      

        if icoarsen==1
            
            dlat_target = dlat_GOES;
            dlon_target = dlon_GOES;
            
%            dlat_target = dlat_AMSRE;            
%            dlon_target = dlon_AMSRE;
            
        %average to the coarser resolution of goes
       
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_target/dlat_UM));
        
       
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_target/dlon_UM));
        
        
        

        var_UM_n5 = reduce_matrix_subsample_mean(var_UM,N,M);
        lwp_UM_n5 = reduce_matrix_subsample_mean(lwp,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        else
            var_UM_n5 = var_UM;
            lwp_UM_n5 = lwp;
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

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
        
        %--- set some options for these particular plot loops
%        set_screening = {'none'};
%        modis_data_plot = 'Map of 2D data from outside driver script';
        i577 = 'MODIS_plot_UW';

        iset_min_clim=0;
        clim_min=0;
        iset_max_clim=1;
        clim_max=200;

        

        

        
     
        
        isave_plot=0;
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'General GCM-style x-axis simple'; %'Dummy data'; %dummy data
%        y_axis_vals = 'UM LWP';
         y_axis_vals = 'General GCM-style';
         datatype = 'gcm_data';
         
switch joint_plot_case
    case 'LWP vs Nd'
        iarea_normalize=1; %Normalise each 2D bin by the area of X and Y bin
       % Probably important when using unevely spaced bins.         
        
        logflag=0;
        dlogflag=0;           
         
        ylabelstr = ['Nd for LWP.GT.' num2str(thresh_LWP_DRIVER) ' g m^{-2}, cm^{-3}'];
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;    
        Xbins = Xbins_DRIVER; ichoose_Xbins=1;          
        Y_driver = var_UM_n5; %Nd
        X_driver = lwp_UM_n5; %LWP
        xlabelstr = 'LWP (g m^{-2})';
        
        iminovr=1;
        mincovOvr = 1e-4;
%        mincovOvr = 1e-2;
        imaxovr=1;
        maxcovOvr = 1e-2;
        
    case 'LWC vs Nd (3D)'
        logflag=0;
        dlogflag=0; 
        iarea_normalize=1; %Normalise each 2D bin by the area of X and Y bin
       % Probably important when using unevely spaced bins. 
        
        ylabelstr = ['Nd (mg^{-1})'];
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;    
        Xbins = Xbins_DRIVER; ichoose_Xbins=1;          
        Y_driver = var_UM_n5; %Nd
        X_driver = lwp_UM_n5; %LWP
        xlabelstr = 'LWC (g kg^{-1})';
        
        %Use the time dimension for the height to make the screening for
        %lat and lon work
        gcm_time_matlab_UM = [1:size(var_UM_n5,1)]';
        
        iminovr=1;
        mincovOvr = 1e-7;
%        mincovOvr = 1e-2;
        imaxovr=1;
        maxcovOvr = 1e-4;   
        maxcovOvr = 1e-3;           
        
        
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
%         close(gcf);
%         waterVapourMay2005
%         close(gcf);
%         
%         %store the PDF data
%          switch pdf_type_driver
%             case 'normal'
%                 ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
%                 xdat_import(idat_driver) = xdat_norm;
%             case 'cumulative'
%                 ydat_import(idat_driver) = ydat_cum;
%                 xdat_import(idat_driver) = xdat_cum;
%          end
%         
%        
%         labs_import(idat_driver).l = labs_UM(idat_UM).l;
%         xlab_import = xlab;
%         ylab_import = ylab;
% %        ioverride_savePDF=1;
% %        save_1D_pdfs;

    if iplot_precip_contours==1
        [LWP,ND]=meshgrid(Xbins,Ybins);
        Precip=24*precip_rate_Wood_2008(LWP,ND);  %mm/day
        cont_precip=10.^[-2:1:1];
        [cs_precip,h_precip] = contour(LWP,ND,Precip,cont_precip,'g','linewidth',2);
        clabel(cs_precip,h_precip,'fontsize',25,'color','g','rotation',0);
    end


        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end

end % for loop

end


savedir = '/home/disk/eos1/d.grosvenor/modis_work/lwp_vs_Nd_plots/';
savename = [savedir 'lwp_vs_nd_2D_histogram'];
%saveas_ps_fig_emf(gcf,[savename],'',0,1);
        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------

% 
% %--- run the file to set up the defaults
% watervap_defaults
% 
% %--- set some options for this particular plot
% graph=0; %graph choice in watervap
% titlenam = ['Nd PDFs for GOES LWP.GT.' num2str(thresh_LWP_DRIVER) ' g m^{-2}'];
% xlab='Droplet Number Concentration (cm^{-3})';
% ylab = ylab_import;
% xlims=0;
% xlimits=[0 100];
% 
% izlim=0;
% zmin=1500;
% zmax=3000;
% 
% lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
% nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
% 
% isave_plot=0;
% 
% %idate_ticks_fix=1;
% %iaxis_square=0; %switch to make axis square
% ichoose_styles=1;
% 
% line_pattern = line_pattern_DRIVER;  line_colour=line_colour_DRIVER; marker_style=marker_style_DRIVER;
% 
% 
% %---  Main script to do plots and save
% DRIVER_lineplot_watervap
% 
%         
        
        
        
%         if isave_plot==1
%             saveas_ps_fig_emf(gcf,[savename],'',0,1);
%             close(gcf);
%         end
   
     

%    xdat_import(idat).x =







                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
