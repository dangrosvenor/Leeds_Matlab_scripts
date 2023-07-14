% PDFs of LWP
% If want to calculate LWP from qL and rho then use the routine
% calc_LWP_multi.m and save as .mat
% rather than calculating here. Then specify the .mat file here.
% For AMSRE at the moment will use the loaded file - fix later.
    %    Using multi_read_amsre_daily to load the amsre file
% For GOES it loads in the requested file using
    % read_GOES_vocals_netcdf_files
        
%% Note that pole_lon = 284 for this case (centred at 76W)    

titlenam_DRIVER = 'LWP PDFs ship'; %Becomes the savename for the figure too

time_period_ship = '13Nov 00:00 - 14Nov 00:00 UTC'; save_str_Werner='13Nov_to_14_Nov';
time_period_ship = '12Nov 06:00 - 14Nov 00:00 UTC'; save_str_Werner='All available minus 6hrs spin-up';
%time_period_ship = 'Daytime 13Nov 08:00 - 20:00'; save_str_Werner='Daytime';
%time_period_ship = 'Nighttime 12Nov 20:00 - 13Nov 08:00'; save_str_Werner='Nighttime';
%time_period_ship = 'Daytime2 13Nov 17:00 - 23:00'; save_str_Werner='Daytime';
%time_period_ship = 'Nighttime2 13Nov 02:00 - 08:00'; save_str_Werner='Nighttime';



i_plot_RHB=1;
iplot_UM=1;
i_calc_UM=1;
iplot_final=1; %Whether to do the final plot, or e.g. just the PDF comparison metric
i_just_calc_Dtot=0; %Flag to say whether to do pdfs for the whole domain (for calc of Werner, Dtot, param).
                    %So, if want to plot the LWP PDFs then set to zero.
include_RWP=0;          
i_old_ship_pdf_linestyles=0;

% Set the values for how many surrounding points to use here.
%i for lat, j for lon
i_inds = -1:1; j_inds = i_inds; %3x3 square with ship in middle
%i_inds = -2:2; j_inds = i_inds; %5x5 square
%        i_inds = -195:20:195; j_inds = i_inds;
%        i_inds = -195:5:195; j_inds = i_inds;
%        i_inds = -4:2:4; j_inds = i_inds;
%i_inds = 0; j_inds = 0; %i_inds;   %for the actual point - is quite noisy...

%% times to do the comparison for (model and ship use the same)
% ---- USING UTC TIME ----
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)
switch time_period_ship
    case '13Nov 00:00 - 14Nov 00:00 UTC'
        time_limit_RHB = [datenum('13-Nov-2008') datenum('14-Nov-2008')]; %Just look at second day for now to allow for model spin-up
    case '12Nov 06:00 - 14Nov 00:00 UTC'
        time_limit_RHB = [datenum('12-Nov-2008 06:00') datenum('14-Nov-2008 00:00')]; %Allow 6 hrs for model spin-up        
    case 'Daytime 13Nov 08:00 - 20:00'
        time_limit_RHB = datenum('13-Nov-2008 14:00') - time_shift + [-6 6]/24; %Daytime - 6hrs either side of 2pm LST on 13th Nov
    case 'Nighttime 12Nov 20:00 - 13Nov 08:00'
        time_limit_RHB = datenum('13-Nov-2008 02:00') - time_shift + [-6 6]/24; %Nighttime - 6hrs either side of 2am LST on 13th Nov
    case 'Daytime2 13Nov 17:00 - 23:00'
        time_limit_RHB = datenum('13-Nov-2008 20:00') + [-3 3]/24; %Daytime - 3hrs either side of 20 UTC 15:12 LST on 13th Nov
    case 'Nighttime2 13Nov 02:00 - 08:00'
        time_limit_RHB = datenum('13-Nov-2008 05:00') + [-3 3]/24; %Nighttime - 6hrs either side of 2am LST on 13th Nov        
end

Ybins_DRIVER=[0:10:50 60:10:300];
Ybins_DRIVER=[0:10:50 60:10:1000];
Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;   %As used for spatial PDFs
        
iasc_desc=1; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass
%iasc_desc=2; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass

%For 13th Nov the AMSRE overpass was at 06:30 UTC for the descending
%and 18:42 for the ascending, so pick 07:00 and 19:00 UM output times.


%Options for how to plot the points surrounding the ship grid box in the
%model - added these to get some idea of the spatial variability.
surrounding_points_option = 'plot separately';
surrounding_points_option = 'combine into one PDF';                       


i_plot_goes = 1;

%% Times get selected here
if iasc_desc==1
    i_plot_goes=1; %Whether to plot GOES data or not
    time_select = datenum('13-Nov-2008 19:00'); %for UM
else
    i_plot_goes=0; %Whether to plot GOES data or not
    time_select = datenum('13-Nov-2008 07:00'); %for UM
end



icoarsen=0;

idat_driver=0;
clear fileUM xdat_import ydat_import
%LAT_val_DRIVER = [-23.5 -16.44]; LON_val_DRIVER = [-85.93 -78.08]; %Smaller region to but out model edges
%LAT_val_DRIVER = [-24.5 -15.44]; LON_val_DRIVER = [-86.93 -77.08]; %GOES regin for UM comparison xkqk 26thOct POC
%LAT_val_DRIVER = [-24.5 -15.44]; LON_val_DRIVER = [-84 -77.08]; %Trying to match AMSRE and GOES domains


%LAT_val_DRIVER = [-21 -16.44]; LON_val_DRIVER = [-82 -78.08]; %Much smaller region in NW corner to avoid boundary issues.
  %This is what had set for AGu... but is this right?? Looks like the western edge of 
  %the domain is actually around -87... so for the smaller region in NW
  %would want -87 to -82....??
  
%LAT_val_DRIVER = [-21 -15.44]; LON_val_DRIVER = [-86.93 -82]; %Revised (for internal seminar) smaller region in NW corner 

%% 12th Nov case
% Domain runs from -78.93 to -73.08 W and -22.70 to -17.28 (based on the edges)
LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER = [-78.93 -73.08];
% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
%LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];


%The location of the ship (want to pick one point only) - 20S, 75W
LAT_val_DRIVER = -20; LON_val_DRIVER = -75;
%LAT_val_DRIVER = -21.4; LON_val_DRIVER = -75.9; %NOT the ship - Place where PDFs match quite well
%LAT_val_DRIVER = -20.8543; LON_val_DRIVER = -76.4911; %NOT the ship - Place where PDFs match quite well

%Pick a sub-region to do the analysis in - avoid the flow from the
%boundary, etc.
% Set as [] to do whole domain
LAT_region_bounds = []; LON_region_bounds = []; % Whole region
%LAT_region_bounds = [-21 -17.5362]; LON_region_bounds = [-78.9232 -73.3366]; %Allowing 0.25 in from the full domain boundary
 



LAT_val_UM = LAT_val_DRIVER; LON_val_UM = LON_val_DRIVER;




% -- For other option setting see inside the loops
pdf_type_driver='normal';
%pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths (also normalises)

%--- Load and process the data



                            
if i_plot_RHB==1

%% ------------------------------
% ------ Ship LWP data --------
% ------------------------------
idat_driver=idat_driver+1;

LAT_val = LAT_val_DRIVER;
LON_val = LON_val_DRIVER;
        
%--- Load and process the data
fileRHB='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/VOCALS_obs_from_IanB/vocals_rhb_lwp10min_v1.nc';
ncRHB=netcdf(fileRHB);
 
time=ncRHB{'jd_jdfrac'}(:); %"julian day + day fraction, 10-min res,time marker at beginning"
time(time>1e30)=NaN;
%convert into matlab time
time_matlab_RHB = datenum(2008,1,1) -1 + time; %(in UTC here rather than local)

lwp_RHB = ncRHB{'lwp'}(:);
lwp_RHB(lwp_RHB<-98.9)=NaN;
lwp_RHB(lwp_RHB>1e30)=NaN;

itime_RHB = find(time_matlab_RHB>=time_limit_RHB(1) & time_matlab_RHB<=time_limit_RHB(2) );

Y_driver = lwp_RHB(itime_RHB);
%xdat_import(idat_driver).x = time_matlab_RHB + time_shift; %
%labs_import(idat_driver).l = 'RHB ship';



         
         
        
% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='RHB';
        gcm_str='RHB';
       
        %month_amsre = goes_month;
        %year_amsre = goes_year;

        
        
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
        
        Ybins=Ybins_DRIVER; ichoose_Ybins=1; 
        ylabelstr='LWP (g m^{-2})';
        
        datatype = 'makeshift';
        
        
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
%        y_axis_vals = 'GOES LWP';
        y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data
        
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
       
        labs_import(idat_driver).l = ['Ship, \mu=' num2str(Y_mean_overall,'%.1f') ];
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;


        line_pattern_DRIVER(idat_driver).p='-';  line_colour_DRIVER(idat_driver).c=[0 0 1]; marker_style_DRIVER(idat_driver).m='*';


else
        %Quick fix to get the linestyles to be the same as without GOES
    idat_driver=idat_driver+1;
    labs_import(idat_driver).l = ' ';
    xdat_import(idat_driver).x=NaN;
    ydat_import(idat_driver).y=NaN;    

end



if iplot_UM==1
    


%% ------------------------------
% ------ UM data --------
% ------------------------------
LAT_val_DRIVER = LAT_val_UM;
LON_val_DRIVER = LON_val_UM;

dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
idat=1;

for i=1:99
    flag{i}='';
end



% - 12th Nov case
%fileUM{idat} = '/xlhg-u/xlhgu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; idat=idat+1;
%fileUM{idat} = '/xlhg-u/xlhgu_LWP_RWP_10min_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = ''; pole_lat=70; pole_lon=284; idat=idat+1;


% -- new runs with increased Nd
fileUM{idat} = '/xmmz-u/xmmzu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
fileUM{idat} = '/xmmz-v/xmmzv_LWP_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
fileUM{idat} = '/xmmz-w/xmmzw_LWP_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;    
fileUM{idat} = '/xmmz-x/xmmzx_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
fileUM{idat} = '/xmmz-n/xmmzn_LWP_.pp.nc.mat'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;


  
if i_just_calc_Dtot==1
   
    fileUM = fileUM(1);  %just do the first one if calcualating Dtot since will just show this in the paper
    
end
        

for idat_UM=1:length(fileUM)
    %idat_driver=idat_driver+1;   %Am looping over a few locations instead
    %below
    
    if i_calc_UM==1
    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    


        

        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        time = time_select;
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
        

         vars_in.var = 'LWP_10min';
%         vars_in.flag = flag{idat_UM};
         vars_in.flag = '';         
%         vars_in.file_lwp =  [dirUM fileUM{idat_UM}]; 
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'_LWP_.pp.nc.mat','_LWP_RWP_.pp.nc')];          
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = [];
    
     [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
     lwp = lwp*1e3;
     
     if include_RWP==1
         time = time_select;
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
        
         vars_in.var = 'RWP_10min';
         vars_in.flag = '';          
%         vars_in.flag = flag{idat_UM};
%         vars_in.file_lwp =  [dirUM fileUM{idat_UM}]; 
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'_LWP_.pp.nc.mat','_LWP_RWP_.pp.nc')];            
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = [];
    
         [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
         rwp = rwp*1e3;
         lwp = lwp + rwp;
         
     end
     
     it_UM = find(gcm_time_matlab_UM>=time_limit_RHB(1) &  gcm_time_matlab_UM<=time_limit_RHB(2));
    
        

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
        
        
        

        lwp_UM_n5 = reduce_matrix_subsample_mean(lwp,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        else
            lwp_UM_n5 = lwp;
        end
        
        [ilat_ship,ilon_ship] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,LAT_val_UM,LON_val_UM,0.1);
        
     if i_just_calc_Dtot==0
     

   

     else

         if length(LAT_region_bounds)==0
             ilat_botL=1; ilon_botL=1;
             ilat_topR=size(gcm_Plat2D_UM,1); ilon_topR=size(gcm_Plat2D_UM,2);
         else
             %Find the bottom left corner of the sub-region selected.
             [ilat_botL,ilon_botL] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,min(LAT_region_bounds),min(LON_region_bounds),0.1);
             %Find the top right corner of the sub-region selected.
             [ilat_topR,ilon_topR] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,max(LAT_region_bounds),max(LON_region_bounds),0.1);
         end
         
         [ilons_all,ilats_all] = meshgrid(ilat_botL:ilat_topR,ilon_botL:ilon_topR);
         
         ilin_all = sub2ind(size(gcm_Plat2D_UM),ilons_all,ilats_all);         
         
         Dtot_all = NaN*ones(size(gcm_Plat2D_UM)); %Make it the size of the full domain for consistency
         Ybins_mid = ( Ybins_DRIVER(2:end)+Ybins_DRIVER(1:end-1) ) / 2;
         
         for ilat=ilat_botL:ilat_topR
             for ilon=ilon_botL:ilon_topR
            
                 qh = ndHistc_run(lwp_UM_n5(it_UM,ilat,ilon), Ybins_DRIVER);
                 %Compare to the ship data          
                 [Dtot_all(ilat,ilon),D]=PDF_distance_Werner_1D(Ybins_mid,ydat_import(1).y,qh',1e6,'cumsum');
             end
             fprintf(1,'\nilat loop %d of %d...',ilat,length(ilat_botL:ilat_topR));
         end
         

         i_inds = []; j_inds = []; %set to avoid the loop below if just calculating Dtot
         
     end        
        
           

        num_surr_points = length(i_inds)*length(j_inds);
        
        iadd_pdf=0;
        idat_driver2=0;
        idat_driver_surr=0;
        for ii=i_inds
            for jj=j_inds
                idat_driver_surr=idat_driver_surr+1   %Counter to store all of the variability lines from surrounding points
                                                %- will combine all of
                                     %these into one PDF (or have option of plotting separately.
                                                
        switch surrounding_points_option
            case 'plot separately'
                Y_driver = lwp_UM_n5(it_UM,ilat_ship+ii,ilon_ship+jj);
                iadd_pdf=1;
                lab_str_extra = [' (' num2str(ii) ',' num2str(jj) ')'];
            case 'combine into one PDF'
                if idat_driver_surr==1
                   Y_driver = []; 
                end
                if idat_driver_surr==num_surr_points  %if is the last one then plot the combined pdf
                    iadd_pdf=1;
                end
                dat = lwp_UM_n5(it_UM,ilat_ship+ii,ilon_ship+jj);                
                Y_driver = cat(1,Y_driver,dat(:)); %Concatenate the data together for one PDF.
                lab_str_extra = '';
        end
        
        datatype = 'makeshift';
        
        
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
         pdf2D_defaults  %defaults for pdf2D_plot_commands
         
        
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
        y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data    
        
        Ybins=Ybins_DRIVER; ichoose_Ybins=1; 
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

        
if iadd_pdf==1  %Depending on whether the separate or combined PDF option is chosen           
        idat_driver2 = idat_driver2 + 1;
        
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
                ydat_import2(idat_driver2) = ydat_norm; %Use the non-cumulative PDF data
                xdat_import2(idat_driver2) = xdat_norm;
            case 'cumulative'
                ydat_import2(idat_driver2) = ydat_cum;
                xdat_import2(idat_driver2) = xdat_cum;
         end
         
         Y_mean_import(idat_driver) = Y_mean_overall;
         X_mean_import(idat_driver) = X_mean_overall;
        
       
        labs_import2(idat_driver2).l = [labs_UM(idat_UM).l lab_str_extra ', \mu=' num2str(Y_mean_overall,'%.1f') ];
        xlab_import2 = xlab;
        ylab_import2 = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;

end



            end
        end
        
    end  %i_calc_UM
        
        %Store them all (done like this for compatability with the option
        %to pick out the min and max below, although prob won't use this)
        for i=1:idat_driver2
            idat_driver=idat_driver+1;
        
            ydat_import(idat_driver) = ydat_import2(i);
            xdat_import(idat_driver) = xdat_import2(i);
            labs_import(idat_driver).l = labs_import2(i).l;    
            
            line_pattern_DRIVER(idat_driver)=line_patternUM(idat_UM);  line_colour_DRIVER(idat_driver)=line_colourUM(idat_UM); marker_style_DRIVER(idat_driver)=marker_styleUM(idat_UM);

        end
        xlab_import = xlab_import2;
        ylab_import = ylab_import2;
        
%Actually, prob not a good idea to highlight the min and max of the PDFs since PDFs need to be taken as a whole really since
%lack of LWP at low values would be compensated by more at higher values in
%one curve, but if pick some parts from one curve and some from another
%then is not really a fair comparison

%         %Pick out the actual closest position (ilat,ilon)
%         idat_driver=idat_driver+1;
%         imid = floor(idat_driver2/2); %Middle position (only use odd numbers in total, e.g. 9 or 25 probably).
%         ydat_import(idat_driver) = ydat_import2(imid);
%         xdat_import(idat_driver) = xdat_import2(imid);
%         xlab_import(idat_driver) = xlab_import2(imid);
%         ylab_import(idat_driver) = ylab_import2(imid);  
%         labs_import(idat_driver).l = [labs_UM(idat_UM).l];
%         
%         ymin=1e9*ones(size(ydat_import2(imid).y)); ymax=zeros(size(ydat_import2(imid).y));
%         %Now find min and max from all the surroudnging points
%         for i=1:idat_driver2
%             ymin = min(ymin,ydat_import2(i).y);
%             ymax = max(ymax,ydat_import2(i).y);            
%         end
%         
%         %Store min and max
%         idat_driver=idat_driver+1;
%         ydat_import(idat_driver).y = ymin
%         xdat_import(idat_driver) = xdat_import2(imid);
%         xlab_import(idat_driver) = xlab_import2(imid);
%         ylab_import(idat_driver) = ylab_import2(imid);  
%         labs_import(idat_driver).l = [labs_UM(idat_UM).l ' min'];
%         
%         %Store min and max
%         idat_driver=idat_driver+1;
%         ydat_import(idat_driver).y = ymax
%         xdat_import(idat_driver) = xdat_import2(imid);
%         xlab_import(idat_driver) = xlab_import2(imid);
%         ylab_import(idat_driver) = ylab_import2(imid);  
%         labs_import(idat_driver).l = [labs_UM(idat_UM).l ' max'];        
        
        
        

end

end

        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------
if iplot_final==1

%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = [titlenam_DRIVER ' for ' time_period_ship];
xlab='Liquid Water Path (g m^{-2})';
ylab = ylab_import;
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=0;

%idate_ticks_fix=1;
%iaxis_square=0; %switch to make axis square

if i_old_ship_pdf_linestyles==1

ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
       istyle=1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; line_widths(istyle).l = 3; istyle=istyle+1;       
       
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; line_widths(istyle).l = 2; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; line_widths(istyle).l = 2; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
       
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; line_widths(istyle).l = 3; istyle=istyle+1;
       
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.7 0.7 0]; marker_style(istyle).m='s'; line_widths(istyle).l = 2; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[0.7 0.7 0]; marker_style(istyle).m='s'; line_widths(istyle).l = 2; istyle=istyle+1;       


else
    ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
    line_pattern = line_pattern_DRIVER;  line_colour=line_colour_DRIVER; marker_style=marker_style_DRIVER;
    
       
end
       
%---  Main script to do plots and save
DRIVER_lineplot_watervap


switch surrounding_points_option
    case 'plot separately'
        uistack(h(6).h,'top');  %Move the 6th line to be on top
end

uistack(h(1).h,'top'); 
set(gca,'xlim',[0 550]);
set(gca,'ylim',[0 0.02]);


        
        
        
        
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
        
end
   
     

%    xdat_import(idat).x =



%% Calculate a score metric for the UM PDFs vs the ship PDF
clear Dtot
if iplot_UM==1 & i_plot_RHB==1
    
    if i_just_calc_Dtot==1

        Dtot = Dtot_all;
        [a,b]=minALL(Dtot);
        lat_minW = gcm_Plat2D_UM(b(1),b(2));
        lon_minW = gcm_Plon2D_UM(b(1),b(2));
                
        i=findstr(fileUM{1},'/');
        savedir_UM_Dtot = [dirUM fileUM{1}(1:i(2))];
        savefile_Dtot = [savedir_UM_Dtot 'Dtot_all_' save_str_Werner '.mat'];        
        save(savefile_Dtot,'Dtot','lat_minW','lon_minW');

        
        

    else

        switch pdf_type_driver
            case 'cumulative'
                error('Cumulative PDF selected - use a normal PDF');
            case 'normal'
                for i=2:length(ydat_import)
                    fprintf(1,'i=%i, ',i);
                    [Dtot(i-1),D]=PDF_distance_Werner_1D(xdat_import(1).x,ydat_import(1).y,ydat_import(i).y,1e6);
                end

        end



    %Example of finding the min Werner value and getting the lat
    %and lon, and the index in the overall array
    [a,b]=min(Dtot);
    [J,I]=ind2sub([sqrt(length(Dtot)) sqrt(length(Dtot))],b);
    %N.B. need to put the j-index first here since the j loop earlier was the inner loop
    %So if Dtot is rearranged into a square matrix called D
    %(equivalent to what the ind2sub above would refer to)
    %then the j values from the loop will fill along D(j,1) first
    % (since for a matrix D(1:n) starts expands along D(:,1) first - i.e.
    % along columns, or along the leftmost index
    lat_minW = gcm_Plat2D_UM(ilat_ship+i_inds(I),ilon_ship+j_inds(J));
    lon_minW = gcm_Plon2D_UM(ilat_ship+i_inds(I),ilon_ship+j_inds(J));
    
    end

end

                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
