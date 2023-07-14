   
%% Note that pole_lon = 284 for this case (centred at 76W)    
day_or_night = 'day';
%day_or_night = 'night';
day_or_night = 'all';

irestrict_region_DRIVER = 1; %whether to do the regional selection when load data to save memory

fprintf(1,'\n*** N.B. - if says cannot load the file, then try closing Matlab and re-starting it (lack of memory I think) ***\n');
fprintf(1,'\n*** Also, will need to change the caxis colorbar to see the colours!!\n');

isave_plot=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

zlevs_file = '/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/zlevs_orig_L70_40';
[N_levs,z_levs,dz_levs] = read_zlevs_UM(zlevs_file);
%set to =1 for Nd
iasc_desc=1; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass
%iasc_desc=2; %Whether to use ascending (=1, daytime) and descending (=2, nighttime) overpass

%For 13th Nov the AMSRE overpass was at 06:30 UTC for the descending
%and 18:42 for the ascending, so pick 07:00 and 19:00 UM output times.

iplot_amsre=0; %set to zero for Nd - will keep code in for other parameters.
i_plot_goes=0;
i_plot_UM = 1;
iplot_radar=0; iload_ship_dat=1;

i_apply_min_dBZ=1; %Whether to apply the height dependent min dBZ that estimated from the ship obs to the model data
prc_val = 99.99; %percentile to use for approximating max dBZ
prc_val = 99.9;

% Choose normalization of PDFs
iuseYnorm_DRIVER=0;  % Each column will add up to frequency of 1.0
iuseXnorm_DRIVER=0;  % Each row will add up to frequency of 1.0
iuse_overall_norm_DRIVER=1; %Had this set to one - but changed to zero for testing for sample size effect.

dir_radar='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/VOCALS_obs_from_IanB/';
% .mat file that specifies the estimated min detectable dBZ from the ship
% obs (using the PDF)
min_dBZ_file = [dir_radar 'min_dBZ.mat'];

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

%LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov

%LAT_val_DRIVER = [-30.5 -6.44]; LON_val_DRIVER = [-95.93 -68.08]; %Wider view - whole GOES domain

LAT_val_DRIVER = [-20.5 -19.5]; LON_val_DRIVER = [-75.5 -74.5]; %1 degree region around the ship (20S, 75W)
res=0.25/16;LAT_val_DRIVER = [-20-res -20+res]; LON_val_DRIVER = [-75-res -75+res]; %Testing smaller regions around the ship (20S, 75W)



Xbins_DRIVER = [-0.01 10.^[log10(5):0.05:log10(300)]]; Xbins_driver = Xbins_DRIVER;
Xbins_DRIVER = [-1005:2:20]; Xbins_driver = Xbins_DRIVER;
%Ybins_DRIVER = [0:10:290 300:50:950 1000:500:6000];   %need to be replicated in pdf2d_plot_commands
%Ybins_DRIVER = [0:10:290 300:50:550];
Ybins_DRIVER = [0:110:1500]; %z_levs(1:28); %[0:20:1500];
%Ybins_DRIVER = [0:10:1500]; %z_levs(1:28); %[0:20:1500];
imax=20;
Ybins_DRIVER = (0.5*(z_levs(1:imax)+z_levs(2:imax+1))); %*** NEEDS to be of size [1 N], not [N 1] !!!
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

%Normalization
iuseYnorm = iuseYnorm_DRIVER;
iuseXnorm = iuseXnorm_DRIVER;
iuse_overall_norm = iuse_overall_norm_DRIVER;

%--- Load and process the data


if iplot_radar==1
%% ------------------------------
% ------ Ship radar data --------
% ------------------------------
LAT_val = LAT_val_UM;
LON_val = LON_val_UM;
%N.B. - the regional screening for this is done in pdf2D_plot_commands
%(case 'UM LWP')


ifile=1;

clear ship_radar_file base_time

% N.B. - base_time does not need to be set here below (is read from the
% file name now)

%Day 317 is 12th Nov, 2008.
ship_radar_file{ifile}='20083170000MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170100MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170200MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170300MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170400MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170500MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170600MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170700MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170800MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083170900MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171000MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171100MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171200MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171300MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;


ship_radar_file{ifile}='20083171400MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171500MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 15:00'); ifile=ifile+1;

ship_radar_file{ifile}='20083171600MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171700MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171800MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083171900MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083172000MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083172100MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083172200MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083172300MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;



%Day 318 is 13th Nov, 2008.
ship_radar_file{ifile}='20083180000MMCRMom.nc'; base_time{ifile} = datenum('13-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180100MMCRMom.nc'; base_time{ifile} = datenum('13-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180200MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180300MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180400MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180500MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180600MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180700MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180800MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083180900MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181000MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181100MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181200MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181300MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;


ship_radar_file{ifile}='20083181400MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181500MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 15:00'); ifile=ifile+1;

ship_radar_file{ifile}='20083181600MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181700MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181800MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083181900MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083182000MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083182100MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083182200MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;
ship_radar_file{ifile}='20083182300MMCRMom.nc'; base_time{ifile} = datenum('12-Nov-2008 14:00'); ifile=ifile+1;


%time_select = datenum('12-Nov-2008 15:00'); %Want to do a 30 min average either side of this

switch day_or_night
    case 'day'
        av_time = [ datenum('13-Nov-2008 17:00') datenum('13-Nov-2008 23:00') ];
    case 'night'
        av_time = [ datenum('13-Nov-2008 02:00')  datenum('13-Nov-2008 08:00') ];
    case 'all'
        %av_time = [datenum('12-Nov-2008 14:30') datenum('12-Nov-2008 15:30')];
        av_time = [datenum('13-Nov-2008 02:00') datenum('13-Nov-2008 23:59')];
        av_time = [datenum('12-Nov-2008 06:00') datenum('13-Nov-2008 23:59')];
end




if iload_ship_dat==1
    
    
    ref_ship=[]; time_ship=[]; sig_to_noise=[];
    
for idat_ship=1:length(ship_radar_file)
    idat_driver=idat_driver+1;
    
    %Calculate the day and time from 
    year_str = ship_radar_file{idat_ship}(1:4);
    day_str = ship_radar_file{idat_ship}(5:7);
    hour_str = ship_radar_file{idat_ship}(8:9);    
    [a,b] = date_from_day_of_year_func(str2num(day_str),str2num(year_str));
    base_time{idat_ship} = b + str2num(hour_str)/24;
    
%    line_pattern_DRIVER(idat_driver)=line_patternUM(idat_UM);  line_colour_DRIVER(idat_driver)=line_colourUM(idat_UM); marker_style_DRIVER(idat_driver)=marker_styleUM(idat_UM);    
    

        nc_ship = netcdf([dir_radar ship_radar_file{idat_ship}]);
        Heights_read = nc_ship{'Heights'}(:); % size is [nmodes Nheights]
        Heights = Heights_read(1,:); %Only the first mode seems to be used.
        %Size [1 120]
        
        ref_ship_read = nc_ship{'Reflectivity'}(:); %size [N_time N_heights]
        ref_ship = cat(1,ref_ship,ref_ship_read);
        
%        min_ref_read = nc_ship{'MinimumDetectableReflectivity'}(:); %These are all NaNs...
%        min_ref = cat(1,min_ref,min_ref_read);
        
        sig_to_noise_read = nc_ship{'SignalToNoiseRatio'}(:);
        sig_to_noise = cat(1,sig_to_noise,sig_to_noise_read);
        
        time_ship_read = nc_ship{'time_offset'}(:); 
        time_ship_read = base_time{idat_ship} + time_ship_read/3600/24;
%        time_base_ship = nc_ship{'base_time'}(:);
        time_ship = cat(1,time_ship,time_ship_read);
        
        

        


end

end


extra_title_info=['Ship radar dBZ ' datestr(av_time(1),'HH:MM') ' to ' datestr(av_time(2),'HH:MM') ' UTC'];

itimes_av = find( time_ship>=av_time(1) & time_ship<=av_time(2) );

ref_av = ref_ship(itimes_av,:);
times_av = time_ship(itimes_av);
time=meanNoNan(times_av,1);
sig_to_noise_av = sig_to_noise(itimes_av,:);

% Use the signal to noise ratio to remove likely clear-sky data.
min_sig=-14;
ref_av(sig_to_noise_av<min_sig)=-1000;
% Set the values in the first 3 levels to the min value since it is just
% ground clutter.
ref_av(:,1:3)=-1000;

     
siz=size(ref_av);     
height_dBZ = repmat(Heights,[siz(1) 1]);
     
      

%         if icoarsen==1
%             
%             dlat_target = dlat_GOES;
%             dlon_target = dlon_GOES;
%             
% %            dlat_target = dlat_AMSRE;            
% %            dlon_target = dlon_AMSRE;
%             
%         %average to the coarser resolution of goes
%        
%         d=diff(gcm_Plat2D_UM,[],1);
%         dlat_UM = meanNoNan(meanNoNan(d,1),1);
%         N = ceil(abs(dlat_target/dlat_UM));
%         
%        
%         d=diff(gcm_Plon2D_UM,[],2);
%         dlon_UM = meanNoNan(meanNoNan(d,1),1);
%         M = ceil(abs(dlon_target/dlon_UM));
%         
%         
%         
% 
%         dBZ_UM_n5 = reduce_matrix_subsample_mean(dBZ,N,M);
%         height_UM_n5 = reduce_matrix_subsample_mean(height_dBZ,N,M);
%         gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
%         gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
%         %Work out the cell edges (as halfway between the centres)
%         [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
%         
%         else
%             dBZ_UM_n5 = dBZ;
%             height_UM_n5 = height_dBZ;
%         end
        
        
% ----- Set various things

          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
%        titlenam_driver = ['LWP for ' time_str ' ' labs_import(idat).l];
%        units_str_plot = 'g m^{-2}';
         
%        mod_data_type='AMSRE';
        gcm_str_select='UM';
        gcm_str='UM';

       
        month_amsre = [1:length(times_av)];
        year_amsre = [1:length(times_av)];

        
        
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

        
        iminovr=1;
        mincovOvr = 0; %1e-4;
%        mincovOvr = 1e-2;
        imaxovr=1;
        maxcovOvr = 0.005; %1;
        

        
        logflag=0;
        dlogflag=0;        
        

                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'General GCM-style x-axis simple2'; %'Dummy data'; %dummy data
%        y_axis_vals = 'UM LWP';
%         y_axis_vals = 'General GCM-style';
        y_axis_vals = 'General y-axis no ilat simple';
        
         datatype = 'gcm_data';
         datatype = 'makeshift';         
         
%        ylabelstr = ['Nd for LWC.GT.' num2str(thresh_LWP_DRIVER) ' g m^{-3}, cm^{-3}'];
        ylabelstr = ['Height (m)'];
%        Ybins = Heights; ichoose_Ybins=1;
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;    
        Xbins = Xbins_DRIVER; ichoose_Xbins=1;          
        Y_driver = height_dBZ; %
        X_driver = ref_av; %
        xlabelstr = 'Radar reflectivity (dBZ)';
        
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
        man_choose_plotTimeHeight_graph2=1;
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


%        [LWP,ND]=meshgrid(Xbins,Ybins);
%        Precip=24*precip_rate_Wood_2008(LWP,ND);  %mm/day
%        cont_precip=10.^[-2:1:1];
%         [cs_precip,h_precip] = contour(LWP,ND,Precip,cont_precip,'g','linewidth',2);
%        clabel(cs_precip,h_precip,'fontsize',25,'color','g','rotation',0);


set(gca,'xlim',[-45 20]);
set(gca,'ylim',[0 1500]);
%caxis([0 0.2]);
caxis([0 0.003]);

% -- Estimate the min detectable dBZ for the ship (vs height)
min_dBZ = NaN*ones([1 length(height)]);
thresh_min_dBZ = 1e-4;
i=find(timesTH(1).t>-100);
for iz=1:length(height)
    ii=find(pdat(1).p(iz,i)>thresh_min_dBZ);
    if length(ii)>0
        min_dBZ(iz) = timesTH(1).t(i(ii(1)));
    end
end

ii = find(X_driver>-45);
prctile_minus45dBZ = prctile(X_driver(ii),prc_val)
prctile_ALL = prctile(X_driver(:),prc_val)
prc_val


save(min_dBZ_file,'min_dBZ','thresh_min_dBZ','pdat','timesTH','height');

        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
%            close(gcf);
        end



end %if iplot_radar










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
        
        [LWP,ND]=meshgrid(Xbins,Ybins);
        Precip=24*precip_rate_Wood_2008(LWP,ND);  %mm/day
        cont_precip=10.^[-2:1:1];
%        [cs_precip,h_precip] = contour(LWP,ND,Precip,cont_precip,'g','linewidth',2);
%        clabel(cs_precip,h_precip,'fontsize',25,'color','g','rotation',0);
        

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

load(min_dBZ_file);



% - 12th Nov case
%fileUM{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc.mat'; fileUM_Nd{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc';labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;

%fileUM{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Nd_fixed_act-10';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
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


%%

%fileUM{idat} = '/xmmz-u/xmmzu_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xmmz-v/xmmzv_Nd_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xmmz-w/xmmzw_Nd_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;%
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;    
fileUM{idat} = '/xmmz-x/xmmzx_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
%fileUM{idat} = '/xmmz-n/xmmzn_Nd_.pp.nc.mat'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;



for idat_UM=1:length(fileUM)
    idat_driver=idat_driver+1;
    
    line_pattern_DRIVER(idat_driver)=line_patternUM(idat_UM);  line_colour_DRIVER(idat_driver)=line_colourUM(idat_UM); marker_style_DRIVER(idat_driver)=marker_styleUM(idat_UM);    
    

    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    




        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
        
%          vars_in.var = 'Nd';
%          vars_in.flag = flag{idat_UM};
% %         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
%          vars_in.file_lwp =  [dirUM fileUM{idat_UM}];          
% %         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
%          vars_in.pole_lat = pole_lat;
%          vars_in.pole_lon = pole_lon;
%          vars_in.time_in = time_select;
%     
%      [var_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
%         time = [ datenum('13-Nov-2008 01:00') : 2/24 : datenum('13-Nov-2008 23:00') ];
         
%For ship PDFs for day and night took 6 hours either side of the 02:00 LST (07 UTC) and 14:00 LST (19 UTC) times that defined the min and max LWP from the satellite timeseries.
%But looking at the model LWP timeseries, the nightimte max is actually at
% 00 LST and the min at 15 LST. Probably sensible to only do 3 horus either side too, since
% otherwise start getting into the nex peak/trough.
% So, will do 21 to 3 LST on 12th/13th and 12 to 18 LST on 13th or :-
% 01:48 to 07:48 UTC (02-08 UTC) on 13th  and   16:48 to 22:48 UTC on 13th
% (17-23 UTC)
% But only have 5 hours of data beyond 19 UTC, so will use 5 hours either side         
        switch day_or_night
            case 'day'
                time = [ datenum('13-Nov-2008 17:00') : 0.5/24 : datenum('13-Nov-2008 23:00') ];
            case 'night'
                time = [ datenum('13-Nov-2008 02:00') : 0.5/24 : datenum('13-Nov-2008 08:00') ];
            case 'all'
                time = [ datenum('12-Nov-2008 06:00') : 0.5/24 : datenum('14-Nov-2008 00:00') ];                      
        end
        
        extra_title_info=['dbZ_vs_height_2D_histogram ' labs_UM(idat_UM).l ', ' datestr(time(1),'dd-mmm-yyyy HH:MM') ' to ' datestr(time(end),'dd-mmm-yyyy HH:MM') ' UTC '];                

         vars_in.var = 'dBZ';
         vars_in.flag = flag{idat_UM};
%         vars_in.file_Nd =  [dirUM fileUM_Nd{idat_UM}]; 
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','dBZ')];          
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
%         vars_in.time_in = time_select; 
%         vars_in.time_in = av_time; %Can choose a range of times? No - can specify multiple ones, but need to explicity state each
         vars_in.time_in = time;
         vars_in.irestrict_region = irestrict_region_DRIVER;
         vars_in.lat_restrict = LAT_val_DRIVER;
         vars_in.lon_restrict = LON_val_DRIVER;         
    
     [dBZ,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
     if length(size(vars_in.time_in))>1 & irestrict_region_DRIVER==0
         siz=size(dBZ);
         dBZ = reshape(dBZ,[siz(1) siz(2) sqrt(siz(3)) sqrt(siz(3))]);
         dBZ = permute(dBZ,[3 4 1 2]);

         %Returns an array of size [NZ NX*NY].
         %So need to replicate the heights to match

         %i_LWP = find(lwp<thresh_LWP_DRIVER);
         %var_UM(i_LWP)=NaN;

         siz=size(dBZ);

         height_dBZ = repmat(z_levs(1:siz(4)),[1 siz(1) siz(2) siz(3)]);
         height_dBZ = permute(height_dBZ,[2 3 4 1]);
     else
         siz=size(dBZ);
         height_dBZ = repmat(z_levs(1:siz(2)),[1 siz(1) siz(3)]);
         height_dBZ = permute(height_dBZ,[2 1 3]);       
     end
     
     if i_apply_min_dBZ==1
         load(min_dBZ_file)
         
         mid_z_min_dBZ = 0.5*(height(1:end-1)+height(2:end));
         
         %Permute to make the height axis the first index to make life
         %easier for setting all values < min_dBZ to -1000
         dBZ = permute(dBZ,[2 1 3]);

        for iz=1:length(min_dBZ)-1 %ignore the last height since this is didn't produce a sensible answer

            %find the nearest matching height - min_dBZ based on bin edges,
            %which were chosen to be the mid-points between model levels,
            %so that the central point of the bin edges is approx the model
            %levels. But is not exactly equal to them
            [minval iz2] = min(abs(z_levs - mid_z_min_dBZ(iz)));      
            N_ship = sum(pdat(1).p,2); %number in each height bin
            if N_ship(iz)==0
               dBZ(iz2,:)=NaN; %ignore levels where was no ship data at all
            elseif isnan(min_dBZ(iz)) %if it's NaN then assume no dBZ is detectable, but that was set to low dBZ
                dBZ(iz2,:)=-1000;
            else
                ilow = find(dBZ(iz2,:,:)<min_dBZ(iz));
                dBZ(iz2,ilow) = -1000;
            end
        end
        
        dBZ = permute(dBZ,[2 1 3]);
        
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
        
        
        

        dBZ_UM_n5 = reduce_matrix_subsample_mean(dBZ,N,M);
        height_UM_n5 = reduce_matrix_subsample_mean(height_dBZ,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        else
            dBZ_UM_n5 = dBZ;
            height_UM_n5 = height_dBZ;
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

        
        iminovr=1;
        mincovOvr = 0; %1e-4;
%        mincovOvr = 1e-2;
        imaxovr=1;
        maxcovOvr = 0.6; %1;
        

        
        logflag=0;
        dlogflag=0;        
        

        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
           
% --- PDF2d options        
        

%        screen_type = 'gcm_screening';

        if irestrict_region_DRIVER ==1
            x_axis_vals = 'General GCM-style x-axis simple2'; %'Dummy data'; %dummy data
            %        y_axis_vals = 'UM LWP';
            %         y_axis_vals = 'General GCM-style';
            y_axis_vals = 'General y-axis no ilat simple';
        else


            %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
            %        x_axis_vals = 'General GCM-style x-axis simple'; %'Dummy data'; %dummy data
            x_axis_vals = 'General reg lat lon screening x-axis simple 4D'; %'Dummy data'; %dummy data
            %        y_axis_vals = 'UM LWP';
            y_axis_vals = 'General GCM-style 4D';
        end
        
        datatype = 'gcm_data';
         
%        ylabelstr = ['Nd for LWC.GT.' num2str(thresh_LWP_DRIVER) ' g m^{-3}, cm^{-3}'];
        ylabelstr = ['Height (m)'];
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;    
        Xbins = Xbins_DRIVER; ichoose_Xbins=1;          
        Y_driver = height_UM_n5; %Nd
        X_driver = dBZ_UM_n5; %LWP
        xlabelstr = 'Radar reflectivity (dBZ)';
        
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
        man_choose_plotTimeHeight_graph2=1;
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


        [LWP,ND]=meshgrid(Xbins,Ybins);
        Precip=24*precip_rate_Wood_2008(LWP,ND);  %mm/day
        cont_precip=10.^[-2:1:1];
%         [cs_precip,h_precip] = contour(LWP,ND,Precip,cont_precip,'g','linewidth',2);
%        clabel(cs_precip,h_precip,'fontsize',25,'color','g','rotation',0);


xlims_dbz = [-45 20];
ylims_dbz = [0 1500];
set(gca,'ylim',ylims_dbz);
set(gca,'xlim',xlims_dbz);
%caxis([0 0.2]);
%caxis([0 0.01]);

savename=[savedir 'dBZ_vs_height_PDF_' labs_UM(idat_UM).l '_' upper(day_or_night)];

        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
%            close(gcf);
        end
        
        i=find(Xbins>=-45);
        frac_within_plot_bounds = sum(sum(qh(:,i)))./sum(qh(:))
        
        qh_norm = qh./sum(qh(:));
        
        
        
        %plot the sums along each axis
        figure
        sum_X_norm = sum(qh_norm(:,1:end-1),1);
        plot(mid_Xbins,sum_X_norm);
        set(gca,'xlim',xlims_dbz);
        xlabel(xlabelstr);
        
        
        cumsum_X = cumsum([0 sum_X_norm(i(1:end-1))]);
        bins_cumsum_X = Xbins(i);
        
     
        
        ii = find(X_driver>-45);
        prctile_minus45dBZ = prctile(X_driver(ii),prc_val)        
        prctile_ALL = prctile(X_driver(:),prc_val)
        prc_val
        

        figure
        %just do for dBZ>-45 (range of plot visible)
        i=find(Xbins>-45);
        plot(sum(qh_norm(1:end-1,i),2),mid_Ybins);
        set(gca,'ylim',ylims_dbz);
        ylabel(ylabelstr);

    
    

end % for loop

end


        
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

max_cont_val=3e-3;

%Find the max height of a given frequency contour
max_cont = max(pdat(1).p(:,10:512),[],2);
imax=find(max_cont>=max_cont_val);
% Report the lower edge of the corresponding bin
max_height_cont_val = z_levs(imax(1))


                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
