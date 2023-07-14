try

%Have added experimental option to not include certain lines in teh legend.
%Just set the labs(i).l for those lines to NaN

savenotes_filedir='/home/disk/eos1/d.grosvenor/notes/';
    
if ~exist('man_choose_water_graph')
    graph = 0; % Case where the data comes from another script (xdat_import, ydat_import, etc.)
    graph=40;  %ice_budget
    graph=46;  %timeseries_dan
    %graph=45;  %av vert wind
    %graph=61; %wind profile


    %graph=57; %emm and LEM comparisons, inc. adiabatic LWC

    %graph=56; %vapour eq model plots

    %graph=3; %SF4 vapour data
    %graph=48; %total water before/after plots
    %graph=482; %mean TTL as function of ice sat MR
    %graph=483; %tot and mean TTL from TwoD slices (2d/3d comparison)
    %graph=484; %tot and vap changes from TwoD (low/high CCN)
    %graph=485; %tot and vap changes from TwoD (3D/2D)
    %graph=4852; %tot and vap changes from TwoD (2D)
    %graph=486; %tot and vap changes from TwoD (1km/2km)
    %graph=487; %vap from several 3d cases for vapour and water
    %graph=4872; %vap from two 3d cases for temp change
    %graph=4873; %min vap from two 3d cases
    %graph=4874; %vap/tot from 3d CCN cases
    %graph=4875; %min vap from two 3d cases (general comparison for any runs)
    %graph=4876; %mean vap and tot change for N 3D runs (27Jul07)



    %graph=488; %vapour plots from TwoD (averaged over the 500km either side of y=0)
    %graph=44; %maxw
    %graph=53; %potemp perturbations (percentile plots to check temps)
    %graph=54; %ice sat mr LEM
    %graph=62; %aerosol distribution
    %graph=63; %vapour ozone
    %graph=64; %DMI ozonesondes
    %graph=65;
    %graph=17; %top adown cumulative (weighted) means of total water, etc.
    %graph=66; %3d min
    %graph=67; %tracerflux
    %graph=68; %eddy heat flux contribution
    %graph=69; %mass flux of updraught
    %graph=70; %max ice no.
    %graph=71; %average speeds for Manchester 10k
    %graph=72; %running time breakdown
    %graph=73; %MAC3 plots
    %graph=74; %WRF max profiles
    %graph=744; %WRF maean profiles
    %graph=75; %WRF wind speed
    %graph=76; %WRF wind dir
    %graph=77; %max LEM Q-wfield profiles
%    graph=78; %general WRF plots
    %graph=79; %WRF constant lat/lon slice
%    graph=80; %WRF lat/lon slice of mean contributions to melt
    %graph=81; %WRF distance/lon vs. variable plots
    %graph=82; %houghton/Smith hydraulic plots
    %graph=83; % WRF microphysics for ice heteorogeneous freezing, Bigg's immersion droplet freezing and contact nucleation
    %graph=84; % Froude number as function of depth (thi) for continuous stratification
%    graph=85; % IN concentrations
    %graph=86; % Profiles different locations from the cross section
    
%    graph=87; % CAS/CIP size distributions
%    graph=88; % Particle separation distributions
%    graph=89; % Aircraft profiles

    %graph=399; %DMI sounding and ice sat and LEM vapour profile

%    graph=90; % ACPIM vs Matlab profiles
%     graph=91; %CAS/hotwire analysis
%     graph=92; %CAS/hotwire analysis - overall means of all flights
%     graph=93; %CAS number analysis vs mode size
%     graph=94;
% graph=95;  %condensation rate (dLWC/dz)
% graph=952;  %condensation rate (dLWC/dz) Rob vs Albrecht
% graph=953; %latent heat temperature variation

%  ***   1D mean  *****
 graph=96;  %1D means from 2D histogram PDF data
 graph = 966;
%  ********************
 
% graph=99; % N points along x from 2D PDF
 
% graph=97;  %1D graph from 2D histogram data

%  ***   1D PDF  *****
 graph = 977; %new 1D PDF from 2D histo data - can choose either axis
%  ******************** 

%graph=98;% 1D PDF of the data produced by pdf2D_plot_commands_L2 (just for x-axis)

%graph=100; % Mode along x from 2D PDF
%graph=101; %nd vs szal from constant latitudes
%graph=102; %dnd/dSZA vs sza from constant latitudes
%graph=103; %Overall mean of dNd/dSZA for each SZA bin
%graph=104;
%graph=105; %1D PDFs for lat/lon cells from mock L3 data
%graph=106; %gaussian Reff or Tau
%graph=107; %plots along MPACE flight track
%graph=108; %profiles from MPACE flight data
%graph=109;  %X means from PDF for different seasons as stored by
%monthly_means_from_plot_global.m
%graph=110;  %N datapoints from PDF for different seasons as stored by monthly_means_from_plot_global.m
%graph=111;  %1D PDFs for different seasons as stored by monthly_means_from_plot_global.m
%graph=112 %gcm 1D pdfs
%graph=113; %gcm dirunal & monthly cycles
%graph=114 %timeseries from values from monthly_means... script - seasonal, monhtly, etc.
%graph=115; %calc gcm variables vs longitude
%graph=116; %plot gcm variables vs longitude
%graph=117; %gcm quick profiles at certain locations.
%graph=118; %Mean no. overpasses vs latitude
%graph=119; %Local solar time of the L3 min SZA - can perhaps tell if it wasn't the 10:30LT swath
%graph=120; %Solar and senosr orbit properties vs longitude
%graph=121; %SZA vs timeed
 % *** case 122 is the one used for the std CTT plots, etc. ***
%graph=122; %e.g. Nd vs SZA for all Re wavelengths. multiple means of 2D pdfs (lines on same plot) - e.g.
%different wavelengths for Re - also does single plots with error bars
%graph=122.1; %plots difference between high and low sza values at each sctt
%graph=123; %compare real L3 with mock L3
%graph=124; %fraction of cold clouds in model vs longitude
%graph=125; %CloudSat cold clouds precip vs total precip
%graph=126; %Terrain height with land fraction marked
%graph=127;  %Rain rate - ascending and descending
%graph=128; %Makes the table of relative differences in Nd from low to high
%SZA and the contibutions from tau and Re
%graph=129 %calculate the N lines (vs SZA) when hold tau constant and then Re constant
%graph=130; %PDFs of tau and Re at low and high SZA
%graph=131; %multi plot of means or PDFs from plotTime 2d pdfs
%graph=132; %pdfs from CALIPSO PDFs
%graph=133; %CF vs homog factor and other 1D means of 2D pdfs
%graph=134; %simple 1D PDFs
%graph=135; %Longitude plot showing the diurnal range with bars for lon bins
%graph=136; %plot of monthly means from monthly_means_from_plot_global - monthly, seasonal, etc.
%graph=137; %plots of max and min sza for a certain latitude
%graph=138; %PDFs of different parameters (e.g. 1.6, 2.1 and 3.7um
%differences between lowSZA and allSZA data)
%graph=139; %Plots of mean % differerence between <65 and allSZA plots vs latitude (1D means for the sides of the 2D plots)
%graph=140; % PDFs of various things from saved .mat files - e.g. Tau from GCMs (COSP) and MODIS at high CF
%graph=141; %Longitude vs Reff plots for POLDER vs MODIS colocated
%graph=142;  %Runs plotTimeheight several times to give multiple lines - PDFs or means
%graph=143;   %plot of various lines saved using save_vars_mat_run
%graph=144; %timeseries plot of F0, hbar, U, N from Fohn cross section
%graph=145; %sza vs local time
%graph=146; % vs longitude comparisons of WRF and aircraft (Fohn paper)
%graph=147;  %1d PDFs from saved xdat, ydat etc. to look at the pre/post mphys LWP issue - a PDF for each model
%graph=148;  %1d PDFs from saved xdat, ydat etc. to look at the pre/post mphys LWP issue - two different PDFS for each model
%graph=149; % Fraction of points for which >95% of LWP is removed by mphys
%graph=150; % Seasonal Nd cycles for multiple years
%graph=151; % Seasonal Nd cycles for 2007, different wavelengths
%graph=152; % UM timeseries
%graph=1522; %UM timseries from Moose NetCDF output
%graph=153; % Heterogeneity parameter comparison for SZA paper - Cahalan gamma_tau vs sigCTT
%graph=154; % Nd vs height for different months
%graph=155; % Aircraft transects at 3000m
%graph=156; % WRF timerseries of 3000m region of aircraft 
%graph=157; % Rothera and WRF timerseries at Rothera
%graph=158; %Plots of RWP contribution from rain for AMSRE
%graph=159; %Plots of precent diffs MODIS minus AMSRE LWP retrieval with precip
%graph=160; %Plots of precent contribution to MODIS LWP from RWP
%graph=161; %qR from the above simple model
%graph=162; %Some profiles of CF etc. from the UM cloud scheme runs

%% end of graph choice


watervap_defaults  %Default flags now done in here - generally should run this before a script that uses man_choose_water_graph


%    clear man_choose_water_graph %clear for next time
end

if ~exist('line_pattern')
    line_pattern(1).p=NaN;
    line_colour(1).c=NaN;
    marker_style(1).m=NaN;
end


if ~exist('savedir')
    savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';
end

if ~exist('i_highlight_path') | i_highlight_path==0
    time_highlight_path=[];
else
    clear i_highlight_path
end
if ~exist('iadd_terrain') | iadd_terrain==0
    iadd_terrain2=0;
else
    iadd_terrain2=1;
    clear iadd_terrain
end
if ~exist('iplot_latlon') | iplot_latlon==0
    iplot_latlon2=0;
else
    iplot_latlon2=1;
    clear iplot_latlon
end


if ~exist('highlight_type')
    highlight_type='fill';
end
    



f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
%  Ms=1/6*pi*rhoS*Ds.^3;
% 
% 
% a=find(times~=0);
% dt(2:length(times))=times(2:end)-times(1:end-1);
% dt(a(1))=dt(a(2));
% a=find(times==0);
% dt(a)=1; %avoi5d divide by zero
% dtt=repmat(dt,[150 1]);
% dtt2=repmat(dtt,[1 1 14]);
% dtt3=permute(dtt2,[1 3 2]);




if ~exist('subplotting'); subplotting=0; end
if ~exist('idir'); idir=1; end
if ~exist('file_type'); file_type=''; end
    
if ~exist('justplot')
    justplot=0;
end

if justplot==0

clear labs xdat ydat diff xpos ypos point_labs
figname='Vapour graph';



%xlab='Water Vapour Mixing Ratio (ppmv)';
%ylab='Height (km)';
     
     
    
%% start of graphs

switch graph       
     case 9999999
        % template to copy for new graph
        
        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';          
        titlenam = ['XXX for ' tstr];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Height (m)';
        xlab= 'WRF ice number concentration (L^{-1})';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        HGT=WRFUserARW(nca(1).nc,'Z',time,ilat(iloc),ilon(iloc));

        idat=idat+1;
        xdat(idat).x=0.005*exp(0.304*(-tc)); %
        ydat(idat).y=HGT; %
        labs(idat).l='';
        
        
    case 0  
        
        figname=titlenam;
        savename=figname;

       %
       start_time=0;
       start_time = datenum('01-Jan-1970')*24; %
       istyle=1;
       idat=0;
       clear end_time start_time2
       for idat2=1:length(xdat_import)
           idat=idat+1;
           
           xdat(idat).x=xdat_import(idat).x; %           
           ydat(idat).y=ydat_import(idat).y;
           labs(idat).l=labs_import(idat).l;
                  
%           line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
           
           end_time(idat) = max(xdat(idat).x);
           start_time2(idat) = min(xdat(idat).x);
           

           ismooth_x(idat)=ismooth_x_import(idat);
           ismooth_y(idat)=ismooth_y_import(idat);           
           
       end      

       
 
       if idate_ticks_fix==1
           date_ticks_fix
       end
     
     
   
              


        
    case 162
        % Profiles of CF, precip frac, etc. for the new precip fraction
        % coding using netCDF files generated from pp files. Did this via a
        % transfer to MASS, but this stage is not necessary since the LWP,
        % RWP, etc. are in the pp files. But will need to copy things to
        % MASS for the bigger runs to avoid filling the ASCI disk.
        
        iload_nc = 1; %flag for whether to laod the netCDF files or not
        
        UM_run = 'xkjk';
        nc_filename = 'xkjki_LWP_RWP_CFs_12.pp.nc';
        it_nc = 1; %time index to use
        
        nc_dir = ['/home/disk/eos1/d.grosvenor/UM/' UM_run '/output/'];
        nc_filepath = [nc_dir nc_filename];
        

        titlenam = ['Various cloud and precip profiles'];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 1];
        
        izlim=1;
        zmin=0;
        zmax=5;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Height (km)';
        xlab= 'Fraction';



        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        if iload_nc==1
        
        %---- Read data if not present
          nc_UM = netcdf([nc_filepath],'nowrite');
          %Data is labbellbed as 'unspecified, 'unspecified_1'  etc. after conversion to
          %NetCDF. These correspond to the same order as in the pp file.
          
          z_um = nc_UM{'hybrid_ht'}(:);
          
          LWP=nc_UM{'unspecified'}(it_nc,:,:,:);
          RWP=nc_UM{'unspecified_1'}(it_nc,:,:,:);
          
          %find position of max RWP to get a location with some rain
          [max_rwp,imax_rwp] = maxALL(RWP);

          %Specify which location to use
          inds = imax_rwp;
          
%           idat=idat+1;
%           labs(idat).l = 'Cloud fraction post mphys';
%           xdat(idat).x = nc_UM{'unspecified_2'}(it_nc,:,inds(1),inds(2));
%           
%           idat=idat+1;
%           labs(idat).l = 'Cloud fraction pre mphys';
%           xdat(idat).x = nc_UM{'unspecified_5'}(it_nc,:,inds(1),inds(2));
% 
%           idat=idat+1;
%           labs(idat).l = 'Precip fraction post mphys';
%           xdat(idat).x = nc_UM{'unspecified_3'}(it_nc,:,inds(1),inds(2));
% 
%           idat=idat+1;
%           labs(idat).l = 'Precip fraction pre mphys';
%           xdat(idat).x = nc_UM{'unspecified_6'}(it_nc,:,inds(1),inds(2));
% 
          idat=idat+1;
          labs(idat).l = 'Precip frac in-cloud post mphys';
          xdat(idat).x = nc_UM{'unspecified_4'}(it_nc,:,inds(1),inds(2));

          idat=idat+1;
          labs(idat).l = 'Precip frac in-cloud pre mphys';
          xdat(idat).x = nc_UM{'unspecified_7'}(it_nc,:,inds(1),inds(2));
          
        end
        
        for idat=1:length(xdat)
            ydat(idat).y = z_um/1e3;
        end
          

          
          
        

        
        
        
    case 161
        % qR contribution for given LWP and Nd
        

        titlenam = ['qR for given cloud LWP and N_d'];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 180];
        
        izlim=0;
        zmin=0;
        zmax=20;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='qR (g m^{-3})';
        xlab= 'LWP_c + RWP (g m^{-2})';



        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        %---- run scrip to get results
        RWP_from_precip_run2
        
        labs = labs_RWP;

        for idat=1:3
            qR=PR_save{idat}/3600./Vfall_save{idat}*1e3;
            ydat(idat).y = qR;
            RWP = RWP_save{idat};
            xdat(idat).x=LWP+RWP;
        end


        
        
    case 160
        % RWP contribution for given LWP and Nd for MODIS retrieval
        
        %---- run scrip to get results
        RWP_from_precip_run
%        RWP_from_precip_run2
        

        titlenam = ['Percentage of RWP detected by MODIS for r_{vR} = ' num2str(rv*1e6) ' \mum'];
%        titlenam = ['Percentage of RWP detected by MODIS'];        
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 180];
        
        izlim=0;
        zmin=-40;
        zmax=0;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Percentage detected (%)';
        xlab= 'LWP_c + RWP (g m^{-2})';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        

        
        labs = labs_RWP;


        
        idat=0;
        idat=idat+1;
        RWP = RWP_01_N50;
        ydat(idat).y = 100 * (LWP_R_modis_01_N50) ./ (RWP);
        xdat(idat).x=LWP+RWP;
        
%         idat=idat+1;
%         RWP = RWP_02_N50;
%         ydat(idat).y = 100 * (LWP_R_modis_02_N50) ./ (RWP);
%         xdat(idat).x=LWP+RWP;
        
        idat=idat+1;
        RWP = RWP_01_N100;
        ydat(idat).y =  100 * (LWP_R_modis_01_N100) ./ (RWP);
        xdat(idat).x=LWP+RWP;
        
%         idat=idat+1;
%         RWP = RWP_02_N100;
%         ydat(idat).y = 100 * (LWP_R_modis_02_N100) ./ (RWP); 
%         xdat(idat).x=LWP+RWP;
        
    case 159
        % RWP contribution for given LWP and Nd for MODIS retrieval
        
        %---- run scrip to get results
%        RWP_from_precip_run
        RWP_from_precip_run2
        
%        titlenam = ['Percentage bias of MODIS relative to AMSR-E due to lack of rain detection by MODIS for r_{vR} = ' num2str(rv*1e6) ' \mum'];
        titlenam = ['Percentage bias of MODIS relative to AMSR-E due to lack of rain detection by MODIS'];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 180];
        
        izlim=1;
        zmin=-20;
        zmax=0;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='PB_{MODIS} (%)';
        xlab= 'LWP_c + RWP (g m^{-2})';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        

        
%         labs = labs_RWP;
% 
%         for idat=1:2
%             xdat(idat).x=LWP; %
%         end
%         
%         idat=0;
%         idat=idat+1;
%         RWP = RWP_01_N50;
%         ydat(idat).y = PB_01_N50;
%         xdat(idat).x=LWP+RWP;
%         
% %         idat=idat+1;
% %         RWP = RWP_02_N50;        
% %         ydat(idat).y = PB_02_N50; 
% %         xdat(idat).x=LWP+RWP;
%         
%         idat=idat+1;
%         RWP = RWP_01_N100;
%         ydat(idat).y = PB_01_N100; 
%         xdat(idat).x=LWP+RWP;
% %         
% %         idat=idat+1;
% %         RWP = RWP_02_N100;
% %         ydat(idat).y = PB_02_N100;   
% %         xdat(idat).x=LWP+RWP;


   labs = labs_RWP;

        for idat=1:3
            RWP = RWP_save{idat};
            ydat(idat).y = PB_save{idat};
            xdat(idat).x=LWP+RWP;
        end
        
    case 158
        % RWP contribution for given LWP and Nd
        

        titlenam = ['Contribution of RWP to total LWP for given cloud LWP and N_d'];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 180];
        
        izlim=0;
        zmin=0;
        zmax=20;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Contribution to total LWP from RWP (%)';
        xlab= 'LWP_c + RWP (g m^{-2})';



        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        %---- run scrip to get results
        RWP_from_precip_run2
        
        labs = labs_RWP;

        for idat=1:3
            RWP = RWP_save{idat};
            ydat(idat).y = 100*RWP./(LWP+RWP);
            xdat(idat).x=LWP+RWP;
        end

        
        
    case 157
        % Rothera timeseries comparisons
        % run read_AWS_data.m first to laod AWS data  - need to find the
        % file first!!
        % and read_rothera_dat to laod the Rothera data

        
        irecalc=0; %Set to zero if don't want to do the lengthy calculation again and just want to replot
        %irecalc='ALL' % = recalc all
%        irecalc='Field only' % = already know ilat ilon        
           %N.B. have saved ilat and ilon values for the flight path
           %specified below in a .mat file - so should only need to do
           %'Field only'
           
           %select variable
        ylab='Mean Sea Level Pressure (hPa)';  ivar_roth = 5; ivar_AWS=6; offset=315;%ivar_roth is the index in the 
            % rothera data file for the variable required.
%        ylab='Wind direction (degrees)';   ivar_roth = 8; ivar_AWS=8; offset=0;

        lat_roth = -67.5667; lon_roth = -68.1333; %Rothera
        lat_aws = -67.01; lon_aws = -61.55; %AWS          
        %Time for WRF
%        time=16; %6th Jan:- 11=6UTC, 13=12UTC, 15=18UTC, 16=21 UTC, 17=00UTC
           %Aircraft took off at 19:20 and started descent at 20:15 (flew
           %at 3000m). So, midpoint of the upper level run was around
           %19:45. Closest to 21 UTC model output.
           
        iaxis_square=0;
                   
        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';          
        titlenam = ['Rothera and WRF timeseries'];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.


        xlab= 'Time (UTC)';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        
% Get WRF lat lon indices
      [ilat_roth,ilon_roth] = getind_latlon_quick(lat2d.var,lon2d.var,lat_roth,lon_roth,0.1);
      [ilat_aws,ilon_aws] = getind_latlon_quick(lat2d.var,lon2d.var,lat_aws,lon_aws,0.1);      
      time_utc_wrf = [0:3:72];
      
        idat=idat+1;
        xdat(idat).x=datenum('05-Jan-2006 00:00')+time_utc_wrf/24; %convert to Matlab time
%        ydat(idat).y = p3000_time;
        labs(idat).l='WRF Rothera';
        switch ylab
            case 'Mean Sea Level Pressure (hPa)'
                ydat(idat).y = 0.01* nc{'PSFC'}([1:25],ilat_roth,ilon_roth);
            case 'Wind direction (degrees)'
                u=WRFUserARW(nc,'U10',[1:25]);
                v=WRFUserARW(nc,'V10',[1:25]);
                u = u.var(:,ilat_roth,ilon_roth);
                v = v.var(:,ilat_roth,ilon_roth);
                ydat(idat).y = wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat(iloc),ilon(iloc),DX,DY);                                
                
        end
        
        idat=idat+1;
        xdat(idat).x=datenum('05-Jan-2006 00:00')+time_utc_wrf/24; %convert to Matlab time       
%        ydat(idat).y = p3000_time;
        labs(idat).l='WRF AWS';
         switch ylab
            case 'Mean Sea Level Pressure (hPa)'
                 ydat(idat).y = 0.01* nc{'PSFC'}([1:25],ilat_aws,ilon_aws);
             case 'Wind direction (degrees)'
                u=WRFUserARW(nc,'U10',[1:25]);
                v=WRFUserARW(nc,'V10',[1:25]);
                u = u.var(:,ilat_aws,ilon_aws);
                v = v.var(:,ilat_aws,ilon_aws);
                ydat(idat).y = wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat(iloc),ilon(iloc),DX,DY);                
         end
         
        
        idat=idat+1;
        xdat(idat).x=datenum('05-Jan-2006 00:00')+time_utc_wrf/24; %convert to Matlab time
%        ydat(idat).y = 0.01* nc{'PSFC'}([1:25],ilat_aws,ilon_aws);
        ydat(idat).y = meanNoNan(p3000_time,1) + offset;
        labs(idat).l=['WRF upper plus ' num2str(offset)];
        
        idat=idat+1;
        xdat(idat).x=time_roth_matlab; %Matlab time
        ydat(idat).y=dat_rothera(ivar_roth,:); %take average over all lons     
        labs(idat).l='Station';        
        
%        idat=idat+1;
%        xdat(idat).x=AWShrs(inds)'/24; %Matlab time
%        ydat(idat).y=AWSdat2(ivar_AWS,inds)'; % 
%        labs(idat).l='AWS';        

switch ylab
    case 'Wind direction (degrees)'
         for idat=1:length(ydat)
             i90 = find(ydat(idat).y<90);
             ydat(idat).y(i90) = ydat(idat).y(i90)+360;
         end
end
                
          [Y2,M2,D2,H2,MM2] = datevec(xdat(idat).x);
          inew_day2 = find(H2==0);
          
          tick_freq = 6; %no. hours between ticks
          hr_start = floor(H2(1)/tick_freq)*tick_freq;

       xtickvals_set = [datenum(Y2(1),M2(1),D2(1),hr_start,MM2(1),0) : 6/24 : xdat(1).x(end)];
     
       [Y,M,D,H,MM] = datevec(xtickvals_set);
       month = datestr(xtickvals_set,'mmm');
       time_str144 = datestr(xtickvals_set,'HH');
        
       for ixtick_lab=1:length(xtickvals_set)
           xticklabs_set{ixtick_lab} = [time_str144(ixtick_lab,:)];
       end
       
       inew_day = find(H==0);
        
       for ixtick_lab=inew_day
           xticklabs_set{ixtick_lab} = [num2str(D(ixtick_lab)) '-' month(ixtick_lab,:)];
       end 
        
    iset_xticks = 1;
    iset_xticklabs=1;
    
    
        
        
    case 156
        % Upper level comparison to aircraft at 3000m (in response to ACPD
        % ref). Timeseries of means over eastern portion of teh leg (away
        % from mountain)
        
        %Change settings in along_flight_track_timeseries.m
        
        irecalc=0; %Set to zero if don't want to do the lengthy calculation again and just want to replot
        %irecalc='ALL' % = recalc all
%        irecalc='Field only' % = already know ilat ilon        
           %N.B. have saved ilat and ilon values for the flight path
           %specified below in a .mat file - so should only need to do
           %'Field only'
           
           %select variable
        ylab='Pressure (hPa)';           ivar_air = 4;
%        ylab='Wind direction (degrees)';   ivar_air = 10;                   
           
        %Time for WRF
%        time=16; %6th Jan:- 11=6UTC, 13=12UTC, 15=18UTC, 16=21 UTC, 17=00UTC
           %Aircraft took off at 19:20 and started descent at 20:15 (flew
           %at 3000m). So, midpoint of the upper level run was around
           %19:45. Closest to 21 UTC model output.
        
           iaxis_square=0;
        
        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';          
        titlenam = ['Upper level aircraft track WRF timeseries'];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.


        xlab= 'Time (UTC)';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
 

% ----------------------------------------
switch irecalc
    case 'ALL'
        along_flight_track_timeseries
end
% ----------------------------------------




        idat=idat+1;
        xdat(idat).x=datenum('05-Jan-2006 00:00')+time_utc_p3000/24; %convert to Matlab time
        ydat(idat).y=meanNoNan(p3000_time,1); %take average over all lons
        labs(idat).l='';
        
               
          [Y2,M2,D2,H2,MM2] = datevec(xdat(idat).x);
          inew_day2 = find(H2==0);
          
          tick_freq = 6; %no. hours between ticks
          hr_start = floor(H2(1)/tick_freq)*tick_freq;

       xtickvals_set = [datenum(Y2(1),M2(1),D2(1),hr_start,MM2(1),0) : 6/24 : xdat(1).x(end)];
     
       [Y,M,D,H,MM] = datevec(xtickvals_set);
       month = datestr(xtickvals_set,'mmm');
       time_str144 = datestr(xtickvals_set,'HH');
        
       for ixtick_lab=1:length(xtickvals_set)
           xticklabs_set{ixtick_lab} = [time_str144(ixtick_lab,:)];
       end
       
       inew_day = find(H==0);
        
       for ixtick_lab=inew_day
           xticklabs_set{ixtick_lab} = [num2str(D(ixtick_lab)) '-' month(ixtick_lab,:)];
       end 
        
    iset_xticks = 1;
    iset_xticklabs=1;
    
    
    %calculate the aircraft value for in-bound bit
    i3000_IN = i3000_2(1):i3000_2(end);
    mean_aircraft_val_IN = meanNoNan(dat_flt19_proc(i3000_IN,ivar_air),1);
    mean_aircraft_time_IN = meanNoNan(time_flt19(i3000_IN,1),1);
    mean_aircraft_time_IN_matlab = datenum('06-Jan-2006 00:00') + mean_aircraft_time_IN/24;
    
    %Ascent out of region
%    i3000_OUT = 717471:797475;
    i3000_OUT = 717471:797170;  
       %Data cuts off at -63.7 lon
    mean_aircraft_val_OUT = meanNoNan(dat_flt19_proc(i3000_OUT,ivar_air),1);
    mean_aircraft_time_OUT = meanNoNan(time_flt19(i3000_OUT,1),1);
    mean_aircraft_time_OUT_matlab = datenum('06-Jan-2006 00:00') + mean_aircraft_time_OUT/24;
    
    %calculate the aircraft value for in-bound bit, but limited to same lon
    %as outbound - not necessary for case C
    i3000_IN2= 184101:238242;
    mean_aircraft_val_IN2 = meanNoNan(dat_flt19(i3000_IN2,ivar_air),1);
    mean_aircraft_time_IN2 = meanNoNan(time_flt19(i3000_IN2,1),1);
    mean_aircraft_time_IN2_matlab = datenum('06-Jan-2006 00:00') + mean_aircraft_time_IN2/24;
    
    
%    plot(mean_aircraft_time_IN_matlab,mean_aircraft_val_IN,'bo','markerfacecolor','b','markersize',10)    
%    plot(mean_aircraft_time_OUT_matlab,mean_aircraft_val_OUT,'bd','markerfacecolor','b','markersize',10)
    
    
        
    case 155
        % Upper level comparison to aircrfat at 3000m (in response to ACPD
        % ref).
        
        irecalc=0; %Set to zero if don't want to do the lengthy calculation again and just want to replot
        %irecalc='ALL' % = recalc all
%        irecalc='Field only' % = already know ilat ilon        
           %N.B. have saved ilat and ilon values for the flight path
           %specified below in a .mat file - so should only need to do
           %'Field only'
        irecalc = 'Use p3000_time';    %if have calculated in along_flight_track

        ylab='Pressure (hPa)'; ivar_air=4;
        ylab='Wind direction (degrees)'; ivar_air=10;

           
        %Time for WRF
        time=16; %6th Jan:- 11=6UTC, 13=12UTC, 15=18UTC, 16=21 UTC, 17=00UTC
           %Aircraft took off at 19:20 and started descent at 20:15 (flew
           %at 3000m). So, midpoint of the upper level run was around
           %19:45. Closest to 21 UTC model output.
           
        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';          
        titlenam = ['Upper level aircraft comparison for ' tstr];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.



        xlab= 'Longitude';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        switch irecalc
            case {'ALL','Field only'}

                %indices for the first 3000m leg
                i3000 = [62921:238242];
                %make some coarse indices for faster calcs
                i3000_2=[i3000(1):100:i3000(end)];

                LATS = dat_flt19(i3000_2,2);
                LONS = dat_flt19(i3000_2,3);

                %Load in previously calculated ilat and ilon indices for the WRF
                %array that correspond to the flight track (LATS and LONS above)
                load_file_Ant='~/mat_files_various/Antarctica/ilat_ilon_3000m.mat';
                load(load_file_Ant,'ilat','ilon');

        end
        
        field_str = 'Pressure'; 
        zfind = 3000; %(m)
        
        switch irecalc
            case 'ALL'
                [field,HGT,iz3000,ilat,ilon] = get_field_along_flight_path(nc,lat2d,lon2d,time,LATS,LONS,zfind,field_str);
            case 'Field only'
                [field,HGT,iz3000,ilat,ilon] = get_field_along_flight_path(nc,lat2d,lon2d,time,LATS,LONS,zfind,field_str,ilat,ilon);
            case 'Use p3000_time';    %if have calculated in along_flight_track  
                iselect = find(times==time);
                field = p3000_time(:,iselect);
        end



        idat=idat+1;
        xdat(idat).x=dat_flt19_proc(i3000,3); %
        ydat(idat).y=dat_flt19_proc(i3000,ivar_air); %
        labs(idat).l='Aircraft';
        
        idat=idat+1;
        xdat(idat).x=dat_flt19(i3000_2,3); %
        ydat(idat).y=field; %
        labs(idat).l='WRF';
        

        
    case 154

        
         %switch between different regions etc. using titlenam:-

        titlenam = ['Nd vs CTH from MODIS individual mockL3 for 60-70S, 160-60W'];
        
        switch titlenam
            case 'Nd vs CTH from MODIS individual mockL3 for 60-70S, 160-60W'
                filename = ['saved_vars_for_Nd_vs_height_60-70S_60_160W'];  %CF80 all CTH, all CTT
                filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
                filename_savevars = [filedir_savevars filename '.mat'];
                
                tags={'Jan','March','Sep'};
        end

        
        load(filename_savevars);

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=[0 100];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
%        nmark=0;

ierror_bars='horiz2';   %whether to plot error bars

        xlab= 'N_d (cm^{-3})';
        ylab= 'Cloud Top Height (km)';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;      
        

        
        %see screen_edits_driver for the screening thresholds
        % and screening_eval_str_driver
        

        logflag=0;
        ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
        
        
%        Nsamp_all = eval(['season_Ndatap_driver_' tag '{1};']); %Ndatap is actually the number of gridboxes with data
        %So it shows the amount of area covered (with at least one time
        %sample for each).
%        max_Nsamp = max(Nsamp_all);
%        thresh_min_samp = max_Nsamp*0.5; %require x% of all the area covered     
        

       
       istyle=1;
       idat=0;
       for idat2=1:length(tags)
           tag = tags{idat2};
           idat=idat+1;
           
           xdat(idat).x=eval(['xdat_' tag '.x;']);
           ydat(idat).y=eval(['ydat_' tag '.y;']);
           labs(idat).l=tag;
           
           errordatL(idat).dat = eval(['errordatL_' tag '.dat;']);
           errordatU(idat).dat = eval(['errordatU_' tag '.dat;']);
           
           max_err = 8;
           icut = find(errordatL(idat).dat>max_err/2);
           ydat(idat).y(icut)=NaN;
           

           
       end                     

%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
       
%       save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],['BASE screening:- ' eval(['screening_eval_str_' tag])] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );
%       write_notes_to_file( save_notes_filepath, ['AND :- ' eval(['screen_edits_' tag '{4}']) ] );
%       write_notes_to_file( save_notes_filepath, ['FOR REGION :- 40-60S, 0-360E'] );
%        write_notes_to_file( save_notes_filepath, ['with mockL3 screenings:-' eval(['thresh_str_multiL2_driver_' tag '{1}'] ); %need to sort this one out

       

        
     case 153
        filename = ['Het_params_comparison_SZA_paper_Cahaan_vs_sigCTT']; %Mean sig_CTT for each
             %gamma_tau bin.
        filename = ['Het_params_comparison_SZA_paper_Cahaan_vs_sigCTT_binned_by_sigCTT'];
             %Mean gamma_tau for each sig_CTT bin.
             
        ierror_bars='horiz2';   %whether to plot error bars
        
        irestrict_Ndatap=0; %flag to say whether to restrict using the no. datapoints 
        thresh_Ndatap=20; %Minimum number of datapoints requried
        

        filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
        filename_savevars = [filedir_savevars filename '.mat'];
        load(filename_savevars);
        
        tags = {'het_low_SZA','het_high_SZA'};
        tags_disp = {'50\leqSZA<55^\circ','SZA\geq75^\circ'};
        
        thresh_str_display = [thresh_str '_' tags{1}];
       
        titlenam = ['Relationship between different heterogeneity parameters ' thresh_str_display];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 0.25];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=eval(['ylab_' tags{1}]);
        xlab=eval(['xlab_' tags{1}]);
        
%         switch filename
%             case 'saved_vars_for_LWP_vs_CF' %Nd
%                 xlab= 'MOD06 liquid cloud fraction';
%         end


        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        

        for idat=1:length(tags)

            %        idat=idat+1;
            ydat(idat).y=eval(['ydat_' tags{idat} '.y']); %
            xdat(idat).x=eval(['xdat_' tags{idat} '.x']); %
            errordatU(idat).dat = eval(['errordatU_' tags{idat} '.dat;']);
            errordatL(idat).dat = eval(['errordatL_' tags{idat} '.dat;']);
            labs(idat).l=tags_disp{idat};
            
            %Remove points based on no. datapoints
            if irestrict_Ndatap==1
                inan=find(eval(['Ndatap_' tags{idat} '<thresh_Ndatap;']));
                ydat(idat).y(inan)=NaN;
            end

        end
        
        
        
        istyle=1; %the counter for the styles we specify

        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='^'; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='^'; istyle=istyle+1;


%           titlenam = [titlenam ' low sensSZA'];
        
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
        

    
        
       save_notes_filepath = write_notes_to_file([savenotes_filedir filename],['Het param comparison for low and high SZA.']);
       write_notes_to_file( save_notes_filepath, ['Low SZA screening:- ' eval(['thresh_str_' tags{1}]) ]);
       write_notes_to_file( save_notes_filepath, ['High SZA screening:- ' eval(['thresh_str_' tags{2}]) ]);       
       write_notes_to_file( save_notes_filepath, ['Title from PlotTime (contains corr coeff) for low SZA :- ' eval(['titlenam_' tags{1}]) ]);
       write_notes_to_file( save_notes_filepath, ['Title from PlotTime (contains corr coeff) for high SZA :- ' eval(['titlenam_' tags{2}]) ]);


    case 152   
        %Use read_UM_timeser in work/Leeds to read in the text files
        
        mfac=1;
        
        tag='LWP';
        tag='LowCF';
        tag='LWP_Adrian';
        
        switch tag
            case {'LWP','LWP_Adrian'}
                titlenam = ['Domain mean LWP timeseries'];
                ylab='LWP (g m^{-2})';
                mfac=1e3;
            case 'LowCF'
                titlenam = ['Domain mean low cloud fraction timeseries'];
                ylab='Cloud Fraction';
        end
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=[0 100];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        nmark=0;

        xlab= 'UTC Time (hours)';



        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;   
        
%        idatetick=1;      datetick_type = 'HH:MM';
        
        %data is ordered for 2006-2010 for CF 0.1-0.4 (cells 1-5), then 2006-2010 for
        %0.4-0.6 (cells 6-10), etc. (then 0.6-0.8 [cells 11-15], then 0.8-1.0 [cells 16-20])
        
        %see screen_edits_driver for the screening thresholds
        % and screening_eval_str_driver
        

        logflag=0;
        ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
        
        clear labs_UM
     
       labs_UM{1} = '4a';
       labs_UM{2} = '3d';
       labs_UM{3} = '3d no cloud scheme';
       labs_UM{4} = '4a 142L';  
       labs_UM{5} = '4a Smag BL'; 
       
       clear labs_UM       
       idat=1;
       labs_UM{idat} = '(s) No cloud scheme, Rain OFF'; idat=idat+1; %s
       labs_UM{idat} = '(d) Rain OFF, RHcrit=0.999999999, minCF=0.01, Sat adjust'; idat=idat+1; %d
       labs_UM{idat} = '(w) Rain OFF, RHcrit=0.999999999, minCF=0.05, Sat adjust'; idat=idat+1; %w
       
       labs_UM{idat} = '(g) Rain ON, seeding, RHcrit=0.999999999, minCF=0.05, Sat adjust OFF'; idat=idat+1; %g
       labs_UM{idat} = '(h) Rain ON, seeding OFF, RHcrit=0.999999999, minCF=0.05, Sat adjust OFF'; idat=idat+1; %h  
       
       labs_UM{idat} = '(i) Rain ON, seeding,. RHcrit=UKV, minCF=0.05, Sat adjust OFF'; idat=idat+1; %i
       labs_UM{idat} = '(k) Rain ON, seeding OFF,. RHcrit=UKV, minCF=0.05, Sat adjust OFF'; idat=idat+1; %k  
       
       
       clear labs_UM
       idat=1;
       labs_UM{idat} = '(xkjks) CASIM, No cloud scheme'; idat=idat+1; %s
       labs_UM{idat} = '(xkjkg) CASIM, Cloud Scheme, RHcrit =0.999'; idat=idat+1; %d
       labs_UM{idat} = '(xkcvo) 3D mphys, Cloud Scheme, UKV RHcrit'; idat=idat+1; %w              
       labs_UM{idat} = '(xkjki) CASIM Cloud scheme, UKV RHcrit'; idat=idat+1; %i
       labs_UM{idat} = '(xkcvv) 3D mphys, Cloud Scheme, RHcrit 0.999'; idat=idat+1; %k      
       
       
       clear labs_UM
       idat=1;       %all based on CASIM Cloud scheme, UKV RHcrit        
       labs_UM{idat} = '(xkjki) 100 cm^{-3}'; idat=idat+1; %i
       labs_UM{idat} = '(xkjkx) 50 cm^{-3}'; idat=idat+1; %k      
       labs_UM{idat} = '(xkjky) 10 cm^{-3}'; idat=idat+1; %s     

       %1000km 26th Oct POC cases
       clear labs_UM
       idat=1;       %all based on CASIM Cloud scheme, UKV RHcrit        
       labs_UM{idat} = '(xkqkf) No cloud scheme 100cm^{-3}'; idat=idat+1; %i
       labs_UM{idat} = '(xkqkg) Cloud scheme 100cm^{-3}'; idat=idat+1; %k      
       labs_UM{idat} = '(xkqkj) Cloud scheme 400cm^{-3}'; idat=idat+1; %k           

       
       % -- 
       %
       start_time=0;
       start_time = datenum('01-Jan-1970')*24; %
       istyle=1;
       idat=0;
       clear end_time
       for idat2=1:length(labs_UM)
           idat=idat+1;
           
%           xdat(idat).x=time_UM{idat2}(:) - time_UM{idat2}(1); %
           xdat(idat).x=(time_UM{idat2}(:) + start_time)/24; %           
           ydat(idat).y=eval(['mfac*' tag '_UM{idat2}(:);']);
           labs(idat).l=labs_UM{idat2};
           
        
           line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
           
           end_time(idat) = xdat(idat).x(end);
       end      

       
       [temp,idat_max] = max(end_time);
       [Y,M,D,H,MM] = datevec(xdat(idat_max).x);
       H1 = H(1);
       H=H1:H1+6;
       istart_temp=find(rem(H,6)==0);
       start_val=H(istart_temp(1));
       Htime = datenum(Y(1),M(1),D(1),start_val,0,0);        
              
       xtickvals_set = [Htime : 6/24 : xdat(idat_max).x(end)];
        
        [Y,M,D,H,MM] = datevec(xtickvals_set);
        month = datestr(xtickvals_set,'mmm');
        time_str144 = datestr(xtickvals_set,'HH');
        
        
        
        for ixtick_lab=1:length(xtickvals_set)
            xticklabs_set{ixtick_lab} = [time_str144(ixtick_lab,:)];
        end
        
       istart_temp=find(H==0);
        for ixtick_lab=[istart_temp:4:length(xtickvals_set)]
            xticklabs_set{ixtick_lab} = [num2str(D(ixtick_lab)) '-' month(ixtick_lab,:)];
        end 
        
     iset_xticks = 1;
     iset_xticklabs=1;
     
     
     
       
%        for idat2=2:2
%            idat=idat+1;
%            
%            xdat(idat).x=[1:12]; %
%            ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
%            labs(idat).l='2.1 \mum';
%            
%            err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
%            %Sampling error is not really useful if only have a few points
%            %since the std dev is very small!
%            N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
%            icut = find(N_samp<thresh_min_samp);
%            ydat(idat).y(icut)=NaN;
%            
%            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.6 0.6 0.6]; marker_style(istyle).m='d'; line_widths(istyle).l = 2; istyle=istyle+1;
%            
% 
%        end
%        
%        for idat2=3:3
%            idat=idat+1;
%            
%            xdat(idat).x=[1:12]; %
%            ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
%            labs(idat).l='3.7 \mum';
%            
%            err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
%            %Sampling error is not really useful if only have a few points
%            %since the std dev is very small!
%            N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
%            icut = find(N_samp<thresh_min_samp);
%            ydat(idat).y(icut)=NaN;
%            
%            line_pattern(istyle).p= '--';  line_colour(istyle).c=[1.0 0.0 0.0]; marker_style(istyle).m='s'; line_widths(istyle).l = 2; istyle=istyle+1;
%            
% 
%        end

       
%       save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],['BASE screening:- ' eval(['screening_eval_str_' tag])] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );
%       write_notes_to_file( save_notes_filepath, ['AND :- ' eval(['screen_edits_' tag '{4}']) ] );
%       write_notes_to_file( save_notes_filepath, ['FOR REGION :- 40-60S, 0-360E'] );
       
     
              
% --- Now plot the means of all the years   

% %CF 0.1-0.4
%        %Gather the ydat data into an array, so can do MeanNoNan
%       dat_CF01 = NaN*ones([12 5]);
%       idat2=0;
%       for iyears_CF01=1:5
%           idat2=idat2+1;
%           dat_CF01(:,idat2) = ydat(iyears_CF01).y;
%       end
%       
%       %Put into plot structure
%       idat=idat+1;
%       xdat(idat).x=[1:12]; %
%       ydat(idat).y=meanNoNan(dat_CF01,2);
%       labs(idat).l='CF 0.1-0.4';
%       line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 3; istyle=istyle+1;
% 
%            
% %CF 0.8-1.0
%        %Gather the ydat data into an array, so can do MeanNoNan
%       dat_CF01 = NaN*ones([12 5]);
%       idat2=0;
%       for iyears_CF01=6:10
%           idat2=idat2+1;
%           dat_CF01(:,idat2) = ydat(iyears_CF01).y;
%       end
%       
%       %Put into plot structure
%       idat=idat+1;
%       xdat(idat).x=[1:12]; %
%       ydat(idat).y=meanNoNan(dat_CF01,2);
%       labs(idat).l='CF 0.8-1.0';
%       line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.6 0.6 0.6]; marker_style(istyle).m='o'; line_widths(istyle).l = 3; istyle=istyle+1;
        
        
    
    case 151
        % Seasonal Nd cycles for 2007, CF>0.8, but all wavelengths
         %switch between different regions etc. using titlenam:-
        titlenam = ['Seasonal cycle of Nd from MODIS mockL3 for 40-60S, 0-360E'];
        titlenam = ['Seasonal cycle of Nd from MODIS individual mockL3 for 60-70S, 160-60W'];
        
        switch titlenam
            case 'Seasonal cycle of Nd from MODIS mockL3 for 40-60S, 0-360E'
                %        filename = ['saved_vars_for_Nd_2007_all_wavelengths_CF80_seasonal_SouthernOcean_mockL3'];  %for the 'MODIS Nd multiple year seasonal cycles'
                filename = ['saved_vars_for_Nd_2007_all_wavelengths_CF99_CTT_268_seasonal_SouthernOcean_mockL3'];  %As above, but for CF>0.99 and CTT>268K
                filedir_savevars = '/home/disk/eos1/d.grosvenor/';
                filename_savevars = [filedir_savevars filename '.mat'];
                tag = 'Nd_2007_wavelengths';

            case 'Seasonal cycle of Nd from MODIS individual mockL3 for 60-70S, 160-60W'           
                filename = ['saved_vars_Nd2007_all_wlengths_CF80_seasonal_60-70S_60-160W_L3_individ_20140425T063158'];  %CF80 all CTH, all CTT
                filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
                filename_savevars = [filedir_savevars filename '.mat'];
                tag = 'Nd_2007_wavelengths';
        end

        

        load(filename_savevars);

        

        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=[0 100];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
%        nmark=0;

        ylab='N_d (cm^{-3})';
        xlab= 'Month';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;      
        

        
        %see screen_edits_driver for the screening thresholds
        % and screening_eval_str_driver
        

        logflag=0;
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
        
        
        Nsamp_all = eval(['season_Ndatap_driver_' tag '{1};']); %Ndatap is actually the number of gridboxes with data
        %So it shows the amount of area covered (with at least one time
        %sample for each).
        max_Nsamp = max(Nsamp_all);
        thresh_min_samp = max_Nsamp*0.5; %require x% of all the area covered     
        

       % -- Data is ordered as 1.6, 2.1 and then 3.7um in cells of
       % season_Ndatap_driver_...  All for CF>0.8
       istyle=1;
       idat=0;
       for idat2=1:1
           idat=idat+1;
           
           xdat(idat).x=[1:12]; %
           ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
           labs(idat).l='1.6 \mum';
           
           err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
           %Sampling error is not really useful if only have a few points
           %since the std dev is very small!
           N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
           icut = find(N_samp<thresh_min_samp);
           ydat(idat).y(icut)=NaN;
           
           line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
           
       end                     
       
       for idat2=2:2
           idat=idat+1;
           
           xdat(idat).x=[1:12]; %
           ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
           labs(idat).l='2.1 \mum';
           
           err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
           %Sampling error is not really useful if only have a few points
           %since the std dev is very small!
           N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
           icut = find(N_samp<thresh_min_samp);
           ydat(idat).y(icut)=NaN;
           
           line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.6 0.6 0.6]; marker_style(istyle).m='d'; line_widths(istyle).l = 2; istyle=istyle+1;
           

       end
       
       for idat2=3:3
           idat=idat+1;
           
           xdat(idat).x=[1:12]; %
           ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
           labs(idat).l='3.7 \mum';
           
           err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
           %Sampling error is not really useful if only have a few points
           %since the std dev is very small!
           N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
           icut = find(N_samp<thresh_min_samp);
           ydat(idat).y(icut)=NaN;
           
           line_pattern(istyle).p= '--';  line_colour(istyle).c=[1.0 0.0 0.0]; marker_style(istyle).m='s'; line_widths(istyle).l = 2; istyle=istyle+1;
           

       end

       
       save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],['BASE screening:- ' eval(['screening_eval_str_' tag])] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );
%       write_notes_to_file( save_notes_filepath, ['AND :- ' eval(['screen_edits_' tag '{4}']) ] );
%       write_notes_to_file( save_notes_filepath, ['FOR REGION :- 40-60S, 0-360E'] );
%        write_notes_to_file( save_notes_filepath, ['with mockL3 screenings:-' eval(['thresh_str_multiL2_driver_' tag '{1}'] ); %need to sort this one out

       
     
              
% --- Now plot the means of all the years   

% %CF 0.1-0.4
%        %Gather the ydat data into an array, so can do MeanNoNan
%       dat_CF01 = NaN*ones([12 5]);
%       idat2=0;
%       for iyears_CF01=1:5
%           idat2=idat2+1;
%           dat_CF01(:,idat2) = ydat(iyears_CF01).y;
%       end
%       
%       %Put into plot structure
%       idat=idat+1;
%       xdat(idat).x=[1:12]; %
%       ydat(idat).y=meanNoNan(dat_CF01,2);
%       labs(idat).l='CF 0.1-0.4';
%       line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 3; istyle=istyle+1;
% 
%            
% %CF 0.8-1.0
%        %Gather the ydat data into an array, so can do MeanNoNan
%       dat_CF01 = NaN*ones([12 5]);
%       idat2=0;
%       for iyears_CF01=6:10
%           idat2=idat2+1;
%           dat_CF01(:,idat2) = ydat(iyears_CF01).y;
%       end
%       
%       %Put into plot structure
%       idat=idat+1;
%       xdat(idat).x=[1:12]; %
%       ydat(idat).y=meanNoNan(dat_CF01,2);
%       labs(idat).l='CF 0.8-1.0';
%       line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.6 0.6 0.6]; marker_style(istyle).m='o'; line_widths(istyle).l = 3; istyle=istyle+1;
%       
      
    
    case 150
        % Seasonal Nd cycles for multiple years
        filename = ['saved_vars_for_Nd_multi_year_seasonal_SouthernOcean_L3'];  %for the 'MODIS Nd multiple year seasonal cycles'
        filedir_savevars = '/home/disk/eos1/d.grosvenor/';
        filename_savevars = [filedir_savevars filename '.mat'];

        load(filename_savevars);
        tag = 'Nd_multi_year_seasonal_SouthernOcean_L3';
        
        titlenam = ['Seasonal cycle of Nd from MODIS L3 for 40-60S, 0-360E'];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=[0 100];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        nmark=0;

        ylab='N_d (cm^{-3})';
        xlab= 'Month';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;      
        
        %data is ordered for 2006-2010 for CF 0.1-0.4 (cells 1-5), then 2006-2010 for
        %0.4-0.6 (cells 6-10), etc. (then 0.6-0.8 [cells 11-15], then 0.8-1.0 [cells 16-20])
        
        %see screen_edits_driver for the screening thresholds
        % and screening_eval_str_driver
        

        logflag=0;
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
        
        
        

       % -- first plot the results for CF = 0.1-0.4
       %  These are the first 5 cells in season_mean_driver
       istyle=1;
       idat=0;
       for idat2=1:5
           idat=idat+1;
           
           xdat(idat).x=[1:12]; %
           ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
           labs(idat).l=NaN;
           
           err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
           %Sampling error is not really useful if only have a few points
           %since the std dev is very small!
           N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
           icut = find(N_samp<200);
           ydat(idat).y(icut)=NaN;
           
           line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 0.5; istyle=istyle+1;
           
       end                     
       
       for idat2=16:20
           idat=idat+1;
           
           xdat(idat).x=[1:12]; %
           ydat(idat).y=eval(['season_mean_driver_' tag '{idat2};']); %
           labs(idat).l=NaN;
           
           err_samp = eval(['season_std_driver_' tag '{idat2} ./ sqrt(season_Ndatap_driver_' tag '{idat2});']); %
           %Sampling error is not really useful if only have a few points
           %since the std dev is very small!
           N_samp = eval(['season_Ndatap_driver_' tag '{idat2};']); %
           icut = find(N_samp<200);
           ydat(idat).y(icut)=NaN;
           
           line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.6 0.6 0.6]; marker_style(istyle).m='o'; line_widths(istyle).l = 0.5; istyle=istyle+1;
           

       end

       
       save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],['BASE screening:- ' eval(['screening_eval_str_' tag])] );
       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );
       write_notes_to_file( save_notes_filepath, ['AND :- ' eval(['screen_edits_' tag '{4}']) ] );
       write_notes_to_file( save_notes_filepath, ['FOR REGION :- 40-60S, 0-360E'] );
       
     
              
% --- Now plot the means of all the years   

%CF 0.1-0.4
       %Gather the ydat data into an array, so can do MeanNoNan
      dat_CF01 = NaN*ones([12 5]);
      idat2=0;
      for iyears_CF01=1:5
          idat2=idat2+1;
          dat_CF01(:,idat2) = ydat(iyears_CF01).y;
      end
      
      %Put into plot structure
      idat=idat+1;
      xdat(idat).x=[1:12]; %
      ydat(idat).y=meanNoNan(dat_CF01,2);
      labs(idat).l='CF 0.1-0.4';
      line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.0 0.0 0.0]; marker_style(istyle).m='o'; line_widths(istyle).l = 3; istyle=istyle+1;

           
%CF 0.8-1.0
       %Gather the ydat data into an array, so can do MeanNoNan
      dat_CF01 = NaN*ones([12 5]);
      idat2=0;
      for iyears_CF01=6:10
          idat2=idat2+1;
          dat_CF01(:,idat2) = ydat(iyears_CF01).y;
      end
      
      %Put into plot structure
      idat=idat+1;
      xdat(idat).x=[1:12]; %
      ydat(idat).y=meanNoNan(dat_CF01,2);
      labs(idat).l='CF 0.8-1.0';
      line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.6 0.6 0.6]; marker_style(istyle).m='o'; line_widths(istyle).l = 3; istyle=istyle+1;
      

      
      
        
     case 149
        % Fraction of points for each LWP bin for which >95% of LWP is removed by mphys
        % Run after pdf2d... with teh correct settings (pre-mphys LWP vs
        % fraction of LWP removed)
                      
        titlenam = ['Fraction of points within each LWP bin for which >95% of LWP is removed by mphys'];
        titlenam = [titlenam ' xlab from 2D PDF- ' xlabelstr ', ylab- ' ylabelstr ', case 149 WaterVap'];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[0 100];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Fraction';
        xlab= 'Pre-mphys LWP (g m^{-2})';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;      
        
        fqh95 = qh(1:end-1,1)./sum(qh(1:end-1,1:end-1),2); %the first column of the 2D PDF of pre-mphys LWP vs precentage of LWP removed
        % corresponds to fraction < -95%
        % The last row is ignored because it is filled with dummy values
        % for plotting purposes
        % Number of points are divided by the total number with that LWP
        

        idat=idat+1;
        xdat(idat).x=mid_Ybins; %
        ydat(idat).y=fqh95; %
        labs(idat).l='';

        logflag=1;
        
        
 case 148
        % plots of LWP and CF PDFs for CAMCLUBBv2 and CAM5 to look at the
        % 'empty' cloud issue
        
%These are the filenames that want to load. One for each model etc. Can repeat same thing for all models.      
       plot148 = {'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1',...
           'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.GT.1'...
           };
       
        plot148 = {'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1',...
           'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.TLWP.LT.1'...
           };
       
       plot148 = {'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1_all_times_and_locations',...
           'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.GT.1_all_times_and_locations'...
           };
       
       plot148 = {'TLWP_PDFs_CAPT_1Dpdfs'...
           };
       
       
       
       
       
       gcm_str148 = 'CAMCLUBBv2_prepostLWP';
       gcm_str148 = 'CAM5_prepostLWP';
                            
       
      
        
       
        
        xlims=0;
        
        izlim=0;
        zmin=1500;
        zmax=3000;
        
        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.
        
       


        
        


    
        figname='';
        for idat=1:length(plot148)
            
            filedir_savevars = '/home/disk/eos1/d.grosvenor/';
            filename_savevars = [filedir_savevars plot148{idat} '.mat'];
            load(filename_savevars); %load all the variables
               
            titlenam = eval(['titlenam_' gcm_str148]);
            figname=[figname titlenam];
            savename=figname;
        
            ylab=eval(['ylab_' gcm_str148]);
            xlab=eval(['xlab_' gcm_str148]);

       
            %Look out in case data is stored in the 2nd xdat for older saves. Was due
            %to a bug in the case 977 plot - fixed now for newer saved data
            xdat(idat).x=eval(['xdat_' gcm_str148 '(1).x;']);
            ydat(idat).y=eval(['ydat_' gcm_str148 '(1).y;']);
            
            switch gcm_str148
                case 'CAMCLUBBv2_prepostLWP'
                    titlenam = 'CAMCLUBBv2';
                case 'CAM5_prepostLWP'
                    titlenam = 'CAM5';                    
                otherwise
                    titlenam = gcm_str148;
            end
            
            switch plot148{idat}
                case {'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1','saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1_all_times_and_locations'}
                    xlims=1;
                    xlimits=[xdat(1).x(2) 4.5];
                    x_axis_type = 'log10_matlab';
                    xlab='Precipitation Rate (mm day^{-1})';
%                    titlenam = 'CF >= 0.05';
                    fsize_title = 18; %Make the fontsize for the title large
                    
                    labs(idat).l = 'LWP<1';
                    titlenam = [titlenam ', CF>=0.05'];
                    
                case {'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.GT.1','saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.GT.1_all_times_and_locations'}
                    xlims=1;
                    xlimits=[xdat(1).x(2) 4.5];                    
                    
                     xlab='Precipitation Rate (mm day^{-1})';
%                    titlenam = 'CF >= 0.05';
                    fsize_title = 18; %Make the fontsize for the title large
                    
                    labs(idat).l = 'LWP>=1';
                    titlenam = [titlenam ', CF>=0.05'];
                    
                case 'saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.TLWP.LT.1'
                    xlims=1;
                    xlimits=[xdat(1).x(2) 4.5];
                    x_axis_type = 'log10_matlab';
                    xlab='Precipitation Rate (mm day^{-1})';
%                    titlenam = 'CF >= 0.05';
                    fsize_title = 18; %Make the fontsize for the title large
                    
                    labs(idat).l = 'TLWP<1';
                    titlenam = [titlenam ', CF>=0.05'];  
                    
                case 'TLWP_PDFs_CAPT_1Dpdfs'
                    xlims=0;
                    xlimits=[0 100];
%                    x_axis_type = 'log10_matlab';
                    xlab='LWP (g m^{-2})';
%                    titlenam = 'CF >= 0.05';
                    fsize_title = 18; %Make the fontsize for the title large
                    
                    labs(idat).l = gcm_str148;
                    titlenam = [titlenam ', CF>=0.05'];                     

                    
            end



        end

       



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane


        
     case 147
        % plots of LWP and CF PDFs for CAMCLUBBv2 and CAM5 to look at the
        % 'empty' cloud issue
        xdat_ind = 1; %the index at which the data is stored in xdat - was a bug in the early saved data meaning
        % that it was stored in the 2nd index for some reason
        
       plot147 = 'saved_vars_for_CF_cumulative_PDFs_for_LWP.LT.1'; xdat_ind=2;
       plot147 = 'saved_vars_for_CF_PDFs_for_LWP.LT.1'; xdat_ind=2;
       plot147 = 'saved_vars_for_LWP_PDFs_for_CF.GT.0.05'; xdat_ind=2; %~/modis_work/plots/1D-PDF LWP GCM grid-box mean vs CF GCM, LAT=-39.58 to 9.42, LON=-139.38 to -50.62 CAMCLUBBv2-prepostLWP-CAMCLUBBv2-prepostLWP for CALIPSO-MHCF.GTE.-0.01.AND.LT.0.3, ocean only
       plot147 = 'saved_vars_for_LWP_PDFs_for_CF.GT.0.05_all_times_and_locations'; xdat_ind=2;
       plot147 = 'saved_vars_for_fraction_points_95percent_LWP_removed_mphys_all_points_all_times_and_locations'; xdat_ind=1;
       
       
       gcm_strs = {'CAMCLUBBv2_prepostLWP',...
           'CAM5_prepostLWP'...
           };
       
       filedir_savevars = '/home/disk/eos1/d.grosvenor/';
       filename_savevars = [filedir_savevars plot147 '.mat'];
       load(filename_savevars); %load all the variables
               
       titlenam = eval(['titlenam_' gcm_strs{1}]);
       ylab=eval(['ylab_' gcm_strs{1}]);
       xlab=eval(['xlab_' gcm_strs{1}]);
        
        figname=titlenam;
        savename=figname;
        
        xlims=0;
        
        izlim=0;
        zmin=1500;
        zmax=3000;
        
        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.
        
        switch plot147
            case 'saved_vars_for_LWP_PDFs_for_CF.GT.0.05'
                xlims=1;
                xlimits=[1 1e3];
                x_axis_type = 'log10_matlab';
                xlab='Grid-box mean LWP (g m^{-2})';
                titlenam = 'CF >= 0.05';
                fsize_title = 18; %Make the fontsize for the title large
            case 'saved_vars_for_fraction_points_95percent_LWP_removed_mphys_all_points_all_times_and_locations'
                xlims=1;
                xlimits=[0 100];
                x_axis_type = 'log10_matlab';
%                xlab='Grid-box mean LWP (g m^{-2})';
                titlenam = 'Fraction of points within each LWP bin for which 95% of LWP was removed by mphys';
                fsize_title = 18; %Make the fontsize for the title large
        end

        title_str_notes = eval(['titlenam_' gcm_strs{1}]);
        thresh_str_notes = eval(['thresh_str_' gcm_strs{1}]);
        save_notes = [thresh_str_notes '\n' title_str_notes]; 
        


    
        
        for idat=1:length(gcm_strs)
            %the data is stored in the 2nd xdat for some reason?? Was a bug
            %somewhere. Have fixed for later saved datta.
            xdat(idat).x=eval(['xdat_' gcm_strs{idat} '(xdat_ind).x;']);
            ydat(idat).y=eval(['ydat_' gcm_strs{idat} '(xdat_ind).y;']);
            
            switch gcm_strs{idat}
                case 'CAMCLUBBv2_prepostLWP'
                    labs(idat).l = 'CAMCLUBBv2';
                case 'CAM5_prepostLWP'
                    labs(idat).l = 'CAM5';                    
                otherwise
                    labs(idat).l=gcm_strs{idat};
            end



        end

       



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane



       
        
        
 case 146
        % 
            
        var146 = 'Horizontal wind speed (cross section component, m s^{-1})';                    
        var146 = 'Potential Temperature (K)';              
        var146 = 'w (m s^{-1})';

        
        titlenam = var146;
        
        figname=['Aircraft vs WRF comparison of ' titlenam];
        savename=figname;

        xlims=1;
        xlimits=[-67 -64];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.
        iaxis_square=0;
        
         idat=0;
        
         iz =findheight_nearest(zz(1).z,3.0); 

        idat=idat+1;
        xdat(3).x=lon_slice; %WRF
        xdat(2).x = lon19av;
        xdat(1).x = lon19;        
        
        switch var146
            case 'Horizontal wind speed (cross section component, m s^{-1})'
                ylab='U (m s^{-1})';              
                ydat(3).y = Ucomp_cross_AP(iz,:);
                ydat(2).y = u19av;    
                ydat(1).y = u19;                  

            case 'Potential Temperature (K)'
                ylab = 'Potential Temperature (K)';  
                ydat(3).y = pot_cross_AP(iz,:);
                ydat(2).y = pot19av;
                ydat(1).y = pot19;                
              
            case 'w (m s^{-1})'
                ylab = 'w (m s^{-1})';  
                ydat(3).y = vertwind_AP(iz,:);   
                ydat(2).y = w19av;
                ydat(1).y = w19;                
        end
        
        xlab= 'Time (UTC)';
        
        %Find regions where the potential temperature vertical gradient is
        %low - more likely to be well-mixed
        zrep = (repmat(zz(1).z,[size(pot_cross_AP,2) 1]))';
        dth_dz = diff(pot_cross_AP,1)./diff(zrep,1);
        dth_slice=dth_dz(iz-1,:);
        dth_slice(lon_slice>xlimits(2))=NaN;
        dth_slice(lon_slice<xlimits(1))=NaN;
        iith=find(dth_slice<3);
%        line([lon_slice(min(iith)) lon_slice(min(iith))],[ylims(1),ylims(2)]);
%        line([lon_slice(max(iith)) lon_slice(max(iith))],[ylims(1),ylims(2)]);

%The potemp gradient never actually becomes zero and so perhaps cannot be
%considered to be truly well-mixed. However, it does get lower in the
%suspected well-mixed region - down to around 2 K/km. Might be issues with
%the coarse vertical spacing preventing explicityly resolved zero dth/dz.


        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
      

        labs(3).l = 'WRF';
        labs(2).l = 'Aircraft smoothed';
        labs(1).l = 'Aircraft raw'; 
        
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
                
        istyle=1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0.7 0.7]; marker_style(istyle).m='o'; line_widths(istyle).l = 1; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='^'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0];
%        marker_style(istyle).m='^'; istyle=istyle+1;

        
%        xtickvals_set = [xdat(1).x(1) : 6/24 : xdat(1).x(end)];
        
%        [Y,M,D,H,MM] = datevec(xtickvals_set);
%        month = datestr(xtickvals_set,'mmm');
%        time_str144 = datestr(xtickvals_set,'HH');
        
%        for ixtick_lab=1:length(xtickvals_set)
%            xticklabs_set{ixtick_lab} = [time_str144(ixtick_lab,:)];
%        end
        
%        for ixtick_lab=[1:4:13]
%            xticklabs_set{ixtick_lab} = [num2str(D(ixtick_lab)) '-' month(ixtick_lab,:)];
%        end 
        
%     iset_xticks = 1;
%     iset_xticklabs=1;
        
         
           
        
 case 145
        % SZA vs local time vertical plot for alongside the Nd vs SZA plot
        iaxis_square=0; %switch to make axis square
        
        titlenam='sza vs time of day';
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        xlab='\theta_{0}';
        ylab= 'Local Time';

        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        %calculate sza
        hours=0:0.2:24;
        times=datenum(2007,6,30) + hours/24;
        sza75=sun_pos(times,72,0); %results between the start and end of the period are very slight 
                                   %main variation is north south
        sza75=sun_pos(times,75,0);



        idat=idat+1;
        xdat(idat).x=sza72; %
        ydat(idat).y=hours; %
        labs(idat).l='72^{o} N';
        
         idat=idat+1;
        xdat(idat).x=sza75; %
        ydat(idat).y=hours; %
        labs(idat).l='75^{o} N';

    case 144
        % Timeseries from Fohn cross sections
        % Run make_cross_section_timeseries to load data and make the
        % timeseries (after having run
        % multiple_cross_sections_for_timeseries)
%            idatetick=1;      datetick_type = 'dd-mmm HH';
            
        var144 = 'Horizontal wind speed (cross section component, m s^{-1})';
%        var144 = 'N (s^{-1})';        
%        var144 = 'Non-dimensional mountain height';    
%         var144 = '2km Wind Direction (degrees)';
%         var144 = '1km Wind Direction (degrees)';         
          var144 = '1 km Relative Humidity (%)';          
          var144 = '2 km Relative Humidity (%)';                    
          
        titlenam = var144;
        
        figname=['Timeseries of ' titlenam];
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        iaxis_square=0;
        
         idat=0;
        

        idat=idat+1;
        xdat(idat).x=matlab_time; %
        
        switch var144
            case 'Horizontal wind speed (cross section component, m s^{-1})'
                ylab='U (m s^{-1})';
                ydat(idat).y = U_mean_timser; %
            case 'N (s^{-1})'   
                ylab ='N (s^{-1})';
                ydat(idat).y = N_Brunt_timser;
            case 'Non-dimensional mountain height'  
                ylab = 'Non-dimensional mountain height';  
                ydat(idat).y = Hbar_timser;
            case '2km Wind Direction (degrees)'
                ylab = '2km Wind Direction (degrees)';  
                ydat(idat).y = WindDir_2km_timser;
            case '1km Wind Direction (degrees)'
                ylab = '1km Wind Direction (degrees)';  
                ydat(idat).y = WindDir_1km_timser;    
            case '1 km Relative Humidity (%)'
                ylab = 'Relative humidity (%)';
                ydat(idat).y = RH_1km_timser; 
            case '2 km Relative Humidity (%)'
                ylab = 'Relative humidity (%)';
                ydat(idat).y = RH_2km_timser;                 
        end
        
        xlab= 'Time (UTC)';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
      

        labs(idat).l='';
        
        xtickvals_set = [xdat(1).x(1) : 6/24 : xdat(1).x(end)];
        
        [Y,M,D,H,MM] = datevec(xtickvals_set);
        month = datestr(xtickvals_set,'mmm');
        time_str144 = datestr(xtickvals_set,'HH');
        
        for ixtick_lab=1:length(xtickvals_set)
            xticklabs_set{ixtick_lab} = [time_str144(ixtick_lab,:)];
        end
        
        for ixtick_lab=[1:4:13]
            xticklabs_set{ixtick_lab} = [num2str(D(ixtick_lab)) '-' month(ixtick_lab,:)];
        end 
        
     iset_xticks = 1;
     iset_xticklabs=1;
        
         
        

        
 case 143
        % load in the .mat file first - for names see below.               
       
        titlenam = ['LWP comparison MODIS and AMSRE vs CF ' thresh_str];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='LWP (grid-box average) (g m^{-2})';
        
        switch filename
            case 'saved_vars_for_LWP_vs_CF' %Nd
                xlab= 'MOD06 liquid cloud fraction';
        end


        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        

        
        idat=idat+1;
        ydat(idat).y=xdat_AMSREvsCF(1).x; %
        xdat(idat).x=ydat_AMSREvsCF(1).y; %
%        errordatU(idat).dat = sqrt(errordatU_lowSZA(1).dat.^2 + errordatU_highSZA(1).dat.^2 );
%        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='AMSRE';
        
        idat=idat+1;
        ydat(idat).y=xdat_MODIS_MOD35_15_percent_vsCF(1).x; %
        xdat(idat).x=ydat_MODIS_MOD35_15_percent_vsCF(1).y; %
        %        errordatU(idat).dat = sqrt(errordatU_lowSZA(1).dat.^2 + errordatU_highSZA(1).dat.^2 );
        %        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='MODIS MOD35, minus 15%';
        
        idat=idat+1;
        ydat(idat).y=xdat_MODIS_MOD35_vsCF(1).x; %
        xdat(idat).x=ydat_MODIS_MOD35_vsCF(1).y; %
        %        errordatU(idat).dat = sqrt(errordatU_lowSZA(1).dat.^2 + errordatU_highSZA(1).dat.^2 );
        %        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='MODIS MOD35';
        
        idat=idat+1;
        ydat(idat).y=xdat_MODIS_MOD06_vsCF(1).x; %
        xdat(idat).x=ydat_MODIS_MOD06_vsCF(1).y; %
        %        errordatU(idat).dat = sqrt(errordatU_lowSZA(1).dat.^2 + errordatU_highSZA(1).dat.^2 );
        %        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='MODIS MOD06';
        
        
        
        
        istyle=1; %the counter for the styles we specify

        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='^'; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='^'; istyle=istyle+1;


%           titlenam = [titlenam ' low sensSZA'];
        
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
        
%        ierror_bars='horiz2';  %no error bars - probably should add some
        

        
        
case 142  %PDFs of MODIS and Nacc Nd (Puijo)
    
    pdf_plot142 = 'MODIS vs DMPS Nd matches';
    pdf_plot142 = 'MODIS vs CDP all data';
    pdf_plot142 = 'MODIS L3 vs CDP all data';    
%    pdf_plot142 = 'days of the year for MODIS vs CDP all data';  
  %For the following use load_VOCALS_POLDER_saved_data.m to load in the
  %.mat files
    pdf_plot142 = 'POLDER Reff comparison vs longitude for VOCALS .mat file';
%    pdf_plot142 = 'POLDER Reff comparison vs longitude for VOCALS .mat file LT 20um re16 re37';
    pdf_plot142 = 'POLDER Reff comparison vs longitude for VOCALS .mat file - no confidence screening';
%    pdf_plot142 = 'Het index vs longitude for VOCALS .mat file with POLDER co-location';
    
iswap_xy_142=0;    
mean_or_pdf = 'pdf';

set_colocate_POLDER=1;
i_plot_diff_from_POLDER=0; %May get overridden below
axis1D = 'x'; %Do the mean over x axis (i.e. mean within each y bin)
       
    
    
    nXpdf_131 = 10;
    nYpdf_131 = 10;
    
    %Some dummy default values
    Xbins_set=[0 1];
    Ybins_set=[0 1];
    Zbins_set=[0 1];
    
    pdf2D_defaults
    
    switch pdf_plot142
        case 'Het index vs longitude for VOCALS .mat file with POLDER co-location'
            
            mean_or_pdf = 'mean';
            iswap_xy_142=0;
            
            axis1D = 'y'; %Do the mean over y instead of x
            
            thresh_CF=[-0.1 1.000001];
%            thresh_CF=[0.799 1.000001];
            
            time_mean_str = {'choose_years'}; years_required_for_mean = {[2007]};
            days_required_for_mean = {[1:366]};                        
              
            xvars = {...
%                'Heterogeneity Parameter Cahalan Optical Depth (Seethala)',...
                 'Longitude'...
                };                 
            N142 = length(xvars);            
            
            time_mean_str142 = repmat(time_mean_str,[1 N142]);
            years_required_for_mean142= repmat(years_required_for_mean,[1 N142]);
            days_required_for_mean142 = repmat(days_required_for_mean,[1 N142]);
            
            yvar_str = {'Longitude'};
            yvar_str = {'Heterogeneity Parameter Cahalan Optical Depth (Seethala)'};     
            
            yvars = repmat(yvar_str,[1 N142]);
            
            
            
            
            inotoverride_screening142(1:99) = {0};
            
            screen_str_rep = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};
            screen_type142 = repmat(screen_str_rep,[1 N142]);
            
            labs_multi_specify = {...
                'Heterogeneity',...           
                };
            
%            Xbins_set=[1:10:365];
%            Ybins_set=[1:10:365];
%            Zbins_set=[0:25:2000];



        case 'POLDER Reff comparison vs longitude for VOCALS .mat file LT 20um re16 re37'
            
            mean_or_pdf = 'mean';
            iswap_xy_142=1;
            thresh_CF=[-0.1 1.000001];
%            thresh_CF=[0.799 1.000001];

            inotoverride_screening142(1:99) = {0};
            
            time_mean_str = {'choose_years'}; years_required_for_mean = {[2007]};
            days_required_for_mean = {[1:366]};                        
              
            xvars = {...
              'MockL3 reff 1.6um timeseries3 re.LT.20um',...
              'MockL3 reff 2.1um timeseries3 re.LT.20um',...              
              'MockL3 reff 3.7um timeseries3 re.LT.20um',...
                'POLDER reff daymean',...
                'MockL3 reff 1.6um timeseries3',...
                'MockL3 reff 2.1um timeseries3',...
                'MockL3 reff 3.7um timeseries3',...  
                };                 
            N142 = length(xvars);            
            
            time_mean_str142 = repmat(time_mean_str,[1 N142]);
            years_required_for_mean142= repmat(years_required_for_mean,[1 N142]);
            days_required_for_mean142 = repmat(days_required_for_mean,[1 N142]);
            
            yvar_str = {'Longitude'};
%            yvar_str = {'Heterogeneity Parameter Cahalan Optical Depth (Seethala)'};            
            
               yvars = repmat(yvar_str,[1 N142]);
            
            POLDER_paper_thresh_vals  %Thresh_reff etc.            
            
            screen_str_rep = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};
            screen_type142 = repmat(screen_str_rep,[1 N142]);
            
                
            
            labs_multi_specify = {...
                'MockL3 reff 1.6um timeseries3 re.LT.20um',...
                'MockL3 reff 2.1um timeseries3 re.LT.20um',...                
                'MockL3 reff 3.7um timeseries3 re.LT.20um',...
                'POLDER'...,
                'MODIS mockL3 1.6\mum',...
                'MODIS mockL3 2.1\mum',...                
                'MODIS mockL3 3.7\mum',...  
                };
            
%            Xbins_set=[1:10:365];
%            Ybins_set=[1:10:365];
%            Zbins_set=[0:25:2000];

            for iset=1:length(xvars)
                iset_gcm_str{iset}=1;
                gcm_str_multi{iset} = 'MODIS';
                switchable_str{iset} = 'MODIS.'; %Nothing as default for MODIS, etc.
            end
            
         %Override some of the defaults   
           gcm_str_multi(1:3) = {'switchable'};
%           gcm_str_multi{6} = 'switchable';
%           gcm_str_multi{7} = 'switchable';           

%           switchable_str{1} = 'mockL3_no_conf.'; %add the dot on the end
           switchable_str(1:3) = {'mockL3_conf_re16_re37_LT_20um.'}; %add the dot on the end  




        case 'POLDER Reff comparison vs longitude for VOCALS .mat file'
            
            datatype='timeseries3'; %LAT=MLAT; LON=MLON;
            
            mean_or_pdf = 'mean';
            iswap_xy_142=1;
            thresh_CF=[-0.1 1.000001];
%            thresh_CF=[0.799 1.000001];

            inotoverride_screening142(1:99) = {0};
            
            time_mean_str = {'choose_years'}; years_required_for_mean = {[2007]};
            days_required_for_mean = {[1:366]};                        
            
            %These become x_axis_vals, which is what decides the case in
            %pdf2D_plot_commands
            xvars = {...
              'MockL3_no_conf reff 1.6um timeseries3'...                
              'MockL3_no_conf reff 2.1um timeseries3'...
              'MockL3_no_conf reff 3.7um timeseries3'...              
              'R_{eff 2.1 \mum} (\mum)',...
                'POLDER reff daymean',...
                'MockL3 reff 1.6um timeseries3',...
                'MockL3 reff 2.1um timeseries3',...
                'MockL3 reff 3.7um timeseries3',...  
                };                 
            N142 = length(xvars);            
            
            time_mean_str142 = repmat(time_mean_str,[1 N142]);
            years_required_for_mean142= repmat(years_required_for_mean,[1 N142]);
            days_required_for_mean142 = repmat(days_required_for_mean,[1 N142]);
            
            yvar_str = {'Longitude'};
%            yvar_str = {'Heterogeneity Parameter Cahalan Optical Depth (Seethala)'};            
            
               yvars = repmat(yvar_str,[1 N142]);
            
            POLDER_paper_thresh_vals  %Thresh_reff etc.            
            
            screen_str_rep = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};
            screen_type142 = repmat(screen_str_rep,[1 N142]);
            
                
            
            labs_multi_specify = {...
                'MODIS mockL3 no conf 1.6\mum',...                
                'MODIS mockL3 no conf 2.1\mum',...
                'MODIS mockL3 no conf 3.7\mum',...                
                'MODIS standard L3 2.1\mum',...
                'POLDER'...,
                'MODIS mockL3 1.6\mum',...
                'MODIS mockL3 2.1\mum',...                
                'MODIS mockL3 3.7\mum',...  
                };
            
%            Xbins_set=[1:10:365];
%            Ybins_set=[1:10:365];
%            Zbins_set=[0:25:2000];

            for iset=1:length(xvars)
                iset_gcm_str{iset}=1;
                gcm_str_multi{iset} = 'MODIS';
                switchable_str{iset} = ''; %Nothing as default for MODIS, etc.
            end
            
         %Override some of the defaults   
           gcm_str_multi(1:3) = {'switchable'};
%           gcm_str_multi{6} = 'switchable';
%           gcm_str_multi{7} = 'switchable';           

        %Make it take these values from the following structures
           switchable_str(1:3) = {'mockL3_no_conf.'}; %add the dot on the end
  

        case 'POLDER Reff comparison vs longitude for VOCALS .mat file - no confidence screening'
            i_plot_diff_from_POLDER=0; %Whether to plot the difference from the POLDER line       
            datatype='timeseries3'; %LAT=MLAT; LON=MLON;
            
            mean_or_pdf = 'mean';
            iswap_xy_142=1;
            thresh_CF=[-0.1 1.000001];
%            thresh_CF=[0.799 1.000001];

            inotoverride_screening142(1:99) = {0};
            
            time_mean_str = {'choose_years'}; years_required_for_mean = {[2007]};            
            time_mean_str = {'choose_years'}; years_required_for_mean = {[2005:2012]};
            days_required_for_mean = {[1:366]};                        

%with standard MODIS            
            xvars = {...
              'MockL3_no_conf reff 1.6um timeseries3'...                
              'MockL3_no_conf reff 2.1um timeseries3'...
              'MockL3_no_conf reff 3.7um timeseries3'...              
              'R_{eff 2.1 \mum} (\mum)',...
              'POLDER reff daymean',...
                };    
%without            
            xvars = {...
              'MockL3_no_conf reff 1.6um timeseries3'...                
              'MockL3_no_conf reff 2.1um timeseries3'...
              'MockL3_no_conf reff 3.7um timeseries3'...              
              'POLDER from structure reff daymean',...
                };
            
            N142 = length(xvars);            
            
            time_mean_str142 = repmat(time_mean_str,[1 N142]);
            years_required_for_mean142= repmat(years_required_for_mean,[1 N142]);
            days_required_for_mean142 = repmat(days_required_for_mean,[1 N142]);
            
            yvar_str = {'Longitude'};
%            yvar_str = {'Heterogeneity Parameter Cahalan Optical Depth (Seethala)'};            
            
               yvars = repmat(yvar_str,[1 N142]);
            
            POLDER_paper_thresh_vals  %Thresh_reff etc.            

            %24th Nov, 2016 - had the following screening set, but not sure
            %why for co-located POLDER data?
            screen_str_rep = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};
            %Changed to no screening
            screen_str_rep = {'none'};            
            screen_type142 = repmat(screen_str_rep,[1 N142]);
            
                
            
            labs_multi_specify = {...
                'MODIS mockL3 no conf 1.6\mum',...                
                'MODIS mockL3 no conf 2.1\mum',...
                'MODIS mockL3 no conf 3.7\mum',...                
                'MODIS standard L3 2.1\mum',...
                'POLDER'...,
                };
            labs_multi_specify = {...
                'MODIS mockL3 no conf 1.6\mum',...                
                'MODIS mockL3 no conf 2.1\mum',...
                'MODIS mockL3 no conf 3.7\mum',...                
                'POLDER'...,
                };            
            
%            Xbins_set=[1:10:365];
%            Ybins_set=[1:10:365];
%            Zbins_set=[0:25:2000];

            %Defaults :-
            for iset=1:length(xvars)
                iset_gcm_str{iset}=1;
                gcm_str_multi{iset} = 'MODIS';
                switchable_str{iset} = ''; %Nothing as default for MODIS, etc.
            end
            
         %Override some of the defaults   
           gcm_str_multi(1:4) = {'switchable'};
%           gcm_str_multi{6} = 'switchable';
%           gcm_str_multi{7} = 'switchable';           

           switchable_str(1:4) = {'mockL3_no_conf.'}; %add the dot on the end
  


            
            
        case 'days of the year for MODIS vs CDP all data'
             xvars = {'CDP Nd all data'...
                'Nd from individual pixels timeseries3'...
                };
            
            xvars = {'days for CDP Nd all data'...
                'days for Nd from individual pixels 3.7\mum timeseries3'...
                };
            
               yvars = {'Dummy data for 1D', ...
                'Dummy data for 1D'};         
            
            inotoverride_screening142 = {0,1};
            
            screen_type142 = {'none',''};
            
            labs_multi_specify = {...
                'CDP all'...,
                'MODIS'};
            
            Xbins_set=[1:10:365];
            Ybins_set=[1:10:365];
            Zbins_set=[0:25:2000];
            
        case 'MODIS vs DMPS Nd matches'

            xvars = {'DMPS Nd estimate'...
                'Nd from individual pixels timeseries3'...
                };

            xvars = {'DMPS Nd estimate'...
                'Nd from individual pixels 3.7\mum timeseries3'...
                };
            
            
            yvars = {'Dummy data for 1D', ...
                'Dummy data for 1D'};
            
            inotoverride_screening142 = {1,1};
            
             labs_multi_specify = {...
               'DMPS'...,
               'MODIS'};
           
           Xbins_set=[0:25:2000];
           Ybins_set=[0:25:2000];
           Zbins_set=[0:25:2000];
            

        case 'MODIS vs CDP all data'
            xvars = {'CDP Nd all data'...
                'Nd from individual pixels timeseries3'...
                };
            
            xvars = {'CDP Nd all data'...
                'Nd from individual pixels 3.7\mum timeseries3'...
                };
            
            yvars = {'Dummy data for 1D', ...
                'Dummy data for 1D'};
            
            inotoverride_screening142 = {0,1};
            
            screen_type142 = {'none',''};
            
            labs_multi_specify = {...
                'CDP all'...,
                'MODIS'};
            
            Xbins_set=[0:25:2000];
            Ybins_set=[0:25:2000];
            Zbins_set=[0:25:2000];
            
             inotnan = find(isnan(Date_Time_Swath.timeseries3)==0);
             daynum_timeseries3 = NaN * ones(size(Date_Time_Swath.timeseries3));
             date_str = datestr(Date_Time_Swath.timeseries3(inotnan));        
             daynum_timeseries3(inotnan) = day_of_year_from_date_func(date_str);            
             days_required_for_mean = [275:305]; time_mean_str = 'Oct';

             
        case 'MODIS L3 vs CDP all data'
            xvars = {'CDP Nd all data'...
                'Nd from grid vals timeseries3'...
                };    
            
            yvars = {'Dummy data for 1D', ...
                'Dummy data for 1D'};            
            
            inotoverride_screening142 = {0,1};
            
            screen_type142 = {'none',''};
            
            labs_multi_specify = {...
                'CDP all'...,
                'MODIS L3'};
            
            Xbins_set=[0:25:2000];
            Ybins_set=[0:25:2000];
            Zbins_set=[0:25:2000];
                              
             days_required_for_mean = [275:305]; time_mean_str = 'Oct';

             

    end
        

        

        
       
        
       
    
    %normalise_by_bins = 0;
    
%%  
    clear xdat_multi ydat_multi labs_multi Npdf meanpdf ihtot_142

        for idat_multi=1:length(xvars)
            
            if iset_gcm_str{idat_multi}==1
                gcm_str = gcm_str_multi{idat_multi};
                
                 eval(['LAT = ' switchable_str{idat_multi} 'LAT;']);
                 eval(['LON = ' switchable_str{idat_multi} 'LON;']);
                 
                switch gcm_str
                    case 'MODIS'
                        
                    otherwise                       
                        eval(['daynum_timeseries3_' gcm_str ' = ' switchable_str{idat_multi} 'daynum_timeseries3;']);
                        eval(['modisyear_timeseries3_' gcm_str ' = ' switchable_str{idat_multi} 'modisyear_timeseries3;']);
                end
            end
            

            colocate_POLDER=set_colocate_POLDER;
            ioverride_pdf_varchoose=1; %override the defaults
            ioverride_pdf=1; %override the defaults  
            inotoverride_screening = inotoverride_screening142{idat_multi}; %don't override the screening if set to one
            
            time_mean_str = time_mean_str142{idat_multi};
            years_required_for_mean= years_required_for_mean142{idat_multi};
            days_required_for_mean = days_required_for_mean142{idat_multi};
            
            
            
            if inotoverride_screening==0
                screen_type = screen_type142{idat_multi};
            end
            
            %flags to allow the direct specification of Xbins, Ybins, Zbins
            ichoose_Xbins=1; Xbins=Xbins_set; 
            ichoose_Ybins=1; Ybins=Ybins_set; 
            ichoose_Zbins=1; Zbins=Zbins_set; 
            
                                    
            i577 = 'MODIS_plot_UW';
%            band_str = band_str_multi{idat_multi};


            x_axis_vals = xvars{idat_multi};
            y_axis_vals = yvars{idat_multi};
            
            
%            screen_type = screen_strs{idat_multi};
%            eval(thresh_val_strs{idat_multi});
%            time_mean_str = time_mean_str_multi{idat_multi};
%            days_required_for_mean = days_multi_specify{idat_multi};
%            SST = SST_multi_specify{idat_multi};
            
            man_choose_plotTimeHeight_graph=1;
            logflag=0;
            dlogflag=0;
            
            noplot=1; %flag to say to just do the calculations and not to plot
            
            nXpdf=nXpdf_131;
            nYpdf=nYpdf_131;
            ndims_hist=2;
            ioverride_time_selection=1;
            ioverride_years_time_screen=1;
            plotTimeHeightVap3
%            close(gcf);

            
           
            %plot the PDF
            man_choose_water_graph=1;
            switch mean_or_pdf
                case 'pdf'
                    graph=977; %1d pdf
                case 'mean'
                    graph=966; %mean of 2D PDF - formerly used case 96 - can change back if things don't work
            end
            
            noplot=0;  iwrite_text_dat=1;
            waterVapourMay2005
%            close(gcf);
%            man_choose_water_graph=1;

   %save some values
            xdat_multi(idat_multi).x = xdat(1).x;
            ydat_multi(idat_multi).y = ydat(1).y';
%            labs_multi(idat_multi).l = xvars{idat_multi};               
            labs_multi(idat_multi).l = labs_multi_specify{idat_multi}; 
            if ~exist('N')
                N=NaN;
            end
            Npdf(idat_multi) = N; %total number of datapoints
            meanpdf(idat_multi) = X_mean_overall; %overall mean value
            ihtot_142{idat_multi} = ihtot; %the indices of the values that were used
            errordatL_multi(idat_multi) =  errordatL;
            errordatU_multi(idat_multi) =  errordatU;
            
        end   
        
        noplot=0;
        
        if iswap_xy_142==0
            xdat=xdat_multi;
            ydat=ydat_multi;
        else
            for idat=1:length(xdat_multi)
                xdat(idat).x = ydat_multi(idat).y;
                ydat(idat).y = xdat_multi(idat).x;
            end
            ylab_save = ylab;
            ylab = xlab;
            xlab = ylab_save;
        end
        
       
        
        labs=labs_multi;
        errordatL = errordatL_multi;
        errordatU = errordatU_multi;
        
        if i_plot_diff_from_POLDER==1
            ipol=0;
            for idat=1:length(xdat_multi)
               if length(strfind(labs(idat).l,'POLDER'))>0
                   ipol=idat;
               end
            end
            for idat=1:length(xdat_multi)
                ydat(idat).y = ydat(idat).y - ydat(ipol).y;
            end
        end        
    
   
    grid on
    title(['time=' num2str(itime)]);
    
%     titlenam = [data_type_cf_cal ' ' norm_str ' frequency for LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(end)) ' LON=' num2str(LON_val(1)) ' to ' num2str(LON_val(end))];
        
        figname=titlenam;
        savename=figname;


        xlims=0;
        xlimits=[0 1000];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ierror_bars='vert2';
%        ierror_bars='none'; 




        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
       istyle=1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
        

% end case 142

%%        
    case 141 %
        dataset = 'Pre Sep 2014';
        dataset = 'Post Sep 2014 - using VOCALS .mat file';
        
        plot_global_maps_defaults

        ichoose_styles=0;
        
        
        idat_141=0;

        
        recalc=1;
        
        
        if recalc==1
            clear xdat_save ydat_save std_dev_save Ndat_save sumsq_save
            
            switch dataset
                case 'Pre Sep 2014'
        
        idat_141=idat_141+1;
        modis_data_plot='POLDER effective radius - using colocated data'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'POLDER';
        
        
        idat_141=idat_141+1;        
        modis_data_plot='Cloud Effective Radius 1.6\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'MODIS 1.6\mum';

        
        idat_141=idat_141+1;
        modis_data_plot='Cloud Effective Radius 2.1\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'MODIS 2.1\mum';

        
        idat_141=idat_141+1;
        modis_data_plot='Cloud Effective Radius 3.7\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'MODIS 3.7\mum';
        
    case 'Post Sep 2014 - using VOCALS .mat file'
        idat_141=idat_141+1;
        modis_data_plot='POLDER effective radius - using individual days from VOCALS save file'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'POLDER';
        
        
        idat_141=idat_141+1;        
        modis_data_plot='Re1.6 mockL3'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        ico_locate_POLDER=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'MODIS 1.6\mum';

        
        idat_141=idat_141+1;
        modis_data_plot='Re2.1 mockL3'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        ico_locate_POLDER=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'MODIS 2.1\mum';

        
        idat_141=idat_141+1;
        modis_data_plot='Re3.7 mockL3'; mod_data_type='timeseries3 lambert';
        ioverride_plotglobal_thresh=1;
        ico_locate_POLDER=1;
        plot_global_maps
        titlenam_save = ['mean values vs longitude ' titlenam];
        man_choose_water_graph=1;
        graph=115;
        watervapourMay2005
        xdat_save(idat_141) = xdat; ydat_save(idat_141) = ydat; std_dev_save(idat_141)=std_dev_gcm; Ndat_save(idat_141) = Ndat_gcm; sumsq_save(idat_141) = sumsq;
        labs_save(idat_141).l = 'MODIS 3.7\mum';
        
            end
                    

        
        end
        
                 
        ylab = 'Effective Radius (\mum)';
        xlab= 'Longitude';
            
        izlim=1;
        zmin=0;
        zmax=20;
        
        clear xdat ydat
nbins=4;
        for idat=1:length(xdat_save)
            ni = length(xdat_save(idat).x)/nbins;
            for i=1:(floor(ni)-1)
                inds = (i-1)*nbins+1:(i-1)*nbins+1+3;
                Ns = Ndat_save(idat).N(inds);
                me = meanNoNan(ydat_save(idat).y(inds).*Ns,2,'sum')./sum(Ns);
                xdat(idat).x(i)=meanNoNan(xdat_save(idat).x(inds),2);
                ydat(idat).y(i)=me;
                N = meanNoNan(Ns,2,'sum');
                Ndat(idat).N(i) = N;
%                sum_std = meanNoNan(Ndat_save(idat).N(inds).*std_dev_save(idat).s(inds),2,'sum');
                
%                std_dev(idat).dat(i) = ( sum_std./Ndat(idat).N(i) );
%using formula variance = N/(N-1) * { mean(x.^2) - (mean(x)).^2 }
                std_dev(idat).dat(i) = sqrt ( N./(N-1) .* ( sum(sumsq_save(idat).y(inds))./N - me.^2 ) );
                errordatU(idat).dat(i) = std_dev(idat).dat(i) ./ sqrt(N);
            end
            labs(idat).l=labs_save(idat).l;
            errordatL(idat).dat = errordatU(idat).dat;
            
        end
        titlenam = titlenam_save;
%        y_axis_type = 'log10_matlab';


lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
figname=titlenam;
savename=[savedir figname];

xlims=0;
xlimits=1000*[0 0.025];


nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        
        ierror_bars='vert2';       
        ichoose_styles=1;
        
        idat_141=1;
        line_pattern(idat_141).p= '-';  line_colour(idat_141).c=[1 0 0]; marker_style(idat_141).m='d';
        idat_141=idat_141+1;
        line_pattern(idat_141).p= '-';  line_colour(idat_141).c=[0 0 0]; marker_style(idat_141).m='o';
        idat_141=idat_141+1;
        line_pattern(idat_141).p= '-';  line_colour(idat_141).c=[0 0 1]; marker_style(idat_141).m='o';
        idat_141=idat_141+1;
        line_pattern(idat_141).p= '-';  line_colour(idat_141).c=[0.6 0.6 0.6]; marker_style(idat_141).m='o';
        idat_141=idat_141+1;
        
        
        
 case 140
        % PDFs of Tau etc. from GCMs (COSP) and MODIS at high CF
        % Uses Xbins_XXX and PDF_XXX, where XXX is MODIS, CAM5 etc.
        % Think I created these using pdf2d routine, then created 1D PDF
        % using watervap case 977 and saved them in a mat file ready to be
        % loaded all at once
        
       pdf_plot_type = 'PDF';
%       pdf_plot_type = 'Cumulative';

       plot_case_140 = 'TLWP PDFs CPT';
       plot_case_140 = 'LWP PDFs CPT';    %For CAM models   
       plot_case_140 = 'LWP PDFs CPT2';   %For AM3  
       plot_case_140 = 'AMSRE clear-sky bias LWP';
       
       iload=1;
       iuse_style_file=1;
       
       lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
       
       switch plot_case_140
           case {'AMSRE clear-sky bias LWP'}
               load_path ='/home/disk/eos1/d.grosvenor/saved_vars/LWP_PDFs_AMSRE_clear-sky_bias_tests_1Dpdfs.mat';
               iuse_style_file=1;
                style_file='Styles_clear_sky_LWP_PDFs_22ndSep2015';
               
           case 'TLWP PDFs CPT'       
               load_path = '/home/disk/eos1/d.grosvenor/saved_vars/TLWP_PDFs_CAPT_1Dpdfs.mat';
               style_file='Styles_TLWP_PDFs_11Oct2014';  %Uses a matlab script to specify the styles for each line
               %The best way to plot seems to be as a normal (i.e. not
               %cumulative) PDF with log axes for both - this then shows up
               %the differences at both ends of the LWP scale.
           case {'LWP PDFs CPT'}
               load_path ='/home/disk/eos1/d.grosvenor/saved_vars/LWP_PDFs_CAPT_1Dpdfs.mat';
               style_file='Styles_TLWP_PDFs_11Oct2014';  %Uses a matlab script to specify the styles for each line  
               
          case {'LWP PDFs CPT2'} %For AM3 use a different file where the AMSRE data has been averaged to 2x2 degrees
               load_path ='/home/disk/eos1/d.grosvenor/saved_vars/LWP_PDFs_CAPT_ocean_only_2deg_1Dpdfs.mat';
               style_file='Styles_TLWP_PDFs_11Oct2014';  %Uses a matlab script to specify the styles for each line  
           
               
       end
       
       if iload==1
           load(load_path);
       end
       



        
       switch plot_case_140
           case 'AMSRE clear-sky bias LWP'
               plot140_temp = {'LWP_day'}; gcm_str140_temp={'AMSRE_time3_0pt01','AMSRE_time3_0pt03','AMSRE_time3_0pt05','AMSRE_time3_0pt07','AMSRE_time3_0pt09','AMSRE_time3_0pt2','AMSRE_time3_1pt01'}; plot140 = repmat(plot140_temp,[1 length(gcm_str140_temp)]); plot_type140='multi_both';
               plot140_temp = {'LWP_day'}; gcm_str140_temp={'AMSRE_time3_0pt03','AMSRE_time3_0pt05','AMSRE_time3_0pt07','AMSRE_time3_0pt09','AMSRE_time3_0pt2'}; plot140 = repmat(plot140_temp,[1 length(gcm_str140_temp)]); plot_type140='multi_both';               
%               plot140_temp = {'LWP_day'}; gcm_str140_temp={'AMSRE_time3_0pt03','AMSRE_time3_0pt05','AMSRE_time3_0pt07','AMSRE_time3_0pt09','AMSRE_time3_0pt2','AMSRE_time3_0pt8'}; plot140 = repmat(plot140_temp,[1 length(gcm_str140_temp)]); plot_type140='multi_both';                              

           case 'TLWP PDFs CPT'
               %Case where want multiple variables and models. Just replicate the
               %variables for eeach model.
               plot140_temp = {'TLWP_day','TLWP_night'}; gcm_str140_temp={'AMSRE_time3','CAM5_prepostLWP','CAMCLUBBv2_prepostLWP'}; plot140 = repmat(plot140_temp,[1 length(gcm_str140_temp)]); plot_type140='multi_both';




           case {'LWP PDFs CPT'}
               plot140_temp = {'LWP_day','LWP_night'}; gcm_str140_temp={'AMSRE_time3','CAM5_prepostLWP','CAM5_CLUBB_COSP','CAMCLUBBv2_prepostLWP'}; plot140 = repmat(plot140_temp,[1 length(gcm_str140_temp)]); plot_type140='multi_both';
               
           case {'LWP PDFs CPT2'}
               plot140_temp = {'LWP_day','LWP_night'}; gcm_str140_temp={'AMSRE_time3','AM3','AM3_CLUBBv1_2deg','AM3_CLUBBv2_COSP_200km'}; plot140 = repmat(plot140_temp,[1 length(gcm_str140_temp)]); plot_type140='multi_both';               
               lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
           otherwise



               plot140 = {'Tau'};
               plot140 = {'LWP'};
               %plot140 = {'Nd'};
               %plot140 = {'CF'};
               %plot140 = {'LWP_COSP'};
               %plot140 = {'Re'};
 %              plot140 = {'LWP_in_cloud'};
 %              plot140 = {'LWP_in_cloud_COSP_MODIS_CF'};
               %plot140 = {'CALISPO_Daytime','MOD06_CF','MOD35_CF','Model_CF','COSP_CALIPSO_CF','COSP_MODIS_CF'}; gcm_str140={'CALIPSO_monthly','MODIS','MODIS',gcm_str gcm_str gcm_str};
               %plot140 = {'MOD35_CF','Model_CF','COSP_CALIPSO_CF','COSP_MODIS_CF'}; gcm_str140={'MODIS',gcm_str gcm_str gcm_str}; plot_type140='multi_variable';
               %plot140 =
               %{'LWP_in_cloud_AMSRE_MOD35','LWP_in_cloud_COSP_MODIS_CF','LWP_in_cloud_model_CF','TLWP_in_cloud_COSP_MODIS_CF'}; gcm_str140={'MODIS',gcm_str gcm_str gcm_str};
%               plot140 = {'Precip_rate','Precip_rate'}; gcm_str140={'CAM5_COSP','CAM5_CLUBBv2_COSP'}; plot_type140='multi_model';



       end

%The code below replicates the model names for each variable as e.g. {model1
%model2 model1 model2 model1 model2} for the case where we have three variables and two
%models
switch plot_type140
    case 'multi_both'   
        Nmodels = length(gcm_str140_temp);
        Nplots_temp = length(plot140_temp);
        clear gcm_str140
        for i=1:Nmodels
            ind0 = (i-1)*Nplots_temp+1;
            inds = ind0:ind0+Nplots_temp-1;
            temp = repmat(gcm_str140_temp(i),[1 Nplots_temp]);
            gcm_str140(inds) = temp;
        end

end


Nvars140 = length(plot140);



ireduce_modis = 0;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-0; %-1 means that all points have markers. Otherwise only plot the number specified.
        



switch plot140{1}   
    case 'TLWP_day'
        titlenam = ['TLWP PDFs for ' LAT_str ' and ' LON_str];
        xlab= 'TLWP (g m^{-2})';
        xlims=1; xlimits=[0 200];
        
        x_axis_type = 'log10_matlab'; y_axis_type = 'log10_matlab'; xlimits=[0 700]; izlim=1;  zmin=5e-5; zmax=5e-2; lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
    case 'LWP_day'
                
        switch plot_case_140
           case 'AMSRE clear-sky bias LWP'               
               titlenam = ['LWP PDFs for ' LAT_str ' and ' LON_str];
               xlab= 'LWP (g m^{-2})';
               xlims=1;

               x_axis_type = ''; y_axis_type = ''; xlimits=[-50 50]; izlim=0;  zmin=5e-5; zmax=5e-2; lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane


            otherwise
                titlenam = ['LWP PDFs for ' LAT_str ' and ' LON_str];
                xlab= 'LWP (g m^{-2})';
                xlims=1; xlimits=[0 200];

                x_axis_type = 'log10_matlab'; y_axis_type = 'log10_matlab'; xlimits=[0 700]; izlim=1;  zmin=5e-5; zmax=5e-2; %lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        end
        
        
    case 'Tau'        
        titlenam = ['Optical Depth PDFs for ' LAT_str ' and ' LON_str];
        xlab= 'COSP Optical Depth';
        pdf_str = 'Tau'; 
 case 'Re'        
        titlenam = ['Effective radius PDFs for ' LAT_str ' and ' LON_str];
        xlab= 'COSP Effective radius (\mum)';
        pdf_str = 'Re';     
         xlims=1;
        xlimits=[3 30];
    case 'LWP'
        titlenam = ['LWP PDFs for ' LAT_str ' and ' LON_str];
         xlab= 'LWP (g m^{-2})';
        pdf_str = 'LWP';         
    case 'Nd'
        titlenam = ['Droplet Concentration PDFs for ' LAT_str ' and ' LON_str];
        xlab= 'Droplet Concentration (cm^{-3})';
        pdf_str = 'Nd';        
    case 'CF'
        titlenam = ['Cloud Fraction PDFs for ' LAT_str ' and ' LON_str];
        xlab= 'Cloud Fraction';
        pdf_str = 'CF';  
     case 'LWP_COSP'
        titlenam = ['LWP COSP (re*1.6) PDFs for ' LAT_str ' and ' LON_str];
         xlab= 'LWP COSP (g m^{-2})';
        pdf_str = 'LWP_COSP';   
     case 'LWP_in_cloud'
        titlenam = ['In-cloud LWP PDFs for ' LAT_str ' and ' LON_str];
         xlab= 'In-cloud LWP (g m^{-2})';
        pdf_str = 'LWP_in_cloud'; 
    case 'LWP_in_cloud_COSP_MODIS_CF'
        titlenam = ['In-cloud LWP PDFs for ' LAT_str ' and ' LON_str ', CF GTE ' num2str(100*CF_gcm_thresh(1),'%2.0f') ' AND LT ' num2str(100*CF_gcm_thresh(2),'%2.0f') ' percent'];
        xlab= 'In-cloud LWP (g m^{-2})';
        pdf_str = 'LWP_in_cloud_COSP_MODIS_CF';

    case 'Model_CF'
        titlenam = ['Cloud Fraction PDFs for ' LAT_str ' and ' LON_str];
         xlab= 'Cloud Fraction';     
    case 'MOD35_CF'
        titlenam = ['Cloud Fraction PDFs for ' LAT_str ' and ' LON_str];
         xlab= 'Cloud Fraction'; 
    case 'LWP_in_cloud_AMSRE_MOD35'
         titlenam = ['In-cloud normalised LWP PDFs for ' LAT_str ' and ' LON_str ', CF GTE ' num2str(100*CF_gcm_thresh(1),'%2.0f') ' AND LT ' num2str(100*CF_gcm_thresh(2),'%2.0f') ' percent'];
         xlab= 'LWP (g m^{-2})';
    otherwise
        titlenam = 'PDF';
        xlab = plot140{1};
        pdf_str = plot140{1};
end

titlenam = [titlenam ' ' time_mean_str];

        figname=titlenam;
        savename=figname;
        
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours

       
        istyle=1;



        idat=0;
        clear labs
        
        for ivar140=1:Nvars140
        
            if Nvars140==1

                if ireduce_modis==1
                    idat=idat+1;
                    data_str = 'MODIS';
                    labs(idat).l='MODIS minus 20%';
                    case140_calc_xdat_ydat
                    xdat(idat).x = xdat(idat).x*0.8;
                    line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                end

                idat=idat+1;
                data_str = 'MODIS';
                pdf_str_bk = pdf_str;
                pdf_str ='LWP_in_cloud';
                data_str = 'AMSRE_MOD35_MODIS';
                %        labs(idat).l='MODIS';
                labs(idat).l='AMSRE daytime ';
                case140_calc_xdat_ydat
                pdf_str = pdf_str_bk;

                %models
                idat=idat+1;
                data_str = 'CAM5';
                data_str = 'CAM5_COSP';
                labs(idat).l='CAM5';
                case140_calc_xdat_ydat


                idat=idat+1;
                data_str = 'CAM5_CLUBB';
                data_str = 'CAM5_CLUBB_COSP';
                labs(idat).l='CAM5-CLUBB';
                case140_calc_xdat_ydat

                idat=idat+1;
                data_str = 'CAM5_CLUBBv2_COSP';
                labs(idat).l='CAM5-CLUBBv2';
                case140_calc_xdat_ydat


                %         idat=idat+1;
                %         data_str = 'AM3';
                %         data_str = 'AM3_COSP';
                %         xdat(idat).x=eval(['0.5*(Xbins_' pdf_str '_' data_str '(1:end-1) + Xbins_' pdf_str '_' data_str '(2:end));']); %
                %         ydat(idat).y=eval(['PDF_' pdf_str '_' data_str '(1:end-1)./sum(PDF_' pdf_str '_' data_str '(1:end-1));']); %
                %         labs(idat).l='AM3';
                %         mean_pdf(idat) = sum(ydat(idat).y.*xdat(idat).x);

                %         idat=idat+1;
                %         data_str = 'POLDER';
                %         xdat(idat).x=eval(['0.5*(Xbins_' pdf_str '_' data_str '(1:end-1) + Xbins_' pdf_str '_' data_str '(2:end));']); %
                %         ydat(idat).y=eval(['PDF_' pdf_str '_' data_str '(1:end-1)./sum(PDF_' pdf_str '_' data_str '(1:end-1));']); %
                %         labs(idat).l='POLDER';
                %         mean_pdf(idat) = sum(ydat(idat).y.*xdat(idat).x);


            else %if Nvars140==1 

                idat=idat+1;
                pdf_str = plot140{ivar140};
                data_str = gcm_str140{ivar140};
                
                switch data_str
                    case 'CAM5_prepostLWP'
                        data_str_nice = 'CAM5';
                    case 'CAMCLUBBv2_prepostLWP'
                        data_str_nice = 'CAMCLUBBv2';   
                    case 'AMSRE_time3'
                        data_str_nice = 'AMSRE';                        
                    otherwise
                        data_str_nice = data_str;
                end
                
                switch plot_type140
                    case 'multi_model'
                        labs(idat).l = remove_character(data_str_nice,'_',' ');
                    case 'multi_variable'
                        labs(idat).l = remove_character(pdf_str,'_',' ');
                    case 'multi_both'
                        labs(idat).l = [remove_character(pdf_str,'_',' ') ' ' remove_character(data_str_nice,'_',' ')];
                end
                
                
% External code to get the required data
                case140_calc_xdat_ydat


            end %if Nvars140==1
        
        if Nvars140==1

            for idat2=1:length(labs)
                switch labs(idat2).l
                    case {'MODIS  ','MODIS','CALIPSO monthly nighttime ','CALIPSO monthly NIGHTTIME ','CALIPSO monthly DAYTIME ','CALIPSO monthly average ',' CLOUDSAT PRECIP  ','ascending CLOUDSAT PRECIP  ','descending CLOUDSAT PRECIP  ','AMSRE daytime ','AMSRE nighttime ','AMSRE average '}
                        line_pattern(istyle).p= '-.';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                        if exist('man_choose_water_graph') & man_choose_water_graph==1
                            switch labs(idat).l
                                case {'ascending CLOUDSAT PRECIP  ',' CLOUDSAT PRECIP  ','descending CLOUDSAT PRECIP  '}
                                    labs(idat2).l = 'CloudSat';
                                case {'AMSRE daytime ','AMSRE nighttime ','AMSRE average '}
                                    labs(idat2).l = 'AMSRE';
                            end
                        end


                    case {'CAM5  CAM5 1deg','CAM5 COSP  CAM5 1deg','CAM5'}
                        labs(idat2).l = 'CAM5 1deg';
                        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                    case {'CAM5 CLUBB  CAMCLUBB','CAM5 CLUBB COSP  ','CAM5-CLUBB'}
                        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'CAM5 CLUBB';
                    case {'CAM5-CLUBBv2'}
                        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'CAM5 CLUBBv2';
                    case 'CAM5  CAM5 2deg'
                        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.6 0.6 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'CAM5 2deg';
                    case {'AM3  2deg','AM3'}
                        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'AM3 2deg';
                    case 'AM3  0pt5deg'
                        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'AM3 0.5deg';
                    case 'AM3 CLUBB  AM3CLUBB'
                        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'AM3 CLUBB';
                    case 'AM3 CLUBBv2'
                        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'AM3 CLUBB';
                    case 'POLDER'
                        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 1]; marker_style(istyle).m='d'; istyle=istyle+1;
                        labs(idat2).l = 'POLDER';

                end

            end

        else

            ichoose_styles=0;


        end %if Nvars140==1




        end


        if iuse_style_file==1
            eval(style_file);  %Uses a matlab script to specify the styles for each line
        end
 
 %quick fix as forgot to normalise by bin widths... fix and remove later
%  for idat=1:length(xdat)
%      dx = mean(diff(xdat(idat).x))'
%      ydat(idat).y = ydat(idat).y / dx;
%  end
 
 
 
 
 
        
       
        
    case 139 
        %run 'Percentage Nd difference ALL vs low SZA - specific days' from
        %plot_global_maps first
%for the tall skinny graph set 
        iaxis_square=0;
        %and after plotting do:-
%        set(gca,'position',[0.3280    0.1100    0.170    0.68141]) 
%         also set(gca,'yticklabels','') to take away ytick labels
  %lso makes it the same height as the plot_global plot (in lat dimension)
  
    
        xlims=1;
        xlimits=[0 20];
        
        ierror_bars='vert2';
        ierror_bars='none';        

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        nmark=0;
        
        titlenam = ['Mean percentage difference in Nd for low vs all SZA for ' time_mean_str];
       
        
        xlab = 'Percentage difference';
        ylab = 'Latitude';
        ylab = 'Miller y co-ordinate'; 
        
        [x,ylim_miller]=m_ll2xy(0,82);
        [x,yticks_miller]=m_ll2xy(0,[-60 -30 0 30 62])
       
        izlim=1;
        zmin=-ylim_miller;
        zmax=ylim_miller;

                
                
        figname=[ylab titlenam];
        savename=[savedir figname];
                
       

%        y_axis_type = 'log10_matlab';


        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
       idat=0;


       idat = idat+1;       
%       ydat(idat).y = MLAT;     
       ydat(idat).y = 5/4*log(tan(pi/4+2/5*MLAT*pi/180)); %miller cylindrical projection- not linear!
       xdat(idat).x = meanNoNan(meanNoNan(dat_modis(:,:,time_inds_average),2),2);       
       labs(idat).l='';
       
       %N.B. or can do e.g. [x,y]=m_ll2xy(0,82), which calculates x and y
       %for the current projection
        
       
        

        
        xlims=0;
        xlimits=[0 13];
        


        
case 138  %
    
    LAT_val = [25 34]; LON_val = [122:129]; %China Sea
    LAT_val = [44:51]; LON_val = [-30:-6]; %SW UK/Atlantic - ocean only
    LAT_val = [40:50]; LON_val = [30:60]; %Central Europe    
    LAT_val = [-43 -40]; LON_val = [110 144]; %west of Tasmania (similar to Boers, QJRMS, 1998)
    LAT_val = [72 75]; LON_val = [-3:48]; %Arctic summer box region
       
    
    xvars = {'Percentage Nd difference allSZA vs lowSZA'...
             'Percentage Nd difference allSZA vs lowSZA'...
             'Percentage Nd difference allSZA vs lowSZA'...
            };
        
        band_str_multi = {'16'...
            '21'...
            '37'...
            };
    
    %normalise_by_bins = 0;

        for idat_multi=1:length(xvars)
        
            ioverride_pdf_varchoose=1; %override the defaults
            ioverride_pdf=1; %override the defaults  
            
            %flags to allow the direct specification of Xbins, Ybins, Zbins
            ichoose_Xbins=ichoose_Xbins_set; Xbins=Xbins_set; 
            ichoose_Ybins=ichoose_Ybins_set; Ybins=Ybins_set; 
            ichoose_Zbins=ichoose_Zbins_set; Zbins=Zbins_set; 
            
                                    
            i577 = 'MODIS_plot_UW';
            band_str = band_str_multi{idat_multi};
            x_axis_vals = xvars{idat_multi};
            screen_type = screen_strs{idat_multi};
            eval(thresh_val_strs{idat_multi});
            time_mean_str = time_mean_str_multi{idat_multi};
            days_required_for_mean = days_multi_specify{idat_multi};
%            SST = SST_multi_specify{idat_multi};
            
            man_choose_plotTimeHeight_graph=1;
            logflag=0;
            dlogflag=0;
            noplot=1; %flag to say to just do the calculations and not to plot
            nXpdf=nXpdf_131;
            nYpdf=N_sza_bins;
            ndims_hist=2;
            ioverride_time_selection=1;
            plotTimeHeightVap3
%            close(gcf);

            
           
            %plot the PDF
            man_choose_water_graph=1;
            switch mean_or_pdf
                case 'pdf'
                    graph=977; %1d pdf
                case 'mean'
                    graph=96; %mean of 2D PDF
            end
            
            noplot=1;
            waterVapourMay2005
%            close(gcf);
            man_choose_water_graph=1;

            
            

            xdat_multi(idat_multi).x = xdat(1).x;
            ydat_multi(idat_multi).y = ydat(1).y';
%            labs_multi(idat_multi).l = xvars{idat_multi};               
            labs_multi(idat_multi).l = labs_multi_specify{idat_multi};   

  
        end    

   
    
   
    grid on
    title(['time=' num2str(itime)]);
    
     titlenam = [data_type_cf_cal ' ' norm_str ' frequency for LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(end)) ' LON=' num2str(LON_val(1)) ' to ' num2str(LON_val(end))];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[-2 10];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.






        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
       istyle=1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
        
        
        



        
case 137 %plots of max vs min sza for a year of Terra data
        
lat_sza = 70;

[tmp,ilat]=min(abs(MLAT-lat_sza));
lat_sza_str = num2str(MLAT(ilat));

    
        xlims=0;
        xlimits=1000*[0 0.025];
        
        ierror_bars='vert2';
        ierror_bars='none';        

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        nmark=0;
        
        titlenam = ['Solar Zenith Angle properties for ' lat_sza_str '^{o}'];
       
        
        xlab = 'Day of year';
        ylab = 'Solar Zenith Angle (degrees)';
       
        izlim=1;
        zmin=0;
        zmax=90;

                
                
        figname=[ylab titlenam];
        savename=[savedir figname];
                
       

%        y_axis_type = 'log10_matlab';


        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
       idat=0;


       idat = idat+1;       
       %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       ydat(idat).y = squeeze(min(Solar_Zenith_Maximum.timeseries3(ilat,:,:),[],2));
   %so instead will plot the max SZA for one location and show how this
   %varies with time. Will give some idea of the range of the lowest SZA
   %available across all lons (excluding the small regions just mentioned).
       ydat(idat).y = squeeze((Solar_Zenith_Maximum.timeseries3(ilat,1,:)));       
       xdat(idat).x = daynum_timeseries3;       
       labs(idat).l='Zonal min of daily max';
        
       
        

        
        xlims=0;
        xlimits=[0 13];
        
        
    
    case 136 %timeseries from values from monthly_means_from_plot_global script       
        

        xlims=0;
        xlimits=1000*[0 0.025];
        
        ierror_bars='vert2';
%        ierror_bars='none';        

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        titlenam = ['Monthly cycle for ' season_LAT_str{1} ' ^{o} N, ' season_LON_str{1} ' ^{o}E for ' thresh_str];
       
        

        xlab = xlab_str;
        ylab = season_datatype{1};
        %This bit just changes the y limits based on the ylab already set
        %in PlotTimeHeightVap3
                switch season_datatype{1}
                    case {'Number of droplets cell values time mean - specific days'}
                        izlim=1;
                        zmin=0;
                        zmax=100;                        
                        ylab = 'Droplet concentration (cm^{-3})';                         
                    case 'Liquid Cloud Fraction time mean - specific days'
                        ylab = 'MODIS Liquid Cloud Fraction';
                        izlim=1;
                        zmin=0;
                        zmax=1;
                    case 'Ice Cloud Fraction time mean - specific days'
                        ylab = 'MODIS Ice Cloud Fraction';
                        izlim=1;
                        zmin=0;
                        zmax=1;
                    case 'Combined Cloud Fraction time mean - specific days'
                        ylab = 'MODIS Combined Cloud Fraction';
                        izlim=1;
                        zmin=0;
                        zmax=1;    
                    case 'Low Cloud Fraction COSP-CALIPSO'
                        ylab = 'CALIPSO Low Cloud Fraction';
                        izlim=1;
                        zmin=0;
                        zmax=1;
                        titlenam = [titlenam ' for ' years_calipso_str];
                    case 'Mid Cloud Fraction COSP-CALIPSO'
                        ylab = 'CALIPSO Mid Cloud Fraction';
                        izlim=1;
                        zmin=0;
                        zmax=1;
                        titlenam = [titlenam ' for ' years_calipso_str];
                    case 'High Cloud Fraction COSP-CALIPSO'
                        ylab = 'CALIPSO High Cloud Fraction';
                        izlim=1;
                        zmin=0;
                        zmax=1;
                        titlenam = [titlenam ' for ' years_calipso_str];                        
                end
                
                 figname=[ylab titlenam];
                 savename=[savedir figname];
                
       

%        y_axis_type = 'log10_matlab';


        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
        for idat=1:1  %length(time_means_multi)


%        xdat(idat).x=X_mean; %
       
        ydat(idat).y = season_mean; %
        switch multi_case
            case 'Bimonthly'
                xdat(idat).x = [1.5:2:11.5];
            otherwise               
                xdat(idat).x = 1:length(ydat(idat).y);
        end
        
%        labs(idat).l=['LON=' LON_str_multi(idat).dat, ',LAT=' LAT_str_multi(idat).dat ' me=' num2str(Xmean_multi(idat).dat,'%3.2f')];
        labs(idat).l='';
        
        error_type = 'Sampling error';
        error_type = 'Spatial std dev of all monthly means';
%        error_type = 'Spatial and temporal std dev sampling error';
                
        switch error_type
            case 'Sampling error'
                sampling_errorU(idat).dat = season_std ./ sqrt(season_Ndatap);
            case 'Spatial std dev of all monthly means'
                sampling_errorU(idat).dat = season_std; %spatial std of the monthly means for all grid points
            case 'Spatial and temporal std dev sampling error'
                sampling_errorU(idat).dat = season_stdoverall ./ sqrt(season_Noverall);
        end
        
        errordatU(idat).dat = sqrt (  sampling_errorU(idat).dat.^2 );
        errordatL(idat).dat = sqrt (  sampling_errorU(idat).dat.^2 );
        
        Ndatapt_thresh = 300;
        irem_pts = find(season_Ndatap<Ndatapt_thresh);
        xdat(idat).x(irem_pts)=NaN;
        
        end
        
        xlims=1;
        xlimits=[0 13];
        
        
        
    case 135
        % bar plot of diurnal ranges for longitude bins
        %Values are stored in the array value_array{ilat,ivar,iline}
        %ilat is the index of the latitude, ivar is the variable and iline
        %is the line of the plot saved (e.g. Cloudsat, Am3, etc)
        % value_array is made in -- saveload_plot_MULTI_lon_transect --
        %  - set icalc_diurnal to =1 and set time_of_day_multi =
        % {'daytime','nighttime'}
        
        ilat_use=1; %index of the latitude to use; 1=30S, 2=20S, 3=10S
        
        var_plot_label = var_plot;
       
       switch var_plot
           case 'low_CF'
               unit_str = '';
               istart=1; 
           case 'Precip_rate'
               unit_str = '(mm hr^{-1})';
               istart=2; %don't bin the precip data as already at the desired 4 deg spacing
           case {'LWP','LWP2','TWP','TLWP'}
               istart=1; %bin the AMSRE LWP data
               
               switch var_plot
                   case 'TWP'
                       var_plot_label = 'TWP';  %the y-label
                   case 'TLWP'
                       var_plot_label = 'TLWP';  %the y-label                       
                   otherwise
                       var_plot_label = 'LWP';  %the y-label
               end
               unit_str = '(g m^{-2})';
               
%                %de Szoeke values, 75-85W, 20S
%                add_points=1;
%                %night (01:30am AMSRE equivalent()
%                xpos(1).x = -80;
%                ypos(1).y = 130;
%                point_labs(1).lab = 'Ship (night)';
%                
%                xpos(2).x = -80;
%                ypos(2).y = 50;
%                point_labs(2).lab = 'Ship (day)';
%                
                %de Szoeke values, 75-85W, 20S
%               add_points=1;
               %night (01:30am AMSRE equivalent()
               xbox(1).x = [-85 -75];
               ybox(1).y = [130 50];
               colbox(1).c = 'k';
               box_type = 'rectangle';
%               point_labs(1).lab = 'Ship (night)';
%               
%                xpos(2).x = -80;
%                ypos(2).y = 50;
%                point_labs(2).lab = 'Ship (day)';

                  script_name='draw_box';
                  if  strcmp(lon_load_multi{ilat_use},'20S') %only if is the 20S plot
                      iexecute_script=1;
                  end
               

               
%               add_points_to_plot(xpos,ypos,point_labs,8,11)

                 highlight_type = 'fill';


               
               
       end
        

       ylab=[remove_character(var_plot_label,'_',' ') ' ' unit_str];
       xlab= 'Longitude';
       

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        
        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.
         lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        iovr_leg_line=1; %flag which plots a '-' line in the legend even if linestyle has been set to 'none'



       
        idat=0;
        
       
        

                                                    
        
        
        lon_bin_edges=[-104   -100   -96   -92   -88   -84   -80   -76  -72];
        lon_bin_mid = 0.5*(lon_bin_edges(1:end-1) + lon_bin_edges(2:end));
        
        values_binned = value_array;
        
        
        
        for idat=istart:length(value_array{1,1,1}) %all lines/models/obs etc.
            for iday=1:2
                lon = value_array{ilat_use,1,iday}(idat).x;
                dat = value_array{ilat_use,1,iday}(idat).y;
                dlon = meanNoNan(diff(lon),2);
                lon_edges = [lon-dlon/2 lon(end)+dlon/2];


                values_binned{ilat_use,1,iday}(idat).y=[];
                for ilon=1:length(lon_bin_edges)-1
                    ibin = find(lon_edges>=lon_bin_edges(ilon) & lon_edges<lon_bin_edges(ilon+1));
                    if length(ibin)>0
                        values_binned{ilat_use,1,iday}(idat).y(ilon) = meanNoNaN(dat(ibin),2);
                    else
                         values_binned{ilat_use,1,iday}(idat).y(ilon) = NaN;
                    end
                end
                
                values_binned{ilat_use,1,iday}(idat).x = lon_bin_mid;

            end



        end
        
         ierror_bars='horiz2';
         ierror_bars='vert2';
         ierror_bars='vert2 filled';         
         
         ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
         istyle=1;
        
         for idat=1:length(values_binned{1,1,1}) %all lines/models/obs etc.


             xdat(idat).x = values_binned{ilat_use,1,1}(idat).x;
             day = values_binned{ilat_use,1,1}(idat).y;
             night = values_binned{ilat_use,1,2}(idat).y;
             ydat(idat).y=0.5*(day+night); %

% use error bars to represent the diurnal range of the data
             errordatU(idat).dat = night - ydat(idat).y;
             errordatL(idat).dat = ydat(idat).y - day;


             labs(idat).l=values_binned{ilat_use,1,1}(idat).lab;
             
             
             switch labs(idat).l
                 case {'AMSRE','CloudSat','MODIS  ','CALIPSO monthly nighttime ','CALIPSO monthly DAYTIME ','CALIPSO monthly average ',' CLOUDSAT PRECIP  ','ascending CLOUDSAT PRECIP  ','descending CLOUDSAT PRECIP  ','AMSRE daytime ','AMSRE nighttime ','AMSRE average '}
                     %line_pattern(istyle).p= '-.';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;                    
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     
                     switch labs(idat).l
                         case {'CALIPSO monthly nighttime ','CALIPSO monthly DAYTIME '}
                             labs(idat).l = 'CALIPSO';                             
                     end
                 case {'CAM5 1deg','CAM5  CAM5 1deg','CAM5 COSP  CAM5 1deg'}
%                     labs(idat2).l = 'CAM5 1deg';
%                     line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                 case {'CAM5 CLUBBv2'}
%                     line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0 0 0.6]; marker_style(istyle).m='d'; istyle=istyle+1;
%                     labs(idat2).l = 'CAM5 CLUBB';
                 case {'CAM5 CLUBB','CAM5 CLUBB  CAMCLUBB','CAM5 CLUBB COSP  '}
%                     line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
                      line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
%                     labs(idat2).l = 'CAM5 CLUBB';
                 case {'CAM5 2deg','CAM5  CAM5 2deg'}
%                     line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.6 0.6 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0.6 0.6 0]; marker_style(istyle).m='d'; istyle=istyle+1;
%                     labs(idat2).l = 'CAM5 2deg';
                 case {'AM3 2deg','AM3  2deg'}
%                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;                     
%                     labs(idat2).l = 'AM3 2deg';
                 case {'AM3  0pt5deg','AM3 0.5deg'}
%                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;                     
%                     labs(idat2).l = 'AM3 0.5deg';
                 case {'AM3 CLUBB  AM3CLUBB','AM3 CLUBB'}
%                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
%                     labs(idat2).l = 'AM3 CLUBB';
                 case {'AM3 CLUBBv2'}
%                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= 'none';  line_colour(istyle).c=[0.4 0.4 0.0]; marker_style(istyle).m='d'; istyle=istyle+1;
%                     labs(idat2).l = 'AM3 CLUBB';

             end


             
             if idat==2 | idat==3
                 ibar=(idat-4); %i.e. -1 or -2
             elseif idat==4 | idat==5
                 ibar = idat-3; %= 1 or 2
             end
             
             if idat>1
%                 ydat(idat).y=ydat(1).y * (1+(idat-1)/3);
                 xdat(idat).x= lon_bin_mid + ibar*0.4;
%                 errordatU(idat).dat = ydat(idat).y*0.2;
%                 errordatL(idat).dat = ydat(idat).y*0.2;
             end
             

         end
        

         titlenam = [var_plot_label ' Longitude Transect at ' [lon_load_multi{ilat_use}] ' for ' time_str_116_2];
         figname=titlenam;
         savename=figname;
        
        
        
        
    case 134
        %plot simple 1D PDFS of tau and Re

        pdf_type = 'screened';
%        pdf_type = 'simple';

        %xlab = ['N_d (cm^{-3}) for ' low_or_high_CF];
        %xlab = ['N_d (cm^{-3})'];
        xlab = ['Homogeneity Parameter timeseries3 using mean W for each pixel'];
        xlab = ['Homogeneity Parameter timeseries3 using mean W for each pixel - screened'];
        xlab = ['Cloud Fraction'];
        %        xlab = ['LWP'];

        %In these scripts plotTimeHeightVap3 is just run in order to do the
        %screening and make the X and Y values with the correct changing to deal
        %with ndhistc. We just use the X and Y values so the bin sizes don't really
        %matter. Can just select to get the desired Xbins or Ybins for the final
        %PDF
        
        ichoose_Ybins=0; %default

        switch pdf_type
            case 'screened'



                thresh_sensZA=[0 90];

                thresh_CF = [0.1 1.0];
                thresh_CF = [0.1 0.8];
                thresh_CF = [0.8 1.0];
                %                thresh_CF = [0.85 1.0];
                %                thresh_CF = [0.6 1.0];
                %                thresh_CF = [0.0 1.0];



                thresh_relAZ = [50 130];
                %                thresh_relAZ = [60 120];
                thresh_relAZ = [0 180];
                thresh_SZA = [50 55];
                %                thresh_SZA = [75 85];
                %                thresh_SZA = [70 85];
                %                thresh_SZA = [0 85];
                thresh_CTT = [273-100 273-5];
                thresh_CTT = [273-5 373];

                %changed this to the homogeneity factor now
                stdW_boundary=1;
                thresh_stdW = [stdW_boundary 1e9];
                thresh_stdW = [0 stdW_boundary];
                thresh_stdW = [0 1e9];


                switch xlab
                    case 'Homogeneity Parameter timeseries3 using mean W for each pixel - screened'
                        

                        y_axis_vals = 'Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins = [0:0.5:50]; ichoose_Ybins=1;
                        %x_axis_vals = 'Homogeneity Parameter timeseries3 using Cahalan log mean W
                        %(pixel level)'; 
                        x_axis_vals=y_axis_vals;
                         
                        case134_set = 'low/high SZA and CF';
                        
                    case 'Cloud Fraction'
                        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins=[0:0.05:1.0 1.00001]; ichoose_Ybins=1;
                        %        y_axis_vals='Cloud Fraction from grid vals timeseries3';
                        %        Ybins=[0.8:0.02:1.0 1.00001]; ichoose_Ybins_set=1;
                        
                        x_axis_vals=y_axis_vals;
                        
                        case134_set = 'allCF low/high SZA';


                end
                
                  %------------------
                   case134_pdf_calcs
                  %------------------




                    case 'simple'

                        switch xlab
                            case 'Homogeneity Parameter timeseries3 using mean W for each pixel'


                        ierror_bars='none';

                        istyle = 1;
                        Npdf_pts = 50;

                        idat2=0;
                        idat2=idat2+1;
                        %dat=homog_time3_W;
                        dat=homog_time3_meanW;
                        dat(dat>1000)=1000;
                        dat=dat(:);
                        %        Xbins=make_PDF_bins(dat,Npdf_pts);
                        Xbins = [0:0.1:20];
                        ydat(idat2).y = ndHistc_run([dat],Xbins);
                        ydat(idat2).y = ydat(idat2).y / sum(ydat(idat2).y);
                        xdat(idat2).x = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );
                        labs(idat2).l = '';

                    case 'LWP'
                        ierror_bars='none';

                        istyle = 1;
                        Npdf_pts = 50;

                        idat2=0;
                        idat2=idat2+1;
                        dat=Cloud_Water_Path_Liquid.timeseries3;
                        %                dat=W_time3;
                        %        dat=homog_time3_meanW(ilat,ilon,itime);
                        %                dat(dat>1000)=1000;

                        Xbins=make_PDF_bins(dat,Npdf_pts);
                        %               Xbins = [0:0.1:20];
                        dat=dat(:);
                        ydat(idat2).y = ndHistc_run([dat],Xbins);
                        ydat(idat2).y = ydat(idat2).y / sum(ydat(idat2).y);
                        xdat(idat2).x = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );
                        labs(idat2).l = '';

                end



        end



        noplot=0;
        ylab = 'Frequency';
        plot_error_bars = 0;


    case 133
        mean_or_pdf = 'mean';
%        mean_or_pdf = 'pdf';
        
        thresh_NP=10;
        N_sza_bins=25;
%      N_sza_bins=12;  

time_mean_str = 'ALL'; %choose time period

%choose type of screening
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';
screen_type ='none';

%set defaults of these, which may get overwritten below
ichoose_Xbins=0; %flags to allow the direct specification of Xbins, Ybins, Zbins
ichoose_Ybins=0;
ichoose_Zbins=0;



        y_axis_vals='Mean SZA timeseries3';
        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0:0.1:1.0 1.00001]; ichoose_Ybins_set=1;
        % %        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0.8:0.02:1.0 1.00001]; ichoose_Ybins_set=1;
        %         y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:10]; ichoose_Ybins_set=1;
        %         y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:25]; ichoose_Ybins_set=1;
        %         y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:2:25]; ichoose_Ybins_set=1;
        % %        Ybins_set=[0:0.1:0.9 0.99 1.00001]; ichoose_Ybins_set=1;
        %         y_axis_vals='Mean CTT timeseries3, y-axis';  ichoose_Ybins_set=0;  %Ybins_set=[0:2:25];
        %         y_axis_vals='Mean CTH timeseries3, y-axis';  ichoose_Ybins_set=0;  %Ybins_set=[0:2:25];
        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0:0.1:1.0 1.00001]; ichoose_Ybins_set=1;
        %CTH using CTT from all available clouds (liquid, ice and undetermined)
        y_axis_vals = 'Cloud Top Height all clouds'; ichoose_Ybins_set=1; Ybins_set=[-2:0.2:10];

        %x_axis_vals = 'Homogeneity Parameter timeseries3 using W from mean tau and Re';
        x_axis_vals = 'Homogeneity Parameter timeseries3 using mean W for each pixel';
        %x_axis_vals = 'Homogeneity Parameter timeseries3 using Cahalan log mean W
        %(pixel level)';
        x_axis_vals = 'CF all'; %ice + liquid + undetermined cloud fraction
        
        
        ichoose_Ybins=ichoose_Ybins_set; Ybins=Ybins_set;

        man_choose_plotTimeHeight_graph=1;
        i577 = 'MODIS_plot_UW';
        logflag=0;
        dlogflag=0;
        noplot=0; %flag to say to just do the calculations and not to plot
        nXpdf=500;
        nYpdf=N_sza_bins;
        ndims_hist=2;
        ioverride_time_selection=1;
        ioverride_pdf_varchoose=1; %override the defaults
        ioverride_pdf=1; %override the defaults
        plotTimeHeightVap3

        man_choose_water_graph=1;
         switch mean_or_pdf
                case 'pdf'
                    graph=977; %1d pdf
                case 'mean'
                    graph=96; %mean of 2D PDF
         end
        noplot=0;
        waterVapourMay2005

        noplot=1; %flag to say to just do the calculations and not to plot





        
        
        
case 132  %CALIPSO CFAD CTH pdfs - run read_cfads_calipso_monthly_night_ipsl first
    
    LAT_val = [25 34]; LON_val = [122:129]; %China Sea
    LAT_val = [44:51]; LON_val = [-30:-6]; %SW UK/Atlantic - ocean only
    LAT_val = [40:50]; LON_val = [30:60]; %Central Europe    
    LAT_val = [-43 -40]; LON_val = [110 144]; %west of Tasmania (similar to Boers, QJRMS, 1998)
    LAT_val = [72 75]; LON_val = [-3:48]; %Arctic summer box region
    LAT_val = [-70 -60]; LON_val = [-20:0]; %Seascape Neumayer region  
    LAT_val = [-50 20]; LON_val = [-160 -60]; %VOCALS
    LAT_val = [-40.5 10.5]; LON_val = [-140 -50]; %VOCALS CAPT (whole map for CPT paper plots)
%    LAT_val = [-40.5 -30.5]; LON_val = [-140 -100]; %SW region where there is high cloud.
%    LAT_val = [-25.5 -15.5]; LON_val = [-80 -70];     
    
    %Select the months to include.
    %Months start Jan 2007 and run for 4 years (48 in total)
    itimes_cfad={[1 2],[6 7 8]}; timstr_cfad={'DJF','JJA'};
    itimes_cfad={[6 18 30 42]}; timstr_cfad={'June'}; %2007-2010
    itimes_cfad={[1 2 12 13 14 24 25 26 36 37 38 48]}; timstr_cfad={'DJF'}; %2007-2010 DJF (for SEASCAPE)    
    %[1 2]; %Jan Feb 2007
    %itime=[6 7 8]; %JJA 2007
    itimes_cfad={[1:48]}; timstr_cfad={'ALL'};
        
    normalise_by_bins = 0;



    for icfad=1:length(itimes_cfad)
        itime = itimes_cfad{icfad};
        plotCFAD_1D_dist_for_region
        
        if normalise_by_bins==1
            bin_widths = diff(cfad_alts_edges_CALIPSO/1e3);
            norm_str='bin normalised';
            ylab='Bin normalised frequency';
        else
            bin_widths=1;
            norm_str='';   
            ylab='Cloud Fraction';
        end
        
        N = sum(cf_cfad);
        N = 1;
        
        xdat(icfad).x = cfad_alts_centre_CALIPSO/1e3;
        ydat(icfad).y = cf_cfad./bin_widths/N; 
        labs(icfad).l=timstr_cfad{icfad};
    end
    
    xcal=xdat;
    ycal=ydat;
        
        
    set(gca,'xlim',[-2 10]);
    grid on
    title(['time=' num2str(itime)]);
    
     titlenam = [data_type_cf_cal ' ' norm_str ' frequency for LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(end)) ' LON=' num2str(LON_val(1)) ' to ' num2str(LON_val(end))];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[-2 10];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.


        xlab= 'Cloud Top Height (km)';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
       istyle=1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
        
        
        
case 131  %run multiple runs of PlotTimeHeightVap3 in order to get N lines for the means
           %of 2d pdfs (e.g. Nd vs solarZA, etc)
           
      clear xvars errorvars
      
      %Some later additions (18th April, 2016) - trying to get this working
      %with MODIS load using load_saved_modis_vars for a regional dataset
      %(restricted region, not global).
      
      datatype='timeseries3'; LAT=MLAT; LON=MLON;
                     
      i_plot_norm=1; %flag to say whether we want to plot a normalised PDF or a PDF of N
      axis1D = 'x';
      axis1D = 'y';
      
      logbin_norm=0;
      
      iadd_data=1; %flag to add another PDF to that calculated from MODIS (e.g. CALIPSO pdf)
      
      ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
      
%      ivary_thresh=1; %flag to say that we want a different thresh_CF etc for each line
       param_vary = 'time';
       param_vary = 'thresh';
       

      plot_error_bars = 0; %is reset below in certain cases
      thresh_NP=10;
%      thresh_NP=50;
      
      N_sza_bins=25;
%      N_sza_bins=12;      
      
      ichoose_Xbins_set=0; Xbins_set=NaN;
      ichoose_Ybins_set=0; Ybins_set=NaN;
      ichoose_Zbins_set=0; Zbins_set=NaN;
      
      mean_or_pdf = 'pdf';
%      mean_or_pdf = 'mean';   
      mean_or_pdf = 'cf_array';  

      minfrac_CF = 0.9; %minimum fraction of the sampled points that had successful cloudy/clear/phase
      %determination (i.e. Npix/Nptot_mockL3 =
      %Cloud_Fraction_Liquid_Pixel_Counts./Cloud_Fraction_Liquid./Total_pix
      %els

      minfrac_NpNd = 0.9; 
      %Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3./Cloud_Fraction_Liquid_Pixel_Counts.timeseries3
       %Fraction of points that remain after all previous filtering for
       %which we have an Nd retrieval. Restriction (4) in the SZA paper.
      
      thresh_NP_Nd = 10; %min no. of pixels required for an Nd, re, tau, etc measurement to count
      %(Uses Cloud_Fraction_Liquid_Pixel_Counts2). For this screening usually
      %thresh_NP is the number of pixels that the swath must have covered -
      %i.e. the total number of pixels available. Only a portion of those 
      

switch mean_or_pdf
    case 'pdf'
        nXpdf_131 = 150;        
    case 'mean'
        nXpdf_131 = 500;
end

%Have saved some of the PDFs here:-
save_pdfs_filename='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs.mat';
%see save_TauRe_attribution_PDF_data for the script to add the data to the
%arrays and save
      
      if ~exist('ioverride_122') | ioverride_122==0

        multi_plot_case_122 = 'Nd vs CTH';
%        multi_plot_case_122 = 'Nd vs CTH_no_subtract';        
%        multi_plot_case_122 = 'Nd vs SST-CTT';
%        multi_plot_case_122 = 'Nd vs SST';
%        multi_plot_case_122 = 'Nd vs CTT';
%        multi_plot_case_122 = 'Nd vs minCTT';        
        multi_plot_case_122 = 'Nd vs CTT all';
%        multi_plot_case_122 = 'Nd vs homog';
%        multi_plot_case_122 = 'Nd vs CTT std_dev';      
%        multi_plot_case_122 = 'Nd vs tau-homog Cahalan';
%        multi_plot_case_122 = 'CF vs CTH';     
%        multi_plot_case_122 = 'CTH vs CTH'; %for CTH pdf
%        multi_plot_case_122 = 'Re vs CTT'; 
%        multi_plot_case_122 = 'Tau vs CTT';
%        multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
%         multi_plot_case_122 = 'Homog vs SZA';
%         multi_plot_case_122 = 'Nd vs SZA';   
%         multi_plot_case_122 = 'CF_MOD35 vs CTH_Ryan MOD35'; %attempting to replicate the time-averaged CF that CALIPSO will
         multi_plot_case_122 = 'CF_MOD35 vs CTH MOD35'; %attempting to replicate the time-averaged CF that CALIPSO will
           %report - when there is no cloud MODIS CTT is NaN and so there
           %is no CTH - but zeros don't get reported by this method. So,
           %making a profile of zero CFs for each height
%         multi_plot_case_122 = 'CF_MOD35 vs CTH_max MOD35'; %using the min CTT and so max CTH
%         multi_plot_case_122 = 'CF_MOD35 vs CTH_max_CTP_hybrid MOD35'; %using the min CTT and so max CTH, but also using 
                                                                       %the min CTP and converting to CTH using std. atmos (see filtering_data_get) 

%         multi_plot_case_122 = 'CF_all vs CTH all'; %plotting CF vs CTH for MODIS is biased since
         %we need a cloud in the scene to get a CTH and then the whole
         %scene is likely to be cloudy. would be good to plot tau vs CTH,
         %but need tau for all cloud phases.
%         multi_plot_case_122 = 'CF liq vs CF liq';
%         multi_plot_case_122 = 'CF ice vs CF ice';         
%         multi_plot_case_122 = 'CF undet vs CF undet';  
%         multi_plot_case_122 = 'CTT all vs CTT all';  
%         multi_plot_case_122 = 'CTH all vs CTH all';           
%         multi_plot_case_122 = 'Tau vs Tau';  
%         multi_plot_case_122 = 'Re vs Re';   
%         multi_plot_case_122 = 'Nd vs Nd';
%         multi_plot_case_122 = 'Nd_factor vs Nd_factor'; %the
                                        %multiplication factor to go from tau^0.5*Re^(-2.5) to N
%         multi_plot_case_122 = 'Percent diff Nd_all and Nd_lowSZA';

      end
              
      switch multi_plot_case_122
          case {'Tau vs Tau','Re vs Re','Nd vs Nd','Nd_factor vs Nd_factor','Percent diff Nd_all and Nd_lowSZA'}
              iswap_xy=0;
              y_axis_vals = 'Dummy data for 1D';     
              axis1D = 'x';        
          case {'CTH all vs CTH all'}
              iswap_xy=1; %flag to swap the defined x and y
              y_axis_vals = 'Cloud Top Height all clouds'; ichoose_Ybins_set=1; Ybins_set=[-2:0.2:10];          
           case {'CTT all vs CTT all'}
              iswap_xy=1; %flag to swap the defined x and y
              y_axis_vals='CTT all clouds';  ichoose_Ybins_set=1;  Ybins_set=[230:1:283];
          case {'Nd vs CTT all'}
              y_axis_vals='CTT all clouds';  ichoose_Ybins_set=1;  Ybins_set=[230:1:283];
          case {'CF_all vs CTH all'}
              y_axis_vals = 'Cloud Top Height all clouds'; ichoose_Ybins_set=1; Ybins_set=[-2:0.2:10];
          case {'CF_MOD35 vs CTH MOD35'}
              y_axis_vals = 'Mean CTH timeseries3 flag zeroCF, y-axis'; ichoose_Ybins_set=1; Ybins_set=[-0.48*6:0.48:20]; %Ybins_set=[-2:0.2:20]; 
              param_vary=''; %Prevent the legend labelling.
          case {'CF_MOD35 vs CTH_Ryan MOD35'}
              y_axis_vals = 'Mean CTH_Ryan timeseries3 flag zeroCF, y-axis'; ichoose_Ybins_set=1; Ybins_set=[-0.48*6:0.48:20]; %Ybins_set=[-2:0.2:20]; 
              param_vary=''; %Prevent the legend labelling.              
          case {'CF_MOD35 vs CTH_max MOD35'}
              y_axis_vals = 'Mean CTH_max timeseries3 flag zeroCF, y-axis'; ichoose_Ybins_set=1; Ybins_set=[-0.48*6:0.48:20]; %Ybins_set=[-2:0.2:20]; 
              param_vary=''; %Prevent the legend labelling.  
          case 'CF_MOD35 vs CTH_max_CTP_hybrid MOD35'  
              y_axis_vals = 'Mean CTH_max_CTP_min timeseries3 flag zeroCF, y-axis'; ichoose_Ybins_set=1; Ybins_set=[-0.48*6:0.48:20]; %Ybins_set=[-2:0.2:20];
              param_vary=''; %Prevent the legend labelling.
          case {'Nd vs CTH','CF vs CTH'}
%              y_axis_vals='Mean CTH timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[-2:0.2:10];
              y_axis_vals='Mean CTH timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[-2:0.5:10];              
          case {'Nd vs CTH_no_subtract'}
              y_axis_vals='Mean CTH_no_subtract timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[-2:0.5:10];              
          case {'Nd vs CTT','Re vs CTT','Tau vs CTT'}
              y_axis_vals='Mean CTT timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[230:1:283];
          case {'Nd vs minCTT'}
              y_axis_vals='Min CTT timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[250:1:283];              
          case {'Nd vs SST-CTT'}
                 y_axis_vals= 'SST - Mean CTT timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[-2:1.5:10];              
          case {'Nd vs SST'}
                 y_axis_vals= 'SST timeseries3, y-axis';  ichoose_Ybins_set=1;  Ybins_set=[-2:1:15];                               
          case {'Nd vs homog'}
              y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:10]; ichoose_Ybins_set=1;
              y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:25]; ichoose_Ybins_set=1;
              y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:2:25]; ichoose_Ybins_set=1;    
              y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:0.5:25 35 45 55 65 100]; ichoose_Ybins_set=1;               
              y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:2:25 35 45 55 65 100]; ichoose_Ybins_set=1;                             
              y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:0.5:100]; ichoose_Ybins_set=1;                   
          case {'Nd vs tau-homog Cahalan'}
              y_axis_vals='Homogeneity Parameter Cahalan Optical Depth (Seethala)'; Ybins_set=[0:0.01:1]; ichoose_Ybins_set=1;              
          case {'Homog vs SZA','Nd vs SZA'}
               y_axis_vals='Mean SZA timeseries3'; Ybins_set=[0:5:90]; ichoose_Ybins_set=1;
          case {'Nd vs CTT std_dev'}
              y_axis_vals='Cloud Top Temp standard deviation, liquid pixels'; Ybins_set=[0:0.1:25]; ichoose_Ybins_set=1;
%              y_axis_vals='Cloud Top Temp standard deviation, liquid pixels'; Ybins_set=[0:0.2:2 2.2:1:25]; ichoose_Ybins_set=1;              
          case {'CF liq vs CF liq'}
              y_axis_vals = 'Cloud Fraction from grid vals timeseries3'; ichoose_Ybins_set=1; Ybins_set=[0:0.01:1];
          case {'CF ice vs CF ice'}  
              iswap_xy=1; %flag to swap the defined x and y
              y_axis_vals = 'CF ice'; ichoose_Ybins_set=1; Ybins_set=[0:0.01:1];              
          case {'CF undet vs CF undet'}  
              iswap_xy=1; %flag to swap the defined x and y
              y_axis_vals = 'CF undet'; ichoose_Ybins_set=1; Ybins_set=[0:0.01:1];                
          case {'CF undet vs CF undet'}  
              iswap_xy=1; %flag to swap the defined x and y
              y_axis_vals = 'CF undet'; ichoose_Ybins_set=1; Ybins_set=[0:0.01:1];                              
          otherwise
              y_axis_vals='Mean SZA timeseries3';
              y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0:0.1:1.0 1.00001]; ichoose_Ybins_set=1;
              %        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0.8:0.02:1.0 1.00001]; ichoose_Ybins_set=1;
          
              %        Ybins_set=[0:0.1:0.9 0.99 1.00001]; ichoose_Ybins_set=1;
      end

      

        
        switch mean_or_pdf
            case 'mean'
                plot_error_bars = 1;
            case 'pdf'
                plot_error_bars = 0;
        end
    
         

    
                                 
                
%                screen_strs={'NP + CF mockL3','NP + CF + MEAN sensZA'};
%                screen_strs={'NP + CF + MEAN sensZA','NP + CF + MEAN sensZA','NP + CF + MEAN sensZA'};                
%                screen_strs={'NP + CF + MEAN sensZA','NP + CF + MEAN sensZA'};
%                screen_strs={'NP + CF + MEAN sensZA + MEAN relAZ','NP + CF + MEAN sensZA + MEAN relAZ'};                
            sens_limit = '41.4';
            sens_limit = '90';  
            


 if ~exist('ioverride_122') | ioverride_122==0
                

                thresh_cf_122 = {[0.0 1.0],[0.0 1.0]}; 
%                thresh_cf_122 = {[0.8 1.0],[0.1 0.8]}; 
%                thresh_cf_122 = {[0.8 1.0],[0.8 1.0]};
%                thresh_cf_122 = {[0.8 1.0],[0.8 1.0],[0.8 1.0]};                
%                thresh_cf_122 = {[0.99 1.0],[0.99 1.0]};  
                thresh_cf_122 = {[0.9999 1.0],[0.9999 1.0]};             
%                thresh_cf_122 = {[0.9999 1.0],[0.9999 1.0]}; thresh_cf_122=[thresh_cf_122 thresh_cf_122 thresh_cf_122];
%                thresh_cf_122 = {[0.1 0.8],[0.1 0.8]};                     
%                thresh_cf_122 = {[0.99 1.0],[0.99 1.0],[0.99 1.0],[0.99 1.0],[0.99 1.0],[0.99 1.0]};                       
%                thresh_cf_122 = {[-0.01 1.0],[-0.01 1.0]};                



 end
                
                thresh_sensZA_122 = {[41.4 90],[41.4 90]};
                thresh_sensZA_122 = {[0 41.4],[0 41.4]};                
             
                thresh_AZ_122 = [50 130];
%                thresh_AZ_122 = [60 120];    
                thresh_AZ_122 = {[0 180],[0 180],[0 180]}; 
                thresh_AZ_122 = {[0 180],[0 180],[0 180]}; thresh_AZ_122=[thresh_AZ_122 thresh_AZ_122 thresh_AZ_122];  
%                thresh_AZ_122 = {[0 180],[0 180],[0 180],[0 180],[0 180],[0 180]};   
                
                thresh_solarZA_122 = [50 55];
%                thresh_solarZA_122 = [75 85];  
%                thresh_solarZA_122 = [70 85];                  
%                thresh_solarZA_122 = [0 85];      
                thresh_solarZA_122 = {[50 55],[75 85]};  
                thresh_solarZA_122 = {[50 55],[75 85]}; thresh_solarZA_122=[thresh_solarZA_122 thresh_solarZA_122 thresh_solarZA_122];  
%                thresh_AZ_122 = {[50 55],[50 55]};   
%                thresh_solarZA_122 = {[75 85],[75 85]};   
%                thresh_solarZA_122 = {[0 85],[0 85]};    
%                thresh_solarZA_122 = {[0 65],[0 65]};    
%                thresh_solarZA_122 = {[0 72],[0 72]};                    
%                thresh_solarZA_122 = {[0 85],[0 85],[0 85]};    
%                thresh_solarZA_122 = {[0 75],[75 85]};                    
%                thresh_solarZA_122 = {[0 65],[65 85]};                                    
                
%                thresh_solarZA_122 = {[0 30],[50 65]};  
%                thresh_solarZA_122 = {[0 65],[75 85]};     
%                thresh_solarZA_122 = {[50 55],[75 85]};     
%                thresh_solarZA_122 = {[50 55],[55 60],[60 65],[65 70],[70 75],[75 85]};     
                
%                thresh_CTT_122 = {[273-100 373],[273-100 373],[273-100 373],[273-100 373],[273-100 373],[273-100 373]};   
                thresh_CTT_122 = {[273-100 373],[273-100 373]};  
%                thresh_CTT_122 = {[260 268],[260 268]};                  
%                thresh_CTT_122 = {[273-100 268],[268 373]};         
%                thresh_CTT_122 = {[273-100 253],[268 373]};    
%                thresh_CTT_122 = {[263 373],[263 373]};   
%                thresh_CTT_122 = {[265 373],[265 373],[265 373]};   
%                thresh_CTT_122 = [273-5 373];    
                thresh_CTT_122 = {[273-5 373],[273-5 373]};   
                thresh_CTT_122 = {[273-5 373],[273-5 373]}; thresh_CTT_122=[thresh_CTT_122 thresh_CTT_122 thresh_CTT_122];  
%                thresh_CTT_122 = {[273 373],[273 373]};                   
%                thresh_CTT_122 = [273-100 273-3];  

                 thresh_CTH_122 = {[-2 20],[-2 20],[-2 20]}; %km
                 thresh_CTH_122 = {[-2 20],[-2 20],[-2 20]}; thresh_CTH_122=[thresh_CTH_122 thresh_CTH_122 thresh_CTH_122];                   
%                 thresh_CTH_122 = {[-2 2],[-2 2],[-2 2],[-2 2],[-2 2],[-2 2]}; %km
%                 thresh_CTH_122 = {[-2 2],[-2 2]}; %km

                thresh_sigCTT_122 = {[0 0.65],[0 0.65],[0 0.65]}; thresh_sigCTT_122=[thresh_sigCTT_122 thresh_sigCTT_122 thresh_sigCTT_122];
                thresh_sigCTT_122 = {[1 1e9],[1 1e9],[1 1e9]}; thresh_sigCTT_122=[thresh_sigCTT_122 thresh_sigCTT_122 thresh_sigCTT_122];
                thresh_sigCTT_122 = {[0 1e9],[0 1e9],[0 1e9]}; thresh_sigCTT_122=[thresh_sigCTT_122 thresh_sigCTT_122 thresh_sigCTT_122];                
%                 thresh_CTH = [-2 2];
                 
                %changed this to the homogeneity factor now
                stdW_boundary=1;
                thresh_stdW_122 = [stdW_boundary 1e9];                
                thresh_stdW_122 = [0 stdW_boundary];
                thresh_stdW_122 = {[0 1e9],[0 1e9],[0 1e9]}; 
                thresh_stdW_122 = {[0 1e9],[0 1e9],[0 1e9]}; thresh_stdW_122=[thresh_stdW_122 thresh_stdW_122 thresh_stdW_122];  
%                thresh_stdW_122 = {[0 1e9],[0 1e9],[0 1e9],[0 1e9],[0 1e9],[0 1e9]};  

thresh_homog_tau = [0 0.7];
thresh_homog_tau = [0 1.1];    

   clear line_pattern line_colour
       istyle=1;
%       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
%       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;

       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;


%        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='s'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='s'; istyle=istyle+1;
%        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%         

                

               
             
%                 screen_strs={...
%                     'none',...
%                     'none',...
%                     };
% 
%                 screen_strs={...
%                     'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
%                     'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
%                     'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
%                     };
                
%                screen_strs={...
%                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
%                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
%                    };
  

                                                               
        
        switch multi_plot_case_122
              case {'Homog vs SZA'}                                      
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


                    
                 xvars={...
                      'Homogeneity Parameter timeseries3 using mean W for each pixel'...
                     ,'Homogeneity Parameter timeseries3 using mean W for each pixel'...                    
                     };
                 
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
                     };
                 
           case {'Tau vs CTT','Tau vs Tau'}                                      
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


                 for ilines=1:length(thresh_cf_122)
                     xvars{ilines}='Tau reduced dataset Re_1.6 Re_3.7';

                 end
                 
                 for ilines=1:length(thresh_cf_122)
                     errorvars{ilines} = 'Optical Depth uncertainty from grid vals timeseries3';                     
                 end  
                 
                 ichoose_Xbins_set=1;  Xbins_set=[0:0.1:2.3]; Xbins_set=10.^(Xbins_set);
                 logbin_norm=1;
                 x_axis_type = 'log10_matlab'; %make xscale = 'log'
                 
            case 'Percent diff Nd_all and Nd_lowSZA'
                                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


xvars = {'Percentage Nd difference 1.6um allSZA vs lowSZA'...
    'Percentage Nd difference 2.1um allSZA vs lowSZA'...
    'Percentage Nd difference 3.7um allSZA vs lowSZA'...
    };
                 
                 for ilines=1:length(thresh_cf_122)
                     errorvars{ilines} = 'Optical Depth uncertainty from grid vals timeseries3';                     
                 end  
                 
%                 ichoose_Xbins_set=1;  Xbins_set=[0:0.1:2.3]; Xbins_set=10.^(Xbins_set);
%                 logbin_norm=1;
%                 x_axis_type = 'log10_matlab'; %make xscale = 'log'
                 
                 
        Xbins_set = [-70:4:70]; ichoose_Xbins_set=1;

           case {'Re vs CTT','Re vs Re'}                                      
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


                    
%                  xvars={...
%                       'Nd from grid vals timeseries3'...
%                      ,'Nd from grid vals timeseries3'...                    
%                      };
%                  
%                  errorvars={...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
%                      };


for ilines=1:length(thresh_cf_122)
    if ~exist('ioverride_122_save')
        xvars{ilines}='R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%        xvars{ilines}='R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%        xvars{ilines}='R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
    else
        xvars{ilines}=xvar_save_str;
    end
end

                 
%                  for ilines=1:2
%                      xvars{ilines}='R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%                  end
%                  for ilines=3:4
%                      xvars{ilines}='R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';                     
%                  end
%                  for ilines=5:6
%                      xvars{ilines}='R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%                  end
                 
                 
                 for ilines=1:length(thresh_cf_122)
%                     errorvars{ilines} = 'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                     
                     errorvars{ilines} = 'Re 2.1 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                     errorvars{ilines} = 'Re 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)';                     
                 end  
                 
                 ichoose_Xbins_set=1;  Xbins_set=[1:0.75:20];
                 
           case {'Nd vs Nd'}                                      
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


                    
%                  xvars={...
%                       'Nd from grid vals timeseries3'...
%                      ,'Nd from grid vals timeseries3'...                    
%                      };
%                  
%                  errorvars={...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
%                      };

                 for ilines=1:length(thresh_cf_122)
                     if ~exist('ioverride_122_save')
                     xvars{ilines}= 'Nd_{1.6} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7';                     
                     xvars{ilines}='Nd from grid vals timeseries3';
                     xvars{ilines}= 'Nd_{3.7} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7';
                     else
                         xvars{ilines} = xvar_save_str;
                     end
                 end
                 
                 for ilines=1:length(thresh_cf_122)
                     errorvars{ilines} = 'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                 end 
             
             ichoose_Xbins_set=1;  Xbins_set=[0:5:500];
             
             choose_Xbins_set=1;  Xbins_set=[0:0.1:3]; Xbins_set=10.^(Xbins_set);
             logbin_norm=1;
             
             x_axis_type = 'log10_matlab'; %make xscale = 'log'
             
             
            case 'Nd_factor vs Nd_factor'
                ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


                    
%                  xvars={...
%                       'Nd from grid vals timeseries3'...
%                      ,'Nd from grid vals timeseries3'...                    
%                      };
%                  
%                  errorvars={...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
%                      };

                 for ilines=1:length(thresh_cf_122)
%                     xvars{ilines}= 'Nd_{1.6} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7';                     
%                     xvars{ilines}='Nd from grid vals timeseries3';
                     xvars{ilines}= 'Nd factor';
                 end
                 
                 for ilines=1:length(thresh_cf_122)
                     errorvars{ilines} = 'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                 end 
             
             
             choose_Xbins_set=0;  Xbins_set=[0:0.1:3]; Xbins_set=10.^(Xbins_set);
             logbin_norm=0;
             
                 
           case {'Nd vs CTT all','Nd vs SZA','Nd vs SST','Nd vs CTH_no_subtract','Nd vs SST-CTT','Nd vs CTH','Nd vs CTT','Nd vs homog','Nd vs CTT std_dev','Nd vs tau-homog Cahalan','Nd vs minCTT'}                                      
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';


                    
%                  xvars={...
%                       'Nd from grid vals timeseries3'...
%                      ,'Nd from grid vals timeseries3'...                    
%                      };
%                  
%                  errorvars={...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                      'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
%                      };

                 for ilines=1:length(thresh_cf_122)
                     xvars{ilines}='Nd from grid vals timeseries3'
                 end
                 
                 for ilines=1:length(thresh_cf_122)
                     errorvars{ilines} = 'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                 end        
                 
            case 'CF vs CTH'
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';
                  
                 xvars={...
                      'Cloud Fraction from grid vals timeseries3'...
                     ,'Cloud Fraction from grid vals timeseries3'...                    
                     };
                 
                 %N.B. - not used for just sampling error
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
                     };
                 
            case {'CF_MOD35 vs CTH MOD35','CF_MOD35 vs CTH_Ryan MOD35','CF_MOD35 vs CTH_max MOD35','CF_MOD35 vs CTH_max_CTP_hybrid MOD35'}
                                ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';
                  
                 xvars={...
                     'MOD35 Cloud Fraction from grid vals timeseries3'...
%                      'CF all'...
%                     ,'CF all'...                    
                     };
                 
                 %N.B. - not used for just sampling error
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
                     };
                 
                 
                 
            case 'CF_all vs CTH all'  
                ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';
                  
                 xvars={...
%                     'MOD35 Cloud Fraction from grid vals timeseries3'...
                      'CF all'...
%                     ,'CF all'...                    
                     };
                 
                 %N.B. - not used for just sampling error
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
                     };
                 
                 
             case {'CF liq vs CF liq','CF ice vs CF ice','CF undet vs CF undet','CTT all vs CTT all','CTH all vs CTH all'}  
                ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';
                  
%                 xvars={...
%                      'Cloud Fraction from grid vals timeseries3'...
%                     ,'Cloud Fraction from grid vals timeseries3'...                    
%                     };
                 
                 xvars={...
                     'Dummy data for 1D'...
                     ,'Dummy data for 1D'...
                     };
                 
                 %N.B. - not used for just sampling error
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
                     };   
                 
              case 'CF ice vs CF ice'  
                ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';
                  
                 xvars={...
                      'CF ice'...
                     ,'CF ice'...                    
                     };
                 
                 %N.B. - not used for just sampling error
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...                    
                     };      
                 
                 
             
        end
                 


for ilines=1:length(thresh_cf_122)
%                screen_strs{ilines}='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';
%                screen_strs{ilines}='NP2 + CF2 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH'; %CF2 uses Fraction of domain in which Nd was available (as used for CTT, etc.)
%                screen_strs{ilines}='NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH'; %CF2 uses Fraction of domain in which Nd was available (as used for CTT, etc.)                
                screen_strs{ilines}='NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH + sigCTT'; %          
%                screen_strs{ilines}='NP4 + CF4 + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH'; %CF2 uses Fraction of domain in which Nd was available (as used for CTT, etc.)                                

%                screen_strs{ilines}='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH + Cahalan_tau_homog';
%                screen_strs{ilines}='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH';                
                screen_strs{ilines}='none';
%                screen_strs{ilines}='MEAN solarZA';
%                screen_strs{ilines}='MEAN solarZA + CF';
%                screen_strs{ilines}='MEAN solarZA + CF + NP';                
end
                
          switch param_vary
              case 'time'
%           if ivary_thresh==0
                    
                    time_mean_str_multi={....
                        ['DJF']...
                        ['MAM']...                        
                        ,['JJA']...
                        };

                    days_multi_specify={....
                        [336:366 1:60]...
                        [61:152]...                          %MAM
                        ,[153:244]...
                        };

                    thresh_cf_str=['[' num2str(thresh_cf_122{1}(1)) ' ' num2str(thresh_cf_122{1}(2)) ']'];
                    thresh_AZ_str=['[' num2str(thresh_AZ_122{1}(1)) ' ' num2str(thresh_AZ_122{1}(2)) ']'];
                    thresh_stdW_str=['[' num2str(thresh_stdW_122{1}(1)) ' ' num2str(thresh_stdW_122{1}(2)) ']'];
                    thresh_SZA_str=['[' num2str(thresh_solarZA_122{1}(1)) ' ' num2str(thresh_solarZA_122{1}(2)) ']'];
                    thresh_CTT_str=['[' num2str(thresh_CTT_122{1}(1)) ' ' num2str(thresh_CTT_122{1}(2)) ']'];
                    thresh_CTH_str=['[' num2str(thresh_CTH_122{1}(1)) ' ' num2str(thresh_CTH_122{1}(2)) ']'];
                    thresh_sigCTT_str=['[' num2str(thresh_sigCTT_122{1}(1)) ' ' num2str(thresh_sigCTT_122{1}(2)) ']'];                    

                    thresh_val_strs={...
                        ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                        ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                        };
                    
                    
                    thresh_val_strs={...
                        ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str '; thresh_sigCTT=' thresh_sigCTT_str ';'],...
                        ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str '; thresh_sigCTT=' thresh_sigCTT_str ';'],...
                        ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str '; thresh_sigCTT=' thresh_sigCTT_str ';'],...
                        };
                    
                    

                case 'thresh'
                
                    clear thresh_val_strs
                    for ithresh=1:length(xvars)
                        
                        if length(thresh_cf_122)<length(xvars)
                            ithresh2=1;
                        else
                            ithresh2=ithresh;
                        end

                        
                        thresh_cf_str=['[' num2str(thresh_cf_122{ithresh2}(1)) ' ' num2str(thresh_cf_122{ithresh2}(2)) ']'];
                        thresh_AZ_str=['[' num2str(thresh_AZ_122{ithresh2}(1)) ' ' num2str(thresh_AZ_122{ithresh2}(2)) ']'];
                        thresh_stdW_str=['[' num2str(thresh_stdW_122{ithresh2}(1)) ' ' num2str(thresh_stdW_122{ithresh2}(2)) ']'];
                        thresh_SZA_str=['[' num2str(thresh_solarZA_122{ithresh2}(1)) ' ' num2str(thresh_solarZA_122{ithresh2}(2)) ']'];
                        thresh_CTT_str=['[' num2str(thresh_CTT_122{ithresh2}(1)) ' ' num2str(thresh_CTT_122{ithresh2}(2)) ']'];
                        thresh_CTH_str=['[' num2str(thresh_CTH_122{ithresh2}(1)) ' ' num2str(thresh_CTH_122{ithresh2}(2)) ']'];
                        thresh_sigCTT_str=['[' num2str(thresh_sigCTT_122{ithresh2}(1)) ' ' num2str(thresh_sigCTT_122{ithresh2}(2)) ']'];
                        thresh_sens_str=['[' num2str(thresh_sensZA_122{ithresh2}(1)) ' ' num2str(thresh_sensZA_122{ithresh2}(2)) ']'];
                        
%                        thresh_val_strs{ithresh} = ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'];
                        thresh_val_strs{ithresh} = ['thresh_CF=' thresh_cf_str '; thresh_sensZA=' thresh_sens_str '; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str '; thresh_sigCTT=' thresh_sigCTT_str ';'];                        
                        
%                        time_mean_str_multi{ithresh} = ['CF=' thresh_cf_str];
                        time_mean_str_multi{ithresh} = ['CTT=' thresh_CTT_str ' K'];  
                        time_mean_str_multi{ithresh} = ['SZA=' thresh_SZA_str ''];  
%                        time_mean_str_multi{ithresh} = xvars{ithresh}(26:30); %band_str2;

                        days_multi_specify{ithresh}=[153:244];
                        days_multi_specify{ithresh}=[1:366];   
%                        days_multi_specify{ithresh}=[336:366 1:60]; 
%                        days_multi_specify{ithresh}=[153:244]; %153:244 = JJA

%solar DJF
%         days_required_for_mean = [354-45:366 1:mod(354+45,366)]; time_mean_str = 'solar DJF';
    %    %straight MAM
%        days_required_for_mean = [79-45:79+45]; time_mean_str = 'solar MAM';
    %    %straight JJA
%        days_multi_specify{ithresh}=[171-46:171+45]; %time_mean_str = 'solar JJA';
    %    %straight SON
%        days_required_for_mean = [263-46:263+45]; time_mean_str = 'solar
%        SON';


                    end
                    
              otherwise
                  for ithresh=1:length(xvars)
                  time_mean_str_multi{ithresh}='';
                  end
                       
                    
                    
%           end
                
                
          end
              
                 
%SW of UK                
%                  SST_multi_specify={....
%                       [273.15+10]...
%                      ,[273.15+18]...                     
%                      };   
                 
%                  %China
%                   SST_multi_specify={....
%                       [273.15+12]...
%                      ,[273.15+27]...                     
%                      };
%                  
%                  for isst=1:length(SST_multi_specify)
%                      labs_multi_specify{isst} =
%                      [time_mean_str_multi{isst} ' SST=' num2str(SST_multi_specify{isst})];
%                  end
                 
%                 for isst=1:length(SST_multi_specify)
%                     labs_multi_specify{isst} = [time_mean_str_multi{isst}];
%                 end
                 
                 labs_multi_specify = time_mean_str_multi;
                 
                 for ilab=1:length(labs_multi_specify)
                     if length(strfind(labs_multi_specify{ilab},'SZA'))>0
                         labs_multi_specify{ilab} = ['SZA=' labs_multi_specify{ilab}(6:7) '-' labs_multi_specify{ilab}(9:10) '^{o}'];
                     end

                 end

                 
                
%              xlab_multi = ['N_d (cm^{-3}) for CF ' num2str(thresh_cf_122(1),'%1.2f') '-' num2str(thresh_cf_122(2),'%1.2f')];                           
              xlab_multi = [' for CF ' num2str(thresh_cf_122{1}(1),'%1.2f') '-' num2str(thresh_cf_122{1}(2),'%1.2f')];                                               
                
                
                  
              
        
        
        clear allSZA_pdfs highSZA_pdf xdat_multi ydat_multi labs_multi errordatU_multi errordatL_multi sampling_errorU sampling_errorL
% -------------------------------------------------------------------------       
%    loop through the multiple runs of plotTimeHeightVap3        
% -------------------------------------------------------------------------

        for idat_multi=1:length(xvars)
        
            ioverride_pdf_varchoose=1; %override the defaults
            ioverride_pdf=1; %override the defaults  
            
            %flags to allow the direct specification of Xbins, Ybins, Zbins
            ichoose_Xbins=ichoose_Xbins_set; Xbins=Xbins_set; 
            ichoose_Ybins=ichoose_Ybins_set; Ybins=Ybins_set; 
            ichoose_Zbins=ichoose_Zbins_set; Zbins=Zbins_set; 
            
                                    
            i577 = 'MODIS_plot_UW';
            x_axis_vals = xvars{idat_multi};
            
            if length(thresh_val_strs)<length(xvars)
                idat_multi2=1;
            else
                idat_multi2=idat_multi;
            end
            screen_type = screen_strs{idat_multi2};
            eval(thresh_val_strs{idat_multi2});
            time_mean_str = time_mean_str_multi{idat_multi2};
            days_required_for_mean = days_multi_specify{idat_multi2};
%            SST = SST_multi_specify{idat_multi};
            
            man_choose_plotTimeHeight_graph=1;
            logflag=0;
            dlogflag=0;
            noplot=1; %flag to say to just do the calculations and not to plot
            nXpdf=nXpdf_131;
            nYpdf=N_sza_bins;
            ndims_hist=2;
            ioverride_time_selection=1;
         %% --- Does calculations based on plotTimeHeightVap3            
            plotTimeHeightVap3
%            close(gcf);

%sampling error = std_dev / sqrt(N)
            sampling_errorU(idat_multi).dat = std_dev_X ./ sqrt(NX_vals);
            sampling_errorL(idat_multi).dat = std_dev_X ./ sqrt(NX_vals);
            
            stdev_multi(idat_multi).dat = std_dev_X;
            
            [meanXvals(idat_multi),Nvals,stdXvals(idat_multi)] = meanNoNan(X,1);
            [meanYvals(idat_multi),Nvals,stdYvals(idat_multi)] = meanNoNan(Y,1);            
            
            Xvals_save(idat_multi).dat = X;
            Yvals_save(idat_multi).dat = Y;
            
            for isza=1:size(qh,1)-1
                allSZA_pdfs{isza}(idat_multi).dat = expand_PDF(mid_Xbins,qh(isza,1:end-1));                
            end
            
           
            
            man_choose_water_graph=1;
            switch mean_or_pdf
                case 'pdf'
                    graph=977; %1d pdf
                case 'mean'
                    graph=96; %mean of 2D PDF
                case 'cf_array'
                    graph=961; %case where populate an array of CF vs height for each datapoint with the observed CF at the observed 
                    %CTH,, but with zeros above it and NaNs below. If have
                    %no CF then is all zeros - trying to match what CALIPSO
                    %3d CF does.
%                    graph=966; axis1D = 'y';
            end
            
            noplot=1;
            waterVapourMay2005
%            close(gcf);

man_choose_water_graph=1;

            
            

            xdat_multi(idat_multi).x = xdat(1).x;
            ydat_multi(idat_multi).y = ydat(1).y';
%            labs_multi(idat_multi).l = xvars{idat_multi};               
            labs_multi(idat_multi).l = labs_multi_specify{idat_multi2};   

            %only remove points if we are doing a mean (since pdfs show no.
            %points anyway).
            switch mean_or_pdf
                case 'mean'
                    iremove_points_122=1;
                    thresh_points_122=2;
                    if iremove_points_122==1
                        iremove_122=find(NX_vals<thresh_points_122);    
                        xdat_multi(idat_multi).x(iremove_122)=NaN;
                    end

            end
            
            
            if plot_error_bars==1
                switch error_type
                    case 'sampling'
                        errordatU_multi(idat_multi).dat = sampling_errorU(idat_multi).dat;
                        errordatL_multi(idat_multi).dat = sampling_errorL(idat_multi).dat;
                        
                    case 'standard deviation'
                        errordatU_multi(idat_multi).dat = stdev_multi(idat_multi).dat;
                        errordatL_multi(idat_multi).dat = stdev_multi(idat_multi).dat;                        
                        
                    case 'bootstrap sampling'
                        ioverride_bootstrap=1;
                        xvars_boot = xvars{idat_multi};  
                        bootstrap_from_2d_PDF
                        errordatU_multi(idat_multi).dat = err_bootstrap;
                        errordatL_multi(idat_multi).dat = err_bootstrap;

                    case {'combined','mean_percent','combined + sampling'}
                        %for the errors we do a 2D pdf of the absolute uncertainty
                        %(x-axis) vs SZA (y). This calculates the sqrt of the  sum of the squares
                        %of all the uncertainties divided by N, which should be the
                        %combined error (if the errrors are assumed random)

                        ioverride_pdf_varchoose=1; %override the defaults
                        ioverride_pdf=1;
                        
                        %flags to allow the direct specification of Xbins, Ybins, Zbins
                        ichoose_Xbins=ichoose_Xbins_set; Xbins=Xbins_set; 
                        ichoose_Ybins=ichoose_Ybins_set; Ybins=Ybins_set;
                        ichoose_Zbins=ichoose_Zbins_set; Zbins=Zbins_set;
                        
                        x_axis_vals = errorvars{idat_multi};
                        screen_type = screen_strs{idat_multi};
                        eval(thresh_val_strs{idat_multi});
                        nXpdf=nXpdf_131;
                        nYpdf=N_sza_bins;
                        man_choose_plotTimeHeight_graph=1;
                        logflag=0;
                        dlogflag=0;
                        noplot=1; %flag to say to just do the calculations and not to plot
                        ndims_hist=2;
                 %% --- Does calculations based on plotTimeHeightVap3
                        plotTimeHeightVap3
                        %                close(gcf);

                        %this is the mean of the sqrt of the sum of the square of
                        %the error
                        
                        switch error_type
                            case 'combined + sampling'
                                 errordatU_multi(idat_multi).dat = sqrt ( rms_X.^2 + sampling_errorU(idat_multi).dat.^2 );
                                 errordatL_multi(idat_multi).dat = sqrt ( rms_X.^2 + sampling_errorU(idat_multi).dat.^2 );
                        
                            case 'combined'
                                errordatU_multi(idat_multi).dat = rms_X; %from PlotTimeHeightVap3
                                errordatL_multi(idat_multi).dat = rms_X;
                            case 'mean_percent'
                                errordatU_multi(idat_multi).dat = X_mean;
                                errordatL_multi(idat_multi).dat = X_mean;   
                        end

                end
                  
                  
              titlenam = [titlenam ', ' remove_character(error_type,'_','-') ' error bars'];

            end
        
        end    
        
        iwrite_text_dat=0;
        noplot=0;
                                  
        xdat = xdat_multi;
        ydat = ydat_multi;
        labs = labs_multi;
        xlab = [xlab xlab_multi];
        
        if plot_error_bars==1
            errordatU = errordatU_multi;
            errordatL = errordatL_multi;
            ierror_bars = ierror_bars_multi;
        end

        
%        titlenam = ['Nd vs SZA for lat=' num2str(MLAT(ilat)) ' and lon=' num2str(MLON(ilon))];
aqua_terra_str='';
if exist('aqua_terra_timeseries3')
    for istr=1:length(aqua_terra_timeseries3)
        aqua_terra_str=[aqua_terra_str aqua_terra_timeseries3{istr} ' '];
    end
end


if iadd_data==1
    for ilab=1:length(labs)
        labs(ilab).l = [labs(ilab).l ' MODIS'];
    end

    xdat(end+1) = xcal; %got this by running read_cfads_calipso_monthly_night_ipsl
    ydat(end+1) = ycal; %and then case 132 of watervap (and copied the data into xcal,ycal)
    labs(end+1).l = 'CALIPSO Cloud Fraction';
    ierror_bars=0;
    line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
    
    xdat(end+1) = xryan; %
    ydat(end+1) = yryan; %
    labs(end+1).l = 'MODIS CTH (Ryan)';
    ierror_bars=0;
    line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='s'; istyle=istyle+1;
    
    xdat(end+1) = xCTH_max; %
    ydat(end+1) = yCTH_max; %
    labs(end+1).l = 'MODIS max CTH';
    ierror_bars=0;
    line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.7]; marker_style(istyle).m='^'; istyle=istyle+1;  
    
    xdat(end+1) = xCTP_max; %Get this by running 'CF_MOD35 vs CTH_max_CTP_hybrid MOD35' in case 131 
    ydat(end+1) = yCTP_max; %and saving xdat as xCTP_max
    labs(end+1).l = 'MODIS max CTP';
    ierror_bars=0;
    line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.7 0.4 0.4]; marker_style(istyle).m='>'; istyle=istyle+1;  
    
%    xlab = 'Cloud top height or Height (km)';
%    ylab = 'Bin normalised frequency or Cloud fraction';    
    
    xlab = 'Height (km)';
    ylab = 'Cloud fraction';        

end

    
        titlenam = [titlenam ' ' aqua_terra_str];

        figname=titlenam;
        savename=figname;

%        xlims=1;
%        xlimits=1000*[0 0.025]; %set in case number 96 (mean) or 977 (pdf)

        izlim=1;
%        zmin=1500;
%        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

%        ylab='SZA'; %set in case number 96
%        xlab='UTC Time';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

%        posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];

    
    
        
        
        
  case 130
        %plot PDFS of tau and Re
        
        ioverride_122=1;
        
        Npdf_pts=40;
        
        
        graph=128;
        noplot=1;
%        man_choose_water_graph=1; waterVapourMay2005

noplot=0;
ichoose_styles=1;
        
        ierror_bars='none';
        
        ilab=1;
        line_labs{ilab} = 'LowCF, Low Sens'; ilab=ilab+1;
        line_labs{ilab} = 'LowCF, High Sens'; ilab=ilab+1;
        line_labs{ilab} = 'HighCF, Low Sens'; ilab=ilab+1;
        line_labs{ilab} = 'HighCF, High Sens'; ilab=ilab+1;  
        
        ilab=1;
        line_cols{ilab} = [0 0 1]; ilab=ilab+1;
        line_cols{ilab} = [0 1 0]; ilab=ilab+1;
        line_cols{ilab} = [1 0 0]; ilab=ilab+1;
        line_cols{ilab} = [0.4 0.4 0.4]; ilab=ilab+1;  
        
        
        
         istyle = 1;
         
        
idat2=0;
        for idat=1:length(Tau_all_pdfs{1})
            idat2=idat2+1;
            a=Tau_all_pdfs{1}(idat).dat;
            Xbins=make_PDF_bins(a,Npdf_pts);
            ydat(idat2).y = ndHistc_run([a],Xbins);
            ydat(idat2).y = ydat(idat2).y / sum(ydat(idat2).y);
            xdat(idat2).x = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );
            labs(idat2).l = [line_labs{idat} ', low SZA'];
            line_pattern(istyle).p= '-';  line_colour(istyle).c=line_cols{idat}; marker_style(istyle).m='d'; istyle=istyle+1;
            
            idat2=idat2+1;            
            a=Tau_all_pdfs{12}(idat).dat;
            Xbins=make_PDF_bins(a,Npdf_pts);
            ydat(idat2).y = ndHistc_run([a],Xbins);
            ydat(idat2).y = ydat(idat2).y / sum(ydat(idat2).y);
            xdat(idat2).x = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );
            labs(idat2).l = [line_labs{idat} ', high SZA'];
            line_pattern(istyle).p= '--';  line_colour(istyle).c=line_cols{idat}; marker_style(istyle).m='d'; istyle=istyle+1;
            
        end
        
        
        

%titlenam = [titlenam Re_or_tau_constant ' CONSTANT, varying scale_fac'];

%xlab = ['N_d (cm^{-3}) for ' low_or_high_CF];
%xlab = ['N_d (cm^{-3})'];
xlab = ['Optical Depth'];

ylab = 'Frequency';


%calculates the fraction of the overall Nd increase (from low to high SZA)
%that the constant Re or constant tau increase accounts for (and writes the
%table)
%make_sza_rel_increases_table_vals



                     
                     plot_error_bars = 0;


%                      istyle = 1;
%                      line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
%                      line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
%                      line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%                      line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%                      line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
%                      line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
% 



        
        
  
                     
        
    case 129
        %calculate the N lines (vs SZA) when hold tau constant and then Re
        %constant
        
        ioverride_122=1;
        

        man_choose_water_graph=1;
        graph=128;
        noplot=1;
        waterVapourMay2005
        
        %saves these lines:-
%          tau_vs_sza (length 4 - low CF+low sensor, low CF+high sensor, etc.)
%          Re_vs_sza_lowCF (rest are length length 6)
%          Re_vs_sza_highCF
%          Nd_vs_sza_lowCF
%          Nd_vs_sza_highCF

low_or_high_CF = 'lowCF';
%low_or_high_CF = 'highCF';
Re_or_tau_constant = 'Re';
%Re_or_tau_constant = 'Tau';

titlenam = [titlenam Re_or_tau_constant ' CONSTANT, varying scale_fac'];

xlab = ['N_d (cm^{-3}) for ' low_or_high_CF];

switch low_or_high_CF
    case 'lowCF'

        xdat = Nd_vs_sza_lowCF;
        ydat = ydat_Nd_lowCF;
        labs = labs_Nd_lowCF;
        errorU = errorU_Nd_lowCF;
        errorL = errorL_Nd_lowCF;
        Lxdat = length(xdat);
        line_pattern = patt_Nd_vs_sza;
        line_colour = line_colour_Nd_vs_sza;
        marker_style = marker_Nd_vs_sza;
        
        for idat=1:length(Nd_vs_sza_lowCF)
            itau = 2-mod(idat,2); %laternating 1s and 2s for low CF
            %will need 3s & 4s for high CF, so use 4-mod(idat,2)

            for isza=1:length(Re_all_pdfs_lowCF)
                %makes a big difference whether we find the mean
                %tau and then apply the ^0.5 factor, or whether we
                %apply the power to each tau and then take the mean
                % *** make sure that (idat) is used for Re and
                % (itau) for Tau for Re_all_pdfs_lowCF{isza}(idat)
                Re_vals(isza) = mean( (Re_all_pdfs_lowCF{isza}(idat).dat).^-2.5 );
                Tau_vals(isza) = mean( (Tau_all_pdfs{isza}(itau).dat).^0.5 );
                %as would have been used in the actual calc of mean
                %Nd
                prod_vals(isza) = mean( Re_all_pdfs_lowCF{isza}(idat).dat.^-2.5 .* Tau_all_pdfs{isza}(itau).dat.^0.5 );
            end
            const = Nd_vs_sza_lowCF(idat).x(1) / prod_vals(1); %const from N1 = k * mean(Tau^0.5 .* Re.^-2.5)
%            scale_fac = Re_vals(1)*Tau_vals(1) / prod_vals(1);
            scale_fac = Nd_vs_sza_lowCF(idat).x' ./ (Re_vals.*Tau_vals);
%            scale_fac(idat).dat = Re_vals.*Tau_vals ./ prod_vals;
            
            %scale_fac_lowCF_save(idat).dat = Re_vals.*Tau_vals ./ prod_vals;
            scale_fac_lowCF_save(idat).dat = scale_fac;
            
            
            switch Re_or_tau_constant
                case 'Re'

%                    xdat(Lxdat+idat).x = (Nd_vs_sza_lowCF(idat).x(1) ./ tau_vs_sza(itau).x(1).^0.5 ) .* tau_vs_sza(itau).x.^0.5 ;
%                    xdat(Lxdat+idat).x = (Nd_vs_sza_lowCF(idat).x(1) ./ vals_vs_sza(1) ) .* vals_vs_sza ;   
                    
                    xdat(Lxdat+idat).x = Re_vals(1) .* Tau_vals .* scale_fac ;  
                    labs(Lxdat+idat).l = [labs(idat).l ' Re constant'];                    
                case 'Tau'
%                    xdat(Lxdat+idat).x = Nd_vs_sza_lowCF(idat).x(1) * (Re_vs_sza_lowCF(idat).x ./ Re_vs_sza_lowCF(idat).x(1) ).^-2.5 ;                    
%                    labs(Lxdat+idat).l = [labs(idat).l ' Tau constant'];
                    

%                    xdat(Lxdat+idat).x = (Nd_vs_sza_lowCF(idat).x(1) ./ tau_vs_sza(itau).x(1).^0.5 ) .* tau_vs_sza(itau).x.^0.5 ;
%                    xdat(Lxdat+idat).x = Nd_vs_sza_lowCF(idat).x(1) ./ vals_vs_sza(1)  .* vals_vs_sza ;      
                    xdat(Lxdat+idat).x = Tau_vals(1) .* Re_vals .* scale_fac ;                          
                    labs(Lxdat+idat).l = [labs(idat).l ' Tau constant'];    
                    
                    
            end
            ydat(Lxdat+idat).y = ydat(idat).y;

            errordatU(Lxdat+idat).dat=zeros(size(xdat(Lxdat+idat).x));
            errordatL(Lxdat+idat).dat=zeros(size(xdat(Lxdat+idat).x));   
        end

    case 'highCF'

        xdat = Nd_vs_sza_highCF;
        ydat = ydat_Nd_highCF;
        labs = labs_Nd_highCF;
        errorU = errorU_Nd_highCF;
        errorL = errorL_Nd_highCF;
        Lxdat = length(xdat);
        line_pattern = patt_Nd_vs_sza;
        line_colour = line_colour_Nd_vs_sza;
        marker_style = marker_Nd_vs_sza;        
        
        for idat=1:length(Nd_vs_sza_highCF)
            itau = 4-mod(idat,2); %laternating 1s and 2s for low CF
            %will need 3s & 4s for high CF, so use 4-mod(idat,2)

            for isza=1:length(Re_all_pdfs_highCF)
                %makes a big difference whether we find the mean
                %tau and then apply the ^0.5 factor, or whether we
                %apply the power to each tau and then take the mean
                % *** make sure that (idat) is used for Re and
                % (itau) for Tau for Re_all_pdfs_highCF{isza}(idat)
                Re_vals(isza) = mean( (Re_all_pdfs_highCF{isza}(idat).dat).^-2.5 );
                Tau_vals(isza) = mean( (Tau_all_pdfs{isza}(itau).dat).^0.5 );
                %as would have been used in the actual calc of mean
                %Nd
                prod_vals(isza) = mean( Re_all_pdfs_highCF{isza}(idat).dat.^-2.5 .* Tau_all_pdfs{isza}(itau).dat.^0.5 );
            end
            const = Nd_vs_sza_highCF(idat).x(1) / prod_vals(1); %const from N1 = k * mean(Tau^0.5 .* Re.^-2.5)
%            scale_fac = Re_vals(1)*Tau_vals(1) / prod_vals(1);
            scale_fac = Nd_vs_sza_highCF(idat).x' ./ (Re_vals.*Tau_vals);
            %scale_fac = Re_vals.*Tau_vals ./ prod_vals(1);
            
%            scale_fac_highCF_save(idat).dat = Re_vals.*Tau_vals ./ prod_vals;
            scale_fac_highCF_save = scale_fac;
            
            switch Re_or_tau_constant
                case 'Re'

%                    xdat(Lxdat+idat).x = (Nd_vs_sza_highCF(idat).x(1) ./ tau_vs_sza(itau).x(1).^0.5 ) .* tau_vs_sza(itau).x.^0.5 ;
%                    xdat(Lxdat+idat).x = (Nd_vs_sza_highCF(idat).x(1) ./ vals_vs_sza(1) ) .* vals_vs_sza ;   
                    
                    xdat(Lxdat+idat).x = Re_vals(1) .* Tau_vals .* scale_fac ;  
                    labs(Lxdat+idat).l = [labs(idat).l ' Re constant'];                    
                case 'Tau'
%                    xdat(Lxdat+idat).x = Nd_vs_sza_highCF(idat).x(1) * (Re_vs_sza_highCF(idat).x ./ Re_vs_sza_highCF(idat).x(1) ).^-2.5 ;                    
%                    labs(Lxdat+idat).l = [labs(idat).l ' Tau constant'];
                    

%                    xdat(Lxdat+idat).x = (Nd_vs_sza_highCF(idat).x(1) ./ tau_vs_sza(itau).x(1).^0.5 ) .* tau_vs_sza(itau).x.^0.5 ;
%                    xdat(Lxdat+idat).x = Nd_vs_sza_highCF(idat).x(1) ./ vals_vs_sza(1)  .* vals_vs_sza ;      
                    xdat(Lxdat+idat).x = Tau_vals(1) .* Re_vals .* scale_fac ;                          
                    labs(Lxdat+idat).l = [labs(idat).l ' Tau constant'];    
                    
                    
            end
            ydat(Lxdat+idat).y = ydat(idat).y;

            errordatU(Lxdat+idat).dat=zeros(size(xdat(Lxdat+idat).x));
            errordatL(Lxdat+idat).dat=zeros(size(xdat(Lxdat+idat).x));   
        end
        


      
        
        
end

%calculates the fraction of the overall Nd increase (from low to high SZA)
%that the constant Re or constant tau increase accounts for (and writes the
%table)
make_sza_rel_increases_table_vals



                     
                     plot_error_bars = 1;


istyle = length(line_pattern)+1;
                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
                     line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
                     line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
                     line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;




        
        
    case 128
        
        %Need rel changes in
        % dN_2.1, tau, dN_dtau, dRe_2.1, dN_dRe_2.1, dN_1.6, dRe_1.6, dN_dRe_1.6, 
        % dN_3.7, dRe_3.7, dN_dRe_3.7, 
        
        clear rel_incs_low_low rel_incs_low_high rel_incs_high_low rel_incs_high_high
        
        iplot_128 = 0; %flag to say whether to plot all the graphs or not
        iplot_128 = 1;
        
        preserve_flags=1;
        
        iwrite_table=0;
        
        izlim_128 = 0;
        zlims_128 = [0 90];
        
        savedir_128 = 'SZA_plots/vs_CF/';
        tag_128 = '_lowSZA_COLD_lowSensZA'; 
        
        scrsz=get(0,'ScreenSize');
        posit=[scrsz];

        

        
        
%% Optical depth, both CF ranges        
        
        ioverride_122=1;
        
%         multi_plot_case_122 = 'Nd, different wavelengths vs SZA';
%        multi_plot_case_122 = 'Re, different wavelengths vs SZA';  
        multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
        
        %   var_plot_122_21mum = 'N_d (cm^{-3})';
        %               var_plot_122_21mum = 'R_{eff} (\mum)';
        var_plot_122_21mum = 'Optical Depth';

        man_choose_water_graph=1;
        graph=122;
        if iplot_128==1
            noplot=0;
        else
            noplot=1;
        end
        ioverride_watervap_newfig=1;
        waterVapourMay2005
        if izlim_128==1
            set(gca,'ylim',zlims_128);
        end
        saveas_ps_fig_emf(gcf,[savedir savedir_128 'tau'],tag_128,0,0);
        
        %save the vs SZA lines for doing the constant Tau varying Re
        %and vice versa lines
        tau_vs_sza = xdat; 
        
        Tau_pdfs = highSZA_pdf;  %highSZA_pdf(i).dat
        Tau_all_pdfs = allSZA_pdfs;  %cell structure of pdfs allSZA_pdfs{isza_bin}(iplot_line).dat
        %where iplot_line=1 is the Re_1.6 low sens, =2 is the Re_1.6 high
        %sens, etc. through to iplot_line=6
        
        
        if iwrite_table==1
            %this produces low CF for low and high sensZA and then high CF for
            %low and high sens
            rel_incs_low_low([2])   = rel_incs([1]);
            rel_incs_low_high([2])   = rel_incs([2]);
            rel_incs_high_low([2])   = rel_incs([3]);
            rel_incs_high_high([2])   = rel_incs([4]);

            clear tau_lowSZA_lowCF
            for ilow=1:2
                tau_lowSZA_lowCF(ilow) = xdat(ilow).x(1);
            end

            clear tau_lowSZA_highCF
            for ilow=3:length(xdat)
                tau_lowSZA_highCF(ilow-2) = xdat(ilow).x(1);
            end

        end
        
        
      
       % *************************************************************************** 
%% Re CF<0.8        

         ioverride_122=1;
        
         multi_plot_case_122 = 'Nd, different wavelengths vs SZA';
%         multi_plot_case_122 = 'Re, different wavelengths vs SZA';  
%        multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
         multi_wavelength_case='Re';

         thresh_cf_122 = [0.1 0.8];
%         thresh_cf_122 = [0.8 1.0];
         thresh_cf_122 = [0.1 1.0];
         
        man_choose_water_graph=1;
        graph=122;
        if iplot_128==1
            noplot=0;
        else
            noplot=1;
        end
        ioverride_watervap_newfig=1;
        waterVapourMay2005
        if izlim_128==1
            set(gca,'ylim',zlims_128);
        end
        saveas_ps_fig_emf(gcf,[savedir savedir_128 'Re_1'],tag_128,0,0);
        
        %save the vs SZA lines for doing the constant Tau varying Re
        %and vice versa lines
        Re_vs_sza_lowCF = xdat; 
        
        if iwrite_table==1

            rel_incs_low_low([4 7 10])   = rel_incs([3 1 5]);
            rel_incs_low_high([4 7 10])   = rel_incs([4 2 6]);

        end

        Re_pdfs_lowCF = highSZA_pdf;  %highSZA_pdf(i).dat 
        Re_all_pdfs_lowCF = allSZA_pdfs; %pdfs from all sza bins
        
        clear Re_lowSZA_highCF
        for ilow=1:length(xdat)
            Re_lowSZA_lowCF(ilow) = xdat(ilow).x(1);
        end
        
        
%% Re CF>0.8        
          ioverride_122=1;
        
         multi_plot_case_122 = 'Nd, different wavelengths vs SZA';
%         multi_plot_case_122 = 'Re, different wavelengths vs SZA';  
%        multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
         multi_wavelength_case='Re';

          thresh_cf_122 = [0.1 0.8];
          thresh_cf_122 = [0.8 1.0];


        man_choose_water_graph=1;
        graph=122;
        if iplot_128==1
            noplot=0;
        else
            noplot=1;
        end        
        ioverride_watervap_newfig=1;
        waterVapourMay2005
        if izlim_128==1
            set(gca,'ylim',zlims_128);
        end        
        saveas_ps_fig_emf(gcf,[savedir savedir_128 'Re_2'],tag_128,0,0);
        
        %save the vs SZA lines for doing the constant Tau varying Re
        %and vice versa lines
        Re_vs_sza_highCF = xdat; 
        
        if iwrite_table==1
            rel_incs_high_low([4 7 10])   = rel_incs([3 1 5]);
            rel_incs_high_high([4 7 10])   = rel_incs([4 2 6]);
        end

        Re_pdfs_highCF = highSZA_pdf;  %highSZA_pdf(i).dat 
        Re_all_pdfs_highCF = allSZA_pdfs; %pdfs from all sza bins
        
        clear Re_lowSZA_highCF
        for ilow=1:length(xdat)
            Re_lowSZA_highCF(ilow) = xdat(ilow).x(1);
        end
        
        
        % *****************************************************************
%% Nd CF<0.8        
        
        
         ioverride_122=1;
        
         multi_plot_case_122 = 'Nd, different wavelengths vs SZA';
%        multi_plot_case_122 = 'Re, different wavelengths vs SZA';  
%        multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
         multi_wavelength_case='Nd';

          thresh_cf_122 = [0.1 0.8];
%          thresh_cf_122 = [0.8 1.0];
          thresh_cf_122 = [0.1 1.0];


        man_choose_water_graph=1;
        graph=122;
        if iplot_128==1
            noplot=0;
        else
            noplot=1;
        end
        ioverride_watervap_newfig=1;
        waterVapourMay2005
        if izlim_128==1
            set(gca,'ylim',zlims_128);
        end        
        saveas_ps_fig_emf(gcf,[savedir savedir_128 'Nd_1'],tag_128,0,0);
        
        %save the vs SZA lines for doing the constant Tau varying Re
        %and vice versa lines
        Nd_vs_sza_lowCF = xdat; 
        ydat_Nd_lowCF = ydat;
        labs_Nd_lowCF = labs;
        errorU_Nd_lowCF = errordatU;
        errorL_Nd_lowCF = errordatL;     
        patt_Nd_vs_sza = line_pattern;
        line_colour_Nd_vs_sza = line_colour;
        marker_Nd_vs_sza = marker_style;

        if iwrite_table==1
            rel_incs_low_low([1 6 9])   = rel_incs([3 1 5]);
            rel_incs_low_high([1 6 9])   = rel_incs([4 2 6]);
        end

        Nd_pdfs_lowCF = highSZA_pdf;  %highSZA_pdf(i).dat
        Nd_all_pdfs_lowCF= allSZA_pdfs; %pdfs from all sza bins
        
        clear Nd_lowSZA_lowCF
        for ilow=1:length(xdat)
            Nd_lowSZA_lowCF(ilow) = xdat(ilow).x(1);
        end
        
        
        
           ioverride_122=1;
%% Nd CF>0.8        
         multi_plot_case_122 = 'Nd, different wavelengths vs SZA';
%        multi_plot_case_122 = 'Re, different wavelengths vs SZA';  
%        multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
         multi_wavelength_case='Nd';
         
%          thresh_cf_122 = [0.1 0.8];
          thresh_cf_122 = [0.8 1.0];

        man_choose_water_graph=1;
        graph=122;
        if iplot_128==1
            noplot=0;
        else
            noplot=1;
        end
        ioverride_watervap_newfig=1;        
        waterVapourMay2005
        if izlim_128==1
            set(gca,'ylim',zlims_128);
        end        
        saveas_ps_fig_emf(gcf,[savedir savedir_128 'Nd_2'],tag_128,0,0);
        
        %save the vs SZA lines for doing the constant Tau varying Re
        %and vice versa lines
        Nd_vs_sza_highCF = xdat; 
        ydat_Nd_highCF = ydat;
        labs_Nd_highCF = labs;
        errorU_Nd_highCF = errordatU;
        errorL_Nd_highCF = errordatL;
        
        if iwrite_table==1
            rel_incs_high_low([1 6 9])   = rel_incs([3 1 5]);
            rel_incs_high_high([1 6 9])   = rel_incs([4 2 6]);
        end

        Nd_pdfs_highCF = highSZA_pdf;  %highSZA_pdf(i).dat 
        Nd_all_pdfs_highCF = allSZA_pdfs; %pdfs from all sza bins
        
        clear Nd_lowSZA_highCF
        for ilow=1:length(xdat)
            Nd_lowSZA_highCF(ilow) = xdat(ilow).x(1);
        end
        
        
        
        
        

        if iwrite_table==1

       % dN due to dtau (keeping Re constant)
       ival=3;   
       N1 = Nd_lowSZA_lowCF(ival);
       N12 = (Tau_pdfs(1).dat./ tau_lowSZA_lowCF(1)).^0.5 * N1;
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_low([3]) = mean( N12 - N1 ) ./ mean((N2 - N1));
       
       ival=4;
       N12 = (Tau_pdfs(2).dat./ tau_lowSZA_lowCF(2)).^0.5 * Nd_lowSZA_lowCF(ival);
       N1 = Nd_lowSZA_lowCF(ival);
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_high([3]) = mean( N12 - Nd_lowSZA_lowCF(ival)) ./ mean((N2 - N1));
       
       ival=3;
       N12 = (Tau_pdfs(3).dat./ tau_lowSZA_highCF(1)).^0.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_low([3]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       ival=4;
       N12 = (Tau_pdfs(4).dat./ tau_lowSZA_highCF(2)).^0.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_high([3]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));

       
       % dN due to dRe_2.1 (keeping Tau constant)
       ival=3;   
       N1 = Nd_lowSZA_lowCF(ival);
       N12 = (Re_pdfs_lowCF(ival).dat./ Re_lowSZA_lowCF(ival)).^-2.5 * N1;
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_low([5]) = mean( N12 - N1 ) ./ mean((N2 - N1));
       
       ival=4;
       N12 = (Re_pdfs_lowCF(ival).dat./ Re_lowSZA_lowCF(ival)).^-2.5 * Nd_lowSZA_lowCF(ival);
       N1 = Nd_lowSZA_lowCF(ival);
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_high([5]) = mean( N12 - Nd_lowSZA_lowCF(ival)) ./ mean((N2 - N1));
       
       ival=3;
       N12 = (Re_pdfs_highCF(ival).dat./ Re_lowSZA_highCF(ival)).^-2.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_low([5]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       ival=4;
       N12 = (Re_pdfs_highCF(ival).dat./ Re_lowSZA_highCF(ival)).^-2.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_high([5]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       
       % dN due to dRe_1.6 (keeping Tau constant)
       ival=1;       
       N12 = (Re_pdfs_lowCF(ival).dat./ Re_lowSZA_lowCF(ival)).^-2.5 * Nd_lowSZA_lowCF(ival);
       N1 = Nd_lowSZA_lowCF(ival);
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_low([8]) = mean( N12 - Nd_lowSZA_lowCF(ival)) ./ mean((N2 - N1));
       
       ival=2;
       N12 = (Re_pdfs_lowCF(ival).dat./ Re_lowSZA_lowCF(ival)).^-2.5 * Nd_lowSZA_lowCF(ival);
       N1 = Nd_lowSZA_lowCF(ival);
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_high([8]) = mean( N12 - Nd_lowSZA_lowCF(ival)) ./ mean((N2 - N1));
       
       ival=1;
       N12 = (Re_pdfs_highCF(ival).dat./ Re_lowSZA_highCF(ival)).^-2.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_low([8]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       ival=2;
       N12 = (Re_pdfs_highCF(ival).dat./ Re_lowSZA_highCF(ival)).^-2.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_high([8]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       
       % dN due to dRe_3.7 (keeping Tau constant)
       ival=5;       
       N12 = (Re_pdfs_lowCF(ival).dat./ Re_lowSZA_lowCF(ival)).^-2.5 * Nd_lowSZA_lowCF(ival);
       N1 = Nd_lowSZA_lowCF(ival);
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_low([11]) = mean( N12 - Nd_lowSZA_lowCF(ival)) ./ mean((N2 - N1));
       
       ival=6;
       N12 = (Re_pdfs_lowCF(ival).dat./ Re_lowSZA_lowCF(ival)).^-2.5 * Nd_lowSZA_lowCF(ival);
       N1 = Nd_lowSZA_lowCF(ival);
       N2 = Nd_pdfs_lowCF(ival).dat;
       rel_incs_low_high([11]) = mean( N12 - Nd_lowSZA_lowCF(ival)) ./ mean((N2 - N1));
       
       ival=5;
       N12 = (Re_pdfs_highCF(ival).dat./ Re_lowSZA_highCF(ival)).^-2.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_low([11]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       ival=6;
       N12 = (Re_pdfs_highCF(ival).dat./ Re_lowSZA_highCF(ival)).^-2.5 * Nd_lowSZA_highCF(ival);
       N1 = Nd_lowSZA_highCF(ival);
       N2 = Nd_pdfs_highCF(ival).dat;
       rel_incs_high_high([11]) = mean( N12 - Nd_lowSZA_highCF(ival)) ./ mean((N2 - N1));
       
       
      %print the table to file in Latex format 
      fid = fopen('/home/disk/eos1/d.grosvenor/matlab/work/MODIS/Nd_attr_table.txt','wt');
                   
      fprintf(fid,'low CF, low sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_low_low(1),rel_incs_low_low(2),rel_incs_low_low(3),rel_incs_low_low(4),rel_incs_low_low(5),rel_incs_low_low(6),rel_incs_low_low(7),rel_incs_low_low(8),rel_incs_low_low(9),rel_incs_low_low(10),rel_incs_low_low(11));
      fprintf(fid,'low CF, high sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_low_high(1),rel_incs_low_high(2),rel_incs_low_high(3),rel_incs_low_high(4),rel_incs_low_high(5),rel_incs_low_high(6),rel_incs_low_high(7),rel_incs_low_high(8),rel_incs_low_high(9),rel_incs_low_high(10),rel_incs_low_high(11));
      fprintf(fid,'high CF, low sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_high_low(1),rel_incs_high_low(2),rel_incs_high_low(3),rel_incs_high_low(4),rel_incs_high_low(5),rel_incs_high_low(6),rel_incs_high_low(7),rel_incs_high_low(8),rel_incs_high_low(9),rel_incs_high_low(10),rel_incs_high_low(11));
      fprintf(fid,'high CF, high sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_high_high(1),rel_incs_high_high(2),rel_incs_high_high(3),rel_incs_high_high(4),rel_incs_high_high(5),rel_incs_high_high(6),rel_incs_high_high(7),rel_incs_high_high(8),rel_incs_high_high(9),rel_incs_high_high(10),rel_incs_high_high(11));
      
      fclose(fid);
      
        end
       
       
      preserve_flags=0;
      clear_flags_watervap
                       
       
return

        iwrite_text_dat=0;
        clear xdat ydat labs

        idat=1;
        xdat(idat).x = xdat_save(1).x;        
        ydat(idat).y = precip_rate_mean;
        labs(idat).l = 'CloudSat Asc desc mean';

        idat=idat+1;
        xdat(idat).x = xdat_save(2).x;
        ydat(idat).y = ydat_save(2).y;
        labs(idat).l = labs_save(2).l;
        
        idat=idat+1;
        xdat(idat).x = xdat_save(3).x;
        ydat(idat).y = ydat_save(3).y;
        labs(idat).l = labs_save(3).l;        

        idat=idat+1;
        xdat(idat).x = xdat_save(1).x;
        ydat(idat).y = precip_rate_asc;
        labs(idat).l = 'CloudSat Ascending';
        
        idat=idat+1;
        xdat(idat).x = xdat_save(1).x;
        ydat(idat).y = precip_rate_desc;
        labs(idat).l = 'CloudSat Descending';        
        
        
        titlenam = ['Precip rates vs longitude at LAT=' num2str(mean(mean_lats))];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        xlab='Longitude';
%        ylab='Precipitation rate (mm hr^{-1})';



        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
        
     
        
    case 127

        
         
%        ioverride_plotglobal_thresh=1;
%        modis_data_plot='Cloud-Sat warm rain precip rate';
        iasc_desc=1; %ascending
%        plot_global_maps;
%        close(gcf);
                       
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data';
        gcm_strs={'CLOUDSAT_PRECIP','CAM5','AM3'};
        use_saved_dat=1;
        var_choose ='Precip rate';
        iscreen_CTT = 0;
        waterVapourMay2005
        close(gcf)

        %store some values
        ydat_save=ydat;
        xdat_save=xdat;
        labs_save = labs;
%        Npoints_gcm_total = Ndat_gcm;
        precip_rate_asc = ydat(1).y;
        
        
%        ioverride_plotglobal_thresh=1;
%        modis_data_plot='Cloud-Sat warm rain precip rate';
        iasc_desc=2; %descending
%        plot_global_maps;
%        close(gcf);
                
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data';
        var_choose ='Precip rate';
        iscreen_CTT = 0;

        waterVapourMay2005
        close(gcf)
        precip_rate_desc = ydat(1).y;

        
%        ioverride_plotglobal_thresh=1;
%        modis_data_plot='Cloud-Sat warm rain precip rate';
        iasc_desc=0; %both
%        plot_global_maps;
%        close(gcf);
        
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data';
        var_choose ='Precip rate';
        iscreen_CTT = 0;
        waterVapourMay2005
        close(gcf)
        precip_rate_mean = ydat(1).y;

        




        iwrite_text_dat=0;
        clear xdat ydat labs

        idat=1;
        xdat(idat).x = xdat_save(1).x;        
        ydat(idat).y = precip_rate_mean;
        labs(idat).l = 'CloudSat Asc desc mean';

        idat=idat+1;
        xdat(idat).x = xdat_save(2).x;
        ydat(idat).y = ydat_save(2).y;
        labs(idat).l = labs_save(2).l;
        
        idat=idat+1;
        xdat(idat).x = xdat_save(3).x;
        ydat(idat).y = ydat_save(3).y;
        labs(idat).l = labs_save(3).l;        

        idat=idat+1;
        xdat(idat).x = xdat_save(1).x;
        ydat(idat).y = precip_rate_asc;
        labs(idat).l = 'CloudSat Ascending';
        
        idat=idat+1;
        xdat(idat).x = xdat_save(1).x;
        ydat(idat).y = precip_rate_desc;
        labs(idat).l = 'CloudSat Descending';        
        
        
        titlenam = ['Precip rates vs longitude at LAT=' num2str(mean(mean_lats))];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        xlab='Longitude';
%        ylab='Precipitation rate (mm hr^{-1})';



        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
        
        
        
    case 126
             
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data';
        var_choose = 'Land Fraction';
        iscreen_CTT = 0;
        %run with no CTT screening to get the total no. points
        waterVapourMay2005
        close(gcf)

        %store the Land fraction and the number of points used
        ydat_save=ydat;
        labs_save = labs;


        
        %run again with CTT screening to get no. warm points
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data'
         var_choose = 'Terrain Height';
        iscreen_CTT = 0;
        %run with no CTT screening to get the total no. points
        waterVapourMay2005
        close(gcf)


%keep xdat and yday - just make text strings for above the lines
iwrite_text_dat=0;
clear text_dat
for idat=1:length(ydat_save)
    for itext=1:length(ydat_save(idat).y)
        str_text=num2str(ydat_save(idat).y(itext),'%1.1f');
        text_dat(idat).text(itext,1:length(str_text)) = str_text;
    end


    if izlim==1
        dz = ( zmax - zmin )/30;
    else
        dz = ( max(ydat(idat).y) - min(ydat(idat).y) ) /30;
    end
    xtext_dat(idat).x = xdat(idat).x;
    ytext_dat(idat).y = ydat(idat).y + dz;

end


ihighlight_points=1;
xdat_highlight(3).x = xdat(1).x;
ydat_highlight(3).y = ydat(1).y;
ydat_highlight(3).y(ydat_save(1).y>0.1)=NaN;

xdat_highlight(4).x = xdat(2).x;
ydat_highlight(4).y = ydat(2).y;
ydat_highlight(4).y(ydat_save(2).y>0.1)=NaN;

labs(1).l=[labs(1).l ' at LAT=' num2str(mean(single_lat_vals(1).val),'%.1f')];
labs(2).l=[labs(2).l ' at LAT=' num2str(mean(single_lat_vals(2).val),'%.1f')];


%  secyA=ydat(1).y
%  secyB=ydat_save(1).y;
% 
%     lab2='Land Fraction';  
%     dual=1;
%     
%     xloc=[1 1 0 0];
%     
%     xdat(3).x=xdat(1).x;
%     ydat(3).y=ydat_save(1).y;
%     labs(3).l='LF';
%     
%     xdat(4).x=xdat(1).x;
%     ydat(4).y=ydat_save(1).y;
%     labs(4).l='LF';
    
     
        
        
        titlenam = ['Terrain Height vs longitude, with land fraction LTE 0.1 marked'];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=1;
        zmin=0;
        zmax=4500;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        xlab='Longitude';
        ylab='Terrain Height (m)';

        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
     
        
     case 125

         
        ioverride_plotglobal_thresh=1;
        modis_data_plot='Cloud-Sat warm rain precip rate';
        plot_global_maps;
        close(gcf);
                       
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data';
        iscreen_CTT = 0;
        %run with no CTT screening to get the total no. points
        waterVapourMay2005
        close(gcf)

        %store the longitude and the number of points used
        xdat_save=xdat;
        labs_save = labs;
%        Npoints_gcm_total = Ndat_gcm;
        precip_rate_warm = ydat(1).y;

        
        ioverride_plotglobal_thresh=1;        
        modis_data_plot='Cloud-Sat all rain precip rate';
        plot_global_maps;
        close(gcf);
        
        %run again with CTT screening to get no. warm points
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data'
        iscreen_CTT = 0;
        %run with no CTT screening to get the total no. points
        waterVapourMay2005
        close(gcf)
        precip_rate_all = ydat(1).y;



        iwrite_text_dat=0;
        clear xdat ydat labs

        idat=1;
        xdat(idat).x = xdat_save(idat).x;
        ydat(idat).y = precip_rate_warm;
        labs(idat).l = 'Warm clouds only';

        xdat(idat+1).x = xdat_save(idat).x;
        ydat(idat+1).y = precip_rate_all;
        labs(idat+1).l = 'All clouds';

        titlenam = ['CloudSat precip rates vs longitude'];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        xlab='Longitude';
        ylab='Precipitation rate (mm hr^{-1})';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
    case 124

        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data'
        iscreen_CTT = 0;
        %run with no CTT screening to get the total no. points
        waterVapourMay2005
        close(gcf)

        %store the longitude and the number of points used
        xdat_save=xdat;
        labs_save = labs;
        Npoints_gcm_total = Ndat_gcm;

        %run again with CTT screening to get no. warm points
        man_choose_water_graph=1;
        graph=115;
        ioverride_longitude_transect_options=1;
        gcm_case = 'Using saved GCM data'
        iscreen_CTT = 1;
        %run with no CTT screening to get the total no. points
        waterVapourMay2005
        close(gcf)
        Npoints_gcm_warm = Ndat_gcm;



        iwrite_text_dat=0;
        clear xdat ydat labs

        idat=2;
        xdat(idat-1).x = xdat_save(idat).x;
        ydat(idat-1).y = (Npoints_gcm_total(idat).N - Npoints_gcm_warm(idat).N) ./ Npoints_gcm_total(idat).N;
        labs(idat-1).l = gcm_strs{idat};

        idat=3;
        xdat(idat-1).x = xdat_save(idat).x;
        ydat(idat-1).y = (Npoints_gcm_total(idat).N - Npoints_gcm_warm(idat).N) ./ Npoints_gcm_total(idat).N;
        labs(idat-1).l = gcm_strs{idat};

        titlenam = ['Cold cloud fraction (of data) vs longitude'];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        xlab='Longitude';
        ylab='Cold cloud fraction of data';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane



        
    case 123
        LAT_val = [72 75]; LON_val = [-3:48]; %Arctic summer region
        
        modis_var_plot = 'Cloud_Fraction_Liquid';
        modis_var_plot = 'Cloud_Effective_Radius_Liquid_Mean';  
        modis_var_plot = 'Cloud_Optical_Thickness_Liquid_Mean';
%        modis_var_plot = 'Solar_Zenith_Minimum';
        
        ilat = find(MLAT_L32007_Arctic>=LAT_val(1) & MLAT_L32007_Arctic<LAT_val(end));
        ilon = find(MLON_L32007_Arctic>=LON_val(1) & MLON_L32007_Arctic<LON_val(end));
        eval_str = [modis_var_plot '_L32007_Arctic.timeseries3(ilat,ilon,:);'];
        dat_modis = eval(eval_str);
        xdat(1).x = dat_modis(:);
        xlab='Actual L3';
     
        
%        d = 171; %day of the year in question
%        day=datenum('01-Jan-2008') + d -1; %in Matlab time
        
%        time = day:1/24:day+1;

%assuming here that the MLAT,MLON is that for the mock L3 and that we load in the
%full global MLAT,MLON
ilat = find(MLAT>=LAT_val(1) & MLAT<LAT_val(end));
ilon = find(MLON>=LON_val(1) & MLON<LON_val(end));
eval_str = [modis_var_plot '_Daily.timeseries3(ilat,ilon,:);'];
dat_modis2 = eval(eval_str);
ydat(1).y = dat_modis2(:);
ylab='Mock L3';

return
        
        labs(1).l = [''];
        
        titlenam = ['Mock L3 vs actual L3 for ' modis_var_plot];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.





        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
 case 122.1
        % load in the .mat file first or set below:-
        
        % --- N.B. - don't forget to set flags in switch filename below for each filename ----
        
        % This is the one for Tau as used in ACPD paper :- 
        filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau'];
        filedir_savevars = '/home/disk/eos1/d.grosvenor/';
        
        %Tau for ACP paper with the top point included as requested by ref#2      
        filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau_for_ACP'];
        filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
        
        %Tau vs gamma_Tau for ACP paper      
        filename = ['saved_vars_for_high_minus_low_sza_vs_gammaTau_Tau_for_ACP'];
        filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';


        filename_savevars = [filedir_savevars filename '.mat'];
        load(filename_savevars);

        %'\tau_{high \theta_0} - \tau_{low \theta_0}' for xlabel -add this
       
        titlenam = ['Difference between low and high SZA vs sctt for ' thresh_str];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        
        %default values that may get overwritten below
        ylab='\sigma_{CTT} (K)';
        itau=0;
        ilow_sens=[1 2 3];

        switch filename
            case 'saved_vars_for_high_minus_low_sza_vs_ctt_Nd' %Nd
                xlab= 'N_{d,high \theta_{0}} - N_{d,low \theta_{0}} (cm^{-3})';
                itau=0;
                ilow_sens=[1 3 5];
            case {'saved_vars_for_high_minus_low_sza_vs_ctt'} %reff
                xlab= 'r_{e,high \theta_{0}} - r_{e,low \theta_{0}} (\mum)';
                itau=0;
                ilow_sens=[1 3 5];
            case 'saved_vars_for_high_minus_low_sza_vs_ctt_Tau'
                xlab= '{\tau}_{high \theta_{0}} - {\tau}_{low \theta_{0}}';
                itau=1;
                ilow_sens=[1 3 5];
            case {'saved_vars_for_high_minus_low_sza_vs_gamma_tau'} %reff
                xlab= 'r_{e,high \theta_{0}} - r_{e,low \theta_{0}} (\mum)';
                itau=0;
                ilow_sens=[1 2 3];  
            case 'saved_vars_for_high_minus_low_sza_vs_ctt_Tau_for_ACP'
                xlab= '\tau_{high \theta_{0}} - \tau_{low \theta_{0}}';
                itau=1;
                ilow_sens=[1 1 1];     
            case 'saved_vars_for_high_minus_low_sza_vs_gammaTau_Tau_for_ACP'
                xlab= '\tau_{high \theta_{0}} - \tau_{low \theta_{0}}';
                itau=1;
                ilow_sens=[1 1 1]; 
                ylab='\gamma_{\tau}';
        end


        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        sens_range = 'low';
%        sens_range = 'high';

switch itau
    case 0
        
        switch sens_range
            case 'low'
        
        idat=idat+1;
        xdat(idat).x=xdat_highSZA(ilow_sens(1)).x - xdat_lowSZA(ilow_sens(1)).x; %
        ydat(idat).y=ydat_highSZA(ilow_sens(1)).y; %
        errordatU(idat).dat = sqrt(errordatU_lowSZA(ilow_sens(1)).dat.^2 + errordatU_highSZA(ilow_sens(1)).dat.^2 );
        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='1.6 {\mu}m, {\theta}=0-41.4';
        
        idat=idat+1;
        xdat(idat).x=xdat_highSZA(ilow_sens(2)).x - xdat_lowSZA(ilow_sens(2)).x; %
        ydat(idat).y=ydat_highSZA(ilow_sens(2)).y; %
       errordatU(idat).dat = sqrt(errordatU_lowSZA(ilow_sens(2)).dat.^2 + errordatU_highSZA(ilow_sens(2)).dat.^2 );
        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='2.1 {\mu}m, {\theta}=0-41.4';
        
        idat=idat+1;
        xdat(idat).x=xdat_highSZA(ilow_sens(3)).x - xdat_lowSZA(ilow_sens(3)).x; %
        ydat(idat).y=ydat_highSZA(ilow_sens(3)).y; %
        errordatU(idat).dat = sqrt(errordatU_lowSZA(ilow_sens(3)).dat.^2 + errordatU_highSZA(ilow_sens(3)).dat.^2 );
        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='3.7 {\mu}m, {\theta}=0-41.4';
        
        
        istyle=1; %the counter for the styles we specify

        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='^'; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
%        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='^'; istyle=istyle+1;


           titlenam = [titlenam ' low sensSZA'];

            case 'high'
                idat=idat+1;
                xdat(idat).x=xdat_highSZA(2).x - xdat_lowSZA(2).x; %
                ydat(idat).y=ydat_highSZA(2).y; %
                errordatU(idat).dat = sqrt(errordatU_lowSZA(2).dat.^2 + errordatU_highSZA(2).dat.^2 );
                errordatL(idat).dat = errordatU(idat).dat;
                labs(idat).l='1.6 {\mu}m, {\theta}>41.4';

                idat=idat+1;
                xdat(idat).x=xdat_highSZA(4).x - xdat_lowSZA(4).x; %
                ydat(idat).y=ydat_highSZA(4).y; %
                errordatU(idat).dat = sqrt(errordatU_lowSZA(4).dat.^2 + errordatU_highSZA(4).dat.^2 );
                errordatL(idat).dat = errordatU(idat).dat;
                labs(idat).l='2.1 {\mu}m, {\theta}>41.4';

                idat=idat+1;
                xdat(idat).x=xdat_highSZA(6).x - xdat_lowSZA(6).x; %
                ydat(idat).y=ydat_highSZA(6).y; %
                errordatU(idat).dat = sqrt(errordatU_lowSZA(6).dat.^2 + errordatU_highSZA(6).dat.^2 );
                errordatL(idat).dat = errordatU(idat).dat;
                labs(idat).l='3.7 {\mu}m, {\theta}>41.4';
        
        
        istyle=1; %the counter for the styles we specify

%        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='^'; istyle=istyle+1;
%        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; istyle=istyle+1;
%        line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
        line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='^'; istyle=istyle+1;


          titlenam = [titlenam ' high sensSZA'];
                
                
        end
        
    case 1
         switch sens_range
            case 'low'
        
        idat=idat+1;
        xdat(idat).x=xdat_highSZA(1).x - xdat_lowSZA(1).x; %
        ydat(idat).y=ydat_highSZA(1).y; %
        errordatU(idat).dat = sqrt(errordatU_lowSZA(1).dat.^2 + errordatU_highSZA(1).dat.^2 );
        errordatL(idat).dat = errordatU(idat).dat;
        labs(idat).l='{\theta}=0-41.4';
        

        
        
        istyle=1; %the counter for the styles we specify

 line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.5 0.5 1]; marker_style(istyle).m='o'; istyle=istyle+1;
% line_pattern(istyle).p= '--'; line_colour(istyle).c=[0.5 0.5 1]; marker_style(istyle).m='^'; istyle=istyle+1;

           titlenam = [titlenam ' low sensSZA'];

            case 'high'
                idat=idat+1;
                xdat(idat).x=xdat_highSZA(2).x - xdat_lowSZA(2).x; %
                ydat(idat).y=ydat_highSZA(2).y; %
                errordatU(idat).dat = sqrt(errordatU_lowSZA(2).dat.^2 + errordatU_highSZA(2).dat.^2 );
                errordatL(idat).dat = errordatU(idat).dat;
                labs(idat).l='{\theta}>41.4';


        
        
        istyle=1; %the counter for the styles we specify

% line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.5 0.5 1]; marker_style(istyle).m='o'; istyle=istyle+1;
 line_pattern(istyle).p= '--'; line_colour(istyle).c=[0.5 0.5 1]; marker_style(istyle).m='^'; istyle=istyle+1;

          titlenam = [titlenam ' high sensSZA'];
                
                
         end
        
         
         
end  %switch itau
       
        ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
        
        ierror_bars='horiz2';
        
        
        


% ---------------------------------------------- case 122 ---------------------------        
        
    case 122  %run multiple runs of PlotTimeHeightVap3 in order to get 3 lines for the means
           %of 2d pdfs (e.g. Nd vs solarZA)
           %calls modisL3_screening_timeseries3  for the screening
           
% --- to plot tau graphs ----
   % set multi_plot_case_122 = 'Nd, just 2.1um vs SZA';
   % and var_plot_122_21mum = 'Optical Depth';
   
% --- to plot re graphs ----
   % multi_plot_case_122 = 'Nd, different wavelengths vs SZA'; 
   % and multi_wavelength_case='Re';
              
        ilow_sens_only=1; %flag to show only low sensZA for the multiple wavelength
                 %cases - e.g. for Nd vs CF plot. Can then also change the
                 %thresh_sensZA_122 thresholds if want all sensZA
        ihigh_sens_only=0; %or set this to show only the high sensZA range
        %It is the second value of thresh_sensZA that determines the
        %boundary
        
        iremove_points_122=1;
        thresh_points_122=10; %originaly set to 10 for ACPD paper
%        thresh_points_122=40;        
        
        ioneRe_only=1;  %set to only plot one Reff value and choose which below
            oneRe = '1.6';
            oneRe = '2.1';
            oneRe = '3.7';
           
      plot_error_bars = 0; %is reset below in certain cases
      thresh_NP=50;
%      thresh_NP=30;                 
      
      minfrac_CF = 0.9; %minimum fraction of the sampled points that had successful cloudy/clear/phase
      %determination (i.e. Npix/Nptot_mockL3 =
      %Cloud_Fraction_Liquid_Pixel_Counts./Cloud_Fraction_Liquid./Total_pixels
      % - restriction (2) as presented in the paper

      minfrac_NpNd = 0.9;        
       %Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3./Cloud_Fraction_Liquid_Pixel_Counts.timeseries3
       %Fraction of points that remain after all previous filtering for
       %which we have an Nd retrieval. Restriction (4) in the SZA paper.
            
      thresh_NP_Nd = 50; %min no. of pixels required for an Nd, re, tau, etc measurement to count
      %(Uses Cloud_Fraction_Liquid_Pixel_Counts2). For this screening
      %usually
      %thresh_NP is the number of pixels that the swath must have covered -
      %i.e. the total number of pixels available. Only a portion of those 
      
      %screening done in -- modisL3_screening_timeseries3  --
      

      
      %for optical depth - make sure that thresh_SZA_range values are
      %strings!
      thresh_SZA_range1 = '[50 55]';
%      thresh_SZA_range1 = '[75 85]';  %this one changes SZA for tau - use the one further below for re/Nd
      
      thresh_SZA_range2 = '[75 85]';
      
%      thresh_SZA_range1 = '[0 90]';
%      thresh_SZA_range2 = '[50 90]';
      
if isstr(thresh_SZA_range1)==0 | isstr(thresh_SZA_range2)==0 
    disp('*** thresh_SZA_range values need to be specified as strings (e.g. = ''[50 55]''***')
    return
end

      N_sza_bins=8;
      
      ichoose_Xbins_set=0; Xbins_set=NaN;
      ichoose_Ybins_set=0; Ybins_set=NaN;
      ichoose_Zbins_set=0; Zbins_set=NaN;
      
    minXbins=-9e99
    maxXbins=-9e99;

    minYbins=-9e99;
    maxYbins=-9e99;
    
    %  For homog factor using meanW for each pixel means the mean of all L2
    %  W values (rather than the mean
% using the L3 mean tau and Re
        y_axis_vals='Mean SZA timeseries3'; ichoose_Ybins_set=0; plotcase_122='normal';
        y_axis_vals='Time of Day'; ichoose_Ybins_set=0; plotcase_122='normal';        
%        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0:0.1:1.0 1.00001]; ichoose_Ybins_set=1; plotcase_122='vs_CF';
%       y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0.8:0.02:1.0 1.00001]; ichoose_Ybins_set=1; plotcase_122='vs_CF';       
%        y_axis_vals='Cloud Fraction from grid vals timeseries3'; Ybins_set=[0.6:0.04:1.0 1.00001]; ichoose_Ybins_set=1;  plotcase_122='vs_CF';              
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:10]; ichoose_Ybins_set=1;  plotcase_122='vs_CF';
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:25]; ichoose_Ybins_set=1; plotcase_122='vs_CF'; 
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:2:25]; ichoose_Ybins_set=1;  plotcase_122='vs_CF';        
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:5:25]; ichoose_Ybins_set=1; plotcase_122='vs_CF';                 
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:1:10 15 20:20:100]; ichoose_Ybins_set=1; plotcase_122='vs_CF';                 
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:2:10 15 20 25 35 55 75]; ichoose_Ybins_set=1; plotcase_122='vs_CF';           
%        y_axis_vals='Homogeneity Parameter timeseries3 using mean W for each pixel'; Ybins_set=[0:2:25 35 50 75 100]; ichoose_Ybins_set=1; plotcase_122='vs_CF';                 
%        Ybins_set=[0:0.1:0.9 0.99 1.00001]; ichoose_Ybins_set=1;
%        y_axis_vals='Mean CTH timeseries3, y-axis'; Ybins_set=[-2:0.5:10]; ichoose_Ybins_set=1; plotcase_122='other';                         
%y_axis_vals='Cloud Top Temp standard deviation, liquid pixels'; Ybins_set=[0:0.4:4]; ichoose_Ybins_set=1; plotcase_122='vs_CF';          
y_axis_vals='Cloud Top Temp standard deviation, liquid pixels'; Ybins_set=[0:0.25:3]; ichoose_Ybins_set=1; plotcase_122='vs_CF';          
%y_axis_vals='Cloud Top Temp standard deviation, liquid pixels'; Ybins_set=[0:0.2:1 2:1:25]; ichoose_Ybins_set=1; plotcase_122='vs_CF';          
% y_axis_vals='Homogeneity Parameter Cahalan Optical Depth (Seethala)'; Ybins_set=[0:0.05:1]; ichoose_Ybins_set=1; plotcase_122='vs_CF';             
% y_axis_vals='Homogeneity Parameter Cahalan Optical Depth (Seethala)'; Ybins_set=[0:0.15:0.75 0.77:0.02:1.01]; ichoose_Ybins_set=1; plotcase_122='vs_CF';  
% y_axis_vals='Homogeneity Parameter Cahalan Optical Depth (Seethala)'; Ybins_set=[0:0.15:0.75 0.77:0.01:1]; ichoose_Ybins_set=1; plotcase_122='vs_CF';   
% y_axis_vals='Homogeneity Parameter Cahalan Optical Depth (Seethala)'; Ybins_set=[0.720001:0.02:1.00001]; ichoose_Ybins_set=1; plotcase_122='vs_CF';  
y_axis_vals='Heterogeneity Parameter Cahalan Optical Depth (Seethala)'; Ybins_set=1-fliplr([0.720001:0.02:1.00001]); ichoose_Ybins_set=1; plotcase_122='vs_CF';  

      if ~exist('ioverride_122') | ioverride_122==0

% ------- choose what to plot here --------------------------------------------------------------          
          
        multi_plot_case_122 = 'Nd, different wavelengths vs SZA'; %also use this for Re now
%        multi_plot_case_122 = 'Re, different wavelengths vs SZA';  
        multi_plot_case_122 = 'Nd, just 2.1um vs SZA'; %also for *** optical depth ***



                
                thresh_cf_122 = [0.1 1.0];
%                thresh_cf_122 = [0.1 0.8];
                thresh_cf_122 = [0.8 1.0]; 
                thresh_cf_122 = [0.99 1.0];                 
%                thresh_cf_122 = [0.85 1.0];                 
%                thresh_cf_122 = [0.6 1.0];                 
%                thresh_cf_122 = [0.0 1.0];

      end
   
                
                thresh_AZ_122 = [50 130];
%                thresh_AZ_122 = [60 120];    
                thresh_AZ_122 = [0 180];     
%                thresh_AZ_122 = [0 80];   
%                thresh_AZ_122 = [100 180];
                thresh_solarZA_122 = [50 55];
%                thresh_solarZA_122 = [75 85];  
%                thresh_solarZA_122 = [50 85];                  
%                thresh_solarZA_122 = [0 90];
%                thresh_solarZA_122 = [0 30];   
%                thresh_solarZA_122 = [50 55];                   

                thresh_CTT_122 = [273-100 373];
%                thresh_CTT_122 = [273-100 273-5];                
                thresh_CTT_122 = [273-5 373]; 
%                thresh_CTT_122 = [273-3 373];   
%                thresh_CTT_122 = [273 373];   
%                thresh_CTT_122 = [273-100 273-3];

                thresh_CTH_122 = [-2 2];
                thresh_CTH_122 = [-20 20];                

%                
                
                %changed this to the homogeneity factor now
                stdW_boundary=1;
                thresh_stdW_122 = [stdW_boundary 1e9];                
                thresh_stdW_122 = [0 stdW_boundary];
                thresh_stdW_122 = [0 1e9]; 
                
                
                %for the different wavelengths case set the higher
                %thresh_sensZA_122 below to be the required boundary if
                %ihigh_sens_only=1. E.g. 
                thresh_sensZA_122 = [0 90]; 
                thresh_sensZA_122 = [0 41.4];  %set so 2nd number is the
%                threshold for sensZA splitting

%                thresh_sensZA_122 = [41.4 90];                  
                
                switch plotcase_122 %this is usully set when we set the y-axis parameter
                    case 'vs_CF'
                        thresh_cf_122 = [0.1 0.8];
                        thresh_cf_122 = [0.99 1.0];    
%                        thresh_cf_122 = [0.9999 1.01];
%                        
%                        thresh_cf_122 = [0.8 1.0];
%                        thresh_cf_122 = [0 1.0]; 
                    case 'normal'
                        thresh_solarZA_122 = [0 85];
                    case 'other'
                end
                
                
                sens_limit = num2str(thresh_sensZA_122(2));
           
        
        switch multi_plot_case_122
%% Tau             
           case 'Nd, just 2.1um vs SZA'
               if ~exist('ioverride_122') | ioverride_122==0
                   var_plot_122_21mum = 'N_d (cm^{-3})';
                   %               var_plot_122_21mum = 'R_{eff} (\mum)';
                   var_plot_122_21mum = 'Optical Depth';
%                   var_plot_122_21mum = 'CTT standard deviation (K)';                   
               end
               
               plot_error_bars = 1;
               
               cf_split = 'no';  %whether to do separate lines for low and high CF or not
                
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';





                   switch var_plot_122_21mum
                       case 'N_d (cm^{-3})'
    

%                xvars={'Nd from grid vals timeseries3','Nd from grid vals timeseries3','Nd from grid vals timeseries3'};
                xvars={'Nd from grid vals timeseries3',...
                    'Nd from grid vals timeseries3'};
                
                 xvars={'Nd from grid vals timeseries3'...
                     };
                
                
                
%                errorvars={'Nd uncertainty from grid vals timeseries3','Nd uncertainty from grid vals timeseries3','Nd uncertainty from grid vals timeseries3'};
                errorvars={'Nd uncertainty from grid vals timeseries3',...
                    'Nd uncertainty from grid vals timeseries3'};
               
                wavelength_str='2.1 \mum, ';
                
    case 'R_{eff} (\mum)'
        %                xvars={'Nd from grid vals timeseries3','Nd from grid vals timeseries3','Nd from grid vals timeseries3'};
%                xvars={'Nd from grid vals timeseries3','Nd from grid vals timeseries3'};
%                errorvars={'Nd uncertainty from grid vals timeseries3','Nd uncertainty from grid vals timeseries3','Nd uncertainty from grid vals timeseries3'};
                errorvars={'Re uncertainty from grid vals timeseries3',...
                           'Re uncertainty from grid vals timeseries3'};
                 
                
                 xvars={'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                        'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                       };
                   
                   wavelength_str='2.1 \mum, ';
                   
    case {'Optical Depth','CTT standard deviation (K)'}
        
        
        %                xvars={'Nd from grid vals timeseries3','Nd from grid vals timeseries3','Nd from grid vals timeseries3'};
%                xvars={'Nd from grid vals timeseries3','Nd from grid vals
%                timeseries3'};
%                errorvars={'Nd uncertainty from grid vals timeseries3','Nd uncertainty from grid vals timeseries3','Nd uncertainty from grid vals timeseries3'};

sens_split = 'yes';  %N.B. for the Reff and Nd case used in the paper set ilow_sens_only=1 instead of this
sens_split = 'no';  %(and change the sensZA limit)

switch sens_split
    case 'yes'
        switch cf_split
            case 'yes'
                switch var_plot_122_21mum
                    case 'Optical Depth'
                        xvars={'Tau reduced dataset Re_1.6 Re_3.7',...
                            'Tau reduced dataset Re_1.6 Re_3.7',...
                            'Tau reduced dataset Re_1.6 Re_3.7',...
                            'Tau reduced dataset Re_1.6 Re_3.7',...
                            };
                    case 'CTT standard deviation (K)'
                        xvars={'Mean CTT timeseries3',...
                            'Cloud Top Temp standard deviation, liquid pixels',...
                            'Cloud Top Temp standard deviation, liquid pixels',...
                            'Cloud Top Temp standard deviation, liquid pixels',...
                            };
                end

            case 'no'
                switch var_plot_122_21mum
                    case 'Optical Depth'
                        xvars={'Tau reduced dataset Re_1.6 Re_3.7',...
                            'Tau reduced dataset Re_1.6 Re_3.7',...
                            };
                    case 'CTT standard deviation (K)'
                        xvars={'Cloud Top Temp standard deviation, liquid pixels',...
                           'Cloud Top Temp standard deviation, liquid pixels',...
                            };

                end
        end

        switch var_plot_122_21mum
            case 'Optical Depth'
                errorvars={'Optical Depth uncertainty from grid vals timeseries3',...
                    'Optical Depth uncertainty from grid vals timeseries3',...
                    'Optical Depth uncertainty from grid vals timeseries3',...
                    'Optical Depth uncertainty from grid vals timeseries3',...
                    };
            case 'CTT standard deviation (K)' %won't bother with measurement uncertainty for this - just 
                % select sampling errors only
%                 errorvars={'Optical Depth uncertainty from grid vals timeseries3',...
%                     'Optical Depth uncertainty from grid vals timeseries3',...
%                     'Optical Depth uncertainty from grid vals timeseries3',...
%                     'Optical Depth uncertainty from grid vals timeseries3',...
%                     };
                

        end

    case 'no'
        switch var_plot_122_21mum
            case 'Optical Depth'
                xvars={'Tau reduced dataset Re_1.6 Re_3.7',...
                    'Tau reduced dataset Re_1.6 Re_3.7',...
                    };

                errorvars={'Optical Depth uncertainty from grid vals timeseries3',...
                    'Optical Depth uncertainty from grid vals timeseries3',...
                    };

            case 'CTT standard deviation (K)'

                xvars={'Cloud Top Temp standard deviation, liquid pixels',...
                    'Cloud Top Temp standard deviation, liquid pixels',...
                    };
% Select sampling errors only.

%                 errorvars={'Optical Depth uncertainty from grid vals timeseries3',...
%                     'Optical Depth uncertainty from grid vals timeseries3',...
%                     };
        end

end

                   

                
                

                   
                   wavelength_str='';
                   
end
                
                %acos(0.5*(1+cos(60*pi/180)))*180/pi = 41.4, so maybe this
                %is the best value (i.e. theta for halfway from cos(0) to
                %cos(60) ).  cos(60)=0.5 so halfway between cos(0)=1 and
                %cos(60) is cos(theta)=0.75, theta=41.4
%                sens_limit = '41.4';
%                sens_limit = '20';                
                
% if ~exist('ioverride_122') | ioverride_122==0
%                 thresh_cf_122 = [0.1 1.0];
% %                thresh_cf_122 = [0.1 0.8];
% %                thresh_cf_122 = [0.8 1.0]; 
% %                thresh_cf_122 = [0.0 1.0];
% 
% end
%                 
%                 thresh_AZ_122 = [50 130];
% %                thresh_AZ_122 = [60 120];    
%                 thresh_AZ_122 = [0 180]; 
%                 
%                 thresh_CTT_122 = [273-100 373];   
%                 thresh_CTT_122 = [273-100 273-5];                   
% %                thresh_CTT_122 = [273-5 373];                   
% %                thresh_CTT_122 = [273-100 273-3];   
%                 
%                 switch plotcase_122
%                     case 'vs_CF'
%                         thresh_solarZA_122 = [50 55];
% %                        thresh_solarZA_122 = [75 85];
%                         thresh_cf_122 = [0.1 1.0];
%                     case 'normal'
%                         thresh_solarZA_122 = [0 85];
%                 end
%                 
%                 stdW_boundary=50;
%                 thresh_stdW_122 = [stdW_boundary 1e9];                
%                 thresh_stdW_122 = [0 stdW_boundary];
%                 thresh_stdW_122 = [0 1e9]; 

thresh_cf_str=['[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']'];
thresh_AZ_str=['[' num2str(thresh_AZ_122(1)) ' ' num2str(thresh_AZ_122(2)) ']'];
thresh_stdW_str=['[' num2str(thresh_stdW_122(1)) ' ' num2str(thresh_stdW_122(2)) ']'];
thresh_SZA_str=['[' num2str(thresh_solarZA_122(1)) ' ' num2str(thresh_solarZA_122(2)) ']'];
thresh_CTT_str=['[' num2str(thresh_CTT_122(1)) ' ' num2str(thresh_CTT_122(2)) ']']; 
thresh_sensZA_str=['[' num2str(thresh_sensZA_122(1)) ' ' num2str(thresh_sensZA_122(2)) ']']; 
thresh_CTH_str=['[' num2str(thresh_CTH_122(1)) ' ' num2str(thresh_CTH_122(2)) ']']; 
                
                
%                screen_strs={'NP + CF mockL3','NP + CF + MEAN sensZA'};
                screen_strs={'NP + CF + MEAN sensZA','NP + CF + MEAN sensZA','NP + CF + MEAN sensZA'};                
                screen_strs={'NP + CF + MEAN sensZA','NP + CF + MEAN sensZA'};
                screen_strs={'NP + CF + MEAN sensZA + MEAN relAZ','NP + CF + MEAN sensZA + MEAN relAZ'};                
                
%                 screen_strs={...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP'...
%                     ,'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP'...
%                     };  
%                 
%                  screen_strs={...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN
%                      solarZA',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
%                     }; 
% % 
%                    screen_strs={...
%                        'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + noice'...
%                        'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + noice'...
%                        'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + noice'...
%                        'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + noice'...
%                        };

 screen_str_single={...
%                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
%                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...                    
                    'NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
%                    'NP4 + CF4 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...                    
%                    'none'
%'NP + CF mockL3'
                    };

                switch sens_split
                    case 'yes'
                        for iscr=1:length(xvars)
                            screen_strs{iscr} = screen_str_single{1};
                        end


                    case 'no'
                        for iscr=1:length(xvars)
                            screen_strs{iscr} = screen_str_single{1};
                        end

                end

                  
                
switch plotcase_122
    
    % ----- vs SZA -------------------------- 
    case 'normal'
%                  thresh_val_strs={...
%                      ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
%                      ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...                                     
%                      ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
%                      ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...                                                         
%                     };
%                 
%                     thresh_val_strs={...
%                      ['thresh_CF=[0.1 0.8]; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
%                      ['thresh_CF=[0.1 0.8]; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...                                     
%                      ['thresh_CF=[0.8 1.0]; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
%                      ['thresh_CF=[0.8 1.0]; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...                                                         
%                     };

switch sens_split
    case 'yes'

        switch cf_split
            case 'yes'
                thresh_val_strs={...
                    ['thresh_CF=[0.1 0.8]; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=[0.1 0.8]; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=[0.8 1.0]; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=[0.8 1.0]; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    };

                labs_multi_specify={....
                    [wavelength_str 'sensZA 0-' sens_limit ', CF 0.1-0.8'],...
                    [wavelength_str 'sensZA GT ' sens_limit ', CF 0.1-0.8'],...
                    [wavelength_str 'sensZA 0-' sens_limit ', CF 0.8-1.0'],...
                    [wavelength_str 'sensZA GT ' sens_limit ', CF 0.8-1.0'],...
                    };
        
            case 'no'                
                thresh_val_strs={...
                    ['thresh_CF=[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    };

                labs_multi_specify={....
                    [wavelength_str 'sensZA 0-' sens_limit],...
                    [wavelength_str 'sensZA GT ' sens_limit],...
                    };
                
        end
        
    case 'no'
        thresh_val_strs={...          
            ['thresh_CF=[0.8 1.0]; thresh_sensZA=' thresh_sensZA_str '; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
            ['thresh_CF=[0.8 1.0]; thresh_sensZA=' thresh_sensZA_str '; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
            };

        labs_multi_specify={....
            [wavelength_str 'sensZA ' num2str(thresh_sensZA_122(1)) '-' num2str(thresh_sensZA_122(2)) ', CF 0.1-0.8'],...
            [wavelength_str 'sensZA ' num2str(thresh_sensZA_122(1)) '-' num2str(thresh_sensZA_122(2)) ', CF 0.8-1.0'],...
            };
end

% ----- vs CF -------------------------- 
    case 'vs_CF'
        switch sens_split
            case 'yes'

                thresh_val_strs={...
                    ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_range1 '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_range1 '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_range2 '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_range2 '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    };

                labs_multi_specify={....
                    [wavelength_str 'sensZA 0-' sens_limit ', SZA ' thresh_SZA_range1],...
                    [wavelength_str 'sensZA GT ' sens_limit ', SZA ' thresh_SZA_range1],...
                    [wavelength_str 'sensZA 0-' sens_limit ', SZA ' thresh_SZA_range2],...
                    [wavelength_str 'sensZA GT ' sens_limit ', SZA ' thresh_SZA_range2],...
                    };


            case 'no'
                thresh_val_strs={...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=' thresh_sensZA_str '; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_range1 '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...                                     
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=' thresh_sensZA_str '; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_range2 '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                    };    
                
                labs_multi_specify={....
                    [wavelength_str 'sensZA ' num2str(thresh_sensZA_122(1)) '-' num2str(thresh_sensZA_122(2)) ', SZA ' thresh_SZA_range1],...
                    [wavelength_str 'sensZA ' num2str(thresh_sensZA_122(1)) '-' num2str(thresh_sensZA_122(2)) ', SZA ' thresh_SZA_range2],...
                    };
                
            case 'use_set_thresh_values'
                
                thresh_val_strs={...
                    ['thresh_CF=[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']; thresh_sensZA=[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...                  
                    };
                
                 labs_multi_specify={....
                    [wavelength_str],...
                    };
                
                
        end
        
        
    end
                                                                                 

                 
                 
                 clear line_pattern line_colour
                 ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
                 istyle=1; %the counter for the styles we specify
               
                 line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.5 0.5 1]; marker_style(istyle).m='o'; istyle=istyle+1;
                 line_pattern(istyle).p= '--'; line_colour(istyle).c=[0.5 0.5 1]; marker_style(istyle).m='^'; istyle=istyle+1;
                 line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.0 0.0 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
                 line_pattern(istyle).p= '--'; line_colour(istyle).c=[0.0 0.0 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
                 
                

                 
                 
%                 xlab_multi = [var_plot_122_21mum ' for CF ' num2str(thresh_cf_122(1),'%1.2f') '-' num2str(thresh_cf_122(2),'%1.2f')];                         
                
 xlab_multi = [var_plot_122_21mum];                                          
                
                
           %Nd or Re - choose below    
           case 'Nd, different wavelengths vs SZA'
                    plot_error_bars = 1;
                    
                    if ~exist('ioverride_122') | ioverride_122==0
                        multi_wavelength_case='Nd';
%                        multi_wavelength_case='Re';
%                        multi_wavelength_case='POLDER';
                    end
                
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
%                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';

switch multi_wavelength_case
    case 'Nd'
                    
                   xvars={...
                      'Nd_{1.6} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'...
                     ,'Nd_{1.6} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'...
                     ,'Nd from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'...
                     ,'Nd from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'...
                     ,'Nd_{3.7} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'...
                     ,'Nd_{3.7} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'...
                     };
                 
                 errorvars={...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     'Nd uncertainty from grid vals timeseries3',...
                     'Nd uncertainty from grid vals timeseries3',...
                     'Nd 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...'Re uncertainty from grid vals timeseries3',...
                     'Nd 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
                     };
                 
    case 'Re'
        
                xvars={...
                    'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    };
                
%N.B. - no instrument uncertainties for 1.6 and 3.7 um wavelengths - using
%the 2.1um ones here
     errorvars={...
            'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
            'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
            'Re uncertainty from grid vals timeseries3',...
            'Re uncertainty from grid vals timeseries3',...
            'Re 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...'Re uncertainty from grid vals timeseries3',...
            'Re 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
            };   
        
         case 'POLDER'
        
             xvars={...
                 'CDR Polder Colocated with MODIS',...
                 };

%             errorvars={...
%                 'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
%                 };
        
        
end
                 
%     cdan(1).c=[1 0 0];  %red
%     cdan(3).c=[0 0.5 0.1]; %dark green
%     cdan(2).c=[0 0 1];  %blue
%     cdan(4).c=[0.5 0.5 0.5];  %grey
%     cdan(5).c=[0 0 0];
%     cdan(6).c=[0.7 0.8 0]; %yellowish
%     cdan(7).c=[1.0 0.7 0.7]; %salmon pink
%     cdan(8).c=[0.7 0.7 1.0]; %light blue
%     cdan(9).c=[0.6 0.6 0.8]; %turqoise?
% 
%         pdan(1).p='-';
%         pdan(2).p='--';
%         pdan(3).p='-.';
    
       clear line_pattern line_colour
       ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
       istyle=1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='^'; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; istyle=istyle+1;
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; istyle=istyle+1;
       line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='^'; istyle=istyle+1;

    % markers(1).m='d';
% markers(2).m='+';
% markers(3).m='o';
% markers(4).m='*';
% markers(5).m='.';
% markers(6).m='x';
% markers(7).m='s';
% markers(8).m='d';
% markers(9).m='^';
% markers(10).m='<';
% markers(11).m='>';
% markers(12).m='p';
% markers(13).m='h';                 
                                 
                
%                screen_strs={'NP + CF mockL3','NP + CF + MEAN sensZA'};
%                screen_strs={'NP + CF + MEAN sensZA','NP + CF + MEAN sensZA','NP + CF + MEAN sensZA'};                
%                screen_strs={'NP + CF + MEAN sensZA','NP + CF + MEAN sensZA'};
%                screen_strs={'NP + CF + MEAN sensZA + MEAN relAZ','NP + CF + MEAN sensZA + MEAN relAZ'};                
%            sens_limit = '41.4';
%            sens_limit = '90';  

 

thresh_cf_str=['[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']'];
thresh_AZ_str=['[' num2str(thresh_AZ_122(1)) ' ' num2str(thresh_AZ_122(2)) ']'];
thresh_stdW_str=['[' num2str(thresh_stdW_122(1)) ' ' num2str(thresh_stdW_122(2)) ']'];
thresh_SZA_str=['[' num2str(thresh_solarZA_122(1)) ' ' num2str(thresh_solarZA_122(2)) ']'];
thresh_CTT_str=['[' num2str(thresh_CTT_122(1)) ' ' num2str(thresh_CTT_122(2)) ']']; 
thresh_CTH_str=['[' num2str(thresh_CTH_122(1)) ' ' num2str(thresh_CTH_122(2)) ']']; 
             
             screen_strs={...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT',...
                    };  
                
                screen_strs={...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT',...
                    };    
                
                screen_strs={...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH',...
                    };
                
                 screen_strs={...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
                    'NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH',...
                    };
                
%                 screen_strs={...
%                     'NP + CF mockL3',...
%                     'NP + CF mockL3',...
%                     'NP + CF mockL3',...
%                     'NP + CF mockL3',...
%                     'NP + CF mockL3',...
%                     'NP + CF mockL3',...
%                     };
                
                
                
%                 screen_strs={...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice',...
%                      'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice',...
%                     };                    
                
             for iscr=1:length(xvars) 
                 screen_strs{iscr}=[...
%                     'NP2 + CF2 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH'...
                     'NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH'...                     
%                     'NP4 + CF4 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH'...                                          
                     ];
             end

                 thresh_val_strs={...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str '; thresh_CTT=' thresh_CTT_str '; thresh_CTH=' thresh_CTH_str ';'],...                                          
                    };
                
                labs_multi_specify={....
                      ['1.6 \mum, sensZA 0-' sens_limit]...
                     ,['1.6 \mum, sensZA GT ' sens_limit]...
                     ,['2.1 \mum, sensZA 0-' sens_limit]...
                     ,['2.1 \mum, sensZA GT ' sens_limit]...
                     ,['3.7 \mum, sensZA 0-' sens_limit]...
                     ,['3.7 \mum, sensZA GT ' sens_limit]...
                     };
                 
                     lp_bk = line_pattern;  
                     lc_bk = line_colour;
                     mk_bk = marker_style;
                     ss_bk = screen_strs;
                     tv_bk = thresh_val_strs;
                     la_bk = labs_multi_specify;
                     xv_bk = xvars;
                     ev_bk = errorvars;
                     
                     switch multi_wavelength_case
                         case 'POLDER'
                         otherwise
                     
                 if ilow_sens_only==1                    
                     clear xvars
                     
                     line_pattern(2).p = lp_bk(3).p; line_pattern(3).p = lp_bk(5).p;
                     line_colour(2).c = lc_bk(3).c; line_colour(3).c = lc_bk(5).c;
                     marker_style(2).m = mk_bk(3).m; marker_style(3).m = mk_bk(5).m;
                     screen_strs{2} = ss_bk{3}; screen_strs{3} = ss_bk{5};
                     thresh_val_strs{2} = tv_bk{3}; thresh_val_strs{3} = tv_bk{5};                     
                     labs_multi_specify{2} = la_bk{3}; labs_multi_specify{3} = la_bk{5};
                     xvars{1} = xv_bk{1}; xvars{2} = xv_bk{3}; xvars{3} = xv_bk{5};
                     errorvars{2} = ev_bk{3}; errorvars{3} = ev_bk{5};   
                     
                 elseif ihigh_sens_only==1                    
                     clear xvars
                     
                     line_pattern(1).p = lp_bk(2).p; line_pattern(2).p = lp_bk(4).p; line_pattern(3).p = lp_bk(6).p;
                     line_colour(1).c = lc_bk(2).c; line_colour(2).c = lc_bk(4).c; line_colour(3).c = lc_bk(6).c;
                     marker_style(1).m = mk_bk(2).m; marker_style(2).m = mk_bk(4).m; marker_style(3).m = mk_bk(6).m;
                     screen_strs{1} = ss_bk{2}; screen_strs{2} = ss_bk{4}; screen_strs{3} = ss_bk{6};
                     thresh_val_strs{1} = tv_bk{2}; thresh_val_strs{2} = tv_bk{4}; thresh_val_strs{3} = tv_bk{6};                     
                     labs_multi_specify{1} = la_bk{2}; labs_multi_specify{2} = la_bk{4}; labs_multi_specify{3} = la_bk{6};
                     xvars{1} = xv_bk{2}; xvars{2} = xv_bk{4}; xvars{3} = xv_bk{5}; xvars{3} = xv_bk{6};
                     errorvars{1} = ev_bk{2}; errorvars{2} = ev_bk{4};  errorvars{3} = ev_bk{6}; 
                     
                 elseif ioneRe_only==1
                     switch oneRe
                         case '1.6'
                             clear xvars

                             line_pattern(1).p = lp_bk(1).p; line_pattern(2).p = lp_bk(2).p;
                             line_colour(1).c = lc_bk(1).c; line_colour(2).c = lc_bk(2).c;
                             marker_style(1).m = mk_bk(1).m; marker_style(2).m = mk_bk(2).m;
                             screen_strs{1} = ss_bk{1}; screen_strs{2} = ss_bk{2};
                             thresh_val_strs{1} = tv_bk{1}; thresh_val_strs{2} = tv_bk{2};
                             labs_multi_specify{1} = la_bk{1}; labs_multi_specify{2} = la_bk{2};
                             xvars{1} = xv_bk{1}; xvars{2} = xv_bk{2};
                             errorvars{1} = ev_bk{1}; errorvars{2} = ev_bk{2};
                             
                         case '2.1'
                             clear xvars

                             line_pattern(1).p = lp_bk(3).p; line_pattern(2).p = lp_bk(4).p;
                             line_colour(1).c = lc_bk(3).c; line_colour(2).c = lc_bk(4).c;
                             marker_style(1).m = mk_bk(3).m; marker_style(2).m = mk_bk(4).m;
                             screen_strs{1} = ss_bk{3}; screen_strs{2} = ss_bk{4};
                             thresh_val_strs{1} = tv_bk{3}; thresh_val_strs{2} = tv_bk{4};
                             labs_multi_specify{1} = la_bk{3}; labs_multi_specify{2} = la_bk{4};
                             xvars{1} = xv_bk{3}; xvars{2} = xv_bk{4};
                             errorvars{1} = ev_bk{3}; errorvars{2} = ev_bk{4};

                         case '3.7'
                             clear xvars

                             line_pattern(1).p = lp_bk(5).p; line_pattern(2).p = lp_bk(6).p;
                             line_colour(1).c = lc_bk(5).c; line_colour(2).c = lc_bk(6).c;
                             marker_style(1).m = mk_bk(5).m; marker_style(2).m = mk_bk(6).m;
                             screen_strs{1} = ss_bk{5}; screen_strs{2} = ss_bk{6};
                             thresh_val_strs{1} = tv_bk{5}; thresh_val_strs{2} = tv_bk{6};
                             labs_multi_specify{1} = la_bk{5}; labs_multi_specify{2} = la_bk{6};
                             xvars{1} = xv_bk{5}; xvars{2} = xv_bk{6};
                             errorvars{1} = ev_bk{5}; errorvars{2} = ev_bk{6};
                             
                             
                     end
                     
                 end
                 
                 
                     end
               
                 switch multi_wavelength_case
                     case 'Nd'
                         xlab_multi = ['N_d (cm^{-3}) for CF ' num2str(thresh_cf_122(1),'%1.2f') '-' num2str(thresh_cf_122(2),'%1.2f')];
                         xlab_multi = ['N_d (cm^{-3})'];
                         
                     case 'Re'
                         xlab_multi = ['R_e (\mum) for CF ' num2str(thresh_cf_122(1),'%1.2f') '-' num2str(thresh_cf_122(2),'%1.2f')];
                         xlab_multi = ['r_e (\mum)'];                         
                     case 'POLDER'
                         xlab_multi = ['POLDER R_e (\mum) for CF ' num2str(thresh_cf_122(1),'%1.2f') '-' num2str(thresh_cf_122(2),'%1.2f')];
                 end



            case 'Re, different wavelengths vs SZA'
                
                re_122_plotcase='normal';
%                re_122_plotcase='CF_var';

                xvars={...
                    'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7',...
                    };
                
%N.B. - no instrument uncertainties for 1.6 and 3.7 um wavelengths - using
%the 2.1um ones here
     errorvars={...
            'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
            'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
            'Re uncertainty from grid vals timeseries3',...
            'Re uncertainty from grid vals timeseries3',...
            'Re 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...'Re uncertainty from grid vals timeseries3',...
            'Re 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)',...
            };


                
                
                plot_error_bars = 1;
                
                 ierror_bars_multi='horiz2';
        %        ierror_bars_multi='vert';
        
                 error_type = 'combined';
                 error_type = 'combined + sampling';  %added in quadrature
%                 error_type = 'bootstrap sampling';
%                 error_type = 'sampling';
%                 error_type = 'mean_percent'; %instrument uncertainty (straight mean of L2 and L3 values)
%                 error_type = 'standard deviation';

                sens_limit = '41.4';
%                sens_limit = '20';           

 if ~exist('ioverride_122') | ioverride_122==0  
     switch re_122_plotcase
         case 'normal'
             thresh_cf_122 = [0.1 1.0];
             thresh_cf_122 = [0.1 0.8];
             thresh_cf_122 = [0.8 1.0];

             thresh_solarZA_122 = [50 55];
             thresh_solarZA_122 = [75 85];
             thresh_solarZA_122 = [0 85];
         case 'CF_var'
             thresh_cf_122 = [0.0 1.0];
             thresh_solarZA_122 = [50 55];
%             thresh_solarZA_122 = [75 85];
     end
 end
                
                thresh_AZ_122 = [50 130];
%                thresh_AZ_122 = [60 120];    
                thresh_AZ_122 = [0 180];     
               
                
                stdW_boundary=50;
                thresh_stdW_122 = [stdW_boundary 1e9];                
                thresh_stdW_122 = [0 stdW_boundary];
                thresh_stdW_122 = [0 1e9]; 

thresh_cf_str=['[' num2str(thresh_cf_122(1)) ' ' num2str(thresh_cf_122(2)) ']'];
thresh_AZ_str=['[' num2str(thresh_AZ_122(1)) ' ' num2str(thresh_AZ_122(2)) ']'];
thresh_stdW_str=['[' num2str(thresh_stdW_122(1)) ' ' num2str(thresh_stdW_122(2)) ']'];
thresh_SZA_str=['[' num2str(thresh_solarZA_122(1)) ' ' num2str(thresh_solarZA_122(2)) ']'];
 
             
             screen_strs={...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
                     'NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA',...
                    };  
                  
   switch re_122_plotcase
         case 'normal'                

                 thresh_val_strs={...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...                     
                    };
                
                labs_multi_specify={....
                      ['1.6 \mum, sensZA 0-' sens_limit]...
                     ,['1.6 \mum, sensZA GT ' sens_limit]...
                     ,['2.1 \mum, sensZA 0-' sens_limit]...
                     ,['2.1 \mum, sensZA GT ' sens_limit]...
                     ,['3.7 \mum, sensZA 0-' sens_limit]...
                     ,['3.7 \mum, sensZA GT ' sens_limit]...
                     };
                 
       case 'CF_var'
           thresh_val_strs={...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[0 ' sens_limit ']; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...
                     ['thresh_CF=' thresh_cf_str '; thresh_sensZA=[' sens_limit ' 90]; thresh_relAZ = ' thresh_AZ_str '; thresh_stdW = ' thresh_stdW_str '; thresh_SZA = ' thresh_SZA_str ';'],...                     
                    };
                
                labs_multi_specify={....
                      ['1.6 \mum, sensZA 0-' sens_limit]...
                     ,['1.6 \mum, sensZA GT ' sens_limit]...
                     ,['2.1 \mum, sensZA 0-' sens_limit]...
                     ,['2.1 \mum, sensZA GT ' sens_limit]...
                     ,['3.7 \mum, sensZA 0-' sens_limit]...
                     ,['3.7 \mum, sensZA GT ' sens_limit]...
                     };
           
           
   end
                 
                 
                 clear line_pattern line_colour marker_style
                 ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
                 istyle=1;
                 line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
                 line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
                 line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
                 line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; istyle=istyle+1;
                 line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
                 line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; istyle=istyle+1;
                    
                
                    
                xlab_multi = ['Re (\mum) for CF ' num2str(thresh_cf_122(1),'%1.2f') '-' num2str(thresh_cf_122(2),'%1.2f')];                             
                
                
                    
        end
                                                       
                        % 0.8<CF<=1.000001
                        %also thresh_sensZA(1)<sensZA<=thresh_sensZA(2)
             
           
                
%                labs_multi_specify=xvars;
                
%                y_axis_vals='Mean SZA timeseries3';
                
 
        
        
        clear allSZA_pdfs highSZA_pdf xdat_multi ydat_multi labs_multi errordatU_multi errordatL_multi sampling_errorU sampling_errorL
% -------------------------------------------------------------------------       
%    loop through the multiple runs of plotTimeHeightVap3        
% -------------------------------------------------------------------------

        for idat_multi=1:length(xvars)
        
            ioverride_pdf_varchoose=1; %override the defaults
            ioverride_pdf=1; %override the defaults  
            
            %flags to allow the direct specification of Xbins, Ybins, Zbins
            ichoose_Xbins=ichoose_Xbins_set; Xbins=Xbins_set; 
            ichoose_Ybins=ichoose_Ybins_set; Ybins=Ybins_set; 
            ichoose_Zbins=ichoose_Zbins_set; Zbins=Zbins_set; 
            
%             i577 = 'MODIS_plot_UW';
%             x_axis_vals = xvars{idat_multi};
%             y_axis_vals_save = y_axis_vals;
%             y_axis_vals = 'Mean SZA timeseries3';
%             screen_type = screen_strs{idat_multi};
%             eval(thresh_val_strs{idat_multi});
%             man_choose_plotTimeHeight_graph=1;
%             logflag=0;
%             dlogflag=0;
%             noplot=1; %flag to say to just do the calculations and not to plot
%             nXpdf=500;
%             nYpdf=12;
%             ndims_hist=2;
%             plotTimeHeightVap3
            
            
            
            i577 = 'MODIS_plot_UW';
            x_axis_vals = xvars{idat_multi};
            screen_type = screen_strs{idat_multi};
            eval(thresh_val_strs{idat_multi});
            man_choose_plotTimeHeight_graph=1;
            logflag=0;
            dlogflag=0;
            noplot=1; %flag to say to just do the calculations and not to plot
            nXpdf=500;
            nYpdf=N_sza_bins;
            ndims_hist=2;
            plotTimeHeightVap3
%            close(gcf);

%sampling error = std_dev / sqrt(N)
            sampling_errorU(idat_multi).dat = std_dev_X ./ sqrt(NX_vals);
            sampling_errorL(idat_multi).dat = std_dev_X ./ sqrt(NX_vals);
            
            stdev_multi(idat_multi).dat = std_dev_X;
            
            %export the final SZA bin PDF by expanding the histogram into
            %its component values (will be some resolution loss due to
            %the values being put into bins - using bin mid-points)
            highSZA_pdf(idat_multi).dat = expand_PDF(mid_Xbins,qh(end-1,1:end-1));
            
            for isza=1:size(qh,1)-1
                allSZA_pdfs{isza}(idat_multi).dat = expand_PDF(mid_Xbins,qh(isza,1:end-1));                
            end
            
           
            
            man_choose_water_graph=1;
            graph=96; %mean of 2D PDF
            noplot=1;
            waterVapourMay2005
%            close(gcf);

            xdat_multi(idat_multi).x = xdat(1).x;
            ydat_multi(idat_multi).y = ydat(1).y';
%            labs_multi(idat_multi).l = xvars{idat_multi};               
            labs_multi(idat_multi).l = labs_multi_specify{idat_multi};   


            if iremove_points_122==1
                iremove_122=find(NX_vals<thresh_points_122);
                xdat_multi(idat_multi).x(iremove_122)=NaN;
                fprintf('\n*** Removing points thresh_points_122=%d\n',thresh_points_122);
            end
            
            
            if plot_error_bars==1
                switch error_type
                    case 'sampling'
                        errordatU_multi(idat_multi).dat = sampling_errorU(idat_multi).dat;
                        errordatL_multi(idat_multi).dat = sampling_errorL(idat_multi).dat;
                        
                    case 'standard deviation'
                        errordatU_multi(idat_multi).dat = stdev_multi(idat_multi).dat;
                        errordatL_multi(idat_multi).dat = stdev_multi(idat_multi).dat;                        
                        
                    case 'bootstrap sampling'
                        ioverride_bootstrap=1;
                        xvars_boot = xvars{idat_multi};  
                        bootstrap_from_2d_PDF
                        errordatU_multi(idat_multi).dat = err_bootstrap;
                        errordatL_multi(idat_multi).dat = err_bootstrap;

                    case {'combined','mean_percent','combined + sampling'}
                        %for the errors we do a 2D pdf of the absolute uncertainty
                        %(x-axis) vs SZA (y). This calculates the sqrt of the  sum of the squares
                        %of all the uncertainties divided by N, which should be the
                        %combined error (if the errrors are assumed random)

                        ioverride_pdf_varchoose=1; %override the defaults
                        ioverride_pdf=1;
                        
                        %flags to allow the direct specification of Xbins, Ybins, Zbins
                        ichoose_Xbins=ichoose_Xbins_set; Xbins=Xbins_set; 
                        ichoose_Ybins=ichoose_Ybins_set; Ybins=Ybins_set;
                        ichoose_Zbins=ichoose_Zbins_set; Zbins=Zbins_set;
                        
                        x_axis_vals = errorvars{idat_multi};
                        screen_type = screen_strs{idat_multi};
                        eval(thresh_val_strs{idat_multi});
                        nXpdf=500;
                        nYpdf=N_sza_bins;
                        man_choose_plotTimeHeight_graph=1;
                        logflag=0;
                        dlogflag=0;
                        noplot=1; %flag to say to just do the calculations and not to plot
                        ndims_hist=2;
                        plotTimeHeightVap3
                        %                close(gcf);

                        %this is the mean of the sqrt of the sum of the square of
                        %the error
                        
                        switch error_type
                            case 'combined + sampling'
                                 errordatU_multi(idat_multi).dat = sqrt ( rms_X.^2 + sampling_errorU(idat_multi).dat.^2 );
                                 errordatL_multi(idat_multi).dat = sqrt ( rms_X.^2 + sampling_errorU(idat_multi).dat.^2 );
                        
                            case 'combined'
                                errordatU_multi(idat_multi).dat = rms_X; %from PlotTimeHeightVap3
                                errordatL_multi(idat_multi).dat = rms_X;
                            case 'mean_percent'
                                errordatU_multi(idat_multi).dat = X_mean;
                                errordatL_multi(idat_multi).dat = X_mean;   
                        end

                end
                  
                  
              titlenam = [titlenam ', ' remove_character(error_type,'_','-') ' error bars'];

            end
        
        end    
        
        iwrite_text_dat=0;
        noplot=0;
                                  
        xdat = xdat_multi;
        ydat = ydat_multi;
        labs = labs_multi;
        xlab = xlab_multi;
        
        if plot_error_bars==1
            errordatU = errordatU_multi;
            errordatL = errordatL_multi;
            ierror_bars = ierror_bars_multi;
        end

        
%        titlenam = ['Nd vs SZA for lat=' num2str(MLAT(ilat)) ' and lon=' num2str(MLON(ilon))];
aqua_terra_str='';
for istr=1:length(aqua_terra_timeseries3)
    aqua_terra_str=[aqua_terra_str aqua_terra_timeseries3{istr} ' '];
end
    
        titlenam = [titlenam ' ' aqua_terra_str];

        figname=titlenam;
        savename=figname;

        xlims=1;
%        xlimits=1000*[0 0.025]; %set in case number 96

        izlim=1;
%        zmin=1500;
%        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

%        ylab='SZA'; %set in case number 96
%        xlab='UTC Time';



        lor=0; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

%        posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];
        
        clear rel_incs
        for idat_inc=1:length(xdat_multi)
        rel_incs(idat_inc) = (xdat_multi(idat_inc).x(end)/xdat_multi(idat_inc).x(1) - 1)*100;
        end
     
            
    case 121
        LAT_find = 90;
        ilat=findheight_nearest(MLAT,LAT_find);
        LON_find = 0;
        ilon=findheight_nearest(MLON,LON_find)
        
        d = 171; %day of the year in question
        day=datenum('01-Jan-2008') + d -1; %in Matlab time
        
        time = day:1/24:day+1;
        
        xdat(1).x = (time - time(1)) *24;
        ydat(1).y = sun_pos(time,MLAT(ilat),MLON(ilon));
        labs(1).l = ['day ' num2str(d)];
        
        titlenam = ['SZA vs time for lat=' num2str(MLAT(ilat)) ' and lon=' num2str(MLON(ilon))];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='SZA';
        xlab='UTC Time';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane



        
    case 120
        LAT_find = 71;
        ilat=findheight_nearest(MLAT,LAT_find);
        d = 171; %day of the year in question
        day=datenum('01-Jan-2008') + d -1; %in Matlab time
                                                      
        titlenam = ['Orbit properties for day=' num2str(d) ' at lat=' num2str(MLAT(ilat))];

        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[-60 60];

        izlim=1;
        zmin=0;
        zmax=80;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Degrees';
        xlab='Longitude';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;        
       
        idat=idat+1;
        xdat(idat).x=MLON; %
        ydat(idat).y=Solar_Zenith_Minimum.timeseries3(ilat,:,d); %
        labs(idat).l='Min SZA';
        
        idat=idat+1;
        xdat(idat).x=MLON; %
        ydat(idat).y=Solar_Zenith_Maximum.timeseries3(ilat,:,d); %
        labs(idat).l='Max SZA';
        
        idat=idat+1;
        xdat(idat).x=MLON; %
        ydat(idat).y=Sensor_Zenith_Minimum.timeseries3(ilat,:,d); %
        labs(idat).l='Min Sensor';
        
        idat=idat+1;
        xdat(idat).x=MLON; %
        ydat(idat).y=Sensor_Zenith_Maximum.timeseries3(ilat,:,d); %
        labs(idat).l='Max Sensor';
        
    case 119
        LAT_find = 51;
        ilat=findheight_nearest(MLAT,LAT_find)
        d = 171; %day of the year in question
        day=datenum('01-Jan-2008') + d -1; %in Matlab time
        
        lon_range = [5 25];
        [ilons1 ilons2] = findheight_nearest(MLON,lon_range(1),lon_range(2));
        ilons=ilons1:ilons2;
        
        clear local_time
        for ilon2=1:length(ilons)
            ilon2
            ilon = ilons(ilon2);
            %the target SZA (min from L3)
            sza_target=Solar_Zenith_Minimum.timeseries3(ilat,ilon,d);
            %finds the UTC time for this SZA
            [time_sza,sza2]=get_time_for_SZA2(MLAT(ilat),MLON(ilon),sza_target,day);
            %convert to local (solar) time based in hours on longitude
            local_time(ilon) = (time_sza + MLON(ilon)/360 - day)*24;            
        end
        
          titlenam = ['Mean No. overpasses'];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Solar time of min SZA';
        xlab='Longitude';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
       
        idat=idat+1;
        xdat(idat).x=MLON; %
        ydat(idat).y=local_time; %
        labs(idat).l='';
        
    case 118
        swath=2330; %km
        R=6378.140;

        lat_swaths = [0:1:80];

        Nswaths = 16*swath ./ (2*pi*R*cos(lat_swaths*pi/180));

        titlenam = ['Mean No. overpasses'];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='No. overpasses';
        xlab= 'Latitude';



        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
       
        idat=idat+1;
        xdat(idat).x=lat_swaths; %
        ydat(idat).y=Nswaths; %
        labs(idat).l='';
        
case 117
        % gcm profiles
        
        
       
        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=1;
        zmin=0;
        zmax=3;
        
        itime=1;
        
        LATS = [-30 -26 -20];
        LONS = [-77.5 -75 -75];
        
        for ilats=1:length(LATS)

            [ilat,ilon] = getind_latlon_quick(gcm_Plat2D_edges',gcm_Plon2D_edges',LATS(ilats),LONS(ilats),0.1);
            lon2=gcm_Plon2D(ilon,ilat); %centre of the box
            lat2=gcm_Plat2D(ilon,ilat);


            nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

            ylab='Height (km)';

            xlab = 'Grid box mean LWC (g m^{-3})';
            xlab = 'Cloud Fraction';
            xlab = 'Temperature (K)';

            switch xlab
                case 'Temperature (K)'
                    xdat(ilats).x = gcm_temp(itime,:,ilon,ilat);
                    labs(ilats).l = ['lat=' num2str(lat2,'%.1f') ',lon=' num2str(lon2,'%.1f')];


                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                    
                case 'Grid box mean LWC (g m^{-3})'
                    xdat(ilats).x = 1e3*gcm_lwc_av(itime,:,ilon,ilat);
                    labs(ilats).l = ['lat=' num2str(lat2,'%.1f') ',lon=' num2str(lon2,'%.1f')];


                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                case 'Cloud Fraction'
                    xdat(ilats).x = gcm_cf(itime,:,ilon,ilat);
                    labs(ilats).l = ['lat=' num2str(lat2,'%.1f') ',lon=' num2str(lon2,'%.1f')];
                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                    xlims=1;
                    xlimits=[0 1];
            end

            ydat(ilats).y = 1e-3*h_full(1,:,ilon,ilat);

        end
        

        

        
        
        titlenam = ['gcm profile of ' xlab ' for time=' num2str(itime) ' ' gcm_str time_UTC_str thresh_str];  
%        titlenam = [titlenam gcm_str time_UTC_str thresh_str]; 
        
        figname=titlenam;
        savename=figname;
        




        
        

        
 case 116
     % search for line_pattern for the specification of linestyles for each
     % model/obs
     
        % longitude transects
        
%        var_plot='Cloud Depth (m)';
%        var_plot='LWP (g m^{-2})';
%        var_plot='low_CF';
%        var_plot='Nd (cm^{-3})';    
%        var_plot='Precip rate (mm hr^{-1})';    

var_plot = var_load;

        add_modis_reduced=0;

        
        switch var_plot
             case {'Reff_COSP','Reff_COSP_allCF','REFFL_maxlayer','REFFL_max_noCF','REFFL_maxliq'}
                unit_str = '(\mum)';
                %                ydat(1).y=Nd_lon.modis;
%                ydat(2).y=Nd_lon.cam5;
%                ydat(3).y=Nd_lon.am3;
                
                izlim=1;
                zmin=5;
                zmax=25;
            
            case 'Precip_rate'
                unit_str = '(mm hr^{-1})';
                %                ydat(1).y=Nd_lon.modis;
%                ydat(2).y=Nd_lon.cam5;
%                ydat(3).y=Nd_lon.am3;
                
                izlim=1;
                zmin=0;
                zmax=0.1;
                
               
                
            case 'Cloud Depth'
                unit_str = '(m)';
                
                izlim=1;
                zmin=150;
                zmax=800;
                

                
            case 'LWP'

                unit_str = '(g m^{-2})';

                izlim=1;
                zmin=0;
                zmax=150;
                
                
                izlim=1;
                zmin=0;
                zmax=250;
                
              case 'LWP2'

                unit_str = 'all CF (g m^{-2})';

                izlim=1;
                zmin=0;
                zmax=150;   
                
            case {'LTS','LTS1000'}

                unit_str = '(K)';

                izlim=1;
                zmin=12;
                zmax=30;     
                
            case 'qv700'

                unit_str = '(g kg^{-1})';

                izlim=1;
                zmin=0;
                zmax=8;                     



            case 'Nd'

                unit_str = '(cm^{-3})';

                izlim=1;
                zmin=0;
                zmax=300;
                
                %add on the in-situ datapoints from Rob
                switch lon_load
                    case '20S';
                        ihighlight_points=1;
                        xdat_highlight(1).x = [-82.5 -77.5 -72.5];
                        ydat_highlight(1).y = [80 145 195];
                end
                
%                 if ihighlight_points==1
%                     for idat=1:length(xdat_highlight)
%                         plot(xdat_highlight(idat).x,ydat_highlight(idat).y,'bo','markersize',16);
%                     end
%                 end
                
                
                

                
            case 'low_CF'
                unit_str = '';
                
                izlim=1;
                zmin=0;
                zmax=1;
                
            otherwise
                disp('*** NEED to set attributes in case 116 of watervap ***');
                unit_str = '';
                
                izlim=0;
                zmin=0;
                zmax=0;
                                              
        end
        
        ylab=[var_plot ' ' unit_str];  
        
        ylab = remove_character(ylab,'_', ' ');
        ichoose_styles=1;
        istyle=1;        

        
        %flag to say whether is nighttime or not - don't plot MODIS if so
        switch time_of_day_multi{itime_multi}
            case 'nighttime'
                inight=1;
                add_modis_reduced=0;
            otherwise
                inight=0;
                %leave add_modis_reduced as set above
        end
        
        
        
        switch add_modis_reduced
            case 1
                idat=1;
                xdat(idat).x = lon_lon(idat).x;
                ydat(idat).y = dat_lon(idat).y*0.8;
                labs(idat).l = remove_character([gcm_load_multi{idat} ' ' am3_dataset_str{idat} ' minus 20%'],'_',' ');
                idat_start=2;
                
                line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;

            otherwise
                idat_start=1;
        end
        
        time_str_116_2='';
        idat2 = idat_start;
        for idat3=idat_start:length(lon_lon)+idat_start-1
            %idat2 is incremented at the end of the loop
            idat = idat3 - idat_start + 1;
            
            labs(idat2).l = remove_character([gcm_load_multi{idat} ' ' am3_dataset_str{idat}],'_',' ');
            
            if length(strfind(labs(idat2).l,'MODIS'))==0 | length(strfind(labs(idat2).l,'POLDER'))==0 | inight==0 
            
            xdat(idat2).x = lon_lon(idat).x;
            ydat(idat2).y = dat_lon(idat).y;
            
            
            
%            time_str_116_2 = cat(2,time_str_116_2,' ',time_str_116{idat});
            

switch labs(idat2).l
    case {'ERAInt  ','MODIS  ','CALIPSO monthly nighttime ','CALIPSO monthly NIGHTTIME ','CALIPSO monthly DAYTIME ','CALIPSO monthly average ',' CLOUDSAT PRECIP  ','ascending CLOUDSAT PRECIP  ','descending CLOUDSAT PRECIP  ','AMSRE daytime ','AMSRE nighttime ','AMSRE average ','POLDER  '}
        line_pattern(istyle).p= '-.';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        if exist('man_choose_water_graph') & man_choose_water_graph==1
            switch labs(idat).l
                case {'ascending CLOUDSAT PRECIP  ',' CLOUDSAT PRECIP  ','descending CLOUDSAT PRECIP  '}
                    labs(idat2).l = 'CloudSat';                    
                case {'AMSRE daytime ','AMSRE nighttime ','AMSRE average '}
                    labs(idat2).l = 'AMSRE';
            end
        end
        
        switch gcm_load_multi{2}
            case 'AMSRE_daytime'
           switch  labs(idat2).l
               case 'MODIS  '
                   istyle=istyle-1;
                    line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
           end
        end

    case {'CAM5  CAM5 1deg','CAM5 COSP  CAM5 1deg'}
        labs(idat2).l = 'CAM5 1deg';
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
    case {'CAM5 CLUBB  CAMCLUBB','CAM5 CLUBB COSP  '}
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'CAM5 CLUBB';
     case {'CAM5 CLUBBv2 COSP  '}
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'CAM5 CLUBBv2';    
    case 'CAM5  CAM5 2deg'
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.6 0.6 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'CAM5 2deg';
    case 'AM3  2deg'
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'AM3 2deg';
    case 'AM3  0pt5deg'
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'AM3 0.5deg';
    case 'AM3 CLUBB  AM3CLUBB'
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'AM3 CLUBB';    
    case 'AM3 CLUBBv2  '
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
        labs(idat2).l = 'AM3 CLUBBv2';            
        
end


              idat2=idat2+1; %increment idat if have added data
            end



        end
        

        if length(time_str_116_2)>1
            time_str_116_2 = time_str_116{2};
        else
            time_str_116_2 = time_str_116{1};
        end
        
        
%         labs(1).l='MODIS 2000-2010';
%         labs(2).l='AM3 2deg 2007-2010';
%         labs(3).l='AM3 0.5deg 2007-2008';
%         labs(4).l='CAM5 2deg 2001';   
%         
%         labs(1).l='MODIS 2007-2010';
%         labs(2).l='CAM5 CLUBB';
%         labs(3).l='AM3 0.5deg 2007-2008';
%         labs(4).l='CAM5 2deg 2001';   
        
        
%            xdat(1).x = lon_lon.modis;        
%            xdat(2).x = lon_lon.cam5;
%            xdat(3).x = lon_lon.am3;
            
%            labs(1).l = 'MODIS';
%            labs(2).l = 'CAM5';
%            labs(3).l = 'AM3';


% **** remove the data that is over land *****

switch lon_load
    case '10S'
        lon_cutoff = -82; %land is at approx 73W at this latitude
    case '20S'        
        lon_cutoff = -73; %land is at approx 73W at this latitude
    case '30S'        
        lon_cutoff = -75; %land is at approx 73W at this latitude
    case '10to30S'
        lon_cutoff = -73; %land is at approx 73W at this latitude   
        lon_cutoff = -80; %land is at approx 73W at this latitude           
        
end

for idat=1:length(xdat)
    icut = find(xdat(idat).x>lon_cutoff);
    xdat(idat).x(icut)=NaN; 
end


                  
        titlenam = [var_plot ' Longitude Transect at ' lon_load ' for ' time_str_116_2];
        if exist('man_choose_water_graph') & man_choose_water_graph==1
            titlenam = [titlenam '_' upper(time_of_day)];
        end
        figname=titlenam;
        savename=figname;

%        titlenam = ['Longitude Transect at ' lon_load];
        

        xlims=0;
        xlimits=1000*[0 0.025];
        
       

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.


        xlab= 'Longitude';



        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

       
        
    case 115 %time_inds_modisL3_timeseries3
        

         slice_2D = 0;
        
%         if icf_low_only==1
%             pstr = ['.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'];
%         else
%             pstr='';
%         end

         pstr='';



        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
%        y_axis_type = 'log10_matlab';
        %LON_val = [-76.5 -72]; LON_val=[-81 -76.5]; LON_val = [-85.5-81];
        %LON_val = [-90 -85.5];

if ~exist('ioverride_lat_vals') |  ioverride_lat_vals ==0
        
        %old lon bins
        LAT_val = [-22 -18; -22 -18; -22 -18; -22 -18]; LON_val = [-76.5 -72; -81 -76.5; -85.5 -81; -90 -85.5];   
        %VOCA ones - plus Klein and Hartmann
        LAT_val = [-22 -18; -22 -18; -22 -18; -20 -10]; LON_val = [-75 -70; -80 -75; -86 -80; -90 -80];
        LAT_val = [-22.74 -18;]; LON_val = [-120 -71.25];
%        LAT_val = [-22.74 -18;]; LON_val = [-120 -56.25];        
%        LAT_val = [-22.74 -18;]; LON_val = [-90 -45.25];        
%        LAT_val = [-22.74 -18]; LON_val = [-76.25 -71.25];
%        LAT_val = [-22.74 -18;]; LON_val = [-104 -71.25]; %for precip - where we start to get cold rain contributions
        %in cloudsat data
%        LAT_val = [-27 -25]; LON_val = [-120 -71.25];  

LAT_val = [-32.74 -28;]; LON_val = [-120 -71.25];
%LAT_val = [-22.74 -18;]; LON_val = [-120 -71.25];
%LAT_val = [-12.74 -8;]; LON_val = [-120 -71.25];

LAT_val = [-32.74 -28;]; LON_val = [-104 -71.25];
LAT_val = [-22.74 -18;]; LON_val = [-104 -71.25];
%LAT_val = [-12.74 -8;]; LON_val = [-104 -71.25];
LAT_val = [-26 -14;]; LON_val = [-120 -71.25];
LAT_val = [-26 1;]; LON_val = [-120 -71.25];

end

%time_inds_modisL3_timeseries3 

time_UTC_str='';
        
                

        irestrict_times=1;
        
        
            ione_lat=0; %flag to say that we only choose one lat index (i.e. only one grid point,
            %and not an average of a few lat indices)
                        
            use_saved_dat=1;
            
            iscreen_CTT = 0;
            iwrite_text_dat=0;
        
        if ~exist('ioverride_longitude_transect_options') | ioverride_longitude_transect_options==0


            
            
            gcm_case = 'CTP';
            gcm_case = 'H';
            gcm_case = 'LWP';
            gcm_case = 'LWP2';            
            gcm_case = 'Nd';
%            gcm_case = 'Pressure_vert_slice';
%            gcm_case = 'LWC_vert_slice (g kg^{-1})';
            gcm_case = 'LTS';
            
 gcm_case = 'Reff_COSP_allCF';
% gcm_case = 'REFFL_maxlayer';
% gcm_case = 'REFFL_maxliq';
% gcm_case = 'REFFL_max_noCF';


%
%            gcm_case = 'Terrain Height';
%            gcm_case = 'Land Fraction';
%            gcm_case = 'CALIPSO Cloud Fraction from CFAD';
%            gcm_case = 'CloudSat dbZ from CFAD';    
%            gcm_case = 'Reff_COSP';   
%            gcm_case = 'Tau_COSP';  
%            gcm_case = 'Nd_COSP';              

%             gcm_case = 'Using saved GCM data';
            
%            var_choose ='Precip rate'; %chose this for CloudSat too
%            var_choose = 'Terrain Height';
%            var_choose = 'Land Fraction';
%            var_choose =  'low CF';
%            var_choose =  'high CF';        
      

            gcm_strs={'MODIS','CAM5','AM3'};
%            gcm_strs={'CALIPSO','CAM5','AM3'};            
            gcm_strs={'CLOUDSAT_PRECIP','CAM5','AM3'};
%            gcm_strs={'CAM5','AM3'};
            gcm_strs={'AM3'};
%            gcm_strs={'CALIPSO'};  
            gcm_strs={'MODIS'};
%            gcm_strs={'CAM5'};
%            gcm_strs={'CAM5_COSP'};
%            gcm_strs={'CAM5_CLUBB'};
%            gcm_strs={'CALIPSO'}; 
%            gcm_strs={'CALIPSO_monthly'};      
%            gcm_strs={'CLOUDSAT_PRECIP'};
%            gcm_strs={'AMSRE'};

                 gcm_time_of_day_select=1; %set to one as a default
            
        end
        
        switch gcm_case
            case 'Using saved GCM data'
                titlenam = ['gcm longitude slice of ' var_choose ' for '];
            otherwise
                titlenam = ['gcm longitude slice of ' gcm_case ' for '];

        end
        
        switch gcm_case
            case {'CALIPSO Cloud Fraction from CFAD','CloudSat dbZ from CFAD'}
                icfad=1;
            otherwise
                icfad=0;
        end

 
        xlab = 'Longitude';

            


        clear mean_val Ndat_gcm single_lat_vals mean_lats mean_LAT_val
        Ndat_gcm=[];
        thresh_str='';
        
        



for idat=1:length(gcm_strs)
    gcm_str = gcm_strs{idat};
         thresh_str2=gcm_strs{idat};
         
         switch gcm_strs{idat}
             case{'POLDER'}
                 Plat2D_edges = eval(['gcm_Plat2D_edges_POLDER;']);
                 Plon2D_edges = eval(['gcm_Plon2D_edges_POLDER']);
                 Plat2D = eval(['gcm_Plat2D_POLDER']);
                 Plon2D = eval(['gcm_Plon2D_POLDER']);      
                 gcm_time_of_day_select=0;
                 gcm_time_UTC = eval(['gcm_time_UTC_POLDER']);
                 time_series_type='POLDER'; %default type of time data
             
             
             case{'AMSRE'}
                 Plat2D_edges = eval(['gcm_Plat2D_edges_AMSRE;']);
                 Plon2D_edges = eval(['gcm_Plon2D_edges_AMSRE']);
                 Plat2D = eval(['gcm_Plat2D_AMSRE']);
                 Plon2D = eval(['gcm_Plon2D_AMSRE']);      
                 gcm_time_of_day_select=0;
                 gcm_time_UTC = eval(['gcm_time_UTC_AMSRE']);
                 time_series_type='AMSRE'; %default type of time data
                 
             case {'CALIPSO','CALIPSO_monthly'}
                 Plat2D_edges = eval(['gcm_Plat2D_edges_CALIPSO_monthly;']);
                 Plon2D_edges = eval(['gcm_Plon2D_edges_CALIPSO_monthly']);
                 Plat2D = eval(['gcm_Plat2D_CALIPSO_monthly']);
                 Plon2D = eval(['gcm_Plon2D_CALIPSO_monthly']);      
                 gcm_time_of_day_select=0;
                 gcm_time_UTC = eval(['gcm_time_UTC_CALIPSO_monthly']);
                 time_series_type='CALIPSO'; %default type of time data
                 
                 
              case {'CALIPSO','CALIPSO_monthly'}
                 Plat2D_edges = eval(['gcm_Plat2D_edges_ERAInt_monthly;']);
                 Plon2D_edges = eval(['gcm_Plon2D_edges_ERAInt_monthly']);
                 Plat2D = eval(['gcm_Plat2D_ERAInt_monthly']);
                 Plon2D = eval(['gcm_Plon2D_ERAInt_monthly']);      
                 gcm_time_of_day_select=0;
                 gcm_time_UTC = eval(['gcm_time_UTC_ERAInt_monthly']);
                 time_series_type='CALIPSO'; %default type of time data    
                 
                 
             case 'MODIS'
                 Plat2D_edges = repmat([LAT+0.5 LAT(end)-0.5],[size(LON,2)+1 1])';
                 Plon2D_edges = repmat([LON-0.5 LON(end)+0.5],[size(LAT,2)+1 1]);
                 Plat2D       = repmat([LAT],[size(LON,2) 1])';
                 Plon2D       = repmat([LON],[size(LAT,2) 1]);       
                 gcm_time_of_day_select=0;
                 
                 if size(Plat2D,1) ~= size(dat_modis,1)
                     fprintf(1,'\n*** PROBLEM - lat/lon size is different to dat_modis - check that the correct dat_modis is present\n');
                     return
                 end
                 if size(Plat2D,2) ~= size(dat_modis,2)
                     fprintf(1,'\n*** PROBLEM - lat/lon size is different to dat_modis - check that the correct dat_modis is present\n');
                     return
                 end
                 time_series_type='daily'; %default type of time data
                 
                 if ~exist('gcm_years_loaded_str')
                     gcm_years_loaded_str = modisyears_str;
                 end
                 
             case 'CLOUDSAT_PRECIP'
                   Plat2D_edges = Plat2D_matt_edges;
                   Plon2D_edges = Plon2D_matt_edges;
%                 Plat2D_edges = repmat([LAT+0.5 LAT(end)-0.5],[size(LON,2)+1 1])';
%                 Plon2D_edges = repmat([LON-0.5 LON(end)+0.5],[size(LAT,2)+1 1]);
                 Plat2D       = Plat2D_matt_centres;
                 Plon2D       = Plon2D_matt_centres;
                 gcm_time_of_day_select=0;
                 time_series_type = 'monthly calipso precip';
                 
                 if size(Plat2D,1) ~= size(dat_modis,2)
                     fprintf(1,'\n*** PROBLEM - lat/lon size is different to dat_modis - check that the correct dat_modis is present\n');
                     return
                 end
                 if size(Plat2D,2) ~= size(dat_modis,3)
                     fprintf(1,'\n*** PROBLEM - lat/lon size is different to dat_modis - check that the correct dat_modis is present\n');
                     return
                 end
                 
                 
             otherwise
                 Plat2D_edges = eval(['gcm_Plat2D_edges_' gcm_strs{idat} ';']);
                 Plon2D_edges = eval(['gcm_Plon2D_edges_' gcm_strs{idat} ';']);
                 Plat2D = eval(['gcm_Plat2D_' gcm_strs{idat} ';']);
                 Plon2D = eval(['gcm_Plon2D_' gcm_strs{idat} ';']);      
                 gcm_time_of_day_select=2;
                 gcm_time_UTC = eval(['gcm_time_UTC_' gcm_strs{idat} ';']);
                 time_series_type='daily'; %default type of time data
         end                                  
                 
         
         daynum_timeseries3 = eval(['daynum_timeseries3_' gcm_strs{idat} ';']);

%% Calculate the longitude range required and store in ilon_vals
         ihtot_max = find(Plat2D_edges>=LAT_val(1) & Plat2D_edges<LAT_val(2) & Plon2D_edges<LON_val(2) );
         [ilats_max,ilons_max]=ind2sub(size(Plat2D_edges),ihtot_max);
         ihtot_min = find(Plat2D_edges>=LAT_val(1) & Plat2D_edges<LAT_val(2) & Plon2D_edges>=LON_val(1) );
         [ilats_min,ilons_min]=ind2sub(size(Plat2D_edges),ihtot_min);
         
         ilon_end = maxALL(ilons_max);
         ilon_start = minALL(ilons_min);
%                thresh_str3=['LAT.GTE.' num2str(LAT_val(idat,1)) ' AND LAT.LTE.' num2str(LAT_val(idat,2)) ' AND LON.GT.' num2str(LON_val(idat,1)) ' AND LON.LTE.' num2str(LON_val(idat,2))];

          ilon_vals=ilon_start:ilon_end;
          
          

time_inds_modisL3_timeseries3
if gcm_time_of_day_select==2
%    time_inds_average=time_inds_average2; %for the approach where a 4D field of zeros is added with NaNs 
    %in places with no data
end

Y_data_all=[];

clear cfad_csection_lons cfad_csection


%% Loop through all the longitude values in ilon_vals
        for ilon_pt=1:length(ilon_vals)
            
            if icfad==1
                if size(time_inds_average2,1)==1
                    %if are just getting a 1D time_inds_average then make
                    %up a 4D array of zeros with NaNs in the places for
                    %the time indices we don't want
                    scf = eval(['size(cf_CFAD_sr_' gcm_strs{idat} ')']);
                    itime_nan = daynum_timeseries3_CALIPSO;
                    itime_nan(time_inds_average)=[];
                    time_inds_av4D = zeros([scf(1) scf(2) scf(3) length(ilon_vals(ilon_pt))]);
                    time_inds_av4D(itime_nan,:,:,:)=NaN;
                else
                    time_inds_av4D = eval(['repmat(time_inds_average2(:,:,ilon_vals(ilon_pt)),[1 1 size(cf_CFAD_dbz_' gcm_strs{idat} ',2)])']);
                    time_inds_av4D = permute(time_inds_av4D,[1 3 2]);
                end
            end

                lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

         % -------  for the satellites use the 2D data in dat_modis (as processed by plot_global_maps) 
         % Are selecting all the lat and time values for a particular longitude index
         % Putting the data into Y_data(time,lats) for one lon index
                switch gcm_strs{idat}
                    %For MODIS set the field here as it is the same each time
                    case {'MODIS','CLOUDSAT PRECIP','AMSRE','POLDER'}
                        %MODIS is ordered (lat,lon,time) -
                        %run plot_global_maps first to do
                        %screening
                        var_str = ['dat_modis'];
                        eval_str = ['Y_data = permute(squeeze(' var_str '(:,ilon_vals(ilon_pt),:)),[2 1]); Nd_str = ''Max Nd in cloud layer'';'];                   
                        
                    case {'CLOUDSAT_PRECIP','CALIPSO','CALIPSO_monthly','ERAInt'}
                        % is ordered (time,lat,lon) -
                        %run plot_global_maps first to do
                        %screening
                        var_str = ['dat_modis'];
                        eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt))); Nd_str = ''Max Nd in cloud layer'';'];                                           
                        time_inds_average = [1:size(dat_modis,1)];
                end

                                 
                
         % -------  Now for the models - aren't using pre-processed dat_modis data here 
         % Are selecting all the lat and time values for a particular longitude index
         % For 2D slices would also select all of the height values
         % Putting the data into Y_data(time,lats) for one lon index
                switch gcm_case
                    
                     case {'qv700'}
                        switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE','ERAInt'}
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_qv700_' gcm_strs{idat}];
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
    
                                 end                                 
                                 eval(eval_str);
                        end
                        

                                ylab='qv700 (g kg^{-1})';
  
                        
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=16;
                        
                        
                        
                    
                    case {'Re37_minus_Re21','Re16_minus_Re21'}
                        switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
%                                         var_str = ['1e-2*gcm_ps_' gcm_strs{idat}]; 
%                                         eval_str = ['Y_data = squeeze( ' var_str2 '(:,:,ilon_vals(ilon_pt))  + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
    
                                 end                                 
                                 eval(eval_str);
                        end
                        
                        switch gcm_case
                            case 'Re37_minus_Re21'
                                ylab='Reff_{3.7} minus Reff_{2.1} (\mum)';
                            case 'Re16_minus_Re21'
                                ylab='Reff_{1.6} minus Reff_{2.1} (\mum)';                                                            
                        end
                        
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=-5;
                        zmax=2;
                        
                    
                    case {'Re16','Re21','Re37'}
                        switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
%                                         var_str = ['1e-2*gcm_ps_' gcm_strs{idat}]; 
%                                         eval_str = ['Y_data = squeeze( ' var_str2 '(:,:,ilon_vals(ilon_pt))  + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
    
                                 end                                 
                                 eval(eval_str);
                        end
                        
                        switch gcm_case
                            case 'Re16'
                                ylab='Effective radius for 1.6\mum (\mum)';
                            case 'Re21'
                                ylab='Effective radius for 2.1\mum (\mum)';                                
                            case 'Re37'
                                ylab='Effective radius for 3.7\mum (\mum)';                                
                        end
                        
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=6;
                        zmax=25;  
                        
                    case 'Sea Level Pressure_vert_slice'
                        switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         %gcm_liq_av is in kg/kg
                                         var_str = ['1e-2*gcm_ps_' gcm_strs{idat}]; 
                                    %Doing a height-lon transect -
                                    %repeating this watervap function for
                                    %each height (zslice_index)
                                          eval_str = ['Y_data = squeeze( ' var_str2 '(:,:,ilon_vals(ilon_pt))  + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
    
                                 end                                 
                                 eval(eval_str);
                        end
                        
                        
                        ylab='Sea Level Pressure (hPa)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=1000;  
                        
case 'Pressure_vert_slice'
                        switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         %gcm_liq_av is in kg/kg
                                         var_str = ['1e-2*gcm_phalf_' gcm_strs{idat}]; 
                                         var_str2 = ['1e-2*gcm_ps_' gcm_strs{idat}];                                          
                                    %Doing a height-lon transect -
                                    %repeating this watervap function for
                                    %each height (zslice_index)
                                          eval_str = ['Y_data = squeeze(squeeze(' var_str '(:,zslice_index,:,ilon_vals(ilon_pt))) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
    
                                 end                                 
                                 eval(eval_str);
                        end
                        
                        
                        ylab='Pressure (hPa)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=1000;                          
                                 
                      case 'LWC_vert_slice (g kg^{-1})'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         %gcm_liq_av is in kg/kg
                                         var_str = ['1e3*gcm_liq_av_' gcm_strs{idat}]; 
                                    %Doing a height-lon transect -
                                    %repeating this watervap function for
                                    %each height (zslice_index)
                                          eval_str = ['Y_data = squeeze(squeeze(' var_str '(:,zslice_index,:,ilon_vals(ilon_pt))) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
    
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
%                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Grid-box mean LWC (g kg^{-1})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=3;                         
                        
                        
                    case 'CTP'
                        Y_data = squeeze(gcm_CTP(:,:,ilon_vals(ilon_pt))/100);
                        ylab='Cloud Top Pressure (hPa)';                          
                    case 'H'
                        Y_data = squeeze( gcm_CTH(:,:,ilon_vals(ilon_pt)) - gcm_CBH(:,:,ilon_vals(ilon_pt)) );
                        ylab='Cloud Depth (m)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=150;
                        zmax=600;
                        
                    case {'LTS','LTS_daily'}
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'ERAInt'}
                                        %set above
                                     otherwise
                                         var_str = ['gcm_LTS_' gcm_strs{idat}];
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Lower Tropospheric Stability (K)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=30; 
                        
                        
                     case 'LTS1000'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'ERAInt'}
                                        %set above
                                     otherwise
                                         var_str = ['gcm_LTS1000_' gcm_strs{idat}];
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Lower Tropospheric Stability based on 1000 hPa (K)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=30;    
                        
                        
                    case 'LWP_COSPCF'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case 'MODIS'
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_lwp_COSPCF_' gcm_strs{idat}];
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='LWP COSP CF filtering (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;    

                    case 'LWP_COSP'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case 'MODIS'
                                        %set above
                                     otherwise
                                         var_str = ['1e3*LWP_COSP_' gcm_strs{idat}];
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='LWP COSP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;    
                  
                        
                    case 'LWP'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case 'MODIS'
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_lwp_minthreshCF_' gcm_strs{idat}];
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Grid-box mean LWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;
                        
                    case 'LWP2'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_lwp_' gcm_strs{idat}];                                        
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Grid-box mean LWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;    
                        
                    case 'TLWP' %grid-box average LWP+RWP for comparison to AMSRE
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_TLWP_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
%                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Grid-box mean TLWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;  

                    case {'LWP_MOD35CF','LWP_MOD06CF'} %MODIS grid-box average LWP using MOD35 CF
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_TLWP_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
%                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Grid-box mean TLWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;  

                        
                        
 case 'TWP' %grid-box average LWP+RWP for comparison to AMSRE
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_TWP_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
%                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Grid-box mean TWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;  
                        
                        
                        
                        
                    case {'REFFL_maxlayer','REFFL_maxliq','REFFL_max_noCF'} %MODIS liq Reff as calculated in gcm_process
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['1e6*gcm_' gcm_case '_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab=[gcm_case ' Effective Radius (\mum)'];
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=30;  
                        
                        
                        
                    case 'Reff_COSP' %COSP MODIS liq Reff
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE','POLDER'}
                                        %set above
                                     otherwise
                                         var_str = ['1e6*liqReCF80_modis_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Effective Radius (\mum)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=30;  
                        
                     case 'Reff_COSP_allCF' %COSP MODIS liq Reff
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE','POLDER'}
                                        %set above
                                     otherwise
                                         var_str = ['1e6*liqRe_allCF_modis_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Effective Radius (\mum)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=30;      
                        

                        
                        
                        
        case 'Tau_COSP' %COSP MODIS liq Reff
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['liqTauCF80_modis_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Optical Depth';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=30;    
                        
        case 'Nd_COSP' %COSP MODIS liq Reff
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case {'MODIS','AMSRE'}
                                        %set above
                                     otherwise
                                         var_str = ['Nd_COSP_CF80_' gcm_strs{idat}];
                                         if slice_2D==1 %if are doing a height-lon transect
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         else %normal longitude line plot
                                             eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                         end
                                 end                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='Droplet Concentration (cm^{-3})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=250;                         
                        
                        
                        
                        
                        
                    case 'CloudSat dbZ from CFAD'
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case 'MODIS'
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_lwp_minthreshCF_' gcm_strs{idat}];
                                         

%here are multiplying the mean dBZ for each height by the CF (=N_CFAD_dbz)
%to give a grid box average dbz. Are also working on 10.^dbz
eval_str = ['Y_data = squeeze( time_inds_av4D + 10.^(0.1*mean_CFAD_dbz_' gcm_strs{idat} '(:,:,:,ilon_vals(ilon_pt))).* cf_CFAD_dbz_' gcm_strs{idat} '(:,:,:,ilon_vals(ilon_pt)));'];
%mean_CFAD_dbz is [time height lat lon]
                                         %eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)));'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;
                        
                        
                    case 'CALIPSO Cloud Fraction from CFAD'
                        icfad=1;
                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case 'MODIS'
                                        %set above
                                     otherwise
                                         var_str = ['1e3*gcm_lwp_minthreshCF_' gcm_strs{idat}];
                                         

%here are multiplying the mean dBZ for each height by the CF (=N_CFAD_dbz)
%to give a grid box average dbz. Are also working on 10.^dbz
eval_str = ['Y_data = squeeze(time_inds_av4D + cf_CFAD_sr_' gcm_strs{idat} '(:,:,:,ilon_vals(ilon_pt)) );'];
%mean_CFAD_dbz is [time height lat lon]
                                         %eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)));'];
                                 end
                                 
                                 eval(eval_str);
                                 
                             case 0
                                 %                        Y_data = gcm_lwp_cloud(ihtot);
                                 Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                         end


                        ylab='';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;
                        
                        
                        
                        
                    case 'Nd'
%                        Y_data = gcm_Nd_meanliq(ihtot);
%                        Y_data = gcm_Nd_maxALL(ihtot);
%                        Y_data = gcm_Nd_max(ihtot);

                         switch use_saved_dat
                             case 1
                                 switch gcm_strs{idat}
                                     case 'MODIS'
                                         %do nothing as set above
                                     otherwise
                                         var_str = ['gcm_Nd_maxlayer_' gcm_strs{idat}];
%                                         var_str = ['gcm_Nd_max_noCF_' gcm_strs{idat}];  
%                                         var_str = ['gcm_Nd_maxliq_' gcm_strs{idat}];                                           
                                         
                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) ); Nd_str = ''Max Nd in cloud layer'';'];

                                 end   
                                 
                              eval(eval_str);
                              
                          case 0
                              Y_data = squeeze(gcm_Nd_maxlayer(:,:,ilon_vals(ilon_pt))); Nd_str = 'Max Nd in cloud layer';
                              %                        Y_data = squeeze(gcm_Nd_max_lwc_cont(:,:,ilon_vals(ilon_pt))); Nd_str = 'Nd at max LWC in top layer';
                      end
                        ylab=['Droplet Concentration (cm^{-3}), ' Nd_str];
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot panel
                        
                         izlim=1;
                         zmin=0;
                         zmax=300;
                         
                    case 'Terrain Height'
                        switch gcm_str
                            case 'AM3'
                                Y_data = (repmat(squeeze(gcm_zsurf(:,ilon_vals(ilon_pt))),[size(gcm_CTP,1) 1]));
                            case {'CAM5','CAM5_CLUBB'}
                                Y_data = squeeze(gcm_surfgeo(:,:,ilon_vals(ilon_pt)))/9.80616;
                        end
                        ylab='Surface Height (m)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                        
%                        thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%            ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1}'];
%                        thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'];        
                        
                          izlim=0;
                          zmin=0;
                          zmax=1;
        
                         thresh_str=[pstr];
                         
                    case 'LWP'
%                        Y_data = gcm_lwp_cloud(ihtot);
                        Y_data = squeeze(1e3*gcm_lwp_minthreshCF(:,:,ilon_vals(ilon_pt)));
                        ylab='Grid-box mean LWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;   
                        
                     case 'Land Fraction'
                        Y_data = (repmat(squeeze(gcm_landmask(:,ilon_vals(ilon_pt))),[1 size(gcm_CTP,1)]))';
                        ylab='Land Fraction';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                        
%                        thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%            ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1}'];
%                        thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'];        
                        
                          izlim=0;
                          zmin=0;
                          zmax=1;
        
                         thresh_str=[pstr];
                         
                         
                    case 'Using saved GCM data'
                        

                        
                        
                        eval_str_screen='';
                          switch use_saved_dat
                             case 1
                                
                                         
                                         switch var_choose
                                             case 'Precip rate'
                                                 inormalise_precip=0;
                                                 if inormalise_precip==1
                                                     normalise_data_get %script

                                                 else
                                                     ylab='Precip rate (mm hr^{-1})';
                                                     normalise_data=1;
                                                 end

                                                 switch gcm_strs{idat}
                                                     case {'MODIS','CLOUDSAT_PRECIP','CALIPSO','CALIPSO_monthly'}
                                                         %set above
                                                     otherwise                                                         
                                                         gcm_precip_norm = eval(['gcm_precT_' gcm_strs{idat} './normalise_data;']);
                                                         %                                                 var_str = ['gcm_precT_' gcm_strs{idat}];
                                                         var_str = 'gcm_precip_norm';
                                                         %convert from m/s to mm/hr
%                                                         eval_str = ['Y_data = 3600*1e3*squeeze(' var_str '(:,:,ilon_vals(ilon_pt)));'];
                                                         eval_str = ['Y_data = 3600*1e3*squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];                                                         
                                                 end
                                                  
                                                 lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                                                 izlim=1;
                                                 zmin=0;
                                                 zmax=0.1;
                                                 
                                             case 'low CF'
                                                 
                                                 switch gcm_strs{idat}
                                                     case {'MODIS','CLOUDSAT_PRECIP','CALIPSO','CALIPSO_monthly'}
                                                         %set above
                                                     otherwise
%                                                         var_str = ['gcm_CF_max_screened_' gcm_strs{idat}];
                                                          var_str = ['cf_isccp_low_' gcm_strs{idat}];
                                                         
                                                         %                                         var_str = ['gcm_CF_maxliq_' gcm_strs{idat}];
                                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                                 end

%                                                 eval(eval_str);

                                                 ylab='Low Cloud Fraction';
                                                 lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                                                 izlim=1;
                                                 zmin=0;
                                                 zmax=1;                                                 
                                                 
                                                 
                                             case 'high CF'
                                                 
                                                 switch gcm_strs{idat}
                                                     case {'MODIS','CLOUDSAT_PRECIP','CALIPSO','CALIPSO_monthly'}
                                                         %set above
                                                     otherwise
%                                                         var_str = ['gcm_CF_max_screened_' gcm_strs{idat}];
                                                          var_str = ['cf_isccp_high_' gcm_strs{idat}];
                                                         
                                                         %                                         var_str = ['gcm_CF_maxliq_' gcm_strs{idat}];
                                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)) + time_inds_average2(:,:,ilon_vals(ilon_pt)) );'];
                                                 end

%                                                 eval(eval_str);

                                                 ylab='High Cloud Fraction';
                                                 lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                                                 izlim=1;
                                                 zmin=0;
                                                 zmax=1;
                                                 
                                                                                           case 'low CF'
                                                 
                                                 switch gcm_strs{idat}
                                                     case {'MODIS','CLOUDSAT_PRECIP','CALIPSO','CALIPSO_monthly'}
                                                         %set above
                                                     otherwise
%                                                         var_str = ['gcm_CF_max_screened_' gcm_strs{idat}];
                                                          var_str = ['cf_isccp_low_' gcm_strs{idat}];
                                                         
                                                         %                                         var_str = ['gcm_CF_maxliq_' gcm_strs{idat}];
                                                         eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)));'];
                                                 end

%                                                 eval(eval_str);

                                                 ylab='Low Cloud Fraction';
                                                 lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                                                 izlim=1;
                                                 zmin=0;
                                                 zmax=1;   
                                                 

                                             case 'Terrain Height'
                                                  var_str = ['gcm_zsurf_' gcm_strs{idat}];
                                                  svar = eval(['size(' var_str ');']);
                                                  if length(svar)==2
                                                      eval_str = ['Y_data = (squeeze(' var_str '(:,ilon_vals(ilon_pt))))'';'];
                                                  else
                                                      eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)));'];                                                      
                                                  end
                                                  ylab='Terrain Height (m)';
                                                  lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                                                  izlim=0;
                                                  zmin=0;
                                                  zmax=0.1;
                                                  
                                                  time_inds_average=1;
                                                  
                                              case 'Land Fraction'
                                                  var_str = ['gcm_landmask_' gcm_strs{idat}];
                                                  svar = eval(['size(' var_str ');']);
                                                  if length(svar)==2
                                                      eval_str = ['Y_data = (squeeze(' var_str '(:,ilon_vals(ilon_pt))))'';'];
                                                  else
                                                      eval_str = ['Y_data = squeeze(' var_str '(:,:,ilon_vals(ilon_pt)));'];                                                      
                                                  end
                                                  ylab='Land Fraction';
                                                  lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                                                  izlim=0;
                                                  zmin=0;
                                                  zmax=0.1;
                                                  
                                                  time_inds_average=1;
                                                  
                                                  
                                         end
                                         
                                        
                                                                                                                         

                                 
                                 eval(eval_str); %also sets up the MODIS/CLOUDSAT data
                                 if iscreen_CTT==1
                                     iscreen = find(~ (screen_data>=gcm_CTT_thresh(1) & screen_data<=gcm_CTT_thresh(2) ) );
                                      eval_str_screen=['Y_data(iscreen) = NaN;'];
                                     var_str2 = ['gcm_CTT_layer_' gcm_strs{idat}];
                                     eval_str2 = ['screen_data = squeeze(' var_str2 '(:,:,ilon_vals(ilon_pt)));'];
                                     eval(eval_str2);
                                     eval(eval_str_screen);
                                 end
                                 
                             case 0                                
                                 Y_data = squeeze( gcm_precT(:,:,ilon_vals(ilon_pt)) );
%                                 Y_data = squeeze( gcm_CF_maxliq(:,:,ilon_vals(ilon_pt)) );                                 
                          end
                         
        
                end  %switch gcm_case
                
                
%% Now finds all of the points that are within the required latitude range (one longitude)
%   Puts into ihtot2 (linear indices)
                ihtot2 = find ( Plat2D_edges(:,ilon_vals(ilon_pt))>=LAT_val(1) & Plat2D_edges(:,ilon_vals(ilon_pt))<LAT_val(2) );
%                thresh_str3=['LAT.GTE.' num2str(LAT_val(ilon_pt,1)) ' AND LAT.LTE.' num2str(LAT_val(ilon_pt,2)) ' AND LON.GT.' num2str(LON_val(ilon_pt,1)) ' AND LON.LTE.' num2str(LON_val(ilon_pt,2))];
                ihtot2 = ihtot2(1:end-1); %added this 24th April, 2013
                     %  -- the last ihtot2 value will be the edge for which
                     %  the datapoint to the east will not lie exclusively
                     %  within the range selected (its edge will be outside
                     %  the range). So we ignore this datapoint.
                     
                lat_vals_all = Plat2D(ihtot2,ilon_vals(ilon_pt));
                mean_lats(ilon_pt) = mean(lat_vals_all);
                
                 if ione_lat==1
                     lat_mid = 0.5* (LAT_val(1) + LAT_val(2));                     
                     [minval imin]=min(abs(lat_vals_all - lat_mid));
                     
%                     ihtot2 = round(median(ihtot2));
                      ihtot2 = ihtot2(imin);                                                               
                     %    ihtot2=14;
                     %ihtot2=11;
%                     if idat==1
%                         ihtot2=112;
%                     end
                      single_lat_vals(idat).val(ilon_pt)  = Plat2D(ihtot2,ilon_vals(ilon_pt));
                 end

%% Now take the ihtot2 linear indices (one lon value) for the 2D x-y plane
%% and replicate over time with the aim of producing linear indices (ihtot)
%% that reference the required XY points at the required times. Then Y_data
%% is reduced to just select this data
% For a 2D slice will want to repeat this process for every height - can
% do this from outside of this - i.e. just run this watervap function multiple times for each
% height required
          if icfad==0
                 
                %these are the lat/lon indices from all latitudes, but only one lon:-
%ihtot2 = find ( Plat2D_edges(:,ilon_vals(ilon_pt))>=LAT_val(1) & Plat2D_edges(:,ilon_vals(ilon_pt))<LAT_val(2) );- need to repeat these for
                %each time and calculate a new index for the 3D matrix
                %lat/lon indices
                IJ = repmat(ihtot2,[1 length(time_inds_average)]);
                %K - replicate all of the time indices
                K = repmat([time_inds_average],[length(ihtot2) 1]);
                ihtot = sub2ind([size(Y_data,1) size(Y_data,2)] , K(:) ,IJ(:));
                
                Y_data = Y_data(ihtot);
                
else
    code_switch = '4D';
    switch code_switch
        case 'old'
                    %these are the lat/lon indices from all latitudes, but only one lon:-
%ihtot2 = find ( Plat2D_edges(:,ilon_vals(ilon_pt))>=LAT_val(1) & Plat2D_edges(:,ilon_vals(ilon_pt))<LAT_val(2) );- need to repeat these for
                %each time and calculate a new index for the 3D matrix
                %lat/lon indices. size(ihtot2)=[N 1]
                sYdata = size(Y_data);
                %for obs = [48 40 90] i.e. [time height lat]
                nTim   = max(size(time_inds_average));
                %when using local time time_inds_aveage is now a 3D array,
                %so use size and not length
                %actually ahve two different cases between the obs and GCMs
                %with obs have a 1D array, but with GCM have 4D
                
                %replicate the lat indices 
                IJ = repmat(ihtot2,[1 nTim size(Y_data,2)]);
                IJ = permute(IJ,[2 3 1]);
                %K - replicate all of the time indices
                %size(time_inds_average) = [1 N]
                K = repmat(time_inds_average,[length(ihtot2) 1 size(Y_data,2)]);
                K = permute(K,[2 3 1]);
      
                %repmat the heights (use all heights)
                H = repmat([1:size(Y_data,2)],[length(ihtot2) 1 nTim]);
                H = permute(H,[3 2 1]);
                ihtot = sub2ind([size(Y_data,1) size(Y_data,2) size(Y_data,3)] , K(:) , H(:), IJ(:));
                
                Y_data = Y_data(ihtot);
                SIJ=size(IJ);
                Y_data = reshape(Y_data,[SIJ(1) SIJ(2) SIJ(3)]);
                
        case '4D'
%these are the lat/lon indices from all latitudes, but only one lon:-
%ihtot2 = find ( Plat2D_edges(:,ilon_vals(ilon_pt))>=LAT_val(1) & Plat2D_edges(:,ilon_vals(ilon_pt))<LAT_val(2) );
%need to replicate these using repmat for EACH ihtot2 value - these are the
%lat values that we want to select from Y_data
%(replicate over each time and height) and calculate a new index for the 3D matrix
%lat/lon indices. size(ihtot2)=[N 1]
                
                %now we have the 4D time_inds_av4D, which contains NaNs at
                %all locations that we want removing.
                sYdata = size(Y_data);
                %for obs = [48 40 90] i.e. [time height lat]
                nTim   = size(time_inds_av4D,1);
                %when using local time time_inds_aveage is now a 3D array,
                %so use size and not length
                %actually ahve two different cases between the obs and GCMs
                %with obs have a 1D array, but with GCM have 4D
                

                IJ = repmat(ihtot2,[1 nTim sYdata(2)]);
                IJ = permute(IJ,[2 3 1]);
                %K - replicate all of the time indices
                %are using all time indices here since we have the 4D array
                %with NaNs in 
                K = repmat([1:nTim],[length(ihtot2) 1 sYdata(2)]);
                K = permute(K,[2 3 1]);
      
                %repmat the heights (use all heights)
                H = repmat([1:sYdata(2)],[length(ihtot2) 1 nTim]);
                H = permute(H,[3 2 1]);
                ihtot = sub2ind([size(Y_data,1) size(Y_data,2) size(Y_data,3)] , K(:) , H(:), IJ(:));
                
                Y_data = Y_data(ihtot);
                SIJ=size(IJ);
                Y_data = reshape(Y_data,[SIJ(1) SIJ(2) SIJ(3)]);
                
    end
    
    
end
                
                
          
%                xdat(1).x(ilon_pt) = mean(Plon2D_edges(ihtot2,ilon_vals(ilon_pt)));  

xdat(idat).x(ilon_pt) = mean(Plon2D(ihtot2,ilon_vals(ilon_pt))); 
lon_edges(ilon_pt) = mean(Plon2D_edges(ihtot2,ilon_vals(ilon_pt)-1)); 
lon_edges(ilon_pt+1) = mean(Plon2D_edges(ihtot2,ilon_vals(ilon_pt)));

if icfad==1
    %if using CFADS then we have the grid box (cloud fraction weighted)
    %values, so we can just do the mean
    %average over latitude
    [xx,Nmean] = meanNoNan(Y_data,3);
    %    xx = xx.*Nmean;  %%this now gives the sum rather than the mean

    %now avearage over time
    [xx2,Nmean] = meanNoNan(xx,1);
    %    cfad_csection(ilon_pt,:) = xx2.*Nmean;

    switch gcm_case
        case 'CloudSat dbZ from CFAD'
            xx2 = 10*log10(xx2);
    end

    cfad_csection(ilon_pt,:) = xx2;

    ydat(idat).y(ilon_pt) = 0;
    Y_data_all=0;

    cfad_csection_lons(ilon_pt) = xdat(idat).x(ilon_pt);

else

    %special case where onnly have one time index in the data
    if size(Y_data,1)==1
        Y_data = Y_data'
    end

    [ydat(idat).y(ilon_pt),Ndat_gcm(idat).N(ilon_pt),std_dev_gcm(idat).s(ilon_pt)] = meanNoNan(Y_data,1); %
    [sumsq(idat).y(ilon_pt)] = meanNoNan(Y_data.^2,1,'sum'); %sum of the squares - can be used for std dev
    %collect all the data together for a final mean at the end
    Y_data_all = [Y_data_all; Y_data];

end
                
               

                end %switch gcm_case
        

        mean_LAT_val(idat) = mean(mean_lats);
        
           %            labs(idat).l=[thresh_str2 ',mean=' num2str(mean_val(idat),'%.1f')];
            labs(idat).l=[thresh_str2 ',me=' num2str(meanNoNan(Y_data_all,1),'%.2f')];

            
        end %ilon_pt
        
        if length(time_UTC_str)>0
            time_UTC_str = [' for ' time_UTC_str ' UTC at LAT=' num2str(mean_LAT_val)];
        else
            time_UTC_str = [' at LAT=' num2str(mean_LAT_val)];
        end

        
%            titlenam = [titlenam ' for ' time_UTC_str ' UTC ' thresh_str '
%            LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(2))];



             if use_saved_dat==1            
                 titlenam = [titlenam 'comparison ' time_UTC_str gcm_years_loaded_str thresh_str];
             else
                 titlenam = [titlenam gcm_str time_UTC_str thresh_str];
             end
             
              if iscreen_CTT==1
                  titlenam = [titlenam ' CTT.GTE.' num2str(gcm_CTT_thresh(1))];
              end
             
            if ione_lat==1
                titlenam=[titlenam ' SINGLE LAT'];
            end
            
            titlenam_orig = titlenam;
            titlenam=[titlenam ' ' time_mean_str '_days '];
            
            

            clear text_dat xtext_dat ytext_dat
            for idat=2:length(Ndat_gcm)
                for itext=1:length(Ndat_gcm(idat).N)
                    str_text=num2str(Ndat_gcm(idat).N(itext));
                    text_dat(idat).text(itext,1:length(str_text)) = str_text;
                end
                if izlim==1
                    dz = ( zmax - zmin )/30;
                else
                    dz = ( max(ydat(idat).y) - min(ydat(idat).y) ) /30;
                end
                xtext_dat(idat).x = xdat(idat).x;
                ytext_dat(idat).y = ydat(idat).y + dz;

            end

switch gcm_case
    case 'Using saved GCM data'
        gcm_case = remove_character(var_choose,' ','_');
        
        switch gcm_str
            case 'CLOUDSAT_PRECIP'
                gcm_case = [gcm_case '_' asc_desc];
            case 'CALIPSO_monthly'
%                gcm_case = [gcm_case '_' calipso_daynight_label];
        end
        
end

        
        figname=titlenam;
        savename=[savedir figname];    
        
        
    case 114 %timeseries from values from monthly_means_from_plot_global script             
        titlenam = ['mean values timeseries ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=1000*[0 0.025];
        

   

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab = season_vals_save(1).NdPDF_xlab;
        xlab=x_lab_seasonal(1).lab;
        
        %This bit just changes the y limits based on the ylab already set
        %in PlotTimeHeightVap3
                switch ylab
                    case {'LWP (g m^{-2}) Painemal','LWP (g m^{-2})'}
                        izlim=1;
                        zmin=0;
                        zmax=150;
                    case {'Cloud Depth (m)','Cloud Depth (m) Painemal'}
                        izlim=1;
                        zmin=150;
                        zmax=400;    
                    case 'N_d (cm^{-3})'
                        izlim=1;
                        zmin=0;
                        zmax=160; 
                    case 'Cloud Fraction'
                        izlim=1;
                        zmin=0;
                        zmax=1;     
                end

%        y_axis_type = 'log10_matlab';


        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        
        for idat=1:length(time_means_multi)


%        xdat(idat).x=X_mean; %
       
        ydat(idat).y = time_means_multi(idat).dat; %
        switch multi_case
            case 'Bimonthly'
                xdat(idat).x = [1.5:2:11.5];
            otherwise               
                xdat(idat).x = 1:length(ydat(idat).y);
        end
        
        labs(idat).l=['LON=' LON_str_multi(idat).dat, ',LAT=' LAT_str_multi(idat).dat ' me=' num2str(Xmean_multi(idat).dat,'%3.2f')];
        
        end
        
       
    
    case 113 %gcm dirunal cycle timeseries
        


        
        if icf_low_only==1
            pstr = ['.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'];
        else
            pstr='';
        end



        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
%        y_axis_type = 'log10_matlab';
        %LON_val = [-76.5 -72]; LON_val=[-81 -76.5]; LON_val = [-85.5-81];
        %LON_val = [-90 -85.5];

        %old lon bins
        LAT_val = [-22 -18; -22 -18; -22 -18; -22 -18]; LON_val = [-76.5 -72; -81 -76.5; -85.5 -81; -90 -85.5];   
        %VOCA ones - plus Klein and Hartmann
        LAT_val = [-22 -18; -22 -18; -22 -18; -20 -10]; LON_val = [-75 -70; -80 -75; -86 -80; -90 -80];
        LAT_val = [-22.74 -18; -22.74 -18; -22.74 -18; -20 -10]; LON_val = [-76.25 -71.25; -81.25 -76.25; -87.25 -81.25; -91.25 -81.25];
%        LAT_val = [-22.74 -18]; LON_val = [-76.25 -71.25];
         
        %think there is a 3-4 hour time difference for LT in Chile e.g. LT = UTC-4
        times_of_day = [0 3 6 9 12 15 18 21]; %UTC
        times_of_day = [15 18 21];  
%        times_of_day = [6 9 12]; %UTC   

time_UTC_str='';
        
        time_scale = 'diurnal';
        time_scale = 'monthly';   
        time_scale = 'bimonthly';
%        time_scale = 'seasonal';           

        irestrict_times=1;


        gcm_case = 'CTP';
        gcm_case = 'H';
%        gcm_case = 'LWP';
%        gcm_case = 'Nd';
%        gcm_case = 'CF';  
        
        


         months=[1:12];   
         
 
        

         switch time_scale
             case 'diurnal'
                 times_plot = times_of_day;
                 xlab = 'Time (UTC)';
                 titlenam = ['gcm diurnal cycle of ' gcm_case ' for '];
             case 'monthly'
                 times_plot = months;
                 xlab = 'Month';
                 titlenam = ['gcm monthly of ' gcm_case ' for '];
             case 'seasonal'
                 times_plot = [12 1 2; 3 4 5; 6 7 8; 9 10 11];
                 times_plot = times_plot';
                 xlab = 'Season';
                 titlenam = ['gcm seasonal of ' gcm_case ' for '];  
                 
             case 'bimonthly'
                 times_plot = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12];
                 times_plot = times_plot';
                 xlab = 'Month';
                 titlenam = ['gcm bimonthly of ' gcm_case ' for '];      
                 
         end

      
     



        clear mean_val
        Ndat_gcm=[];
        thresh_str='';


        for idat=1:size(LON_val,1);
            Y_data_all=[];
            for itime_di = 1:size(times_plot,2)
                
                switch time_scale
                    case 'diurnal'
                        iday_gcm2D_inds = [find(gcm_time_UTC==times_plot(itime_di))];
                    case {'monthly','seasonal','bimonthly'}
                        switch time_scale
                            case {'seasonal','bimonthly'}
                                iday_gcm2D_inds=[];
                                for imonth=1:size(times_plot,1)
                                    iday_gcm2D_inds2 = find(gcm_month==times_plot(imonth,itime_di));
                                    iday_gcm2D_inds=[iday_gcm2D_inds; iday_gcm2D_inds2];
                                end
                            case 'monthly'
                                iday_gcm2D_inds = [find(gcm_month==times_plot(itime_di))];
                        end
                        
                        
                        if irestrict_times == 1
                            time_inds_average2=[];
                            time_UTC_str=' ';
                            for itime_mean = 1:length(times_of_day)
                                %find the indices for a particular hour
                                % = times_of_day(itime_mean) from the
                                % indices already filtered for month/season
                                itimes_find = find(gcm_time_UTC(iday_gcm2D_inds)==times_of_day(itime_mean));
                                %combine in one array
                                time_inds_average2 = [time_inds_average2; iday_gcm2D_inds(itimes_find)];
                                time_UTC_str=[time_UTC_str num2str(times_of_day(itime_mean))  ','];
                            end
                            
                            time_UTC_str(end)=[];
                            iday_gcm2D_inds = time_inds_average2;



                        end                                        
                end
                


                ihtot2 = find ( gcm_Plat2D_edges>=LAT_val(idat,1) & gcm_Plat2D_edges<LAT_val(idat,2) & gcm_Plon2D_edges>=LON_val(idat,1) & gcm_Plon2D_edges<LON_val(idat,2) );
%                thresh_str3=['LAT.GTE.' num2str(LAT_val(idat,1)) ' AND LAT.LTE.' num2str(LAT_val(idat,2)) ' AND LON.GT.' num2str(LON_val(idat,1)) ' AND LON.LTE.' num2str(LON_val(idat,2))];

                [ilons,ilats]=ind2sub(size(gcm_Plat2D_edges),ihtot2);
                ihtot2_addlon = sub2ind(size(gcm_Plat2D_edges),ilons+1,ilats);
                ihtot2_addlat = sub2ind(size(gcm_Plat2D_edges),ilons,ilats-1);
%                thresh_str2=['LON=' num2str(minALL(gcm_Plon2D_edges(ihtot2)),'%.2f') ' to ' num2str(maxALL(gcm_Plon2D_edges(ihtot2))+dlon(1),'%.2f') ',LAT=' num2str(minALL(gcm_Plat2D_edges(ihtot2)),'%.2f') ' to ' num2str(maxALL(gcm_Plat2D_edges(ihtot2))+dlat(1),'%.2f')];
                thresh_str2=['LON=' num2str(minALL(gcm_Plon2D_edges(ihtot2)),'%.2f') ' to ' num2str(maxALL(gcm_Plon2D_edges(ihtot2_addlon)),'%.2f') ',LAT=' num2str(minALL(gcm_Plat2D_edges(ihtot2)),'%.2f') ' to ' num2str(maxALL(gcm_Plat2D_edges(ihtot2_addlat)),'%.2f')];                
                
                %these are the lat/lon indices - need to repeat these for
                %each
                %time and calculate a new index for the 3D matrix
                %lat/lon indices
                IJ = repmat(ihtot2',[length(iday_gcm2D_inds) 1]);
                %K - replicate all of the time indices
                K = repmat([iday_gcm2D_inds],[1 length(ihtot2)]);
                ihtot = sub2ind([size(gcm_CTP,1) size(gcm_CTP,2)*size(gcm_CTP,3)] , K(:) ,IJ(:));

                lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane


                switch gcm_case
                    case 'CTP'
                        Y_data = gcm_CTP(ihtot)/100;
                        ylab='Cloud Top Pressure (hPa)';

                    case 'H'
                        Y_data = (gcm_CTH(ihtot)-gcm_CBH(ihtot));
                        ylab='Cloud Depth (m)';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=150;
                        zmax=400;
                    case 'LWP'
%                        Y_data = gcm_lwp_cloud(ihtot);
                        Y_data = 1e3*gcm_lwp_minthreshCF(ihtot);
                        ylab='Grid-box mean LWP (g m^{-2})';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=1;
                        zmin=0;
                        zmax=150;
                    case 'Nd'
%                        Y_data = gcm_Nd_meanliq(ihtot);
%                        Y_data = gcm_Nd_maxALL(ihtot);
%                        Y_data = gcm_Nd_max(ihtot);
%                        Y_data = gcm_Nd_max_screen(ihtot); Nd_str = 'Max after screening';
                        Y_data = gcm_Nd_max_lwc_cont(ihtot); Nd_str = 'Nd at max LWC in top layer';
                        ylab=['Droplet Concentration (cm^{-3}), ' Nd_str];
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot panel
                        
                         izlim=1;
                         zmin=0;
                         zmax=160;
                    case 'CF'
%                        Y_data = gcm_CF_max(ihtot);
                        Y_data = gcm_CF_max_screened(ihtot);
                        ylab='Max Cloud Fraction';
                        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                        
%                        thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%            ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1}'];
%                        thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'];        
                        
                          izlim=1;
                          zmin=0;
                          zmax=1;
        
                         thresh_str=[pstr];
                end
                
                
                switch time_scale
                    case 'seasonal'
                        xdat(idat).x = 1:size(times_plot,2);
                    case 'bimonthly'
                        xdat(idat).x = 0.5*(times_plot(1,:)+times_plot(2,:));
                    otherwise
                        xdat(idat).x(itime_di) = times_plot(itime_di);
                end
                
                [ydat(idat).y(itime_di),Ndat_gcm(idat).N(itime_di)] = meanNoNan(Y_data,1); %
                Y_data_all = [Y_data_all; Y_data];
                
            end
            
            



            %            labs(idat).l=[thresh_str2 ',mean=' num2str(mean_val(idat),'%.1f')];
            labs(idat).l=[thresh_str2 ',me=' num2str(meanNoNan(Y_data_all,1),'%.2f')];


        end
        
        if length(time_UTC_str)>0
            time_UTC_str = [' for ' time_UTC_str ' UTC '];
        end
        

        
%            titlenam = [titlenam ' for ' time_UTC_str ' UTC ' thresh_str ' LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(2))];
            titlenam = [titlenam gcm_str time_UTC_str thresh_str];            

        
        figname=titlenam;
        savename=[savedir figname];    
        
    case 112 %gcm 1D pdfs


        thresh_str=['CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];




        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        y_axis_type = 'log10_matlab';
        %LON_val = [-76.5 -72]; LON_val=[-81 -76.5]; LON_val = [-85.5-81];
        %LON_val = [-90 -85.5];

       
     
        
        
        diurnal=0;
        
        if diurnal==1
            LAT_val = [-22 -18]; LON_val = [-76.5 -72];  %one location
        %        for diurnal PDFS 
        else
             LAT_val = [-22 -18]; LON_val = [-76.5 -72; -81 -76.5; -85.5 -81; -90 -85.5];
        end
        
        
        %think there is a 3-4 hour time difference for LT in Chile e.g. LT = UTC-4
        times_of_day = [0 3 6 9 12 15 18 21]; %UTC
%        times_of_day = [15 18 21];  
        

        %for diurnal PDFs replicate the longitude for each time required
        if diurnal==1
            LON_val = repmat(LON_val,[length(times_of_day) 1]);
        end


        gcm_case = 'CTP';
        gcm_case = 'H';
        gcm_case = 'LWP';
        gcm_case = 'Nd';   
        gcm_case = 'CF';           

        clear mean_val

        for idat=1:size(LON_val,1);
            
            if diurnal==1
                iday_gcm2D_inds = [find(gcm_time_UTC==times_of_day(idat))];
            else
                iday_gcm2D_inds=[];
                for igcm_time=1:length(times_of_day)
                    iday_gcm2D_inds2 = [find(gcm_time_UTC==times_of_day(igcm_time))];
                    iday_gcm2D_inds = [iday_gcm2D_inds; iday_gcm2D_inds2];
                end
                iday_gcm2D_inds = sort(iday_gcm2D_inds);
            end
            
  

            ihtot = find ( Plat>=LAT_val(1) & Plat<=LAT_val(2) & Plon>LON_val(idat,1) & Plon<=LON_val(idat,2) ); thresh_str3=['LAT.GTE.' num2str(LAT_val(1)) ' AND LAT.LTE.' num2str(LAT_val(2)) ' AND LON.GT.' num2str(LON_val(idat,1)) ' AND LON.LTE.' num2str(LON_val(idat,2))]; thresh_str2=['LON=' num2str(LON_val(idat,1)) ' to ' num2str(LON_val(idat,2))];
            %these are the lat/lon indices - need to repeat these for each
            %time and calculate a new index for the 3D matrix
            IJ = repmat(ihtot',[length(iday_gcm2D_inds) 1]);
            %K all of the time indices replicated
            K = repmat([iday_gcm2D_inds],[1 length(ihtot)]);
            ihtot = sub2ind([size(gcm_CTP,1) size(gcm_CTP,2)*size(gcm_CTP,3)] , K(:) ,IJ(:));

            ylab = 'Frequency';
            
            if diurnal==1
                thresh_str2 = [num2str(times_of_day(idat)) ' UTC'];
            end

            lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane


            switch gcm_case
                case 'CTP'
                    X_bins = [400:50:1050];
                    Y_data = gcm_CTP(ihtot)/100;
                    xlab='Cloud Top Pressure (hPa)';

                case 'H'
                    X_bins = [0:50:2000];
                    Y_data = (gcm_CTH(ihtot)-gcm_CBH(ihtot));
                    xlab='Cloud Depth (m)';
                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                    
                case 'LWP'
                    X_bins = [0:10:550];
                    Y_data = gcm_lwp_cloud(ihtot);
                    xlab='Cloud mean LWP (g m^{-2})';       
                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot panel
                case 'Nd'
                    X_bins = [0:5:550];
                    Y_data = gcm_Nd_meanliq(ihtot);
                    xlab='Droplet Concentration (cm^{-3})';       
                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot panel  
               case 'CF'
                    X_bins = [-0.02:0.02:1.02];
                    Y_data = gcm_CF_max(ihtot);
                    xlab='Cloud Fraction';  
                    lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                             
            end
            
            
            

            xdat(idat).x = 0.5 * (X_bins(1:end-1) + X_bins(2:end));
            Y_data(isnan(Y_data))=''; %remove NaN values
            ydat(idat).y = ndhistc(Y_data,X_bins); %
            mean_val(idat) = mean(Y_data);
            labs(idat).l=[thresh_str2 ',me=' num2str(mean_val(idat),'%.1f')];



        end
        
        titlenam = ['gcm PDF ' xlab thresh_str];
        
        if diurnal==0
            titlenam = [titlenam ' LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(2))];
        else
            titlenam = [titlenam ' LAT=' num2str(LAT_val(1)) ' to ' num2str(LAT_val(2)) ', LON=' num2str(LON_val(1,1)) ' to ' num2str(LON_val(1,2)) ];
        end
        figname=titlenam;
        savename=[savedir figname];

        
 case 111
     %1D Nd PDFs from different seasons from 2D histograms
     
     plot_seasonal_or_location = 'seasonal';
     plot_seasonal_or_location = 'location';  %iseason=4;
              

     titlenam = ['1D PDF ' tit(1).tit];

        xlims=1;
        xlimits=[0 600];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.
        
          xlab = season_vals_save(1).NdPDF_xlab;
%        xlab=xlabelstr;
        ylab='N data points';
        
        switch plot_seasonal_or_location
            case 'seasonal'


                ylab=season_vals_save(1).pdf_ylab;

                %        y_axis_type = 'log10_matlab';
                nmark=-1;

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;

                for idat=1:length(season_vals_save)


                    %        xdat(idat).x=X_mean; %
                    xdat(idat).x = season_vals_save(idat).NdPDF_xbins;
                    ydat(idat).y = season_vals_save(idat).NdPDF; %
                    labs(idat).l=[season_vals_save(idat).pdflab ' me=' num2str(seasonal_meanX(idat),'%3.1f')];

                end

            case 'location'
                ylab=season_vals_save(1).pdf_ylab;

                %        y_axis_type = 'log10_matlab';
                nmark=-1;

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;

                for idat=1:length(location_vals_save)


                    %        xdat(idat).x=X_mean; %
                    xdat(idat).x = location_vals_save{idat}(iseason).NdPDF_xbins;
                    ydat(idat).y = location_vals_save{idat}(iseason).NdPDF; %
%                    labs(idat).l=[location_vals_save{idat}(iseason).pdflab ' me=' num2str(seasonal_meanX(idat),'%3.1f')];
                    labs(idat).l=[location_vals_save{idat}(1).lon_lab ' ^{o}E'];                    

                end
                
                titlenam = [location_vals_save{idat}(iseason).pdflab ' ' titlenam];


        end
        
  
        
        figname=titlenam;
        savename=[savedir figname];
        
        
        
     case 110
     %plots from different seasons of X_means from 2D histograms
              

        titlenam = ['No, data points ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab = season_vals_save(1).ylabelstr;
%        xlab=xlabelstr;
        xlab='N data points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;


        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Ndatap;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='DJF';
        
        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Ndatap;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='MAM';
        
        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Ndatap;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='JJA';
        
        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Ndatap;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='SON';
   
    
    case 109
     %plots from different seasons of X_means from 2D histograms
              

        titlenam = ['Nd means for ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab = season_vals_save(1).ylabelstr;
%        xlab=xlabelstr;
        xlab='Mean Droplet Concentration (cm^{-3})';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;


        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Xmean;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='DJF';
        
        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Xmean;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='MAM';
        
        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Xmean;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='JJA';
        
        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x = season_vals_save(idat).Xmean;
        ydat(idat).y = season_vals_save(idat).Yvals; %
        labs(idat).l='SON';
        
        
        
    case 108
        
xlab='Nd (cm^{-3})';
%        xlab='Reff (\mum)';        
%        xlab='Phase';      
%        xlab='Ice number concentration (L^{-1})';
%        xlab='Height (m)';
%xlab='Liquid Water Content (g m^{-3})';
%xlab='Ice Water Content (g m^{-3})';
%xlab='Distance from start of profile (km)';
%xlab='Altitude (m)';

% *************************** timeseries flag *****************************
itimeseries = 0; %flag to make it plot a timeseries instead of a profile
% *************************************************************************


profile_set = 'Ascent and 2nd over Oliktok';
profile_set = 'First ascent profile over Oliktok';
%profile_set = 'Partial descent profile over Oliktok';
  profile_set = 'Subsequent profiles over Oliktok';
%  profile_set = 'Subsequent profiles over Oliktok - profile 1';
%  profile_set = 'Subsequent profiles over Oliktok - profile 2';
  profile_set = 'Subsequent profiles over Oliktok - profile 3';  
%  profile_set = 'Porpoise profile on transit to Barrow (1)';
%  profile_set = 'Porpoise profiles on transit to Barrow (2)';
%  profile_set = 'Porpoise profiles on transit to Barrow (3)';
%  profile_set = 'Porpoise profiles on transit to Barrow (4)';
%  profile_set = 'Porpoise profiles on transit to Barrow (5)';
%  profile_set = 'Porpoise profiles on transit to Barrow (6)';
%  profile_set = 'ALL';



%set the type of y-axis
ylab=['Altitude (m)'];
%ylab=['Normalised Height'];

CT_distance = [-150 0]; %metres below and above the height of the max LWC location for averaging

switch profile_set
    case 'Ascent and 2nd over Oliktok'
        
        time_profiles=[20.3 20.42;20.53 20.63];
%                time_profiles=[20.53 20.63];
        CB_height = [292.75; 1000];
        CB_height = [350; 9.9e9];    
        
    case 'First ascent profile over Oliktok'
        time_profiles=[20.336 20.4];
%                time_profiles=[20.53 20.63];
        CB_height = [350]; 
%        CB_height = [320];         
        
        iCT_mean = [20:26];
        
     case 'Partial descent profile over Oliktok'
        time_profiles=[20.53 20.63];
        
        CB_height = [350]; 
        
        
    case 'Subsequent profiles over Oliktok'
        
        time_profiles=[20.86 21.135];
        CB_height = [430;];  
        
    case 'Subsequent profiles over Oliktok - profile 1'
        time_profiles=[20.875 20.965];
        CB_height = [490;];
        
    case 'Subsequent profiles over Oliktok - profile 2'
        time_profiles=[20.965 21.04];
        CB_height = [490;];    
        
    case 'Subsequent profiles over Oliktok - profile 3'
        time_profiles=[21.052 21.135];
        CB_height = [490;]; 
        
    case 'Porpoise profile on transit to Barrow (1)'
        
        time_profiles=[21.135 21.33;21.33 21.435];
        CB_height = [430; 580];  
        
     case 'Porpoise profiles on transit to Barrow (2)'
        
        time_profiles=[21.435 21.4865;21.4865 21.545];
        CB_height = [575;650];      
        
    case 'Porpoise profiles on transit to Barrow (3)'
        time_profiles=[21.545 21.6;21.6 21.67];
        CB_height = [550;700];

    case 'Porpoise profiles on transit to Barrow (4)'
        time_profiles=[21.718 21.77;21.77 21.85];
        CB_height = [450;600];
        
    case 'Porpoise profiles on transit to Barrow (5)'
        time_profiles=[21.85 21.98;21.98 22.045];
        CB_height = [700; 850]; 
        
    case 'Porpoise profiles on transit to Barrow (6)'
        time_profiles=[22.045 22.12;];
        CB_height = [575]; 
        
     case 'ALL'
        time_profiles=[mpace_time(1) mpace_time(end);];
        CB_height = [575];
        lwidth=0.1;
        
end
        
clear inds xline time_profs
for itime=1:size(time_profiles,1);
    ii = find(mpace_time>=time_profiles(itime,1) & mpace_time<=time_profiles(itime,2));
    inds(itime,1:length(ii)) = ii;

    mean_lon(itime) = mean(mpace_lon(ii));
    lon_range(itime,:) = [min(mpace_lon(ii)) max(mpace_lon(ii))];
     time_profs(itime).t=mpace_time(ii);

    if CB_height(itime)<9e9
        iCB=findheight_nearest(mpace_alt,CB_height(itime));
        tCB=273.15+mpace_temp(iCB);
        pCB=100*mpace_press(iCB);
        %             pressure=100*[mpace_press(iCB):-1:min(mpace_press(ii))];
        pressure=100*mpace_press(ii);


        for iad=1:length(pressure)
            [xline(itime).x(iad),Tad(iad)]=adLWC_PaulLawson_simple(pCB,tCB,pressure(iad));
        end

        yline(itime).y = mpace_alt(ii);               
      
    else
        xline(itime).x=[];
        yline(itime).y=[];
    end

end
%         inds=inds(:);
%         inds(inds==0)=[];
        
        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=1;
        zmin=0;
        zmax=1800;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        X=0.5*(LONS(1:end-1)+LONS(2:end));
        Y=0.5*(LATS(1:end-1)+LATS(2:end));
        
        thresh_str='';
        %                ihtot = 1:length(mpace_lat_mapped);
        ihtot=[];
        
%[mpace_time,mpace_phase,mpace_temp,mpace_height,mpace_cwc,mpace_lwc,mpace_iwc,mpace_rew,mpace_rei,mpace_Nd,mpace_Ni]         
                

       

                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];

                tol=0.1;      
               
               idat2=-1;
               for idat=1:size(time_profiles,1)
                   
                   idat2=idat2+2;
            
                   ii=inds(idat,:);
                   ii(ii==0)=[];
                   
                   switch xlab
                       case 'Liquid Water Content (g m^{-3})'
                           xdat(idat2).x = mpace_lwc(ii);
                           
                           xlims=1;
                           xlimits=[0 0.6];
                           
                       case 'Ice Water Content (g m^{-3})'
                           xdat(idat2).x = mpace_iwc(ii);
                           
                       case 'Nd (cm^{-3})'
                           xdat(idat2).x = mpace_Nd(ii);
                           xline(idat).x=[]; %switch off the adibatic LWC profile
                           
                             xlims=1;
                             xlimits=[0 150];
                             
                       case 'Distance from start of profile (km)'
                           start_lat = mpace_lat(ii(1))*ones(size(mpace_lat(ii)));
                           start_lon = mpace_lon(ii(1))*ones(size(mpace_lat(ii)));
                           xdat(idat2).x =distlatlon(start_lat,start_lon,mpace_lat(ii),mpace_lon(ii));
                           xdat(idat2).x(1) = 0; %set first value to zero as gives imaginary answer
                           xline(idat).x=[]; %switch off the adibatic LWC profile
                           
                           
                            xlims=1;
                            xlimits=[0 max(xdat(idat2).x)*1.15];
                           
                       case 'Altitude (m)';
                            xdat(idat2).x = mpace_alt(ii);
                           xline(idat).x=[]; %switch off the adibatic LWC profile
                           
                             xlims=0;
                             xlimits=[0 150];
                             
                       case 'Reff (\mum)';
                            xdat(idat2).x = mpace_rew(ii);
                           xline(idat).x=[]; %switch off the adibatic LWC profile
                           
                             xlims=1;
                             xlimits=[0 20];      
                   end
                   
                   if itimeseries==1 | strcmp(profile_set,'ALL')==1
                           xline(idat).x=[]; %switch off the adibatic LWC profile
                   end
                    
                    
                xdat(idat2).x(ihtot) = NaN;
                
                switch ylab
                    case 'Altitude (m)'
                        ydat(idat2).y = mpace_alt(ii);
                    case 'Normalised Height'
                        ydat(idat2).y = mpace_height(ii);
                end
                
                lwc = mpace_cwc(ii);
                [maxval,iymax] = max(lwc);
                maxheight = ydat(idat2).y(iymax);
                
                
                
                 iCT_save(idat2).dat = find(ydat(idat2).y>=maxheight+CT_distance(1) & ydat(idat2).y<=maxheight+CT_distance(2));
                
%   calculating the overall profile mean and the cloud top mean (using defined indices above)
%                labs(idat2).l=['LON=' num2str(mean_lon(idat))];
                mean_vals(idat2).me = meanNoNan(xdat(idat2).x,2);                                
                mean_vals(idat2).me2 = meanNoNan(xdat(idat2).x(iCT_save(idat2).dat),2);
%                mean_vals(idat2).me2=99;
                
                labs(idat2).l=['LON=' num2str(lon_range(idat,1),'%.2f') ' to ' num2str(lon_range(idat,2),'%.2f') ' me=' num2str(mean_vals(idat2).me) ' me2=' num2str(mean_vals(idat2).me2) ' for ' num2str(CT_distance(1)) ' to ' num2str(CT_distance(2)) ' m'];                
                nmark(idat2)=-1;
                
                if length(xline(idat).x)>0

                    xdat(idat2+1).x = xline(idat).x;
                    ydat(idat2+1).y = yline(idat).y;
                    labs(idat2+1).l=['Adiabatic'];
                    nmark(idat2+1)=0;
                    ismooth_y(idat2+1)=0;

                    

                else

                    idat2=idat2-1;

                end
                

                smooth_mode='mean';
                Nsmooth_window = 6;
                
               end
                
%                thresh_str=[thresh_str 'within ' num2str(tol*100) ' %'];                
                               
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
                
                

         title_str = 'MPACE flight profiles, ';
         param_str = xlab;

         if itimeseries==1
             %swap the y-axis to the x-axis
             for idat=1:length(xdat)
                 ydat(idat).y = xdat(idat).x;
                 xdat(idat).x = time_profs(idat).t
             end
             
              izlim=xlims;              
              zmin=xlimits(1);
              zmax=xlimits(2);
        
              xlims=0;
%              xlimits=1000*[0 0.025];


               ylab = xlab;
               xlab = 'Time (UTC)';
               
               lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
               
               title_str = 'MPACE timeseries, ';
               param_str = ylab;
             
         end

             
         titlenam = [title_str profile_set ' ' param_str];
         titlenam=remove_character(titlenam,'_','-');

         figname=titlenam;
         savename=[savedir figname];

         for idat=1:length(xdat)

             if ismooth_y(idat)==1
                 labs(idat).l=[labs(idat).l ' smoothed ' num2str(Nsmooth_window) ' steps'];
                 titlenam=[titlenam ' smoothed ' num2str(Nsmooth_window) ' steps'];
             end

         end
         
         
        
    case 107
        
        ylab='Nd (cm^{-3})';
%        ylab='Reff (\mum)';        
%        ylab='Phase';
%        ylab='Cloud Water Content (g m^{-3})';
%        ylab='Optical Depth';
%         ylab = 'Condensation rate (g m^{-4})';  
%         ylab = 'Cloud Top Tempeature (K)';
%        ylab='Ice Water Content (g m^{-3})';
%        ylab='Ice number concentration (L^{-1})';
%        ylab='Height (m)';
%        ylab='Height ratio within cloud (m)';
%         ylab = 'Wind speed (m/s)';
%         ylab='Wind dir (degrees)';
%         ylab='Cloud Fraction';
%         ylab='Homogeneity Factor Wood';
%         ylab='Homogeneity Factor Cahalan';
         ylab='Nd uncertainty (%)';
%         ylab='Nd mean comparison (cm^{-3})';


map_flight_track_onto_L2
%now should use lat_new and lon_new as the new mpace coords


        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        X=0.5*(LONS(1:end-1)+LONS(2:end));
        Y=0.5*(LATS(1:end-1)+LATS(2:end));
        
        thresh_str=thresh_str_mock_L3;
        %                ihtot = 1:length(mpace_lat_mapped);
        ihtot=[];
        
%[mpace_time,mpace_phase,mpace_temp,mpace_height,mpace_cwc,mpace_lwc,mpace_iwc,mpace_rew,mpace_rei,mpace_Nd,mpace_Ni]         
                

        switch ylab
            case 'Height ratio within cloud (m)'
                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];

                idat=idat+1;

                ibase=find(mpace_height>-tol & mpace_height<tol);                
                xdat(idat).x = mpace_lon;
                ydat(idat).y = mpace_height;
                labs(idat).l='Height ratio within cloud';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
                
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
                
                
            case 'Height (m)'  %cloud base and top heights
              
                
                % 0=no cloud, 1=ice, 2=mixed, 3=water
                thresh_H=0.75;
                thresh_H=-0.2;
                
%                dat2D = meanNoNan(W_all,3);

%                ihtot = find(~ (mpace_height<=thresh_H) ); thresh_str = ['for H.LTE.' num2str(thresh_H)]; 
%                ihtot = find(~ (mpace_height>=thresh_H) ); thresh_str = ['for H.GTE.' num2str(thresh_H)]; 

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];

                mpace_cb_and_ct_height_calc
                %tol=0.1;   %tol now set in the function above   
                
                idat=idat+1;
                xdat(idat).x=cb_lons;               
                ydat(idat).y = cloud_base_height;
                labs(idat).l='Cloud base';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
                idat=idat+1;             
                xdat(idat).x=ct_lons;            
                ydat(idat).y = cloud_top_height;
                labs(idat).l='Cloud top';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
                thresh_str=[thresh_str 'within ' num2str(tol*100) ' % '];   
                

                
                %satellite cloud thickness + a baseline height
                
%                dat2D = meanNoNan(H_all,3);
                dat2D = H_meanTauReff;
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;                
                ydat(idat).y = interp2(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped);
                labs(idat).l=['Satellite Cloud Depth + aircraft cloud base height'];
                
                %                 add_base_height = 600;
%                 ydat(idat).y = interp2(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped) + add_base_height;
%                 labs(idat).l=['Satellite Cloud Depth + ' num2str(add_base_height) 'm'];


                % add the aircraft cloud base height to the satellite cloud depth                   
                xdat(idat).x=xdat(idat).x(ibase);
                ydat(idat).y=ydat(idat).y(ibase) + cloud_base_height;
                
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;                
                ydat(idat).y = interp2(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped);
                labs(idat).l=['Satellite Cloud Depth'];

                
               

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
                
                
                case 'Ice number concentration (L^{-1})'
                % 0=no cloud, 1=ice, 2=mixed, 3=water
                thresh_H=0.75;
                thresh_H=-0.2;
                
%                dat2D = meanNoNan(W_all,3);

                ihtot = find(~ (mpace_height<=thresh_H) ); thresh_str = ['for H.LTE.' num2str(thresh_H)]; 
%                ihtot = find(~ (mpace_height>=thresh_H) ); thresh_str = ['for H.GTE.' num2str(thresh_H)]; 

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];
                              
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = mpace_Ni;
                labs(idat).l='Aircraft';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
                
               

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
                
            case 'Ice Water Content (g m^{-3})'
                % 0=no cloud, 1=ice, 2=mixed, 3=water
                thresh_H=0.9;
                
%                dat2D = meanNoNan(W_all,3);
                
%                ihtot = find(~ (mpace_height>=thresh_H) ); thresh_str = ['for H.GTE' num2str(thresh_H)]; 

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];
                              
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = mpace_iwc;
                labs(idat).l='Aircraft';
                ismooth_y(idat)=1;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
               

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            
            case 'Optical Depth'

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];                                             

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                thresh_H=[-0.90 99]; %lwidth=0; %0.0;
            thresh_H=[0.75 1.1]; lwidth=0; %remove the line
                ihtot = find(~ (mpace_height>=thresh_H(1) & mpace_height<=thresh_H(2)) ); thresh_str = [' for H.GTE' num2str(thresh_H(1)) ' and H.LTE' num2str(thresh_H(2))];
                  
                
%                cw_aircraft=8e-4;  thresh_str=[thresh_str 'using constant cw=' num2str(cw_aircraft)];
%                cw_aircraft=interp2(X,Y,meanNoNan(cw_all,3),mpace_lon_mapped,mpace_lat_mapped)*1000; thresh_str=[thresh_str 'using varying cw']
                cw_aircraft=0.8*1e3*map_flight_onto_sat(X,Y,cw_meanTauReff,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                
                W_aircraft = 0.5 * mpace_cwc.^2 ./ cw_aircraft /1000;
              
                
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                rhoL=1000;
                ydat(idat).y = W_aircraft.* 9/5 ./(mpace_rew/1e6) / rhoL;
                labs(idat).l='Aircraft Tau (using aircraft Re, cw*0.8)';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                rhoL=1000;
                Re_sat = 1e-6*map_flight_onto_sat(X,Y,meanRe_mockL3,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                ydat(idat).y = W_aircraft.* 9/5 ./Re_sat / rhoL;
                labs(idat).l='Aircraft Tau (using satellite Re, cw*0.8)';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;

                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,meanTau_mockL3,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['Satellite Tau'];
                
                
                
            case 'Cloud Top Tempeature (K)'

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];                                             

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;                
                ydat(idat).y = interp2(X,Y,meanCTT_mockL3,mpace_lon_mapped,mpace_lat_mapped);
                labs(idat).l='Cloud Top Temperature';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;

         
                
                
              case 'Condensation rate (g m^{-4})'

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];                                             

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
                 %cw_aircraft is the cw from the satellite CTT interpolated onto aircraft lon grid              
%                cw_aircraft=interp2(X,Y,meanNoNan(cw_all,3),mpace_lon_mapped,mpace_lat_mapped)*1000; %thresh_str=[thresh_str 'using varying cw']
                cw_aircraft=interp2(X,Y,meanNoNan(cw_meanTauReff,3),mpace_lon_mapped,mpace_lat_mapped)*1000; %thresh_str=[thresh_str 'using varying cw']                
               
                % W_aircraft = 0.5 * mpace_cwc.^2 ./ cw_aircraft /1000;
                
              %calculate teh aircraft cloud base and heights (approx and poor
              %resolution)
                mpace_cb_and_ct_height_calc
                %tol=0.1;   %tol now set in the function above   
                
              
                
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                rhoL=1000;
                ydat(idat).y = cw_aircraft;
                labs(idat).l='cw';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                rhoL=1000;
                ydat(idat).y = cw_aircraft;
                labs(idat).l='cw';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
             
                
            
            case 'Cloud Water Content (g m^{-3})'  
                % 0=no cloud, 1=ice, 2=mixed, 3=water
                thresh_H=0.9;
                
%                dat2D = meanNoNan(W_all,3); %mean of all individual W values
                dat2D = W_meanTauReff; %value using the mean tau and reff for the lat/lon cell
                %makes very little difference to the LWP - perhaps more of
                %a difference for Nd
                
                
%                ihtot = find(~ (mpace_height>=thresh_H) ); thresh_str = ['for H.GTE' num2str(thresh_H)]; 

                xlab=['Longitude'];
                
                idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];

               cw_aircraft = 1e3*map_flight_onto_sat(X,Y,cw_meanTauReff,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
               
                % W_aircraft = 0.5 * mpace_cwc.^2 ./ cw_aircraft /1000;
                
              %calculate teh aircraft cloud base and heights (approx and poor
              %resolution)
                mpace_cb_and_ct_height_calc
                %tol=0.1;   %tol now set in the function above   
                
                %now bin the longitude values from the aircraft cloud base
                %into regular lon bins
                LON_bins = [-157:0.5:-148];
                [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals,meanvals_ignore_zeros,std_dev,mean_error]=bin_data2(cb_lons,cloud_base_height,LON_bins,cloud_base_height);
                x_int = mid_points;
                y_int = meanvals;
                x_int(isnan(meanvals))=[];
                y_int(isnan(meanvals))=[];                
                cb_int = interp1(x_int,y_int,mpace_lon);
                
               [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals,meanvals_ignore_zeros,std_dev,mean_error]=bin_data2(ct_lons,cloud_top_height,LON_bins,cloud_top_height);
                x_int = mid_points;
                y_int = meanvals;
                x_int(isnan(meanvals))=[];
                y_int(isnan(meanvals))=[];                
                ct_int = interp1(x_int,y_int,mpace_lon);      


                              
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = mpace_cwc;
                labs(idat).l='Aircraft';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;
                
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 W_aircraft = 0.5 * mpace_cwc.^2 / 8e-4 /1000;
%                 ydat(idat).y =  W_aircraft *2;
%                 labs(idat).l='Aircraft LWP (from cw assumption) x2 (kg/m2)';
%                 ismooth_y(idat)=0;
%                 smooth_mode='mean';
%                 Nsmooth_window = 6;
%                 
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 rhoL=1000;
%                 Re_sat = 1e-6*interp2(X,Y,meanRe_mockL3,mpace_lon_mapped,mpace_lat_mapped);
%                 ydat(idat).y = W_aircraft.* 9/5 ./Re_sat / rhoL  /50;
%                 labs(idat).l='Aircraft Tau (from cw assumption) /50';
%                 ismooth_y(idat)=0;
%                 smooth_mode='mean';
%                 Nsmooth_window = 6;
%                 
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = 2*interp2(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped);
%                 labs(idat).l=['Satellite LWPx2'];
%                 
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = interp2(X,Y,meanTau_mockL3,mpace_lon_mapped,mpace_lat_mapped)/50;
%                 labs(idat).l=['Satellite Tau/50'];
    
re_case = '1.6 \mum';
re_case = '2.1 \mum';
re_case = '3.7 \mum';

switch re_case
    case '1.6 \mum'

        %                dat2D = meanNoNan(H_all,3);
        dat2D = H_meanTauReff_16;
        H_1D = map_flight_onto_sat(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
        thresh_str=re_case;
    case '2.1 \mum'
        dat2D = H_meanTauReff;
        H_1D = map_flight_onto_sat(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
        thresh_str=re_case;
    case '3.7 \mum'
        dat2D = H_meanTauReff_37;
        H_1D = map_flight_onto_sat(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
        thresh_str=re_case;
end
                
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;                
%                 ydat(idat).y = H_1D/1000;
%                 labs(idat).l=['Satellite H/1000'];
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                cw_plot = 8e-7;
%                cw_plot = cw; %cw should have been set in MODIS_N_H_calc is in kg/m4 so *1000
%                ydat(idat).y = H_1D*cw_plot*1000;    
                ydat(idat).y = 0.8*H_1D.*cw_aircraft;                    
%                labs(idat).l=['Satellite LWC (cw*H), cw=' num2str(cw_plot*1000) ' g/m^{4}'];
                labs(idat).l=['Satellite LWC (0.8*cw*H), CTT dependent cw, g/m^{4}'];                
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;               
                ydat(idat).y = 0.8*cw_aircraft.*(ct_int-cb_int);
                labs(idat).l=['cw*aircraft H*0.8'];
                
                

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
               
                
                
            case 'Phase'  
                % 0=no cloud, 1=ice, 2=mixed, 3=water
                thresh_H=0.75;
                
%                ihtot = find(~ (mpace_height>=thresh_H) ); thresh_str = ['for H.GTE' num2str(thresh_H)]; 

                xlab=['Longitude'];
                idat=0;
                
              
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = mpace_phase;
                labs(idat).l='Aircraft';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
               
                
            case  'Nd (cm^{-3})'

               

                xlab=['Longitude'];
                
                dat2D_pixel=meanNd_mockL3; %mean of all individual Nd values within lat/lon cell
                dat2D=N_meanTauReff; %makes mean of tau and reff for cell and use this for Nd
%                dat2D_16=meanNd_mockL3_16; %mean of all individual Nd values within lat/lon cell
                dat2D_16=N_meanTauReff_16; %makes mean of tau and reff for cell
%                dat2D_37=meanNd_mockL3_37; %mean of all individual Nd values within lat/lon cell
                dat2D_37=N_meanTauReff_37; %makes mean of tau and reff for cell
                
                
                 re_case = '1.6 \mum';
%                re_case = '2.1 \mum';
%                re_case = '3.7 \mum';   
                
                thresh_str = [re_case ' '];
                
            
            thresh_H=[-0.90 99]; %lwidth=0; %0.0;
            thresh_H=[0.75 1.1]; lwidth=0; %remove the line
%            thresh_H=[0.9 1.1]; lwidth=0; %remove the line            

            LWC_thresh = 0.1;
                

                
                %make an Reff vs lon array for the aircraft obs
                %use the values near cloud top
                thresh_H2=0.95; %0.0;
                
                ihtot = find(~ (mpace_height>=thresh_H2) ); 
                x_reff = mpace_lon;
                %make other values NaN
                x_reff(ihtot) = NaN;
                reff_lon = mpace_rew;
                
                %remove the very small point as is probably an erroneous
                %measurement. And remove the zeros too.
                ihtot = find(reff_lon<8); 
                x_reff(ihtot) = NaN;
                icut = find(isnan(x_reff)==1);
                x_reff(icut)=[];
                reff_lon(icut)=[];
                
                reff_lon2 = interp1(x_reff,reff_lon,mpace_lon,'linear','extrap');
                
            %optical depth from satellite                    
                
            tau_mpace = map_flight_onto_sat(X,Y,meanTau_mockL3,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
            CTT_mpace = map_flight_onto_sat(X,Y,meanCTT_mockL3,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                       
            [mpace_Nd_satRe,H_mpace,W_mpace,k_mapce,Q_mpace,cw_mpace]=MODIS_N_H_func(tau_mpace,reff_lon2*1e-6,Wflag,NaN,CTT_mpace);
                
            


                ihtot = find(~ (mpace_height>=thresh_H(1) & mpace_height<=thresh_H(2) & mpace_cwc>LWC_thresh) ); thresh_str = [thresh_str 'for H.GTE' num2str(thresh_H(1)) ' and H.LTE' num2str(thresh_H(2)) ' and LWC.GTE.' num2str(LWC_thresh) ' g m^{-3}']; 

                
% %                tol=0.05; tol2=0.05; lwidth=0;
%                 tol=2; tol2=2.1; lwidth=0;                
% %                tol=99;  tol2=99;
%                 ihtot = find( mpace_height<=1-tol | mpace_height>=1+tol2 ); thresh_str = [thresh_str 'within ' num2str(tol*100) '% of cloud top']; 
%                 
                

idat=0;


                idat=idat+1;                
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = mpace_Nd;
                ismooth_y(idat)=0;
                labs(idat).l='Aircraft';
                

                
                switch re_case
                    case '1.6 \mum'

                        idat=idat+1;
                        xdat(idat).x = mpace_lon;
                        xdat(idat).x(ihtot) = NaN;
                        ydat(idat).y = map_flight_onto_sat(X,Y,dat2D_16,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                        labs(idat).l=['Satellite 1.6 \mum'];

                    case '2.1 \mum';

                        idat=idat+1;
                        xdat(idat).x = mpace_lon;
                        xdat(idat).x(ihtot) = NaN;
                        ydat(idat).y = map_flight_onto_sat(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                        labs(idat).l=['Satellite 2.1 \mum'];

                    case '3.7 \mum';

                        idat=idat+1;
                        xdat(idat).x = mpace_lon;
                        xdat(idat).x(ihtot) = NaN;
                        ydat(idat).y = map_flight_onto_sat(X,Y,dat2D_37,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                        labs(idat).l=['Satellite 3.7 \mum'];

                end
%                 
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = map_flight_onto_sat(X,Y,dat2D_pixel,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
%                 labs(idat).l=['Average of pixel values 2.1 \mum'];

                
%                 idat=idat+1;                
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = mpace_Nd_satRe;
%                 ismooth_y(idat)=0;
%                 labs(idat).l='Aircraft Re, Sat Tau';
%                 
%                 idat=idat+1;                
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = mpace_Nd_satRe*sqrt(0.8);
%                 ismooth_y(idat)=0;
%                 labs(idat).l='sqrt(0.8)*Aircraft Re, Sat Tau';                
                



                smooth_mode='mean';
                Nsmooth_window = 6;

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case  'Reff (\mum)'
                

                Re_plot ='all';
                Re_plot = '1.6 \mum';
%                Re_plot = '2.1 \mum';
%                Re_plot = '3.7 \mum';

                thresh_str = Re_plot;



                thresh_H=[-0.90 99]; %lwidth=0; %0.0;
                thresh_H=[0.90 1.1]; lwidth=0; %remove the lines    
                LWC_thresh = 0.1;
                
                ihtot = find(~ (mpace_height>=thresh_H(1) & mpace_height<=thresh_H(2) & mpace_cwc>LWC_thresh) ); thresh_str = [thresh_str 'for H.GTE' num2str(thresh_H(1)) ' and H.LTE' num2str(thresh_H(2)) ' and LWC.GTE.' num2str(LWC_thresh) ' g m^{-3}']; 
                
                thresh_str=[thresh_str ' NO Re removal'];

                xlab=['Longitude'];
idat=0;



                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = mpace_rew;
                labs(idat).l='Aircraft';
                ismooth_y(idat)=0;
                smooth_mode='mean';
                Nsmooth_window = 6;

switch Re_plot
    case {'1.6 \mum','all'}

                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;                
                ydat(idat).y = map_flight_onto_sat(X,Y,meanRe_mockL3_16,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['Satellite 1.6 \mum'];
    case {'2.1 \mum','all'}
%                 
                idat=idat+1;                
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,meanRe_mockL3,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['Satellite 2.1 \mum'];
                
    case {'3.7 \mum','all'}
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,meanRe_mockL3_37,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['Satellite 3.7 \mum'];

end
                


                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                
            case  'Wind speed (m/s)'

                thresh_H=[0 0.8];
                
                ihtot = find(~ (mpace_height>=thresh_H(1) & mpace_height<thresh_H(2)) ); thresh_str = ['for H.GTE' num2str(thresh_H(1)) 'AND.H.LT' num2str(thresh_H(2))]; 

                xlab=['Longitude'];

                xdat(1).x = mpace_lon;
                xdat(1).x(ihtot) = NaN;
                ydat(1).y = mpace_wind;
                labs(1).l=['Wind speed'];

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                                
            case  'Wind dir (degrees)'
             
                thresh_H=[0 0.8];
%                thresh_H=[-99 99];
                               
                ihtot = find(~ (mpace_height>=thresh_H(1) & mpace_height<thresh_H(2)) ); thresh_str = ['for H.GTE' num2str(thresh_H(1)) 'AND.H.LT' num2str(thresh_H(2))]; 

                xlab=['Longitude'];

                xdat(1).x = mpace_lon;
                xdat(1).x(ihtot) = NaN;
                ydat(1).y = mpace_winddir;
                labs(1).l=['Wind dir'];

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                                
                
                
            case  'Cloud Fraction'

               

                xlab=['Longitude'];
                
                idat=0;
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = interp2(X,Y,CF_mockL3,mpace_lon_mapped,mpace_lat_mapped);
                labs(idat).l=['Satellite CF'];
                
            

                
              
                



                smooth_mode='mean';
                Nsmooth_window = 0;

                lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                
            case  'Homogeneity Factor Wood'
                xlab=['Longitude'];
                thresh_str = mockL3_box_type_str;                
                
                Hom = (meanW_mockL3 ./ stdW_mockL3).^2;

                idat=0;

                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,Hom,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['LWP Homogeneity factor, Wood (2006)'];
                
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = interp2(X,Y,(logW_mockL3 ./ meanW_mockL3),mpace_lon_mapped,mpace_lat_mapped);
%                 ydat(idat).y = min(ydat(1).y) + (ydat(2).y-min(ydat(2).y))*(max(ydat(1).y)-min(ydat(1).y));
%                 labs(idat).l=['LWP Homogeneity factor, Cahalan (1994)'];
                
                                             
                smooth_mode='mean';
                Nsmooth_window = 0;

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                                
            case  'Homogeneity Factor Cahalan'

                xlab=['Longitude, LATstep=' num2str(LAT_step)];
                

                Hom = (logW_mockL3 ./ meanW_mockL3);

                idat=0;

                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,Hom,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['LWP Homogeneity factor, Cahalan (1994)'];

                                             
                smooth_mode='mean';
                Nsmooth_window = 0;

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                  
                
            case 'Nd uncertainty (%)'
                xlab=['Longitude, LATstep=' num2str(LAT_step) ' No screening'];

                idat=0;

                dat2D = meanNd_un_mockL3;
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['Average of pixel values'];
                
                dat2D = sqrt( (0.5*meanTau_un_mockL3).^2 + (-5/2*meanRe_un_mockL3).^2 );
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = map_flight_onto_sat(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace);
                labs(idat).l=['Using mean Tau & Re un'];
                
                

                                             
                smooth_mode='mean';
                Nsmooth_window = 0;

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                  
              case 'Nd mean comparison (cm^{-3})'
                xlab=['Longitude, LATstep=' num2str(LAT_step) ' No screening'];
                

               

                idat=0;

                dat2D = meanNd_mockL3; %mean taken after Nd calculated
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = interp2(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped);
                labs(idat).l=['Average of pixel values'];
                
                
                dat2D=N_meanTauReff; %makes mean of tau and reff for cell and use this for Nd
                
                idat=idat+1;
                xdat(idat).x = mpace_lon;
                xdat(idat).x(ihtot) = NaN;
                ydat(idat).y = interp2(X,Y,dat2D,mpace_lon_mapped,mpace_lat_mapped);
                labs(idat).l=['Using mean Tau & Re'];
                
%                 
%                 dat2D=N(row_L2_inds,col_L2_inds); %1km pixel values                                
%                 X2=Plat;
%                 Y2=Plon;
%                 
%                 idat=idat+1;
%                 xdat(idat).x = mpace_lon;
%                 xdat(idat).x(ihtot) = NaN;
%                 ydat(idat).y = interp2(X2,Y2,dat2D,mpace_lon_mapped,mpace_lat_mapped);
%                 labs(idat).l=['1km pixel values'];
                
                

                                             
                smooth_mode='mean';
                Nsmooth_window = 0;

                lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                      
                
                
         end
         
         
         titlenam = [ylab ' from ' aq_terr_str ' ' datestr(minALL(scantime_matlab),'dd mmm HHMM') ' along MPACE flight path ' thresh_str];
         titlenam=remove_character(titlenam,'_','-');
         
         figname=titlenam;
         savename=[savedir figname];
         
         for idat=1:length(xdat)

             if ismooth_y(idat)==1
                 labs(idat).l=[labs(idat).l ' smoothed ' num2str(Nsmooth_window) ' steps'];
%                 titlenam=[titlenam ' smoothed ' num2str(Nsmooth_window) ' steps'];
             end

         end

        
        
    case 106
        %Nd distribution from gaussian Reff and Tau distributions

        %create data from a normal PDF
        Ntot=1e4; %actual number will be slightly different due to rounding
        Nd_bins2=[0:10:1000];
        
        var_for_norm = 'Reff';
        var_for_norm = 'Tau';
        
        switch var_for_norm
            case 'Reff'
                bins = [0:0.01:30]; %

                mu_vals=[12.5 16.5]; %Reff
                sigma_vals=[3];
                sigma_vals=[1];
                other_vals=[15 25]; %Tau vals

            case 'Tau'
                bins = [0:0.01:30]; %
                mu_vals=[15 25]; %of tau
                sigma_vals=[7.5];
%                sigma_vals=[15];
                other_vals=[12.5 16.5]; %Reff vals

        end
        
        
        %%%
        
        xlims=0;
        xlimits=[0 1000];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Frequency';        
        xlab= 'Nd (cm^{-3})';
      
        

        idat=0;
for imu=1:length(mu_vals)
    for isigma=1:length(sigma_vals)
        for iother=1:length(other_vals)
        idat=idat+1;
        
                     
        sigma=sigma_vals(isigma);
        mu=mu_vals(imu);
        val=other_vals(iother);



        fpdf=normpdf(bins,mu,sigma);

        Npoints = round(Ntot*fpdf);
        Ntot_actual = sum(Npoints);
        
        dat_norm = NaN*ones([1 Ntot_actual]);
        iloc=1;
        for i=1:length(Npoints)
            dat_norm(iloc:iloc+Npoints(i)-1)=bins(i);
            iloc=iloc+Npoints(i);
        end
        
        switch var_for_norm
            case 'Reff'
                [N,H,W,k,Q,cw]=MODIS_N_H_func(val,dat_norm*1e-6,'calc',0);
            case 'Tau'
                [N,H,W,k,Q,cw]=MODIS_N_H_func(dat_norm,val*1e-6,'calc',0);
        end
        
        
        
        N(N>Nd_bins2(end))=Nd_bins2(end-1);
        N(isnan(N))=Nd_bins2(end-1);
        Npdf=ndHistc(squeeze(N)', Nd_bins2);
        
      
        
          xdat(idat).x = 0.5*(Nd_bins2(1:end-1)+Nd_bins2(2:end));
          ydat(idat).y = Npdf/Ntot_actual;
          
            [maxval,imax]=max(Npdf);
            mode_val = xdat(idat).x(imax);
            mean_val = mean(N);
            
            mean_val/mode_val
        
        
          labs(idat).l=['Mu=' num2str(mu) ',sig=' num2str(sigma) ',val=' num2str(val) ',mode=' num2str(mode_val) ',mean=' num2str(mean_val)];
          
        end
    end
end

            
            xlims=1;
            xlimits=[0 500];






        lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        titlenam = [xlab ' PDF for ' var_for_norm ' gaussian'];
        
        figname=titlenam;
%        savename=remove_character(figname,'\','');
%        savename=remove_character(savename,'/','');        

%%%
        
        




    case 105
        %1D PDFs for lat/lon cells
        
        
        

        xlims=0;
        xlimits=[0 1000];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Frequency';
        
        xlab= 'Nd (cm^{-3})';
%        xlab= 'Tau';
%        xlab= 'LWP (g m^{-2})';
        xlab= 'Re (\mum)';  i_diff_re=1; %flag to plot the different wavelengths instead of diff positions
%        xlab = 'Re uncertainty (%)';
%        xlab = 'Tau uncertainty (%)';
%       xlab = 'Nd uncertainty (%)';        
        
        LAT_val=[73 72 71 72];
        LON_val=[-160 -153 -160 -163];
        
        LAT_val=[71 71 71];
        LON_val=[-153 -153 -153];
        
        if i_diff_re==0
            nplot=length(LAT_val);
        else
            nplot=3; %for the 3 different wavelengths
        end
        
        clear ilats ilons
        
        for idat=1:nplot

            ilats(idat)=findheight_nearest(LATS,LAT_val(idat));
            ilons(idat)=findheight_nearest(LONS,LON_val(idat));
            
            ilat=ilats(idat);
            ilon=ilons(idat);

            switch xlab
                case 'Nd (cm^{-3})'
                    Xbins = Nd_bins; %

                    ydat(idat).y=NdPDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                    labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];
                    
                     xlims=1;
                     xlimits=[0 600];
                    
                case 'Tau'
                    Xbins = Tau_bins; %

                    ydat(idat).y=TauPDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                    labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];  
                    
                    ismooth_x(idat)=0;
                    ismooth_y(idat)=1;                    
                    smooth_mode='sum';
                    smooth_mode='mean';
                    Nsmooth_window = 30;
                    
                case 'LWP (g m^{-2})'
                     Xbins = W_bins; %

                    ydat(idat).y=WPDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                    labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];  
                    
                     xlims=1;
                     xlimits=[0 1000];
                     

                    ismooth_y(idat)=1;                    
                    smooth_mode='sum';
                    smooth_mode='mean';
                    Nsmooth_window = 6;
                    
                case 'Re (\mum)'
                     Xbins = Re_bins; %
                     
                     if i_diff_re==0

                         ydat(idat).y=RePDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                         labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];
                         ismooth_y(idat)=1;

                     else
                         switch idat
                             case 1
                                 ydat(idat).y=RePDF_mockL3_16(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                                 labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon)) ', 1.6 \mum'];
                                 ismooth_y(idat)=1;


                             case 2
                                 ydat(idat).y=RePDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                                 labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon)) ', 2.1 \mum'];
                                 ismooth_y(idat)=1;

                             case 3
                                 ydat(idat).y=RePDF_mockL3_37(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                                 labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon)) ', 3.7 \mum'];
                                 ismooth_y(idat)=1;
                         end


                     end
                    
                     xlims=1;
                     xlimits=[0 25];
                     

                  
                    smooth_mode='sum';
                    smooth_mode='mean';
                    Nsmooth_window = 6;
                     
                case 'Re uncertainty (%)'
                    Xbins = Re_un_bins; %

                    ydat(idat).y=Re_un_PDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                    labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];  
                    
                     xlims=1;
                     xlimits=[0 20];
                     
                 case 'Tau uncertainty (%)'
                    Xbins = Tau_un_bins; %

                    ydat(idat).y=Tau_un_PDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                    labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];  
                    
                     xlims=1;
                     xlimits=[0 150];
                     
                    ismooth_y(idat)=1;                    
                    smooth_mode='sum';
                    smooth_mode='mean';
                    Nsmooth_window = 6;
                    
                     
                case 'Nd uncertainty (%)'
                    Xbins = Nd_un_bins; %

                    ydat(idat).y=Nd_un_PDF_mockL3(ilat,ilon,:)./Np_mockL3(ilat,ilon); %
                    labs(idat).l=[num2str(LATS(ilat)) ', ' num2str(LONS(ilon))];  
                    
                     xlims=1;
                     xlimits=[10 80];
            end



            xdat(idat).x = 0.5 * ( Xbins(1:end-1)+Xbins(2:end) ); %
            ydat(idat).y = squeeze(ydat(idat).y);

        end




        lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        titlenam = [xlab ' PDF for lat/lon cells'];
        
        figname=titlenam;
%        savename=remove_character(figname,'\','');
%        savename=remove_character(savename,'/','');        

        
        
    case 104
         %plots of the overall mean (for each SZA bin) dNd/dSZA vs SZA for constant latitude points

%        run sza_multiple_lats_process
              
days=[9:16];

        titlenam = ['N positive dNd dSZA vs SZA for days ' num2str(days(1)) ' to ' num2str(days(end)) ' at constant latitude for ' tit(1).tit];        
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=[0 350];
        
        izlim=1;
        zmin=0;
        zmax=90;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

%        ylab=ylabelstr;
%        xlab=xlabelstr;
%        xlab='N points';

ylab = 'Max SZA';
xlab = 'N datapoints with dNd/dSZA positive';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

%        idat=0;
%        idat=idat+1;

    ipos=1;    

    
    dat_grad=grad_sza(:,:,days);   
    iposit=find(dat_grad>0);
    inan=find(isnan(dat_grad)==1);
    temp_dat = ones(size(dat_grad));
    temp_dat(inan)=0;
    Ntot = sum(sum(temp_dat,2),3);
    
    xdat(ipos).x=zeros(size(dat_grad));
    xdat(ipos).x(iposit)=1;
    xdat(ipos).x=sum(sum(xdat(ipos).x,2),3);
    xdat(ipos).x=sum(sum(xdat(ipos).x,2),3)./Ntot;    
%    xdat(ipos).x=Ntot;
        
%    xdat(ipos).x=prod(size(iposit));
%        xdat(ipos).x=grad_sza(:,idat,2);
        
        ydat(ipos).y=mid_sza(:,1,1); %
%        mid_lats(idat)=0.5*(LAT_sza(1,1)+LAT_sza(2,1));
        %this is the mid point of the bin - NOTE is not LAT_width wide
        %e.g. if LAT_width=2 then the actual bin width is 3 degrees        
        labs(ipos).l='';
        
        ismooth_x=1;
        smooth_mode='sum';
        smooth_mode='mean';        
        Nsmooth_window = 30;
        
        
                        

        
        
        
        
        
         case 103
%plots of the overall mean (for each SZA bin) dNd/dSZA vs SZA for constant latitude points

%        run sza_multiple_lats_process

  xlab = 'Mean Normalised dN_d/dSZA';
%  xlab = 'N datapoints';
    
              
days=[339:365 1:6];
days=[339:355];
days=[356:365 1:6];

        titlenam = ['Overall mean dNd dSZA vs SZA for days ' num2str(days(1)) ' to ' num2str(days(end)) ' at constant latitude for ' tit(1).tit];        
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=[-0.3 0.3];
        xlimits=[-2 2];
        
        izlim=1;
        zmin=0;
        zmax=90;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

%        ylab=ylabelstr;
%        xlab=xlabelstr;
%        xlab='N points';

ylab = 'Max SZA';





        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

%        idat=0;
%        idat=idat+1;

    ipos=1;    
    
    %normalise the gradients by the mean value over each latitude profile
    mean_for_each_lat = meanNoNan(Nd_sza(:,:,:),1);
    mean_for_each_lat_3D = repmat(mean_for_each_lat,[1 1 size(Nd_sza,1)-1]);
    mean_for_each_lat_3D = permute(mean_for_each_lat_3D,[3 1 2]);
    
    mean_for_each_lat_3D2 = repmat(mean_for_each_lat,[1 1 size(Nd_sza,1)]);
    mean_for_each_lat_3D2 = permute(mean_for_each_lat_3D2,[3 1 2]);
    
    mean_for_each_lat = meanNoNan(Nd_rebin(:,:,:),1);
    mean_for_each_lat_3D_rebin = repmat(mean_for_each_lat,[1 1 size(Nd_rebin,1)-1]);
    mean_for_each_lat_3D_rebin = permute(mean_for_each_lat_3D_rebin,[3 1 2]);
    
    sza_thresh=40;
    
    dsza=diff(mid_sza(1:2));
    sza1D=[mid_sza-dsza/2 mid_sza(end)+dsza/2];
    sza3D=squeeze( repmat(sza1D,[1 1 size(Nd_sza,2) size(Nd_sza,3)]) );
    
    ithresh=find(sza1D<sza_thresh);
    factor1D = 10 * (sza1D-sza_thresh)/ (90-sza_thresh) ;
    factor1D = 5 * exp ( ( (sza1D-sza_thresh)/ (90-sza_thresh) ).^2 ) ;    
%    factor1D = 1 * ones(size(sza1D));
    
    factor1D(ithresh)=0;
    
    factor3D=squeeze( repmat(factor1D,[1 1 size(Nd_sza,2) size(Nd_sza,3)]) );
    Nd_sza_test = Nd_sza.*factor3D;
    Nd_sza_test = mean_for_each_lat_3D2.*factor3D;    
    
    grad_sza_test = diff(Nd_sza_test,1)./diff(sza3D,1);   
    

    
    %%%%%%
    thresh_Np=3;
    Nd_sza_thresh = Nd_sza;
%    Nd_sza_thresh = Nd_rebin;
    
    ithresh_Np = find(Np_sza<thresh_Np);
%    ithresh_Np = find(Np_rebin<thresh_Np);
    
    Nd_sza_thresh(ithresh_Np)=NaN;  %remove data with only a few datapoints for the gradient calc

    sza3D_rebin=squeeze( repmat(new_sza_bins,[1 1 size(Nd_rebin,2) size(Nd_rebin,3)]) );
    
%    grad_sza = diff(Nd_sza_thresh,1)./diff(sza3D_rebin,1);
    
        grad_sza = diff(Nd_sza_thresh,1)./diff(sza3D,1);

   
    
    

%    grad_sza_norm = grad_sza./mean_for_each_lat_3D_rebin;    
    grad_sza_norm = grad_sza./mean_for_each_lat_3D;        
    %test version with artificial factor applied
%    grad_sza_norm = grad_sza_test./mean_for_each_lat_3D;    
    


switch xlab
    case 'Mean Normalised dN_d/dSZA'
        xdat(ipos).x=meanNoNan( meanNoNan(grad_sza_norm(:,:,[days]),3) , 2);

        %    xlab = 'Mean Unnormalised dN_d/dSZA';
        %    xdat(ipos).x=meanNoNan( meanNoNan(grad_sza(:,:,[days]),3) , 2);

        ydat(ipos).y=mid_sza; %
%        ydat(ipos).y=0.5*(new_sza_bins(1:end-1) + new_sza_bins(2:end) );
        
        xlims=1;
        xlimits=[-0.3 0.3];
        xlimits=[-2 2];
        xlimits=[-0.5 0.5];

    case 'N datapoints'
        Np_sza_thresh = Np_sza;
        Np_sza_thresh(ithresh_Np)=NaN;
        [mean_nnan,ndat]=meanNoNan(Np_sza_thresh,3);
        sum_nnan = mean_nnan.*ndat;
        [mean_nnan,ndat]=meanNoNan(sum_nnan,2);
        xdat(ipos).x= mean_nnan.*ndat;
        
        
        
        ydat(ipos).y=sza1D;

end
    
    
%        xdat(ipos).x=grad_sza(:,idat,2);
        
        

        %this is the mid point of the bin - NOTE is not LAT_width wide
        %e.g. if LAT_width=2 then the actual bin width is 3 degrees        
        labs(ipos).l='';
        
        ismooth_x=1;
        Nsmooth_window = 10;
        
                        

    
 
    
        
     case 102
%plots of dNd/dSZA vs SZA for constant latitude points

%        run sza_multiple_lats_process
              
days=1;

        titlenam = ['dNd dSZA vs SZA for days ' num2str(days(1)) ' to ' num2str(days(end)) ' at constant latitude for ' tit(1).tit];        
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=[0 350];
        
        izlim=1;
        zmin=0;
        zmax=90;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
        xlab=xlabelstr;
%        xlab='N points';

ylab = 'Max SZA';
xlab = 'dNd/dSZA';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

%        idat=0;
%        idat=idat+1;

    ipos=0;
    
    idats=length(LAT_vals_sza)-40:1:length(LAT_vals_sza);
    idats=1:length(LAT_vals_sza);
    for idat=1:length(idats)
        ipos=ipos+1;

        %        xdat(idat).x=X_mean; %
        xdat(ipos).x=meanNoNan(grad_sza(:,idats(idat),days),3);
%        xdat(ipos).x=grad_sza(:,idat,2);
        
        ydat(ipos).y=mid_sza(:,idat,1); %
        mid_lats(idat)=0.5*(LAT_sza(1,idats(idat))+LAT_sza(2,idats(idat)));
        
        labs(ipos).l=num2str( mid_lats(idat) );
        %this is the mid point of the bin - NOTE is not LAT_width wide
        %e.g. if LAT_width=2 then the actual bin width is 3 degrees
        
        if length(xdat(ipos).x(isnan(xdat(ipos).x))) == length(xdat(ipos).x)
            xdat(ipos)=[];
            ydat(ipos)=[];
            labs(ipos)=[];
            ipos=ipos-1;
        end
                        
    end
    
    
    if LAT_sza(2,1)<0
        xdat2=xdat;
        ydat2=ydat;
        labs2=labs;
                 
        for idat=1:length(xdat)            
            
            xdat(length(xdat)-idat+1) = xdat2(idat);              
            ydat(length(ydat)-idat+1) = ydat2(idat);            
            labs(length(labs)-idat+1) = labs2(idat);

        end


    end
    
    
        
    case 101
%plots of Nd vs SZA for constant latitude points

%        run sza_multiple_lats_process
              

        titlenam = ['Nd vs SZA at constant latitude for ' tit(1).tit];        
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=1;
        xlimits=[0 350];
        
        izlim=1;
        zmin=0;
        zmax=80;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
        xlab=xlabelstr;
%        xlab='N points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

%        idat=0;
%        idat=idat+1;

    ipos=0;
    for idat=1:length(LAT_vals_sza)
        ipos=ipos+1;

        %        xdat(idat).x=X_mean; %
        xdat(ipos).x=Nd_sza(:,idat);
        ydat(ipos).y=SZA_sza(:,idat); %
        labs(ipos).l=num2str( 0.5*(LAT_sza(1,idat)+LAT_sza(2,idat)) );
        %this is the mid point of the bin - NOTE is not LAT_width wide
        %e.g. if LAT_width=2 then the actual bin width is 3 degrees
        
        if length(xdat(ipos).x(isnan(xdat(ipos).x))) == length(xdat(ipos).x)
            xdat(ipos)=[];
            ydat(ipos)=[];
            labs(ipos)=[];
            ipos=ipos-1;
        end
                        
    end
    
    
    if LAT_sza(2,1)<0
        xdat2=xdat;
        ydat2=ydat;
        labs2=labs;
                 
        for idat=1:length(xdat)            
            
            xdat(length(xdat)-idat+1) = xdat2(idat);              
            ydat(length(ydat)-idat+1) = ydat2(idat);            
            labs(length(labs)-idat+1) = labs2(idat);

        end


    end
    
    

        
    case 100
%mode along one dimension of a 2D histogram
              

        titlenam = ['Mode for ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
        xlab=[xlabelstr ' mode'];
%        xlab='N points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;


        idat=idat+1;
%        xdat(idat).x=X_mean; %
%        xdat(idat).x=sum(qh,2);
        xdat(idat).x=X_mode;
        
        ydat(idat).y=Ybins; %
        labs(idat).l='';        
        
    case 99
%plot of N points along one dimension of a 2D histogram
              

        titlenam = ['N points for ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
%        xlab=xlabelstr;
        xlab='N points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;


        idat=idat+1;
%        xdat(idat).x=X_mean; %
        xdat(idat).x=sum(qh(1:end-1,1:end-1),2);
        ydat(idat).y=mid_Ybins; %
        labs(idat).l='';
        

     case 98
         % 1D PDF of the data produced by pdf2D_plot_commands_L2 (just for
         % x-axis)

         dX=(Xbins(end)-Xbins(1))/nXpdf;
         Xbins(1)=Xbins(1)-dX/1e6;
         Xbins(end)=Xbins(end)+dX/1e6;

         %the y values that we want to average are just the counts - each pixel
         %counts as one point to add to the number distribution
         Y=ones(size(X(:)));

         %the x-values will be Nd, etc. that we want to bin into
         %this function ignores the NaN data
         [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals,meanvals_ignore_zeros]=bin_data(X(:),Y,Xbins);

         sizX=size(X);
         [maxval,imean]=max(sizX);
         meanX=meanNoNan(X,imean);
         
         titlenam = ['Normalised frequency distribution for ' thresh_str extra_title_info '. Mean=' num2str(meanX)];

         figname=[titlenam ' ' xlabelstr];
         savename=[savedir figname];

         xlims=0;
         xlimits=1000*[0 0.025];

         izlim=0;
         zmin=1500;
         zmax=3000;

         nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

         ylab='Normalised Frequency';
         xlab= xlabelstr;



         lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

         idat=0;


         idat=idat+1;
         xdat(idat).x=mid_points; %
         ydat(idat).y=sum_vals./sum(sum_vals(~isnan(sum_vals))); %
         labs(idat).l='';
         

    ismooth_y(idat)=0;
    Nsmooth_window=5;
    smooth_mode='mean';
         edit 

   case 977
       %1D PDF from the 2D PDF data

       %iuse_overall_norm=1;

       if ~exist('logbin_norm')
           logbin_norm=0;
       end
       
       if ~exist('axis1D')
           axis1D = 'x';
           %axis1D = 'y';
       end
       


       ikeep_xabove_zero=0; %flag to stop x being below zero (usually if want a log scale).


       if ~exist('i_plot_norm')
           i_plot_norm=1;
       end

       if ~exist('i_div_bin_widths')
           i_div_bin_widths=1;
           %    i_div_bin_widths=0;
       end

       %clear i_div_bin_widths; i_div_bin_widths=1

       %x_axis_type = 'log10_matlab'; %might also want to set this
       if ~exist('ioverride_pdf')
           
           pdf_type = 'normal';
           %       pdf_type = 'cumulative';
           %       pdf_type = 'cumulative * number'; %I.e. the total value of the quantity for each bin
       end

       cumulative=0;

        savename = ['1D-PDF ' tit(1).tit];

       figname=savename;
       switch screen_type
           case 'none'
               titlenam=[figname ' ' xlabelstr ' ' ylabelstr];
           case 'gcm_screening'
               titlenam=[figname ' ' thresh_str];
       end

       xlims=0;
       xlimits=1000*[0 0.025];

       izlim=0;
       zmin=1500;
       zmax=3000;

       nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.


       lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane





       idat=0;


%Will take the approach of calculating both the normal and cumulative and
%then selecting the desired one at the end. This is so that can outside
%scripts can choose to save either one without having to set pdf_type
       clear ydat_norm xdat_norm xdat_cum ydat_cum

       idat=idat+1;
       switch axis1D %the axis to have as the base of the 1D distribution
           case 'y'
               xlab=ylabelstr;

               ydat_norm(idat).y=(sum(qh(1:end-1,1:end-1),2))'; %transpose otherwise is a different size to bin_widths
               xdat_norm(idat).x=mid_Ybins; %
               



                       
                       

% -- Calculate cumulative frequency

                       N = sum(ydat_norm(idat).y);
                       ytemp = cumsum( ydat_norm(idat).y )/N;
                       ydat_cum(idat).y = cat(1,[0 ytemp]);
                       %bin edges
                       xdat_cum(idat).x = Ybins;
                       
                       
                       

% -- Calculate cumulative VOLUME frequency - Tomek thinks this is frequency
% times rain rate and thena cumulative plot of that

%                        N = sum(ydat_norm(idat).y);
%                        ytemp=cumsum(ydat_norm.y/N.*mid_Xbins); %cumulative sum of the frequency * rain rate (but plot against original
%                                                                 %rain rate
%                        ytemp = ytemp/ytemp(end); %Scale to 1 at the end of the cumulative sum
%                        ydat_cum_vol(idat).y = cat(1,[0 ytemp]); %add a zero to the start
%                        %bin edges
%                        xdat_cum_vol(idat).x = Ybins;
%                        

  


                       %                ydat(idat).y=(sum(qh(1:end-1,1:end-1),2))'; %transpose otherwise is a different size to bin_widths
                       %                xdat(idat).x=mid_Ybins; %
               

                       
                       if logbin_norm==1
                           bin_widths_norm = diff(log10(Ybins));
                       else
                           bin_widths_norm = diff(Ybins);
                       end
                       
               switch pdf_type
                   case 'normal'
                       ylab='N datapoints';
                       
                   case 'cumulative'
                       ylab='Cumulative Frequency'; 
                       
                       bins = Xbins; %Just in case are used elsewhere for legacy
                       bin_widths = diff(bins);
                       
                   case 'cumulative * number'
                       ylab='Cumulative Volume Frequency'; 
                       
                       bins = Xbins; %Just in case are used elsewhere for legacy
                       bin_widths = diff(bins);    

               end



           case 'x'
               xlab=xlabelstr;
               ydat_norm(idat).y=sum(qh(1:end-1,1:end-1),1);
               xdat_norm(idat).x=mid_Xbins; 
               
               
               N = sum(ydat_norm(idat).y);
               ytemp = cumsum( ydat_norm(idat).y )/N;
               ydat_cum(idat).y = cat(1,[0 ytemp]);
               %bin edges
               xdat_cum(idat).x = Xbins;






               
               if logbin_norm==1
                   bin_widths_norm = diff(log10(Xbins));
               else
                   bin_widths_norm = diff(Xbins);
               end
                       
               switch pdf_type
                   case 'normal'    
                       ylab='N datapoints';
                       
                   case 'cumulative'
                       bins = Ybins;
                       bin_widths = diff(bins); %Just in case are used elsewhere for legacy
                       
                       ylab='Cumulative Frequency';
                       
                   case 'cumulative * number'
                       bins = Ybins;
                       bin_widths = diff(bins); %Just in case are used elsewhere for legacy
                       
                       ylab='Cumulative Volume Frequency';    
               end


end
       


       %for bin normalised n_i/N = f_i * dx_i, where f_i is the y-axis and
       %N=sum(ni)
       %so f_i = n_i/N/dx_i. Need to calculate N before dividing by bin widths
       
       %Calculate normalised PDFs if required (not for cumulative plots)

       ydat_bk = ydat_norm(idat).y; 
       N = sum(ydat_bk);

       labs(idat).l='';
       
       %Dividing by the bin widths (unnormalised to give n_i from area in
       %each bar
       if i_div_bin_widths==1
           ydat_norm(idat).y = ydat_norm(idat).y ./ bin_widths_norm;
       end
          
       %Normalise if required (total area sums to one)
       if i_plot_norm==1
           ydat_norm(idat).y = ydat_norm(idat).y / N;
       end

       switch pdf_type
           case 'normal'

               if i_plot_norm==1
                   ylab='Normalised Frequency';
               end

               if i_div_bin_widths==1
                   ylab='Bin Normalised Frequency';

                   if logbin_norm==1
                       ylab='df/dlog(x)';
                   end
               end                         

       end
       
       
       
       
       switch pdf_type
           case 'cumulative'
               xdat = xdat_cum;
               ydat = ydat_cum;
           case 'cumulative * number'
               xdat = xdat_cum_vol;
               ydat = ydat_cum_vol;    

           case 'normal'
               xdat = xdat_norm;
               ydat = ydat_norm;
               bin_widths = bin_widths_norm;

       end
       

       if ikeep_xabove_zero==1
           xdat(idat).x(xdat(idat).x<0) = 1e-3;
       end



    case 97
        %plot of mean values along one dimension of a 2D histogram

        titlenam = ['% days for CF=1 for ' tit(1).tit];

        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];

        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
        xlab='% of days';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;


        idat=idat+1;
        xdat(idat).x=100*pdf_normX(:,end-1); %is end-1 because a row of zeros was added for pcolor plotting
        ydat(idat).y=Ybins; %
        labs(idat).l='';
        
    case 966
%plot of mean values along one dimension of a 2D histogram
%Update of case 96 where can plot mean of either axis


       if ~exist('axis1D')
           axis1D = 'x';
           axis1D = 'y';
       end
       
       thresh_Ndatap=100;
       thresh_Ndatap=30;
%       thresh_Ndatap=2;
       %        thresh_Ndatap=100;


       
       switch axis1D
           case 'x'                             
               titlenam = ['X mean for ' tit(1).tit];
      %        titlenam = ['N points for ' tit(1).tit];        
           case 'y'
               titlenam = ['Y mean for ' tit(1).tit];
       end
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=1;
        
        switch xlabelstr
            case '\sigma_{CTT} (K)'
                xlimits = [0 3];
            case {'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 1.6 \mum} (\mum)','R_{eff 2.1 \mum} (\mum)','R_{eff 3.7 \mum} (\mum)'}
                xlimits=[10 15]; %Reff
            case {'N_d (cm^{-3})','N_{d Re 1.6} (cm^{-3})','N_{d Re 3.7} (cm^{-3})'}
                xlimits=[60 240];
                xlimits=[40 240];  
                xlimits=[0 200];                  
            case 'Mean Optical Depth'
                xlimits=[0 50]; %Tau
            case 'Mean Optical Depth ^{1/2}'
                 xlimits=[0 100].^0.5; %Tau
            case 'Reff ^{-5/2} (\mum ^{-2.5})'
                xlimits=[20 5].^(-5/2); %Tau
            case 'R_e uncertainty (%)'
                xlimits=[0 30];                
             case 'Optical Depth Uncertainty (%)'
                xlimits=[0 50];
             case 'Optical Depth Uncertainty'
                xlimits=[0 50];  
            case 'Cloud Fraction from grid vals timeseries3'
                xlimits=[0 1];
             case 'Cloud Fraction'
                xlimits=[0 1];   
            otherwise
                xlims=0;
                
        end
        

        
        
        switch xlabelstr
            case 'Mean SZA timeseries3'

                izlim=1;
                zmin=40; zmax=85;
                %        zmin=0; zmax=65;

            case 'Cloud Fraction from grid vals timeseries3'
                izlim=1;
                zmin=0; zmax=1;
        end
        

        switch axis1D
            case 'y'   %plotting mean of each X bin on yaxis
                izlim_save=izlim;
                zmin_save =[zmin zmax];
                
                izlim = xlims;
                zmin = xlimits(1);
                zmax = xlimits(2);
                
                xlims = izlim_save;
                xlimits = zmin_save;
                
                ylab=ylabelstr;
                xlab=xlabelstr;
                
                
                idat=1;
                ydat(idat).y = Y_mean; %
                
                Ndatap = sum(qh(1:end-1,1:end-1),1); %total number of points in each X bin

                %        xdat(idat).x=sum(qh,2);
                xdat(idat).x=mid_Xbins; %
                labs(idat).l='';                
                
            case 'x'  %plotting the binining variable on the y-axis and the means on the x-axis
                ylab=ylabelstr;
                xlab=xlabelstr;
                
                idat=1;
                xdat(idat).x = X_mean; %
                
                Ndatap = sum(qh(1:end-1,1:end-1),2); %total number of points in each Y bin

                %        xdat(idat).x=sum(qh,2);
                ydat(idat).y=mid_Ybins; %
                labs(idat).l='';
        
        
                
        end
        
        %Will remove points that have fewer datapoints then thresh_Ndatap
        %So if set to 1 then it will include points with one
        %datapoint
                
        xdat(idat).x(Ndatap<thresh_Ndatap)=NaN;
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];


        

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.


%        xlab='N points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        

        




        
        if ~exist('man_choose_water_graph') | man_choose_water_graph==0
            iwrite_text_dat=1;
        end
        
        clear text_dat
        for idat=1:length(Ndatap)
            str_text=num2str(Ndatap(idat));
            text_dat(1).text(idat,1:length(str_text)) = str_text;
        end
        
        if xlims==1
            dx = ( xlimits(end) - xlimits(1) )/30;
        else
            dx = ( max(xdat(1).x) - min(xdat(1).x) ) /30;
        end

        xtext_dat(1).x = xdat(1).x + dx;
        ytext_dat(1).y = ydat(1).y;
        
        switch axis1D
            case 'x'
                ierror_bars='horiz2';
                errordatU(1).dat = std_dev_X ./ sqrt(NX_vals);
                errordatL(1).dat = std_dev_X ./ sqrt(NX_vals);

            case 'y'
                ierror_bars='vert2';
                errordatU(1).dat = std_dev_Y ./ sqrt(NY_vals);
                errordatL(1).dat = std_dev_Y ./ sqrt(NY_vals);
        end
        
        
        
    case 96
%plot of mean values along one dimension of a 2D histogram
              
        titlenam = ['X mean for ' tit(1).tit];
%        titlenam = ['N points for ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=1;
        
        switch xlabelstr
            case '\sigma_{CTT} (K)'
                xlimits = [0 3];
            case {'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 1.6 \mum} (\mum)','R_{eff 2.1 \mum} (\mum)','R_{eff 3.7 \mum} (\mum)'}
                xlimits=[10 15]; %Reff
            case {'N_d (cm^{-3})','N_{d Re 1.6} (cm^{-3})','N_{d Re 3.7} (cm^{-3})'}
                xlimits=[60 240];
                xlimits=[40 240];  
                xlimits=[0 200];                  
            case 'Mean Optical Depth'
                xlimits=[0 50]; %Tau
            case 'Mean Optical Depth ^{1/2}'
                 xlimits=[0 100].^0.5; %Tau
            case 'Reff ^{-5/2} (\mum ^{-2.5})'
                xlimits=[20 5].^(-5/2); %Tau
            case 'R_e uncertainty (%)'
                xlimits=[0 30];                
             case 'Optical Depth Uncertainty (%)'
                xlimits=[0 50];
             case 'Optical Depth Uncertainty'
                xlimits=[0 50];  
            case 'Cloud Fraction from grid vals timeseries3'
                xlimits=[0 1];
             case 'Cloud Fraction'
                xlimits=[0 1];   
            otherwise
                xlims=0;
                
        end
        
        switch xlabelstr
            case 'Mean SZA timeseries3'

                izlim=1;
                zmin=40; zmax=85;
                %        zmin=0; zmax=65;

            case 'Cloud Fraction from grid vals timeseries3'
                izlim=1;
                zmin=0; zmax=1;
        end
        

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
        xlab=xlabelstr;
%        xlab='N points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        Ndatap = sum(qh(1:end-1,1:end-1),2); %total number of points in each Y bin
        



        idat=idat+1;
        xdat(idat).x=X_mean; %
        
%Will remove points that have fewer datapoints then the threshold below
%So if set to 1 then it will include points with one datapoint
        thresh_Ndatap=100;
        thresh_Ndatap=30;
        thresh_Ndatap=2;        
%        thresh_Ndatap=100;    
        
        xdat(idat).x(Ndatap<thresh_Ndatap)=NaN;
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];
        
%        xdat(idat).x=sum(qh,2);
        ydat(idat).y=mid_Ybins; %
        labs(idat).l='';
        
        if ~exist('man_choose_water_graph') | man_choose_water_graph==0
            iwrite_text_dat=1;
        end
        
        clear text_dat
        for idat=1:length(Ndatap)
            str_text=num2str(Ndatap(idat));
            text_dat(1).text(idat,1:length(str_text)) = str_text;
        end
        
        if xlims==1
            dx = ( xlimits(end) - xlimits(1) )/30;
        else
            dx = ( max(xdat(1).x) - min(xdat(1).x) ) /30;
        end

        xtext_dat(1).x = xdat(1).x + dx;
        ytext_dat(1).y = ydat(1).y;
        
        ierror_bars='horiz2';
        errordatU(1).dat = std_dev_X ./ sqrt(NX_vals);
        errordatL(1).dat = std_dev_X ./ sqrt(NX_vals);
        
        
    case 961
        %This is a modfied copy of case 96.
%construct a cloud frequency array where if the CF was zero then we have an
%array of zeros. If there was cloud then the fraction occurs at one height
%with all CFs above zero and all CFs below NaN to simulate CALIPSO
              
        titlenam = ['X mean for ' tit(1).tit];
%        titlenam = ['N points for ' tit(1).tit];        
        
        figname=titlenam;
        savename=[savedir figname];

        xlims=1;
        
        switch xlabelstr
            case '\sigma_{CTT} (K)'
                xlimits = [0 3];
            case {'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7','R_{eff 1.6 \mum} (\mum)','R_{eff 2.1 \mum} (\mum)','R_{eff 3.7 \mum} (\mum)'}
                xlimits=[10 15]; %Reff
            case {'N_d (cm^{-3})','N_{d Re 1.6} (cm^{-3})','N_{d Re 3.7} (cm^{-3})'}
                xlimits=[60 240];
                xlimits=[40 240];  
                xlimits=[0 200];                  
            case 'Mean Optical Depth'
                xlimits=[0 50]; %Tau
            case 'Mean Optical Depth ^{1/2}'
                 xlimits=[0 100].^0.5; %Tau
            case 'Reff ^{-5/2} (\mum ^{-2.5})'
                xlimits=[20 5].^(-5/2); %Tau
            case 'R_e uncertainty (%)'
                xlimits=[0 30];                
             case 'Optical Depth Uncertainty (%)'
                xlimits=[0 50];
             case 'Optical Depth Uncertainty'
                xlimits=[0 50];  
            case 'Cloud Fraction from grid vals timeseries3'
                xlimits=[0 1];
             case 'Cloud Fraction'
                xlimits=[0 1];   
            otherwise
                xlims=0;
                
        end
        
        switch xlabelstr
            case 'Mean SZA timeseries3'

                izlim=1;
                zmin=40; zmax=85;
                %        zmin=0; zmax=65;

            case 'Cloud Fraction from grid vals timeseries3'
                izlim=1;
                zmin=0; zmax=1;
        end
        

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab=ylabelstr;
        xlab=xlabelstr;
%        xlab='N points';



        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        Ndatap = sum(qh(1:end-1,1:end-1),2); %total number of points in each Y bin
        

%Construct CF array. Will have length(Ybins) (Y=height) for each datapoint
%within the histogram. So will be of size [length(Ybins) * Ntot]

        Ntot = sum(Ndatap,1);
%        CF_array = zeros([Ntot length(Ybins)]); %actually, this causes memory issues - jstu take the sum as go along and keep track of N
        N_vs_CTH = zeros([1 length(Ybins)-1]);
        CFsum_vs_CTH = zeros([1 length(Ybins)-1]);        
        
%        iNtot=0;
        for iY=1:size(qh,1)-1
            for iX=1:size(qh,2)-1
%                iNtot=iNtot+1;
                %populate the height of the bin with the cloud fraction in
                %Xbins - the rest of the column will be zeros since this is
                %what was initialized with.
%                CFarray(iNtot,iY)=Xbins(iX);

                Nbin = qh(iY,iX);
                CFsum_vs_CTH(iY)=CFsum_vs_CTH(iY) + Nbin.*Xbins(iX);
                %Only add one to number if we are above the height in
                %question (akin to setting all below to NaN)
                N_vs_CTH(iY:end) = N_vs_CTH(iY:end) + Nbin; 
%                N_vs_CTH(iY) = N_vs_CTH(iY) + Nbin;   %in case just want to count the actual height (akin to all other heights =NaN)              
                
                %Set all altiudes below this to NaN
%                CFarray(iNtot,1:iY-1)=NaN;                
                %For the points that had zero CF and so no CTH will have
                %set those to low CTH that will be out of range of the
                %Ybins and so shoudl remain as all zeros
                
            end
            
        end
        

        idat=idat+1;
%        xdat(idat).x=meanNoNaN(CFarray,1); %Now average all of those CFs
        N_vs_CTH(N_vs_CTH==0)=NaN; %avoid divide by zero
        xdat(idat).x = CFsum_vs_CTH ./ N_vs_CTH; %Now average all of those CFs  
      
        
%Will remove points that have fewer datapoints then the threshold below
%So if set to 1 then it will include points with one datapoint
        thresh_Ndatap=100;
        thresh_Ndatap=30;
        thresh_Ndatap=2;        
%        thresh_Ndatap=100;    
        
        xdat(idat).x(Ndatap<thresh_Ndatap)=NaN;
        titlenam=[titlenam ' with N.LT.' num2str(thresh_Ndatap) ' removed'];
        
        xdat(idat).x = xdat(idat).x(2:end)'; %Ignore the first bin since this is what the CTH was set to when have low CF
        
%        xdat(idat).x=sum(qh,2);
        ydat(idat).y=mid_Ybins(2:end); %
        labs(idat).l='';
        
        if ~exist('man_choose_water_graph') | man_choose_water_graph==0
            iwrite_text_dat=1;
        end
        
        clear text_dat
        for idat=1:length(Ndatap)
            str_text=num2str(Ndatap(idat));
            text_dat(1).text(idat,1:length(str_text)) = str_text;
        end
        
        if xlims==1
            dx = ( xlimits(end) - xlimits(1) )/30;
        else
            dx = ( max(xdat(1).x) - min(xdat(1).x) ) /30;
        end

        xtext_dat(1).x = xdat(1).x + dx;
        ytext_dat(1).y = ydat(1).y;
        
        ierror_bars='horiz2';
        errordatU(1).dat = std_dev_X ./ sqrt(NX_vals);
        errordatL(1).dat = std_dev_X ./ sqrt(NX_vals);
        
        
        xdat_swap = xdat;
        ydat_swap = ydat;
        
        for idat=1:length(xdat)
        xdat(idat).x = ydat_swap(idat).y;
        ydat(idat).y = xdat_swap(idat).x;
        end
        
        
 
case 952
        % condensation rate (dLWC/dz)
        
      fsize=15;
        titlenam = ['Condensation rate'];
        
        figname=titlenam;
        savename=figname;


        
        izlim=0;
        zmin=0;
        zmax=8e-3;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='dLWC/dz (g m^{-4})';
        xlab= 'Temperature (K)';


        T=[270:0.2:313];
        T=[230:0.2:313];  
        
        xlims=1;
        xlimits=[270 300];
        xlimits=[T(1) T(end)];
        
        P_plot=[1000 700 500];

        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        idat=0;
        
       for iP=1:length(P_plot)
        P=P_plot(iP);
        Parr=ones(size(T))*P*100; %convert to Pa
        [dqldz,dqldz_dan,dLWCdT,dqldz_dan2,dqldz_dan3,dqldz_dan4,dqdz_num,dqldz_dan5,dqldz_dan6,dqldz_dan7,dqldz_dan_goff,dqldz_Lc,qs_flag,dqldz_sami]=adlwcgm2_rob(T,Parr);
% 
         idat=idat+1;
         xdat(idat).x=T; %
         ydat(idat).y=dqldz_Lc; %
         labs(idat).l=['Rob constant L, P=' num2str(P) ' hPa'];
        
        idat=idat+1;
        xdat(idat).x=T; %
        ydat(idat).y=dqldz_dan7; %
        labs(idat).l=['Albrecht dq_dT, P=' num2str(P) ' hPa'];
        
        idat=idat+1;
        xdat(idat).x=T; %
        ydat(idat).y=dqldz_dan5; %
        labs(idat).l=['Albrecht, P=' num2str(P) ' hPa'];
%         
%         idat=idat+1;
%         xdat(idat).x=T; %
%         ydat(idat).y=dqldz_dan6; %
%         labs(idat).l=['Albrecht__num, P=' num2str(P) ' hPa'];
        
%          idat=idat+1;
%         xdat(idat).x=T; %
%         ydat(idat).y=dqldz_dan_goff; %
%         labs(idat).l=['Albrecht dq_dT ' qs_flag ', P=' num2str(P) ' hPa'];

        idat=idat+1;
        xdat(idat).x=T; %
        ydat(idat).y=dqldz_sami; %
        labs(idat).l=['Sami, P=' num2str(P) ' hPa'];
        
       end
        

        
        case 953
        % condensation rate (dLWC/dz)
        
      
        titlenam = ['Latent Heat Temperature Variation'];
        
        figname=titlenam;
        savename=figname;


        
        izlim=0;
        zmin=0;
        zmax=8e-3;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='% change from T=273 K';
        xlab= 'Temperature (K)';


        T=[270:0.2:300];
        T=[230:0.2:300];  
        
        xlims=1;
        xlimits=[270 300];
        xlimits=[T(1) T(end)];
        
        P_plot=[1000];

        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        idat=0;
        
        LT=latent_heat_temp_variation(T);
        [minval,i0]=min(abs(T-273.15));
        
        
       for iP=1:length(P_plot)
%         P=P_plot(iP);
%         Parr=ones(size(T))*P*100; %convert to Pa

         idat=idat+1;
         xdat(idat).x=T; %
         ydat(idat).y=(LT-LT(i0))./LT(i0)*100; %
         labs(idat).l=['L(T)'];
        

       end
   
 
        
 case 95
        % condensation rate (dLWC/dz)
        

        titlenam = ['Condensation rate'];
        
        figname=titlenam;
        savename=figname;

        xlims=1;
        xlimits=[270 300];
        xlimits=[200 300];
        
        izlim=1;
        zmin=0;
        zmax=8e-3;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='dLWC/dz (g m^{-4})';
        xlab= 'Temperature (K)';


        T=[270:0.2:300];
        T=[230:0.2:300];        

        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        idat=0;
        
        idat=idat+1;
        P=1010;
        Parr=ones(size(T))*P*100; %convert to Pa
        [dqldz,dqldz2,dLWCdT]=adlwcgm2_rob(T,Parr);

        xdat(idat).x=T; %
        ydat(idat).y=0.8*dqldz; %
        labs(idat).l=['P=' num2str(P) ' hPa, 80% of ad'];
        
        %%%
        
        idat=idat+1;
        P=860;
        Parr=ones(size(T))*P*100; %convert to Pa;
        [dqldz,dqldz2,dLWCdT]=adlwcgm2_rob(T,Parr);

        xdat(idat).x=T; %
        ydat(idat).y=0.8*dqldz; %
        labs(idat).l=['P=' num2str(P) ' hPa, 80% of ad'];
        
        %%%
        
        idat=idat+1;
        P=710;
        Parr=ones(size(T))*P*100; %convert to Pa;
        [dqldz,dqldz2,dLWCdT]=adlwcgm2_rob(T,Parr);

        xdat(idat).x=T; %
        ydat(idat).y=0.8*dqldz; %
        labs(idat).l=['P=' num2str(P) ' hPa, 80% of ad'];
        
         %%%
        
        idat=idat+1;
        P=710;
        Parr=ones(size(T))*P*100; %convert to Pa;
        rho =density(Parr,T);        
        cw=1000.*exp(-20.41+0.02774*T); %from Ralph Bennartz


        xdat(idat).x=T; %
        ydat(idat).y=0.8*cw; %
        labs(idat).l=['P=' num2str(P) ' hPa, Ralph, 80% of ad'];
        
           %%%
        
        idat=idat+1;
        P=710;
        Parr=ones(size(T))*P*100; %convert to Pa;
        rho =density(Parr,T);        
        cw=1000.*rho.*exp(-20.41+0.02774*T); %from Ralph Bennartz


        xdat(idat).x=T; %
        ydat(idat).y=0.8*cw; %
        labs(idat).l=['P=' num2str(P) ' hPa, Ralph*rho, 80% of ad'];
        
        
case 94  %run case 'CAS hotwire matches' in scatter_plot first
        % CAS number analysis
        
        iaxis_square=0; %switch to make axis square
    
        
                
        xlims=1;
        xlimits=[-2 35];
        
        izlim=0;
        zmin=0;
        zmax=4;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        xlab= 'Mode CAS diameter (\mum)';
        
        
     

        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;                
        
        ithresh=1;        
        
        
        type_plot_num='median';
%        type_plot_num='number';
        
        switch type_plot_num
            case 'median'
                titlenam = ['CAS numbers vs mode, concentration threshold ' num2str(ithresh)];
                ylab='Meidan total droplet number (cm^{-3})';
            case 'number'
                titlenam = ['CAS number datapoints vs mode, concentration threshold ' num2str(ithresh)];
                ylab='Number datapoints for meidan total droplet number';
        end
        
        
        
        
        for idat=1:length(lwc_ratio_dat)
             xdat(idat).x=0.5*(CAS_bins(1:end-1) + CAS_bins(2:end)); %
             labs(idat).l=['Flight ' num2str(lwc_ratio_dat(idat).flight)];        
             switch ithresh
                 case 1
                      switch type_plot_num
                         case 'median'
                             ydat(idat).y=lwc_ratio_dat(idat).median_ratio;
                         case 'number'
                             ydat(idat).y=lwc_ratio_dat(idat).nvals;
                     end
                     figname=['Conc. threshold (1) ' titlenam]; 
                 case 2
                     switch type_plot_num
                         case 'median'
                             ydat(idat).y=lwc_ratio_dat(idat).median_ratio2;
                         case 'number'
                             ydat(idat).y=lwc_ratio_dat(idat).nvals2;
                     end
                     figname=['Conc. threshold (2) ' titlenam];                    
                 case 3
                     switch type_plot_num
                         case 'median'
                             ydat(idat).y=lwc_ratio_dat(idat).median_ratio3;
                         case 'number'
                             ydat(idat).y=lwc_ratio_dat(idat).nvals3;
                     end
                     figname=['Conc. threshold (3) ' titlenam];  
                 case 4
                      switch type_plot_num
                         case 'median'
                             ydat(idat).y=lwc_ratio_dat(idat).median_ratio4;
                         case 'number'
                             ydat(idat).y=lwc_ratio_dat(idat).nvals4;
                     end
                     figname=['Conc. threshold (4) ' titlenam]; 
                 case 5
                      switch type_plot_num
                         case 'median'
                             ydat(idat).y=lwc_ratio_dat(idat).median_ratio5;
                         case 'number'
                             ydat(idat).y=lwc_ratio_dat(idat).nvals5;
                      end
                     figname=['Conc. threshold (5) ' titlenam];   
                 case 6
                     switch type_plot_num
                         case 'median'
                             ydat(idat).y=lwc_ratio_dat(idat).median_ratio6;
                         case 'number'
                             ydat(idat).y=lwc_ratio_dat(idat).nvals6;
                     end
                     figname=['Conc. threshold (6) ' titlenam];   

                     
             end

             
        end
        
        
        
           
        
        savename=figname;
        
        
        
 case 93  %run case 'CAS hotwire matches' in scatter_plot first
        % CAS number analysis - overall means
            
        xlims=0;
        xlimits=1000*[0 0.025];
               
        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        
        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio;
            ratio_array(ilwc,:) = ydata;
            Nratio_array(ilwc,:) = lwc_ratio_dat(ilwc).nvals;
        end

        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio2;
            ratio_array2(ilwc,:) = ydata;
            Nratio_array2(ilwc,:) = lwc_ratio_dat(ilwc).nvals2;            
        end

        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio3;
            ratio_array3(ilwc,:) = ydata;
            Nratio_array3(ilwc,:) = lwc_ratio_dat(ilwc).nvals3;            
        end
        
        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio4;
            ratio_array4(ilwc,:) = ydata;
            Nratio_array4(ilwc,:) = lwc_ratio_dat(ilwc).nvals4;            
        end
        
        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio5;
            ratio_array5(ilwc,:) = ydata;
            Nratio_array5(ilwc,:) = lwc_ratio_dat(ilwc).nvals5;            
        end
        
        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio6;
            ratio_array6(ilwc,:) = ydata;
            Nratio_array6(ilwc,:) = lwc_ratio_dat(ilwc).nvals6;            
        end

         

        xlab= 'Mode CAS diameter (\mum)';
        
        ylab='Mean of meidan ratio CAS to hotwire';
        ylab='Mean no. datapoints for meidan ratio CAS to hotwire';
%        ylab='Mean of meidan ratio CAS to hotwire - campaign in 2 halves';

        titlenam = [ylab ' for all flights'];
        
        figname=titlenam;
        savename=figname;
       
        
        switch ylab
            case 'Mean of meidan ratio CAS to hotwire'
                izlim=0;
                zmin=0;
                zmax=4;
                lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array,1);
                labs(idat).l=['ALL points'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array2,1);
                labs(idat).l=['>25% of max'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array3,1);
                labs(idat).l=['>50% of max'];
                
                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array4,1);
                labs(idat).l=['>65% of max'];
                
                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array5,1);
                labs(idat).l=['>75% of max'];
                
                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array6,1);
                labs(idat).l=['>85% of max'];

            case 'Mean no. datapoints for meidan ratio CAS to hotwire'
                izlim=0;
                zmin=0;
                zmax=4;
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                

                idat=0;

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array,1);
                labs(idat).l=['ALL points'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array2,1);
                labs(idat).l=['>25% of max'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array3,1);
                labs(idat).l=['>50% of max'];
                
                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array4,1);
                labs(idat).l=['>65% of max'];
                
                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array5,1);
                labs(idat).l=['>75% of max'];
                
                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array6,1);
                labs(idat).l=['>85% of max'];
                
            case 'Mean of meidan ratio CAS to hotwire - campaign in 2 halves'
                izlim=1;
                zmin=0;
                zmax=4;
                lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;
                
                
                ithresh=3;      
                
                titlenam = ['CAS to hotwire ratio median, lwc threshold ' num2str(ithresh)];
                
                switch ithresh
                    case 1
                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array(1:7,:),1);
                        labs(idat).l=['108 and before'];

                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array(8:end,:),1);
                        labs(idat).l=['Post 108'];
                        
                        figname=['LWC threshold (1) ' titlenam];  
                        
                    case 2
                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array2(1:7,:),1);
                        labs(idat).l=['108 and before'];

                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array2(8:end,:),1);
                        labs(idat).l=['Post 108'];
                        
                        figname=['LWC threshold (2) ' titlenam];  
                        
                    case 3                       
                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array3(1:7,:),1);
                        labs(idat).l=['108 and before'];

                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array3(8:end,:),1);
                        labs(idat).l=['Post 108'];
                        
                        figname=['LWC threshold (3) ' titlenam];  
                
                end

        
        end
        
        
        
   case 92  %run case 'CAS hotwire matches' in scatter_plot first
        % CAS/hotwire analysis - overall means
        
    
    

        xlims=0;
        xlimits=1000*[0 0.025];
        
       

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        
        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio;
            ratio_array(ilwc,:) = ydata;
            Nratio_array(ilwc,:) = lwc_ratio_dat(ilwc).nvals;
        end

        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio2;
            ratio_array2(ilwc,:) = ydata;
            Nratio_array2(ilwc,:) = lwc_ratio_dat(ilwc).nvals2;            
        end

        for ilwc=1:length(lwc_ratio_dat)
            ydata=lwc_ratio_dat(ilwc).median_ratio3;
            ratio_array3(ilwc,:) = ydata;
            Nratio_array3(ilwc,:) = lwc_ratio_dat(ilwc).nvals3;            
        end

         

        xlab= 'Mode CAS diameter (\mum)';
        
        ylab='Mean of meidan ratio CAS to hotwire';
        ylab='Mean no. datapoints for meidan ratio CAS to hotwire';
        ylab='Mean of meidan ratio CAS to hotwire - campaign in 2 halves';

        titlenam = [ylab ' for all flights'];
        
        figname=titlenam;
        savename=figname;
       
        
        switch ylab
            case 'Mean of meidan ratio CAS to hotwire'
                izlim=1;
                zmin=0;
                zmax=4;
                lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array,1);
                labs(idat).l=['CAS 0.02, HW 0.05'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array2,1);
                labs(idat).l=['CAS 0.02, HW 0.1'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(ratio_array3,1);
                labs(idat).l=['CAS 0.1, HW 0.1'];

            case 'Mean no. datapoints for meidan ratio CAS to hotwire'
                izlim=0;
                zmin=0;
                zmax=4;
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                

                idat=0;

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array,1);
                labs(idat).l=['CAS 0.02, HW 0.05'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array2,1);
                labs(idat).l=['CAS 0.02, HW 0.1'];

                idat=idat+1;
                xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                ydat(idat).y=meanNoNan2(Nratio_array3,1);
                labs(idat).l=['CAS 0.1, HW 0.1'];
                
            case 'Mean of meidan ratio CAS to hotwire - campaign in 2 halves'
                izlim=1;
                zmin=0;
                zmax=4;
                lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;
                
                
                ithresh=3;      
                
                titlenam = ['CAS to hotwire ratio median, lwc threshold ' num2str(ithresh)];
                
                switch ithresh
                    case 1
                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array(1:7,:),1);
                        labs(idat).l=['108 and before'];

                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array(8:end,:),1);
                        labs(idat).l=['Post 108'];
                        
                        figname=['LWC threshold (1) ' titlenam];  
                        
                    case 2
                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array2(1:7,:),1);
                        labs(idat).l=['108 and before'];

                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array2(8:end,:),1);
                        labs(idat).l=['Post 108'];
                        
                        figname=['LWC threshold (2) ' titlenam];  
                        
                    case 3                       
                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array3(1:7,:),1);
                        labs(idat).l=['108 and before'];

                        idat=idat+1;
                        xdat(idat).x=0.5*(CAS_bins(1:end-1)+CAS_bins(2:end)); %
                        ydat(idat).y=meanNoNan2(ratio_array3(8:end,:),1);
                        labs(idat).l=['Post 108'];
                        
                        figname=['LWC threshold (3) ' titlenam];  
                
                end

        
        end
             
        
     
        
        
        
        
    case 91  %run case 'CAS hotwire matches' in scatter_plot first
        % CAS/hotwire analysis
        
        iaxis_square=0; %switch to make axis square
    
        
                
        xlims=1;
        xlimits=[-2 35];
        
        izlim=1;
        zmin=0;
        zmax=4;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        xlab= 'Mode CAS diameter (\mum)';
        
        
        ylab='Meidan ratio CAS to hotwire';

        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;                
        
        ithresh=2;        
        titlenam = ['CAS to hotwire ratio median, lwc threshold ' num2str(ithresh)];
        
        for idat=1:length(lwc_ratio_dat)
             xdat(idat).x=0.5*(CAS_bins(1:end-1) + CAS_bins(2:end)); %
             labs(idat).l=['Flight ' num2str(lwc_ratio_dat(idat).flight)];        
             switch ithresh
                 case 1
                     ydat(idat).y=lwc_ratio_dat(idat).median_ratio;
                     figname=['LWC threshold (1) ' titlenam];
                 case 2
                     ydat(idat).y=lwc_ratio_dat(idat).median_ratio2;
                     figname=['LWC threshold (2) ' titlenam];                     
                 case 3
                     ydat(idat).y=lwc_ratio_dat(idat).median_ratio3;
                     figname=['LWC threshold (3) ' titlenam];                     
             end

             
        end
        
        savename=figname;
        
        
        
        
    case 90
        % ACPIM profiles
        
%        tstr=Times(time,:);
%        iund=findstr('_',tstr);
%        tstr(iund)=' ';          
        
        xlims=0;
        
        izlim=1;
        zmin=1500;
        zmax=5000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        ylab='Height (m)';        
        
        xlab= 'Cloud water Mixing Ratio (g m^{-3})';
%        xlab= 'Droplet Concentration (cm^{-3})';
%        xlab= 'Potential temperature (K)';        
%        xlab= 'Water vapour (g kg^{-1})';   
%        xlab= 'RH';   
%        xlab='Pressure (hPa)';


        DT=0.5; %ACPIM timestep - should make Acpim output this
        time_plot=79; %seconds
        time_plot=100.5; %seconds
%        time_plot=0; %seconds
        
        %Matlab model index
        it2=findheight_nearest(time,time_plot);
        
        %ACPIM index
        it=time_plot/DT +1; %assuming that it=1 is time=0
        if size(acpim.qc,2)==1 %if a parcel run
            it = [1:size(acpim.qc,1)];
            zInit = 2455; %initial height of parcel (m) if a parcel run
        %if doing a parcel run then change the Matlab wave model
        %indices to be all the time indices for the right starting height
            it2 = [1:size(LWC_humi_ALL,2)];            
            iz_wave = findheight_nearest(YY(:,1),zInit);
        else
            iz_wave=[1:size(YY,1)];
            zInit=0;
        end
        
        
        
        %ACPIM density array
        rho=density(acpim.P,acpim.T);
        
        
        
        idat=0;
        
        switch xlab
            case 'Cloud water Mixing Ratio (g m^{-3})'
                idat=idat+1;                
                xdat(idat).x=1000*acpim.qc(it,:).*rho(it,:); %convert to g/m3 from kg/kg
%                xdat(idat).x=1000*acpim.qc(it,:); %convert to g/m3 from kg/kg                
                ydat(idat).y=acpim.Z(it,:)+zInit; %
                labs(idat).l='ACPIM';                               
                
                idat=idat+1;
%                xdat(idat).x=LWChumi_reg(:,it2); %
%                xdat(idat).x=max(LWChumi_reg(:,:),[],2); %  
%                ydat(idat).y=z;
                
                xdat(idat).x=LWC_humi_ALL(iz_wave,it2); %                
                ydat(idat).y=YY(iz_wave,it2); %
                labs(idat).l='Wave model humicap';
                
%                 idat=idat+1;                
%                 xdat(idat).x=LWC_fp_ALL(iz_wave,it2); %                
%                 ydat(idat).y=YY(iz_wave,it2); %
%                 labs(idat).l='Wave model humicap';
                

                %aircraft LWC - max LWC with height profiles
                idat=idat+1;
                %aircraft CAS lwc timeseries
                set_column_numbers_for_flight_data; %gives col_alt, etc.
                indsCAS=1:length(CAS_time_all);
                lwc_air = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all(indsCAS));
                z_air = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))';
                %bin the lwc into height bins
                zbins = [0:25:5000];
                [meanvals,maxvals,max_inds,mid_bins]=bin_data(z_air,lwc_air,zbins);
                maxvals(maxvals>3)=0; %remove the dodgy point
                
                xdat(idat).x = maxvals;
                ydat(idat).y = mid_bins;
                labs(idat).l='Hotwire LWC max';
                
                xlims=1;
                xlimits=[-0.2 0.8];
                
            case 'Droplet Concentration (cm^{-3})'
                idat=idat+1;
                xdat(idat).x=1e-6*acpim.conc2(it,:).*rho(it,:); %convert to #/cm3 from #/kg
                ydat(idat).y=acpim.Z(it,:)+zInit; %
                labs(idat).l='ACPIM';                               
                               
%                 %aircraft LWC - max LWC with height profiles
%                 idat=idat+1;
%                 %aircraft CAS lwc timeseries
%                 set_column_numbers_for_flight_data; %gives col_alt, etc.
%                 indsCAS=1:length(CAS_time_all);
%                 lwc_air = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all(indsCAS));
%                 z_air = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))';
%                 %bin the lwc into height bins
%                 zbins = [0:25:5000];
%                 [meanvals,maxvals,max_inds,mid_bins]=bin_data(z_air,lwc_air,zbins);
%                 maxvals(maxvals>3)=0; %remove the dodgy point
                
%                 xdat(idat).x = maxvals;
%                 ydat(idat).y = mid_bins;
%                 labs(idat).l='Hotwire LWC max';

                    xlims=0;
                    xlimits=[0 0.8];
                    
            case 'Potential temperature (K)'
                idat=idat+1;
                xdat(idat).x=acpim.T(it,:).*(1e5./acpim.P(it,:)).^0.286; %
                ydat(idat).y=acpim.Z(it,:)+zInit; %
                labs(idat).l='ACPIM';    
                
                idat=idat+1;
%                xdat(idat).x=pot_reg(:,it2); %
                 xdat(idat).x=pot_wave(iz_wave,it2); %
%                ydat(idat).y=z; %
                ydat(idat).y=YY(iz_wave,it2); %
                labs(idat).l='Wave model';
                
            case 'Water vapour (g kg^{-1})'
                idat=idat+1;
                qsat_acpim = SatVapPress(acpim.T(it,:),'buck2','liq',acpim.P(it,:),1)/f;
                xdat(idat).x=acpim.R(it,:).*qsat_acpim; %convert to #/cm3 from #/kg
                ydat(idat).y=acpim.Z(it,:)+zInit; %
                labs(idat).l='ACPIM';    
                
                idat=idat+1;
%                xdat(idat).x=qv_humi_reg(:,it2); %
%                ydat(idat).y=z; %

                xdat(idat).x=qv_wave_humi(iz_wave,it2);
                ydat(idat).y=YY(iz_wave,it2); %
                
                labs(idat).l='Wave model';
                
             case 'RH'
                idat=idat+1;                
                xdat(idat).x=acpim.R(it,:); %
                ydat(idat).y=acpim.Z(it,:)+zInit; %
                labs(idat).l='ACPIM';    
                
                idat=idat+1;
%                xdat(idat).x=qv_humi_reg(:,it2); %
%                ydat(idat).y=z; %
                qsat_wave = SatVapPress(Tad_humi_ALL(iz_wave,it2)+273.15,'buck2','liq',PY(iz_wave,it2),1)/f;
%                qsat_wave = SatVapPress(Tad_humi_ALL(iz_wave,it2)+273.15,'goff','liq',PY(iz_wave,it2),1)/f;                
                xdat(idat).x=qv_wave_humi(iz_wave,it2) ./ qsat_wave;
                ydat(idat).y=YY(iz_wave,it2); %
                
                labs(idat).l='Wave model';    
                
            case 'Pressure (hPa)'
                idat=idat+1;
                xdat(idat).x=100*acpim.P(it,:); %convert to #/cm3 from #/kg
                ydat(idat).y=acpim.Z(it,:)+zInit; %
                labs(idat).l='ACPIM';    
                
                idat=idat+1;
%                xdat(idat).x=qv_humi_reg(:,it2); %
%                ydat(idat).y=z; %

                xdat(idat).x=100*PY(iz_wave,it2);
                ydat(idat).y=YY(iz_wave,it2); %
                
                labs(idat).l='Wave model';
                
        end
   
        titlenam = ['ACPIM plots at t=' num2str(time_plot) ' s'];
        figname=[titlenam '-' xlab];
        savename=figname;
        
case 89
        % profiles from aircraft        
        
       
        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.
        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
        ylab='Height (m)';
%        ylab='Potential Temperature (K)';
%        ylab='Equivalent Potential Temperature (K)';
        
        eval_str=['dat_flt=dat_flt' flight_no ';'] %run the CAS read first or manually set flight_no 
        eval(eval_str);
        
        if ~exist('man_choose_case89_xlab')

            xlab= 'Temperature (^{o}C)';
%            xlab= 'Potential temperature (K)';       
            xlab='Vapour mixing ratio (kg kg^{-1}) from frost point hygrometer';
            xlab='Vapour mixing ratio (kg kg^{-1}) from humicap';       
%             xlab= 'Pressure (hPa)';
             xlab= 'Wind speed (m s^{-1})';
            xlab= 'Wind direction (degrees)';
%            xlab= 'Total ice number (L ^{-1})';    
%            xlab='Mean ice size (microns)';
%            xlab='Total CIP mass Jonny (mg m^{-3})';
%            xlab='Total CIP number Jonny (cm^{-3})';
%            xlab='Ice+small CIP number Jonny (L^{-1})';
            
            %set the times to plot from the timeseries for the profile
            %set to [] for all times.
            times=[14.65 14.835; 15.285 16.097]/24;
            %        times=[14.7075 14.835];
            %        times=[14.81 14.835];
            times=[14.05 14.67]/24;
            times=[13.25 14]/24;
            times=[13+35/60 14+10/60]/24;
            times=[14+0/60 14+20/60]/24;
            times=[21.6772 22.0281]/24; %flight 102 end profile
            times=[19.6861 20.3949]/24; %flight 102 start profile            
            times=[19.6861 20.3949;21.6772 22.0281]/24; %flight 102 start profile        
            
            
            times=[];



            ihighlight_cloud=0;
        else
            clear man_choose_case89_xlab
        end

        
       
            
            
            tstr=date_str;
            iund=findstr('_',tstr);
            tstr(iund)=' ';
            titlenam = ['Profile of ' xlab ' at ' tstr ' for flight ' flight_no];

            figname=titlenam;
            savename=figname;

        
clear indsALL
        if length(times)==0
            indsALL(1).i=1:length(dat_flt(:,1));
            times=0; %so that size(times,1)=1 (as used later)
        else
            for i_inds=1:size(times,1)
                [time_0,time_1]=findheight_nearest(dat_flt(:,1)/1000/3600/24,times(i_inds,1),times(i_inds,2));
                indsALL(i_inds).i=time_0:time_1;
            end
        end
        
        if size(dat_flt,2)==15
            %for flt_19
            col_temp=6;
            col_alt=11;
            col_lat=2;
            col_lon=3;
            col_press=4;
            col_wind=9;
            col_winddir=10;
        else
            %for Feb2010 flights
            col_temp=5;
            col_alt=12;
            col_lat=2;
            col_lon=3;
            col_press=6;
            col_wind=9;
            col_winddir=10;
            col_frostpoint_hygro=7;
            col_frostpoint_humi=8;
            col_airspeed=4;
        end

                idat=0;
                
        for idat2=1:size(times,1)
            
            switch ylab
                case 'Potential Temperature (K)'
                    idat=idat+1;
                    inds=indsALL(idat).i;
                    T=dat_flt(inds,col_temp)+273.15; %
                    P=dat_flt(inds,col_press); %
                    potemp = T.*(1000./P).^0.286;
                    ydat(idat).y=potemp;   
                    
                    idat=idat-1;  %reset for later use  
                    
                case 'Equivalent Potential Temperature (K)'
                    idat=idat+1;
                    inds=indsALL(idat).i;
                    T=dat_flt(inds,col_temp)+273.15; %
                    P=100*dat_flt(inds,col_press); %
                    potemp = T.*(1000./P).^0.286;
                    
                    eval_str = ['qv=qv_flt' flight_no '_fp(inds);'];
                    eval_str = ['qv=qv_flt' flight_no '_humi(inds);'];                    
                    eval(eval_str);                    
                    
                    equiv = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                    ydat(idat).y=equiv;  
                    
                    idat=idat-1;   
                otherwise %use height
                    idat=idat+1;
                    ydat(idat).y=dat_flt(inds,col_alt);
                    idat=idat-1;
                    
            end

        switch xlab
            case 'Temperature (^{o}C)'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=dat_flt(inds,col_temp); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(2),1)/1000/3600) ' UTC'];
                
            case 'Potential temperature (K)'
%                for idat=1:size(times,1)
                    idat=idat+1;
                    inds=indsALL(idat).i;
                    T=dat_flt(inds,col_temp)+273.15; %
                    P=dat_flt(inds,col_press); %
                    xdat(idat).x = T.*(1000./P).^0.286;
%                    ydat(idat).y=dat_flt(inds,col_alt); %
                    labs(idat).l=[datestr(dat_flt(inds(1),1)/1000/3600,13) ' - ' datestr(dat_flt(inds(end),1)/1000/3600,13) ' UTC'];
 %               end
                
%                 idat=idat+1; inds=indsALL(idat).i;
%                 T=dat_flt(inds,col_temp)+273.15; %
%                 P=dat_flt(inds,col_press); %
%                 xdat(idat).x = T.*(1000./P).^0.286;
%                 ydat(idat).y=dat_flt(inds,col_alt); %
%                 labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(2),1)/1000/3600) ' UTC'];
%                 
                xlims=0;
                xlimits=[270 285];
                xlimits=[271 280];
                
            case 'Vapour mixing ratio (kg kg^{-1}) from frost point hygrometer'
                idat=idat+1; inds=indsALL(idat).i;
                eval_str = ['xdat(idat).x=qv_flt' flight_no '_fp(inds);'];
                eval(eval_str);
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
              
                
            case 'Vapour mixing ratio (kg kg^{-1}) from humicap'
                idat=idat+1; inds=indsALL(idat).i;
                eval_str = ['xdat(idat).x=qv_flt' flight_no '_humi(inds);'];
                eval(eval_str);
%                ydat(idat).y=dat_flt(inds,col_alt); %
                

                                
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
                
             case 'Pressure (hPa)'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=dat_flt(inds,col_press); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                
             case 'Wind speed (m s^{-1})'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=dat_flt(inds,col_wind); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                
                xlims=1;
                xlimits=[0 15];
                
                
             case 'Wind direction (degrees)'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=dat_flt(inds,col_winddir)+180; %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'Total ice number (L ^{-1})'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=interp1(CIP_time_Jonny/3600,1000*ice_no_CIP_Dan,time_flt(inds)); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            
            case 'Mean ice size (microns)'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=interp1(CIP_time_Jonny/3600,mean_ice_size,time_flt(inds)); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'Total CIP mass Jonny (mg m^{-3})'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=1000*interp1(CIP_time_Jonny(1:end-1)/3600,total_mass_Jonny,time_flt(inds)); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            case 'Total CIP number Jonny (cm^{-3})'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=interp1(CIP_time_Jonny(1:end)/3600,ice_no_tot_CIP_Dan,time_flt(inds)); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            case 'Ice+small CIP number Jonny (L^{-1})'
                idat=idat+1; inds=indsALL(idat).i;
                xdat(idat).x=1000*interp1(CIP_time_Jonny(1:end)/3600,sum(ice_PSD+small_PSD,1),time_flt(inds)); %
%                ydat(idat).y=dat_flt(inds,col_alt); %
                labs(idat).l=[num2str(dat_flt(inds(1),1)/1000/3600) ' - ' num2str(dat_flt(inds(end),1)/1000/3600) ' UTC'];
                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                xlims=1;
                xlimits=[0 0.02*1000];
                                
        end
        
        labs(idat).l=[datestr(dat_flt(inds(1),1)/1000/3600/24,13) ' - ' datestr(dat_flt(inds(end),1)/1000/3600/24,13) ' UTC'];
        

    end %for idat2=1:size(times,1)
        
                idat=length(xdat);

                if ihighlight_cloud==1
                    highlight_cloud_method = 'CAS LWC';

                    xdat(3:idat+1)=xdat(2:idat);
                    ydat(3:idat+1)=ydat(2:idat);
                    labs(3:idat+1)=labs(2:idat);
                    idat=1;

                    switch highlight_cloud_method

                        case 'CAS LWC'
                            inds=indsALL(1).i;
                            times=dat_flt(inds,1)/1000;                            
                            indsCAS=find(CAS_time_all>=times(1) & CAS_time_all<=times(end));
                            lwc_threshold=0.01;
                            icloud = find( LWC_dist_cas(indsCAS) > lwc_threshold );
                           
%                            xdat(idat+1).x = NaN*ones(size(times));
                            xdat(idat+1).x = NaN*ones(size(indsCAS)); 
                            ydat(idat+1).y = NaN*ones(size(indsCAS));                             
                            xdat(idat+1).x(icloud) = interp1(times,xdat(1).x,CAS_time_all(indsCAS(icloud)));
                            ydat(idat+1).y(icloud) = interp1(times,ydat(1).y,CAS_time_all(indsCAS(icloud)));
%                            ydat(idat+1).y = dat_flt(inds,col_alt);
                            labs(idat+1).l=['Cloud (CAS>' num2str(lwc_threshold) ' g m^{-3})'];
                    end



                    nmark=zeros(size(xdat));
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                end
        
        
        
        
        
        
        
        
        if izlim==0
            izlim=1;
            zmin=0;
            zmax=max(ydat(1).y);
        end
        
case 88
        %
        
        tstr=date_str;
        iund=findstr('_',tstr);
        tstr(iund)=' ';          
        titlenam = ['Particle separation distribution for ' tstr ' for flight ' flight_no];
        
        figname=titlenam;
        savename=figname;
        
        ixtick_relabel_log=0;
        x_axis_type='';
%        x_axis_type='log10_matlab';
        

        
        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='Log10 of particle separation N';   
%        ylab='Particle separation N';           
%        xlab='???';
        
        instrument_sd_all={'CAS'};
%        instrument_sd_all={'CIP'};
%        instrument_sd_all={'CAS','CIP'};



        lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;    
                          
            
        time_plot=14.695;
        time_plot=14.698;
        time_plot=14.707;
        time_plot=14.7098;
        time_plot=16.6812;
        time_plot=17.3595;
        time_plot=17.4330;
        time_plot=17.48157;
        time_plot=17.582; %nebuliser not working
        time_plot=17.6215; %nebuliser
        time_plot=17.3212; %no beads 2
        time_plot=[19.13 19.33]; %flight 104 19:00
        time_plot=[20.0425 20.0975]; %flight 104 20:00
        time_plot=20.09;
        time_plot=[21.0012 21.109];
        time_plot=[19.3667 19.6333]; %flight l05
        time_plot=[18.2928 18.3216]; %flight l23 - CAS/CIP overlap attempt - good match
        time_plot=[19.846 19.854]; %flight l23 - CAS/CIP overlap attempt2
%        time_plot=[19.86 19.88]; %flight l23 - CAS/CIP overlap - CAS counts test
        time_plot=[18.794 18.804]; %flight l22 - CAS/CIP overlap - CAS counts test
        time_plot=[18.93 18.965]; %flight l22 - CAS/CIP overlap - CAS counts test2        
        time_plot=[18.262 18.324]; %flight l20 - CAS/CIP overlap - CAS counts test  
        time_plot=[18.262 18.273]; %flight l20 - CAS/CIP overlap - CAS counts test 
%        time_plot=[16.472 16.628]; %flight l13 - CAS/CIP overlap - CAS counts test         
%        time_plot=[16.421 16.4325]; %flight l13 - cloud "spike" before the last main broader spike                         
%        time_plot=[19.7984 19.806]; %flight l05 - CAS/CIP overlap - CAS counts test   
        time_plot=[20.7894 20.7902]; %flight 102       
        
%        time_plot=[CAS_time_all(1) CAS_time_all(end-1)]/3600;
        
        savename = [savename ' ' num2str(time_plot) ' ' ylab];
        
        
        clear itime
        
    for idat=1:length(instrument_sd_all)
            instrument_sd=instrument_sd_all{idat};
            
        
        for itime_plot=1:length(time_plot)
            switch instrument_sd
                case 'CAS'
                    itime(itime_plot)=findheight_nearest(CAS_time_all/3600,time_plot(itime_plot));
                case 'CIP'
                    itime(itime_plot)=findheight_nearest(CIP_time_all/3600,time_plot(itime_plot));

            end
        end
        
        time_inds = itime(1):itime(end);
        
        ydat(idat).y = mean(CAS_psep_all(time_inds,:),1);
        xdat(idat).x=[1:size(CAS_psep_all,2)]; %
                  
            switch ylab
                case 'Particle separation N'
                    iytick_relabel_log=0;          
                    y_axis_type='';
                    
                case 'Log10 of particle separation N'    
                    iytick_relabel_log=1;
                    y_axis_type='log10_matlab';
            end
            
            switch length(time_plot)
                case 2
                    labs(idat).l=['Mean ' datestr(time_plot(1)/24,15) '-' datestr(time_plot(2)/24,15) ' UTC ' instrument_sd];  
                otherwise
                    labs(idat).l=[num2str(time_plot) ' UTC ' instrument_sd];  
            end
            
    end        
        
case 87
        %
        if exist('tstr')
            tstr=date_str;
            iund=findstr('_',tstr);
            tstr(iund)=' ';
        else
            tstr='';
            flight_no='';
            date_str='';
        end
        
        ixtick_relabel_log=1;
%         x_axis_type='log10';
%         x_axis_type='log10_matlab';
%         x_axis_type='';

        if length(x_axis_type)==0
            ixtick_relabel_log=0;
        end
        
        
        
        xlims=0;
        xlimits=1000*[0 0.025];
        
        iadd_line=0;
       
        
        izlim=0;
        zmin=1500;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

xlab='Diameter (\mum)';

if ~exist('ichoose_times_size_dists')
    
        iytick_relabel_log=1;       
        y_axis_type='log10_matlab';
        y_axis_type='';
        
        x_axis_type='log10_matlab';
        x_axis_type='';
                
        

        ylab='N (cm^{-3})';
%        ylab='dN/dlogD (cm^{-3} \mum^{-1})';
%        ylab='dN/dD (cm^{-3} \mum^{-1})';
        

        
        instrument_sd_all={'BAS CAS'};
%        instrument_sd_all={'MAN CAS Karl'};   
%        instrument_sd_all={'CDP Karl','MAN CAS Karl'};           
%        instrument_sd_all={'BAS CAS'};        
%        instrument_sd_all={'MAN CAS','BAS CAS'}; 
%        instrument_sd_all={'CAS 29th June 2010'}; %comparison to Manchester CAS
        instrument_sd_all={'CIP'}; %uses CIP_counts_all
%        instrument_sd_all={'CAS','CIP'};
%        instrument_sd_all={'CAS','CAS back','CIP'}; 
%        instrument_sd_all={'CAS','CAS back','CIP_Jonny','CIP','CIP_Ice'};         
        instrument_sd_all={'CIP_Ice_Jonny'}; %uses ice_PSD
%        instrument_sd_all={'Welas'};
%        instrument_sd_all={'FSSP'};        

        instrument_ratio=0;

end

if length(y_axis_type)==0
        iytick_relabel_log=0;
end
        


if instrument_ratio==1
    titlenam = ['Ratio of particle size dists for ' tstr ' for flight ' flight_no];
else
    titlenam = ['Particle size distribution for ' tstr ' for flight ' flight_no];
end
        
        figname=titlenam;
        savename=figname;



        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        lor=1;
        
        idat=0;    
        
         if ~exist('icas_count')
                icas_count=1;
            end

            %%% these are set in plotTimeHeightVap3.m - might want to override here
            % or to make sure they are the same if requried
            if ~exist('cut_off_size')
                cut_off_size=1; %size (microns) below which to ignore counts for particle concentration
            end
            if ~exist('air_speed_type')
                %air_speed_type = 'aircraft';
                %air_speed_type = 'constant 60m/s';
                air_speed_type = 'CIP probe';
            end
            if ~exist('TAS') | ~exist('CIP_time')
                air_speed_type = 'constant';
            else
                if length(TAS)~=length(CIP_time)
                    air_speed_type = 'constant';
                else
                    air_speed_type = 'CIP probe';
                end
            end
            

    %%%    ------------------------------------   %%%
            
           %get the sample volume and total concentrations, plus air speed if required.           
%             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
%                 ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip]...
%                 =cas_sample_volume_and_stats(dat_flt,CAS_time_all...
%                 ,CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all...
%                 ,air_speed_type,cut_off_size,TAS_all,20);   
            
            

if ~exist('ichoose_times_size_dists')
            
        time_plot=14.695;
        time_plot=14.698;
        time_plot=14.707;
        time_plot=14.7098;
        time_plot=16.6812;
        time_plot=17.3595;
        time_plot=17.4330;
        time_plot=17.48157;
        time_plot=17.582; %nebuliser not working
        time_plot=17.6215; %nebuliser
        time_plot=17.3212; %no beads 2
        time_plot=[19.13 19.33]; %flight 104 19:00
        time_plot=[20.0425 20.0975]; %flight 104 20:00
        time_plot=20.09;
        time_plot=[21.0012 21.109];
        time_plot=[19.3667 19.6333]; %flight l05
%         time_plot=[18.2928 18.3216]; %flight l23 - CAS/CIP overlap attempt - good match
%         time_plot=[19.846 19.854]; %flight l23 - CAS/CIP overlap attempt2
% %        time_plot=[19.86 19.88]; %flight l23 - CAS/CIP overlap - CAS counts test
%         time_plot=[18.794 18.804]; %flight l22 - CAS/CIP overlap - CAS counts test
%         time_plot=[18.93 18.965]; %flight l22 - CAS/CIP overlap - CAS counts test2        
%         time_plot=[18.262 18.324]; %flight l20 - CAS/CIP overlap - CAS counts test  
% %        time_plot=[18.262 18.273]; %flight l20 - CAS/CIP overlap - CAS counts test 
% %         time_plot=[16.472 16.628]; %flight l13 - CAS/CIP overlap - CAS counts test         
% %         time_plot=[16.421 16.4325]; %flight l13 - cloud "spike" before the last main broader spike                         
% %        time_plot=[19.7984 19.806]; %flight l05 - CAS/CIP overlap - CAS counts test                 
%         time_plot=[19.35 19.48]; %flight l20 - 19.35 - 19.49 period (19:21 - 19:29) 
% %        time_plot=[19.335+idir*0.005 19.335+(idir+1)*0.005]; %flight l20 - 19.35 - 19.49 period (19:21 - 19:29) 
% %        time_plot=[19.335+idir*0.005]; %flight l20 - 19.35 - 19.49 period (19:21 - 19:29) 
%         time_plot=[19.1104 19.1132];
% %        time_plot=[20.296 20.298];  
%         time_plot=[13.1308 13.1328]; 
%         time_plot=[12.6952 12.6974];     
%         time_plot=[17.3595 17.3603]; %15 micron 
%         time_plot=[17.4330 17.4338]; %30 micron    
%         time_plot=[17.4816 17.4852]; %Lycopodium            
%         time_plot=[17.6176 17.6257]; %Nebulizer
%         time_plot=[17.3203 17.3236]; %No particles   
         time_plot=[16.5412 16.5428]; %31/01/10  15 microns 1
%         time_plot=[16.6154 16.6172]; %31/01/10  15 microns 2
%         time_plot=[16.9284 16.93]; %31/01/10  15 microns 3       
%         time_plot=[16.708 16.7097]; %31/01/10  30 microns 
%         time_plot=[16.733 16.739]; %31/01/10  empty
%         time_plot=[16.4945 16.4963]; %31/01/10  lipo       
%         time_plot=[16.8534 16.86]; %31/01/10  nebu    
%         time_plot=[15.4 15.9]; %flight 99 - CAS/CIP overlap - CAS counts test         
%         time_plot=[14.42 14.9]; %flight 100 - CAS/CIP overlap - CAS counts test  
%         time_plot=[13.72 14.1]; %flight 100 - CAS/CIP overlap - CAS counts test  
%         time_plot=[21.1348 21.1384]; %flight 102 - CAS/CIP overlap - CAS counts test          
%         time_plot=[19.9 20.34]; %flight 102 - CAS/CIP overlap - CAS counts test    
%         time_plot=[21 22]; %flight 108 - CAS/CIP overlap - CAS counts test   
          time_plot=[13+56/60+1/3600 13+57/60+40/3600]; %flight 100 - second droplet period on descent
          time_plot=[13+52/60+46/3600 13+53/60+31/3600]; %flight 100 - second droplet period on ascent
          time_plot=[14+13/60+0/3600 14+19/60+0/3600];
%          time_plot=[welas_dat.time_of_day(1)/3600 welas_dat.time_of_day(end)/3600];
          time_plot=[20+20/60 20+40/60]; %leg 1, flight 102 approx 60-160 km along track
%          time_plot=[20+40/60 21+00/60]; %leg 2, flight 102 approx 155-240 km along track          
          time_plot=[21+00/60 21+20/60]; %leg 3, flight 102 approx 155-240 km along track          
          time_plot=[21+30/60 21+42/60]; %leg 4, flight 102 approx 155-240 km along track   
          time_plot=[20+19.4/60 20+23/60]; %leg 1, first out of cloud portion     
%          time_plot=[20.425 20.445]; %leg 1, first out of cloud portion               
          time_plot=[9.852 10.188; 9.852 10.188]; %leg 1, first out of cloud portion                     
          time_plot=[9.3960 13.6824; 9.3960 13.6824]; %leg 1, first out of cloud portion  
          time_plot=[20.6328 20.6472]; %flight 102 high droplet concs at 20:38 or 146.7 km
%          time_plot=[20.5824 20.5944]; %flight 102 slightly lower droplet concentrations 20:35 or 130.5 km
          time_plot=[13 15.5]; %flight 101 whole flight for CIP dist
          time_plot=[13.5 15.5]; %flight 100 whole flight for CIP dist
          
%          time_plot=[];          
else
    clear ichoose_times_size_dists
end

CAS_LWC_cut_off_sizes=[0 50];
         
        str_times=num2str(time_plot);
        str_times=str_times(:)';
        savename = [savename ' ' str_times ' ' ylab];
        
%        log_bins = log10(CAS_bins);
        
        clear itime
        
    for idat=1:length(instrument_sd_all)
            instrument_sd=instrument_sd_all{idat};
            if size(time_plot,1)>=idat
                str_times_single=num2str(time_plot(idat,:));
            else
                str_times_single='';
            end
%            str_times_single=str_times_single;
            
            switch instrument_sd
                case 'MAN CAS'
                   data_particle = data_CAS_PACS(41:41+29,:)';
                   time_timeseries = data_CAS_PACS(1,:);

                [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_PACS',data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                    
                case 'Welas'
                    data_particle = 1e-6*welas_dat.conc1'; %*1e-6 to convert to cm^-3                                        
                case 'FSSP'
                    disp('*** WARNING - applying FSSP flow speed correction factor ***');
                   data_particle = data_FSSP(23:42,:)'/3.2906;
%            sample_volume_FSSP = 100*airspeed_constant * 0.24e-2; %0.24e-2 is the CAS laser area in cm^2
            
                    time_timeseries = time_FSSP*3600;

                    %just do at first to calculate the sample volume applied
                    [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                            ,CAS_total_number_cutoff ...                                    
                            ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                            ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                            =cas_sample_volume_and_stats2...
                            (0,time_timeseries,...
                            bins_FSSP',data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                            ,CAS_LWC_cut_off_sizes,airspeed_constant);

                        %now compensate for the applied sample volume
                        %since we already have concentrations for FSSP
                    data_particle_scaled = data_particle.*sample_volume_CAS; 
                        %this only applies for the products of cas_sample_volume_and_stats2
                        %for plotting size distributions, etc. just use data_particle
                    
                    %and recalculate the other stuff   
             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_FSSP',data_particle_scaled,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                
                %but want to use data_particle for plotting distributions, etc
                    
                case 'BAS CAS'
                    time_timeseries = CAS_time_all;
                    data_particle = CAS_counts_all;
            
%                     [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
%                         ,CAS_total_number_cutoff ...                                    
%                         ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
%                                 ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
%                         =cas_sample_volume_and_stats2...
%                         (dat_flt,time_timeseries,...
%                         CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,airspeed_constant); 
                    
                    
                    [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (dat_flt,time_timeseries,...
                    CAS_bins,data_particle,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                
                case  'MAN CAS Karl'
                
                    time_timeseries = CAS_time_Karl; %seconds from 0 UTC
                    data_particle = CAS_per_cc;
                    sample_volume_CAS = ones(size(data_particle));
                    
                 case  'CDP Karl'
                
                    time_timeseries = CDP_time_Karl; %seconds from 0 UTC
                    data_particle = CDP_per_cc;
                    sample_volume_CDP = ones(size(data_particle));
                 
                    
            end
                                
            
         if size(time_plot,1)>=idat
            time_plot_idat=time_plot(idat,:);
         else 
            time_plot(idat,:)=[time_timeseries(1) time_timeseries(end)];
            time_plot_idat=time_plot(idat,:)/3600;
         end
         
         for itime_plot=1:size(time_plot_idat,2)
            switch instrument_sd
                case {'CAS','CAS back'}
                    itime(itime_plot)=findheight_nearest(CAS_time_all/3600,time_plot_idat(itime_plot));
                case 'CIP'
                    itime(itime_plot)=findheight_nearest(CIP_time_all/3600,time_plot_idat(itime_plot));
                case {'CIP_Jonny','CIP_Ice_Jonny'}
                    itime(itime_plot)=findheight_nearest(CIP_time_Jonny/3600,time_plot_idat(itime_plot));                    
                case 'MAN CAS'
                    itime(itime_plot)=findheight_nearest(data_CAS_PACS(1,:)'/3600,time_plot_idat(itime_plot));
                case 'Welas'
                    itime(itime_plot)=findheight_nearest(welas_dat.time_of_day'/3600,time_plot_idat(itime_plot));    
                case 'FSSP'
                    itime(itime_plot)=findheight_nearest(time_FSSP,time_plot_idat(itime_plot));
                case {'BAS CAS'}
                    itime(itime_plot)=findheight_nearest(CAS_time_all/3600,time_plot_idat(itime_plot));    
                case {'MAN CAS Karl'}
                    itime(itime_plot)=findheight_nearest(CAS_time_Karl/3600,time_plot_idat(itime_plot));                        
                case {'CDP Karl'}
                    itime(itime_plot)=findheight_nearest(CDP_time_Karl/3600,time_plot_idat(itime_plot));                                            

            end
         end
        
        
        time_inds = itime(1):itime(end);
        
        icalc_Nlimits=0;
        
        if icalc_Nlimits==1
        
%calculte the number of particles within a size limit        
        
         Nlim_sizes=[0.6 1.04]; %calculate the number in between these limits                
         Nlim_sizes=[0.6 10]; %calculate the number in between these limits                         
                Ndrops = mean( data_particle(time_inds,1:end)./sample_volume_CAS(time_inds,1:end) )';
                if Nlim_sizes(1)==0
                    ilims_lower=1;
                else
                    ilims01 = find(CAS_bins>=Nlim_sizes(1));
                    ilims_lower = ilims01(1)+1; %add one as the first number bin is for D<0.61
                end
                ilims02 = find(CAS_bins<=Nlim_sizes(2));                
                ilims_upper = ilims02(end); %bin N+1 is for sizeN to sizeN+1. Want N for 
                                            %sizes<sizeN so is bin N that we require
                Nlims = sum(Ndrops(ilims_lower:ilims_upper));

                %example of displaying the number info
%               text(0.84,1.5,['Total 0.61-1.03 \mum = ' num2str(Nlims) ' cm^{-3}'],'fontsize',14)

%%%%% end of size limit particle number calculation %%%%%%

        end

                  
            switch ylab
                case 'dN/dlogD (cm^{-3} \mum^{-1})'
                    
                    
                    
                    switch instrument_sd
                        case 'MAN CAS'
                                                 
                            log_bins = log10(bins_PACS');
                            dlogD=repmat(diff(log_bins),[1 size(data_particle,1)])';
                            ydat(idat).y = mean(data_particle(time_inds,2:end)./dlogD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);                                                       
                            xdat(idat).x=(bins_PACS(2:end)'+bins_PACS(1:end-1)')/2;  %
                            

                        case 'CIP'
                            CIP_bins2 = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries
                            log_bins = log10([CIP_bins2; CIP_bins2(end)+25]);
                            dlogD=repmat([diff(log_bins)],[1 size(CIP_counts_all,1)])';
                            ydat(idat).y = mean(CIP_counts_all(time_inds,:)./dlogD(time_inds,:)./sample_volume_CIP(time_inds,:),1);
                            
                            xdat(idat).x=(CIP_bins); %CIP bins are already the mid-points
                            lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                                                            
                            
                         case 'CIP_Ice_Jonny'
                            CIP_bins2 = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries
                            log_bins = log10([CIP_bins2; CIP_bins2(end)+25]);                            
                            dlogD=repmat([diff(log_bins)],[1 size(ice_PSD,2)]);
                            ydat(idat).y = mean(ice_PSD(:,time_inds)./dlogD(:,time_inds),2);
                            xdat(idat).x=(CIP_bins); %CIP bins are already the mid-points
%                            lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                            fsize=10;                              
                            
                        case 'Welas'


                            sample_volume_CAS = ones(size(data_particle)); 
                            refrac_factor=1.2375; %the value used for glass beads for the fog instrument
                            %from the straight line between the first and last point.
                            refrac_factor=1.1;
                            refrac_factor=1;
                            
                            orig_bins = welas_dat.size(:,1)/refrac_factor; %is the size form the calibration files
                            orig_bins = water2u_105;
                            %the radius or diameter? What is the conversion for droplets from latex beads?
                            
                            
                            log_bins = log10(orig_bins);
                            dlogD = repmat(diff(log_bins)',[size(data_particle,1) 1]);

                            ydat(idat).y = mean(data_particle(time_inds,2:end)./dlogD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);                                                       
                            xdat(idat).x =(orig_bins(2:end)'+orig_bins(1:end-1)')/2;  %
                            
                            
                            %set nfilter to one for no smoothing
                            nfilter=1; bfilter=ones([1 nfilter])*1/nfilter;
                            ydat(idat).y=filter(bfilter,1,ydat(idat).y);                            
                            xdat(idat).x(1:nfilter-1)=[]; ydat(idat).y(1:nfilter-1)=[];
                            
                            
                            
                        case 'FSSP'

                            sample_volume_CAS = ones(size(data_particle)); 
                            log_bins = log10(bins_FSSP);
                            dlogD = repmat(diff(log_bins),[size(data_particle,1) 1]);

                            ydat(idat).y = mean(data_particle(time_inds,2:end)./dlogD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);                                                       
                            xdat(idat).x=(bins_FSSP(2:end)'+bins_FSSP(1:end-1)')/2;  %
                            
                       case 'MAN CAS Karl'    
                            log_bins = log10(CAS_bins_Karl);
                            dlogD = repmat(diff(log_bins),[size(data_particle,1) 1]);                            
                            ydat(idat).y = mean(data_particle(time_inds,2:end)./dlogD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);
                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CAS_bins_Karl(2:end)+CAS_bins_Karl(1:end-1))/2; %
                            
                            
%                             ydat(idat).y = mean( data_particle(time_inds,2:end)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here
% 
%                             %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
%                             xdat(idat).x=(bins_PACS(2:end)'+bins_PACS(1:end-1)')/2;  %
                            
                        case 'CDP Karl'
                            log_bins = log10(CDP_bins_Karl);
                            dlogD = repmat(diff(log_bins),[size(data_particle,1) 1]);                            
                            ydat(idat).y = mean(data_particle(time_inds,1:end)./dlogD(time_inds,:)./sample_volume_CDP(time_inds,1:end),1);
                            xdat(idat).x=(CDP_bins_Karl(2:end)+CDP_bins_Karl(1:end-1))/2; %    
                            
                            %for the CDP we only have one count per bin (per mid-point) whereas with CAS
                            %there is an extra count for the counts less than the first bin edge
                            
                        

                           
                        otherwise   
                            log_bins = log10(CAS_bins);
                            dlogD=repmat(diff(log_bins),[1 size(CAS_counts_all,1)])';
                            ydat(idat).y = mean(CAS_counts_all(time_inds,2:end)./dlogD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);
                            xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2; %                        
                    end
                    
                    
                    
                    xlims=0;
                    xlimits=[0 50];                    
%                    xlimits=[20 50];
               
                    
                case 'dN/dD (cm^{-3} \mum^{-1})'
                    switch instrument_sd
                        case 'CAS'
                            dD=repmat(diff(CAS_bins),[1 size(CAS_counts_all,1)])';
                            ydat(idat).y = (1/1)*mean(CAS_counts_all(time_inds,2:end)./dD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);
                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2; %
                        case 'CAS back'
                            dD=repmat(diff(CAS_bins_back),[1 size(CAS_back_all,1)])';
                            ydat(idat).y = (1/1)*mean(CAS_back_all(time_inds,2:end)./dD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);
                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CAS_bins_back(2:end)+CAS_bins_back(1:end-1))/2; %                            
                        case 'CIP'
                            CIP_bins2 = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries
                            dD=repmat([diff(CIP_bins2); 25],[1 size(CIP_counts_all,1)])';
                            ydat(idat).y = mean(CIP_counts_all(time_inds,:)./dD(time_inds,:)./sample_volume_CIP(time_inds,:),1);
                            xdat(idat).x=(CIP_bins); %CIP bins are already the mid-points
                            lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                            
                        case 'CIP_Jonny'
                            CIP_bins2 = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries
                            dD=repmat([diff(CIP_bins2); 25],[1 size(CIP_PSD_Jonny,1)])';
                            ydat(idat).y = mean(CIP_PSD_Jonny(time_inds,:)./dD(time_inds,:),1);
                            xdat(idat).x=(CIP_bins); %CIP bins are already the mid-points
                            lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                            fsize=10;
                        case 'CIP_Ice_Jonny'
                            CIP_bins2 = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries
                            dD=repmat([diff(CIP_bins2); 25],[1 size(ice_PSD,2)]);
                            ydat(idat).y = mean(ice_PSD(:,time_inds)./dD(:,time_inds),2);
                            xdat(idat).x=(CIP_bins); %CIP bins are already the mid-points
%                            lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                            fsize=10;    
                            
                        case 'MAN CAS Karl'    
                            dD=repmat(diff(CAS_bins_Karl),[1 size(data_particle,1)])';
                            ydat(idat).y = (1/1)*mean(data_particle(time_inds,2:end)./dD(time_inds,:)./sample_volume_CAS(time_inds,2:end),1);
                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CAS_bins_Karl(2:end)+CAS_bins_Karl(1:end-1))/2; %
                            
                            
%                             ydat(idat).y = mean( data_particle(time_inds,2:end)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here
% 
%                             %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
%                             xdat(idat).x=(bins_PACS(2:end)'+bins_PACS(1:end-1)')/2;  %
                            
                        case 'MAN CDP Karl'
                            dD=repmat(diff(CDP_bins_Karl),[1 size(data_particle,1)])';
                            ydat(idat).y = (1/1)*mean(data_particle(time_inds,2:end)./dD(time_inds,:)./sample_volume_CDP(time_inds,2:end),1);
                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CDP_bins_Karl(2:end)+CDP_bins_Karl(1:end-1))/2; %
                            
                            
                            
                            
%                            ydat(idat).y = mean( data_particle(time_inds,2:end)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here

                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
 %                           xdat(idat).x=(CAS_bins(2:end)'+ CAS_bins(1:end-1)')/2;  %    
                            
                    end
                    
                    
                    xlims=1;
                    xlimits=[0.6 1000];
  
                case 'N (cm^{-3})'
                    switch instrument_sd
                        case 'MAN CAS'
                                                 
                     
                            ydat(idat).y = mean( data_particle(time_inds,2:end)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here

                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(bins_PACS(2:end)'+bins_PACS(1:end-1)')/2;  %
                            
                        case 'MAN CAS Karl'
                            ydat(idat).y = mean( data_particle(time_inds,2:end)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here

                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CAS_bins(2:end)'+ CAS_bins(1:end-1)')/2;  %
                            

                        case 'CDP Karl'
                            ydat(idat).y = mean( data_particle(time_inds,1:end)./sample_volume_CDP(time_inds,1:end) )'; %note the transpose ' here

                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(CDP_bins_Karl(2:end)'+ CDP_bins_Karl(1:end-1)')/2;  %                            
                            
                         case 'CIP_Ice_Jonny'
                            CIP_bins2 = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries                            
                            ydat(idat).y = mean(ice_PSD(:,time_inds),2);
                            xdat(idat).x=(CIP_bins); %CIP bins are already the mid-points
%                            lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                            fsize=10;   
                            
                            
                        otherwise                            
                            ydat(idat).y = mean( CAS_counts_all(time_inds,1:end)./sample_volume_CAS(time_inds,1:end) )'; %note the transpose ' here
                            xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2;  %    
                            xdat(idat).x=[0 xdat(idat).x'];
                    end
                    
%                    ydat(idat).y = mean( CAS_counts_all(time_inds,1:end-1)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here                    
%                    xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2;  %
                    
                    
                    labs(idat).l=[str_times_single ' UTC'];   
                    
                    
                    xlims=0;
                    xlimits=[0 50];
%                    xlimits=[0 10];                      
                    
                case 'LWC (g m^{-3})'
                    switch instrument_sd
                        case 'MAN CAS'
                                                 
                     
                            ydat(idat).y = mean(LWC_size_dist(time_inds,:) )'; %note the transpose ' here

                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=(bins_PACS(2:end)'+bins_PACS(1:end-1)')/2;  %
                           
                        otherwise                            
                            ydat(idat).y = mean( LWC_size_dist(time_inds,:))'; %note the transpose ' here
                            xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2;  %                           
                    end
                    
%                    ydat(idat).y = mean( CAS_counts_all(time_inds,1:end-1)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here                    
%                    xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2;  %
                    
                    
                    labs(idat).l=[str_times_single ' UTC'];   
                    
                    
                    xlims=1;
                    xlimits=[0 50];
                    
                    
                case 'dLWC/dlogD (g m^{-3} \mum^{-1}))'
                    switch instrument_sd
                        case 'MAN CAS'
                                                 
                            log_bins = log10(bins_PACS');
                            dlogD=repmat(diff(log_bins),[1 size(data_particle,1)])';
                            
                            ydat(idat).y = mean(LWC_size_dist(time_inds,:) ./(dlogD(time_inds,:)) )'; %note the transpose ' here

                            %                    xdat(idat).x=(log_bins(2:end)+log_bins(1:end-1))/2; %
                            xdat(idat).x=((bins_PACS(2:end)'+bins_PACS(1:end-1)')/2);  %
                           
                        otherwise      
                            log_bins = log10(CAS_bins);
                            dlogD=repmat(diff(log_bins),[1 size(CAS_counts_all,1)])';
                            
                            ydat(idat).y = mean( LWC_size_dist(time_inds,:)./(dlogD(time_inds,:)) )'; %note the transpose ' here
                            xdat(idat).x= ((CAS_bins(2:end)+CAS_bins(1:end-1))/2);  %                           
                    end
                    
%                    ydat(idat).y = mean( CAS_counts_all(time_inds,1:end-1)./sample_volume_CAS(time_inds,2:end) )'; %note the transpose ' here                    
%                    xdat(idat).x=(CAS_bins(2:end)+CAS_bins(1:end-1))/2;  %
                    
                    
                    labs(idat).l=[str_times_single ' UTC'];   
                    
                    
                    xlims=1;
                    xlimits=[0 50];
                    
            end
            
            if size(time_plot,1)==1
                idat_str=1;
            else
                idat_str=idat;
            end
            
            if length(time_plot)>1
                    labs(idat).l=[datestr(time_plot(idat_str,1)/24,13) '-' datestr(time_plot(idat_str,2)/24,13) ' ' instrument_sd];  
            else
                    labs(idat).l=[datestr(time_plot(1)/24,13) ' UTC ' instrument_sd];  
            end
            
            
            if instrument_ratio==1
                ratio_dat(idat).dat=ydat(idat).y;

                if idat==2
                    ydat(1).y=ratio_dat(2).dat./ratio_dat(1).dat;
                    ydat(2)=[];
                    xdat(2)=[];
                    labs(2)=[];
                    labs(1).l='Ratio';
                    titlenam=['Ratio ' instrument_sd_all{2} ' divided by ' instrument_sd_all{1} ' for ' datestr((time_plot(2,1)-1)/24,13) '-' datestr((time_plot(2,2)-1)/24,13) ' and ' datestr(time_plot(1,1)/24,13) '-' datestr(time_plot(1,2)/24,13)];
%                    y_axis_type='';
                    xlims=1;
                    xlimits=[0 50];
                    ylab='Ratio';
                end
            end
                    
                                
            
    end
    
    if iadd_line==1
        min_line=9e9;
        max_line=-9e9;
        for iline=1:length(xdat)
            min_line = min(min_line,min(ydat(iline).y));
            max_line = max(max_line,max(ydat(iline).y));
        end
        min_line = ydat(1).y(9)/1e1;
        max_line = ydat(1).y(9)*1e1;        
        addlineX=[xdat(1).x(9);xdat(1).x(9)];
        addlineY=[min_line;max_line];
        
        min_line = ydat(1).y(21)/1e1;
        max_line = ydat(1).y(21)*1e1;        
        addlineX(:,2)=[xdat(1).x(21);xdat(1).x(21)];
        addlineY(:,2)=[min_line;max_line];
    end

        
      
            
        
    case 86
        % profiles from the cross section - NEED TO HAVE RUN plotTime... first for cross section
        
        time_prof=0;
%        time_prof=6;

if time_prof==0
     tstr=Times(time,:);
else
     tstr=Times(time_prof,:);
end
        
       
        iund=findstr('_',tstr);
        tstr(iund)=' ';          
        titlenam = ['Potential temperature cross section profile for ' tstr];
        titlenam = ['Wind speed cross section profile for ' tstr];
%        titlenam = ['N profile for ' tstr];
%        titlenam = ['L profile for ' tstr];

        
        figname=titlenam;
        savename=figname;

        
        izlim=0;
        zmin=0.55;
        zmax=2.5;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        
%        HGT=zz(1).z;

if time_prof==0
        HGT=zz(1).z;
        dist=timesTH(1).t;
else
        HGT=XY_pot_cross_data(1).Y_cross;
        dist=XY_pot_cross_data(1).X_cross;
end
        

        
        switch titlenam
            case ['L profile for ' tstr]
                pot=pot_cross_15UTC;
                U_cross = U_cross_15UTC;
                
                z=1000*repmat(HGT',[1 size(pot,2)]);
                DAT = sqrt( 9.81 ./ pot(2:end,:) .* diff(pot,1)./diff(z,1) ) ./ U_cross(2:end,:); %N 
                HGT = HGT(2:end);
                
                ylab='Height (km)';
                xlab= 'L (m^{-1})';
                xlims=0;
                xlimits=[271 295];

                
            case ['N profile for ' tstr]
                pot=pot_cross_15UTC;
                z=1000*repmat(HGT',[1 size(pot,2)]);
                DAT = sqrt( 9.81 ./ pot(2:end,:) .* diff(pot,1)./diff(z,1) ); %N 
                HGT = HGT(2:end);
                
                ylab='Height (km)';
                xlab= 'N (s^{-1})';
                xlims=0;
                xlimits=[271 295];

                                
            case ['Potential temperature cross section profile for ' tstr]
                ylab='Height (km)';
                xlab= 'Potential temperature (K)';
                xlims=1;
                xlimits=[271 295];
                
                switch time_prof
                    case 6
                        DAT=pot_cross_15UTC; %time_prof=6
                    case 11
                        DAT=pot_cross;  %time_prof=11
                    case 0
                        DAT=pdat(1).p;
                end
                


            case ['Wind speed cross section profile for ' tstr]
                xlab= 'Wind speed (m s^{-1})';
                xlims=1;
                xlimits=[0 12];
                
                switch time_prof
                    case 6
                        DAT=U_cross_15UTC; %time_prof=6
                    case 11
                        DAT=u_cross;   %time_prof=11
                    case 0
                        DAT=pdat(1).p;
                end
                

        end


        lor=2; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        

        



        
        ipos=1;
        idat=idat+1;
        xdat(idat).x=DAT(:,ipos); %
        ydat(idat).y=HGT; %
        labs(idat).l=[num2str(dist(ipos),'%.0f') ' km'];
        
%         ipos=3;
%         idat=idat+1;
%         xdat(idat).x=DAT(:,ipos); %
%         ydat(idat).y=HGT; %
%         labs(idat).l=[num2str(dist(ipos),'%.0f') ' km'];
%         
%         ipos=20;
%         idat=idat+1;
%         xdat(idat).x=DAT(:,ipos); %
%         ydat(idat).y=HGT; %
%         labs(idat).l=[num2str(dist(ipos),'%.0f') ' km'];
%         
        ipos=76;
        idat=idat+1;
        xdat(idat).x=DAT(:,ipos); %
        ydat(idat).y=HGT; %
        labs(idat).l=[num2str(dist(ipos),'%.0f') ' km'];  
        
        ipos=113;
        idat=idat+1;
        xdat(idat).x=DAT(:,ipos); %
        ydat(idat).y=HGT; %
        labs(idat).l=[num2str(dist(ipos),'%.0f') ' km'];  
        
    case 85
        % Ice Nucleation schemes
        
%        tstr=Times(time,:);
%        iund=findstr('_',tstr);
%        tstr(iund)=' ';   
        IN_type = 'IN concentrations';
%        IN_type = 'Ice procuction rates';
        
        iydir=-1;
        
        switch IN_type
            case 'IN concentrations'
                titlenam = ['Ice nuclei concentations'];

                figname=titlenam;
                savename=figname;

                xlims=0;
                xlimits=[-18 0];

                izlim=1;
                zmin=-20;
                zmax=-5;

                nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

                xlab='Heterogeneous IN number cocentration (L^{-1})';
                ylab= 'Temperature (^{o}C)';



                lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;

                tc=[-30:0];

                idat=idat+1;
                xdat(idat).x=0.005*exp(0.304*(-tc)); %per L
                ydat(idat).y=tc; %
                labs(idat).l='WRF (Coooper)';

                %         idat=idat+1;
                %         qv=satvappress(273+tc,'goff','liq');
                %         qvi=satvappress(273+tc,'goff','ice');
                %         xdat(idat).x=exp(-6.69+12.96*qv./qvi)*1e-3; %IS WRONG - qv should be vapour mixing ratio vapour pressure
                %         %and not the saturation vapour pressure
                %         ydat(idat).y=tc; %
                %         labs(idat).l='Meyers';

                idat=idat+1;
                xdat(idat).x=0.01*exp(-0.6*tc)*1e-3; %per L
                ydat(idat).y=tc; %p
                labs(idat).l='Fletcher';
                
                a=0.0000594;
                b=3.33;
                c=0.0264;
                d=0.0033;
                TK=tc+273.15;
                
                naer=0.1; %number of aerosol between 0.5 and 1 micron
                idat=idat+1;
                xdat(idat).x=a*(273.16-TK).^b.*naer.^(c.*(273.16-TK)+d);
                ydat(idat).y=tc; %p
                labs(idat).l=['DeMott Naer= ' num2str(naer)];
                
                naer=0.3; %number of aerosol between 0.5 and 1 micron
                idat=idat+1;
                xdat(idat).x=a*(273.16-TK).^b.*naer.^(c.*(273.16-TK)+d);
                ydat(idat).y=tc; %p
                labs(idat).l=['DeMott Naer= ' num2str(naer)];
                
                naer=3; %number of aerosol between 0.5 and 1 micron
                idat=idat+1;
                xdat(idat).x=a*(273.16-TK).^b.*naer.^(c.*(273.16-TK)+d);
                ydat(idat).y=tc; %p
                labs(idat).l=['DeMott Naer= ' num2str(naer)];
                
               


                
                

            case 'Ice procuction rates' %e.g. Bigg's and contact that produce a rate of ice number formation

                

               

                xlims=0;
                xlimits=[-18 0];

                izlim=1;
                zmin=-20;
                zmax=-5;

                nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

                xlab='Ice number production rate (L^{-1} s^{-1})';
                ylab= 'Temperature (^{o}C)';

                lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                idat=0;

                tc=[-30:0];
                P=650; %
                P=600; %            
                QC3D=0.3; %g/m3
                [final_contact_rate, rate_bigg] = microphysics_WRF_contact_and_Biggs(tc,P,QC3D);
                
                titlenam = ['Ice number formation rates for P=' num2str(P) ' mb, QC=' num2str(QC3D) ' g m^{-3}'];
                figname=titlenam;
                savename=figname;

                idat=idat+1;
                xdat(idat).x=rate_bigg; %per L per sec
                ydat(idat).y=tc; %
                labs(idat).l='WRF Biggs';

                idat=idat+1;
                xdat(idat).x=final_contact_rate; %per L per sec
                ydat(idat).y=tc; %p
                labs(idat).l='Contact freezing';
                
                idat=idat+1;
                xdat(idat).x=final_contact_rate+rate_bigg; %per L per sec
                ydat(idat).y=tc; %p
                labs(idat).l='Biggs + contact freezing';



        end
        
    case 84
        % Froude number for continuous stratification
                  
        titlenam = ['Froude vs L*thi number for constant stratification'];
        
        figname=titlenam;
        savename=figname;

        xlims=0;
        xlimits=1000*[0 0.025];
        
        izlim=1;
        zmin=0;
        zmax=2;


        nmark=0; %-1 means that all points have markers. Otherwise only plot the number specified.

        ylab='F';
        xlab='0.5*L*thi';


        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        idat=0;
        
        L=0.01/10;
        
        hhat=1;
        del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
        H0=( hhat - del_hat + acos(hhat./del_hat) )/L;
        hm=hhat/L;
                    
                    
%        H0=3*pi/2/L;
%        h=1/L;
        
%        x=[0:0.01:2*pi]; %=L*thi
        x=[0:0.01:L*H0]; %=L*thi
        h=[0:hm/100:hm];
        Ld=x+L.*(h-H0); %=L*del from del=H0+del-h

        idat=idat+1;
        xdat(idat).x=0.5*x; %0.5*L*thi
        ydat(idat).y=sqrt( (1-Ld.*sin(x))./(1-cos(x)) ); %F as calculated from the formulation in Durran 
        labs(idat).l='F';
                            
                            
                            
                           
        
    case 83
        % WRF microphysics for ice heteorogeneous freezing, Bigg's immersion droplet freezing and contact nucleation
        time=15;
        
        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';  
        
        titlenam = ['WRF Morrison microphysics for ' tstr];
        figname=titlenam;
        savename=figname;
        
        
        
                    xlims=0;
                    xlimits=1000*[0 0.025];
                    
                    nmark=-1;
                    
                    iloc=1;
                            
                            ylab='Height (m)';
                            xlab= 'WRF ice number concentration (L^{-1})';
                                      
                            izlim=1;
                            zmin=000;
                            zmax=3000;
                            
                            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                            
                            idat=0;
                            HGT=WRFUserARW(nca(1).nc,'Z',time,ilat(iloc),ilon(iloc));
                                                                                                 
                            
                            T3D = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                            P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;                           
                            RHO=density(P,T3D);

                            QC3D=nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc)); %cloud MR kg/kg
                            NDCNST=250; %currently set constant (250 per cc) in version 3.0.1.1 of WRF
                            NC3D=NDCNST.*1.E6./RHO; %convert to #/kg

                            tc=T3D-273.15;      
                            
                            %%% contact nucleation
                            idat=idat+1;
                            xdat(idat).x=0.005*exp(0.304*(-tc)); %per L
                            ydat(idat).y=HGT; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='Heterogeneous IN';
                            
                            
                            
                            idat=idat+1;
                            RIN=0.1e-6;
                            CONS37=4.*pi*1.38e-23/(6.*pi*RIN);
                            RHOW=997;
                            CONS26=pi/6.*RHOW;                            
                            
                            
                            DUM2 = 1.496e-6*T3D.^1.5./(T3D+120);
                            MU = DUM2./RHO;
                            NACNT=exp(-2.80+0.262*(-tc))*1000; %per m3
                            DUM = 7.37.*T3D./(288.*10.*P)/100;
                            DAP = CONS37.*T3D.*(1.+DUM/RIN)./MU;
                            
                            DUM3 = P./(287.15*T3D);
                            PGAM=0.0005714*(NC3D./1.E6./DUM3)+0.2714;
                            PGAM=1./(PGAM.^2)-1;
                            PGAM=max(PGAM,2);
                            PGAM=min(PGAM,10);

                            CDIST1 = NC3D./gamma(PGAM+1);
                            LAMC=(CONS26.*NC3D.*GAMMA(PGAM+4)./(QC3D.*GAMMA(PGAM+1))).^(1/3);
                            
                            % LAMMIN, 60 MICRON DIAMETER
                            % LAMMAX, 1 MICRON
                            
                            LAMMIN = (PGAM+1.)/60.e-6;
                            LAMMAX = (PGAM+1.)/1.e-6;

                            iless=find(LAMC<LAMMIN);
                            LAMC(iless) = LAMMIN(iless);
                            NC3D(iless)= exp(3.*log(LAMC(iless))+log(QC3D(iless))+    ...
                                log(GAMMA(PGAM(iless)+1.))-log(GAMMA(PGAM(iless)+4.)))/CONS26;
                            
                            imore=find(LAMC>LAMMAX);
                            LAMC(imore) = LAMMAX(imore);
                            NC3D(imore) = exp(3.*log(LAMC(imore))+log(QC3D(imore))+  ...
                                log(GAMMA(PGAM(imore)+1.))-log(GAMMA(PGAM(imore)+4.)))/CONS26;                            


%NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
%                    GAMMA(PGAM(K)+2.)/                         &
%                    LAMC(K)
                            
                            %this is a rate in #/kg/s
                            rate_contact=2.*pi*DAP.*NACNT.*CDIST1.*gamma(PGAM+2)./LAMC;
                            DT=0.5*60;
                            DT=500;
                            xdat(idat).x = rate_contact*DT.*RHO/1000; %work out for arbitrary time
                            ydat(idat).y=HGT; 
                            labs(idat).l='Contact nucleation';

                            AIMM=0.66;                            
                            BIMM=100;
                            CONS40=pi/6.*BIMM;                           
                            
                            %%% immersion nucleation (of cloud droplets to form ice) - Bigg's
                            idat=idat+1;
                            rate_bigg=CONS40*exp(log(CDIST1)+log(gamma(PGAM+4))-3*log(LAMC)).*exp(AIMM*(273.15-T3D));
                            xdat(idat).x=rate_bigg*DT.*RHO/1000;
                            ydat(idat).y=HGT; 
                            labs(idat).l='Bigg''s immersion nucleation';
                            %% N.B. doesn't depend on droplet number as cancels out in rate_bigg formula
                            %% (from formulae for PGAM, CDIST1 and LAMC)
                            
                            %% aside - mass of ice frozen from Biggs does depend on droplet size
                            CONS39=pi*pi/36.*RHOW*BIMM;
                            MNUCCC = CONS39*exp(log(CDIST1)+log(gamma(7.+PGAM))-6.*log(LAMC)).*exp(AIMM*(273.15-T3D));

        
                                
                                
                               
    case 82
        
                   
                    titlenam = 'Houghton/Smith';
        
                    xlims=0;
                    xlimits=1000*[0 0.025];
                                      
                    gd=1;
                    U=20;
                    N=0.02;

                    hhat=[0:0.01:1];

                    %continuous stratification case (Fig. 15b)
                    del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
                    H0_hat_strat=hhat - del_hat + acos(hhat./del_hat);
                    H0_crit=U/N*H0_hat_strat; %
                    F0_crit=1 ./ H0_hat_strat; %=U/NH=F0=1/H0_hat for Fig 15b


                    %single layer case (Fig. 15a)

                    H0_crit2 = [0:1:1500000];  %N.B. - require very large H0 values to reach high h/H0 values (to approach the limit
                    %towards h/H0=1,F0=0
                    F0_crit2=U./sqrt(gd.*H0_crit2);
                    F0_crit3 = U/N./H0_crit2;
                    h_crit2 = H0_crit2 .* ( 1 + 0.5*F0_crit2.^2 - 1.5*F0_crit2.^(2/3) );
                    h_crit3 = H0_crit2 .* ( 1 + 0.5*F0_crit3.^2 - 1.5*F0_crit3.^(2/3) );
                    hhat_Houghton=h_crit3*N/U;
                    
                    plot_case='stratified';
                    plot_case='h/H0,F0 space';
%                    plot_case='h/H0,F0 space clean';
%                    plot_case='h vs H0';                    
%                    plot_case='h vs H0 dual';                                        
                    
                    switch plot_case
                        case 'stratified'

                            xdat(1).x=hhat_Houghton;
                            ydat(1).y=F0_crit3; %using Smith formula for F0 but with F0=U/(N*H0) - acheived by approximating
                            %g' with N^2*H0 where H0=dz in N^2=g/theta * dtheta/dz
                            labs(1).l='Smith using N';

                            xdat(2).x=hhat;
                            ydat(2).y=F0_crit; %plot the original formulation for continuous stratification
                            % are different but.... if multiply F0 by a factor of 0.62...
                            labs(2).l='Original Smith';

                            xdat(3).x=hhat_Houghton;
                            ydat(3).y=F0_crit3*0.62; %they overlay each other almost exactly. Why 0.62??? - epsilon? =0.622
                            %is independant of N and U
                            labs(3).l='Using N and *0.62';
                            
                            xlab='F0=U/(N*H0)';
                            ylab='hN/U';


                        case 'h/H0,F0 space'
                            idat=0;
                            
                            idat=idat+1;
                            xdat(idat).x=h_crit2./H0_crit2;
                            ydat(idat).y=F0_crit2; %using Smith formula for F0 for an inversion F0=U/sqrt(g'H0)
                            labs(idat).l='Original inversion case';
                            
                            %now give range of h' values and calculate H0 using (20). Note F0=1/H0' since H0' = H0*N/U
                            idat=idat+1;                            
                            xdat(idat).x=2*hhat./H0_hat_strat;
                            ydat(idat).y=sqrt(2)*F0_crit; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='Stratified formula with Smith&Sun transforms';
                            % i.e. using Heff instead of H0 and F0=U/(N*Heff)=sqrt(2)U/(N*H0)

                            %        plot(1/(0.62^2)*hhat./H0_hat_strat,1/0.62*F0_crit,'k--');
                            
                            idat=idat+1;
                            xdat(idat).x=hhat./H0_hat_strat/0.63;
                            ydat(idat).y=F0_crit/0.63;
                            labs(idat).l='Stratified formula with 0.63 transforms';
                            
%                             idat=idat+1;                            
%                             xdat(idat).x=hhat./H0_hat_strat/0.63;
%                             ydat(idat).y=F0_crit/0.63;
%                             labs(idat).l='Stratified formula with 0.63,0.86 transforms';
                            
                            idat=idat+1;                            
                            xdat(idat).x=hhat./H0_hat_strat;
                            ydat(idat).y=1./(sqrt(2)*sin(0.5*H0_hat_strat));
                            labs(idat).l='Stratified formula using Durran F0';
                            
                            idat=idat+1;                            
                            xdat(idat).x=hhat./H0_hat_strat;
                            ydat(idat).y=F0_crit;
                            labs(idat).l='Stratified formula using Durran F0';
                            
                            
                            
%                             idat=idat+1;                            
%                             xdat(idat).x=2*hhat./H0_hat_strat;
%                             ydat(idat).y=1./(sqrt(2)*sin(0.5*H0_hat_strat));
%                             labs(idat).l='Stratified formula using Durran F0 and Heff';
                                                
                            x_start=0.5/0.63; %=0.7937 approx h/Heff for Antarctic Peninsula (with h=1500, H0=3000)
                            y_start=0.18/0.63; %=0.2857 - using N and U for AP case, and H0=3000 (U/N = 0.18*3000 = 540)
                            b=[0:1500]; %values of effective blocking to try

%                             idat=idat+1;                            
%                             xdat(idat).x=(1500-b)./(0.63*(3000-b)); %increasing values of blocking height
%                             ydat(idat).y=540./(3000-b)/0.63;
%                             labs(idat).l='Bee line for AP';


    %%% jump boundary
%                             idat=idat+1;                            
%                             [F,h_h0]=mountain_Houghton_solve_cr_0_line(3500,0.5);
%                             xdat(idat).x=h_h0;
%                             ydat(idat).y=F;
%                             labs(idat).l='Jump boundary';                            
                                                                                    
                            ylab='F_0';
                            xlab='h/H_0';
                            
                            xlims=1;
                            xlimits=[0 1];
                                                        
                            izlim=1;
                            zmin=0;
                            zmax=1;
                            
                        case 'h vs H0 dual'
                            Ha_hat=3*pi/6;
                            Ha_hat=1;
                            
                            Ha_hat=4;
                            
                            dmin=1.01*atan(1/Ha_hat);  %minimum d(=Hb_hat-Ha_hat) allowed from eqn (28) of Smith and Sun
                            Hb_start=Ha_hat+dmin;
                            Hb_hat=[Hb_start:(9*pi/6-Hb_start)/100:9*pi/6];

                            
                            idat=0;
                            
                            idat=idat+1;
                            for iHb=1:length(Hb_hat);
                                [xdat(idat).x(iHb),dA,dB]=Smith_dual_layer_find_hpeak(Ha_hat,Hb_hat(iHb),1);
                            end
                            ydat(idat).y=Hb_hat;  %calculated above
                            labs(idat).l='Stratified exact solution';
                            
                            %now calculate effective F0 and plug into eq (16) for 0.63 transform
                            idat=idat+1;
                            D=0.5;
                            Heff=Ha_hat + D * (Hb_hat-Ha_hat);
                            r=(Hb_hat-Ha_hat)/Ha_hat;
                            Feff=sqrt( (1+r).^2./(Hb_hat.^2.*r.*(1+D*r)) );
%                            Feff=F0_crit3;
                            h_crit2 = Heff .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(idat).x=h_crit2;
                            ydat(idat).y=Hb_hat; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='Heff transform into inversion formula';
                            
                            
                            %attempt at using the 0.63 idea where the change of potential temperature for g' is taken as only that over Heff
                            %so that g'=g*dtheta/theta=N^2*d*0.63 (or another factor instead of 0.63)
                            idat=idat+1;
                            D=0.63;
                            Heff=Ha_hat + D * (Hb_hat-Ha_hat);
                            Feff=sqrt( 1./(Heff.*(Heff-Ha_hat)) );
%                            Feff=F0_crit3;
                            h_crit2 = Heff .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(idat).x=h_crit2;
                            ydat(idat).y=Hb_hat; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='Heff transform into inversion formula';
                            
                            
                            
%continuosly stratified effective F0 into eq (16) for 0.63 transform
                            idat=idat+1;
                            Feff=1/0.63 * U/N./H0_crit2;
%                            Feff=F0_crit3;
                            h_crit2 = H0_crit2*0.63 .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(idat).x=h_crit2*N/U;
                            ydat(idat).y=H0_crit2*N/U; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='0.63 transform into inversion formula';
                            
%continuosly stratified effective F0 into eq (16) for 0.63 transform using Hb_hat
                            idat=idat+1;
                            D=0.63;
                            Feff=1/D * 1./Hb_hat;
%                            Feff=F0_crit3;
                            h_crit2 = Hb_hat*D .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(idat).x=h_crit2;
                            ydat(idat).y=Hb_hat; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='0.63 transform into inversion formula using Hb_hat';                            
                            
                            
                            %Smith-Sun transform
                            idat=idat+1;
                            Feff=sqrt(2)*U/N./H0_crit2;
                            h_crit2 =  H0_crit2*0.5 .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(idat).x=h_crit2*N/U;
                            ydat(idat).y=H0_crit2*N/U; %plot the original formulation for continuous stratification with transformed variables
                            labs(idat).l='Smith-Sun transform into inversion formula';
                            % i.e. using Heff instead of H0 and F0=U/(N*Heff)=sqrt(2)U/(N*H0)

                            %        plot(1/(0.62^2)*hhat./H0_hat_strat,1/0.62*F0_crit,'k--');

%                            xdat(3).x=hhat./H0_hat_strat/0.63;
%                            ydat(3).y=F0_crit/0.63;
%                            labs(3).l='Stratified formula with 0.63 transforms';
                            
                            ylab='Hhat';
                            xlab='hhat';
                            
                            xlims=1;
                            xlimits=[0 1];
                                                        
                            izlim=1;
                            zmin=0;
                            zmax=5;
                            
                            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                       
                            case 'h vs H0'
                            xdat(1).x=hhat;
                            ydat(1).y=H0_hat_strat;  %calculated above
                            labs(1).l='Stratified exact solution';
                            
                            %now calculate effective F0 and plug into eq (16) for 0.63 transform
                            Feff=1/0.63 * U/N./H0_crit2;
%                            Feff=F0_crit3;
                            h_crit2 = H0_crit2*0.63 .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(2).x=h_crit2*N/U;
                            ydat(2).y=H0_crit2*N/U; %plot the original formulation for continuous stratification with transformed variables
                            labs(2).l='0.63 transform into inversion formula';
                            
                            %Smith-Sun transform
                            Feff=sqrt(2)*U/N./H0_crit2;
                            h_crit2 = 0.5 * H0_crit2 .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(3).x=h_crit2*N/U;
                            ydat(3).y=H0_crit2*N/U; %plot the original formulation for continuous stratification with transformed variables
                            labs(3).l='Smith-Sun transform into inversion formula';
                            % i.e. using Heff instead of H0 and F0=U/(N*Heff)=sqrt(2)U/(N*H0)

                            %        plot(1/(0.62^2)*hhat./H0_hat_strat,1/0.62*F0_crit,'k--');

%                            xdat(3).x=hhat./H0_hat_strat/0.63;
%                            ydat(3).y=F0_crit/0.63;
%                            labs(3).l='Stratified formula with 0.63 transforms';
                            
                            ylab='Hhat';
                            xlab='hhat';
                            
                            xlims=1;
                            xlimits=[0 1];
                                                        
                            izlim=1;
                            zmin=0;
                            zmax=5;
                            
                            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                            
                       case 'h/H0,F0 space clean'
                           
                            xdat(1).x=h_crit2./H0_crit2;
                            ydat(1).y=F0_crit2; %using Smith formula for F0 for an inversion F0=U/sqrt(g'H0)
                            labs(1).l='Critical boundary';  
                            
                            xdat(1).x(end+1)=1;
                            ydat(1).y(end+1)=0;
                            
                            [F,h_h0]=mountain_Houghton_solve_cr_0_line(3500,0.5);
                            xdat(2).x=h_h0;
                            ydat(2).y=F;
                            labs(2).l='Jump boundary';    
                            
                            xdat(3).x=hhat./H0_hat_strat;
                            ydat(3).y=1./(sqrt(2)*sin(0.5*H0_hat_strat)); %using Durran's formula for F0
                            labs(3).l='Durran F'; 
                                                                                    
                            ylab='F0';
                            xlab='h/H0';
                            
                            xlims=1;
                            xlimits=[0 1.3];
                                                        
                            izlim=1;
                            zmin=0;
                            zmax=2.4;
                            
%%% adding labelled points onto graph - Durran, etc.      
                            add_points=1;
                            ipoint=0;
%Miller and Durran, 1991, Fig. 5.
                            N=0.02094;
                            U=20;
                            L=N/U;
                            hm=1146.1;
                            
                            %fzero(@mountain_Smith_hhat_for_blocking,[0.01 1],[],L,h_mountain,Haim)
                            %gives an effective mountain height for a given actual height of Haim
                            %with the assumption of as much blocking as is needed
                            %have used this below to find effective values for different streamlines
                            %assuming that the flow can select whatever mountain height it likes                                                         
                            
                            ipoint=ipoint+1;
                            H0=3.6249e+003;
                            xpos(ipoint).x=571.0775/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='M4.2km';
                            
%Durran, 1987 Part 2, Fig 1a,b and c - applying 0.63 transforms
                            N=0.01047;
                            U=20;
                            L=N/U;
                            H0=6.0e3; %divinding streamline is approx 5 km
                        %(a) Lh=0.3
                            ipoint=ipoint+1;
                            xpos(ipoint).x=0.3/L/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='D2a';
                        %(a) Lh=0.4
                            ipoint=ipoint+1;                        
                            xpos(ipoint).x=0.4/L/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='D2b';
                        %(a) Lh=0.5
                            ipoint=ipoint+1;                        
                            xpos(ipoint).x=0.5/L/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='D2c';    
                            
%Durran, 1987 Part 2, Fig 1a,b and c - NOT applying 0.63 transforms - using Durran F0
                            N=0.01047;
                            U=20;
                            L=N/U;
                            H0=6.0e3; %divinding streamline is approx 5 km
                            F0=1/(sqrt(2)*sin(0.5*L*H0));
                            
                        %(a) Lh=0.3
                            ipoint=ipoint+1;
                            xpos(ipoint).x=0.3/L/H0;
                            ypos(ipoint).y=F0;
                            point_labs(ipoint).lab='D2a';
                        %(a) Lh=0.4
                            ipoint=ipoint+1;                        
                            xpos(ipoint).x=0.4/L/H0;
                            ypos(ipoint).y=F0;
                            point_labs(ipoint).lab='D2b';
                        %(a) Lh=0.5
                            ipoint=ipoint+1;                        
                            xpos(ipoint).x=0.5/L/H0;
                            ypos(ipoint).y=F0;
                            point_labs(ipoint).lab='D2c';                                
                            
%Durran, 1987 Part1, Fig. 16
                            N=0.01;
                            U=10;
                            L=N/U;
                            H0=6e3;
                            hm=1400;
                            
                            ipoint=ipoint+1;                            
                            xpos(ipoint).x=hm/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='D1';
                            
                            %fzero(@mountain_Smith_hhat_for_blocking,[0.01 1],[],L,h_mountain,Haim)
                            %gives an effective mountain height for a given actual height of Haim
                            %with the assumption of as much blocking as is needed
                            %have used this below to find effective values for different streamlines
                            %assuming that the flow can select whatever mountain height it likes                                                         
                            
                            H0=4.4861e+003;
                            ipoint=ipoint+1;                            
                            xpos(ipoint).x=800/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='D5km';
                         
                            H0=2.8539e+003;
                            ipoint=ipoint+1;                            
                            xpos(ipoint).x=253.8/H0/0.63;
                            ypos(ipoint).y=1/L/H0/0.63;
                            point_labs(ipoint).lab='D4km';
                            
                            
                            

                    end
                    
                    
        case 81    %plot of aircraft distance vs. wind speed to show variation
        
        [asL23_1 asL23_2]= findheight(time_flt19,20.76,21.55); %L-shaped sections of level flight ignoring when wind speed went wrong (too high)
        [asL23_3 asL23_4]= findheight(time_flt19,21.06,21.23);
        [asL23_5 asL23_6]= findheight(time_flt19,21.6,21.96);
        inds_L3=asL23_4-2500:asL23_2+3500;   %3rd level L-shaped flight sections                
        inds_L3=asL23_5:asL23_6+0;   %3rd level L-shaped flight sections                

    xvar='dist';
%    xvar='lon';
    switch xvar
        case 'lon'           
            xlab = 'Longitude';
            xdat(1).x = dat_flt19(inds_L3,3);
        case 'dist'
            np = 100; %number of points to split the flight path into
            inds = round([inds_L3(1):(inds_L3(end)-inds_L3(1))/np:inds_L3(end)]);
            LAT_plot = dat_flt19(inds,2);
            LON_plot = dat_flt19(inds,3);

            [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_plot,LON_plot,0.1);

            x_vals_flight = (ilon-1)*dx_grid;
            y_vals_flight = (ilat-1)*dy_grid;
                       
            dist_flight = cumsum(sqrt((diff(x_vals_flight)).^2 + (diff(y_vals_flight)).^2)); %find distances between each point and the last and sum cumulatively
            time_dist_flight=time_flt19(inds(2:end));
            
            xdat(1).x = interp1(time_dist_flight,dist_flight,time_flt19(inds_L3),'','extrap'); 
            xlab = 'Distance (km)';
    end    
    
    plotcase='wind';

    
iref = findheight(LAT,-65.58); %Larsen B point (approx centre of Larsen B)

    
    switch plotcase
        case 'wind'    
            
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            
            
            idat=1;
            
            ydat(idat).y = dat_flt19(inds_L3,9);                  
            labs(idat).l = ['L3'];
            labs(idat).l = ['L4'];
            

            
            
%            ylab = 'Mean melt rate contribution (mm day^{-1})';
%            figname=['Mean melt rate contributions for ' filestr];  
            
            ylab = 'Wind Speed (m s^{-1})';
            figname=['Wind speed vs. distance or Lon for ' filestr]; 
            
            izlim=0;
            zmin=-8;
            zmax=2;
              
        
              

            
    end
    
               
    
    titlenam=figname;
    savename=[figname];
    
    xlims=1;
    xlimits=[min(xdat(1).x) max(xdat(1).x)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    case 80    %plot against latitude
%run calc_melt_tot first then, which will also run heat_fluxes_mean_along_latitude.m

    latlon='lon';
    switch latlon
        case 'lat'           
            xlab = 'Longitude';
            xdat(1).x = LON;
        case 'lon'
            xlab = 'Latitude';
            xdat(1).x = LAT;
    end    
    
    plotcase='normal';
%    plotcase='relative_diffs';
%    plotcase='cloud';
%    plotcase='RH';
%    plotcase='n_tot';

offset_using_ref_point = 'yes'; %for melt, sw and lw take away the values at the ref location
%offset_using_ref_point = 'no'; %or don't if no
    
iref = findheight_nearest(LAT,-65.58); %Larsen B point (approx centre of Larsen B)

    
    switch plotcase
        case 'normal'    
            
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            
            
            idat=1;
            
            switch offset_using_ref_point
                case 'yes'
                    ydat(idat).y = melt_dat - melt_dat(iref);
                case 'no'
                    ydat(idat).y = melt_dat;
            end
            xdat(idat).x = LAT;            
            labs(idat).l = ['Melt - ' num2str(melt_dat(iref),'%.2f')]; idat=idat+1;
            
            switch offset_using_ref_point
                case 'yes'
                    ydat(idat).y = sw_dat - sw_dat(iref);
                case 'no'
                    ydat(idat).y = sw_dat;
            end
            xdat(idat).x = LAT;            
            labs(idat).l = ['SW - ' num2str(sw_dat(iref),'%.2f')]; idat=idat+1;
            
            switch offset_using_ref_point
                case 'yes'
                    ydat(idat).y = lw_dat - lw_dat(iref);
                case 'no'
                    ydat(idat).y = lw_dat;
            end
            xdat(idat).x = LAT;            
            labs(idat).l = ['LW + ' num2str(-lw_dat(iref),'%.2f')]; idat=idat+1;
            
            ydat(idat).y = sh_dat; 
            xdat(idat).x = LAT;            
            labs(idat).l = 'SH'; idat=idat+1;
            
            ydat(idat).y = lh_dat; 
            xdat(idat).x = LAT;
            labs(idat).l = 'LH'; idat=idat+1;
            
            ydat(idat).y = sh_dat+lh_dat; 
            xdat(idat).x = LAT;
            labs(idat).l = 'SH+LH'; idat=idat+1;
            
            ydat(idat).y = grd_dat; 
            xdat(idat).x = LAT;
            labs(idat).l = 'GRD'; idat=idat+1;
            
            
%            ylab = 'Mean melt rate contribution (mm day^{-1})';
%            figname=['Mean melt rate contributions for ' filestr];  
            
            ylab = 'Melt contribution (mm)';
            figname=['Melt contributions for ' filestr]; 
            
            izlim=0;
            zmin=-8;
            zmax=2;
              
        case 'relative_diffs'                  

            idat=1;            
            
            ydat(idat).y = (sw_dat-sw_dat(iref))./(melt_dat-melt_dat(iref)); 
            xdat(idat).x = LAT;            
            labs(idat).l = 'SW'; idat=idat+1;
            
            ydat(idat).y = (lw_dat-lw_dat(iref))./(melt_dat-melt_dat(iref));  
            xdat(idat).x = LAT;            
            labs(idat).l = 'LW'; idat=idat+1;
            
            ydat(idat).y = (sh_dat-sh_dat(iref))./(melt_dat-melt_dat(iref)); 
            xdat(idat).x = LAT;            
            labs(idat).l = 'SH'; idat=idat+1;
            
            ydat(idat).y = (lh_dat-lh_dat(iref))./(melt_dat-melt_dat(iref)); 
            xdat(idat).x = LAT;
            labs(idat).l = 'LH'; idat=idat+1;
            
            shlh_dat = sh_dat + lh_dat;
            ydat(idat).y = (shlh_dat-shlh_dat(iref))./(melt_dat-melt_dat(iref)); 
            xdat(idat).x = LAT;
            labs(idat).l = 'SH+LH'; idat=idat+1;
            
            ydat(idat).y = (grd_dat-grd_dat(iref))./(melt_dat-melt_dat(iref)); 
            xdat(idat).x = LAT;
            labs(idat).l = 'GRD'; idat=idat+1;
            
            
            ylab = 'Relative contribution to mean melt rate (mm day^{-1})';
            figname=['Relative contribution to mean melt rate for ' filestr];    
            
            izlim=1;
            zmin=-1;
            zmax=1.5;
            
        case 'cloud'    
            
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            
            
            idat=1;
            
            ydat(idat).y = 1000*cond_dat; %convert to g/m2 of condensate
            xdat(idat).x = LAT;            
            labs(idat).l = ['Cond']; idat=idat+1;                      
            
            
            ylab = 'Total condensate (g m^{-2})';
            figname=['Latitude mean total condensate for ' filestr];   
            
        case 'RH'    
            
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            
            
            idat=1;
            
            ydat(idat).y = rh_dat; %convert to g/m2 of condensate
            xdat(idat).x = LAT;            
            labs(idat).l = ['RH']; idat=idat+1;                      
            
            
            ylab = 'Relative humidity';
            figname=['Latitude mean relative humidity for ' filestr];   
            
         case 'n_tot'    
            
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            
            
            idat=1;
            
            ydat(idat).y = n_tot_dat; %convert to g/m2 of condensate
            xdat(idat).x = LAT;            
            labs(idat).l = ['Days']; idat=idat+1;                      
            
            
            ylab = 'Mean no. days';
            figname=['Latitude mean number of melting days for ' filestr];      
              

            
    end
    
    
    ichoose_styles=1;  % '-'=solid, '--'=dashed, ':'=dotted, '-.'=dash-dot
istyle=1;    
       line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '--';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '-'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='^'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '-.';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;
        line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0.7 0.7]; marker_style(istyle).m='o'; line_widths(istyle).l = lwidth; istyle=istyle+1;        
    
               
    
    titlenam=figname;
    savename=[figname];
    
    xlims=1;
    xlimits=[xdat(1).x(1) xdat(1).x(end)];
    xlimits=[-69 xdat(1).x(end)];    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 79
        
    time=14;
    
    tstr=Times(time,:);
    iund=findstr('_',tstr);
    tstr(iund)=' '; 
    
    
    xlims=0;
    
    latlon='lon';
    switch latlon
        case 'lat'
            lat=-67.5;
            icon_inds = get_inds_constant_lat(lat,lat2d,lon2d); %get indices for a constant latitude slice
            xlab = 'Longitude';
            labs(1).l = ['LAT = ' num2str(lat)];
            xdat(1).x = lon2d.var(icon_inds);
        case 'lon'
            
            lon=-61.25;
%            lon=-62.5;
            icon_inds = get_inds_constant_lat(lon,lon2d,lat2d); %for lon just swap around input of lat and lon
            xlab = 'Latitude';
            labs(1).l = ['LON = ' num2str(lon)];
            xdat(1).x = lat2d.var(icon_inds);
    end
    
    variable='SWDOWN';
    variable='Melt';
%    variable='10m Wind Speed';
%    variable='n-level Wind Speed';
    variable='Terrain';    
%     variable='TSK';   
%     variable='TSLB';   
%     variable='TSLB';
%     variable='SNOWDEN';       
    
    switch variable
        case 'Melt'
            LW=nc{'GLW'}(time,:);
            SW=nc{'SWDOWN'}(time,:); %downwelling SW
            SH=-nc{'HFX'}(time,:);
            LH=-nc{'LH'}(time,:); %negative as WRF convention is that these are fluxes into air
            ALBEDO=nc{'ALBEDO'}(time,:);
            EMISS=nc{'EMISS'}(time,:);
            GRDFLX=nc{'GRDFLX'}(time,:); %sign convention is that this is flux into the surface layer
            TSK=nc{'TSK'}(time,:);
            %    ALBEDO=0.78;
            %sign convention is postive means energy going into ground
            SW_UP=ALBEDO.*SW; %think albedo means this for WRF
            SW_NET=SW-SW_UP;

            %EMISS=0.98;
            LW_UP=5.67e-8.*TSK.^4; %Boltzmann W/m2 using skin temperature (need emissivity?)
            LW_NET=EMISS.*(LW-LW_UP);


             MNET=LW_NET + SW_NET + LH+SH + GRDFLX;
%            MNET=LW_NET;
            %        MNET=LH+SH;

            
            melt_var='MNET';            
%            melt_var='SW_NET';
%            melt_var='LW_NET';
%            melt_var='SH+LH';
            melt_var='GRDFLX';
%            melt_var='SH';
%            melt_var='LH';
%            melt_var='SH+LH+GRDFLX';            
            
            switch melt_var
                case 'MNET'
                    slice_dat = MNET;
                    zmin=0; 
                    zmax=200;
                    izlim=1;
                    figname=['Melt flux at ' tstr ' for ' filestr];
                case 'SW_NET'
                    slice_dat = SW_NET;
                    zmin=230; 
                    zmax=280;
                    izlim=1;
                    figname=['SW net flux at ' tstr ' for ' filestr];
                case 'LW_NET'
                    slice_dat = LW_NET;
                    zmin=230; 
                    zmax=280;
                    izlim=0;
                    figname=['LW net flux at ' tstr ' for ' filestr];
                case 'SH+LH'
                    slice_dat = LH+SH;
                    zmin=-40; 
                    zmax=10;
                    izlim=1;
                    figname=['Sensible + latent heat flux at ' tstr ' for ' filestr];  
                case 'GRDFLX'
                    slice_dat = GRDFLX;
                    zmin=-40; 
                    zmax=10;
                    izlim=0;
                    figname=['Ground heat flux at ' tstr ' for ' filestr];     
                case 'SH'
                    slice_dat = SH;
                    zmin=-30; 
                    zmax=20;
                    izlim=1;
                    figname=['Sensible heat flux at ' tstr ' for ' filestr];    
                 case 'LH'
                    slice_dat = LH;
                    zmin=-30; 
                    zmax=0;
                    izlim=1;
                    figname=['Latent heat flux at ' tstr ' for ' filestr];     
                    
                case 'SH+LH+GRDFLX'
                    slice_dat = SH+LH+GRDFLX;
                    zmin=-90; 
                    zmax=-20;
                    izlim=1;
                    figname=['Sensible + latent + ground heat flux at ' tstr ' for ' filestr];  
                    
                
                    
            end

            
            ylab = 'Flux (W m^{-2})';
            

            

            

        case 'SWDOWN'

            slice_dat = nc{'SWDOWN'}(time,:);
            ylab = 'Flux (W m^{-2})';
            figname=['Shortwave downwelling radiation at ' tstr ' for ' filestr];
            
    case '10m Wind Speed'            
       u10 = nc{'U10'}(time,:); %10 m winds
       v10 = nc{'V10'}(time,:); %        
       slice_dat = sqrt(u10.^2+v10.^2);
        
       ylab = 'Wind Speed (m s^{-1})';
       figname=['10m wind speed at ' tstr ' for ' filestr];
       
    case 'n-level Wind Speed'    
       ih_wrf=1;
       u = 0.5* (nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) ); %2d wind at one height
       v = 0.5* (nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) ); %2d wind at one height 
       slice_dat = sqrt(u.^2+v.^2);
        
       ylab = 'Wind Speed (m s^{-1})';
       figname=['Wind speed at model level ' num2str(ih_wrf) ' at ' tstr ' for ' filestr];   
       
    case 'Terrain'

            slice_dat = nc{'HGT'}(time,:);
            ylab = 'Terrain height (m)';
            figname=['Terrain height at ' tstr ' for ' filestr];
            
      case 'TSK'

            slice_dat = nc{'TSK'}(time,:)-273.15;
            ylab = 'Skin temperature (^{o}C)';
            figname=['Skin temperature at ' tstr ' for ' filestr];      

            
       case 'TSLB'

            slice_dat = nc{'TSLB'}(time,1,:)-273.15;
            ylab = 'Soil temperature (^{o}C)';
            figname=['Soil temperature at ' tstr ' for ' filestr];      
     
        case 'SNOWDEN'

            SNOWH=nc{'SNOWH'}(time,:); %snow depth m - varies by a massive amount - something wrong?
            SNOW=nc{'SNOW'}(time,:); %
            
            SNODEN=SNOW./SNOWH; %snow density in kg/m3 - varies quite a lot as time goes on
            SNODEN(SNODEN>400)=400; %limited to 400 kg/m3 in the code
            
            slice_dat = SNODEN;
            ylab = 'Snow density (kg m^{-3})';
            figname=['Snow density at ' tstr ' for ' filestr];      
              

    end
    
    
    
    ydat(1).y = slice_dat(icon_inds);
    
    
    
    titlenam=figname;
    savename=[figname labs(1).l];
    
    xlims=1;
    xlimits=[xdat(1).x(1) xdat(1).x(end)];
    

case 78
    %different WRF plots

    dual=0;
    lor=-1;
    
    no_sort=0; %flag to stop sorting in height
       
    %num of markers - set to -1 for all data points to have a marker
    nmark=[-1 -1 -1 -1 -1 -1]; %if give an array then these apply for the different lines
    nmark=0; %set to zero for no markers
    
    Nlocs=0;
    
    incep=0;
    
    i_multi_wrf=1; %flag to say that want to plot from more than one WRF output file
    ilabel_rundir=0; %flag to say to label the legend according to the file directory
    %othewise uses the A,B,C etc. location labels.
    i_paper_labels=0; %flag to say whether to write the paper-friendly titles or the full titles with date, time and simulation name 
    i_label_time=0; %flag to label lat/lon (=0) or not(=1)
    no_title=0; %switch off the title

    
    %%%%%%%%%%%%%%  choose things to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aircraft_comp = 1; %flag to say if want the aircraft comparison
    roth_comp = 0;  %flag to say if want Rothera sounding plotted
    iwrf_profs = 1; %flag to say whether to plot wrf profiles or not
    %probably don't want to change this - try setting aircraft_comp=1
     % ascent_str='034'; %'4' is the L-shaped segments - variation of L-shaped segs is minimal
     ascent_str='03'; %0 and 3 are the first and last ascent
%     ascent_str='1'; %upwind side
    % ascent_str='03'; %5 for L-shaped
    %ascent_str='6'; %mini descent for the last ascent

    iplot_aircraft_locs_model=1; %flag to say whether we want the aircraft ascent and descent 
    %locations to be plotted for the model (aircraft locs specified in LAT_extra shortly)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear location_lab
    
    %order of these labels is: aircraft ascent/descent, etc.; specified LAT/LON values; specified extra_x/y values
    i_override_loc_lab=1; %set to one to use the labels below - make sure there are enough labels
    location_lab(1).l = 'A';
    location_lab(2).l = 'B';  %'M'
    location_lab(3).l = 'C';
    location_lab(4).l = 'D';
    location_lab(5).l = 'WJ';
    location_lab(6).l = 'WK';
    
    
    
%    for iloc_lab=1:10
%        location_lab(iloc_lab).l = 'J';
%    end 
    
%    location_lab(2).l = '-68, -75';
%    location_lab(3).l = '-67, -75';
    
    
    

 
    izlim=1;
    zmin=0;
    zmax=12000;
    zmax=3000;
%    zmax=9000;
        



    %	xlims=0;
    %	xlimits=[980 995];
    %	xlimits=[0 1.5e5];


    logflag=0;
 %   iydir = -1;



 switch file_type
     case 'wrfout'
         time_array=[8 9 10 11 12 17 18 19]; %6th Jan: 11=06, 12=09, 13=12, 14=15, 15=18, 16=21 UTC
         time_array=[10 11 12 12];
        time_array=[13];
%        time_array=[14];
%        time_array=[15];
%         time_array=idir;
     case 'met_em'
         time_array=1;
 end
 
 


        var='pressure';
    var='temperature';
%    var='equiv_potemp';
%    var='potemp';
%     var='vapour';
%    var='ice';
%    var='RH';
     var='wind speed';
%     var='wind speed component';
%     var='cloud';
     var='wind dir';
%     var='Froude'; 
%      var='density gradient';
%      var='dpot/dz';
%      var='westerly wind'         
%       var='stratification integral'
%        var='Scorer parameter';
%        var='ridge_height';
%        var='Houghton Smith formulae';
        
%    LAT=[-67.5702 -67.1420];
%    LON=[-68.1297 -61.6650];


    ylab='Height (m)';

%stores locations along the flight track so that the same locations can be plotted for the model or for 
%plotting on the map (stored in LAT_extra)
    if aircraft_comp == 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%   points to plot on the flight track based on the time along the flight track
        times_flight_loc = [19.84565 20.056 20.12 20.2755]; %20.12 is the time where the max was seen for aircraft data - others are just points
        %along the aircraft track
        times_flight_loc = [20.3755 21.9815 22.0205]; %20.3755 is the time where the max was seen for aircraft data (index=269899)
        %21.9815 is the time of second wind peak on the mini descent before the last ascent
        times_flight_loc = [20.3755 22.0205]; 

%        times_flight_loc = [22.0054]; %the time of the lowest point of the mini-descent
        
        clear it_flt;
        for iflt=1:length(times_flight_loc)
            it_flt(iflt) = findheight(time_flt19,times_flight_loc(iflt));
        end
        LAT_extra = dat_flt19(it_flt,2)';
        LON_extra = dat_flt19(it_flt,3)';

    end
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%% can set LAT LON here %%%%%%%%%%%
% LAT_extra = 24.08;
% LON_extra = 38.06;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
LAT_extra2 =[];  %the locations specified by extra_x and extra_y will be stored in LAT/LON_extra2
LON_extra2 =[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%   extra points to plot - give as extra_x and y km

% extra_x = [575]; %for 12UTC, 6th Jan, ncep polar
% extra_y = [351];
 
% extra_x = [556]; %for 03UTC, 7th Jan, ncep polar
% extra_y = [355];
 
% extra_x = [570]; %for 12UTC, 6th Jan, ecmwf
% extra_y = [355];

%  extra_x = [600]; %random 
%  extra_y = [400];
  
%  extra_x = [275]; %lefthand side of the equiv cross sections
%  extra_y = [380];

%extra_x = [605 800]; %for 18UTC, 6th Jan, ecmwf_ml_0.5
%extra_y = [275 200];
%  
%  extra_x = [634]; %for 21UTC, 6th Jan, ncep_3dom_nudging
%  extra_y = [271]; 
% % 
%  extra_x = [640]; %for 21UTC, 6th Jan, ecmwf_ml_0.5_nudging - 67, 62.2
%  extra_y = [320];
%  
   extra_x = [587 654.48 625 525 675 592 685]; %for 12UTC, 6th Jan, ecmwf_ml_0.5_nudging
   extra_y = [325 296.87 120 225 325 222 80];
 
%%%%%%%%%%%%%%%%%     17th Feb, 2010    %%%%%%%%%%%%%%%%
   extra_x = [557 512]; %for 12UTC, 6th Jan, ecmwf_ml_0.5_nudging
   extra_y = [307 240]; %points E and F - will relabel C and D as of 27th Oct.
   
%   extra_x = [682 626]; %for 21UTC, 6th Jan, ecmwf_ml_0.5_nudging
%   extra_y = [270 277]; %middle of jet beyond the edge of the ice shelf - points C and D
   
%   extra_x = [595]; %for 21UTC, 6th Jan, ecmwf_ml_0.5_nudging
%   extra_y = [282]; %middle of jet beyond the edge of the ice shelf - points C and D

%   extra_x = [605]; %06/04/10 - For 15 UTC - the point in the jet centre 
%   extra_y = [300]; %

%   extra_x = [583 550]; %13/09/10 - For the comparison to analysis at 12 UTC on 6th
%   extra_y = [298 260]; %

%   extra_x = [557 512 530]; %13/09/10 - For the comparison to analysis at 12 UTC on 6th
%   extra_y = [307 240 265]; %
%   Nlocs=[2 1];

%    extra_x = [646 620];
%    extra_y = [298 332];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   extra_x = [459 557 512]; %01/03/10 - 12 UTC on 6th - looking at a point from down south
%   extra_y = [146 307 240]; %to investigave the southerly wind properties
   
%   extra_x = [472 605 565]; %01/03/10 - 15 UTC on 6th - looking at a point from down south
%   extra_y = [135 300 312]; %to investigave the southerly wind properties

%   extra_x = [482 605 600]; %01/03/10 - 21 UTC on 6th - looking at a point from down south
%   extra_y = [212 300 220]; %to investigave the southerly wind properties


  %       extra_x = [extra_x 624.37]; %Larsen AWS location (domain 3)
  %       extra_y = [extra_y 283.125];
  
%   
%   extra_x = [587 654.48 625 525 675 592 685 570 690 100]; %for 12UTC, 6th Jan, ecmwf_ml_0.5_nudging
%   extra_y = [325 296.87 120 225 325 222 80 275 10 400];

%extra_x = [680];  %21 UTC 6th Jan ecmwf_ml_0.5_nudging case - strongest wind region to the east of the aircraft location
%extra_y = [270];  %N.B. using the new (proper) DX and DY grid sizes. (15th Feb, 2010)

%    extra_x = [568 700]; %for 12 & 21 UTC, 6th Jan, ecmwf_ml_wps_snow
%    extra_y = [292 291];
    
%    extra_x = [626 590];   % 633    % 580];   %651]; %for 18UTC, 6th Jan, ncepR2_seaice
%    extra_y = [255 254];   % 275    % 274];   %281];
    
%    extra_x = [522.75 485.5];   % 514]; %Avery 7th December case B,C and D
%    extra_y = [407 407];  % 418];
    
%    extra_x = [479]; %Avery 1st Dec case
%    extra_y = [396.5];
  
%    extra_x = [522.75];
%    extra_y = [407];
    
%    extra_x = [485.5];
%    extra_y = [407];

%    extra_x=[300 350];
%    extra_y=[400 350];
    
    
%    [x_cross,y_cross]=find_xy_position_along_cross_section_line(x_line,y_line,100);
    
%    extra_x=[x_line(1) x_cross 1];
%    extra_y=[y_line(1) y_cross 450];
    
    
    
%    extra_x=[250 235.5 10 450];
%    extra_y=[375 413 435 360];

%    extra_x=[200];
%    extra_y=[500];

    
 %   extra_x=[290];
 %   extra_y=[430];

%  extra_x=[];
%  extra_y=[];  


  
% 

% extra_x = [587]; %for 12UTC, 6th Jan, ecmwf_ml_0.5_nudging
% extra_y = [325];

%  extra_x = [200 300]; %extra points to the westside of the peninsula
%  extra_y = [600 500];






% check whether any have been added to the list of LAT/LON values too

 nlat = size(lat2d.var,1);
 nlon = size(lat2d.var,2);
 
 i_grid = DX * [0:nlon-1];
 j_grid = DY * [0:nlat-1];
 
 
 %now calculate the lat/lon of the extra_x and extra_y points and put them in LAT_extra2
 for iflt=1:length(extra_x)
     i_extra = findheight_nearest(i_grid,extra_x(iflt));
     j_extra = findheight_nearest(j_grid,extra_y(iflt));
     LAT_extra2(iflt) = lat2d.var(j_extra,i_extra);
     LON_extra2(iflt) = lon2d.var(j_extra,i_extra);
 end
 
 


 
       
     
     LAT=[-67.55 -67.62 -67.55]; %first ascent from Rothera
     LON=[-68.1 -67.8 -67.5];
     
%     LAT=[-67.55 -67.55 -67.55]; %Rothera base - put in twice so can show two different times for the same location
%     LON=[-68.1 -68.1 -68.1];
     
%      LAT=[-68.2031 -68.2031 -68.2031]; %to the west and south of Rothera - Marg. Bay
%      LON=[-68.6920 -68.6920 -68.6920];
%      
%      LAT=[-68.34 -68.34 -68.34 -68.34 -68.34];  %point H used in the print outs
%      LON=[-69.01 -69.01 -69.01 -69.01 -69.01];
%      
%      
% %     LAT=[-68.1746 -68.1746 -68.1746 -68.1746 -68.1746];  %new point further north than H - approx latitude of I
% %     LON=[-69.0425 -69.0425 -69.0425 -69.0425 -69.0425];
% 
%      LAT=[-67.7 -67.7];  %newer point at the lefft side of the cross section throug the Fohn 
%      LON=[-70 -70];
%      
%       LAT=[-68.5573 -68.3340 -68.3327];  %even newer point at the lefft side of the new NE pointing cross section throug the Fohn 
%       LON=[-68.8593 -70.5348 -67.7867];
%       
%       LAT=[-68.5573 -68.5573 -68.5573 -68.5573 -68.5573 -68.5573 -68.5573 -68.5573];  %even newer point at the lefft side of the new NE pointing cross section throug the Fohn 
%       LON=[-68.8593 -68.8593 -68.8593 -68.8593 -68.8593 -68.8593 -68.8593 -68.8593];
%       
%       LAT=[-68.3186 -68.3186 -68.3186 -68.3186 -68.3186 -68.3186 -68.3186 -68.3186];  %even newer point at the lefft side of the new NE pointing cross section throug the Fohn 
%       LON=[-70.5566 -70.5566 -70.5566 -70.5566 -70.5566 -70.5566 -70.5566 -70.5566];
%       
%       LAT=[-68.3186 -68 -67 -70 -68 -70];
%       LON=[-70.5566 -75 -75 -83 -83 -86];
%       
%       LAT=[-68.7255];
%       LON=[-76.4423];


      LAT=[-67.55 -67.62 -67.55 -68.3186 -68.1 -68.7 -67.6 -66.2];  %4th Nov, 2010
      LON=[-68.1 -67.8 -67.5 -70.5566 -71.3 -76.4 -66.3 -68.5];  %points as seen for the upwind plots in Antarctica notes Dec08.doc
      
      LAT=[-68.1 -68.7 -66.2 -67.3];  %4th Nov, 2010
      LON=[-71.3 -76.4 -68.5 -72.9];  %points as seen for the upwind plots in Antarctica notes Dec08.doc
      
     
%      LAT=[-68.3186 -68.3186 -68.3186 -67];  %5th August, 2009
%      LON=[-70.5566 -70.5566 -70.5566 -75];
      
%      LAT=[-68.1 -68.1];
%      LON=[-70 -70];
     
     
%     LAT=[-68.1746 -68.1746 -68.1746 -68.1746 -68.1746];  %close to the peninsula for Froude analysis - looks like have low level blocking
%     LON=[-67.7088 -67.7088 -67.7088 -67.7088 -67.7088];  %location I
     
%     LAT=[-68.3131 -68.3131 -68.3131]; %to the west and south of Rothera - Marg. Bay
%     LON=[-69.2269 -69.2269 -69.2269];
     
%     LAT=[-68.2367 -68.2367 -68.2367]; %to the west and south of Rothera - Marg. Bay
%     LON=[-69.3400 -69.3400 -69.3400];
     
%     LAT=[-77.6479 -77.6479 -77.6479]; %over the continent at the base of the Peninsula (d02)
%     LON=[-93.4854 -93.4854 -93.4854];

%     LAT=[-68.34 -68.1746]; %locations H and I 
%     LON=[-69.01 -67.7088]; %to get labels right in plotTime...

%     LAT=[-66.8367 -66.9043 -67.1643 -66.3952 -66.8750]; %Avery Plateau 1st Dec 1995
%     LON=[-65.4891 -65.5770 -66.0512 -64.4811 -65.6243]; 
%     
%     LAT=[-66.6967];
%     LON=[-65.6505];
%     
%     LAT=[-66.8367]; %Avery Plateau 1st Dec 1995
%     LON=[-65.4891]; 

%        LAT=[-67.01]; %Larsen AWS
%        LON=[-61.55];
 
% ********************   Don't forget to set/unset these

    LAT=[];
    LON=[];

% *******************  

    
%    LAT=[-66.9043]; %nearby spot to Avery Plateau location where conentrations are larger
%    LON=[-65.5770];


%add the location of the aircraft ascent/descent if requested
if iplot_aircraft_locs_model==1  
    LAT = [LAT_extra LAT];
    LON = [LON_extra LON];
end

%now add the requested extra_x, extra_y locations for model plot
LAT=[LAT LAT_extra2];
LON=[LON LON_extra2];     
     
   
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%   points within the max flow   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 imax_flow=0;
 if imax_flow==1
 %find the location of the maximum wind speed on level 4 at around the lon of the aircraft max   
     lon_fli=-61.7; %latitude of the airdcraft at time of max speed
%     lon_fli=-63;
%     lon_fli=-60.85;
     dfli = 0.05;
     ilon_fli = find(lon2d.var>lon_fli-dfli & lon2d.var<lon_fli+dfli & lat2d.var>-68 & lat2d.var<-66);
     
     
     max_sp=0;

         u=nc{'U'}(time,4,:,:);
         v=nc{'V'}(time,4,:,:);
         u = 0.5 * ( u(:,1:end-1) + u(:,2:end) );
         v = 0.5 * ( v(1:end-1,:) + v(2:end,:) );
         
         sp=sqrt(u(ilon_fli).^2+v(ilon_fli).^2);
         [max_sp imax] = max(sp);
         
         [lat_max lon_max]=ind2sub(size(lon2d.var),ilon_fli(imax));
         LAT = [LAT lat2d.var(lat_max,lon_max)];
         LON = [LON lon2d.var(lat_max,lon_max)];

 end

    %1500m point on the descent is at -61.8, -67.15 (i.e. very similar to
    %the first aircraft trajectory point)
    
    %Nlocs determines how many locations to plot for each of the runs specified in rundir.dir
    %E.g. might want 2 locations for file 1 and two for file 2 - then put Nlocs=[2 1]
    %If only have one run then Nlocs is set to length(LAT)
    if sum(Nlocs)~=length(LAT) | length(LAT)==0
        Nlocs=[length(LAT)]; %switch off the extra user added locations
    end


    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);

    i=0;
    
%%%%%%%%% add aircraft profiles? Locations of ascents/descents described here %%%%%%%%%    
    if aircraft_comp == 1
        
        nmark(1)=0;

        [as_1 as_2]= findheight(time_flt19,19.332,19.53);  %first ascent over peninsula ('1')
        as_inds = [as_1:as_2];
        
        [as2_1 as2_2]= findheight(time_flt19,22.005,22.146); %final ascent back over peninsula  ('2')
        as2_inds = [as2_1:as2_2];

        [as3_1 as3_2]= findheight(time_flt19,20.744,20.78);
        [as4_1 as4_2]= findheight(time_flt19,21.176,21.201);
        [as5_1 as5_2]= findheight(time_flt19,21.552,21.577);
        asL_inds = [as3_1:as3_2 as4_1:as4_2 as5_1:as5_2];     %ascents during L-shaped sections ('4')
        
        [asL23_1 asL23_2]= findheight(time_flt19,20.76,21.55); %L-shaped sections of level flight ignoring when wind speed went wrong (too high)
        [asL23_3 asL23_4]= findheight(time_flt19,21.06,21.23);
        asL23_inds=[asL23_1:asL23_3 asL23_4:asL23_2]; %these two periods are then the 2nd and 3rd level flight sections

        

        [ds_1 ds_2]= findheight(time_flt19,20.24,20.4);  %first descent after peninsula ('0')
        ds_inds = [ds_1:ds_2];
        
        [asL_mini_descent_inds_1 asL_mini_descent_inds_2]= findheight(time_flt19,21.959,22.0052);
        asL_mini_descent_inds = [asL_mini_descent_inds_1:asL_mini_descent_inds_2];
        
        
        
        
        for iflight=1:length(ascent_str)
            i=i+1;
            

            
            if strfind(ascent_str(iflight),'1')
                inds = as_inds;
                labs(i).l = 'Aircraft ascent';
            elseif strfind(ascent_str(iflight),'3')
                inds = as2_inds;
                labs(i).l = 'Aircraft ascent2';
            elseif strfind(ascent_str(iflight),'4')
                inds = asL_inds;
                labs(i).l = 'Aircraft ascent-L';    
            elseif strfind(ascent_str(iflight),'0')
                inds = ds_inds;
                labs(i).l = 'Aircraft descent';
            elseif strfind(ascent_str(iflight),'5')
                inds = asL23_inds;
                labs(i).l = 'Aircraft L2 & L3'; 
            elseif strfind(ascent_str(iflight),'6')
                inds = asL_mini_descent_inds;
                labs(i).l = 'Aircraft mini descent';     
            end
            
            ydat(i).y = dat_flt19(inds,11) - 65;            
            disp('******  minus 65  ********');

            switch var
                case 'pressure'
              %      xlab='Pressure difference, WRF minus aircraft (hPa)';

              %      wrf_int = interp1(ydat(i).y,xdat(1).x,dat(ds_inds,11),[],'extrap');
              %      xdat(i).x = wrf_int - dat(ds_inds,4);
              %      ydat(i).y = dat(ds_inds,11);
                    
                    
                    xlab='Pressure (mb)';
                    
                    xdat(i).x = dat_flt19(inds,4);
                                                            



                    xlims=1;
                    xlimits=[975 990];
                    %		xlimits=[800 995];

                    izlim=1;
                    zmin=0;
                    zmax=120;
%                    zmax=3000;

                    %		logflag=1;

                    if iwrf_profs==0
                        figname=['Pressure profile '];
                    end

                case 'temperature'
                    xlab='Temperature (^{o}C)';
                    xdat(i).x = dat_flt19(inds,6); %use ds_inds for descent
                    


                    izlim=1;
                    zmin=0;
%                    zmax=3000;

                    xlims=1;
                    xlimits=[-15 6];
                    
                    if iwrf_profs==0
                        figname=['Temperature profile '];
                    end
                    
                case 'potemp'
                    xlab='Potential Temperature (K)';
                    T=dat_flt19(inds,6)+273.15;
                    P=dat_flt19(inds,4)*100;
                                        
                    xdat(i).x = T.*(1000e2./P).^0.286; %use ds_inds for descent
                    


                    izlim=1;
                    zmin=0;
                    zmax=3000;

                    xlims=0;
                    xlimits=[-15 -25];
                    
                    if iwrf_profs==0
                        figname=['Temperature profile '];
                    end    
                    
                    
                case 'equiv_potemp'
                    xlab='Equivalent potential Temperature (K)';
                    T=dat_flt19(inds,6)+273.15; %K
                    P=dat_flt19(inds,4)*100; %Pa                                            
                    qv = qv2_flt19(inds); %kg/kg
                    
                    xdat(i).x = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                                        
                   
                    


                    izlim=1;
                    zmin=0;
                    zmax=3000;

                    xlims=0;
                    xlimits=[-15 -25];
                    
                    if iwrf_profs==0
                        figname=['Temperature profile '];
                    end      
                    
                case 'vapour'
                    xlab='Vapour mixing ratio (g kg^{-1})';
                    xdat(i).x = 1000*qv2_flt19(inds);
                    

                    izlim=1;
                    zmin=0;
                    zmax=3000;

                    xlims=1;
                    xlimits=[0 5];
                    
                    if iwrf_profs==0
                        figname=['Vapour profile '];
                    end
                    
                case 'RH'
                    xlab='RH (%)';
                    T = dat_flt19(inds,6)+273.15;
                    P = dat_flt19(inds,4)*100;
                    qsat = SatVapPress(T,'goff','liq',P,1) / f;
                    xdat(i).x = qv2_flt19(inds)./qsat *100;
                    

                    izlim=1;
                    zmin=0;
                    zmax=3000;

                    xlims=0;
                    xlimits=[0 5];
                    
                    if iwrf_profs==0
                        figname=['RH profile '];
                    end    
                    
                case 'wind speed'
                    xlab='Wind speed (m s^{-1})';

                    xdat(i).x = dat_flt19(inds,9);
                    
                    

                    nfilter=20; bfilter=ones([1 nfilter])*1/nfilter;
                    xdat(i).x=filter(bfilter,1,xdat(i).x);
                    xdat(i).x(end-nfilter+1:end)=[]; ydat(i).y(end-nfilter+1:end)=[];
                    xdat(i).x(1:nfilter)=[]; ydat(i).y(1:nfilter)=[];    

                    izlim=1;
                    zmin=0;
                    zmax=3000;

                    xlims=1;
                    xlimits=[0 20];

                    if iwrf_profs==0
                        figname=['Wind speed profile '];
                    end
                    
                    
                    
                case 'wind dir'
                    xlab='Wind direction (degrees)';
                    
                    xdat(i).x = dat_flt19(inds,10)+180;
                    
                    
                    %wind aircraft data is very noisy so smooth over 20 points
                    nfilter=10; bfilter=ones([1 nfilter])*1/nfilter;                                     
                    xdat(i).x=filter(bfilter,1,xdat(i).x);
                    xdat(i).x(end-nfilter+1:end)=[]; ydat(i).y(end-nfilter+1:end)=[];
                    xdat(i).x(1:nfilter)=[]; ydat(i).y(1:nfilter)=[];    

                    izlim=1;
                    zmin=0;
                    zmax=3000;

                    xlims=1;
                    xlimits=[150 300];
                    
                    if iwrf_profs==0
                        figname=['Wind direction profile '];
                    end
            end

            [ydat(i).y,I] = sort(ydat(i).y);
            xdat(i).x = xdat(i).x(I);
            
        end
    
    end
    
    i_previous=i;
    nplots_previous=i;
    
    
    if iwrf_profs==1

      Nwrf_files=length(rundir);
    ipos_locs=1;
     for ifile_wrf=1:Nwrf_files
         
         if i_multi_wrf==1
                dire2 = [dire(ifile_wrf).dir rundir(ifile_wrf).dir];
                cd(dire2);
         end
            

            


        iloc_inds=ipos_locs:ipos_locs+Nlocs(ifile_wrf)-1;
        ipos_locs=ipos_locs+Nlocs(ifile_wrf);
        
        for iloc2=1:length(iloc_inds)    
            iloc=iloc_inds(iloc2);
            
            i=i+1;
            
            
            if i_multi_wrf==1 & ilabel_rundir==1
%                for ilab=1:length(rundir)
                    location_lab(iloc).l = rundir(ifile_wrf).dir;
                    location_lab(iloc).l = remove_character(location_lab(iloc).l,'_',' ');
%                end
            end


            
            %if are using more than one WRF file then cd to the correct directory
            
            if length(time_array)==1
                i_previous=i-1;
                time=time_array;

                if i_multi_wrf==1
                    if is_met_em(ifile_wrf)==1
                        Times=(nca(ifile_wrf).nc{'Times'}(:))';
                    else
                        Times=(nca(ifile_wrf).nc{'Times'}(:));
                    end
                else
                    if is_met_em(1)==1
                        Times=(nca(1).nc{'Times'}(:))';
                    else
                        Times=(nca(1).nc{'Times'}(:));
                    end
                end

            else
                time=time_array(iloc);
                if is_met_em(1)==1
                    Times=(nca(i-i_previous).nc{'Times'}(:))';
                else
                    Times=(nca(1).nc{'Times'}(:));
                end

            end

            


  
            
            

            if i_multi_wrf==1 & is_met_em(ifile_wrf)==1
                tstr=Times(:)';
            elseif i_multi_wrf==0 & is_met_em(1)==1
                tstr=Times(:)';
            else
                tstr=Times(time,:);
            end
                    

            iund=findstr('_',tstr);
            tstr(iund)=' ';  
            
            [year,month,day,hour,mins,sec,month_text]=WRF_time_strings(tstr);

            if (i_multi_wrf==1 & is_met_em(ifile_wrf)==1) | (i_multi_wrf==0 & is_met_em(1)==1)
                if prod(size(nc{'GHT'})) > 0
                    ydat(i).y = nca(ifile_wrf).nc{'GHT'}(1,:,ilat(iloc),ilon(iloc));
                    incep=1;
                else
                    incep=0;

                    pres=nca(ifile_wrf).nc{'PRES'}(1,:,ilat(iloc),ilon(iloc)); %BUT should the first pressure and temperature be counted - know
                    % that the pressure is the surface pressure - what does that mean for the temperature - skintemp?
                    pres=pres(2:end);   %don't want the first level - at least for ecmwf ml data
                    psfc=nca(ifile_wrf).nc{'PSFC'}(:); psfc=psfc(ilat(iloc),ilon(iloc));
                    pmsl=nca(ifile_wrf).nc{'PMSL'}(:); pmsl=pmsl(ilat(iloc),ilon(iloc));
                    temp=nca(ifile_wrf).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                    temp=temp(2:end);

                    %Height of the terrain as input from the analysis (underlying) model (coarse resolution interpolated onto wrf grid).
                    % (PSFC - this is the pressure at the height of the analysis terrain
                    % and not the pressure at the height of the high resolution terrain
                    %as taken from the geogrid program and given in HGT_M - tested this by checking the PSFC field in met_em
                    %files and wrfout and they are consistent with this).
                    %*** Actually, looking at the first level of PRES and comparing it to PSFC it seems the two are very close with both positive
                    %and negative differences. the PSFC lines look less smooth. PRES seems to follow the terrain more closely in terms of this
                    %smoothness anf the peaks and troughs correspond more closely - so this may be the pressure to choose to
                    %correspond to the height of SOILHGT
                    %NOTE that SOILHGT can be negative
                    %difference between choosing PSFC and PRES level one is about 80-100 m in the troposphere

                    %Here the hydrostatic equation is solved starting at the surface pressure, which is known to be at the height
                    %of the analysis terrain. The temperature vs pressure profile is described by PRES and TT and the required temperature
                    %at a given pressure (e.g. PSFC to start with) is found by interpolation in the hydrostatic2 function

                    %NOTE, might want to use this for finding pressure corrections for NCEP analysis files
                    %                     soilhgt = nc{'GHT'}(1,:,ilat(iloc),ilon(iloc)); %height of the second (=1000 hPa) pressure level
                    %                     soilhgt = soilhgt(2);
                    %                     [Y,I]=sort(pres,'descend'); %sort the data in order to take into account pres(2)(=surface) being lower than first level of 1000 hPa
                    %                     pres=pres(I);
                    %                     temp=temp(I);


                    if prod(size(nca(ifile_wrf).nc{'SOILHGT'})) > 0
                        soilhgt=nca(ifile_wrf).nc{'SOILHGT'}(:); soilhgt=soilhgt(ilat(iloc),ilon(iloc));  %height of PSFC level
                        pres2=pres;
                        temp2=temp;

                    else
                        pres2 = [pmsl pres];
                        skinT = temp(1)*(pmsl/pres(1))^0.286; %temperature assuming a constant potential temperature between msl and the surface point
                        %skinT = nc{'SKINTEMP'}(:); %using skin temperature as sea level temp - could be inaccurate?
                        %skinT = skinT(ilat(iloc),ilon(iloc));
                        %    skinT=skinT+5; %to test temperature sensitivity - does not seem that senstitive - 5 degree increase led to only 0.14 hPa increase in
                        %in estimated pressure of a height below soilhgt (tried 91 m wheras soilhgt was 111m - difference was even less for lower altitudes).
                        temp2 = [skinT temp];
                        %                soilhgt=nc{'SOILHGT'}(:); soilhgt=soilhgt(ilat(iloc),ilon(iloc));  %height of PSFC level
                        soilhgt=0;
                    end




                    PSPAN=[pres2(1) pres2(end)]; %pressure range for integration
                    [P,hp] = ODE45(@hydrostatic2,PSPAN,soilhgt,[],pres2,temp2); %enter the initial height for the given the first value of PSPAN
                    %note that some of the pressure levels in the PRES array will be below the

                    ydat(i).y = interp1(P,hp,pres);

                    %                     pres=nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    %                     ydat(i).y = find_height_from_p_ant_d03(pres);

                end
            else
                ydat(i).y = WRFUserARW(nca(1).nc,'Z',time,ilat(iloc),ilon(iloc));
            end

            %	ydat(i).y = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));
            switch var
                case 'pressure'
                    figname=['Pressure profile at ' tstr ' for ' filestr];
                    
                    if is_met_em(1)
                        xdat(i).x = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    else                        
                        xdat(i).x = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));
                    end
                   % iydir = -1; %reverse the direction of the pressure axis so is right way around
                   
               case 'equiv_potemp'
                    figname=['Equivalent potential temperature profile at ' tstr ' for ' filestr];
                    
                    if is_met_em(1)
%                        xdat(i).x = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    else        
                        potemp = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                        P = nc{'P'}(time,:,ilat(iloc),ilon(iloc)) + nc{'PB'}(time,:,ilat(iloc),ilon(iloc));
                        T = potemp ./ ( (1e5./P).^0.286 );
                        qv = nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                        xdat(i).x = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';

                    end
                    
                    xlims=0;
                    xlimits=[275 290];
                    
                   
                   case 'potemp'
                    figname=['Potential temperature profile at ' tstr ' for ' filestr];
                    xlab='Potential temperature (K)';
                    if is_met_em(1)
%                        xdat(i).x = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    else        
                        xdat(i).x = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;                        
                    end
                    
                    xlims=1;
                    xlimits=[270 295];
%                    xlimits=[270 305];
%                    xlimits=[270 325];
                    
                case 'temperature'

                    figname=['Temperature profile at ' tstr ' for ' filestr];
                    xlab='Temperature (^{o}C)';

                    if is_met_em(ifile_wrf)==1
                        xdat(i).x = nca(ifile_wrf).nc{'TT'}(1,:,ilat(iloc),ilon(iloc)) - 273.15;
                        if incep==0
                            xdat(i).x = xdat(i).x(2:end);
                        end
                    else
                        xdat(i).x = WRFUserARW(nca(ifile_wrf).nc,'tc',time,ilat(iloc),ilon(iloc));
                        xdat(i).x = [get_wrf_point_surface(nca(ifile_wrf).nc,'T2',time,ilat(iloc),ilon(iloc))-273.15 xdat(i).x];
                        terr_level = nca(ifile_wrf).nc{'HGT'}(:,ilat(iloc),ilon(iloc));
                        terr_level=terr_level(time);
                        ydat(i).y = [terr_level+2 ydat(i).y]; %add air temp at 2 m
                    end
                    
                    xlims=1;
                    xlimits=[-15 6];

                case 'vapour'
                    figname=['Vapour profile at ' tstr ' for ' filestr];
                    xlab='Vapour mixing ratio (g kg^{-1})';
                    
                    if is_met_em(1)
                        rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                        T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                        P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                        qsat = satvappress(T,'goff','liq',P,1)/f;
                        xdat(i).x = 1000 * rh/100 .* qsat;
                    else
                        xdat(i).x = 1000*nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                    end
                    
                 case 'ice'

                    xlab='Ice number concentration (L^{-1})';
%                    xlab='Snow number concentration (L^{-1})';
%                    xlab='Graupel number concentration (L^{-1})';
%                    xlab='Rain number concentration (L^{-1})';                    
%                    xlab='Cloud mixing ratio (g kg^{-1})';
%                    xlab='Ice mixing ratio (g kg^{-1})';
%                    xlab='Snow mixing ratio (g kg^{-1})';                    
%                    xlab='Graupel mixing ratio (g kg^{-1})'; 
%                    xlab='Total condensate mixing ratio (g kg^{-1})'; 
%                    xlab='Total number concentration (L^{-1})'; 
%                    xlab='Water supersaturation (%)'; 
%                    xlab='Ice supersaturation (%)';                     
                    figname=[xlab '  C ' tstr ' for ' filestr];
                                       
                    if is_met_em(1)
                        rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                        T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                        P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                        qsat = satvappress(T,'goff','liq',P,1)/f;
                        xdat(i).x = 1000 * rh/100 .* qsat;
                    else
                        
                            T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                            P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;
                            rho=density(P,T);
                            
                            %numbers from WRF are in #/kg so multiply by the density to get #/m3 and then divide by 1000 to get #/L
                            
                        switch xlab
                            case 'Ice number concentration (L^{-1})'                        
                                xdat(i).x = nc{'QNICE'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                            case 'Snow number concentration (L^{-1})'    
                                xdat(i).x = nc{'QNSNOW'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;   
                            case 'Graupel number concentration (L^{-1})'    
                                xdat(i).x = nc{'QNGRAUPEL'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000; 
                            case 'Rain number concentration (L^{-1})'
                                xdat(i).x = nc{'QNRAIN'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000; 
                            case 'Ice mixing ratio (g kg^{-1})'
                                xdat(i).x = 1000*nc{'QICE'}(time,:,ilat(iloc),ilon(iloc)); 
                            case 'Snow mixing ratio (g kg^{-1})'    
                                xdat(i).x = 1000*nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));   
                            case 'Graupel mixing ratio (g kg^{-1})'                                   
                                xdat(i).x = 1000*nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                            case 'Cloud mixing ratio (g kg^{-1})'    
                                xdat(i).x = 1000*nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));  
                            case 'Total condensate mixing ratio (g kg^{-1})'
                                xdat(i).x = 1000*( nc{'QICE'}(time,:,ilat(iloc),ilon(iloc))+nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc))...
                                    +nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc))+nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc))...
                                    +nc{'QRAIN'}(time,:,ilat(iloc),ilon(iloc)) );
                            case 'Total number concentration (L^{-1})'                                    
                                xdat(i).x = ( nc{'QNICE'}(time,:,ilat(iloc),ilon(iloc))+nc{'QNSNOW'}(time,:,ilat(iloc),ilon(iloc))...
                                    +nc{'QNGRAUPEL'}(time,:,ilat(iloc),ilon(iloc))+nc{'QNRAIN'}(time,:,ilat(iloc),ilon(iloc)) ).*rho/1000;
                            case 'Water supersaturation (%)'
                                qv=nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                qvs=satvappress(T,'goff','liq',P,1)/f;
                                xdat(i).x=100*(qv./qvs-1);
                            case 'Ice supersaturation (%)'
                                qv=nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                qvs=satvappress(T,'goff','ice',P,1)/f;
                                xdat(i).x=100*(qv./qvs-1);                                
                                
                        end





                    end   
                    
                    xlims=0;
                    xlimits=[-20 1];
                    
                 case 'RH'
                    figname=['RH profile at ' tstr ' for ' filestr];
                    xlab='Relative humidity (%)';
                    
                    if is_met_em(1)
                        xdat(i).x = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                    else                        
                        qv = nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                        T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                        P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;
                        qsat = satvappress(T,'goff','liq',P,1)/f;
                        xdat(i).x = qv./qsat *100;
                    end   

                case 'wind speed'
                     
                    if is_met_em(ifile_wrf)                        
                        u=0.5*(nca(ifile_wrf).nc{'UU'}(1,:,ilat(iloc),ilon(iloc))+nca(ifile_wrf).nc{'UU'}(1,:,ilat(iloc),ilon(iloc)+1));
                        v=0.5*(nca(ifile_wrf).nc{'VV'}(1,:,ilat(iloc),ilon(iloc))+nca(ifile_wrf).nc{'VV'}(1,:,ilat(iloc)+1,ilon(iloc)+1));                        
%                        v=nca(ifile_wrf).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else
                        u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                        v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                    end

                    xdat(i).x= sqrt( u.^2 + v.^2 );

                    figname=['Wind speed profile at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Wind speed profiles'];
                    end
                    xlab='Wind speed (m s^{-1})';
                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end
                    
                case 'wind speed component'                   
                   
                   xlims=0;
                   xlimits=[0 2.5];

                   
                    if is_met_em(1)
                        u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                        v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else
                        if i_multi_wrf==1
                            dire2 = [dire(iloc).dir rundir(iloc).dir];
                            cd(dire2);
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                        else
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                        end
                    end

                    sp = sqrt( u.^2 + v.^2 );
                    
                    
            %%%%%%%%%%%% get wind dir so can calculate easterly component
                    jnorth = ilat(iloc) + 10;
                    lons_north = lon2d.var(jnorth,:);
                    [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );
                    
                    %angle of the local north line relative to the grid
                    thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) ); 
                    
                    clear dir
                    for iuv=1:length(u)

                        theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                        if u(iuv)==0 & v(iuv)==0
                            dir(iuv) = 0;
                        elseif u(iuv)>=0 & v(iuv)>=0
                            dir(iuv) = theta2;
                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<=0 & v(iuv)<=0
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<0 & v(iuv)>0
                            dir(iuv) = 360 + theta2; %theta2 is negative
                        end




                    end
                                       
                    dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north

                    
                    %%%%% easterly component 
                    xdat(i).x = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)                                                            


                    figname=['Westerly wind speed profile at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Westerly wind speed profiles'];
                    end
                    xlab='Westerly wind speed'; %dimensionless
                                                           
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end         
                    

                case 'wind dir'
                    xlab='Wind direction (degrees)';
                    figname=['Wind dir profile at ' tstr ' for ' filestr];
                    
                    if is_met_em(ifile_wrf)
                        u=nca(ifile_wrf).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                        v=nca(ifile_wrf).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else                       
                        u=WRFUserARW(nca(ifile_wrf).nc,'u',time,ilat(iloc),ilon(iloc));
                        v=WRFUserARW(nca(ifile_wrf).nc,'v',time,ilat(iloc),ilon(iloc));
                    end                                        
                    
                    xdat(i).x=wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat(iloc),ilon(iloc),DX,DY);
                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end
                    
                    xlims=1;
                    xlimits=[130 325];
                    
                case 'cloud'
                    cloud=nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                    xdat(i).x=1000*cloud;


                    figname=['Cloud mixing ratio profile at ' tstr ' for ' filestr];
                    xlab='Total condensed water mixing ratio (g kg^{-1})';
                    
               case 'Froude'                   
                   H=1.5e3; %set height of mountain
%                   H=ydat(i).y(1);  %or use the height of the topograpy at the location in question
                   
                   dH=500; %layer depth for the sub/super critical type of Froude number above the mountain
                            %Tom suggested at least 1000m
                   
                   xlims=0;
                   xlimits=[0 2.5];
                   
                   %%% for Froude number need the wind speed and direction (to calc the wind speed
                   %%% perpendicular to the mountain) and the potemp
                   
                    if is_met_em(1)
                        u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                        v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else
                        if i_multi_wrf==1
                            dire2 = [dire(iloc).dir rundir(iloc).dir];
                            cd(dire2);
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                        else
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                        end
                    end

                    sp = sqrt( u.^2 + v.^2 );
                    
                    
            %%%%%%%%%%%% get wind dir so can calculate easterly component
                    jnorth = ilat(iloc) + 10;
                    lons_north = lon2d.var(jnorth,:);
                    [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );
                    
                    %angle of the local north line relative to the grid
                    thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) ); 
                    
                    clear dir
                    for iuv=1:length(u)

                        theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                        if u(iuv)==0 & v(iuv)==0
                            dir(iuv) = 0;
                        elseif u(iuv)>=0 & v(iuv)>=0
                            dir(iuv) = theta2;
                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<=0 & v(iuv)<=0
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<0 & v(iuv)>0
                            dir(iuv) = 360 + theta2; %theta2 is negative
                        end




                    end
                                       
                    dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north
                    
                    i_wind_component=1;
                    if i_wind_component==1
                        %%%%% component in THI direction component
                        THI=70*pi/180; %e.g. 45 = southwesterly
                        spE = sp.*cos(dir-THI); %get the wind direcion in the THI component 
                        info = 'using wind component';
                    else
                        spE = sp; %OR ignore direction and just use the magnitude
                        info='';
                    end
                    
                    %%%% potemp
                    if is_met_em(1)
%                        pot = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    else        
                        pot = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;                        
                    end
                    
       %%%%%%%%calculate Froude = 1/h_hat = U / (H*sqrt(g/pot * dpot/dz) )
%                    xdat(i).x = spE(2:end)./ (H.* sqrt(9.81./pot(2:end) .*(diff(pot)./diff(ydat(i).y)) ) );

                    Fr_type='below';  %using (H-z) but using the local N value
                    Fr_type='below_mean';  %using (H-z) and the average N and U value for the heights from z to H                   
%                    Fr_type='below_layers'; 
%                    Fr_type='one_layer';                     
%                    Fr_type='sub/super';
%                    Fr_type='sub/super one layer';                    
%                    Fr_type='sub/super Holton';
                    switch Fr_type
                        case 'below'
                            xdat(i).x = spE(2:end)./ ((H-ydat(i).y(2:end)).* sqrt(9.81./pot(2:end) .*(diff(pot)./diff(ydat(i).y)) ) );
                            if i==1
                                H_froude=ydat(1).y;
                            end
                            ydat(i).y = ydat(i).y(2:end);
                            % Froude of 0.5 corresponds to h_hat (effective barrier height) of 2
                            % so Fr of less than 0.5 suggests that flow will be blocked (van Lipzig paper, 2008)

                            ibel_0 = find(xdat(i).x<0);
                            xdat(i).x(ibel_0)=[];
                            ydat(i).y(ibel_0)=[];
                            
                         case 'below_mean'
                            iH = findheight(ydat(i).y,H) ;
                            for ih=2:iH
                                U = mean(spE(ih:iH));
                                U = spE(iH);
                                xdat(i).x(ih) = U./ (  (H-ydat(i).y(ih)).* sqrt( 9.81./mean(pot(ih:iH)) .*(diff(pot([ih iH]))./diff(ydat(i).y([ih iH]))) )  );
                            end

                            if i==1
                                H_froude=ydat(1).y;
                            end
                            ydat(i).y = ydat(i).y(2:iH);                            

                            ibel_0 = find(xdat(i).x<0);
                            xdat(i).x(ibel_0)=[];
                            ydat(i).y(ibel_0)=[];   
                            
                        case 'one_layer'                            
                                iH = findheight(ydat(i).y,H) ;                                                        
                                U = spE(iH);
%                                pot_layer = mean( pot(1:iH) );   
                                pot_layer =pot(iH);                                   
                                dpot_dz = (pot(iH) - pot(1))  /(ydat(i).y(iH) - ydat(i).y(1));
                                
                            xdat(i).x = U/ ( H* sqrt(9.81./pot_layer*dpot_dz) );                                                           
%                            xdat(i).x = 2*pi*60e3*U/ ( 60e3* sqrt(9.81./pot_layer*dpot_dz) );                                                           
                            ydat(i).y = H;

                            
                        case 'below_layers'
                            %doing dH m layers for this starting at lowest level                            
                            dH=600;
                            layers = [ydat(i).y(1) ydat(i).y(1)+dH:H-dH H];
                            [ilayers,I] = unique( findheight(ydat(i).y,layers) );                                                        
                            layers=layers(I);
                            
                            clear U pot_layer z_layer
                            for ifr=2:length(ilayers)
                                U(ifr-1) = mean(  spE( ilayers(ifr-1) : ilayers(ifr) )  );
                                pot_layer(ifr-1) = mean(  pot( ilayers(ifr-1) : ilayers(ifr) )  );
                                z_layer(ifr-1) = (ydat(i).y(ilayers(ifr-1)) + ydat(i).y(ilayers(ifr)) )/2; %centre of layer
                            end                            

                            
                            dpot_dz = diff(pot_layer)./diff(z_layer);
                            xdat(i).x = - U(2:end)./ ( (z_layer(2:end)-H).* sqrt(9.81./pot_layer(2:end) .*dpot_dz) );                                                           
                            ydat(i).y = z_layer(2:end);

                        case 'sub/super'
                            %need 1 km layers for this starting at mountain top
                            %for Fr<1 have subcritical flow where modes may propagate against the flow
                            %for Fr>1 flow is supercritical and no modes may propagate upstream (see Orr,2008)
                            layers2 = [H:dH:ydat(i).y(end)];
                            ilayers = [findheight(ydat(i).y,layers)];
                            layers = ydat(i).y(ilayers);
                            
                            for ifr=2:length(ilayers)
                                U = mean(  spE( ilayers(ifr-1) : ilayers(ifr) )  );
                                pot_layer = mean(  pot( ilayers(ifr-1) : ilayers(ifr) )  );                                
%                                dpot_dz = (pot(ilayers(ifr)) - pot(ilayers(ifr-1))  ) ./ (layers(ifr) - layers(ifr-1));
                                dpot_dz = ( pot(ilayers(ifr)) - pot(ilayers(ifr-1)) ) ./ (layers(ifr) - layers(ifr-1));
                                xdat(i).x(ifr-1) = U./ ( (layers(ifr)-H).* sqrt(9.81./pot_layer .*dpot_dz) );                               
                            end                            
                            % super/sub critical test for above the mountain
                             ydat(i).y = layers(2:end);
                             
                          case 'sub/super one layer'
                            %need 1 km layers for this starting at mountain top
                            %for Fr<1 have subcritical flow where modes may propagate against the flow
                            %for Fr>1 flow is supercritical and no modes may propagate upstream (see Orr,2008)
                            iH=findheight_nearest(ydat(i).y,H)-1;
                            layer_top = ydat(i).y(iH)+dH; %layer top is just above mountain plus dH
                            ilayer = [findheight_nearest(ydat(i).y,layer_top)];
                            if  ydat(i).y(ilayer)-ydat(i).y(iH) < dH 
                                ilayer=ilayer+1;
                            end
                            
                            %get pressure in Pa
                            if is_met_em(1)
                                P = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc));
                            else
                                P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc))*100;                                
                            end 
                            %get temperature in K
                            if is_met_em(1)
                                T = nca(i-i_previous).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));                                
                            else
                                T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc))+273.15;                                
                            end
                            
                            rho = density(P,T);
                                                        
                                U = mean(  spE( iH:ilayer ) );
%                                pot_layer = mean(  pot( iH:ilayer )  );                                
                                pot_layer = pot(ilayer);    
 %                               pot_layer = rho(ilayer);           
                                
                                layer_depth = ydat(i).y(ilayer)-ydat(i).y(iH);
                                
                                dpot_dz = (pot(ilayer) - pot(iH)  ) / layer_depth;                                
%                                dpot_dz = -(rho(ilayer) - rho(iH)  ) / layer_depth;                                
                                
                                N = sqrt(9.81/pot_layer * dpot_dz);

                                             
                            % super/sub critical test for above the mountain
                             ydat(i).y = ydat(i).y(ilayer);   
                             xdat(i).x = U/ ( (ydat(i).y-H) * N );                               
                             
                         case 'sub/super Holton';    
                             %need 1 km layers for this starting at mountain top
                            %for Fr<1 have subcritical flow where modes may propagate against the flow
                            %for Fr>1 flow is supercritical and no modes may propagate upstream (see Orr,2008)
                            
                            if dH==0 %use all layers from WRF
                                layers = [H ydat(i).y];
                                ilayers = [1:length(ilayers)];
                            else %use layers of depth dH
                                layers = [H:dH:ydat(i).y(end)];
                                ilayers = findheight(ydat(i).y,layers);
                            
                            end
                            
                            %get pressure in Pa
                            if is_met_em(1)
                                P = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc));
                            else
                                P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc))*100;                                
                            end 
                            %get temperature in K
                            if is_met_em(1)
                                T = nca(i-i_previous).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));                                
                            else
                                T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc))+273.15;                                
                            end
                            
                            rho = density(P,T);
                            
                            for ifr=2:length(ilayers)-1
                                U = mean(  spE( ilayers(ifr) : ilayers(ifr+1) )  );
                                rho_layer = mean(  rho( ilayers(ifr) : ilayers(ifr+1) ) );
                                drho = mean(  rho( ilayers(ifr-1) : ilayers(ifr) )  ) - rho_layer;   %%?? is this the right drho??                            
                                xdat(i).x(ifr-1) = U.^2 ./ ( 9.81*(layers(ifr)-H) * drho/rho_layer ); %Fr=g*H*drho/rho1 from Holton                               
                                %think H is the (mean) height above the mountain of the layer
                            end                            
                            % super/sub critical test for above the mountain
                             ydat(i).y = layers(2:end-1);
                    end

                    figname=['Froude number ' info ' profile at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Froude number ' info ' profiles'];
                    end
                    xlab='Froude number'; %dimensionless
                    
                    

                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end     
                    
            case 'density gradient'                                     
                
                xlims=0;
                xlimits=[0 2.5];

                %get pressure in Pa
                if is_met_em(1)
                    P = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc));
                else
                    P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc))*100;
                end
                %get temperature in K
                if is_met_em(1)
                    T = nca(i-i_previous).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                else
                    T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc))+273.15;
                end

                %density gradient
                xdat(i).x = diff(density(P,T))./diff(ydat(i).y);
                ydat(i).y = ydat(i).y(2:end);
                
                %density
%                xdat(i).x = density(P,T);




                figname=['Density gradient profile at ' tstr ' for ' filestr];
                if i_paper_labels==1
                    figname=['Density gradient profiles'];
                end
                xlab='Density gradient (kg m^{-4})'; 


                if incep==0 & is_met_em(1)==1
                    xdat(i).x = xdat(i).x(2:end);
                end
                    
            case 'dpot/dz'
                   
                    xlims=0;
                    xlimits=1000*[0 0.025];
                    
                    %%%% potemp
                    if is_met_em(1)
%                        pot = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    else        
                        pot = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;                        
                    end
                    
                    
                    xdat(i).x = 1000*(diff(pot)./diff(ydat(i).y));
                    ydat(i).y = ydat(i).y(2:end);
                    

                    

                    figname=['Potential temperature gradient profile at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Potential temperature gradient profiles'];
                    end
                    xlab='Potential temperature gradient (K km^{-1})'; %
                    
                    

                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end    
                    
            case 'westerly wind'                  
                   
                    if is_met_em(1)
                        u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                        v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else
                        if i_multi_wrf==1
                            dire2 = [dire(iloc).dir rundir(iloc).dir];
                            cd(dire2);
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                        else
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                        end
                    end

                    sp = sqrt( u.^2 + v.^2 );
                    
                    
            %%%%%%%%%%%% get wind dir so can calculate easterly component
                    jnorth = ilat(iloc) + 10;
                    lons_north = lon2d.var(jnorth,:);
                    [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );
                    
                    %angle of the local north line relative to the grid
                    thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) ); 
                    
                    clear dir
                    for iuv=1:length(u)

                        theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                        if u(iuv)==0 & v(iuv)==0
                            dir(iuv) = 0;
                        elseif u(iuv)>=0 & v(iuv)>=0
                            dir(iuv) = theta2;
                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<=0 & v(iuv)<=0
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<0 & v(iuv)>0
                            dir(iuv) = 360 + theta2; %theta2 is negative
                        end




                    end
                                       
                    dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north

                    
                    %%%%% easterly component 
                    xdat(i).x = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)                                                            
                    
                            
            
                                        

                    figname=['Westerly wind component profile at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Westerly wind component profiles'];
                    end
                    xlab='Westerly wind component (m s^{-1})'; %
                    
                    

                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end   
                    
                case 'stratification integral'
                   
                    xlims=0;
                    xlimits=[0 0.025];
                    
                %get pressure in Pa
                if is_met_em(1)
                    P = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc));
                else
                    P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc))*100;
                end
                %get temperature in K
                if is_met_em(1)
                    T = nca(i-i_previous).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                else
                    T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc))+273.15;
                end
                
                rho=density(P,T);
                

               if is_met_em(1)
                        u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                        v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else
                        if i_multi_wrf==1
                            dire2 = [dire(iloc).dir rundir(iloc).dir];
                            cd(dire2);
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                        else
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                        end
                    end

                    sp = sqrt( u.^2 + v.^2 );
                    
                    
            %%%%%%%%%%%% get wind dir so can calculate easterly component
                    jnorth = ilat(iloc) + 10;
                    lons_north = lon2d.var(jnorth,:);
                    [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );
                    
                    %angle of the local north line relative to the grid
                    thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) ); 
                    
                    clear dir
                    for iuv=1:length(u)

                        theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                        if u(iuv)==0 & v(iuv)==0
                            dir(iuv) = 0;
                        elseif u(iuv)>=0 & v(iuv)>=0
                            dir(iuv) = theta2;
                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<=0 & v(iuv)<=0
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<0 & v(iuv)>0
                            dir(iuv) = 360 + theta2; %theta2 is negative
                        end




                    end
                                       
                    dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north

                    
                    %%%%% easterly component 
                   spE = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)
                    
                    iH=findheight(ydat(i).y,H);
                    rho2=mean(rho(1:iH));
                    rho2=rho(iH);
                    U = spE(iH);                    
                    
                    for ih=2:iH
                        xdat(i).x(ih-1)=9.81*sum( (H-ydat(i).y(ih:iH)).*-diff(rho(ih-1:iH)) ) - 0.5*rho2*spE(ih)^2;
%                        xdat(i).x(ih-1)=9.81*sum( (H-ydat(i).y(ih:iH)).*-diff(rho(ih-1:iH)) ) - 0.5*rho(1)*10^2;
                    end
                    
                    ydat(i).y=ydat(i).y(2:iH);

                    

                    figname=['Integral profile at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Integral profiles'];
                    end
                    xlab='Integral'; %dimensionless
                    
                    

                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end    
                    
                    
            case 'Scorer parameter'  
                
                 if is_met_em(1)
%                        pot = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                else        
                        pot = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;                        
                 end
                    
                   
                    if is_met_em(1)
                        u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                        v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                    else
                        if i_multi_wrf==1
                            dire2 = [dire(iloc).dir rundir(iloc).dir];
                            cd(dire2);
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                        else
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                        end
                    end

                    sp = sqrt( u.^2 + v.^2 );
                    
                    
            %%%%%%%%%%%% get wind dir so can calculate component
                    jnorth = ilat(iloc) + 10;
                    lons_north = lon2d.var(jnorth,:);
                    [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );
                    
                    %angle of the local north line relative to the grid
                    thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) ); 
                    
                    clear dir
                    for iuv=1:length(u)

                        theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                        if u(iuv)==0 & v(iuv)==0
                            dir(iuv) = 0;
                        elseif u(iuv)>=0 & v(iuv)>=0
                            dir(iuv) = theta2;
                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<=0 & v(iuv)<=0
                            dir(iuv) = 180 + theta2;
                        elseif u(iuv)<0 & v(iuv)>0
                            dir(iuv) = 360 + theta2; %theta2 is negative
                        end




                    end
                                       
                    dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north

                    
                    %%%%% component in certain direction
                    
                    dir_comp=pi/2; %the direction of the component we want taken from north - so for the westerly component set = pi/2
%                    dir_comp=0; %for northerly set dir_comp=0
                    U = sp.*cos(dir_comp-dir); %get the wind direcion in the dir_comp component                                                                                
                    U=sp;
                    
                    
                    Z=ydat(i).y;
                    
                    %z_fine=[Z(1):10:Z(end)];
                    %U_fine=interp1(Z,U,z_fine);
                    %pot_fine=interp1(Z,pot,z_fine);
                    z_fine=Z;
                    U_fine=U;
                    pot_fine=pot;
                    
                    
                    N2 = 9.81 ./ pot_fine(2:end) .* diff(pot_fine)./diff(z_fine) ; %N squared 
                    dUdz = diff(U_fine)./diff(z_fine);
                    dU2dz2 = diff(dUdz)./diff(z_fine(2:end));
                    
                    scorer_term='both';                        
%                    scorer_term='second';
%                    scorer_term='first';
%                    scorer_term = 'wind_shear'
%                    scorer_term='gradient';   
%                    scorer_term='Richardson'; 
                    scorer_term='N'; 
%                    scorer_term='N/U'; 

                    switch scorer_term
                        case 'first'
%%%% just the first term                                        
                    xdat(i).x = N2(1:end)./U_fine(2:end).^2;
                    ydat(i).y = z_fine(2:end);
                    scorer_label='First term of Scorer parameter';
                    units='(m^{-2})';
                    
                        case 'second'
%%%% just the second term                    
                    xdat(i).x =  - dU2dz2./U_fine(2:end-1) ;                    
                    ydat(i).y = z_fine(2:end-1);
                    scorer_label='Second term of Scorer parameter';
                    units='(m^{-2})';
                    
                        case 'both'                    
%%%% both terms                    
                    xdat(i).x = N2(1:end-1)./U_fine(2:end-1).^2 - dU2dz2./U_fine(2:end-1) ;                    
                    ydat(i).y = z_fine(2:end-1);
                    scorer_label='Scorer parameter';
                    units='(m^{-2})';
                                        
                    case 'wind_shear'                    
%%%% both terms                    
                    xdat(i).x = dUdz;                    
                    ydat(i).y = z_fine(2:end);
                    scorer_label='Wind shear';
                    units='(s^{-1})';
                    
                    case 'gradient'                    
%%%% both terms                    
                     scorer = N2(1:end-1)./U_fine(2:end-1).^2 - dU2dz2./U_fine(2:end-1) ;        
                     xdat(i).x = diff(scorer) ./ diff(z_fine(2:end-1));
                    ydat(i).y = z_fine(3:end-1);
                    scorer_label='Scorer parameter';
                    units='(m^{-3})';
                    
                    case 'Richardson'                    
%%%% both terms                    

                    dUdz= diff(u)./diff(z_fine);
                    dVdz= diff(v)./diff(z_fine);
                    
                    xdat(i).x = N2./( (dUdz).^2 + (dVdz).^2);                       
                    ydat(i).y = z_fine(2:end);
                    scorer_label='Richardson number';
                    units=''; %no units
                    
                    xlims=1;
                    xlimits=[0 2];
                    
                    
                    case 'N'                                 
                  
                    xdat(i).x = sqrt(N2);
                    ydat(i).y = z_fine(2:end);
                    scorer_label='N';
                    units='s^{-1}'; %no units
                    
                    xlims=0;
                    xlimits=[0 2];
                    
                    case 'N/U'                                 
                  
                    xdat(i).x = sqrt(N2)./U_fine(2:end);
                    ydat(i).y = z_fine(2:end);
                    scorer_label='N/U';
                    units='m^{-1}'; %no units
                    
                    xlims=0;
                    xlimits=[0 2];
                    
                    end
                    
            
                                        

                    figname=[scorer_label ' at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=[scorer_label ' profiles'];
                    end
                    xlab=[scorer_label ' ' units]; %
                                                           
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end    
                    
                     
                    
            
                                        

                    figname=[scorer_label ' at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=[scorer_label ' profiles'];
                    end
                    xlab=[scorer_label ' ' units]; %
                                                           
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end    
                    
                    
                case 'ridge_height'
                   %%% run heat_fluxes_mean_along_latitude.m to get the max ridge height
                    xlims=1;
                    xlimits=[-69 -67];                    
                    
                    xdat(i).x = LAT_ridge;
                    ydat(i).y = peak_vs_lat;
                                       
                    figname=['Max ridge height at ' tstr ' for ' filestr];
                    if i_paper_labels==1
                        figname=['Max ridge height'];
                    end
                    xlab='Latitude'; %
                    
                    

                    
                    if incep==0 & is_met_em(1)==1
                        xdat(i).x = xdat(i).x(2:end);
                    end     
                    
                    no_sort=1;
                    
                case 'Houghton Smith formulae'
                   
                    xlims=0;
                    xlimits=1000*[0 0.025];
                                      
                    gd=1;
                    U=20;
                    N=0.02;

                    hhat=[0:0.01:1];

                    %continuous stratification case (Fig. 15b)
                    del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
                    H0_hat_strat=hhat - del_hat + acos(hhat./del_hat);
                    H0_crit=U/N*H0_hat_strat; %
                    F0_crit=1 ./ H0_hat_strat; %=U/NH=F0=1/H0_hat for Fig 15b


                    %single layer case (Fig. 15a)

                    H0_crit2 = [0:1:7000];  %N.B. - require very large H0 values to reach high h/H0 values (to approach the limit
                    %towards h/H0=1,F0=0
                    F0_crit2=U./sqrt(gd.*H0_crit2);
                    F0_crit3 = U/N./H0_crit2;
                    h_crit2 = H0_crit2 .* ( 1 + 0.5*F0_crit2.^2 - 1.5*F0_crit2.^(2/3) );
                    h_crit3 = H0_crit2 .* ( 1 + 0.5*F0_crit3.^2 - 1.5*F0_crit3.^(2/3) );
                    hhat_Houghton=h_crit3*N/U;
                    
                    plot_case='stratified';
                    plot_case='h/H0,F0 space';
                    plot_case='h vs H0';                    
                    
                    switch plot_case
                        case 'stratified'

                            xdat(1).x=hhat_Houghton;
                            ydat(1).y=F0_crit3; %using Smith formula for F0 but with F0=U/(N*H0) - acheived by approximating
                            %g' with N^2*H0 where H0=dz in N^2=g/theta * dtheta/dz
                            labs(1).l='Smith using N';

                            xdat(2).x=hhat;
                            ydat(2).y=F0_crit; %plot the original formulation for continuous stratification
                            % are different but.... if multiply F0 by a factor of 0.62...
                            labs(2).l='Original Smith';

                            xdat(3).x=hhat_Houghton;
                            ydat(3).y=F0_crit3*0.62; %they overlay each other almost exactly. Why 0.62??? - epsilon? =0.622
                            %is independant of N and U
                            labs(3).l='Using N and *0.62';
                            
                            xlab='F0=U/(N*H0)';
                            ylab='hN/U';


                        case 'h/H0,F0 space'
                            xdat(1).x=h_crit2./H0_crit2;
                            ydat(1).y=F0_crit2; %using Smith formula for F0 for an inversion F0=U/sqrt(g'H0)
                            labs(1).l='Original inversion case';
                            
                            %now give range of h' values and calculate H0 using (20). Note F0=1/H0' since H0' = H0*N/U
                            xdat(2).x=2*hhat./H0_hat_strat;
                            ydat(2).y=sqrt(2)*F0_crit; %plot the original formulation for continuous stratification with transformed variables
                            labs(2).l='Stratified formula with Smith&Sun transforms';
                            % i.e. using Heff instead of H0 and F0=U/(N*Heff)=sqrt(2)U/(N*H0)

                            %        plot(1/(0.62^2)*hhat./H0_hat_strat,1/0.62*F0_crit,'k--');

                            xdat(3).x=hhat./H0_hat_strat/0.63;
                            ydat(3).y=F0_crit/0.63;
                            labs(3).l='Stratified formula with 0.63 transforms';
                            
                            xdat(4).x=hhat./H0_hat_strat/0.63;
                            ydat(4).y=F0_crit/0.63;
                            labs(4).l='Stratified formula with 0.63,0.86 transforms';
                            
                            xdat(5).x=hhat./H0_hat_strat;
                            ydat(5).y=1./(sqrt(2)*sin(0.5*H0_hat_strat));
                            labs(5).l='Stratified formula using Durran F0';
                            
                            xdat(6).x=2*hhat./H0_hat_strat;
                            ydat(6).y=1./(sqrt(2)*sin(0.5*H0_hat_strat));
                            labs(6).l='Stratified formula using Durran F0 and Heff';
                            
                            
                            ylab='F0';
                            xlab='h/H0';
                            
                            xlims=1;
                            xlimits=[0 1];
                                                        
                            izlim=1;
                            zmin=0;
                            zmax=1;
                            
                        case 'h vs H0'
                            xdat(1).x=hhat;
                            ydat(1).y=H0_crit*N/U; 
                            labs(1).l='Stratified exact solution';
                            
                            %now calculate effective F0 and plug into eq (16)
                            Feff=1/0.63 * U/N./H0_crit2;
                            h_crit2 = H0_crit2 .* ( 1 + 0.5*Feff.^2 - 1.5*Feff.^(2/3) );
                            
                            xdat(2).x=h_crit2*N/U;
                            ydat(2).y=H0_crit2*N/U; %plot the original formulation for continuous stratification with transformed variables
                            labs(2).l='0.63 transform into inversion formula';
                            % i.e. using Heff instead of H0 and F0=U/(N*Heff)=sqrt(2)U/(N*H0)

                            %        plot(1/(0.62^2)*hhat./H0_hat_strat,1/0.62*F0_crit,'k--');

%                            xdat(3).x=hhat./H0_hat_strat/0.63;
%                            ydat(3).y=F0_crit/0.63;
%                            labs(3).l='Stratified formula with 0.63 transforms';
                            
                            ylab='Hhat';
                            xlab='hhat';
                            
                            xlims=1;
                            xlimits=[0 1];
                                                        
                            izlim=1;
                            zmin=0;
                            zmax=5;

                    end
                    nmark=0;
                                      

                    figname=['Houghton plots'];
                    if i_paper_labels==1
                        figname=['Houghton/Smith curves'];
                    end
                    
                    
                    
                    
                    ascent_str='';
  
                    
                    


            end  
            
        
            
            if no_sort==0;
                [ydat(i).y I]=sort(ydat(i).y);
                xdat(i).x = xdat(i).x(I);
            end
            
            if ydat(i).y(1)<0 & is_met_em(1)==1 & incep==1
                ydat(i).y=ydat(i).y(2:end);
                xdat(i).x=xdat(i).x(2:end);
            end

         
        
            abc=['ABCDEFGHIJKLM'];

            if length(strfind(ascent_str,'1'))>0
                if i_label_time==1
                    if i_override_loc_lab==1
                        labs(i).l=[location_lab(i-nplots_previous).l ' ' day ' ' month_text ' ' hour ':' mins];  %simple label more suitable for papers
                    else
                        location_lab(i-nplots_previous).l = abc(i-nplots_previous);
                        labs(i).l=[location_lab(i-nplots_previous).l ' ' day ' ' month_text ' ' hour ':' mins];  %simple label more suitable for papers
                    end
                else
                    if i_override_loc_lab==1
                        labs(i).l=[location_lab(i-nplots_previous).l];  %simple label more suitable for papers
                    else
                        location_lab(i-nplots_previous).l = abc(i-nplots_previous);
                        labs(i).l=['W' location_lab(i-nplots_previous).l ' ' num2str(lat2d.var(ilat(iloc),ilon(iloc)),3) ' , ' num2str(lon2d.var(ilat(iloc),ilon(iloc)),3) ];
                    end
                end
            elseif length(strfind(ascent_str,'0'))>0 | length(strfind(ascent_str,'3'))>0 | length(strfind(ascent_str,'4') )>0 | length(strfind(ascent_str,'5') )>0
                if i_label_time==1
                    if i_override_loc_lab==1
                        labs(i).l=[location_lab(i-nplots_previous).l ' ' day ' ' month_text ' ' hour ':' mins];
                    else
                        location_lab(i-nplots_previous).l = abc(i-nplots_previous);
                        labs(i).l=[location_lab(i-nplots_previous).l ' ' day ' ' month_text ' ' hour ':' mins];
                    end
                else
                    if i_override_loc_lab==1
                        labs(i).l=[location_lab(i-nplots_previous).l ' ' num2str(lat2d.var(ilat(iloc),ilon(iloc)),3) ' , ' num2str(lon2d.var(ilat(iloc),ilon(iloc)),3)];  %simple label more suitable for papers
                    else
                        location_lab(i-nplots_previous).l = abc(i-nplots_previous);
                        labs(i).l=['' location_lab(i-nplots_previous).l ' ' num2str(lat2d.var(ilat(iloc),ilon(iloc)),3) ' , ' num2str(lon2d.var(ilat(iloc),ilon(iloc)),3) ];
                    end
                end
            end
            
       
        end 
       
     end
%        labs(i).l=abc(i);
    end
    
    
    
    if roth_comp == 1
        
%fields are
    % 1) YEAR 
    % 2) MONTH 
    % 3) DAY 
    % 4) HOUR 
    % 5) PRESSURE 
    % 6) HEIGHT 
    % 7) TEMPERATURE 
    % 8) DEWPOINT 
    % 9) WIND_DIRECTION 
    %10) WIND_SPEED
    
    
    inds_11 = find(dat_roth(4,:)==11);
    inds_12 = find(dat_roth(4,:)==12);
    labs(i+1).l = 'Rothera sounding 05 Jan 12:00';
    labs(i+2).l = 'Rothera sounding 06 Jan 11:00';


        switch var
            case 'pressure diff'
                xlab='Pressure difference, WRF minus aircraft (hPa)';

                wrf_int = interp1(ydat(i).y,xdat(1).x,dat(ds_inds,11),[],'extrap');
                xdat(i).x = wrf_int - dat(ds_inds,4);
                ydat(i).y = dat(ds_inds,11);
                

                xlims=0;
                xlimits=[965 991];
                %		xlimits=[800 995];

                izlim=1;
                zmin=0;
                zmax=3000;

                %		logflag=1;
                if iwrf_profs==0
                    figname=['Pressure profile '];
                end
                
             case 'pressure'
                xlab='Pressure (mb)';

                xdat(i+1).x = dat_roth(5,inds_11); 
                ydat(i+1).y = dat_roth(6,inds_11);
                

                xlims=0;
                xlimits=[965 991];
                %		xlimits=[800 995];

                izlim=1;
                zmin=0;
                zmax=3000;

                %		logflag=1;
                if iwrf_profs==0
                    figname=['Pressure profile '];
                end   
                
            case 'temperature'
                xlab='Temperature (^{o}C)';                                
                xdat(i+1).x = dat_roth(7,inds_12); 
                ydat(i+1).y = dat_roth(6,inds_12);
                
                xdat(i+2).x = dat_roth(7,inds_11); 
                ydat(i+2).y = dat_roth(6,inds_11);

                
                izlim=1;
                zmin=0;
                zmax=3000;
                
                xlims=1;
                xlimits=[-14 2];
                xlimits=[-15 5];
                
                if iwrf_profs==0
                    figname=['Temperature profile '];
                end
                
            case 'vapour'
                xlab='Vapour mixing ratio (g kg^{-1})';
                
                Tdew = dat_roth(8,inds_12);
                xdat(i+1).x = 1000*vap_from_Tdew(Tdew+273.15,100*dat_roth(5,inds_12));
                ydat(i+1).y = dat_roth(6,inds_12);
                
                Tdew = dat_roth(8,inds_11);
                xdat(i+2).x = 1000*vap_from_Tdew(Tdew+273.15,100*dat_roth(5,inds_11));
                ydat(i+2).y = dat_roth(6,inds_11);

                izlim=1;
                zmin=0;
                zmax=3000;
                
                xlims=1;
                xlimits=[0 20];            
                
                izlim=1;
                zmin=0;
                zmax=3000;
                
                xlims=1;
                xlimits=[0 5];
                
                if iwrf_profs==0
                    figname=['Vapour profile '];
                end
                
            case 'wind speed'
                xlab='Wind speed (m s^{-1})';

                xdat(i+1).x = dat_roth(10,inds_12); 
                ydat(i+1).y = dat_roth(6,inds_12);
                
                xdat(i+2).x = dat_roth(10,inds_11); 
                ydat(i+2).y = dat_roth(6,inds_11);

                izlim=1;
                zmin=0;
                zmax=3000;
                
                xlims=1;
                xlimits=[0 20];
                
                if iwrf_profs==0
                    figname=['Wind speed profile '];
                end
                
            case 'wind dir'
                xlab='Wind direction (degrees)';
                xdat(i+1).x = dat_roth(9,inds_12); 
                ydat(i+1).y = dat_roth(6,inds_12);
                
                xdat(i+2).x = dat_roth(9,inds_11); 
                ydat(i+2).y = dat_roth(6,inds_11);

                izlim=1;
                zmin=0;
                zmax=3000;
                
                xlims=0;
                xlimits=[0 20];
        end
    
        
    i=i+2;
        
    end    

    





    %	xdat(i+1).x = dat2(5,1:12);  %data in hPa
    %	ydat(i+1).y = dat2(6,1:12); %height
    %	labs(i+1).l='Rothera station data';






    
    if is_met_em(1) & i_paper_labels==0
        savename = ['analysis ' figname];
        titlenam = [figname ' (analysis)'];
    else
        savename=figname;
        titlenam = figname;        
    end
    
    if no_title==1
        titlenam='';
    end
    
%    savename = [savename ' ' as_ds_str];
    
    
    



case 77
        %run av_updraught first
    
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
    
%     qstr=1; %vapour
     qstr=2; %liq
%     qstr=3; %rain
%     qstr=4; %snow
%     qstr=5; %graupel
%     qstr=6; %ice
%     qstr=7; %ice NC
%     qstr=8; %graupel NC
%     qstr=9; %snow NC
%     qstr=10; %tracer
%     qstr=11; %graupel density

if qstr<7
    units=' (g kg^{-1})';
    units=' (g m^{-3})';

else
    units=' (kg^{-1})';
end

    itemp=1; %flag to say that want to use the temperature as the vertical coordinate
            
    xlab=['Max Q0' num2str(qstr) units];
	ylab='Height (km)';
    
    if itemp==1;	
        ylab='Temperature (^{o}C)';
        iydir=-1; %reverse direction of y axis
        izlim=1;
        zmin=-45; 
        zmax=15;
    end
    
%    figname=['Max'];
    figname=xlab;
    savename=figname;
    titlenam=savename;
    
    logflag=0;
 
    idats=[1 2];
	for iidat=1:length(idats)
       idat=idats(iidat); 
       
       if qstr<7; 
           factor=1e3.*GridDan(idat).RHO; 
       else
           factor=1;
       end %if want to convert g/m3
        
        [Y,I]=max(ThreeDDan(idat).Q(2:end,2:end,:,qstr),[],1); %find max and index of max for 1st dimension
        I=squeeze(I);  %size=e.g. [97 150]
        [Y2,I2]=max(Y,[],2);  %find max and index of these maxes (2nd dimension of ThreeD)
        I2=squeeze(I2);  %size= [150 1]
                        
        xdat(iidat).x = factor.*squeeze(Y2); %max values stored in Y2      
        ydat(iidat).y = (GridDan(idat).Z + add_ground_height)/1e3 ;     
        
        if itemp==1
               
                clear P_cloud th_cloud T_cloud
                for izqmax=1:length(Y2)
                    th_cloud(izqmax) = ThreeDDan(idat).TH1( I(I2(izqmax),izqmax) , I2(izqmax), izqmax ) ;                    
                    [T_cloud(izqmax), P_cloud(izqmax)] = temp_from_press_and_th( GridDan(idat), th_cloud(izqmax) , ...
                        ThreeDDan(idat).P( I(I2(izqmax),izqmax) , I2(izqmax), izqmax )  , izqmax);
                    %using the actual 3D temperature at the grid point of the max value
                    %I contains the index of the max in 1st dim and I2 for 2nd one
                    
               % T=tempLES(GridDan(idat))-273.15;
                end
                ydat(iidat).y = T_cloud-273.15;
        end            
        
        labs(iidat).l = runName(idat).nam;        	
	end		
    
         f=1e6*28.97/18;
    
    
        T=tempLES(GridDan(idat))-273.15;
        
        iz_th=2;
%        th_base = GridDan(1).THREF(iz_th) + GridDan(1).OLTHBAR(iz_th);
        th_base = th_cloud(iz_th) + GridDan(1).THREF(iz_th);
        
        
        pdat=[GridDan(1).PREFN(2):-10:200]; %set up pressure grid        
        Tdry=th_base./(1e5./pdat).^0.286; %dry adiabat from th_base
   
		
        

       
  icb_temp=1;     
     
     if icb_temp==0   %tries to work out cloud base from moisture in BL
        qv_base=GridDan(1).OLQBAR(20,1); 
        qsat=satvappress(Tdry,'goff','liq',pdat,1)/f; %saturation mixing ratio during dry adiabatic rise
        icb_pa = findheight(qsat,qv_base);
        cb_pa=pdat(icb_pa);
        cb_k = Tdry(icb_pa);
    else           %OR specify cb temp here and works out cb pressure based on dry adiabat 
       cb_k=8+273.15; %temp of cloud base in degrees  
       idry=findheight(Tdry,cb_k);
       cb_pa = pdat(idry);
%       cb_pa = pdat(idry)*1.3;
       idry2 = findheight(T_cloud,cb_k);  %or use the pressure at the LEM gridpoint where the max LWC was and the temperature = cb_k
       cb_pa = P_cloud(idry2);              
     end
     
   %OR base on cloud base according to LEM LWC field  
     ilwc = find(xdat(1).x>0.01); %find first height in LEM grid where LWC starts to form and use this as cloud base
     cb_pa = P_cloud(ilwc(1));
     cb_k = T_cloud(ilwc(1));
     
        
        icb=findheight(GridDan(1).PREFN,cb_pa);
%        icb=findheight(T,cb_pa);

    clear gt_k gt_pa ct_k
           
        gt_pa=GridDan(1).PREFN(icb:end);
        gt_k=T(icb:end)+273.15;  %note: routinte doesn't really use gt_k data
        ct_k=T(end)+273.15;
        
        clear adiabatic_lwc adjusted_lwc adjusted_slwc adtemp
        for iad=1:length(gt_k)
          [adiabatic_lwc(iad),adtemp(iad)] ...
            = adLWC_PaulLawson_simple ( cb_pa,cb_k,gt_pa(iad) ); %returns adiabatic LWC in g/m3
        end
        
        rho = gt_pa'.*28.97e-3/8.3144./(adtemp+273.15); %useful for converting from g/m3 to g/kg
        xdat(iidat+1).x = 0.5*adiabatic_lwc;
        ydat(iidat+1).y = adtemp;   %NOTE - need to plot against the moist adiabatic temperature since the aircraft measurements are
                                    %taken IN CLOUD - so should be different to ambient temperature!
        labs(iidat+1).l = 'half adiabatic LWC';
%         ydat(iidat+1).y = (GridDan(idat).Z(2:end) + add_ground_height)/1e3 ; 
         
%          if itemp==1  %plot LEM LWC as function of adiabatic temperature - actually just plot with actual grid point LEM temperatrure
%                 %T=tempLES(GridDan(idat))-273.15;
%                 %ydat(iidat+1).y = T(icb:end);   
%              
%                 for iidat=1:length(idats)
%                     idat=idats(iidat); 
%                     for iad=1:length(adtemp)   
%                         iiad=findheight(GridDan(idat).PREFN,pdat
%                     
%          end
         
         
         
         
             
        
         xlims=0;
         xlimits=[0 2.8];
        lor=-1;

	%     
         
         nmark=0;
    
case 76
    %WRF wind profile direction plots
    
    dual=0;
    lor=2;
    
    izlim=1;
    zmin=0;
    zmax=1600;
    	
	xlims=0;
%	xlimits=[0 3600];
%	xlimits=[0 1.5e5];
    
    
    logflag=0;
    
    time=16;
    tstr=Times(time,:);
    iund=findstr('_',tstr);
    tstr(iund)=' ';
    
    for i=1:length(ilat)	
		ydat(i).y=WRFUserARW(nc,'Z',time,ilat(i),ilon(i));
		u=WRFUserARW(nc,'umeta',time,ilat(i),ilon(i));
        v=WRFUserARW(nc,'vmeta',time,ilat(i),ilon(i));
        
        [th,r]=cart2pol(u,v); %th : minus indicates angles anti-clockwise from 0 to 180
        
%        xdat(i).x= mod(360 + th*180/pi, 360);
        xdat(i).x= (360 + th*180/pi);
        
        labs(i).l=[num2str(LAT(i)) ' , ' num2str(LON(i)) ];    
    end

        
   
    xlab='Direction (degrees)';    

    
    figname=['Wind direction profile at ' tstr];
    
	ylab='Height (km)';
	titlenam=figname;
    savename=figname;  
    
    
    
    

case 75
    %WRF wind profile plots
    
    dual=0;
    lor=-1;
    
    izlim=0;
    zmin=0;
    zmax=1600;
    	
	xlims=0;
%	xlimits=[0 3600];
%	xlimits=[0 1.5e5];
    
    
    logflag=0;
    
    time=36;  %15,16,17 = 18,21,0 UTC
    tstr=Times(time,:);
    iund=findstr('_',tstr);
    tstr(iund)=' ';
    
    for i=1:length(ilat)	
		ydat(i).y=WRFUserARW(nc,'Z',time,ilat(i),ilon(i));
		u=WRFUserARW(nc,'u',time,ilat(i),ilon(i));;
        v=WRFUserARW(nc,'v',time,ilat(i),ilon(i));
        
        xdat(i).x= sqrt( u.^2 + v.^2 );
        
        labs(i).l=[num2str(LAT(i)) ' , ' num2str(LON(i)) ];    
    end

        
   
    xlab='Speed (m/s)';    

    
	xdat(i+1).x = dat2(10,1:12)*0.5144444; %presumably data is in knots
	ydat(i+1).y = dat2(6,1:12);
	labs(i+1).l='Rothera station data';



    figname=['Wind speed profile at ' tstr];
    
	ylab='Height (km)';
	titlenam=figname;
    savename=figname;    

case 744
    %WRF mean profile plots
    
    dual=0;
    lor=1;
    
    izlim=0;
    zmin=13e3;
    zmax=21e3;
    	
	xlims=0;
	xlimits=[-10 40];
%	xlimits=[0 1.5e5];
    
    
    logflag=0;
    
    time=idir;
    time=23;
    
    tstr=Times(time,:);
    iund=findstr('_',tstr);
    tstr(iund)=' ';
    
    hm='total_water';
%    hm='ice';
%     hm='temp';
%     hm='iceno';

    switch hm
    case 'total_water'
        xlab='Total water mixing ratio (ppmv)';   
        factor=f;
        figname=['Mean total water mixing ratio profile at ' tstr];
        izlim=1;
        zmin=13e3;
        zmax=21e3;
    	xlims=1;
		xlimits=[-10 40];
        
    case 'ice'
        xlab='Ice mixing ratio (kg kg^{-1})';   
        xlab='Ice mixing ratio (ppmv)';   
        factor=f;
        figname=['Max ice mixing ratio profile at ' tstr];
        izlim=1;
        zmin=13e3;
        zmax=21e3;
    	xlims=1;
		xlimits=[-10 10];
        
    case 'temp'
        xlab='Temp (^{o}C)';    
        factor=1;
        figname=['Temperature profile at ' tstr];
        
    case 'iceno'
        xlab='Ice number concentration (kg^{-1})';   
        factor=1;
        figname=['Max ice number profile at ' tstr];
        
    end
    
    
    
%    
    
    
    
    
    for i=1:length(LAT)	       
        switch hm
        case 'total_water'
            xdat(i).x=factor.*mean(mean(nc{'QICE'}(time,:,:,:)+nc{'QVAPOR'}(time,:,:,:),2),3); %produces an [nz nx ny] sized array
        case 'ice'
            [Y,I]=max(nc{'QICE'}(time,:,:,:),[],2); %produces an [nz nx ny] sized array
        end                

        
        z_wrf=WRFUserARW(nc,'Z',time);
        ydat(i).y=squeeze(mean(mean(z_wrf.var,2),3));
        

		    

        %labs(i).l=[num2str(LAT(i)) ' , ' num2str(LON(i)) ];    
        labs(i).l='Mean';
    end

        
   
   
    
	ylab='Height (km)';
	titlenam=figname;
    savename=figname;
    
case 74
    %WRF max profile plots
    
    dual=0;
    lor=1;
    
    izlim=0;
    zmin=13e3;
    zmax=21e3;
    	
	xlims=0;
	xlimits=[-10 40];
%	xlimits=[0 1.5e5];
    
    
    logflag=0;
    
    time=idir;
    time=20;
    
    tstr=Times(time,:);
    iund=findstr('_',tstr);
    tstr(iund)=' ';
    
    hm='total_water';
%    hm='ice';
%     hm='temp';
%     hm='iceno';

    switch hm
    case 'total_water'
        xlab='Total water mixing ratio (ppmv)';   
        factor=f;
        figname=['Max total water mixing ratio profile at ' tstr];
        izlim=0;
        zmin=13e3;
        zmax=21e3;
        
        zmin=300;
        zmax=400;
        
    	xlims=1;
		xlimits=[-10 100];
        
    case 'ice'
        xlab='Ice mixing ratio (kg kg^{-1})';   
        xlab='Ice mixing ratio (ppmv)';   
        factor=f;
        figname=['Max ice mixing ratio profile at ' tstr];
        izlim=1;
        zmin=13e3;
        zmax=21e3;
    	xlims=1;
		xlimits=[-10 10];
        
    case 'temp'
        xlab='Temp (^{o}C)';    
        factor=1;
        figname=['Temperature profile at ' tstr];
        
    case 'iceno'
        xlab='Ice number concentration (kg^{-1})';   
        factor=1;
        figname=['Max ice number profile at ' tstr];
        
    end
    
    
    
%    
    
    
    
    
    for i=1:length(LAT)	       
        switch hm
        case 'total_water'
            [Y,I]=max(nc{'QICE'}(time,:,:,:)+nc{'QSNOW'}(time,:,:,:)+nc{'QGRAUP'}(time,:,:,:)+nc{'QVAPOR'}(time,:,:,:),[],2); %produces an [nz nx ny] sized array
        case 'ice'
            [Y,I]=max(nc{'QICE'}(time,:,:,:),[],2); %produces an [nz nx ny] sized array
        end
        
        %now I is a vector of size [nz nz]
        I=squeeze(I);
        [Y2,I2]=max(Y,[],3);  %find max and index of these maxes (3rd dimension of 3D array)
        I2=squeeze(I2);  %size= [nz 1]
        
        xdat(i).x = factor.*squeeze(Y2); 
        
        
        
        for iz=1:length(Y2)
           % [xdat(i).x(iz) imax]=maxALL(nc{'QNICE'}(time,iz,:,:));
%            [xdat(i).x(iz) imax]=maxALL(nc{'QICE'}(time,iz,:,:));
            
%		    ydat(1).y(iz)=WRFUserARW(nc,'Z',time,I(iz,I2(iz)),I2(iz),iz);
            ydat(1).y(iz)=WRFUserARW(nc,'th',time,I(iz,I2(iz)),I2(iz),iz);

%		xdat(i).x=WRFUserARW(nc,'tc',time,ilat(i),ilon(i));
		end
		    

        %labs(i).l=[num2str(LAT(i)) ' , ' num2str(LON(i)) ];    
        labs(i).l='Max';
    end

        
   
   
    
	ylab='Height (km)';
	ylab='Potential temperature (K)';

	titlenam=figname;
    savename=figname;
    
    
   
    case 73
    %MAC3 plots
    
    dual=0;
    lor=-1;
    
    izlim=0;
    zmin=0;
    zmax=16;
    	
	xlims=0;
	xlimits=[0 3600];
%	xlimits=[0 1.5e5];
    
    
    logflag=0;
    
dz=120; %vertical grid spacing
it=3;

hm='ni';
hm='nd';
hm='ql';

times=[5:5:60];

for itimes=1:length(times)
   % comm=['xdat(' num2str(itimes) ').x=max(mac3(' num2str(itimes) ').' hm ',[],1);'];
    comm=['xdat(' num2str(itimes) ').x=mean(mac3(' num2str(itimes) ').' hm ',1);'];
    
    eval(comm);
    labs(itimes).l=[num2str(times(itimes)) ' mins'];
end
    
%     comm=['xdat(2).x=max(mac3(4).' hm ',[],1);'];
%     eval(comm);
%     labs(2).l='20 mins';
%     
%     comm=['xdat(3).x=max(mac3(5).' hm ',[],1);'];
%     eval(comm);
%     labs(3).l='25 mins';
    
    
    for ix=1:length(xdat)
        ydat(ix).y=[dz:dz:126*dz]/1000;
	end
        
    if strcmp(hm(1),'q')==1    
       	xlab='Mixing Ratio (g kg^{-1})';    
    else
       	xlab='Number concentration (cm^{-3})'; 
        for ix=1:length(xdat)
            xdat(ix).x=xdat(ix).x / 1e6;
        end
    end
    
    figname=[hm ' max profile'];
    figname=[hm ' mean profile'];
    
	ylab='Height (km)';
	titlenam=figname;
    savename=figname;
    
    
 
    
    
    
    case 72
    %running times
    
    dual=0;
    lor=2;
    
    izlim=0;
    zmin=14;
    zmax=22;
    	
	xlims=1;
	xlimits=[0 3600];
%	xlimits=[0 1.5e5];
    
    
    logflag=1;
    


    ydat(1).y=[0 13+59/60 22+51/60 36+41/60 45+9/60];
    labs(1).l='Dan';
    
    ydat(2).y=[0 13+21/60 21+57/60 34+58/60 42+51/60];
    labs(2).l='Grant';
    
    ydat(3).y=[0 14+40/60 23+58/60 37+40/60 45+19/60];
    labs(3).l='Dave W';
    
    ydat(4).y=[0 16+09/60 26+43/60 44+16/60 55+13/60];
    labs(4).l='Dave T';
    
    ydat(5).y=[0 15+50/60 26+45/60 43+24/60 52+39/60];
    labs(5).l='Jonny';
    
    ydat(6).y=[0 10+51/60 18+49/60 31+21/60 38+49/60];
    labs(6).l='Paul';    

    xkm=[0 3 5 8 10];
    for ix=1:6
        xdat(ix).x=ydat(ix).y*60;
        ydat(ix).y=xkm;
	end
        
    
   	xlab='Time (mins)';    
    figname='Time breakdown';
    
	ylab='Distance (km)';
	titlenam=figname;
    savename=figname;
    
    
    case 71
    %running times
    
    dual=0;
    lor=2;
    
    nmark=-1;
    
    izlim=0;
    zmin=14;
    zmax=22;
    	
	xlims=0;
	xlimits=[4.95 5.05];
	xlimits=[0 1.5e5];
    
    
    logflag=0;

    clear times
    
year=2009;
switch year
    case 2008
    times(1).y=[0 13+59/60 22+51/60 36+41/60 45+9/60];
    labs(1).l='Dan';
    
    times(2).y=[0 13+21/60 21+57/60 34+58/60 42+51/60];
    labs(2).l='Grant';
    
    times(3).y=[0 14+40/60 23+58/60 37+40/60 45+19/60];
    labs(3).l='Dave W';
    
    times(4).y=[0 16+09/60 26+43/60 44+16/60 55+13/60];
    labs(4).l='Dave T';
    
    times(5).y=[0 15+50/60 26+45/60 43+24/60 52+39/60];
    labs(5).l='Jonny';
    
    times(6).y=[0 10+51/60 18+49/60 31+21/60 38+49/60];
    labs(6).l='Paul';    
    
    times(7).y=[0 15+55/60 29+04/60 49+31/60 60+38/60];
    labs(7).l='Liz';    
    
    times(8).y=[0 18+25/60 30+39/60 51+37/60 64+4/60];
    labs(8).l='Joz';  
    
    xkm=[0 3 5 8 10];
    
    case 2009
        xkm=[0 3 5 7.5 10];
        
        idat=0;
        
        idat=idat+1;        
        labs(idat).l='Haile';
        times(idat).y=[0 00+08+08/60 00+13+30/60 00+20+35/60 00+27+39/60];
        
        idat=idat+1;
        times(idat).y=[0 11+05/60 18+32/60 28+26/60 38+03/60]; 
        labs(idat).l='Paul'; 
        
        idat=idat+1;        
        labs(idat).l='Grant';
        times(idat).y=[0 00+12+10/60 00+21+04/60 00+32+31/60 00+43+19/60];
        
        idat=idat+1;        
        labs(idat).l='Dave W';
        times(idat).y=[0 00+14+15/60 00+23+40/60 00+35+28/60 00+46+29/60];    %4570   2329 1500 2171 1392
        
        idat=idat+1;        
        labs(idat).l='Ian';
        times(idat).y=[0 00+15+31/60 00+24+51/60 00+36+40/60 00+48+14/60];   %15792   3485 2191 3181 1988
                
        idat=idat+1;
        times(idat).y=[0 15+47/60 26+08/60 39+20/60 52+22/60];
        labs(idat).l='Chris'; 
        
        idat=idat+1;        
        labs(idat).l='Cat';
        times(idat).y=[0 00+16+54/60 00+28+02/60 00+42+26/60 00+56+08/60];
                        
        idat=idat+1;        
        labs(idat).l='Dave T';
        times(idat).y=[0  00+16+21/60 00+28+43/60 00+44+23/60 60+04/60];                
        
        idat=idat+1;        
        labs(idat).l='Hugo';
        times(idat).y=[0 00+19+16/60 00+32+28/60 00+48+59/60 60+04+54/60]; 
                        
        idat=idat+1;        
        labs(idat).l='Rachel';
        times(idat).y=[0 00+19+21/60 00+32+59/60 00+51+06/60 60+08+51/60]; 
        
        
        
        
        
        
        
end

    
    for ix=1:length(times)
%         xdat(ix).x=xkm;
%         speeds=diff(xkm)./diff(ydat(ix).y);
%         ydat(ix).y=[speeds(1) speeds]*60; %*60 to convert from km/min to km/hr

         xdat(ix).x=xkm(2:end);
         speeds=diff(xkm)./diff(times(ix).y);
         ydat(ix).y=[speeds]*60; %*60 to convert from km/min to km/hr
         
         mins=floor(times(ix).y(end));
         if mins>=60
             hours = floor(mins/60);
             if mins-hours*60<10
                 mins_str = [num2str(hours) ':0' num2str(mins-hours*60)];
             else
                 mins_str = [num2str(hours) ':' num2str(mins-hours*60)];
             end
         else
             hours=0;
             mins_str = num2str(mins);
         end
         
         
         
         secs=60*(times(ix).y(end)-mins);
         if secs<10
             secs_str = ['0' num2str(secs)];
         else
             secs_str = [num2str(secs)];
         end
         
         labs(ix).l = [labs(ix).l ' ' mins_str ':' secs_str];
    end
    
    xlims=1;
	xlimits=[0 10];
        
    
   	ylab='Average Speed (km hr ^{-1})';    
%   	ylab='Average Speed (km min^{-1})';    

    figname='Average running speeds';
    
	xlab='Distance (km)';
	titlenam=figname;
    savename=figname;
    
    
    case 70
    %min vap for 3d case
    
    dual=0;
    lor=4;
    
    izlim=1;
    zmin=14;
    zmax=22;
    	
	xlims=1;
	xlimits=[4.95 5.05];
	xlimits=[0 1.5e5];
    
    
    logflag=0;
    

    xdat(1).x=squeeze(max(max(ThreeDDan(1).Q(2:end-1,2:end-1,:))));
    ydat(1).y=Grid.Z/1000 + add_ground_height;
    labs(1).l='3D'; 
    
   	ylab='Height (km)';
    
    figname='Max ice NC';
    
	xlab='Max ice number concentration (kg^{-1})';
	titlenam=figname;
    savename=figname;
    
    case 69
        %run av_updraught first
    
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
 
    xlab=['Time UTC'];
	ylab='Mass Flux (kg s^{-1} m^{-1})';
    
    figname=['Updraught mass flux'];
    savename=figname;
    titlenam=savename;
    
    logflag=0;
 
    idats=[1:5];
	for iidat=1:length(idats)
       idat=idats(iidat); 
                
        ydat(iidat).y = 1000*wflux_diag(idat).w;
        xdat(iidat).x = GridDan(idat).t(1:length(ydat(iidat).y)) + 3 ;         
        labs(iidat).l=runName(idat).nam;        	
	end		
         xlims=1;
         xlimits=[GridDan(1).t(1) GridDan(1).t(60)]+3;
	%     
         
         nmark=0;
         
    case 68
    
    time1=19.75;
    time2=23.75;
    t1=findheight(GridDan(1).t+3,time1);
    t2=findheight(GridDan(1).t+3,time2);
    
    t1=1;
    t2=20;
    
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
 
    xlab=['Eddy Heat Flux Contribution'];
	ylab='Height (km)';
    
    figname=['Eddy Heat flux'];
    savename=figname;
    titlenam=savename;
    
    logflag=0;
    
    error=0;
    [isg,error]=getDGcol('WTHSG',dgstrDan(1).dg,error); %subgrid
    %[icn,error]=getDGcol('WTHCN',dgstrDan(1).dg,error); %resolved
    [icn,error]=getDGcol('WTHAD',dgstrDan(1).dg,error); %resolved
    [iwth,error]=getDGcol('ALL_WTH',dgstrDan(1).dg,error); %resolved
    [iwthsg,error]=getDGcol('ALL_WTHSG',dgstrDan(1).dg,error); %resolved
    [ivw,error]=getDGcol('VW',dgstrDan(1).dg,error); %resolved
    [ivwsg,error]=getDGcol('VWSG',dgstrDan(1).dg,error); %resolved
      
    [iwthA,error]=getDGcol('ALd_A',dgstrDan(1).dg,error); %resolved
    
    [idisr,error]=getDGcol('DISR',dgstrDan(1).dg,error); %resolved
    
    aind=280;%ALu_A
    aind=284; %ACC_A
    aind=iwthA;
    aind=[];
    
   % aind=[];
    
%    [iqsg,error]=getDGcol('WQ01SG',dgstrDan(1).dg,error); %subgrid
%    [iadsg,error]=getDGcol('WQ01AD',dgstrDan(1).dg,error); %resolved
    
    datind=[isg icn]; %sum of subgrid and resolved
  %  datind=[idisr]; 
  %  datind=[iqsg iadsg]; %sum of subgrid and resolved turbulent flux of vapour
  %  datind=[ivw ivwsg]; 
    
    idats=[1 2];
    
    ih=3;
 
    
	for iidat=1:length(idats)
       idat=idats(iidat); 
        
          ih2=length(GridDan(idat).Z);
          
         if length(aind)==1            
            area=TimeAvDan(idat).DGAV(ih:ih2,aind);
            area(area==0)=1;
        else
            area=1;
        end
        
%        xdat(iidat).x= - 1/300*TotMassBudgetALL(GridDan(idat),sum(TimeAvDan(idat).DGAV(:,datind),2),GridDan(idat).t,ih-1,ih2);
        
        xdat(iidat).x= sum(TimeAvDan(idat).DGAV(ih:ih2,datind),2)/npess2(iidat);
        
       % xdat(iidat).x = mean(sum(icediagsALL(idat).i(:,t1:t2,[datind]),3)./area,2)/npess2(idat); %ALu_W. Dividing by no. processors

        ydat(iidat).y = (GridDan(idat).Z(ih:ih2)+620)/1000; 
        
        labs(iidat).l=runName(idat).nam;
        	
	end
	

	
         xlims=0;
         xlimits=[-30 50];
	%     
         
         nmark=0;
         
    case 67
    
    time1=19.75;
    time2=23.75;
    t1=findheight(GridDan(1).t+3,time1);
    t2=findheight(GridDan(1).t+3,time2);
    
    t1=1;
    t2=20;
    
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
 
    xlab=['Tracer Flux'];
	ylab='Height (km)';
    
    figname=['Tracer flux'];
    savename=figname;
    titlenam=savename;
    
    logflag=0;
    
    aind=280;%ALu_A
    aind=284; %ACC_A
    
    aind=[];
    
    datind=151; %ALL_Q10
    datind=157; %ALu_WQ10
    datind=139; %ALL_WQ10
    datind=[145 139]; %ALL_WQSG10
    
    idats=[1 2];
    
	for iidat=1:length(idats)
       idat=idats(iidat); 

         if length(aind)==1            
            area=icediagsALL(idat).i(:,t1:t2,aind)/npess2(idat);
            area(area==0)=1;
        else
            area=1;
        end

        xdat(iidat).x = mean(sum(icediagsALL(idat).i(:,t1:t2,[datind]),3)./area,2)/npess2(idat); %ALu_W. Dividing by no. processors

        ydat(iidat).y = (GridDan(idat).Z+620)/1000; 
        
        labs(iidat).l=runName(idat).nam;
        	
	end
	

	
         xlims=0;
         xlimits=[-30 50];
	%     
         
         nmark=0;
         
    case 66
    %min vap for 3d case
    
    dual=0;
    lor=4;
    
    izlim=1;
    zmin=14;
    zmax=22;
    	
	xlims=0;
	xlimits=[4.95 5.05];
	xlimits=[-0.02 0.02];
    
    
    logflag=0;
    
    meanvap=f*squeeze(mean(mean(ThreeD.Q(2:end-1,:,240,1),1),2));

    xdat(1).x=f*squeeze(min(min(ThreeD.Q(2:end-1,:,:,1),[],1),[],2)) - meanvap;
    ydat(1).y=Grid.Z/1000;
    labs(1).l='Min Vap'; 
    
    xdat(2).x=f*squeeze(max(max(ThreeD.Q(2:end-1,:,:,1),[],1),[],2)) - meanvap;
    ydat(2).y=Grid.Z/1000;
    labs(2).l='Max Vap'; 
    
    xdat(3).x=f*squeeze(mean(mean(ThreeD.Q(2:end-1,:,:,1),1),2)) - meanvap;
    ydat(3).y=Grid.Z/1000;
    labs(3).l='Mean Vap'; 

   
	ylab='Height (km)';
    
    figname='3D min';
    
	xlab='Min Vap (ppmv)';
	titlenam=figname;
    savename=figname;
    
    case 65
    %ozone DMI all plots
    %get from loadvapdata
    
    dual=0;
    lor=4;
    
    izlim=1;
    zmin=0;
    zmax=18;
    	
	xlims=0;
	xlimits=[0 20];
    
    
    logflag=0;
 

    xdat(1).x=Grid.OLQBAR(:,7);
    ydat(1).y=Grid.Z/1000;
    labs(1).l='Hazy';
    
    xdat(2).x=Grid.OLQBAR(:,8);
    ydat(2).y=Grid.Z/1000;
    labs(2).l='Hazy 2';

   
	ylab='Height (km)';
    
    figname='Hazel''s plot';
    
	xlab='Ice No. Conc (#/kg)';
	titlenam=figname;
    savename=figname;
    
case 64
    %ozone DMI all plots
    %get from loadvapdata
    
    dual=1;
    lor=4;
    
    zmin=11;
    zmax=30;
    
    izlim=1;
	
	xlims=1;
	xlimits=[0 20];
    
    
    logflag=0;
 
for idat=1:10
    xdat(idat).x=data(idat).dmi(8,:);
    ydat(idat).y=data(idat).dmi(1,:)/1000;
end
    
labs(1).l='flight_10_24Feb';
labs(2).l='flight_1_10Feb';
labs(3).l='flight_2_12Feb';
labs(4).l='flight_3_13Feb';
labs(5).l='flight_4_16Feb';
labs(6).l='flight_5_17Feb';
labs(7).l='flight_6_19Feb';
labs(8).l='flight_7_21Feb';
labs(9).l='flight_8_21Feb';
labs(10).l='flight_9_23Feb';
    
  
 %   xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='DMI ozonesondes';
    
	xlab='Ozone Mixing Ratio (ppmv)';
	titlenam=figname;
    savename=figname;
   
    
    
case 63
    %SF4 vapour and ozone plots
    %run readOzoneProfileSSS
    %or get from loadvapdata
    %run readsdla
    
    z=GridDan(1).Z+620;
    secyA=z/1000;
    secyB=z/1000;
    lab2='';  
    dual=2;
    
    xloc=[1 1 0 0];
    
    
    
    lor=4;
    
    zmin=11;
    zmax=30;
    
    izlim=1;
	
	xlims=1;
	xlimits=[0 20];
    
    
    logflag=0;
    
    xdat(1).x=sdla(1).s(7,:);
    ydat(1).y=sdla(1).s(2,:)/1000;
    labs(1).l='SDLA 21:57-00:48 UTC';
    
    xdat(2).x=f*dmi(1).p(:,10);
    ydat(2).y=dmi(1).p(:,1)/1000;
    labs(2).l='DMI 20:15-22:08 UTC';
    
     xdat(3).x=dirac2(10,:);
     ydat(3).y=dirac2(4,:)/1000;
     labs(3).l='SSS Ozone';
     
     xdat(4).x=data(10).dmi(8,:);
     ydat(4).y=data(10).dmi(1,:)/1000;
     labs(4).l='DMI Ozone';
    
  
    xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='24th Feb vapour and ice saturation mixing ratios';
    titlenam='Ozone Mixing Ratio (ppmv)';
    savename=figname;

    
    [zvap,I,J]=unique(ydat(2).y);
    vap=xdat(2).x(I);
    
    [zoz I J]=unique(ydat(4).y);
    
    vap2=interp1(zvap,vap,zoz);
    oz=xdat(4).x(I);
    
    
case 62 %aerosol distribution
    idir=1;
    
    izlim=0;

    
	for idat=1:length(Sc)
        xdat(idat).x = 1e6*exp(logD(idat).d(2:end));
        ydat(idat).y = dN(idat).n/1e6;
	end
    
   labs(1).l='Nuclei mode';
   labs(2).l='Accumulation mode';
   
   
	ylab='dN (dlogD)^{-1} (cm^{-3})';
    xlab=['Diameter (microns)'];
    
    
    logflag=1;
    
    
     xlims=0;
	 lor=1;
     nmark=0;
     
     savename=['Aerosol distribution'];
     titlenam=savename;
     
     
    
 case 61
    idir=1;
    
    izlim=0;

	[iz,iz2]=findheight(GridDan(1).Z+620,14e3,20e3);
	[t1]=findheight(GridDan(1).t+3,23.5);
	t2=62; %final dump
	
	xdat(1).x=GridDan(idir).VBAR;
    ydat(1).y = GridDan(1).Z/1000 + 0.62;
    
   labs(1).l='250m res';
   
	ylab='Height (km)';
    xlab=['Horizontal Wind Speed (m s^{-1})'];
    
    
    logflag=0;
    
    
     xlims=0;
	 lor=1;
     nmark=0;
     
     savename=['Horizontal wind speed profile'];
     titlenam=savename;
     
    case 60
    idir=1;
    
    izlim=1;
    zmin=0;
    zmax=18;
    
	for idat=1:3
        xdat(idat).x=width(1).w(:,idat+7)/1000;
        ydat(idat).y = GridDan(1).Z/1000 + 0.62;
        
        t1=GridDan(1).t(idat+7)+3;
        
        mins=(t1-floor(t1))*60;
        minstr=num2str(mins,'%2.0f');
    
        hrs=mod(floor(t1),24);
        hrstr=num2str(hrs,'%2.0f');
        if mins==0; minstr='00';end
        if hrs==0; hrstr='00';end
        labs(idat).l=[hrstr ':' minstr];
        
	end
    
    z2=findheight(GridDan(1).Z+620,2e3);
    xdat(4).x=interp1([2 12.5],[7 20],(GridDan(1).Z(z2:end)+620)/1000,'linear','extrap');
    labs(4).l='Linear variation'
    ydat(4).y = GridDan(1).Z(z2:end)/1000 + 0.62;
    
	ylab='Height (km)';
    xlab=['Updraught Width (km)'];
    
    
    logflag=0;
    
    
     xlims=0;
	 lor=4;
     nmark=0;
     
     savename=['Updraught width'];
     titlenam=savename;
     
    case 59
    idir=1;
    
    izlim=0;

	[izmin,izmax]=findheight(GridDan(1).Z+620,14e3,29e3);
	[t1]=findheight(GridDan(1).t+3,23.5);
	t2=62; %final dump
    
    dumprange=[t1:t2];

    ad_calcs4timeseries;
    
    pisub=f*sum(icediag(idir).i(izmin:izmax,dumprange,[24 25 27]),3);
    
	xdat(1).x=300*sum(icead(:,dumprange),2);
    xdat(2).x=300*sum(fallrate,2);
    xdat(3).x=300*sum(fallrate+icead(:,dumprange),2);
    xdat(4).x=300*sum(pisub,2);
    xdat(5).x=300*sum(microicerate+pisub,2);
    

    
	for idat=1:5
        ydat(idat).y = GridDan(1).Z(izmin:izmax)/1000 + 0.62;
	end
    
   labs(1).l='ad';
   labs(2).l='fall';
   labs(3).l='sum';
   labs(4).l='pisub';
   labs(5).l='microphysics';
   
	ylab='Height (km)';
    xlab=['Ice mixing ratio (ppmv)'];
    
    
    logflag=0;
    
    
     xlims=0;
	 lor=1;
     nmark=0;
     
     savename=['Sum of fall speed source of ice from 17-20 km and from 23.5-1.67 UTC'];
     savename=['Total ice from fall speed flux from 23.5-1.67 UTC'];     
     titlenam=savename;
     
    case 58
    idir=2;
    
    izlim=0;

	[iz,iz2]=findheight(GridDan(1).Z+620,14e3,20e3);
	[t1]=findheight(GridDan(1).t+3,23.5);
	t2=62; %final dump
	
	microicerate=f*sum(icediagsALL(idir).i(iz:iz2,t1:t2,31:33),3)/npes; %ice mixing ratio source rate
	xdat(1).x=300*sum(microicerate,2);

    
	for idat=1:1
        ydat(idat).y = GridDan(1).Z(iz:iz2)/1000 + 0.62;
	end
    
   labs(1).l='250m res';
   
	ylab='Height (km)';
    xlab=['Ice mixing ratio (ppmv)'];
    
    
    logflag=0;
    
    
     xlims=0;
	 lor=1;
     nmark=0;
     
     savename=['Sum of microphysical source of ice from 17-20 km and from 23.5-1.67 UTC'];
     savename=['Total ice formed microphysically from 23.5-1.67 UTC'];     
     titlenam=savename;


    case 57
    idir=1;
    
    izlim=0;
    iydir=-1; %reverse ydir as is pressure
%    zmin=15;
%    zmax=20;
    
   
	ylab='Pressure (hPa)';
    
    
    logflag=0;
    
    i55='total';
    i55='tp';
    i55='lwc';
    i55='lwc2';
    i55='inc';
    %i55='tracer';
  %  i55='tracer_2';
    
%     i55='therm_pos';
      i55='inc_emm';
      i55='ncw_emm';
%     i55='lwc_emm';
 %   i55='lwc_3';
%     i55='iwc_emm';
%i55='rwc';
%i55='lem_lwcad';
%i55='lem_Ndrops_warm_rain';
%i55='lem_Tprofile';
%i55='lem_microdiag';
%i55='lem_lwc';

    
     xlims=0;
	 lor=1;
     nmark=0;
    
     %run EMM_ACCpressure.m first

    
    switch i55
     case 'lem_lwc'
    	%	[tp,h]=plot_tephi(-40+273,50+273,153,5000);
        prate='Q02';

        dgarea='ACu';
       % dgarea='ALL';
		
        figname=['LWC contents vs Height for dump ' num2str(fnall)];
        xlab=['LWC (g kg^{-1})'];
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
         
    idirs=[1 3];
    groundheights=[1000 1000 1000 1000 620]; %%%% NOTE - make sure to set these properly %%%%%%%
     for idat=1:length(idirs)  
        
        idir=idirs(idat);
        groundheight=groundheights(idir);
		%groundheight=620;
		

        
        error=0;
        [dgfind,error]=getDGcol([dgarea '_A'],dgstrDan(idir).dg,error);
        area=TimeAvDan(idir).DGAV(:,dgfind(1));
        i0=find(area==0);
        area(i0)=1e99; %make zero areas=1e99 to avoid divide by zero (diags will be zero anyway)
        [dgfind,error]=getDGcol([dgarea '_' prate],dgstrDan(idir).dg,error);
        
        xdat(idat).x = 1000 * TimeAvDan(idir).DGAV(:,dgfind(1))./area*npess2(idir);
      %  xdat(idat).x = LWCacc*1000;
        ydat(idat).y = (GridDan(idir).Z+groundheight) / 1000;
        labs(idat).l = runName(idir).nam; 
            
     end
     
      for idat=[1:length(idirs)]+2 
        
        idir=idirs(idat-length(idirs));
        groundheight=groundheights(idir);
		
		nz=length(GridDan(idir).Z);
		zdat=GridDan(idir).Z(2:nz)+groundheight;
		pdat=GridDan(idir).PREFN(2:nz);
		qdat=GridDan(idir).OLQBAR(2:nz,1);
%		qdat(1)=qdat(2);
		tdat=tempLES(GridDan(idir));
		tdat=tdat(2:end);
		
        f=1e6*28.97/18;
		qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
		
        [CAPE,CIN,HLCL,TLCL,PLCL]=calc_cape(pdat,tdat,qdat,qsat,zdat);
        
        ipos=findheight(zdat,HLCL);     
        th_start=GridDan(idir).THREF(ipos);
        
		[tad,th_grid,p_grid]=moist_adiabat2(th_start,pdat(ipos),pdat(end)); %temperaure, potemp and pressure during moist saturated rise
        zdat2=interp1(pdat,zdat,p_grid); %find corresponding altitudes

            qsat2=satvappress(tad,'goff','liq',p_grid,1)/f; %the saturation mixing ratio during the adiabatic ascent - i.e. the vapour value
                                                            %of the parcel as are assuming it's at saturation all the way up
  
            labs(idat).l=runName(idir).nam; 
            xdat(idat).x=1000*(qsat(ipos)-qsat2); %the adiabatic LWC is the initial vapour value of the parcel at LCL minus  
            ydat(idat).y=zdat2/1000;                    %value maintained during ascent (=saturation mixing ratio)
            
     end
     
     
     
    case 'lem_microdiag'
        
    idirs=[1:5];
    groundheights=[1000 1000 1000 1000 620]; %%%% NOTE - make sure to set these properly %%%%%%%
     for idat=1:length(idirs)  
        
        idir=idirs(idat);
        groundheight=groundheights(idir);
		%groundheight=620;
		
	%	[tp,h]=plot_tephi(-40+273,50+273,153,5000);
        prate='PRAUT';
       % prate='PIHAL';
         prate='PRAUT';
%         prate='PGACW';
%         prate='PRAUT';
%         prate='PRACW';
%         prate='PIPRM';
     % prate='Q07';

        dgarea='ACu';
       % dgarea='ALL';
		
        figname=[prate ' vs Height for dump ' num2str(fnall)];
        xlab=['Micophysical Process Rate (kg kg^{-1} s^{-1})'];
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        error=0;
        [dgfind,error]=getDGcol([dgarea '_A'],dgstrDan(idir).dg,error);
        area=TimeAvDan(idir).DGAV(:,dgfind(1));
        i0=find(area==0);
        area(i0)=1e99; %make zero areas=1e99 to avoid divide by zero (diags will be zero anyway)
        [dgfind,error]=getDGcol([dgarea '_' prate],dgstrDan(idir).dg,error);
        
        xdat(idat).x = TimeAvDan(idir).DGAV(:,dgfind(1))./area*npess2(idir);
        ydat(idat).y = (GridDan(idir).Z+groundheight) / 1000;
        labs(idat).l = runName(idir).nam; 
            
     end
     
    case 'lem_Tprofile'
        
    idirs=[1 4 5];
    groundheights=[1000 1000 1000 1000 620]; %%%% NOTE - make sure to set these properly %%%%%%%
     for idat=1:length(idirs)  
        
        idir=idirs(idat);
        groundheight=groundheights(idir);
		%groundheight=620;
		
	%	[tp,h]=plot_tephi(-40+273,50+273,153,5000);
		
		nz=length(GridDan(idir).Z);
		zdat=GridDan(idir).Z(2:nz)+groundheight;
		tdat=tempLES(GridDan(idir));
		tdat=tdat(2:end);
		
        figname=['Temperature vs Height'];
        xlab=['Temperature (o^{C})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        xdat(idat).x = tdat-273.15;
        ydat(idat).y = zdat/1000;
        labs(idat).l = runName(idir).nam; 
            
     end
     

   
         
    case 'lem_Ndrops_warm_rain'
        
    idirs=[1 4 5];
    idirs=[1:5];
    groundheights=[1000 1000 1000 1000 620]; %%%% NOTE - make sure to set these properly %%%%%%%
     for idat=1:length(idirs)  
        
        idir=idirs(idat);
        groundheight=groundheights(idir);
		%groundheight=620;
		
	%	[tp,h]=plot_tephi(-40+273,50+273,153,5000);
		
		nz=length(GridDan(idir).Z);
		zdat=GridDan(idir).Z(2:nz)+groundheight;
		pdat=GridDan(idir).PREFN(2:nz);
		qdat=GridDan(idir).OLQBAR(2:nz,1);
%		qdat(1)=qdat(2);
		tdat=tempLES(GridDan(idir));
		tdat=tdat(2:end);
		
        f=1e6*28.97/18;
		qsat=satvappress(tdat,'goff','liq',pdat,1)/f;

        [CAPE,CIN,HLCL,TLCL,PLCL]=calc_cape(pdat,tdat,qdat,qsat,zdat);
        
        ipos=findheight(zdat,HLCL);     
        th_start=GridDan(idir).THREF(ipos);
        
		[tad,th_grid,p_grid]=moist_adiabat2(th_start,pdat(ipos),pdat(end));
        zdat2=interp1(pdat,zdat,p_grid);
        
        figname=['No. conc. vs Height'];
        xlab=['Number conc (cm^{-3})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        m=4/3*pi*(10e-6)^3*1000; %mass of a 20 micron diameter droplet

            qsat2=satvappress(tad,'goff','liq',p_grid,1)/f;
  
           labs(idat).l=runName(idir).nam; 
           lwc=qsat(ipos)-qsat2;
           
               [dgfind,error]=getDGcol([dgarea '_A'],dgstrDan(idir).dg,error);
                area=TimeAvDan(idir).DGAV(:,dgfind(1));
                i0=find(area==0);
                area(i0)=1e99; %make zero areas=1e99 to avoid divide by zero (diags will be zero anyway)
                [dgfind,error]=getDGcol([dgarea '_' prate],dgstrDan(idir).dg,error);
                
                lwc = TimeAvDan(idir).DGAV(:,dgfind(1))./area*npess2(idir);
                lwc = interp1(pdat,lwc(2:end),p_grid);
           
           rho=interp1(pdat,GridDan(idir).RHON(2:nz),p_grid);
           rhobase=GridDan(idir).RHON(ipos);           
           xdat(idat).x = lwc / m * 1e-6 * rhobase;
           ydat(idat).y=zdat2/1000;
            
         end
         
    case 'lem_lwcad'
        
    idirs=[1 4];
    groundheights=[1000 1000 1000 1000 620]; %%%% NOTE - make sure to set these properly %%%%%%%
     for idat=1:length(idirs)  
        
        idir=idirs(idat);
        groundheight=groundheights(idir);
		%groundheight=620;
		
	%	[tp,h]=plot_tephi(-40+273,50+273,153,5000);
		
		nz=length(GridDan(idir).Z);
		zdat=GridDan(idir).Z(2:nz)+groundheight;
		pdat=GridDan(idir).PREFN(2:nz);
		qdat=GridDan(idir).OLQBAR(2:nz,1);
%		qdat(1)=qdat(2);
		tdat=tempLES(GridDan(idir));
		tdat=tdat(2:end);
		
        f=1e6*28.97/18;
		qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
		
		%isup=find(qdat./qsat > 1);
		
		
        %qdat(isup)=qsat(isup);
		%qdat(isup(1))=qsat(isup(1))*0.98;
		
		
%		[CAPE,CIN,HLCL,TLCL,PLCL]=plot_tephi_data2(tdat,pdat,qdat,qsat,zdat);
        [CAPE,CIN,HLCL,TLCL,PLCL]=calc_cape(pdat,tdat,qdat,qsat,zdat);
        
        ipos=findheight(zdat,HLCL);     
        th_start=GridDan(idir).THREF(ipos);
        
		[tad,th_grid,p_grid]=moist_adiabat2(th_start,pdat(ipos),pdat(end));
        zdat2=interp1(pdat,zdat,p_grid);
        
        figname=['LWC vs Height'];
        xlab=['Mixing Ratio (g kg^{-1})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;

            qsat2=satvappress(tad,'goff','liq',p_grid,1)/f;
  
            labs(idat).l=runName(idir).nam; 
            xdat(idat).x=1000*(qsat(ipos)-qsat2);
            ydat(idat).y=zdat2/1000;
            
         end
        
    case 'rwc'
        figname=['RWC vs Height'];
        xlab=['Concentration (kg m^{-3})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        
  
    
%    xdat(idat+1).x = IWCacc;
%    ydat(idat+1).y = (GridDan(1).Z+620)/1000;         
%    labs(idat+1).l='LEM mean IWC'; 
        lor=4;             

        
        
        
    
    case 'total'
        figname=['Depleted total water points'];
        xlab=['Sum of deficit below 5 ppmv of points with tot water LT 5 ppmv (ppmv)'];
        xdat(1).x = dq_tot(idir).d(:,itdehyd,2)*length(GridDan(idir).Y1);
        labs(1).l='Total';
        
    case 'tp'
         ipos=25;
	th_start=Tacc(ipos)*(1e5/pacc(ipos))^0.286;
	adval=adiabat(end,2);
	[tad,th_grid,p_grid]=moist_adiabat2(th_start,pacc(ipos),adval);
     
    th_start=(14.57+273.15)*(1e5/745.9e2)^0.286;
	[tad,th_grid,p_grid]=moist_adiabat2(th_start,pacc(ipos),adiabat(end,2));
    
        figname=['Temperature vs Pressure'];
        xlab=['Temperature (K)'];
        
        xdat(1).x = Tacc;
        ydat(1).y = pacc;         
        labs(1).l='LEM'; 

        xdat(2).x = adiabat(:,3);
        ydat(2).y = adiabat(:,2);         
        labs(2).l='EMM'; 
    
        

        xdat(3).x = tad;
        ydat(3).y = p_grid;         
        labs(3).l=['Sat ad from LEM CB, index=' num2str(ipos)]; 

        
         xlims=1;
         xlimits=[200 300];
         
         lor=3;
        
        
    case 'lwc'
        figname=['LWC vs Pressure'];
        xlab=['Concentration (kg m^{-3})'];
        xdat(1).x = RAINacc;
        ydat(1).y = pacc;         
        labs(1).l='LEM LWC+rain'; 

        xdat(5).x = adiabat(:,4);
        ydat(5).y = adiabat(:,2);         
        labs(5).l='EMM'; 
        
        
        xdat(2).x=1e-3*meanselect(emmdat(1).lwc(:,1:60,1),'dat>1');
        xdat(2).x=1e-3*emmdat(1).lwc(:,20,1);
       % xdat(2).x = 1e-3*max(emmdat(1).lwc(:,1:60,1),[],2);
        ydat(2).y = interp1(rhenv(:,3),rhenv(:,2),Temm);         
        ydat(2).y = interp1(adiabat(:,3)-1,adiabat(:,2),Temm+273.15);  %is minus one since the temperature used
                % in the time-height plots is the adiabatic temperature minus tdiff1(=1). Here are interpolating
                % to get pressure of the temps in LWC time-height array from the adiabatic temp
        labs(2).l='EMM LWC'; 
        
        
        
%         ei=SatVapPress(adiabat(:,3),'goff','liq'); %Pa vapour pressure of EMM adiabat
%         sat=0.622*ei./(adiabat(:,2)-ei); %corresponding MR
%         rhoemm=pemm.*28.97e-3/8.3144./tad;
%         ih=findheight(adiabat(:,2),7.75e4); %pressure where LEM liquid starts to form in dump 7 

        ei=SatVapPress(tad,'goff','liq'); %Pa vapour pressure of EMM adiabat
        sat=0.622*ei./(p_grid-ei); %corresponding MR
        rhoemm=p_grid.*28.97e-3/8.3144./tad;
        ih=25; %pressure where LEM liquid starts to form in dump 7 
        
%        xdat(3).x = (sat(ih)-sat).*rhoemm; %LWC content assuming vapour input of saturation at this point
%        ydat(3).y = p_grid;         
%        labs(3).l='EMM new CB';   
        
        xdat(3).x = CONDacc; %LWC content assuming vapour input of saturation at this point
        ydat(3).y = pacc;         
        labs(3).l='LEM all condensate';   

         xdat(4).x = LWCacc;
        ydat(4).y = pacc;         
        labs(4).l='LEM LWC'; 
        
        lor=4;
    
    case 'lwc2'
        figname=['LWC vs Pressure'];
        xlab=['Concentration (kg m^{-3})'];
        xdat(1).x = RAINacc;
        ydat(1).y = pacc;         
        labs(1).l='LEM LWC+rain'; 

        xdat(5).x = adiabat(:,4);
        ydat(5).y = adiabat(:,2);         
        labs(5).l='EMM'; 
        
        iemm=11;
        
        xdat(2).x=1e-3*meanselect(emmdat(iemm).lwc(:,50:60,1),'dat>1');
        xdat(2).x=1e-3*mean(emmdat(iemm).lwc(:,1:50,1),2);
       % xdat(2).x = 1e-3*max(emmdat(1).lwc(:,1:60,1),[],2);
        ydat(2).y = interp1(rhenv(:,3),rhenv(:,2),Temm);         
        ydat(2).y = interp1(adiabat(:,3)-1,adiabat(:,2),Temm+273.15);  %is minus one since the temperature used
                % in the time-height plots is the adiabatic temperature minus tdiff1(=1). Here are interpolating
                % to get pressure of the temps in LWC time-height array from the adiabatic temp
        labs(2).l='EMM LWC'; 
        
        
        
%         ei=SatVapPress(adiabat(:,3),'goff','liq'); %Pa vapour pressure of EMM adiabat
%         sat=0.622*ei./(adiabat(:,2)-ei); %corresponding MR
%         rhoemm=pemm.*28.97e-3/8.3144./tad;
%         ih=findheight(adiabat(:,2),7.75e4); %pressure where LEM liquid starts to form in dump 7 

    %    ei=SatVapPress(tad,'goff','liq'); %Pa vapour pressure of EMM adiabat
%        sat=0.622*ei./(p_grid-ei); %corresponding MR
 %       rhoemm=p_grid.*28.97e-3/8.3144./tad;
 %       ih=25; %pressure where LEM liquid starts to form in dump 7 
        
%        xdat(3).x = (sat(ih)-sat).*rhoemm; %LWC content assuming vapour input of saturation at this point
%        ydat(3).y = p_grid;         
%        labs(3).l='EMM new CB';   
        
        xdat(3).x = CONDacc2-RAINacc2; %LWC content assuming vapour input of saturation at this point
        ydat(3).y = pacc2;         
        labs(3).l='LEM iwc';   

        xdat(4).x = LWCacc;
        ydat(4).y = pacc;         
        labs(4).l='LEM LWC'; 
        
        
        xdat(5).x =max(TwoD.Q(:,:,2),[],2).*GridDan(1).RHON;

        ydat(5).y = pacc;         
        labs(5).l='LEM max LWC'; 
                
        lor=4;
    
    case 'tracer'
        figname=['Tracer vs Height'];
        xlab=['Concentration (kg m^{-3})'];
    	ylab='Height (km)';

        
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=17;
        
        ixlim=1;
        xlimits=[0 1];
        
        
        
        xdat(1).x = TRACERacc; %tracer value taken from a dump with some criteria for e.g. cloudy points only
        ydat(1).y = (GridDan(1).Z+620)/1000;         
        labs(1).l='LEM tracer'; 
        
        itime=50; %time index for EMM profile to use
        izcb=findheight(GridDan(1).Z+620,vec(1).z(1)*1000); %find height index of EMM cloud base in LEM height grid
        sc=TRACERacc(izcb)/maxALL(emmdat(1).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base
        while isnan(sc)
            izcb=izcb+1;
            sc=TRACERacc(izcb)/maxALL(emmdat(1).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base
        end
           
        xdat(2).x = sc*emmdat(1).qup(:,[50],1);
        ydat(2).y = vec(1).z;         
        labs(2).l='EMM tracer'; 
        
%         
%         xdat(2).x=1e-3*meanselect(emmdat(1).lwc(:,1:60,1),'dat>1');
%         xdat(2).x=1e-3*emmdat(1).lwc(:,20,1);
%        % xdat(2).x = 1e-3*max(emmdat(1).lwc(:,1:60,1),[],2);
%         ydat(2).y = interp1(rhenv(:,3),rhenv(:,2),Temm);         
%         ydat(2).y = interp1(adiabat(:,3)-1,adiabat(:,2),Temm+273.15);  %is minus one since the temperature used
%                 % in the time-height plots is the adiabatic temperature minus tdiff1(=1). Here are interpolating
%                 % to get pressure of the temps in LWC time-height array from the adiabatic temp
%         labs(2).l='EMM LWC'; 
        
       lor=1;    
        
    case 'inc'
        figname=['INC vs Height'];
        xlab=['Number Concentration (m^{-3})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        
        xdat(1).x = INCacc;
        ydat(1).y = (GridDan(1).Z+620)/1000;         
        labs(1).l='LEM mean INC'; 
        
        xdat(1).x=mean(emmdat(2).inczt(:,30:50,1),2);
        ydat(1).y=vec(1).z;

%        xdat(5).x = adiabat(:,4);
%        ydat(5).y = adiabat(:,2);         
%        labs(5).l='EMM'; 
        
        
%        xdat(2).x=1e-3*meanselect(emmdat(1).lwc(:,1:60,1),'dat>1');
        xdat(2).x=mean(emmdat(1).inczt(:,30:50,1),2);
       % xdat(2).x = 1e-3*max(emmdat(1).lwc(:,1:60,1),[],2);
    %    ydat(2).y = interp1(rhenv(:,3),rhenv(:,2),Temm);         
%        ydat(2).y = interp1(adiabat(:,3)-1,adiabat(:,2),Temm+273.15);  %is minus one since the temperature used
                % in the time-height plots is the adiabatic temperature minus tdiff1(=1). Here are interpolating
                % to get pressure of the temps in LWC time-height array from the adiabatic temp
        ydat(2).y=vec(1).z;
        labs(2).l='EMM INC'; 
        
        xdat(3).x = GridDan(1).RHON.*max(sum(TwoD.Q(:,:,7:9),3),[],2);
%        xdat(3).x = INCmaxacc_alltim(1).dat(:,10);
        ydat(3).y = (GridDan(1).Z+620)/1000;
        labs(3).l='LEM max INC'; 
        
        xdat(3).x=mean(emmdat(3).inczt(:,30:50,1),2);
        ydat(3).y=vec(1).z;
        
        
        
%         ei=SatVapPress(adiabat(:,3),'goff','liq'); %Pa vapour pressure of EMM adiabat
%         sat=0.622*ei./(adiabat(:,2)-ei); %corresponding MR
%         rhoemm=pemm.*28.97e-3/8.3144./tad;
%         ih=findheight(adiabat(:,2),7.75e4); %pressure where LEM liquid starts to form in dump 7 

%        ei=SatVapPress(tad,'goff','liq'); %Pa vapour pressure of EMM adiabat
%        sat=0.622*ei./(p_grid-ei); %corresponding MR
%        rhoemm=p_grid.*28.97e-3/8.3144./tad;
%        ih=25; %pressure where LEM liquid starts to form in dump 7 
        
%        xdat(3).x = (sat(ih)-sat).*rhoemm; %LWC content assuming vapour input of saturation at this point
%        ydat(3).y = p_grid;         
%        labs(3).l='EMM new CB';   
        
%        xdat(3).x = CONDacc2-RAINacc2; %LWC content assuming vapour input of saturation at this point
%        ydat(3).y = pacc2;         
%        labs(3).l='LEM iwc';   

%         xdat(4).x = LWCacc;
%        ydat(4).y = pacc;         
%        labs(4).l='LEM LWC'; 
        
        lor=4;  

        
    case 'inc_emm'
        max_mean='Max';
        max_mean='Mean';
        
        tinds=[40:50];
        figname=[max_mean ' INC vs Height for ' num2str(tinds(1)) ' to ' num2str(tinds(end))];
        savename=figname;
%        xlab=['Number Concentration (m^{-3})'];
        xlab=['Number Concentration (kg^{-1})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
                
       for idat=1:length(vec)   
			M = 28*1.67E-27;
			k = 1.38E-23;
			G = 9.81;        
            
			p = emmdat(idat).adiabat(:,2);
			t = emmdat(idat).adiabat(:,3);
			rho=p.*M./k./t;   
            rho2=interp1(emmdat(idat).adiabat(:,1),rho,vec(idir).z*1000); %interpolate onto the vec grid
            
            labs(idat).l=run_name_emm{idat}; 
            switch max_mean
            case 'Mean'                
                xdat(idat).x=mean(emmdat(idat).inczt(:,tinds,1),2)./rho2;
            case 'Max'
                xdat(idat).x=max(emmdat(idat).inczt(:,tinds,1),[],2)./rho2;            
            end
            ydat(idat).y=vec(1).z;
        end
        
	%     xdat(idat+1).x = INCacc;
	%     ydat(idat+1).y = (GridDan(1).Z+620)/1000;         
	%     labs(idat+1).l='LEM mean INC'; 
            lor=4;         
        
    case 'therm_pos'
        figname=['Thermal position vs Height'];
        xlab=['Change in position (km)'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        
        xdat(1).x = thermalwidth(:,3);
        ydat(1).y = thermalwidth(:,2);         
        labs(1).l='EMM thermal'; 

     case 'ncw_emm'
        max_mean='Max';
        max_mean='Mean';
        
        tinds=[40:50];
        figname=[max_mean ' NCW vs Height ' num2str(tinds(1)) ' to ' num2str(tinds(end))];
        
        xlab=['Number Concentration (m^{-3})'];
        xlab=['Number Concentration (kg^{-1})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
    M = 28*1.67E-27;
	k = 1.38E-23;
	G = 9.81;        
        
   for idat=1:length(emmdat) 
		p = emmdat(idat).adiabat(:,2);
		t = emmdat(idat).adiabat(:,3);
		rho=p.*M./k./t;   
        rho2=interp1(emmdat(idat).adiabat(:,1),rho,vec(idir).z*1000); %interpolate onto the vec grid
        rho2=1;
        labs(idat).l=run_name_emm{idat};
        switch max_mean
        case 'Max'
            xdat(idat).x=max(emmdat(idat).ncw(:,30:50,1),[],2)./rho2;
        case 'Mean'
            xdat(idat).x=mean(emmdat(idat).ncw(:,tinds,1),2)./rho2; %convert to #/kg by division by rho
        end
        ydat(idat).y=vec(1).z;
    end
        
        lor=1;         
        
        
           
     case 'lwc_emm'
        figname=['LWC vs Height'];
        xlab=['Number Concentration (m^{-3})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        
   for idat=1:3    
        labs(idat).l=run_name_emm{idat}; 
        xdat(idat).x=1e-3*mean(emmdat(idat).lwc(:,30:50,1),2);
        ydat(idat).y=vec(1).z;
    end
        
        lor=1;  
        
    case 'lwc_3'
        
        figname=['LWC vs Height'];
        xlab=['Concentration (kg m^{-3})'];
        

        idats=[2 10 11];
         for i=1:length(idats)    
            idat=idats(i);
            labs(i).l=[run_name_emm{idat} 'LWC']; 
            xdat(i).x=1e-3*mean(emmdat(idat).lwc(:,30:40,1),2);
            ydat(i).y=vec(idat).z;
        end
        
        xdat(i+1).x = LWCacc;
        ydat(i+1).y = GridDan(1).Z/1000 +0.62;         
        labs(i+1).l='LEM LWC';   
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        
        
    case 'tracer_2'
        
        figname=['tracer vs height'];
        xlab=['Concentration (kg m^{-3})'];
        
        itime=50;
        itime=44;
    
    ivec=10;
        
        izcb=findheight(GridDan(1).Z+620,vec(ivec).z(1)*1000); %find height index of EMM cloud base in LEM height grid
        
        maxTracer=max(TwoD.Q(:,:,10),[],2);
        
        idats=[2 10 11];
        
         for i=1:length(idats)
            idat=idats(i);
       %     sc=TRACERacc(izcb)/maxALL(emmdat(idat).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base
%            while isnan(sc)
%                izcb=izcb+1;
%                sc=TRACERacc(izcb)/maxALL(emmdat(idat).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base
%            end
            

%            sc=maxTracer(izcb)/maxALL(emmdat(idat).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base
            
%            sc=TRACERacc3(1).dat(izcb,10)/maxALL(emmdat(idat).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base

            sc=TRACERmax(1).dat(izcb,10)/maxALL(emmdat(idat).qup(:,itime,1)); %scale factor so that EMM tracer equals LEM at EMM cloud base
            
        %    xdat(i).x = sc*mean(emmdat(idat).qup(:,[50:60],1),2);
            
            xdat(i).x = sc*mean(emmdat(idat).qup(:,[40:44],1),2);
        
            labs(i).l=[run_name_emm{idat}]; 
            ydat(i).y=vec(idats(i)).z;
        end
%         
%        xdat(i+1).x = TRACERacc3(1).dat(:,10); %tracer value taken from a dump with some criteria for e.g. cloudy points only
%        ydat(i+1).y = (GridDan(1).Z+620)/1000;         
%        labs(i+1).l='LEM tracer'; 
        
    for idat2=5:12
%        xdat(idat2-4+i).x = TRACERacc3(1).dat(:,idat2); %tracer value taken from a dump with some criteria for e.g. cloudy points only
        xdat(idat2-4+i).x = TRACERmax(1).dat(:,idat2); %tracer value taken from a dump with some criteria for e.g. cloudy points only        
        ydat(idat2-4+i).y = (GridDan(1).Z+620)/1000;         
        labs(idat2-4+i).l=num2str(idat2); 
    end

%         xdat(idat+1).x = maxTracer %tracer value taken from a dump with some criteria for e.g. cloudy points only
%         ydat(idat+1).y = (GridDan(1).Z+620)/1000;         
%         labs(idat+1).l='LEM max tracer'; 
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
    case 'iwc_emm'
        figname=['IWC vs Height'];
        xlab=['Concentration (kg m^{-3})'];
        
        ylab='Height (km)';
        iydir=1; %normal direction as using height
        
        izlim=1;
        zmin=0;
        zmax=20;
        
        
   for idat=1:length(vec)    
        labs(idat).l=run_name_emm{idat}; 
        xdat(idat).x=1e-3*mean(emmdat(idat).iwczt(:,50:60,1),2);
        ydat(idat).y=vec(1).z;
    end
    
    xdat(idat+1).x = IWCacc;
    ydat(idat+1).y = (GridDan(1).Z+620)/1000;         
    labs(idat+1).l='LEM mean IWC'; 
        lor=4;             

        
        
        
    end

    
%    xdat(2).x = f*sum(icediagsALL(3).i(:,1,[37:42]),3)/npess2(3);    
  %  xdat(2).x = sat;

    
	%	labs(2).l='Vapour mixing ratio (w/out alteration)';
    %    labs(2).l=['Ice saturation mixing ratio'];
	savename=[figname];
    titlenam=figname;
        
         
   
         
    case 56
    idir=1;
    
    zmin=15;
    zmax=20;
    
   
	ylab='Height (km)';
    
    
    logflag=0;
    
    i55='total';
    i55='nn';
    i55='mean';
  %  i55='ascent_rate';
    
     xlims=0;
	 lor=1;
    
    switch i55
    case 'total'
        figname=['Depleted total water points'];
        xlab=['Sum of deficit below 5 ppmv of points with tot water LT 5 ppmv (ppmv)'];
        xdat(1).x = dq_tot(idir).d(:,itdehyd,2)*length(GridDan(idir).Y1);
        labs(1).l='Total';
    case 'nn'
        figname=['Depleted total water points'];
        xlab=['Distance covered by points w/ tot water LT 5 ppmv (km)'];
        xdat(1).x = nn(idir).n(:,itdehyd,2)*(GridDan(1).Y1(2)-GridDan(1).Y1(1))/1000;
        labs(1).l='Distance';
    case 'mean'
        figname=['Simple lifting model vapour mixing ratio after ' num2str(tend*30) ' days'];
        xlab=['Mixing Ratio (ppmv)'];
        for iq=1:13
            xdat(iq).x = qq(:,iq);
            labs(iq).l=[num2str(Nevs(iq))]; 
            ydat(iq).y = zzf; 
            
        end
         xlims=1;
%         xlimits=[5.5 6.7];
         xlimits=[4.75 5.04];
         
         lor=2;
         
    case 'ascent_rate'     
        figname=['Simple Lifting Model Ascent Rate'];
        xlab=['Ascent Rate (km month^{-1})'];
        xdat(1).x = ww;
        labs(1).l=['Simple model']; 
        ydat(1).y = zzf; 
        
        zw(1)=14;
		wref(1)=-0.08;
		zw(2)=16.6;
		wref(2)=0.22;
		zw(3)=16.8;
		wref(3)=0.1; %0.1
		zw(4)=17;
		wref(4)=0.05;
		zw(5)=19;
		wref(5)=-0.01;
        
        wref=wref*1e-5 * 3600*24*30 * wfactor; %convert from cm/s to km/month
        
        xdat(2).x = wref;
        labs(2).l=['Original profile']; 
        ydat(2).y = zw; 
            
         xlims=0;
         xlimits=[4.4 5.2];
         
    end
    
%    xdat(2).x = f*sum(icediagsALL(3).i(:,1,[37:42]),3)/npess2(3);    
  %  xdat(2).x = sat;

    
	%	labs(2).l='Vapour mixing ratio (w/out alteration)';
    %    labs(2).l=['Ice saturation mixing ratio'];
	savename=[figname '_t=' num2str(itdehyd)];
    titlenam=figname;
        
         
         nmark=0;
         
   case 55
    idir=1;
    
    zmin=14;
    zmax=20;
    
   
	ylab='Height (km)';
    
    savename=figname;
    
    logflag=0;
    
    i55='total';
    i55='nn';
    i55='mean';
    
     xlims=0;
	 lor=1;
    
    switch i55
    case 'total'
        figname=['Depleted total water points'];
        xlab=['Sum of deficit below 5 ppmv of points with tot water LT 5 ppmv (ppmv)'];
        xdat(1).x = dq_tot(idir).d(:,itdehyd,2)*length(GridDan(idir).Y1);
        labs(1).l='Total';
    case 'nn'
        figname=['Depleted total water points'];
        xlab=['Distance covered by points w/ tot water LT 5 ppmv (km)'];
        xdat(1).x = nn(idir).n(:,itdehyd,2)*(GridDan(1).Y1(2)-GridDan(1).Y1(1))/1000;
        labs(1).l='Distance';
    case 'mean'
        figname=['Depleted total water points'];
        xlab=['Mean of points w/ tot water LT 5 ppmv (ppmv)'];
        xdat(1).x = 5 - dq_tot(idir).d(:,itdehyd,2)*length(GridDan(idir).Y1)./nn(idir).n(:,itdehyd,2);
        labs(1).l='Mean';
         xlims=1;
         xlimits=[3.6 5.2];
         lor=2;
    end
    
%    xdat(2).x = f*sum(icediagsALL(3).i(:,1,[37:42]),3)/npess2(3);    
  %  xdat(2).x = sat;

    
	for idat=1:1
        ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
	end
	
		
	%	labs(2).l='Vapour mixing ratio (w/out alteration)';
    %    labs(2).l=['Ice saturation mixing ratio'];
	
        
         
         nmark=0;
         

         
   case 54
    idir=1;
    
    zmin=6;
    zmax=19;
    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
    
    
    
    
    figname=['Initial vapour and ice saturation mixing ratios'];
    savename=figname;
    
    logflag=1;
    
    T=tempLES(GridDan(idir)); %K
	ei=SatVapPress(T,'goff','ice'); %Pa
	P=GridDan(idir).PREFN; %Pa
	sat=f*0.622*ei./(P-ei);

    xdat(1).x = f*sum(icediagsALL(1).i(:,1,[37:42]),3)/npess2(idir);
%    xdat(2).x = f*sum(icediagsALL(3).i(:,1,[37:42]),3)/npess2(3);    
    xdat(2).x = sat;

    
	for idat=1:2
        ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
	end
	
		labs(1).l='Vapour mixing ratio (model input)';
	%	labs(2).l='Vapour mixing ratio (w/out alteration)';
        labs(2).l=['Ice saturation mixing ratio'];
	
         xlims=1;
         xlimits=([4.8 5.8]);
         xlimits=([4.5 1000]);
         
         nmark=0;
         
         lor=1;
         
   case 53
    idir=1;
    
    zmin=0;
    zmax=3;
    
    xlab=['Potential Temperature (K)'];
	ylab='Height (km)';
       
    figname=['Potential tempertaure percentiles (each step represents one grid point)'];

    savename=figname;
    
    logflag=0;
    
    pr=100/length(GridDan(1).Y1); %min prctile for one gridpoint
    npr=5;
    
    thref=repmat(GridDan(1).THREF,[1 length(GridDan(1).Y1)]);
    th=thref+TwoD.TH1;
    thav=mean(thref+TwoD.TH1,2);
    
    prcs=[50 100:-pr:100-npr*pr];
    for iprc=1:length(prcs)
        xdat(iprc).x = prctile(th(2:end,:)',prcs(iprc));
        ydat(iprc).y = GridDan(idir).Z(2:end)/1000-0.62 + add_ground_height; 
     %   labs(iprc).l=[num2str( (100-prcs(iprc))/100*(GridDan(1).Y1(end)-GridDan(1).Y1(1) )/1000, '%3.2f'  ) ' km'];
%        labs(iprc).l=[num2str( prcs(iprc) ) ' th percentile'];
        labs(iprc).l=num2str(iprc-1);
    end
    labs(1).l='median';

	
         xlims=0;
         xlimits=([300 332]);
             
         nmark=0;
         
         lor=1;
         
         titlenam=figname;
         
   case 52
    zmin=10;
    zmax=22;
    
    xlab=['Flux of cloud ice (kg m^{-2} s^{-1})'];
	ylab='Height (km)';
    
    figname=['Flux of cloud ice (kg m^{-2} s^{-1})'];
    savename=figname;
    
    logflag=0;
    
	for idat=1:4
        ydat(idat).y = GridDan(1).Z/1000 + 0.62;
	end

    labs(1).l=runName(1).nam;
    labs(2).l=runName(2).nam;
    labs(3).l='Control ice, CCN=960 cm^{-3} updraught';
    labs(4).l='Control updraught, CCN=960 cm^{-3} ice';
    
    xdat(1).x = mean(icediagsALL(1).i(:,dumprange,87),2).*mean(icediagsALL(1).i(:,dumprange,137),2).*GridDan(1).RHON(:)/npess2(1),GridDan(1).Z(1:200)+620
    xdat(2).x = mean(icediagsALL(2).i(:,dumprange,87),2).*mean(icediagsALL(2).i(:,dumprange,137),2).*GridDan(1).RHON(:)/npess2(1),GridDan(1).Z(1:200)+620
    xdat(3).x = mean(icediagsALL(1).i(:,dumprange,87),2).*mean(icediagsALL(2).i(:,dumprange,137),2).*GridDan(1).RHON(:)/npess2(1),GridDan(1).Z(1:200)+620
    xdat(4).x = mean(icediagsALL(2).i(:,dumprange,87),2).*mean(icediagsALL(1).i(:,dumprange,137),2).*GridDan(1).RHON(:)/npess2(1),GridDan(1).Z(1:200)+620

         xlims=0;
         xlimits=([-1e7 1e7]);
         
         nmark=0;
         
         lor=2;
         
   case 51
    zmin=14.5;
    zmax=22;
    
    xlab=['Change in Mixing Ratio (ppmv)'];
	ylab='Height (km)';
    
    figname=['Change in Mixing Ratios'];
    savename=figname;
    
    logflag=0;
    
	for idat=1:2
        ydat(idat).y = GridDan(1).Z/1000 + 0.62;
	end

    labs(1).l='Change in Vapour';
    labs(2).l='Change in Ice';
    
    init=repmat(sum(icediagsALL(idir).i(:,1,[40:42]),3),[1 length(dumprange)]);
    changeice=f*( sum(icediagsALL(idir).i(:,dumprange,[40:42]),3) - init )/npes;
    
    init=repmat(sum(icediagsALL(idir).i(:,1,[37]),3),[1 length(dumprange)]);
    changevap=f*( sum(icediagsALL(idir).i(:,dumprange,[37]),3) - init )/npes;
    

    xdat(1).x = changevap(:,62);
    xdat(2).x = changeice(:,62);
    
%		labs(1).l='Time = 01:10 UTC';
%		labs(1).l='CCN = 240 cm^{-3}';
%		labs(2).l='CCN = 960 cm^{-3}';
	
         xlims=0;
         xlimits=([-1e7 1e7]);
         
         nmark=0;
         
         lor=2;
         
    case 50
        
    iz=132; %16.15km
    
    xlab=['Distance (km)'];
	ylab='LNB (km)';
    
    figname=['LNB_16.15km'];
    savename=figname;
    
    logflag=0;
    
	for idat=1:1
        xdat(idat).x = GridDan(idat).Y1/1000;
        ydat(idat).y = lnb2d(iz,:);
        labs(idat).l=runName(idat).nam;
	end
	
%		labs(1).l='Time = 01:10 UTC';
%		labs(1).l='CCN = 240 cm^{-3}';
%		labs(2).l='CCN = 960 cm^{-3}';
	
         xlims=0;
         xlimits=([-1e7 1e7]);
         
         nmark=0;
         
         lor=2;
         
    case 49
    t1=25.1667;
    t1=22.5;
    it=findheight(GridDan(idir).t+3,t1);
        
    zmin=14;
    zmax=22;
    
    xlab=['Total Water Mixing Ratio (ppmv)'];
	ylab='Height (km)';
    
    figname=['TotalWater'];
    savename=figname;
    
    logflag=0;
    
   
	for idat=1:4
        idir=idat;
        ad_calcs4timeseries;
        xdat(idat).x = f*cumsum(topdown(izmin:izmax,1:it),2)/npess2(idir);
        ydat(idat).y = GridDan(idir).Z(izmin:izmax)/1000 + add_ground_height; 
        labs(idat).l=runName(idat).nam;
	end
	
%		labs(1).l='Time = 01:10 UTC';
%		labs(1).l='CCN = 240 cm^{-3}';
%		labs(2).l='CCN = 960 cm^{-3}';
	
         xlims=1;
         xlimits=([-1e7 1e7]);
         
         nmark=0;
         
         lor=2;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 48
    idir=3;
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
    logflag=0;
    
	
    
%     xdat(1).x = f*sum(icediagsALL(1).i(:,1,[37:42]),3)/npess2(1);
% %     xdat(2).x = f*sum(icediagsALL(1).i(:,it1,[37:42]),3)/npess2(idir);   %f*sum(icediagsALL(idir+1).i(:,it,[37:42]),3)/npess2(idir+1);
% %      xdat(3).x = f*sum(icediagsALL(1).i(:,it2,[37:42]),3)/npess2(idir);
%      
%       xdat(2).x = f*sum(icediagsALL(1).i(:,it1,[37:42]),3)/npess2(1);
%       xdat(3).x = f*sum(icediagsALL(2).i(:,it1,[37:42]),3)/npess2(2);
%       xdat(4).x = f*sum(icediagsALL(3).i(:,it1,[37:42]),3)/npess2(3);
% 
%     
% 	
% 	
%     
%     
% 	
% 		labs(1).l='Initial';
%         
% %         labs(2).l='Final';
% %         labs(2).l=t1str;        
% % 
% %         labs(3).l='Final CCN 960 cm^{-3}';
% %         labs(3).l=t2str;
%         
%         
%         labs(2).l=[t1str ' 250m res tot'];
% 		labs(3).l=[t1str ' 500m res tot'];
% 		labs(4).l=[t1str ' 1km res tot'];
%         
%         
%      xdat(5).x = f*sum(icediagsALL(1).i(:,it1,[37]),3)/npess2(1);
%       xdat(6).x = f*sum(icediagsALL(2).i(:,it1,[37]),3)/npess2(2);
%       xdat(7).x = f*sum(icediagsALL(3).i(:,it1,[37]),3)/npess2(3);
% 
%     
% 
%     
%      
%         labs(5).l=[t1str ' 250m res vap'];
% 		labs(6).l=[t1str ' 500m res vap'];
% 		labs(7).l=[t1str ' 1km res vap'];
    

      xdat(1).x = f*sum(icediagsALL(1).i(:,1,[37:42]),3)/npess2(1); %initial total water 
      
      xdat(2).x = f*sum(icediagsALL(1).i(:,it1,[37:42]),3)/npess2(1);   %total water for 1st case
      xdat(3).x = f*sum(icediagsALL(idir).i(:,it1,[37:42]),3)/npess2(idir);

    
	
	
    
    
	
		labs(1).l='Initial';
        labs(2).l='1 km tot';  %runName(1).nam;    %[t1str ' Control total water'];
		labs(3).l='2 km tot';  %runName(idir).nam;   %[t1str ' CCN 960 cm^{-3} tot'];

        labs(2).l='1 km tot';  %runName(1).nam;    %[t1str ' Control total water'];
		labs(3).l=[runName(idir).nam ' tot'];   %[t1str ' CCN 960 cm^{-3} tot'];

        
     xdat(4).x = f*sum(icediagsALL(1).i(:,it1,[37]),3)/npess2(1);
     xdat(5).x = f*sum(icediagsALL(idir).i(:,it1,[37]),3)/npess2(idir);

        labs(4).l='1 km vap'; %runName(1).nam;    %[t1str ' Control total water'];
		labs(5).l='2 km vap'; %runName(idir).nam;   %[t1str ' CCN 960 cm^{-3} tot'];
        
        labs(4).l='1 km vap';  %runName(1).nam;    %[t1str ' Control total water'];
		labs(5).l=[runName(idir).nam ' vap'];   %[t1str ' CCN 960 cm^{-3} tot'];
        
    
    T=tempLES(GridDan(idir)); %K
	ei=SatVapPress(T,'goff','ice'); %Pa
	P=GridDan(idir).PREFN; %Pa
	sat=f*0.622*ei./(P-ei);
    
  %  xdat(6).x = sat;
  %  labs(6).l=['Ice sat MR'];
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    
            
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
            
            rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air1=flipud(air1);
            air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air2=flipud(air2);
   
            i0=find(xdat(idat).x>6.5);            
%            i0=find(xdat(idat).x>8.5);            

%            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
            
             
%            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
            tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
            tot(idat).t=flipud(tot(idat).t);
            
             
		end
        
        
        
         xlims=0;
         xlimits=([4.8 5.4]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
       %  xlimits=([0 20]);

         %xlimits=[-0.1 0.4];
         

icum=1;  %flag to say to do normal means and not TTL cumulative means.
if icum==1
         
         
            clear xdat ydat labs
                  %  zmin=15;
                  %5  zmax=18.7+0.62;
                  
                  Y1=diff(GridDan(1).Y1([1 end]));                  
                  Y2=diff(GridDan(2).Y1([1 end]));
                  
                  %only for 3d
                  fact=Y2*1000/Y1^2; %conversion factor to multiply 2-d result by for fair comparison (assumes a 1 km length in 3-D)
                  
%                  fact=Y2*1000/Y1; %factor assuming that 2d covers full 300 km of 3d domain

                  
                xlims=0;  
                xdat(1).x=tot(2).t * diff( GridDan(1).Y1([1 end]) ) ./air1;  %tot
                xdat(1).x=tot(4).t * diff( GridDan(1).Y1([1 end]) ) ./air1;  %vapour

                
                %%%  2d
                xdat(2).x=tot(3).t * diff( GridDan(2).Y1([1 end]) ) ./air2;   %tot
                xdat(2).x=tot(5).t * diff( GridDan(2).Y1([1 end]) ) ./air2;   %vapour
                
                %%%% 3d
%                xdat(2).x=tot(3).t * diff( GridDan(2).Y1([1 end]) ) ./air2 * fact; %%%%  only use fact for 3d  %%%%%%%%

                xdat(1).x(end+1)=0;
                xdat(2).x(end+1)=0;
                
                ydat(1).y=GridDan(1).Z(ih1:ih2+1)/1000 + add_ground_height; 
                ydat(2).y=GridDan(2).Z(ih1:ih2+1)/1000 + add_ground_height; 
                
                labs(1).l = '1 km';
                labs(2).l = '2 km';                                
        
                labs(1).l = '3d';
                labs(2).l = '2d (scaled assuming 1 km in 3d)';                                

                labs(1).l = '1 km';
                labs(2).l = runName(2).nam;      
                
                
                xlab='Cumulative total water reduction (kg m^{-1})';
                xlab='Mean total water MR change from 17 km downwards (ppmv)';
                
                        figname=['Cumulative total water'];
                        savename=figname;

end                
    
% 	
%          
         nmark=0;
         
         lor=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 487
    idir=1;
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    
    
            
      
    D=150e3;
    
    iifeature=30; %i index of feature
    ijfeature=15; %j index of feature 
    
     xfeature=0e3; %i index of feature
     yfeature=-75e3; %j index of feature 
   

%    xfeature=-100e3; %i index of feature
%    yfeature=150e3; %j index of feature 
%    yfeature=-125e3; %j index of feature 

    %here y and x refer to as appear on wrap_slice plots (first index in ThreeD refers to y axis
    
    
    [yinds,xinds]=find_inds_feature(GridDan(1),xfeature,yfeature,D);

    
    [iav iav2]=findheight( GridDan(idir).Y1, GridDan(idir).Y1(1)+D/2 , GridDan(idir).Y1(end)-D/2  );
    
    [iav iav2]=findheight( GridDan(idir).Y1, GridDan(idir).Y1(1)+D/2 , GridDan(idir).Y1(end)-D/2  );
    
    inds=[2:iav iav2:length(GridDan(idir).Y1)-1];
    
    [sx sy]=size(ThreeDDan(1).Q);
    [sxB syB]=size(ThreeDDan(idir).Q);
    
    inds=[1:76];
  

    init = f*sum(icediagsALL(1).i(:,1,[37]),3)/npess2(1); %initial total water profile

	idat=1; %don't change this
%             


    idir=1;
      for iz=1:length(GridDan(1).Z)

          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
%          xinds=[152-36:152 1:36]; yinds=[152-36:152 1:36];
         % xinds=[1:152]; yinds=[1:152];
          
          xdat(idat).x(iz,1) = f*  squeeze( mean( mean( ThreeDDan(idir).Q(yinds(2:end-1),xinds(2:end-1),iz) ) )  )  - init(iz); 
      
      end
      labs(idat).l=[runName(idir).nam ' vap']; 
      idat=idat+1;

      
      idir=2;
        for iz=1:length(GridDan(1).Z)
          idir=2;
          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
          xdat(idat).x(iz,1) = f*  squeeze( mean( mean( ThreeDDan(idir).Q(yinds(2:end-1),xinds(2:end-1),iz) ) )  )  - init(iz); 
      
        end
    
%      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(1).Q(yinds,xinds,:) ) )  )  - init;     
      labs(idat).l=[runName(idir).nam ' vap']; 
      idat=idat+1;

      idir=3;
      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) ) )  )  - init; 
      labs(idat).l=[runName(idir).nam ' vap'];
      idat=idat+1;
      
  ice=1;  
  switch ice
  case 1
      idir=1;
      for iz=1:length(GridDan(idir).Z)
%          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz) + totice_44(idir).dat(2:end-1,2:end-1,iz));
          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
          xdat(idat).x(iz,1) = f*  squeeze( mean( mean( ThreeDDan(idir).Q(yinds(2:end-1),xinds(2:end-1),iz) + totice_44(idir).dat(yinds(2:end-1),xinds(2:end-1),iz) ) )  )  - init(iz); 
      
      end
      labs(idat).l=[runName(idir).nam ' tot']; 
      idat=idat+1;
      
      
      idir=2
      for iz=1:length(GridDan(idir).Z)
%          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz) + totice_44(idir).dat(2:end-1,2:end-1,iz));
          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
          xdat(idat).x(iz,1) = f*  squeeze( mean( mean( ThreeDDan(idir).Q(yinds(2:end-1),xinds(2:end-1),iz) + totice_44(idir).dat(yinds(2:end-1),xinds(2:end-1),iz) ) )  )  - init(iz); 
      
      end
      labs(idat).l=[runName(idir).nam ' tot']; 
      idat=idat+1;
      
      idir=3;
      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) + totice_44(idir).dat(2:end-1,2:end-1,:) ) )  )  - init; 
      labs(idat).l=[runName(idir).nam ' tot']; 
      idat=idat+1;
      
  end %switch ice
  
  
      
%       xdat(idat).x = f*  squeeze( mean( mean( ThreeD(1).Q ) )  )  - init; 
%       labs(idat).l=['3D vap']; 
%       idat=idat+1;
% %       
%       xdat(idat).x = f*  squeeze( mean( mean( ThreeD(1).Q ) )  )  - init; 
%       labs(idat).l=['3D vap']; 
%       idat=idat+1;
%       
%       xdat(idat).x = f*squeeze(mean( mean(tot_water44 ) )) - init; 
%       labs(idat).l=['3D tot']; 
%       idat=idat+1;
      
      
     
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
%    zmax=20.7+0.62;

    zmax=18+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    ih2=findheight(GridDan(1).Z+620,17.05e3);
    ih1=findheight(GridDan(1).Z+620,16.05e3);
  %  ih2=findheight(GridDan(1).Z+620,16.05e3);
  %  ih1=findheight(GridDan(1).Z+620,15.05e3);
    
            
%         for idat=1:length(xdat)
%             ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
%             
%             rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
%             dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
%             
%             air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
%             air1=flipud(air1);
%             air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
%             air2=flipud(air2);
%    
%             i0=find(xdat(idat).x>6.5);            
% %            i0=find(xdat(idat).x>8.5);            
% 
% %            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
%             
%              
% %            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
%             tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
%             tot(idat).t=flipud(tot(idat).t);
%             
%             tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
%             
%              
% 		end
        
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
            
            rho=GridDan(idir).RHON(ih1:ih2); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=mean( dz.*rho );
   
           % i0=find(xdat(idat).x>6.5);            

            tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
            
             
		end
        
        
        
         xlims=1;
         xlimits=([-0.5 4]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 4876 %general 3D vapour and total water plot comparison and averaging
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    idir=1;
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    
    
            
      
   

    init = f*sum(icediagsALL(1).i(:,1,[37]),3)/npess2(1); %initial total water profile

	idat=1; %don't change this
%             

    idirs=[1:4];
    Ldir=length(idirs);
    idat=1;
    
    vap=2;
    switch vap
    case 1
      
      for idir2=1:Ldir;
          idir=idirs(idir2);
          xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) ) )  )  - init; 
          labs(idat).l=[runName(idir).nam ' vap'];
          area(idat).dat = diff( GridDan(idir).X1([end 1]) ) * diff( GridDan(idir).Y1([end 1]) ); %area of domain in m^2
          idat=idat+1;          
      end
      
      case 2
      
      for idir2=1:Ldir;
          idir=idirs(idir2);
          xdat(idat).x = f*  squeeze( mean( mean( sum(ThreeDDan(idir).Q(2:end-1,2:end-1,:,1),4) ) )  )  - init; 
          labs(idat).l=[runName(idir).nam ' vap'];
          area(idat).dat = diff( GridDan(idir).X1([end 1]) ) * diff( GridDan(idir).Y1([end 1]) ); %area of domain in m^2
          idat=idat+1;          
      end
      
  end
      
  ice=2;  
  switch ice
  case 1
      
      for idir2=1:Ldir;                  
          idir=idirs(idir2);
          xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) + totice_44(idir).dat(2:end-1,2:end-1,:) ) )  )  - init; 
          labs(idat).l=[runName(idir).nam ' tot']; 
          area(idat).dat = diff( GridDan(idir).X1([end 1]) ) * diff( GridDan(idir).Y1([end 1]) ); %area of domain in m^2
          idat=idat+1;      
      end   
      
  case 2
      for idir2=1:Ldir;
          idir=idirs(idir2);
          xdat(idat).x = f*  squeeze( mean( mean( sum(ThreeDDan(idir).Q(2:end-1,2:end-1,:,1:4),4) ) )  )  - init; 
          labs(idat).l=[runName(idir).nam ' tot'];
          area(idat).dat = diff( GridDan(idir).X1([end 1]) ) * diff( GridDan(idir).Y1([end 1]) ); %area of domain in m^2
          idat=idat+1;          
      end
      
      
  end %switch ice                           
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
%    zmax=20.7+0.62;

    zmax=20+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
 %   ih2=findheight(GridDan(1).Z+620,17e3);
 %   ih1=findheight(GridDan(1).Z+620,15.05e3);
 %   ih1=findheight(GridDan(1).Z+620,14e3);
%    ih2=findheight(GridDan(1).Z+620,17.05e3);
    ih1=findheight(GridDan(1).Z+620,16.05e3);
%    ih2=findheight(GridDan(1).Z+620,16.05e3);
    ih1=findheight(GridDan(1).Z+620,15.9e3);
%    ih2=length(GridDan(1).Z); %index at top of the domain
    
 %   ih1=135; %index for 380 K level in Grid.THREF (actually 379.6 K)
    
            
%         for idat=1:length(xdat)
%             ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
%             
%             rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
%             dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
%             
%             air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
%             air1=flipud(air1);
%             air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
%             air2=flipud(air2);
%    
%             i0=find(xdat(idat).x>6.5);            
% %            i0=find(xdat(idat).x>8.5);            
% 
% %            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
%             
%              
% %            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
%             tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
%             tot(idat).t=flipud(tot(idat).t);
%             
%             tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
%             
%              
% 		end
        
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
            
            rho=GridDan(idir).RHON(ih1:ih2); 
            dz=diff(GridDan(idir).Z(ih1-1:ih2));
            
            air1=mean( dz.*rho );
   
            tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
            
            %total mass increase at each level (vapour or total water)
            mass(idat).dat=sum(xdat(idat).x(ih1:ih2)/f.*area(idat).dat.*dz.*rho);  %total mass increase in kg
            
             
		end
        
        
        
         xlims=1;
         xlimits=([-0.5 4]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 4874 %high low CCN 3D
    idir=1;
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    
    
            
      
    D=150e3;
    
    iifeature=30; %i index of feature
    ijfeature=15; %j index of feature 
    
     xfeature=0e3; %i index of feature
     yfeature=-75e3; %j index of feature 
   

%    xfeature=-100e3; %i index of feature
%    yfeature=150e3; %j index of feature 
%    yfeature=-125e3; %j index of feature 

    %here y and x refer to as appear on wrap_slice plots (first index in ThreeD refers to y axis
    
    
    [yinds,xinds]=find_inds_feature(GridDan(1),xfeature,yfeature,D);

    
    [iav iav2]=findheight( GridDan(idir).Y1, GridDan(idir).Y1(1)+D/2 , GridDan(idir).Y1(end)-D/2  );
    
    [iav iav2]=findheight( GridDan(idir).Y1, GridDan(idir).Y1(1)+D/2 , GridDan(idir).Y1(end)-D/2  );
    
    inds=[2:iav iav2:length(GridDan(idir).Y1)-1];
    
    [sx sy]=size(ThreeDDan(1).Q);
    [sxB syB]=size(ThreeDDan(idir).Q);
    
    inds=[1:76];
  

    init = f*sum(icediagsALL(1).i(:,1,[37]),3)/npess2(1); %initial total water profile

	idat=1; %don't change this
%             

    idir=1;
      for iz=1:length(GridDan(1).Z)

          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
%          xinds=[152-36:152 1:36]; yinds=[152-36:152 1:36];
         % xinds=[1:152]; yinds=[1:152];
         xinds=1:length(GridDan(idir).X1);
         yinds=1:length(GridDan(idir).Y1);
          
        xdat(idat).x(iz,1) = f*  squeeze( mean( mean( ThreeDDan(idir).Q(yinds(2:end-1),xinds(2:end-1),iz) ) )  )  - init(iz); 
      
      end
      labs(idat).l=[runName(idir).nam ' vap']; 
      idat=idat+1;     

      
      idir=2;
      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) ) )  )  - init; 
      labs(idat).l=[runName(idir).nam ' vap'];
      idat=idat+1;
      
      idir=3;
      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) ) )  )  - init; 
      labs(idat).l=[runName(idir).nam ' vap'];
      idat=idat+1;
      
  ice=1;  
  switch ice
  case 1
      idir=1;
      for iz=1:length(GridDan(idir).Z)
%          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz) + totice_44(idir).dat(2:end-1,2:end-1,iz));
          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
          xinds=1:length(GridDan(idir).X1);
          yinds=1:length(GridDan(idir).Y1);
          
          xdat(idat).x(iz,1) = f*  squeeze( mean( mean( ThreeDDan(idir).Q(yinds(2:end-1),xinds(2:end-1),iz) + totice_44(idir).dat(yinds(2:end-1),xinds(2:end-1),iz) ) )  )  - init(iz); 
      
      end
      labs(idat).l=[runName(idir).nam ' tot']; 
      idat=idat+1;
                  
      
      idir=2;
      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) + totice_44(idir).dat(2:end-1,2:end-1,:) ) )  )  - init; 
      labs(idat).l=[runName(idir).nam ' tot']; 
      idat=idat+1;
      
      idir=3;
      xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(idir).Q(2:end-1,2:end-1,:) + totice_44(idir).dat(2:end-1,2:end-1,:) ) )  )  - init; 
      labs(idat).l=[runName(idir).nam ' tot']; 
      idat=idat+1;
      
  end %switch ice
  
  
      
%       xdat(idat).x = f*  squeeze( mean( mean( ThreeD(1).Q ) )  )  - init; 
%       labs(idat).l=['3D vap']; 
%       idat=idat+1;
% %       
%       xdat(idat).x = f*  squeeze( mean( mean( ThreeD(1).Q ) )  )  - init; 
%       labs(idat).l=['3D vap']; 
%       idat=idat+1;
%       
%       xdat(idat).x = f*squeeze(mean( mean(tot_water44 ) )) - init; 
%       labs(idat).l=['3D tot']; 
%       idat=idat+1;
      
      
     
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
%    zmax=20.7+0.62;

    zmax=18+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    ih2=findheight(GridDan(1).Z+620,17.05e3);
    ih1=findheight(GridDan(1).Z+620,16.05e3);
    ih2=findheight(GridDan(1).Z+620,16.05e3);
    ih1=findheight(GridDan(1).Z+620,15.05e3);
    
            
%         for idat=1:length(xdat)
%             ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
%             
%             rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
%             dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
%             
%             air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
%             air1=flipud(air1);
%             air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
%             air2=flipud(air2);
%    
%             i0=find(xdat(idat).x>6.5);            
% %            i0=find(xdat(idat).x>8.5);            
% 
% %            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
%             
%              
% %            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
%             tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
%             tot(idat).t=flipud(tot(idat).t);
%             
%             tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
%             
%              
% 		end
        
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
            
            rho=GridDan(idir).RHON(ih1:ih2); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=mean( dz.*rho );
   
           % i0=find(xdat(idat).x>6.5);            

            tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
            
             
		end
        
        
        
         xlims=1;
         xlimits=([-0.5 4]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 4872
    idir=1;
       
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
   
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC    
    
    xlab=['Temperature change (K)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Temp change'];

    savename=figname;
    
 
    logflag=0;
    
    
        

	idat=1; %don't change this
    
%       xdat(idat).x = f*  squeeze( mean( mean( ThreeDDan(1).Q(yinds,xinds,:) ) )  )  - init; 
%       labs(idat).l=runName(1).nam; 
%       idat=idat+1;
%     
%    init=repmat( icediagsALL(1).i(:,1,246), [1 size(icediagsALL(1).i,2)]  );
    
P=GridDan(1).PREFN;
	init=icediagsALL(1).i(:,1,246) ./(1e5./P).^0.286 ;      
%     
%       xdat(idat).x = icediagsALL(1).i(:,44,246)./(1e5./P).^0.286 - init; 
%       labs(idat).l=[runName(1).nam 'icdiags']; 
%       idat=idat+1;
      
      
      
      
%       xdat(idat).x = squeeze(mean(mean(T,2),3)) - init; 
%       labs(idat).l=[runName(1).nam 'ThreeD']; 
%       idat=idat+1;


      
      D=300e3;
    %  D=150e3;
      idir=1;
      [T]=temp_from_press_and_th(GridDan(idir),ThreeDDan(idir).TH1,ThreeDDan(idir).P);
      for iz=1:length(GridDan(idir).Z)
%          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz) + totice_44(idir).dat(2:end-1,2:end-1,iz));
          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz));
          med=median(median( ThreeDDan(idir).Q(2:end-1,2:end-1,iz,1) ));
          
          
          xfeature=GridDan(idir).X1(mi(2));
          yfeature=GridDan(idir).Y1(mi(1));                    
          [yinds,xinds]=find_inds_feature(GridDan(idir),xfeature,yfeature,D);
          
         % xdat(idat).x(iz,1) = squeeze( mean( mean( T(iz,yinds(2:end-1)-1,xinds(2:end-1)-1)  ) )  )  - init(iz); 
          %xdat(idat).x(iz,1) = squeeze(  T( iz , mi(1),mi(2) ) - init(iz) ); 
          Tiz=T(iz,:,:);
          ivap=find( sum( ThreeDDan(idir).Q(2:end-1,2:end-1,iz,:) , 4) < med*1.05 ); %sum for total water
          xdat(idat).x(iz,1) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
      
      end
      labs(idat).l=[runName(idir).nam]; 
      idat=idat+1;            
      
%       idir=2;
%       [T]=temp_from_press_and_th(GridDan(idir),ThreeDDan(idir).TH1,ThreeDDan(idir).P);    
%       xdat(idat).x = squeeze( mean( mean( T,3 ),2 )  )  - init; 
%       labs(idat).l=[runName(idir).nam]; 
%       idat=idat+1;
%       
%       idir=3;
%       [T]=temp_from_press_and_th(GridDan(idir),ThreeDDan(idir).TH1,ThreeDDan(idir).P);    
%       xdat(idat).x = squeeze( mean( mean( T,3 ),2 )  )  - init; 
%       labs(idat).l=[runName(idir).nam]; 
%       idat=idat+1;
%       
%             aind=107; %ACC_A 
% 			area=TimeAv.DGAV(:,107);
%             area(area==0)=1;
%             
%             
%       init=icediagsALL(idir).i(:,1,383);  
%       idir=1;
%       xdat(idat).x = icediagsALL(idir).i(:,10,384)./area - init; 
%       %xdat(idat).x=area;
%       labs(idat).l=[runName(idir).nam]; 
%       idat=idat+1;

%       
%       xdat(idat).x = xdat(2).x-xdat(1).x; 
%       labs(idat).l='diff'; 
%       idat=idat+1;
%         
    clear diff tot
    
    zmin=15;
    zmin=14;
    zmin=10;
    zmax=16.7+0.62;
    zmax=20.7+0.62;

    %zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17.1e3);
%    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,16.1e3);
    
            
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
            
            rho=GridDan(idir).RHON(ih1:ih2); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=mean( dz.*rho );
   
           % i0=find(xdat(idat).x>6.5);            

            tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
            
             
		end
        
        
        
         xlims=0;
         xlimits=([-1 5]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 4875
    idir=1;        
    

	ylab='Height (km)';
      %  ylab='';        
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];

    savename=figname;
    
 
    logflag=0;        
        

	idat=1; %don't change this    
           
    idirs=[1:3];
    
    mean_min='mean';
%    mean_min='min';
    mean_min='mean2';
%    mean_min='max';
%    mean_min='mean_region';

		dX=75e3;
        dY=diff(GridDan(idir).Y1(1:2));
        dz=diff(GridDan(idir).Z);
        dz=[dz(1); dz];
        dV=GridDan(idir).RHON.*dY*dY.*dz;
        
    
    
    xlab=['Mean vapour (ppmv)'];
         figname=['Mean vapour'];
         
         
    switch mean_min
    case 'min' 
        
         for idat=1:length(idirs)
                idir=idirs(idat);
                xdat(idat).x = f*squeeze(min(min(ThreeDDan(idir).Q(2:end-1,2:end-1,:),[],1),[],2)) ; 
                labs(idat).l=[runName(idir).nam]; 
         end
         
         xlab=['Min vapour (ppmv)'];
         figname=['Min vapour'];
      
    case 'mean'
        
        for idat=1:length(idirs)
            idir=idirs(idat);
            xdat(idat).x = f*squeeze(mean(mean(ThreeDDan(idir).Q(2:end-1,2:end-1,:),1),2)) ; 
            labs(idat).l=[runName(idir).nam]; 
        end
        
         xlab=['Mean vapour (ppmv)'];
         figname=['Mean vapour'];
         
     case 'mean2'
        mean_inds=[1:4]; %indices for q-fields to do means over
        dV=dY*dY
        
        for idat=1:length(idirs)
            idir=idirs(idat);
            xdat(idat).x = dV*squeeze(sum(sum(sum(f*ThreeDDan(idir).Q(2:end-1,2:end-1,:,mean_inds)-5,4),1),2)) /f ; 
            labs(idat).l=[runName(idir).nam]; 
            
            xlab=['Sum ice (kg)'];
            figname=['Sum ice'];
        end
        
        
	case 'max'
        mean_inds=[1:4]; %indices for q-fields to do means over
        for idat=1:length(idirs)
            idir=idirs(idat);
            xdat(idat).x = f*squeeze(max(max(max(ThreeDDan(idir).Q(2:end-1,2:end-1,:,mean_inds),[],4),[],1),[],2)) ; 
            labs(idat).l=[runName(idir).nam]; 
            
            xlab=['Max ice (ppmv)'];
            figname=['Max ice'];
        end   
        
     case 'mean_region'
        mean_inds=[1]; %indices for q-fields to do means over
        
        
        
        for idat=1:length(idirs)
            idir=idirs(idat);
            xdat(idat).x = f*squeeze(mean(mean(sum(ThreeDDan(idir).Q(inds,2:end-1,:,mean_inds),4),1),2)) -5 ; 
            labs(idat).l=[runName(idir).nam]; 
            
            xlab=['Mean ice region (ppmv)'];
            figname=['Mean ice region'];
        end   
         
        
    end
      
        
    clear diff tot
    
    zmin=15;
    zmin=14;
    zmax=16.7+0.62;
    zmax=20.7+0.62;
    zmin=0;

    %zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17.1e3);
%    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,16.1e3);
    
            
         for idat=1:length(xdat)
             ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
%             
%             rho=GridDan(idir).RHON(ih1:ih2); %convert to kg/km3 as xdat in g/kg km 
%             dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
%             
%             air1=mean( dz.*rho );
%    
%            % i0=find(xdat(idat).x>6.5);            
% 
%             tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
%             
%              
        end
        
        
        
         xlims=0;
         xlimits=([-1 5]);
         xlimits=([4.5 6.2]);
         xlimits=([4.5 7]);
       % xlimits=([4.8 5.4]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 4873
    idir=1;        
    

	ylab='Height (km)';
      %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];

    savename=figname;
    
 
    logflag=0;        
        

	idat=1; %don't change this    
           
    idirs=[1:3];
    
    mean_min='mean';
    mean_min='min';
    
    switch mean_min
    case 'min' 
        
         for idat=1:length(idirs)
                idir=idirs(idat);
                xdat(idat).x = f*squeeze(min(min(ThreeDDan(idir).Q(2:end-1,2:end-1,:),[],1),[],2)) ; 
                labs(idat).l=[runName(idir).nam]; 
         end
         
         xlab=['Min vapour (ppmv)'];
         figname=['Min vapour'];
      
    case 'mean'
        
        for idat=1:length(idirs)
            idir=idirs(idat);
            xdat(idat).x = f*squeeze(mean(mean(ThreeDDan(idir).Q(2:end-1,2:end-1,:),1),2)) ; 
            labs(idat).l=[runName(idir).nam]; 
        end
        
         xlab=['Mean vapour (ppmv)'];
         figname=['Mean vapour'];
        
    end
      
        
    clear diff tot
    
    zmin=15;
    zmin=14;
    zmax=16.7+0.62;
    zmax=20.7+0.62;

    %zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17.1e3);
%    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,16.1e3);
    
            
         for idat=1:length(xdat)
             ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
%             
%             rho=GridDan(idir).RHON(ih1:ih2); %convert to kg/km3 as xdat in g/kg km 
%             dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
%             
%             air1=mean( dz.*rho );
%    
%            % i0=find(xdat(idat).x>6.5);            
% 
%             tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
%             
%              
        end
        
        
        
         xlims=1;
         xlimits=([-1 5]);
         xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        xlimits=([4.8 5.4]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 485
    idir=1;
    
%need to load in 2D case, 3D case and then do loadvapdata for icediagsALL and GridDan    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    
    
    D=300e3;
    [iav iav2]=findheight( GridDan(2).Y1, GridDan(2).Y1(1)+D/2 , GridDan(2).Y1(end)-D/2  );	

      init = f*sum(icediagsALL(2).i(:,1,[37:42]),3)/npess2(2); %initial total water profile

        idat=1;	             
      
            
      xdat(idat).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1),3) ,2) - init; 
      labs(idat).l=['2D vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1:6),3) ,2) - init; 
      labs(idat).l=['2D tot']; 
      idat=idat+1;
      
      xdat(idat).x = 0.22 * (  f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1:6),3) ,2) - init  ); 
      labs(idat).l=['2D tot scaled']; 
      idat=idat+1;
      
      xdat(idat).x = f*  squeeze( mean( mean( ThreeD(1).Q ) )  )  - init; 
      labs(idat).l=['3D vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*squeeze(mean( mean(tot_water44 ) )) - init; 
      labs(idat).l=['3D tot']; 
      idat=idat+1;
      
      
     
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
    %zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    
            
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
            
            rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air1=flipud(air1);
            air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air2=flipud(air2);
   
            i0=find(xdat(idat).x>6.5);            
%            i0=find(xdat(idat).x>8.5);            

%            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
            
             
%            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
            tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
            tot(idat).t=flipud(tot(idat).t);
            
             
		end
        
        
        
         xlims=1;
         xlimits=([-1 5]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
 
 case 4852    

    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    
    
    D=50e3;
    [iav iav2]=findheight( GridDan(idir).Y1, GridDan(idir).Y1(idir)+D/2 , GridDan(idir).Y1(end)-D/2  );	

      init = f*sum(icediagsALL(1).i(:,2,[37:42]),3)/npess2(1); %initial total water profile

      idat=1;	                   
      
      idir=2;
      xdat(idat).x = f*mean( sum(TwoDDan(idir).Q(:,[1:iav iav2:length(GridDan(idir).Y1)] ,1),3) ,2) - init; 
      ydat(idat).y=(GridDan(idir).Z+620)/1000;
      labs(idat).l=['1e-3 vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(idir).Q(:,[1:iav iav2:length(GridDan(idir).Y1)] ,1:6),3) ,2) - init; 
      ydat(idat).y=(GridDan(idir).Z+620)/1000;      
      labs(idat).l=['1e-3 tot']; 
      idat=idat+1;
      
      idir=1;
      [iav iav2]=findheight( GridDan(idir).Y1, GridDan(idir).Y1(idir)+D/2 , GridDan(idir).Y1(end)-D/2  );	
      
      xdat(idat).x = f*mean( sum(TwoDDan(idir).Q(:,[1:iav iav2:length(GridDan(idir).Y1)] ,1),3) ,2) - init; 
      ydat(idat).y=(GridDan(idir).Z+620)/1000;
      labs(idat).l=['2D vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(idir).Q(:,[1:iav iav2:length(GridDan(idir).Y1)] ,1:6),3) ,2) - init; 
      ydat(idat).y=(GridDan(idir).Z+620)/1000;      
      labs(idat).l=['2D tot']; 
      idat=idat+1;
      
      izlim=1;
       zmin=15;
% %    zmin=14.5;
     zmax=16.7+0.62;
       
         xlims=1;
         xlimits=([-1 5]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         xlimits=[-1 0.1];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
 

    case 486
    idir=1;
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    D=300e3;
    [iav iav2]=findheight( GridDan(1).Y1, GridDan(1).Y1(1)+D/2 , GridDan(1).Y1(end)-D/2  );	

      init = f*sum(icediagsALL(2).i(:,1,[37:42]),3)/npess2(2); %initial total water 

        idat=1;
	
        
      xdat(idat).x = f*mean( sum(TwoDDan(1).Q(:,[1:iav iav2:length(GridDan(1).Y1)] ,1),3) ,2) - init; 
      labs(idat).l=['1 km res vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(1).Q(:,[1:iav iav2:length(GridDan(1).Y1)] ,1:6),3) ,2) - init; 
      labs(idat).l=['1 km res tot']; 
      idat=idat+1;
      
       [iav iav2]=findheight( GridDan(2).Y1, GridDan(2).Y1(1)+D/2 , GridDan(2).Y1(end)-D/2  );	
      
      xdat(idat).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1),3) ,2) - init; 
      labs(idat).l=['2 km res vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1:6),3) ,2) - init; 
      labs(idat).l=['2 km res tot']; 
      idat=idat+1;
      
      
     
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
    %zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    
            
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
            
            rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air1=flipud(air1);
            air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air2=flipud(air2);
   
            i0=find(xdat(idat).x>6.5);            
%            i0=find(xdat(idat).x>8.5);            

%            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
            
             
%            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
            tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
            tot(idat).t=flipud(tot(idat).t);
            
             
		end
        
        
        
         xlims=1;
         xlimits=([-1 5]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
         

           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

    case 484
    idir=1;
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
 
    logflag=0;
    D=300e3;
    [iav iav2]=findheight( GridDan(2).Y1, GridDan(2).Y1(1)+D/2 , GridDan(2).Y1(end)-D/2  );	

      init = f*sum(icediagsALL(2).i(:,1,[37:42]),3)/npess2(2); %initial total water 

        idat=1;
	
        
      xdat(idat).x = f*mean( sum(TwoDDan(1).Q(:,[1:iav iav2:length(GridDan(1).Y1)] ,1),3) ,2) - init; 
      labs(idat).l=['Control vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(1).Q(:,[1:iav iav2:length(GridDan(1).Y1)] ,1:6),3) ,2) - init; 
      labs(idat).l=['Control tot']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1),3) ,2) - init; 
      labs(idat).l=['CCN=960 cm^{-3} vap']; 
      idat=idat+1;
      
      xdat(idat).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1:6),3) ,2) - init; 
      labs(idat).l=['CCN=960 cm^{-3} tot']; 
      idat=idat+1;
      
      
     
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
    %zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    
            
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
            
            rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air1=flipud(air1);
            air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air2=flipud(air2);
   
            i0=find(xdat(idat).x>6.5);            
%            i0=find(xdat(idat).x>8.5);            

%            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
            
             
%            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
            tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
            tot(idat).t=flipud(tot(idat).t);
            
             
		end
        
        
        
         xlims=1;
         xlimits=([-1 5]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
        % xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
         

           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case 483
    idir=1;
    
    
    %t1=23.5;
%    t1=GridDan(1).t(62)+3;
    
    t1=23.75;
    t1=24.6;
    
  %  t1=25.167; %01:10 UTC
    t1=22.667;
    t1=21.08;
    
    t1=23.46;
    
    
    
    
    
    
    it1=findheight(GridDan(idir).t+3,t1);
%    it1=size(icediagsALL(idir).i,2);
    
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    t2=25.1667;
%    t2=25;
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
    %it=62; % 00:50 UTC
    it=66; % 01:10 UTC

    

    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
  %  ylab='';
    
    
    
    titlenam='';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
        figname=['Total Water and Vapour'];

    savename=figname;
    
    meantot=mean(mean(tot_water44)); %mean at each height
	meantot2=repmat(meantot,[size(tot_water44,1) size(tot_water44,2) 1]);

   [maxtot]=max(mean( abs(tot_water44-meantot2) )); %max of the means along the slices (slices as in 2D orientation)
                                    %for each height
   [ratio]=squeeze(mean(mean( abs(tot_water44-meantot2) ))./maxtot); %ratio of horiz mean to slice with highest mean
                                  %to represent ratio of 2d with likely extrapolated mean over same area
  [maxtot imax]=max(max( abs(tot_water44-meantot2) ));
  imax=squeeze(imax);
  for km=1:length(imax)
      mean_at_maxs(km)=mean( abs(tot_water44(:,imax(km),km)-meantot2(:,imax(km),km)) );
  end
  [ratio2]=squeeze(mean(mean( abs(tot_water44-meantot2) )))./mean_at_maxs'; %ratio of horiz mean to slice with highest mean
                                  %to represent ratio of 2d with likely extrapolated mean over same area
   
    logflag=0;
    D=diff( GridDan(1).Y1([1 end]));
    [iav iav2]=findheight( GridDan(2).Y1, GridDan(2).Y1(1)+D/2 , GridDan(2).Y1(end)-D/2  );	

      xdat(1).x = f*sum(icediagsALL(2).i(:,1,[37:42]),3)/npess2(2); %initial total water 
      
%      xdat(2).x = f*mean( sum(TwoDDan(1).Q(:,:,1:6),3) ,2); %total water for 1st case
	

      xdat(2).x = f*mean( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1:6),3) ,2); 
      
      ixpos=152;
      xdat(3).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case
      labs(3).l=[runName(1).nam ' ' num2str(ixpos)]; 
      
      ixpos=132;
      xdat(4).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case
      labs(4).l=[runName(1).nam ' ' num2str(ixpos)]; 
      
      
	  ixpos=122;
      xdat(5).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case
      labs(5).l=[runName(1).nam ' ' num2str(ixpos)]; 
      
      
	  ixpos=11;
      xdat(6).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case
      labs(6).l=[runName(1).nam ' ' num2str(ixpos)]; 
      

% 	  ixpos=30;
       xdat(7).x = (  squeeze(f*mean(mean(tot_water44(2:end-1,:,:) )))   ); %total water for 1st case    
       labs(7).l=[runName(1).nam ' mean']; 
%       
% 		ixpos=40;
%       xdat(8).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case    
%       labs(8).l=[runName(1).nam ' ' num2str(ixpos)];       
%       
% 		ixpos=100;
%       xdat(9).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case    
%       labs(9).l=[runName(1).nam ' ' num2str(ixpos)];       
% 
 		ixpos=75;
       xdat(8).x = (f*mean( squeeze(tot_water44(2:end-1,ixpos,:)) ))'; %total water for 1st case    
       labs(8).l=[runName(1).nam ' ' num2str(ixpos)];    
       
       xdat(9).x = ( squeeze( f*max(max( tot_water44(2:end-1,:,:)) )) ); %total water for 1st case    
       labs(9).l=[runName(1).nam ' max']; 
       
       xdat(10).x = f*min( sum(TwoDDan(2).Q(:,[1:iav iav2:length(GridDan(2).Y1)] ,1:6),3) ,[],2);
       labs(10).l=[runName(2).nam ' max']; 
      
	
		labs(1).l='Initial';
       % labs(2).l='1 km tot';  %runName(1).nam;    %[t1str ' Control total water'];
	%	labs(3).l='2 km tot';  %runName(idir).nam;   %[t1str ' CCN 960 cm^{-3} tot'];

        labs(2).l=[runName(2).nam ' tot'];  %runName(1).nam;    %[t1str ' Control total water'];
		%labs(3).l=[runName(2).nam ' tot'];   %[t1str ' CCN 960 cm^{-3} tot'];

        
    T=tempLES(GridDan(idir)); %K
	ei=SatVapPress(T,'goff','ice'); %Pa
	P=GridDan(idir).PREFN; %Pa
	sat=f*0.622*ei./(P-ei);
    
    Lx=length(xdat);
    xdat(Lx+1).x = sat;
    labs(Lx+1).l=['Ice sat MR'];
    
     xdat(Lx+2).x = f*sum(icediagsALL(1).i(:,1,[37:42]),3); %initial total water 
     labs(Lx+2).l='Initial 3d';
        
    clear diff tot
    
    zmin=15;
%    zmin=14.5;
    zmax=16.7+0.62;
    zmax=19+0.62;
    
    ih2=findheight(GridDan(1).Z+620,28e3);
    ih2=250;
    
    ih2=findheight(GridDan(1).Z+620,17e3);
    ih1=findheight(GridDan(1).Z+620,15.335e3);
    ih1=findheight(GridDan(1).Z+620,14e3);
    
            
        for idat=1:length(xdat)
            ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
            
            rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=cumsum(flipud( diff( GridDan(1).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air1=flipud(air1);
            air2=cumsum(flipud( diff( GridDan(2).Y1([1 end]) ).*dz.*rho(ih1:ih2)  ) );
            air2=flipud(air2);
   
            i0=find(xdat(idat).x>6.5);            
%            i0=find(xdat(idat).x>8.5);            

%            xdat(idat).x(i0)=xdat(1).x(i0); %make them the same so diff=0 for tot > 6.5 ppmv
            
             
%            tot(idat)=sum( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  );
            tot(idat).t=cumsum(flipud( ( xdat(idat).x(ih1:ih2)-xdat(1).x(ih1:ih2) ).*dz.*rho(ih1:ih2)  ) );
            tot(idat).t=flipud(tot(idat).t);
            
             
		end
        
        
        
         xlims=1;
         xlimits=([4.8 5.4]);
       %  xlimits=([4.5 6.2]);
       %  xlimits=([4.5 7]);
         xlimits=([4 14]);

         %xlimits=[-0.1 0.4];
         

icum=0;  %flag to say to do normal means and not TTL cumulative means.
if icum==1
         
         
            clear xdat ydat labs
                  %  zmin=15;
                  %5  zmax=18.7+0.62;
                  
                  Y1=diff(GridDan(1).Y1([1 end]));                  
                  Y2=diff(GridDan(2).Y1([1 end]));
                  
                  %only for 3d
                  fact=Y2*1000/Y1^2; %conversion factor to multiply 2-d result by for fair comparison (assumes a 1 km length in 3-D)
                  
%                  fact=Y2*1000/Y1; %factor assuming that 2d covers full 300 km of 3d domain

                  
                xlims=0;  
                xdat(1).x=tot(2).t * diff( GridDan(1).Y1([1 end]) ) ./air1;  %tot
                xdat(1).x=tot(4).t * diff( GridDan(1).Y1([1 end]) ) ./air1;  %vapour

                
                %%%  2d
                xdat(2).x=tot(3).t * diff( GridDan(2).Y1([1 end]) ) ./air2;   %tot
                xdat(2).x=tot(5).t * diff( GridDan(2).Y1([1 end]) ) ./air2;   %vapour
                
                %%%% 3d
%                xdat(2).x=tot(3).t * diff( GridDan(2).Y1([1 end]) ) ./air2 * fact; %%%%  only use fact for 3d  %%%%%%%%

                xdat(1).x(end+1)=0;
                xdat(2).x(end+1)=0;
                
                ydat(1).y=GridDan(1).Z(ih1:ih2+1)/1000 + add_ground_height; 
                ydat(2).y=GridDan(2).Z(ih1:ih2+1)/1000 + add_ground_height; 
                
                labs(1).l = '1 km';
                labs(2).l = '2 km';                                
        
                labs(1).l = '3d';
                labs(2).l = '2d (scaled assuming 1 km in 3d)';                                

                labs(1).l = '1 km';
                labs(2).l = runName(2).nam;      
                
                
                xlab='Cumulative total water reduction (kg m^{-1})';
                xlab='Mean total water MR change from 17 km downwards (ppmv)';
                
                        figname=['Cumulative total water'];
                        savename=figname;

end              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         
    case 482    

    
    xlab=['Ice saturation mixing ratio (ppmv)'];
	ylab='Mean TTL total water mixing ratio (ppmv)';
            
    titlenam='';
    figname=['Mean TTL function of ice sat'];
    savename=figname;
    
    logflag=0;
    
    clear diff tot
    
    TTL_avTOT_funcof_IcesatMR; %works out mean TTL from 15 km to 17 km for various ice sat MRs
    
    ydat(1).y=total(1).t;
    ydat(2).y=total(2).t;
    
	
	xdat(1).x=ppmvs; 
	xdat(2).x=ppmvs; 
	
	labs(1).l = '1 km resolution';
	labs(2).l = '2 km resolution';                                
	
	
	nmark=0;
	
	lor=2;
    
    izlim=0;

         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 488
    
    zmin=14;
    zmax=20.5;
    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
    
    
    
    titlenam='Vapour profiles';
%    figname=['Total Water and Vapour at Time = ' num2str(mod(t1,24),'%2.2f') ' UTC'];
    figname=['Water Vapour'];

    savename=figname;
    
    logflag=0;
    
	idirs=[1:4];
    
    for idat=1:length(idirs)
        idir=idirs(idat);
        
      [i1,i2]=findheight(GridDan(idir).Y1,GridDan(idir).Y1(end)-500e3,GridDan(idir).Y1(1)+500e3);

      xdat(idat).x = f*mean(TwoDDan(idir).Q(:,[i1:length(GridDan(idir).Y1) 1:i2],1),2);
	  labs(idat).l=[runName(idir).nam];
      ydat(idat).y=GridDan(idir).Z/1000+0.62;
        
  end
                            	
         xlims=1;
         xlimits=([4.8 5.8]);
         xlimits=([4.5 6.2]);
         xlimits=([4.5 7]);
         xlimits=([4.8 5.1]);

         %xlimits=[-0.1 0.4];
         
         nmark=0;
         
         lor=2; %1=right, 2=left
         
    case 47
        
    zmin=12;
    zmax=21.3;
    
    xlab=['Mixing Ratio (ppmv)'];
	ylab='Height (km)';
    
    figname=['IceSat_initVap'];
    savename=figname;
    
    logflag=0;
    
    T=TempLes(GridDan(idir));
    P=GridDan(idir).PREFN;
    xdat(1).x = SatVapPress(T,'lem','ice',P,1);    
    xdat(2).x = f*icediagsALL(idir).i(:,1,[37])/npess2(idir);
	for idat=1:2
        ydat(idat).y = GridDan(idir).Z/1000 + add_ground_height; 
	end
	
        labs(1).l='Initial Ice Sat Mixing Ratio';
		labs(2).l='Initial Vap Mixing Ratio';
	
         xlims=1;
         xlimits=([3 20]);
         
         nmark=0;
         
         
	case 46
            timeseries_dan
    
    case 45
    
    time1=19.75;
  %  time1=20;
    
    time2=23.46;
    
    time2=20.5;
    
    t1=findheight(GridDan(1).t+3,time1);
    t2=findheight(GridDan(1).t+3,time2);
    
    H=18e3;
    H=20e3;    
    ih=findheight(GridDan(1).Z+620,H);
    
%    t1=1;
%    t2=21;
    
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
 

    
    xlab=['Total mass flow rate (kg s^{-1} m^{-1})'];
    
    
    
  %  xlab=['Total mass flow rate in updraughts larger than 1 m s^{-1} (kg/s)'];
    
  %  xlab=['Updraught diameter (for updraughts larger than 1 m s^{-1} )(km)'];

	ylab='Height (km)';
    
   
    
    
    titlenam='';
    
    logflag=0;
    
    aind=280;%ALu_A
    aind=285; %ACu_A
%    aind=284; %ACC_A
    aind=283; %W>1_A

aind=[];   %if want to ignore area and just get total
    
    datind=137; %ALu_W
 %   datind=47;  %ALu_WQ02
 %   datind=198; %ACC_Q04
 %   datind=193; %ACC_WQ01
 %   datind=223; %ACC_Q01
%    datind=194; %ACC_WQ

plotcase='ACu_W';
plotcase='W>1_W';
%plotcase='ACu_ice_flux';
plotcase='Ice_mass';

switch plotcase
case 'W>1_W'    
    datind=303; %W>1_W
    xlab=['Average vertical wind speed in updraughts larger than 1 m s^{-1}'];
    aind=[];
    aind=283; %ACu_A
    
case 'ACu_W'    
    datind=302; %ACu_W
    xlab=['Average cloudy updraught (m s^{-1})']; 
    aind=[];
    aind=285; %ACu_A
    
case 'ACu_ice_flux'    
    datind=[49:51]; %ACu_W
    xlab=['Ice mass flux (kg m^{-2} s^{-1})'];     
    aind=[];
    aind=285; %ACu_A
    
case 'Ice_mass'    
    datind=[226:228]; %ACC_Q
    xlab=['Ice mixing ratio (kg kg^{-1})'];     
    aind=[];
    aind=284; %ACC_A    
end


savename=xlab;
    
%     figname=['Average vertical wind speed in updraughts larger than 1 m s^{-1}'];    
%     figname=['Total mass flow rate in updraughts larger than 1 m s^{-1} (kg s^{-1})'];
    
    figname=xlab;
    
    idats=[1 2 3];
    
    clear diff
    
	for iidat=1:length(idats)
       idat=idats(iidat); 

         if length(aind)==1            
            area=icediagsALL(idat).i(1:ih,t1:t2,aind)/npess2(idat);
            area(area==0)=1;
        else
            area=1;
        end
        
        
%        xdat(idat).x = mean(icediagsALL(idat).i(:,t1:t2,[137]),2)/npess2(idat); %dividing by no. processors
        xdat(iidat).x = mean(sum(icediagsALL(idat).i(1:ih,t1:t2,[datind]),3)./area,2)/npess2(idat); %ALu_W. Dividing by no. processors
        
%        xdat(idat).x = xdat(idat).x .* length(GridDan(idat).Y1) .* GridDan(idat).RHON(1:ih) * diff( GridDan(idat).Y1(1:2) ) ; %to convert into mass flow per unit in 3rd dim
           %multiply to length(Y) since is averaged over whole domain. dy*rho gives the mass of each element (already have dz from W)
           
     %   xdat(idat).x = mean(icediagsALL(idat).i(:,t1:t2,[1]),2);
        ydat(iidat).y = (GridDan(idat).Z(1:ih)+620)/1000; 
        
        i3d=0;
        if i3d==1
            if idat==1
                    xdat(idat).x = 1000 * xdat(idat).x .* length(GridDan(idat).Y1) .* GridDan(idat).RHON(1:ih) * diff( GridDan(idat).Y1(1:2) ) ; %to convert into mass flow per unit in 3rd dim
            else                    %this factor of 1000 is to simulate the 2d case occupying 1km (1000m) of 3rd dim 
                    xdat(idat).x = xdat(idat).x .* length(GridDan(idat).Y1) * length(GridDan(idat).X1) .* GridDan(idat).RHON(1:ih) * diff( GridDan(idat).Y1(1:2) )^2 ; %to convert into mass flow per unit in 3rd dim
            end
        end
        
        i3dA=0;
        if i3dA==1
            if idat==1
                    xdat(idat).x = length(GridDan(idat).Y1) .* diff( GridDan(idat).Y1(1:2) ) * mean( icediagsALL(idat).i(1:ih,t1:t2,aind),2 )/npess2(idat) /1000;
            else                    %this factor of 1000 is to simulate the 2d case occupying 1km (1000m) of 3rd dim 
                    xdat(idat).x = (length(GridDan(idat).Y1) .* diff( GridDan(idat).Y1(1:2) ) )^2 * mean( icediagsALL(idat).i(1:ih,t1:t2,aind),2 )/npess2(idat) /1000;
                    xdat(idat).x = sqrt(xdat(idat).x / pi )/1000;
            end
        end
        
        labs(iidat).l=runName(idat).nam;
        	
	end
	
%         labs(1).l='250m res';
% 		labs(2).l='500m res';
% 		labs(3).l='1km res';
		%labs(4).l='250m res 13.4';

	
         xlims=0;
         xlimits=[-30 50];
	%     
         
         nmark=0;
    
    
    case 44
    
    time1=19.75;
    time2=23.5;
    time2=21.25;
    time2=20.4;
    
    t1=findheight(GridDan(1).t+3,time1);
    t2=findheight(GridDan(1).t+3,time2);
    
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
 
    xlab=['Max Vertical Wind Speed (m/s)'];
	ylab='Height (km)';
    
    figname=['MaxW wind speed' direc(idir).dir];
    
	for idat=1:4
        xdat(idat).x = max(MaxW(idat).w(:,t1:t2),[],2);
        xdat(idat).x = mean(MaxW(idat).w(:,t1:t2),2);
        
        ydat(idat).y = (GridDan(idat).Z+620)/1000; 
        labs(idat).l = runName(idat).nam;	
	end
	
%         labs(1).l='250m res';
% 		labs(2).l='500m res';
% 		labs(3).l='1km res';
% 		labs(4).l='250m res 13.4';

	
         xlims=0;
         xlimits=[-30 50];
	%     
         
         nmark=0;
         
         savename=figname;
         
    
    case 43
    logflag=0;
    
    izlim=0;
    %z=GridDan(idir).Z;
 
    xlab=['Horizontal Wind Speed (m/s)'];
	ylab='Height (km)';
    
    figname=['V wind speed' direc(idir).dir];
    
	for idat=1:2
        xdat(idat).x = GridDan(idat).VBAR;
        ydat(idat).y = (GridDan(idat).Z+620)/1000; 
        	
	end
	
        labs(1).l='24/02';
        labs(2).l='13/02';

	
         xlims=0;
         xlimits=[-30 50];
	%     
         
         nmark=0;
         
         
         
    case 42
    logflag=2; %for y axis
    
    izlim=0;
    %z=GridDan(idir).Z;
 
    xlab=['Vertical Velocity (m/s)'];
%	ylab='log10 of normalised dn/dw (where n is no. grid points)';

	ylab='Relative Frequency';

    figname=['w distribution'];
    
    %%%%%% note make sure dw (size of vertical velocity bin width) is correct ! %%%%
    dw=0.5;
    
	for idat=1:4
        xdat(idat).x = wb(idat).w-dw/2;
        ydat(idat).y=(wd(idat).w*dw); %data in wd(idat).w is df/dw (f=relative frequency) so divide by dw to get just the rel. freq.
	end
	
	labs(1).l='250m res';
	labs(2).l='500m res';
	labs(3).l='1km res';
	labs(4).l='250m res 13.4';
	
     xlims=1;
     xlimits=[-1 47];
	%     
     nmark=-1;
     
     savename=[figname '_' ylab];
     
    case 41
        ice_budget_test
	case 40
        ice_budget
    case 39
        budget_tim
	case 38 
        ice_sat_mr;
    
    case 37
    logflag=0;
    
    %z=GridDan(idir).Z;
 
    xlab=['Fraction of mass received'];
	ylab='Height (m)';
    
    figname=['dm fraction' direc(idir).dir];
    

    xdat(1).x = dqz2/sumdqz2;
    ydat(1).y=z(z0:zend); %620m already added on in fourevents2.m
    labs(1).l='Fraction';
    
% %     ydat(2).y=z(len:length(z)-1)/1000;
% %     labs(2).l='Vapour';
%     
%     secyA=z/1000 + 0.62;
%     secyB=GridDan(idir).PREFN/100;
%     lab2='Pressure (hPa)';  
%     dual=1;
%     
     xlims=1;
     xlimits=[0 1];
%     
     zmin=15.6e3;  %(km)
     zmax=17.2e3;
     
     nmark=-1;
     
    case 36
    logflag=0;
    
    %z=GridDan(idir).Z;
 
    xlab=['Time UTC'];
	ylab='Mass Decrease (kg/m^2)';
    
    figname=['tot dm' direc(idir).dir];
    

    xdat(1).x = GridDan(idir).t(1:size(mtim,2)) + 3;
    ydat(1).y=mtim;
    labs(1).l='Tot mass decrease due to dry points';
    
% %     ydat(2).y=z(len:length(z)-1)/1000;
% %     labs(2).l='Vapour';
%     
%     secyA=z/1000 + 0.62;
%     secyB=GridDan(idir).PREFN/100;
%     lab2='Pressure (hPa)';  
%     dual=1;
%     
     xlims=1;
     xlimits=[19.5 23];
%     
     zmin=0;  %(km)
     zmax=max(ydat(1).y)*1.1;
     
     nmark=-1;
    

    
    lor=1;
    
    case 35
    logflag=0;
    
    %z=GridDan(idir).Z;
 
    xlab=['dN/dt - Number of overshoot events per month'];
	ylab='Min vapour in 15.8-17km region (ppmv)';
    
    figname=['min q vs dn/dt ' direc(idir).dir];
    

    xdat(1).x = Nevs;
    ydat(1).y=qmin_spread_mass;
    labs(1).l='Min ppmv';
    
% %     ydat(2).y=z(len:length(z)-1)/1000;
% %     labs(2).l='Vapour';
%     
%     secyA=z/1000 + 0.62;
%     secyB=GridDan(idir).PREFN/100;
%     lab2='Pressure (hPa)';  
%     dual=1;
%     
%     xlims=1;
%     %xlimits=[-1.4 0.4];
%     xlimits=[-80 -65];
%     
     zmin=3;  %(km)
     zmax=5;
    

    
    lor=1;
    
	case 34
    Ad_source
    
    case 33
    logflag=0;
    
    xlab='Mixing Ratio Gain After Overshoot (ppmv)';
	ylab='Height (km)';
    
    figname='Net Mass Gain';
    
    xlims=1;
    xlimits=[-1 1];
    
    zmin=14;  %(km)
    zmax=21;   %z(end)/1000;
    

    z=GridDan(idir).Z;
    
    
    dz=GridDan(idir).Z(2:end)-GridDan(idir).Z(1:end-1); 
    rho=GridDan(idir).RHO;
    
    t1=3;
    t2=17; %index for start and end times of difference plot
    
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
       xdat(1).x = f*( sum(icediag4(idir).i(:,t2,[35:42]),3) - sum(icediag4(idir).i(:,t1,[35:42]),3) ); %total water change
       xdat(2).x = f*( sum(icediag4(idir).i(:,t2,[35:36]),3) - sum(icediag4(idir).i(:,t1,[35:36]),3) ); %difference in end and start vapour diags   
    
    ydat(1).y=z/1000;
    labs(1).l='Total Water';
    
    ydat(2).y=z/1000;
    labs(2).l='Vapour';
    
%     ydat(3).y=Grid.Z(2:len+1)/1000;
%     labs(3).l='Microphysics';
%     
%     ydat(4).y=Grid.Z(2:len+1)/1000;
%     labs(4).l='Fall Speed Flux';
%     
%     ydat(5).y=Grid.Z(2:len+1)/1000;
%     labs(5).l='Ice Flux';
%     
%     ydat(6).y=Grid.Z(2:len+1)/1000;
%     labs(6).l='Upwards Ice Flux';
%     
%     ydat(7).y=Grid.Z(2:len+1)/1000;
%     labs(7).l='Downwards Ice Flux';
%     
%     ydat(8).y=Grid.Z(2:len+1)/1000;
%     labs(8).l='Net Vapour Flux';
%     
%     ydat(9).y=Grid.Z(2:len+1)/1000;
%     labs(9).l='Ice Fall + Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=GridDan(idir).PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
       
    case 32
    logflag=0;
    
    z=GridDan(idir).Z;
 
    xlab=['Temp (^oC)'];
	ylab='Height (km)';
    
    figname=['Temp ' direc(idir).dir];
    

    xdat(1).x = tempLES(GridDan(idir)) - 273.15;
    ydat(1).y=z/1000 + 0.62;
    labs(1).l='Temperature';
    
%     ydat(2).y=z(len:length(z)-1)/1000;
%     labs(2).l='Vapour';
    
    secyA=z/1000 + 0.62;
    secyB=GridDan(idir).PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    xlims=1;
    %xlimits=[-1.4 0.4];
    xlimits=[-80 -65];
    
    zmin=14.5;  %(km)
    zmax=20;
    

    
    lor=1;
    
    case 31
    logflag=0;
    
       %idir=4;
    
%     len=length(icediag4(idir).i(:,1,35))-1;
    z=GridDan(idir).Z;
    len=findheight(z,15.05e3)-1;
    %len=length(z)-1;
    
    hstr=num2str( round2(z(len)/1000,1) );
    xlab=['Cumulative Net Gain up from ' hstr 'km (kg/m^2)'];
	ylab='Height (km)';
    
    figname=['Cumulative net gain up from ' hstr 'km ' direc(idir).dir];
    
 
    dz=z(2:end)-z(1:end-1); %finds d/dz of flux
    rho=GridDan(idir).RHO;
    
    for i=len:length(z)-1
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                 
        %airmass=sum( rho(len:i) .*dz(len:i) );
        airmass=1;
        
        %converted into mixing ratio
        xdat(1).x(i-len+1) = sum( rho(len:i) .*dz(len:i) .* ( sum(icediag4(idir).i(len:i,end,[35:42]),3) - sum(icediag4(idir).i(len:i,3,[35:42]),3) )  ) ./ airmass;
        xdat(2).x(i-len+1) = sum( rho(len:i) .*dz(len:i) .* ( sum(icediag4(idir).i(len:i,end,[35:36]),3) - sum(icediag4(idir).i(len:i,3,[35:36]),3) )  ) ./ airmass;
        
        %mass gain
        %xdat(1).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediag3(1).i(len-i+1:len,end,[2 5 8 11]),3) - sum(icediag3(1).i(len-i+1:len,1,[2 5 8 11]),3) )  ) ;
        %xdat(2).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ) ; %difference in end and start vapour diags    
    end
    
    
    ydat(1).y=z(len:length(z)-1)/1000;
    labs(1).l='Total Water';
    
    ydat(2).y=z(len:length(z)-1)/1000;
    labs(2).l='Vapour';
    
    secyA=z/1000;
    secyB=GridDan(idir).PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    %xlims=1;
    %xlimits=[-1.4 0.4];
    xlimits=[-2e-4 2e-4];
    
    zmin=14.5;  %(km)
    zmax=20;
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    
    
    
    case 30
    logflag=0;
    
    xlab='Mean mixing ratio (ppmv)';
	ylab='Height (km)';
    
    figname='Mean mixing ratio';
    
    
    
    zmin=13;  %(km)
    zmax=19;

    idir=4;    
        xdat(1).x = f*squeeze(sum(icediag4(idir).i(:,2,[35 36]),3));
        xdat(2).x = f*squeeze(sum(icediag4(idir).i(:,86,[35 36]),3));
        xdat(3).x = f*squeeze(sum(icediag4(idir).i(:,2,[35:42]),3));
        xdat(4).x = f*squeeze(sum(icediag4(idir).i(:,86,[35:42]),3));
        
        
    
    
      
    for ji=1:4
        ydat(ji).y=GridDan(idir).Z(1:end)/1000;
        secyA=GridDan(idir).Z/1000;
        secyB=GridDan(idir).PREFN/100;
    end
    
    
    
    lab2='Pressure (hPa)';

    labs(1).l='Start Vapour';
    labs(2).l='End Vapour';
    labs(3).l='Start Total Water';
    labs(4).l='End Total Water';
    
    
    %ydat(4).y=Grid.PREFN(1:len)/100;
    
    
    dual=1;
    
    xlims=1;
    xlimits=[3.5 8];
    
    case 29
    logflag=0;
    
    ylab=['Max Ice Concentration (kg^{-1})'];
	xlab='Time (UTC)';
    
    figname=['Max Ice NC'];
    
   
    
    
    
    j=1;
    ydat(j).y=SerDan(j).SER(:,44);
    xdat(j).x=19.75+SerDan(j).SER(:,1)./3600;
    labs(j).l='25km Domain';
    
    j=2;
    ydat(j).y=SerDan(j).SER(:,44);
    xdat(j).x=19.75+SerDan(j).SER(:,1)./3600;
    labs(j).l='Large Damping Layer'; 
    xlims=1;
    %xlimits=[-1.4 0.4];
    xlimits=[time(1)+3 time(end)+3];
    
    j=3;
    ydat(j).y=SerDan(j).SER(:,44);
    xdat(j).x=19.75+SerDan(j).SER(:,1)./3600;
    labs(j).l='CCN = 960 cm^{-1}';
    
    j=4;
    ydat(j).y=SerDan(j).SER(:,44);
    xdat(j).x=19.75+SerDan(j).SER(:,1)./3600;
    labs(j).l='Small Damping Layer';
    
    
    
    zmin=0;  
    zmax=5;
    
    izlim=0; %flag to stop manual scaling of the y axis
    ixtime=1;
    
    case 28
    logflag=0;
    
    ylab=['Min Water Vapour (ppmv)'];
	xlab='Time (UTC)';
    
    figname=['Min Vap'];
    
   
    
    
    
    j=1;
    ydat(j).y=f*min(vapprc(j).v(:,:,1),[],1);
    xdat(j).x=3+time(1:length(ydat(j).y));
    labs(j).l='25km Domain';
    
    j=2;
    ydat(j).y=f*min(vapprc(j).v(:,:,1),[],1);
    xdat(j).x=3+time(1:length(ydat(j).y));
    labs(j).l='Large Damping Layer'; 
    xlims=1;
    %xlimits=[-1.4 0.4];
    xlimits=[time(1)+3 time(end)+3];
    
    j=3;
    ydat(j).y=f*min(vapprc(j).v(:,:,1),[],1);
    xdat(j).x=3+time(1:length(ydat(j).y));
    labs(j).l='CCN = 960 cm^{-1}';
    
    j=4;
    ydat(j).y=f*min(vapprc(4).v(:,:,1),[],1);
    xdat(j).x=3+time(1:length(ydat(j).y));
    labs(j).l='Small Damping Layer';
    
    
    
    zmin=0;  
    zmax=5;
    
    %izlim=0; %flag to stop manual scaling of the y axis
    ixtime=1;
    
    case 27
    logflag=0;
    
    
    
    
    len=length(icediag(1).i(:,1,11))-1;
    len=findheight(z,15.05e3)-1;
    %len=length(z)-1;
    
    hstr=num2str( round2(z(len)/1000,1) );
    xlab=['Cumulative Net Gain up from ' hstr 'km (ppmv)'];
	ylab='Height (km)';
    
    figname=['Cumulative net gain up from ' hstr 'km'];
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); %finds d/dz of flux
    rho=Grid.RHO;
    
    for i=len:length(z)-1
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                 
        airmass=sum( rho(len:i) .*dz(len:i) );
        %airmass=1;
        
        %converted into mixing ratio
        xdat(1).x(i-len+1) = f* sum( rho(len:i) .*dz(len:i) .* ( sum(icediag3(1).i(len:i,end,[2 5 8 11]),3) - sum(icediag3(1).i(len:i,1,[2 5 8 11]),3) )  ) ./ airmass;
        xdat(2).x(i-len+1) = f* sum( rho(len:i) .*dz(len:i) .*( icediag3(1).i(len:i,end,11) - icediag3(1).i(len:i,1,11) ) ) ./ airmass; %difference in end and start vapour diags    
        
        %mass gain
        %xdat(1).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediag3(1).i(len-i+1:len,end,[2 5 8 11]),3) - sum(icediag3(1).i(len-i+1:len,1,[2 5 8 11]),3) )  ) ;
        %xdat(2).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ) ; %difference in end and start vapour diags    
    end
    
    
    ydat(1).y=Grid.Z(len:length(z)-1)/1000;
    labs(1).l='Total Water';
    
    ydat(2).y=Grid.Z(len:length(z)-1)/1000;
    labs(2).l='Vapour';
    
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    %xlims=1;
    %xlimits=[-1.4 0.4];
    xlimits=[-2e-4 2e-4];
    
    zmin=14.5;  %(km)
    zmax=20;
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=2;
    
    
    case 26
    logflag=0;
    
    xlab='Net Mass Gain (kg/m^2)';
	ylab='Height (km)';
    
    figname='Net Mass Gain';
    
    %xlims=1;
    xlimits=[-3e-5 5e-5];
    
    zmin=14.5;  %(km)
    zmax=21;   %z(end)/1000;
    
    len=length(m(1).m);
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); 
    rho=Grid.RHO;
    
    
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
       xdat(1).x = f*( sum(icediag3(1).i(2:len+1,end,[2 5 8 11]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8 11]),3) ); %total water change
       xdat(2).x = f*( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        %xdat(1).x =  ( sum(icediag3(1).i(2:len+1,end,[2 5 8]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8]),3) ); %total water change
        %xdat(2).x = ( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        xdat(3).x = f*( ( m(11).m(1:len)  ) )./ (rho(1:len) .*dz(1:len)); %microphysics from ALL_DQ01                             
        
        %xdat(4).x = f*icediag3(1).i(1:len,1,11),1;
        %xdat(5).x = f*icediag3(1).i(1:len,88,11),1;
        
        xdat(4).x =  f* ( m(9).m(1:len)  ) ./ (rho(1:len) .*dz(1:len)); %fall speed flux                               
        xdat(5).x = - f* ( m(38).m(1:len) ) ./ (rho(1:len) .*dz(1:len)); %flux of ice+snow+graupel
        xdat(6).x = - f* ( m(18+5).m(1:len) + m(18+7).m(1:len) + m(18+9).m(1:len) + m(18+11).m(1:len) + m(18+13).m(1:len) + m(18+15).m(1:len) ) ./ (rho(1:len) .*dz(1:len)); %up 
        xdat(7).x = - f* ( m(19+5).m(1:len) + m(19+7).m(1:len) + m(19+9).m(1:len) + m(19+11).m(1:len) + m(19+13).m(1:len) + m(19+15).m(1:len) ) ./ (rho(1:len) .*dz(1:len)); %down   
        xdat(8).x = - f* ( m(1).m(1:len) ) ./ (rho(1:len) .*dz(1:len));
        xdat(9).x = xdat(4).x + xdat(5).x;

    
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        %xdat(1).x(len-i+1) = xdat(1).x(len-i+1) + sumice ;    
       
    
    
    
    ydat(1).y=Grid.Z(2:len+1)/1000;
    labs(1).l='Total Water';
    
    ydat(2).y=Grid.Z(2:len+1)/1000;
    labs(2).l='Vapour';
    
    ydat(3).y=Grid.Z(2:len+1)/1000;
    labs(3).l='Microphysics';
    
    ydat(4).y=Grid.Z(2:len+1)/1000;
    labs(4).l='Fall Speed Flux';
    
    ydat(5).y=Grid.Z(2:len+1)/1000;
    labs(5).l='Ice Flux';
    
    ydat(6).y=Grid.Z(2:len+1)/1000;
    labs(6).l='Upwards Ice Flux';
    
    ydat(7).y=Grid.Z(2:len+1)/1000;
    labs(7).l='Downwards Ice Flux';
    
    ydat(8).y=Grid.Z(2:len+1)/1000;
    labs(8).l='Net Vapour Flux';
    
    ydat(9).y=Grid.Z(2:len+1)/1000;
    labs(9).l='Ice Fall + Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    case 25
    logflag=0;
    
    xlab='Net Mass Gain (kg/m^2)';
	ylab='Height (km)';
    
    figname='Net Mass Gain';
    
    xlims=1;
    xlimits=[-3e-5 5e-5];
    
    zmin=14.5;  %(km)
    zmax=21;   %z(end)/1000;
    
    len=length(m(1).m);
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); 
    rho=Grid.RHO;
    
    
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
       xdat(1).x = rho(1:len) .*dz(1:len) .* ( sum(icediag3(1).i(2:len+1,end,[2 5 8 11]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8 11]),3) ); %total water change
       xdat(2).x = rho(1:len).*dz(1:len).*( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        %xdat(1).x =  ( sum(icediag3(1).i(2:len+1,end,[2 5 8]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8]),3) ); %total water change
        %xdat(2).x = ( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        xdat(3).x = ( ( m(11).m(1:len)  ) ); %microphysics from ALL_DQ01                             
        xdat(4).x =  ( m(9).m(1:len)  ) ; %fall speed flux                               
        xdat(5).x = - ( m(38).m(1:len) ); %flux of ice+snow+graupel
        xdat(6).x = - ( m(18+5).m(1:len) + m(18+7).m(1:len) + m(18+9).m(1:len) + m(18+11).m(1:len) + m(18+13).m(1:len) + m(18+15).m(1:len) ); %up 
        xdat(7).x = - ( m(19+5).m(1:len) + m(19+7).m(1:len) + m(19+9).m(1:len) + m(19+11).m(1:len) + m(19+13).m(1:len) + m(19+15).m(1:len) ); %down   
        xdat(8).x = - ( m(1).m(1:len) );
        %xdat(9).x = - ( m(12).m + m(13).m + m(14).m );

    
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        %xdat(1).x(len-i+1) = xdat(1).x(len-i+1) + sumice ;    
       
    
    
    
    ydat(1).y=Grid.Z(2:len+1)/1000;
    labs(1).l='Total Water';
    
    ydat(2).y=Grid.Z(2:len+1)/1000;
    labs(2).l='Vapour';
    
    ydat(3).y=Grid.Z(2:len+1)/1000;
    labs(3).l='Microphysics';
    
    ydat(4).y=Grid.Z(2:len+1)/1000;
    labs(4).l='Fall Speed Flux';
    
    ydat(5).y=Grid.Z(2:len+1)/1000;
    labs(5).l='Ice Flux';
    
    ydat(6).y=Grid.Z(2:len+1)/1000;
    labs(6).l='Upwards Ice Flux';
    
    ydat(7).y=Grid.Z(2:len+1)/1000;
    labs(7).l='Downwards Ice Flux';
    
    ydat(8).y=Grid.Z(2:len+1)/1000;
    labs(8).l='Net Vapour Flux';
    
    %ydat(9).y=Grid.Z(2:len+1)/1000;
    %labs(9).l='Net Water Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
case 24
    logflag=0;
    
    xlab='Net Mass Gain (kg/m^2)';
	ylab='Height (km)';
    
    figname='Net Mass Gain';
    
    xlims=1;
    xlimits=[-1e-5 1.5e-5];
    
    zmin=15.8;  %(km)
    zmax=20;   %z(end)/1000;
    
    len=length(m(1).m);
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); 
    rho=Grid.RHO;
    
    
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
       xdat(1).x = rho(1:len) .*dz(1:len) .* ( sum(icediag3(1).i(2:len+1,end,[2 5 8]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8]),3) ); %total water change
       xdat(2).x = rho(1:len).*dz(1:len).*( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        %xdat(1).x =  ( sum(icediag3(1).i(2:len+1,end,[2 5 8]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8]),3) ); %total water change
        %xdat(2).x = ( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        xdat(3).x = ( ( m(11).m(1:len)  ) ); %microphysics from ALL_DQ01                             
        %xdat(4).x =  ( xdat(1).x - xdat(2).x + xdat(3).x  ) ; %fall speed flux                               
        xdat(4).x =  -m(1).m(1:len) + m(11).m(1:len) ; %vapour budget from flux & microphysics                               

        %xdat(5).x =  xdat(4).x .* (  m(9).m(1:len) ./ ( m(38).m(1:len) + m(9).m(1:len) )  )   ; %flux of ice+snow+graupel
        %xdat(5).x = - (m(3).m(1:len) + m(4).m(1:len) + m(5).m(1:len) );
        xdat(5).x =  -m(38).m - m(11).m + m(9).m(1:len); %ice budget from flux
        xdat(6).x =  ( m(9).m(1:len)  ) ; %fall speed flux 
        xdat(7).x =  m(9).m(1:len) - m(38).m(1:len) ; %down
        xdat(8).x =  - m(1).m(1:len) ; %vapour flux gain
        
        %xdat(7).x = - ( m(18+5).m(1:len) + m(18+7).m(1:len) + m(18+9).m(1:len) + m(18+11).m(1:len) + m(18+13).m(1:len) + m(18+15).m(1:len) );
        %xdat(8).x = - ( m(19+5).m(1:len) + m(19+7).m(1:len) + m(19+9).m(1:len) + m(19+11).m(1:len) + m(19+13).m(1:len) + m(19+15).m(1:len) );
        %xdat(9).x = xdat(7).x + xdat(8).x;
    
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        %xdat(1).x(len-i+1) = xdat(1).x(len-i+1) + sumice ;    
       
    
    
    
    ydat(1).y=Grid.Z(2:len+1)/1000;
    labs(1).l='Ice';
    
    ydat(2).y=Grid.Z(2:len+1)/1000;
    labs(2).l='Vapour';
    
    ydat(3).y=Grid.Z(2:len+1)/1000;
    labs(3).l='Microphysics';
    
    ydat(4).y=Grid.Z(2:len+1)/1000;
    labs(4).l='Vapour - flux and microphysics';
    
    ydat(5).y=Grid.Z(2:len+1)/1000;
    labs(5).l='Ice - flux, fall speed & microphysics';
    
    ydat(6).y=Grid.Z(2:len+1)/1000;
    labs(6).l='Fall speed Flux';
    
     ydat(7).y=Grid.Z(2:len+1)/1000;
     labs(7).l='Ice - flux & fall speed flux';
%     
     ydat(8).y=Grid.Z(2:len+1)/1000;
     labs(8).l='Vapour Flux';
%     
%     ydat(9).y=Grid.Z(2:len+1)/1000;
%     labs(9).l='Net Water Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    
case 23
    logflag=0;
    
    xlab='Net Mass Gain (kg/m^2)';
	ylab='Height (km)';
    
    figname='Net Mass Gain';
    
    xlims=1;
    xlimits=[-3e-5 5e-5];
    
    zmin=14.5;  %(km)
    zmax=21;   %z(end)/1000;
    
    len=length(m(1).m);
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); 
    rho=Grid.RHO;
    
    
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
       xdat(1).x = rho(1:len) .*dz(1:len) .* ( sum(icediag3(1).i(2:len+1,end,[2 5 8]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8]),3) ); %total water change
       xdat(2).x = rho(1:len).*dz(1:len).*( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        %xdat(1).x =  ( sum(icediag3(1).i(2:len+1,end,[2 5 8]),3) - sum(icediag3(1).i(2:len+1,1,[2 5 8]),3) ); %total water change
        %xdat(2).x = ( icediag3(1).i(2:len+1,end,11) - icediag3(1).i(2:len+1,1,11) ) ; %difference in end and start vapour diags                                
        xdat(3).x = ( ( m(11).m(1:len)  ) ); %microphysics from ALL_DQ01                             
        xdat(4).x =  ( m(9).m(1:len)  ) ; %fall speed flux                               
        xdat(5).x = - ( m(38).m(1:len) ); %flux of ice+snow+graupel
        xdat(6).x = - ( m(19).m(1:len)  + m(21).m(1:len) ); %up   
        xdat(7).x = - (  m(20).m(1:len) + m(22).m(1:len) ); %down
        xdat(8).x = - ( m(1).m(1:len) );
        %xdat(9).x = - ( m(12).m + m(13).m + m(14).m );

    
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        %xdat(1).x(len-i+1) = xdat(1).x(len-i+1) + sumice ;    
       
    
    
    
    ydat(1).y=Grid.Z(2:len+1)/1000;
    labs(1).l='Total Water';
    
    ydat(2).y=Grid.Z(2:len+1)/1000;
    labs(2).l='Vapour';
    
    ydat(3).y=Grid.Z(2:len+1)/1000;
    labs(3).l='Microphysics';
    
    ydat(4).y=Grid.Z(2:len+1)/1000;
    labs(4).l='Fall Speed Flux';
    
    ydat(5).y=Grid.Z(2:len+1)/1000;
    labs(5).l='Ice Flux';
    
    ydat(6).y=Grid.Z(2:len+1)/1000;
    labs(6).l='Upwards Vapour Flux';
    
    ydat(7).y=Grid.Z(2:len+1)/1000;
    labs(7).l='Downwards Vapour Flux';
    
    ydat(8).y=Grid.Z(2:len+1)/1000;
    labs(8).l='Net Vapour Flux';
    
    %ydat(9).y=Grid.Z(2:len+1)/1000;
    %labs(9).l='Net Water Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    
    
    
    case 22
    logflag=0;
    
    xlab='Net Gain in Vapour (kg/m^2)';
	ylab='Height (km)';
    
    figname='Vapour Gain End-Start';
    
    
    
    zmin=0;  %(km)
    zmax=z(end)/1000;
    
    len=length(icediag3(1).i(:,1,11))-1;
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); 
    rho=Grid.RHO;
    
    
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        xdat(1).x = rho(1:len) .*dz(1:len) .* ( sum(icediag3(1).i(1:len,end,[2 5 8 11]),3) - sum(icediag3(1).i(1:len,1,[2 5 8 11]),3) );
        xdat(2).x = (rho(1:len).*dz(1:len).*( icediag3(1).i(1:len,end,11) - icediag3(1).i(1:len,1,11) ) ); %difference in end and start vapour diags                                
        xdat(3).x = ( ( m(11).m(1:len)  ) ); %difference in end and start vapour diags                                
        xdat(4).x =  ( m(9).m(1:len)  ) ; %difference in end and start vapour diags                                
        xdat(5).x = - ( ( m(10).m(1:len)  ) ); %difference in end and start vapour diags                                
        xdat(6).x = - ( ( m(1).m(1:len)  ) ); %difference in end and start vapour diags                                
    
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        %xdat(1).x(len-i+1) = xdat(1).x(len-i+1) + sumice ;    
       
    
    ydat(2).y=Grid.Z(1:len)/1000;
    labs(2).l='Vapour';
    
    ydat(1).y=Grid.Z(1:len)/1000;
    labs(1).l='Total Water';
    
    ydat(3).y=Grid.Z(1:len)/1000;
    labs(3).l='Microphysics';
    
    ydat(4).y=Grid.Z(1:len)/1000;
    labs(4).l='Fall Speed Flux';
    
    ydat(5).y=Grid.Z(1:len)/1000;
    labs(5).l='Ice Flux';
    
    ydat(6).y=Grid.Z(1:len)/1000;
    labs(6).l='Vapour Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    %xlims=1;
    xlimits=[-2e-4 6e-4];
    
    case 21
    logflag=0;
    
    xlab='Cumulative Net Gain in Vapour (kg/m^2)';
	ylab='Height (km)';
    
    figname='Vapour Gain End-Start';
    
    
    
    zmin=0;  %(km)
    zmax=z(end)/1000;
    
    len=length(icediag3(1).i(:,1,11))-1;
    
    dz=Grid.Z(2:end)-Grid.Z(1:end-1); 
    rho=Grid.RHO;
    
    for i=1:len
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( sum(icediag3(1).i(len-i+1:len,end,[2 5 8 11]),3) - sum(icediag3(1).i(len-i+1:len,1,[2 5 8 11]),3) ) );
        xdat(2).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        xdat(3).x(len-i+1) = sum( ( m(11).m(len-i+1:len)  ) ); %difference in end and start vapour diags                                
        xdat(4).x(len-i+1) = sum( ( m(9).m(len-i+1:len)  ) ); %difference in end and start vapour diags                                
        xdat(5).x(len-i+1) = - sum( ( m(10).m(len-i+1:len)  ) ); %difference in end and start vapour diags                                
        xdat(6).x(len-i+1) = - sum( ( m(1).m(len-i+1:len)  ) ); %difference in end and start vapour diags                                
        xdat(7).x(len-i+1) = - sum( ( m(19).m(len-i+1:len) + m(21).m(len-i+1:len)  ) ); 
        xdat(8).x(len-i+1) = - sum( ( m(20).m(len-i+1:len) + m(22).m(len-i+1:len) ) ); 
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                
        %xdat(1).x(len-i+1) = xdat(1).x(len-i+1) + sumice ;    
    end
    
    ydat(2).y=Grid.Z(1:len)/1000;
    labs(2).l='Vapour';
    
    ydat(1).y=Grid.Z(1:len)/1000;
    labs(1).l='Total Water';
    
    ydat(3).y=Grid.Z(1:len)/1000;
    labs(3).l='Microphysics';
    
    ydat(4).y=Grid.Z(1:len)/1000;
    labs(4).l='Fall Speed Flux';
    
    ydat(5).y=Grid.Z(1:len)/1000;
    labs(5).l='Ice Flux';
    
    ydat(6).y=Grid.Z(1:len)/1000;
    labs(6).l='Vapour Flux';
    
    ydat(7).y=Grid.Z(1:len)/1000;
    labs(7).l='Upwards Vapour Flux';
    
    ydat(8).y=Grid.Z(1:len)/1000;
    labs(8).l='Downwards Vapour Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    %xlims=1;
    xlimits=[-2e-4 6e-4];
    
    
    case 20
    logflag=0;
    
        
    len=findheight(z,18.2e3);
    len2=findheight(z,30e3);
    
    
    xlab='Cumulative mean mixing ratio down from 19km (ppmv)';
	ylab='Height (km)';
    
    figname='Cumulative mean mixing ratio down from 19km ';
    
    
    
    zmin=13;  %(km)
    zmax=19;

    len=length(z);
    
    for i=1:len
        
        
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        xdat(1).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,2)));
        %xdat(2).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,2)));
        
        xdat(2).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        %xdat(4).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        
        %m2= f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        
        %xdat(1).x(len-i+1)=( m2-xdat(1).x(len-i+1) ) / mean(vap(1).v(len-i+1,:,2));
        
        %xdat(2).x(len2-i+1) = f*mean(mean(vap(1).v(len2-i+1:len2,:,2)));
        %xdat(2).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        %m2= f*mean(mean(vap(1).v(len-i+1:len2,:,88)))  ;
        
        %xdat(2).x(len2-i+1)=( m2-xdat(2).x(len2-i+1) ) / mean(vap(1).v(len2-i+1,:,2)) ;
    end
    
      
    
    ydat(1).y=z(1:len)/1000;
    
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';

    labs(1).l='Start';
    
    ydat(2).y=z(1:len)/1000;
    %ydat(4).y=Grid.PREFN(1:len)/100;
    labs(2).l='End';
    
    dual=1;
    
    xlims=1;
    xlimits=[4.5 6];
    
    case 19
    logflag=0;
    
    
    
    figname='Cumulative mean mixing ratio down from 19km ';
    
    
    
    zmin=13.5;  %(km)
    zmax=30;
    
    %xlims=1;
    xlimits=[4.4 6.3];
    
    %len=findheight(z,z(end));
    len=findheight(z,19.5e3);
    len=length(z);
    
    hstr=num2str(round2(z(len)/1000,1));
    xlab=['Mean mixing ratio down from ',hstr,'km (ppmv)'];
	ylab='Height (km)';
    
    secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';
    dual=1;
    
    for i=1:len
        
        
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        %xdat(1).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,2)));
        %xdat(2).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        
        xdat(1).x(len-i+1) = f*mean(icediag3(1).i(len-i+1:len,1,11),1);
        xdat(2).x(len-i+1) = f*mean(icediag3(1).i(len-i+1:len,88,11),1);
        
        %m2= f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        
        %xdat(1).x(len-i+1)=( m2-xdat(1).x(len-i+1) ) / mean(vap(1).v(len-i+1,:,2));
        
        %xdat(2).x(len2-i+1) = f*mean(mean(vap(1).v(len2-i+1:len2,:,2)));
        %xdat(2).x(len-i+1) = f*mean(mean(vap(1).v(len-i+1:len,:,88)));
        %m2= f*mean(mean(vap(1).v(len-i+1:len2,:,88)))  ;
        
        %xdat(2).x(len2-i+1)=( m2-xdat(2).x(len2-i+1) ) / mean(vap(1).v(len2-i+1,:,2)) ;
    end
    
    xdat(1).x=xdat(2).x-xdat(1).x;
    xdat(2)=[];
    
    ydat(1).y=z(1:len)/1000;
    %ydat(1).y=Grid.PREFN(1:len)/100;

    labs(1).l='Start';
    
    %ydat(2).y=z(1:len)/1000;
    %ydat(2).y=Grid.PREFN(1:len)/100;
    %labs(2).l='End';
 
    
    case 18
    logflag=0;
    
    xlab='Mean Vapour (ppmv)';
	ylab='Height (km)';
    
    figname='Vapour Gain End-Start';
      
    zmin=13.5;  %(km)
    zmax=21;
    
    xlims=1;
    xlimits=[3 8];
    
    j=1;
    %xdat(j).x=f*mean(vap(1).v(:,:,2),2); 
    xdat(j).x=f*icediag3(1).i(:,1,11);
    ydat(j).y=z/1000;
    labs(j).l='Mean Vapour at Start (ppmv)';
    
    j=2;
    %xdat(j).x=f*mean(vap(1).v(:,:,88),2);   
    xdat(j).x=f*icediag3(1).i(:,88,11);
    ydat(j).y=z/1000;
    labs(j).l='Mean Vapour at End (ppmv)';
    
    
%     j=3;
%     xdat(j).x=f*mean(vap(1).v(:,:,10),2);   
%     ydat(j).y=z/1000;
%     labs(j).l='Mean Vapour at Dump 10 (ppmv)';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';

	secyA=z/1000;
    secyB=Grid.PREFN/100;
    lab2='Pressure (hPa)';
    dual=1;
    
    
    lor=2;
    
    case 17
    logflag=0;
    
    t1=19.75;
    t2=25.167;
    [it1,it2,t1str,t2str]=time_strings(t1,t2,GridDan(idir).t+3)
       
    z=GridDan(idir).Z+620;
    
    
    %len=length(icediag4(idir).i(:,1,35))-1;
    len=findheight(z,18.5e3)+1;
    %len=length(z)-1;
    
    hstr=num2str( round2(z(len)/1000,1) );
%    xlab=['Averaged Net Gain Down from ' hstr 'km (kg/m^{2})'];
    xlab=['Averaged Net Gain Down from ' hstr 'km (ppmv)'];

	ylab='Height (km)';
    
    
    
    figname=['Cumulative net gain down from ' hstr 'km'];
    
    dz=GridDan(idir).Z(2:end)-GridDan(idir).Z(1:end-1); %finds d/dz of flux
    rho=GridDan(idir).RHON;
    
dirs=[1 2];  
  for idat=1:length(dirs)   
		idir=dirs(idat);
      
    for i=1:len
        %me=mean(icediag3(1).i(len-i+1:len,end,11));
        %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                 
        
        airmass=sum( rho(len-i+1:len) .*dz(len-i+1:len) );
        %airmass=1;

        %converted into mixing ratio
        xdat(idat).x(len-i+1) = f * sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediagsALL(idir).i(len-i+1:len,it2,[37:42]),3) - sum(icediagsALL(idir).i(len-i+1:len,1,[37:42]),3) )  ) ./ airmass;
       % xdat(2).x(len-i+1) = sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediag4(idir).i(len-i+1:len,end,[35:36]),3) - sum(icediag4(idir).i(len-i+1:len,3,[35:36]),3) )  ) ./ airmass;
      %  xdat(1).x(len-i+1) = f * (sum(icediagsALL(idir).i(len-i+1,it2,[37:42]),3) - sum(icediagsALL(idir).i(len-i+1,1,[37:42]),3) );

        
        %mass gain
        %xdat(1).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediag3(1).i(len-i+1:len,end,[2 5 8 11]),3) - sum(icediag3(1).i(len-i+1:len,1,[2 5 8 11]),3) )  ) ;
        %xdat(2).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ) ; %difference in end and start vapour diags    
    end
        
    ydat(idat).y=z(1:len)/1000;
    labs(idat).l=runName(idir).nam;
    
  end
    

    
    secyA=z/1000;
    secyB=GridDan(idir).PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;
    
    xlims=0;
    %xlimits=[-1.4 0.4];
    xlimits=[-400/f 400/f];
    
    zmin=15;  %(km)
    %zmax=z(end)/1000;
    zmax=19;
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    case 16
    logflag=0;
    
    xlab='Net Gain in Vapour (kg/m^2)';
	ylab='Height (km)';
    
    figname='Vapour Gain End-Start';
    
    zmin=0;  %(km)
    zmax=30;
    
    
    
    xdat(1).x =f*( icediag3(1).i(:,end,11) - icediag3(1).i(:,1,11) ); %difference in end and start vapour diags                                
    ydat(1).y=Grid.Z/1000;
    labs(1).l='Net Vapour Gain Profiles';
    
    xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
    [iz,iz2]=findheight(Grid.Z,0e3,30e3);
    ydat(2).y=Grid.Z(iz+1:iz2)/1000;
    labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=2;
    
    
    case 15
    logflag=0;
    
    xlab='Net Gain in Vapour (kg/m^2)';
	ylab='Height (km)';
    
    figname='Vapour Gain';
    
    zmin=0;  %(km)
    zmax=30;
    
    
    
    xdat(1).x =f*( -m(1).m + m(11).m ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
    [iz,iz2]=findheight(Grid.Z,0e3,30e3);
    ydat(1).y=Grid.Z(iz+1:iz2)/1000;
    labs(1).l='Net Vapour Gain';
    
    lor=2;
    
    case 14
    
    logflag=1;
    
    xlab='Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    
    tt=53;
    figname=strcat('Profiles around clouds - High Updraught Case - Time=',num2str(timesTH(tt-49)),' LT','-dump ',num2str(tt) );
    
    j=1;
    xdat(j).x=sdla(7,:);
    ydat(j).y=sdla(2,:)/1000;
    labs(j).l='SDLA SF4 descent 18:57-21:48 LT';
    
    j=2;
    xdat(j).x=f*pr(1).p(:,10);
    ydat(j).y=pr(1).p(:,1)/1000;
    labs(j).l='DMI 17:15-19:08 LT';
    
    j=3;
    xx=257;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
    

    j=4;
    xx=267;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
    
    j=5;
    xx=247;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
    
     j=6;
    xx=260;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
    
    j=7;
    xx=254;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
    
     j=8;
    xx=251;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
    
    j=9;
    xx=248;
    xdat(j).x=f*vap(2).v(:,xx,tt);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l=strcat('x =',' ',num2str(Grid.Y1(xx)/1000),'km');
  
    
    case 13
    logflag=0;
    
    xlab='Total Water Mass from Flux in 4.5 hrs (kg/m^2)';
	ylab='Height (km)';
    
    figname='Total Water from Flux';
    
    izmin=97;  %115;
    izmax=136;
    dumprange=[50:104];
    
    i=2;
    pdat(i).p =   squeeze(Fluxdiag(i).dg(izmin:izmax,6,dumprange) +Fluxdiag(i).dg(izmin:izmax,6+14,dumprange) - Falldiag(i).dg(izmin:izmax,6,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,4,dumprange) +Fluxdiag(i).dg(izmin:izmax,4+14,dumprange) - Falldiag(i).dg(izmin:izmax,4,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,5,dumprange) +Fluxdiag(i).dg(izmin:izmax,5+14,dumprange) - Falldiag(i).dg(izmin:izmax,5,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,1,dumprange) +Fluxdiag(i).dg(izmin:izmax,1+14,dumprange) - Falldiag(i).dg(izmin:izmax,1,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,2,dumprange) +Fluxdiag(i).dg(izmin:izmax,2+14,dumprange) - Falldiag(i).dg(izmin:izmax,2,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,3,dumprange) +Fluxdiag(i).dg(izmin:izmax,3+14,dumprange) - Falldiag(i).dg(izmin:izmax,3,dumprange));



    
            
    j=1;
    
    %diff=pdat(2).p(2:end,:)-pdat(2).p(1:end-1,:);
    %xdat(j).x=sum(diff(:,1:length(dumprange)),2)*300*length(dumprange);
    xdat(j).x=sum(pdat(2).p(:,1:length(dumprange)),2)*300*length(dumprange);
    ydat(j).y=Grid.Z(izmin:izmax)/1000;
    labs(j).l='High Updraught Case';
    
    zmin=Grid.Z(izmin);
    
    lor=2;
    
    
	case 12
    logflag=1;
    
    xlab='Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='Beginning and end mean profiles';
    
    
    
    j=1;
    xdat(j).x=mean(vap(2).v(:,:,50),2);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l='High Updraught Beginning Mean Vapour';
    

    j=2;
    xdat(j).x=mean(vap(2).v(:,:,end),2);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l='High Updraught End Mean Vapour';
    
%     j=3;
%     xdat(j).x=pcents_icemr(2).p(end,2:end,end);
%     ydat(j).y=Grid.Z(2:end)/1000;
%     labs(j).l='High Updraught End Max Ice Sat';
%     
%     j=4;
%     xdat(j).x=f*mean(vap(2).v(:,:,50),2);
%     ydat(j).y=Grid.Z(1:end)/1000;
%     labs(j).l='High Updraught Beginning Vapour';
    
    case 11
	logflag=1;
    
    xlab='Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='Beginning and end MAX and MIN profiles';
    
    
    
    
    xdat(1).x=pcents_icemr(2).p(50,2:end,1);
    ydat(1).y=Grid.Z(2:end)/1000;
    labs(1).l='High Updraught Beginning Ice Sat';
    

    j=2;
    xdat(j).x=pcents_icemr(2).p(end,2:end,1);
    ydat(j).y=Grid.Z(2:end)/1000;
    labs(j).l='High Updraught End Min Ice Sat';
    
    j=3;
    xdat(j).x=pcents_icemr(2).p(end,2:end,end);
    ydat(j).y=Grid.Z(2:end)/1000;
    labs(j).l='High Updraught End Max Ice Sat';
    
    j=4;
    xdat(j).x=f*mean(vap(2).v(:,:,50),2);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l='High Updraught Beginning Vapour';
    

 
    
    case 10
          
        
    for i=1:2
        xdat((i-1)*2+1).x=f*max(TwoDDan(i).Q(:,2:end,1),[],2);
        ydat((i-1)*2+1).y=GridDan(i).Z/1000;
     
        xdat((i-1)*2+2).x=f*min(TwoDDan(i).Q(:,2:end,1),[],2);
        ydat((i-1)*2+2).y=GridDan(i).Z/1000;
        
    end
    
    T=tempLES(GridDan(1)); %K
    ei=SatVapPress(T,'goff','ice'); %Pa
    P=GridDan(1).PREFN; %Pa
    
    xdat(5).x=f*0.622*ei./(P-ei);
    ydat(5).y=GridDan(1).Z/1000;
    
    T=tempLES(GridDan(2)); %K
    ei=SatVapPress(T,'goff','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    xdat(6).x=f*0.622*ei./(P-ei);
    ydat(6).y=GridDan(2).Z/1000;
    
    figname='Max/Min Vapour + Ice Sat';
        %ydat(1).y=(sumPosDep-sumPosUnDep)./sumPosDep;    
    
	labs(1).l='Max for Low Updraught Case';
   	labs(2).l='Min for Low Updraught Case';
	labs(3).l='Max for High Updraught Case';
	labs(4).l='Min for High Updraught Case';
    labs(5).l='Ice Sat Mixing Ratio Low Updraught';
    labs(6).l='Ice Sat Mixing Ratio High Updraught';
	
	xlab='Water Vapour Mixing Ratio (ppmv)';
    %xlab='Aerosol mass (kg)';
	ylab='Height (km)';
    
    
    logflag=1;
    lor=1;
    gridon=1;
    
    %set(gca,'xlim',[9e-9 2e-6]); do this after


    
    case 9
        
  
     xlab='Water Vapour Mixing Ratio (ppmv)';
     ylab='Height (km)';
     logflag=1;
     lor=1;
     
%      pcplot3(1,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'Low Updraught'},f);
%      
%      pcplot3(2,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'High Updraught'},f);
    
    pcs=[0 25 50 75 100];    

    [xdat,ydat,labs]=getpctdat(GridDan(2),f*TwoDDan(2).Q(:,2:end,1),pcs,' high updraught');    

    j=length(pcs)+1;
    
    
     xdat(j).x=f*pr(1).p(:,10);
     ydat(j).y=pr(1).p(:,1)/1000;
     labs(j).l='Bauru Sounding 17:15';
     
     xdat(j+1).x=f*pr(2).p(:,10);
     ydat(j+1).y=pr(2).p(:,1)/1000;
     labs(j+1).l='Campo Grande Sounding 9am';

     xdat(j+2).x=origVap;
     ydat(j+2).y=GridDan(1).Z/1000;
     labs(j+2).l='Original Sounding';

    case 8
    
    logflag=1;
    figname='Mean Percentiles Vapour Graph 50-104';
    xlab='Water Vapour Mixing Ratio (ppmv)';
    
    for i=50:size(vap(1).v,3)
        
        pcs=[0 25 50 75 100];
        
        for j=1:length(pcs)
            pcents(2).p(i,:,j)=prctile(f*vap(2).v(:,:,i)',pcs(j));
        end
    end
    
    for j=1:length(pcs)
        xdat(j).x=mean(pcents(2).p(:,:,j),1);
        ydat(j).y=GridDan(1).Z/1000;
    end
        
    
    lp=length(pcs);
    for j=1:lp
    
        if pcs(j)==0
            labstr='Min';
        elseif pcs(j)==100
            labstr='Max';
        else
            labstr=strcat( num2str(pcs(j)) , 'th percentile' );
        end
    labs(j).l=strcat(labstr,' for ',' high updraught case');
    end
    
    j=length(pcs)+1;
    xdat(j).x=origVap;
    ydat(j).y=GridDan(1).Z/1000;
    labs(j).l='Original Sounding';
    
     
    case 7
        
    logflag=1;

    xlab='Water Vapour Mixing Ratio (ppmv)';
    ylab='Height (km)';
    
%     for i=50:size(vap(2).v,3)
%         
%         pcs=[0 25 50 75 100];
%         
%         for j=1:length(pcs)
%             pc(i,:,j)=prctile(f*vap(2).v(:,:,i)',pcs(j));
%         end
%     end
    
    pcs=[0 25 50 75 100];

    for j=1:length(pcs)
        xdat(j).x=mean(pcents(1).p(50:104,:,j),1);
        ydat(j).y=620/1000+GridDan(1).Z/1000;
    end
        
    
    lp=length(pcs);
    for j=1:lp
    
        if pcs(j)==0
            labstr='Min';
        elseif pcs(j)==100
            labstr='Max';
        else
            labstr=strcat( num2str(pcs(j)) , 'th percentile' );
        end
    labs(j).l=strcat(labstr,' for ',' high updraught case');
    end
    
    j=length(pcs)+1;
    xdat(j).x=origorig;
    ydat(j).y=620/1000+GridDan(1).Z/1000;
    labs(j).l='Dump 1';
    
    j=length(pcs)+2;
    xdat(j).x=min(pcents(1).p(50:104,:,1),[],1); %minimum value over all times
    ydat(j).y=620/1000+GridDan(1).Z/1000;
    labs(j).l='Overall Min';
    
    j=length(pcs)+3;
    xdat(j).x=max(pcents(1).p(50:104,:,5),[],1); %max value over all times
    ydat(j).y=620/1000+GridDan(1).Z/1000;
    labs(j).l='Overall Max';
    

    
    
    case 6
        
    logflag=1;

    xlab='Water Vapour Mixing Ratio (ppmv)';
    
   
        xdat(1).x=f*mean(TwoDDan(1).Q(:,:,1),2);
        ydat(1).y=GridDan(1).Z/1000;
    
        labs(1).l='Start Water Vapour Profile';
        
        figname='Start Vapour Profile';
    
        
    case 5
     xlab='Water Vapour Mixing Ratio (ppmv)';
     ylab='Height (km)';
     logflag=1;
     lor=1;
     
%      pcplot3(1,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'Low Updraught'},f);
%      
%      pcplot3(2,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'High Updraught'},f);
    
    %pcs=[0 25 50 75 100];    

    %[xdat,ydat,labs]=getpctdat(GridDan(2),f*TwoDDan(2).Q(:,2:end,1),pcs,' high updraught');    

    stats={'Corumba','CG','SP','Curitiba','Bauru'};

    for j=1:length(pro)
       
     xdat(j).x=f*pro(j).p(:,10,1);
     ydat(j).y=pro(j).p(:,1,1)/1000;
     labs(j).l=stats{j};
     
    end
    
    j=length(pro)+1;
%     xdat(j).x=f*mean(vap(1).v(:,:,50),2);
    xdat(j).x=origorig;
    ydat(j).y=620/1000 + GridDan(1).Z/1000; %Bauru 620m above msl
    labs(j).l='dump 1';
    
   
    
    case 4
        
    logflag=1;

    xlab='Water Vapour Mixing Ratio (ppmv)';
    
   
        
    figname='Dump 50 percentiles';
    
    pcs=[0 25 50 75 100];    

    [xdat,ydat,labs]=getpctdat(GridDan(2),f*vap(1).v(:,:,50),pcs,' high updraught');   

    
    
    
    
%     j=length(pro)+1;
%     xdat(j).x=f*mean(vap(1).v(:,:,50),2);
%     ydat(j).y=GridDan(1).Z/1000;
%     labs(j).l='dump 50';

case 399
    %SF4 vapour plots
    %run readsound and set dmi=pr 
    %and get LEM data from loadvapdata
    
    lor=3;
    
    zmin=11;
    zmax=21;
    
    izlim=1;
	
	xlims=1;
	xlimits=[0 20];
    
    
    logflag=1;
    
    xdat(1).x=f*dmi(1).p(:,10);
    ydat(1).y=dmi(1).p(:,1)/1000;
    labs(1).l='Original Sounding 20:15-22:08 UTC';
    
    T=dmi(1).p(:,3)+273.15; %K
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=dmi(1).p(:,2)*100; %Pa
    
    xdat(2).x=f*0.622*ei./(P-ei);
    ydat(2).y=dmi(1).p(:,1)/1000;
    labs(2).l='Ice Saturation';  
    
     xdat(3).x=f*GridDan(1).OLQBAR(:,1);
     ydat(3).y=GridDan(1).Z/1000 + 0.62;
     labs(3).l='LEM water vapour';
    
  
    xlab='Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='24th Feb vapour and ice saturation mixing ratios';
    titlenam=figname;
    savename=figname;
    
    
case 3
    %SF4 vapour plots
    %run readsound and set dmi=pr, readSDLA, readallSAW
    %or get from loadvapdata
    
    lor=3;
    
    zmin=11;
    zmax=21;
    
    izlim=1;
	
	xlims=0;
	xlimits=[0 20];
    
    
    zmin=1;
    zmax=21;
    
    
    logflag=0;
    
    xdat(1).x=sdla(1).s(7,:);
    ydat(1).y=sdla(1).s(2,:)/1000;
    labs(1).l='SDLA 21:57-00:48 UTC';
    
   
    
    
    T=sdla(1).s(6,:)+273.15; %K
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=sdla(1).s(5,:)*100; %Pa
    
    xdat(2).x=f*0.622*ei./(P-ei);
    ydat(2).y=sdla(1).s(2,:)/1000;
    labs(2).l='SDLA Ice Sat';
    
    
%     T=tempLES(GridDan(2)); %K
%     ei=SatVapPress(T,'buck2','ice'); %Pa
%     P=GridDan(2).PREFN; %Pa
%     
%     xdat(3).x=f*0.622*ei./(P-ei);
%     ydat(3).y=GridDan(2).Z/1000;
%     labs(3).l='Ice Sat Mixing Ratio High Updraught';
    
  
    
%     
%     T=tempLES(GridDan(2)); %K
%     ei=SatVapPress(T,'goff','ice'); %Pa
%     P=GridDan(2).PREFN; %Pa
%     
%     xdat(4).x=f*0.622*ei./(P-ei);
    

    xdat(3).x=f*dmi(1).p(:,10);
    ydat(3).y=dmi(1).p(:,1)/1000;
    labs(3).l='DMI 20:15-22:08 UTC';
    
    T=dmi(1).p(:,3)+273.15; %K
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=dmi(1).p(:,2)*100; %Pa
    
    xdat(4).x=f*0.622*ei./(P-ei);
    ydat(4).y=dmi(1).p(:,1)/1000;
    labs(4).l='DMI Ice Sat';
    
%     xdat(5).x=saw(1).s(4).sss(7,:);
%     ydat(5).y=saw(1).s(4).sss(3,:)/1000;
%     labs(5).l='SAW water vapour';
    
    %xdat(6).x=f*origorig;
    %ydat(6).y=620/1000 + GridDan(1).Z/1000; %Bauru 620m above msl
    %labs(6).l='LEM start profile';
    
     xdat(5).x=f*GridDan(1).OLQBAR(:,1);
     ydat(5).y=GridDan(1).Z/1000 + 0.62;
     labs(5).l='LEM water vapour';
    
  
    xlab='Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='24th Feb vapour and ice saturation mixing ratios';
    titlenam=figname;
    savename=figname;
    
    
case 2
    logflag=1;
    
    xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='LEM anvil and outside profiles dump 85';
    
    x=3;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=1;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    x=20;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=2;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    x=40;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=3;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    x=60;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=4;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    

case 1
	logflag=1;
    
    xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='Beginning and end min profiles';
    
    
    xdat(1).x=min(vap(2).v(:,:,50),[],2);
    ydat(1).y=Grid.Z/1000;
    labs(1).l='High Updraught Beginning';
    
    j=2;
    xdat(j).x=min(vap(2).v(:,:,end),[],2);
    ydat(j).y=Grid.Z/1000;
    labs(j).l='High Updraught End';
    

end


end

if izlim==0
    zmin='';
    zmax='';
end

if noplot==0
	if subplotting==0    
        if ~exist('ioverride_watervap_newfig') | ioverride_watervap_newfig==0
            scrsz=get(0,'ScreenSize');
            
            if ~exist('ifull_screen') | ifull_screen==0
                posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];
            else
                clear ifull_screen
                posit=[1 -scrsz(2)*0.95 scrsz(3) scrsz(4)];
            end
            
            hf=figure('name',figname,'Position',posit,'color','w');
            
        end
        
        set(gcf,'papersize',[1 9]);

            
%        fsize=18;
%        fsize=12;
%        fsize=26;
%        fsize=30;  

        ixlab=1;
    else  %if subplotting==0
        if ~exist('idirs')
            idirs = [1:xsub*ysub];
        end
        if nsub==1
            %posit is the screen dimensions
            posit(4)=length(idirs)*posit(4)/2.3;
            posit(4)=length(idirs)*posit(4)/1;            
			hf=figure('name',figname,'Position',posit);
        end
        subplot(xsub,ysub,nsub);
        fsize=12;
        fsize=18;
        if nsub==length(idirs)
            ixlab=1;
        else
            ixlab=1; %Overridding this for now
        end
    end
    
    
    for idat=1:length(xdat)
        if ismooth_x(idat)==1
            [ydat(idat).y,xdat(idat).x] = window_average(ydat(idat).y,xdat(idat).x,Nsmooth_window,smooth_mode);
        elseif ismooth_y(idat)==1
            [xdat(idat).x,ydat(idat).y] = window_average(xdat(idat).x,ydat(idat).y,Nsmooth_window,smooth_mode);
        end
            
    end
           

zline=0; %to plot a solid line at x=0
if iplot_3D==0
    [h,ax,ax2]=plotXY6(xdat,ydat,labs,nmark,lwidth,lor,logflag,xlab,ylab,[zmin zmax],...
    zline,dual,secyA,secyB,lab2,fsize,ixlab,ixdir,iydir,xloc,time_highlight_path,highlight_type,ierror_bars,errordatU,errordatL,marksize,ichoose_styles,line_pattern,line_colour,marker_style,line_widths,iovr_leg_line);
%function [H1,ax2]=plotXY3(xdat,ydat,labs,nmark,lwidth,leglor,logflag,xlab,ylab,ylims,zline) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers
%put nmark as -1 for markers for all points
% LEGEND(...,Pos) places the legend in the specified
%     location:
%         0 = Automatic "best" placement (least conflict with data)
%         1 = Upper right-hand corner (default)
%         2 = Upper left-hand corner
%         3 = Lower left-hand corner
%         4 = Lower right-hand corner
%        -1 = To the right of the plot
else
    [h,ax,ax2]=plotXY6_3D(xdat,yydat,ydat,labs,nmark,lwidth,lor,logflag,xlab,ylab,yylab,[zmin zmax],...
    zline,dual,secyA,secyB,lab2,fsize,ixlab,ixdir,iydir,xloc,time_highlight_path,highlight_type,iovr_leg_line);


    if exist('iset_3D_view_properties')==1 & iset_3D_view_properties==1
        set(gca,'View',View);
        set(gca,'CameraPosition', CameraPosition);
        set(gca,'CameraPositionMode', CameraPositionMode);
        set(gca,'CameraTarget', CameraTarget);
        set(gca,'CameraTargetMode', CameraTargetMode);
        set(gca,'CameraUpVector', CameraUpVector);
        set(gca,'CameraUpVectorMode', CameraUpVectorMode);
        set(gca,'CameraViewAngle', CameraViewAngle);
        set(gca,'CameraViewAngleMode', CameraViewAngleMode);
        
    end
        
        
end

set(gca, 'layer', 'top'); %this puts the tick marks on top and so over the 
%timeseries highlighting box - make sure to do this before grid on

YL=textwrap({ylab},40);
ylabel(YL);

if xlims==1
    set(gca,'xlim',xlimits);
end



if gridon==1
    grid on;
end



if (ixtime==1)
    xx=get(gca,'xticklabels');
    xx=str2num(xx);
    xx=num2str(mod(xx,24));
    set(gca,'xticklabels',xx);
end



if idirlabel==1
    xlims=get(gca,'xlim');
    ylims=get(gca,'ylim');
    text(xlims(1),ylims(1)-(ylims(2)-ylims(1))/13,direcDan(idir).dir);
end

if add_points==1
    add_points_to_plot(xpos,ypos,point_labs,8,11)
end

if iaxis_square==1;
    if dual==2
        axes(ax2);
        axis square
        axes(ax);
        axis square
    end
        
    axis square
elseif iplot_3D==1
    dasp=daspect;
    daspect([1 1 dasp(3)]); 
end
    
switch x_axis_type        
    case {'log10_matlab'} %log10_matlab is when the 'yscale' is set to 'log' using matlab built in feature
        set(gca,'xscale','log');
end

switch y_axis_type        
    case {'log10_matlab'} %log10_matlab is when the 'yscale' is set to 'log' using matlab built in feature
        set(gca,'yscale','log');
end

if iadd_line==1
    line(addlineX,addlineY);
end

if ixtick_relabel_log==1
    ticklabels = get(gca,'xticklabel');
    tick_nums = str2num(ticklabels);
    
    switch x_axis_type        
        case {'log10','log10_matlab'}
            new_tick_nums = 10.^tick_nums;            
    end
    
    clear new_xtick_text
    
    for itick=1:length(new_tick_nums)
        xtick_text_i = num2str(new_tick_nums(itick));
        new_xtick_text(itick,1:length(xtick_text_i)) = xtick_text_i;
    end
    
    set(gca,'xticklabel',new_xtick_text);
    
end

if iytick_relabel_log==1
    ticklabels = get(gca,'yticklabel');
    tick_nums = str2num(ticklabels);
    
    switch y_axis_type        
        case {'log10','log10_matlab'}
            new_tick_nums = 10.^tick_nums;            
    end
    
    clear new_ytick_text
    
    for itick=1:length(new_tick_nums)
        ytick_text_i = num2str(new_tick_nums(itick));
        new_ytick_text(itick,1:length(ytick_text_i)) = ytick_text_i;
    end
    
    set(gca,'yticklabel',new_ytick_text);
    
end

if iset_xticks == 1
    set(gca,'xtick',xtickvals_set);
end
if iset_xticklabs == 1
    set(gca,'xticklabel',xticklabs_set);
end



iadd_Temperature_yaxis=0;
if iadd_Temperature_yaxis==1
    ax1=gca;
    ylims_1 = get(ax1,'ylim');
    dylim_1=ylims_1(2)-ylims_1(1);

    ax1_pos=get(ax1,'Position');
    ax2_pos=[ax1_pos(1)-ax1_pos(3)*0.14 0.1100 ax1_pos(3)/1000 0.8150];
    ax2 = axes('Position',ax2_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k');

    Z=dat_flt(:,col_alt);
    T=dat_flt(:,col_temp);
    P=dat_flt(:,col_press);    

    zgrid=[0:10:max(Z)];
    for iz=1:length(zgrid)-1
        ii=find( Z>zgrid(iz) & Z<zgrid(iz+1) );
        tgrid(iz)=mean(T(ii));
        pgrid(iz)=mean(P(ii));
        if length(ii)==0
            tgrid(iz)=NaN;
        end
    end

    zgrid=zgrid(2:end);

    set(ax2,'ylim',[floor(min(tgrid)) ceil(max(tgrid))]);
    set(ax2,'ydir','reverse');
    ylims_2 = get(ax2,'ylim');
    dylim_2=ylims_2(2)-ylims_2(1);

    yticks=get(ax2,'ytick');
    yticks=fliplr(yticks);


    clear ytickstr
    for jticks=1:length(yticks)
        z_tick = dylim_1*(yticks(jticks)-ylims_2(1))/dylim_2;
        new_tick = interp1(zgrid,tgrid,z_tick);
        if ~isnan(new_tick)
%            te=num2str( sigfig(new_tick,2+ceil(log10(abs(new_tick)))) );            
            te=num2str( sigfig(new_tick,2) );  
            ytickstr(jticks,1:length(te))=te;
        end

    end

    set(ax2,'yticklabel',ytickstr);









end



if i_set_dateticks==1
    set(gca,'xtick',date_ticks);     
end

 if idatetick==1; %flag to say the want the xaxis in proper time format rather than decimal time
        %specify the type with datetick_type (see help datetick)
%        datetickzoom('x',datetick_type,'keeplimits','keepticks');
        datetickzoom('x',datetick_type,'keeplimits');        
%        datetickzoom('x',15,'keeplimits');                
        %    zoomadaptivedateticks('on');
        %note need to plot time in days after day 0 of year 0
 end
    
fontsize_figure(gcf,gca,18);  

if (ititle==1)
%    title(titlenam,'fontsize',18,'verticalalignment','middle');
    titlenamw=textwrap({titlenam},50);
    htit=title(titlenamw,'fontsize',fsize_title); %,'verticalalignment','baseline');
    
    if dual==2
        ptit=get(htit,'position');
        set(htit,'position',[ptit(1) ptit(2)*1.015]);
    end
%        title(titlenam,'fontsize',18,'verticalalignment','top');

%    title(titlenam,'fontsize',18);
end


if iexecute_script==1
    eval(script_name);
end

if iadd_nums_above==1
    add_numbers_above_timeseries  %add numbers above the points of a timeseries
end

else
%    disp('****** WARNING noplot set to one - not plotting anything!! ********')
end

figname=remove_character(figname,'/','_');
figname=remove_character(figname,'\','_');
savename=[savedir figname];


if iwrite_text_dat==1
    for idat=1:length(xtext_dat)
        for itext=1:length(text_dat(idat).text)
            text(xtext_dat(idat).x(itext),ytext_dat(idat).y(itext),text_dat(idat).text(itext,:))
        end
    end
end

if ihighlight_points==1
     for idat=1:length(xdat_highlight)
%        plot(xdat_highlight(idat).x,ydat_highlight(idat).y,'bo','markersize',16);
        plot(xdat_highlight(idat).x,ydat_highlight(idat).y,'ro','markersize',12,'markerfacecolor','r');        
     end    
end

      

if ~exist('preserve_flags') | preserve_flags==0
   clear_flags_watervap
end
catch watervap_error
    clear_flags_watervap    
    rethrow(watervap_error)
end

