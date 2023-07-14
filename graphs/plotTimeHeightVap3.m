try %Catching errors so that flags are reset even if the program
    %aborts due to an error. Commands after the catch statement are executed
    %if an error now occurs. The flags are then reset and the error is
    %"rethrown"

    %for 2D PDFs search for -- iuseYnorm -- flag for normalisation along Y (or similarly X)
    %... or ndhistc  - main command is qh = ndHistc_run([X Y], Xbins, Ybin);
    
    
    
plotcase=50;
plotcase=65; %emm
plotcase=577;

idpcolor=0;
fsize_tit = 10; %if want a long title not for final figures
fsize_tit = 8; %larger font more suitable for figures


if exist('mincovOvr') 
    mincovOvr_save=mincovOvr; %save these to reset to later as are changed here if logflag==1
else
    mincovOvr_save=0;
end
if exist('maxcovOvr')
    maxcovOvr_save=maxcovOvr;
else
    maxcovOvr_save=1;
end
    

%clear h min iplot hback
clear min iplot hback

iabc=1; %flag to tell it to add (a), (b), (c), etc. to runName(i).nam

if ~exist('comp')
comp='uni';
comp='UWchallenger';

end

iadd_cont=0; %flag to say whether to draw an extra contour over the color plot (at value cont_val)
idiff_cont_dat=0;
iovride_conts=0;
no_set_ctick_str=1; %flag to prevent the colorbar tick labels from being chagned - causaes the values on the
%colorbar to be changed (to arong ones) at certain fontsizes!!! - set to 1
%unless being very careful to watch out for this!
isquare=0; %flag to make axes square
ione_to_one_line=0; %draw a one-to-one line
iadd_terrain=0;
iplot_latlon=0;
iadd_wind_quivers = 0;
iadd_overlay = 0; %flag to say that want to plot an x,y line on top defined by x_overlay,y_overlay
overlay_style=['''k-''']; %set this as a default use triple ''' to get it to = 'k-'
iadd_streamline_z0=0; %flag to add a streamline overlay (using iadd_overlay)
icont_neg=0;
idraw_streamlines=0;
right_side_extra_bits=1; %switch to put the extreme values in blocks on right hand side of the plot
dan_test=0; %test bit of script, 28th Oct 2009
iytick_relabel=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
pcolor_shading='interp'; %default shading to use for pcolor (only when icont==0)
plot_type='filled contour'; %default type of plot
colorbar_location = 'EastOutside';
iytick_override=0; %flag to override the default tick mark locations
i_only_major_log_ticks=1; %flag to say that only want to label e.g. 1, 10 ,100, etc on log10 scales
y_axis_type=''; %default normal y-axis (the feature for the x-axis has not been implemented yet)
idatetick=0; %flag to say whether the x-axis should be madeto a time format axis
ixtickmode_auto=0; %flag to set xtickmode to 'auto' - otherwise some of the dateticks are missing
%seems to do funny things though....
tick_direction='out';  %makes ticks point outwards from the edge of the plot
colorbar_loc='vert';
%iplot_mean_XY = 'x'; %this is set in pdf2D_plot_commands.m - plots the means of either x
%or y as dots

if ~exist('iplotselect'); iplotselect=0; end
if ~exist('bigcbar'); bigcbar=0; end
if ~exist('lememm'); lememm=0; end
if ~exist('iabc'); iabc=0; end
if ~exist('comp'); comp='uni'; end
if ~exist('icont_extra'); icont_extra=0; end
if ~exist('i3d'); i3d=0; end
if ~exist('wrap2d'); wrap2d=0; end

if iplotselect==1
    plotcase=plotcases(iplot);
end

clear zz timesTH diff dat_overlay

if ~exist('idir'); idir=1; end
if ~exist('ieps'); ieps=0; end
if ~exist('isamescale'); isamescale=0; end
if ~exist('subplotting'); subplotting=0; end
if ~exist('onexlabel'); onexlabel=0; end      %flag to make it so that only the bottom plot has the xlabel on it
if ~exist('iplot'); iplot=1; end 
if ~exist('ikeep_X_above_zero'); ikeep_X_above_zero=0; end 
if ~exist('ipost_plotTime_commands'); ipost_plotTime_commands=0; end 



if exist('npess2')
    if plotcase~=65
        npes=npess2(idir);
    end
end

%plots time height plots

ilem=1;
icont=1;

if ieps==1
    fsize=12;
elseif subplotting==1
    fsize=18;
else
    switch comp
        case 'laptop'
            fsize=14; %26 best for paper figures when more than one plot is being used
            fsize=18;
        case 'UWchallenger'
            fsize=18;
            fsize=12;
            fsize=16;            
    end
end

isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

icolmap=0; %flag to set to another colormap (defined by cmap)

iutc=1; %flag for time axis to be labelled as UTC (otherwise is labelled Local Time)
%add_ground_height=0.62; %height to add to the vertical axis to account for level of ground abv msl.
%add_ground_height=0; %height to add to the vertical axis to account for level of ground abv msl.
%add_ground_height=1.0; %height to add to the vertical axis to account for level of ground abv msl.

try
add_ground_height=add_ground_heights(idir).h;
catch
end

if exist('add_ground_height')==1
%    fprintf('\n**** add_ground_height = %f ********\n',add_ground_height);
end



minZ=0e3;
maxZ=25e3;  %19000;
maxtr=1.0;
%timesTH=Time;
hrstartles=18.67;
dumprange=[1:13];
%dumprange=[1:49];
%dumprange=[1:39];
%dumprange=[1:72];
%dumprange=[1:62]; %62
%dumprange=[1:14];
dumprange=[1:78];

%dumprange=[1:45];
dumprange=[1:21];
%dumprange=[1:40];

%dumprange=[1:44];
%dumprange=[1:15];

vectorf=0;

timesTH(1).t=(dumprange-1)*300/3600 + hrstartles;

jmax=5; %max no. plots on one screen

a1=1;
a2=2; %values for subplot(ai,b,i)

izmin=2;
izmax=2;
f=1e6*28.97/18;


notsame=0; %flag to plot each plot individually in terms of colour scale
offset=0;
clines=0; % when zero removes black contour lines
nocbar=0;
sig=2; %no. sig figs for contour values
nplots2d=1;
itimestamp=0;
manclab=0; %flag to set contour labels manually with mouse
idirstamp=0; %flag to write the directory of the results in the corner

normcbar=0;

dumpint=300; %dump interval in seconds

i2d=0;
izovr=0; %flag to say are setting own z axis
itimelab=0; %flag to say that x axis is time and so should be restricted to <24
figlab='2d contour plot';
iylim=0; %flag to override the setting of the y part of axis (e.g. so can have axis larger than that plotted)
ixlim=0;

nplots2d2d=1;

ncont=25;
clab=1; %flag to label contours
i=1;

phase=1;


clear max timestxt pdat minC maxC



% for i=1:length(times)
%     te=num2str(times(i),3);
% 	timestxt(i,1:length(te))=te;
% end
% tit(1).tit='Max of Low Level Tracer';
% tit(2).tit='Max of Low Level Tracer'

if ~exist('man_choose_plotTimeHeight_graph')
    logflag=0;
    dlogflag=0;
    noplot=0; %flag to say to just do the calculations and not to plot
    iminovr=zeros([1 10]);
    imaxovr=zeros([1 10]);
else
%   clear man_choose_plotTimeHeight_graph
end




switch(plotcase)
    
    
case 65
    emmTimH
case 64
    hms
case 63
    vap_potemp
case 62
    lem_min_temp
case 61
    mpc_min_temp
    
case 60
    mpc_tot_satmr
    
case 59
    fact=1e6*28.97/18;
    logflag=0;
    tit(1).tit=['Equilibrium Vapour Mixing Ratio for dehyd. below ' num2str(h0,'%2.1f')... 
            ' km (ppmv), nclouds= ' num2str(nclouds) ' x-third= ' num2str(x2) ' km'];
    %  tit(1).tit=['Equilibrium Vapour Mixing Ratio (ppmv)'];
    
    savename=tit(1).tit;
    
    nplots2d=1;
    
    clines=0; %makes black contour lines appear
    clab=0;
    
    i2d=3; %tells it are labelling own x axis
    xlabelstr='dN/dt - Number of overshoots per month';
    
    minZ=15.8e3;
    maxZ=17e3;
    ncont=15;
    
    imaxovr=0;
    maxcovOvr=5.6;
    
    iminovr=0;
    mincovOvr=5.6;
    
    sig=3;
    
    
    notsame=1;
    
    izovr=1; %flag to set own z axis
    
    ncont=25;
    
case 58
    dq_dehyd_10thSep2005
    
case 57   %*************************    2-D plots of LEM fields with ice no. *******************************
    fact=1e6*28.97/18;
    logflag=0;
    %tit(1).tit='Height Dependent Tracer';
    onexlabel=1;
    
    switch i57
    case 1
        tit(1).tit='Tot Water Mixing Ratio (ppmv)';
        %            tit(1).tit='Total Condensate (kg/kg)';
    case 2 
        tit(1).tit='Vapour Mixing Ratio (ppmv)';
    end
    
    
    tit(2).tit='Ice Number Concentration (mg^{-1})';
    
    nplots2d=1;
    
    clines=0; %makes black contour lines appear
    clab=0;
    
    i2d=2; %tells it to label x axis in km
    
    minZ=15e3;
    maxZ=22e3;
    ncont=15;
    
    %    minZ=0.1e3;
    %    maxZ=23e3;
    
    
    
    z=GridDan(idir).Z;
    
    notsame=1;
    
    itimestamp=1;
    
    imaxovr=[1 0];
    maxcovOvr=8;
    iminovr=[1 0];
    %mincovOvr=2.5;
    mincovOvr=1.0;
    
%%      Choose i577 plot case below.
case 577
    
    idirstamp=1;
    hrange=8;
    
    switch hrange
    case 1
        minZ=14e3;
        maxZ=23e3;
    case 2
        minZ=0e3;
        maxZ=53e3;
    case 3
        minZ=0.2e3;
        maxZ=19e3;
    case 4
        minZ=15e3;
        maxZ=17e3;    
    case 5
        minZ=13e3;
        maxZ=19e3; 
    case 6
        minZ=14e3;
        maxZ=30e3;  
    case 7
        minZ=0e3;
        maxZ=30e3;      
    case 8
        minZ=0.2e3;
        maxZ=18e3;   
    case 9
        minZ=0e3;
        maxZ=15e3;           
    case 10
        minZ=0e3;
        maxZ=18e3;  
    case 11
        minZ=0e3;
        maxZ=4e3;          
    end
    %ncont=15;
    
    
    fact=1e6*28.97/18;
    itimestamp=1;
    
    if exist('SER')
        time1=SER(end,1)/3600 + 19.75;
        time1=SER(end,1)/3600;
        
        % time1=GridDan(1).t(it)+3;
        
        mins=(time1-floor(time1))*60;
        if round(mins)==60; mins=0; time1=time1+1; end
        minstr=num2str(mins,'%2.0f');
        
        hrs=mod(floor(time1),24);
        hrstr=num2str(hrs,'%2.0f');
        if strcmp(minstr,'0'); minstr='00';end
        if length(hrstr)==1; hrstr=['0' hrstr];end
        if length(minstr)==1;minstr=['0' minstr];end
        timlab=[hrstr ':' minstr];
        
        
    end
    
    if ~exist('man_choose_plotTimeHeight_graph')
            i577='potemp';  %%%% vertical cross section of potemp - for e.g. 3D and 2D plot (Fig. 8 in paper) first of all 
            % load in the ThreeD.TH1 data for the 3D case. Then run wrap slice with iwrap=0 to take the vert slice and put it
            % in TwoDDan(1).TH2. Then load in the 2D data for the 2D case. Then run mulitsaveplot with i3d=1 and i2d=0 (ensures
            % that is wrapped) for the first run and i2d=1, i3d=0 for the second.
            
%     i577='lowtracer';
    %        i577='totwater';
%           i577='vapour';
     i577='general';
%    i577='si'; %supersat wrt ice
%      i577='temppert';
    % % %    i577='rhopert577';
    % %   %  i577='ozone';
    % %  % i577='hydbal';
    % %   i577='icesatMR';
    % %   i577='cdensity';
    %   
    %  
    % % i577='dpdz';
    %  % i577='rhog';
    % 
    %  i577='rad';
    % %i577='lnb';
   % i577='vertvel';
    %i577='inc';
    % 
    % %i577='ARM_radar';
    % 
%     i577='radar';
%    % %i577='w_3d';
%    i577='vap_3d';    %horizontal cross section
%     i577='vap_3d_vert';  %vertical cross sections using data in slice (obtain slice using wrap_slice.m)
    % 
    % %i577='ecmwf_surf';
%     i577 = 'wrf_wind2d';

%   -------- usual WRF plot ------------
     i577 = 'wrf_plot';
%   -----------------------------------

%    i577 = 'wrf_radar_vert';
%  i577 = 'wrf_vert_cross_section';
%  i577 = 'Particle size dist vs time';
%     i577 = 'MODIS_plot';

%   -------- usual UW plot ------------
      i577 = 'MODIS_plot_UW'; 
%   -----------------------------------

%      i577 = 'GCM_longitude_transect_UW';      
    
    else
%        clear man_choose_plotTimeHeight_graph
    end
    
%    i3d=0;
    
    

    
    clines=1; %makes black contour lines appear
    clab=0;
    
    ilem=0;
    
    
    i2d=2; %tells it to label x axis in km
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch i577
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'GCM_longitude_transect_UW'
        iadd_terrain=0;
        iplot_latlon=0;
        iadd_wind_quivers=0;
        iadd_flight_path=0;
        idraw_streamlines=0;
        clab=0;
        manclab=0;
        clines=0;
        icont=0; %flag to say whether want a contour or a colour plot
        pcolor_shading='flat';
        %pcolor_shading='interp';
        plot_type='pcolor';
%        plot_type='filled contour';
        

%%%% NOTE - these may be overrided in the individual plot cases below %%%%%
        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
        dlogflag=0;
        
        mincovOvr = 1e-2;  %usual one for CAS plots
        maxcovOvr = 10e3;

        mincovOvr = 1e-5;  
        maxcovOvr = 1;

        imaxovr=1;
        iminovr=1;

%        mincovOvr = 0;
%        maxcovOvr = 1000; %

        iovride_conts=0;
        conts_ovr=[270:2:310];
        conts_ovr=[270:1:300];
        %        conts_ovr=[0:0.5:19];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
        y_axis_type='log10';
        y_axis_type='log10_matlab'; %matlab's built in log scale - probably the best one to choose
        %as it does the minor tick marks between the major log numbers and puts ticks at appropriate points
        %also sets iytick_relabel to one as Matlab labels the ticks as 10^0 and 10^1 (rather than 1 and 10).
        y_axis_type='';
        %        y_axis_type='Psep';      

        iylim=0; %note is set below in certain cases
        ylims=[0 25];       

        izovr=2;
        i2d=3;

        idirstamp=0;
        itimestamp=0;
        sig=3;
        

        lon_csection_case='CFAD longitude cross section'; %run waterVapourMay2005 first
        lon_csection_case='Model field longitude cross section'; %run make_z_lon_slice first (CPT GCM model)
        
%        short_plot_name=lon_csection_case;
        short_plot_name=[titlenam] %from waterVapourMay2005

        %%% set some default values, so don't need to bother below
        %%% every time
        logflag=0; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
        dlogflag=0;

        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
        maxcovOvr = 1;

        imaxovr=0;
        iminovr=0;
                
                
        switch gcm_case


            case 'CloudSat dbZ from CFAD'
                xlabelstr='Longitude';
                ylabelstr='Height (km)';

                logflag=0; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                dlogflag=0;
                
                imaxovr=1;
                iminovr=1;

                mincovOvr = -60; %for normalised plots, one is the max (although unlikely to reach it)
                maxcovOvr = 5;

        tit(1).tit=[short_plot_name];
                
                
            case 'CALIPSO Cloud Fraction from CFAD'
                xlabelstr='Longitude';
                ylabelstr='Height (km)';

                logflag=0; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                dlogflag=0;
                
                imaxovr=1;
                iminovr=1;

                mincovOvr = 0; %for normalised plots, one is the max (although unlikely to reach it)
                maxcovOvr = 0.5;
                
                srdbz_str = eval(['num2str(sr_thresh_' gcm_str ')']);
                
                short_plot_name =[short_plot_name 'for SR thresh of ' srdbz_str]; 
                tit(1).tit = short_plot_name;

                

        end
                

        

     
        
        

    case 'MODIS_plot_UW'
        
        MODIS_plot='2D plot';
%        MODIS_plot='global flat map';    



        
        iadd_terrain=0;
        iplot_latlon=0;
        iadd_wind_quivers=0;
        iadd_flight_path=0;
        idraw_streamlines=0;
        clab=0;
        manclab=0;
        clines=0;
        icont=0; %flag to say whether want a contour or a colour plot
        pcolor_shading='flat';
        %pcolor_shading='interp';
        plot_type='pcolor';
%        plot_type='filled contour';
        

%%%% NOTE - these may be overrided in the individual plot cases below %%%%%
if ~exist('man_choose_plotTimeHeight_graph')
        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
        dlogflag=0;

        
        mincovOvr = 1e-2;  %usual one for CAS plots
        maxcovOvr = 10e3;

        mincovOvr = 1e-4;  
        maxcovOvr = 1;

        imaxovr=1;
        iminovr=1;
        
end

%        mincovOvr = 0;
%        maxcovOvr = 1000; %

        iovride_conts=0;
        conts_ovr=[270:2:310];
        conts_ovr=[270:1:300];
        %        conts_ovr=[0:0.5:19];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
        y_axis_type='log10';
        y_axis_type='log10_matlab'; %matlab's built in log scale - probably the best one to choose
        %as it does the minor tick marks between the major log numbers and puts ticks at appropriate points
        %also sets iytick_relabel to one as Matlab labels the ticks as 10^0 and 10^1 (rather than 1 and 10).
        y_axis_type='';
        %        y_axis_type='Psep';      

        iylim=0; %note is set below in certain cases
        ylims=[0 25];       

        izovr=2;
        i2d=3;

        idirstamp=0;
        itimestamp=0;
        sig=3;
        
        
        switch MODIS_plot
            case 'global flat map'
                short_plot_name=[MODIS_varname];
                short_plot_name=remove_character(short_plot_name,'_',' ');
                xlabelstr='Longitude';
                ylabelstr='Latitude';

            case '2D plot'
                modis2d_case = 'CF35 vs CF06'; 
                modis2d_case = 'Nd with W calculated vs with MODIS W';
                modis2d_case = 'Nd from histogram vs using mean tau&reff';
                modis2d_case = 'W as funciton of tau&reff';
                modis2d_case = 'Nd vs scattering angle';
                modis2d_case = 'Nd vs cloud fraction';
                modis2d_case = 'Nd vs total no. pixels';
%                modis2d_case = 'Nd vs tau';
                modis2d_case = 'std(Nd) vs.';
%                modis2d_case = 'Total No. pixels vs tau';
%                modis2d_case = 'No. pixels vs CF';
                 modis2d_case = 'vs. LAT';
                 modis2d_case = 'general XY 2D PDF';
                 
%                 modis2d_case = 'Daily max SZA for one location vs time';
%                 modis2d_case = 'Zonal max SZA vs time';
%                 modis2d_case = 'Zonal min SZA vs time';    
%                 modis2d_case = 'Zonal mean of daily max SZA vs time';
%                 modis2d_case = 'Zonal std dev of daily max SZA vs time';
%                 modis2d_case = 'Zonal mean of daily max SZA minus one std dev vs time';
%                 modis2d_case = 'Zonal mean of daily min SZA minus one std
%                 dev vs time';
%                 modis2d_case = '10^{th} percentile of daily max SZA vs time';
%                 modis2d_case = '10^{th} percentile of daily min SZA vs time';                 
%                 modis2d_case = '90^{th} percentile of daily max SZA vs time';
%                 modis2d_case = '90^{th} percentile of daily min SZA vs time';                       
%                 modis2d_case = 'Percent of longitudes with minSZA lte 65 and maxSZA gt 65 vs time';    
%                 modis2d_case = 'Percentage Nd difference of allSZA vs lowSZA vs time and lat';    
%                 modis2d_case = 'Mean percentage Nd difference of allSZA vs lowSZA vs time and lat';                     
%                 modis2d_case = 'Percentage Re difference of allSZA vs lowSZA vs time and lat';
%                 modis2d_case = 'Mean percentage Re difference of allSZA vs lowSZA vs time and lat';    
%                 modis2d_case = 'Mean percentage Tau difference of allSZA vs lowSZA vs time and lat';    
                 
                 tit(1).tit = modis2d_case;
                 
                 short_plot_name='';
  
                
                 
                 
                
                
                %%% set some default values, so don't need to bother below
                %%% every time
                if ~exist('man_choose_plotTimeHeight_graph')
                    logflag=0; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                    dlogflag=0;
                    
                    mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                    maxcovOvr = 1;

                    imaxovr=0;
                    iminovr=0;
                end
                

                
                
                switch modis2d_case
                    case 'std(Nd) vs.'
                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;
                        
                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                    case 'Nd vs tau'
                        ylabelstr='Optical Depth';
                        xlabelstr='N_d from means (cm^{-3})';                                                

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;
                        
                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                        
                    case 'Nd vs no. pixels'
                        ylabelstr='No. pixels';
                        xlabelstr='N_d from means (cm^{-3})';                                                

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;
                        
                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                    case 'Nd vs cloud fraction'
                        ylabelstr='Cloud Fraction';
                        xlabelstr='N_d from means (cm^{-3})';                                                

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;
                        
                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                     case 'Nd vs scattering angle'
                        ylabelstr='Cloud Fraction';
                        xlabelstr='N_d from means (cm^{-3})';                                                

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;
                        
                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                        
                    case 'W as funciton of tau&reff'
                        
                        %                short_plot_name=remove_character(short_plot_name,'_',' ');
                        xlabelstr='Optical depth';
                        ylabelstr='Effective radius (\mum)';                                                

                        logflag=0; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;

                        mincovOvr = 1e-2;  %usual one for CAS plots
                        maxcovOvr = 10e3;

                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 20;

                        imaxovr=0;
                        iminovr=0;
                        
                    case 'CF35 vs CF06'


                        %                short_plot_name=remove_character(short_plot_name,'_',' ');
                        xlabelstr='Cloud Fraction (MOD35)';
                        ylabelstr='Cloud Fraction (MOD06)';

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;

                        mincovOvr = 1e-2;  %usual one for CAS plots
                        maxcovOvr = 10e3;

                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;

                    case 'Nd with W calculated vs with MODIS W'
                        xlabelstr = 'N_d, W calculated';
                        ylabelstr = 'N_d, MODIS W';

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                        dlogflag=0;

                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                         iylim=1; %note is set below in certain cases
                         ylims=[0 400];
                         
                         ixlim=1; %note is set below in certain cases
                         xlims=[0 300];
                         
                     case 'Nd from histogram vs using mean tau&reff'
                        xlabelstr = 'N_d, histogram (cm^{-3})';
                        ylabelstr = 'N_d, means (cm^{-3})';

                        logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so that the min value is not zero
                        dlogflag=0;

                        mincovOvr = 1e-5; %for normalised plots, one is the max (although unlikely to reach it)
                        maxcovOvr = 1;

                        imaxovr=1;
                        iminovr=1;
                        
                         iylim=0; %note is set below in certain cases
                         ylims=[0 800];
                         
                         ixlim=0; %note is set below in certain cases
                         xlims=[0 800];    
                         
                         


                end
                
                
        end
        
        
        
        
        
        
    case 'Particle size dist vs time'
                        

if ~exist('iset_plot_and_instrument')
    var_plot='dN/dlogD (cm^{-3} \mum^{-1})';
%    var_plot='LWC size dist (g cm^{-3})';
    %var_plot='dN/dD (cm^{-3} \mum^{-1})';
%    var_plot='N (cm^{-3})';
    %var_plot='log10 of particle separation';
    %var_plot='Backscatter dN/dlogD (cm^{-3} \mum^{-1})';
    %var_plot='dM/dlogD (\mug m^{-3} \mum^{-1})';
    %var_plot='M (\mug m^{-3})';

    instrument='CAS';
%    instrument='CAS no CIP'; %BAS CAS, but with no CIP (so using constant airspeed 
%          is probably the only difference - mainly for the labstudies)
    %instrument='CIP';
    %instrument='CDP';
    instrument='CAS MAN';
%    instrument='Welas';
%    instrument='FSSP';

else
    clear iset_plot_and_instrument
end

idatetick=0; %flag to say the want the xaxis in proper time format rather than decimal time
datetick_type=13; %specify the type with datetick_type (see help datetick) 15= HH:MM:SS
%datetick_type=15; %specify the type with datetick_type (see help datetick) 15= HH:MM

%NOTE - when using datetick (possibly when not too) some of the xticks are missing unless the
%ylims are set so that the colour plot doesn't go all the way to the edge (i.e. leave a gap
%at the bottom and top of the plot by setting appropriate ylim values).






iadd_terrain=0;
iplot_latlon=0;
iadd_wind_quivers=0;
iadd_flight_path=0;
idraw_streamlines=0;
clab=0;
manclab=0;
clines=0;
icont=0; %flag to say whether want a contour or a colour plot
pcolor_shading='flat';
%pcolor_shading='interp';
plot_type='pcolor';
%plot_type='filled contour';
           
%air_speed_type = 'aircraft'; %from aircraft 
%air_speed_type = 'constant 60m/s'; %assume constant 60 m/s
%air_speed_type = 'CIP probe'; %airspeed measured by CIP/CAS - corrected for airflow
%through the instrument tube using static and dynamic pressures

%%% choose contour range %%%
if ~exist('iset_colour_limits') %if flag is undefined
    
imaxovr=1;
iminovr=1;

mincovOvr = 0;
maxcovOvr = 0.1; %N for volcano

mincovOvr = 0; %dN/dlogD volcano
maxcovOvr = 5;

mincovOvr = 0; %M (microgrammes per m3) volcano
maxcovOvr = 9000;

mincovOvr = 0; %dM/logD (microgrammes per m3) volcano
maxcovOvr = 6e5;

%mincovOvr = 270;
%maxcovOvr = 307;

%mincovOvr = -2;
%maxcovOvr = 20;



logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
dlogflag=0;


mincovOvr = 1e-2;  %usual one for CAS plots
maxcovOvr = 10e3;


% mincovOvr = 1e-2;
% maxcovOvr = 200;
% 
% mincovOvr = 1e-2;
% maxcovOvr = 10000;
% 
% mincovOvr = 1e-2;
% maxcovOvr = 1000;

%mincovOvr = 0;
%maxcovOvr = 0.05;

% mincovOvr = 0;
% maxcovOvr = 0.9;
% 
% mincovOvr = 0;
% maxcovOvr = 0.32;
% 
% mincovOvr = 0;
% maxcovOvr = 0.8;

%mincovOvr = 0;
%maxcovOvr = 4000;

%mincovOvr = 0;
%maxcovOvr = 370;

else
    clear iset_colour_limits %clear the flag for subsequent runs
end

iovride_conts=0;
conts_ovr=[270:2:310];
conts_ovr=[270:1:300];        
%        conts_ovr=[0:0.5:19];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
        y_axis_type='log10';
        y_axis_type='log10_matlab'; %matlab's built in log scale - probably the best one to choose
        %as it does the minor tick marks between the major log numbers and puts ticks at appropriate points
        %also sets iytick_relabel to one as Matlab labels the ticks as 10^0 and 10^1 (rather than 1 and 10).
        y_axis_type='';
%        y_axis_type='Psep';
        

%automatically sets the right flight data here for the airspeed calc based on the CAS
%fileName
%    dat_flt_str = ['dat_flt' flight_no]; %make a string that is the name of the required flight data

%put the data in dat_flt array
%if no_flight_data==0
%    eval(['dat_flt=' dat_flt_str ';']);
%end
    
    
    
    
    
        iylim=0; %note is set below in certain cases
        ylims=[0 25];
        
        
if logflag==1 %if doing log10 of the colour scale then best to do pcolor plot
    plot_type=''; %as otherwise it tends to crash
end
%plot_type='filled contour'; %
            


        

switch instrument
    case {'CAS','CAS no CIP'}
        if idatetick==1
            time_plot2D = (CAS_time_all)/3600/24;
        else
            time_plot2D = (CAS_time_all)/3600;
        end
        
        if strcmp(instrument,'CAS no CIP')==1
            if idatetick==1
                time_plot2D = time_plot2D - 1/24; %to account for Tom LC's computer working on BST
            else
                time_plot2D = time_plot2D - 1;
            end
        end
        
    
        
        switch y_axis_type
            case 'log10'
                iytick_override=1; %set this to re-label the tick marks as the actual diameters
                tick_locs = log10([0.5 0.6 0.7 0.8 0.9 1:10 20 30 40]);
                iylim=0;
                ylims=[tick_locs(1) tick_locs(end)];
                iytick_relabel=1; %re-label the y axis so that shows the actual diameter rather than logD
            case {'log10_matlab'}
                iytick_relabel=1; %re-label the y axis so that shows the actual diameter rather than logD
                iylim=1;
                ylims=[0.5 50];                
        end

        %change the x limits to move the axis slightly to the right in order to see the tick marks.
        ixlim=1;
        xlims=[time_plot2D(1) - (time_plot2D(end) - time_plot2D(1) )/100 time_plot2D(end)];
        
    case 'CIP'
        if idatetick==1
            time_plot2D = (CIP_time_all)/3600/24;
        else
            time_plot2D = (CIP_time_all)/3600;
        end
        
        iylim=1;
        ylims=[0 1550];

        iytick_relabel=0; %re-label the y axis so that shows the actual diameter rather than logD
        y_axis_type='';
        
        ixlim=1;
        xlims=[time_plot2D(1) - (time_plot2D(end) - time_plot2D(1) )/100 time_plot2D(end)];
        
   case 'CDP'
        if idatetick==1
            time_plot2D = CDP_time_all'/3600/24;
        else
            time_plot2D = CDP_time_all'/3600;
        end
        
        iylim=1;
        ylims=[0 52];

        iytick_relabel=0; %re-label the y axis so that shows the actual diameter rather than logD
        y_axis_type='';
        
        ixlim=1;
        xlims=[time_plot2D(1) - (time_plot2D(end) - time_plot2D(1) )/100 time_plot2D(end)]; 
        
        
   case 'CAS MAN' %the Manchester CAS
        if idatetick==1
            time_plot2D = data_CAS_PACS(1,:)'/3600/24;
        else
            time_plot2D = data_CAS_PACS(1,:)'/3600;
        end
        
        switch y_axis_type
            case 'log10'
                iytick_override=1; %set this to re-label the tick marks as the actual diameters
                tick_locs = log10([0.5 0.6 0.7 0.8 0.9 1:10 20 30 40]);
                iylim=0;
                ylims=[tick_locs(1) tick_locs(end)];
                iytick_relabel=1; %re-label the y axis so that shows the actual diameter rather than logD
            case {'log10_matlab'}
                iytick_relabel=1; %re-label the y axis so that shows the actual diameter rather than logD
                iylim=1;
                ylims=[0.5 50];                
        end

        %change the x limits to move the axis slightly to the right in order to see the tick marks.
        ixlim=1;
        xlims=[time_plot2D(1) - (time_plot2D(end) - time_plot2D(1) )/100 time_plot2D(end)];  
        
        
     case 'Welas'
        if idatetick==1
            time_plot2D = welas_dat.time_of_day/3600/24;
        else
            time_plot2D = welas_dat.time_of_day/3600;
        end
        
        iylim=1;
        ylims=[min(welas_dat.size(:,1)) max(welas_dat.size(:,1))];

        iytick_relabel=0; %re-label the y axis so that shows the actual diameter rather than logD
        y_axis_type='';
        
        ixlim=1;
        xlims=[time_plot2D(1) - (time_plot2D(end) - time_plot2D(1) )/100 time_plot2D(end)]; 
    
     case 'FSSP'
        if idatetick==1
            time_plot2D = time_FSSP/24;
        else
            time_plot2D = time_FSSP;
        end
        
        iylim=1;
        ylims=[min(bins_FSSP) max(bins_FSSP)];

        iytick_relabel=1; %re-label the y axis so that shows the actual diameter rather than logD
%        y_axis_type='';
        
        ixlim=1;
        xlims=[time_plot2D(1) - (time_plot2D(end) - time_plot2D(1) )/100 time_plot2D(end)]; 
        
end

%        xlims=[19.3 19.55]/24;
%        xlims=[18.8 19.2]/24;

switch y_axis_type
    case 'log10'
%        ylabelstr='log_{10} D (microns)';
        ylabelstr='Diameter (\mum)';
    case '';        
        ylabelstr='Diameter (\mum)';
    case 'Psep'
        ylabelstr='?';
    otherwise
        ylabelstr='Diameter (\mum)';
end

   xlabelstr='Time (UTC)';
   
   
                
            
%        end
        


switch instrument
    case {'CAS','CIP'};
        short_plot_name=[var_plot ' for flight ' flight_no ' on ' date_str ' at ' CAS_start_time_str ' (' instrument ')'];         
    otherwise
        short_plot_name=[var_plot ' (' instrument ')'];   % 'for ' date_str ' at ' CAS_start_time_str ' (' instrument ')'];tit(1).tit=[short_plot_name];         
end

    

tit(1).tit=[short_plot_name];

        
        


        
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
     
        
        sig=3;
        
%%%%%%%% **********************************************************************        

 case 'wrf_vert_cross_section'


hor_vert=0;   %choice of how to choose cross section - normal cross section where choose the start and end points
%hor_vert=4; %along mountain crest (max mountain height location vs latitude)
%hor_vert=3; %to do a cross section along a streamline

switch hor_vert
    case {0,1,2,3}             
        x_axis_type='dist';
    case 4
        x_axis_type='lat';
end

x_axis_type='lon';
%x_axis_type='lat'; %slice along a constant 
%x_axis_type='dist';


ilon_line=0; %flag to say that want a slice along a constant longitude (set hor_vert=0)
if ilon_line==1
        x_axis_type='lat'; %slice along a constant 
        LON=-66; %the longitude required
end







iadd_streamline_z0=0; %add the streamline used for a streamline cross section
iadd_effective_surface=0;

%% variable for vertical cross section
 if ~exist('man_choose_plotTimeHeight_graph') | man_choose_plotTimeHeight_graph == 0
     
var_plot='Equivalent potential temperature (K)';
var_plot='Specific humidity (g kg^{-1})';
var_plot='Relative humidity (%)';
var_plot='Potential temperature (K)';
%var_plot='Temperature (^{o}C)';
%var_plot='Wind speed (m s^{-1})';
%var_plot='Wind speed perpendicular component (m s^{-1})';
%var_plot='Vertical wind speed (m s^{-1})';
%var_plot='Component (UVW) wind speed (m s^{-1})'; %Horiz wind speed ( sqrt(u^2+v^2+w^2) in orientation of the slice)
%var_plot='Component horizontal (UV) wind speed (m s^{-1})'; %Horiz wind component (u & v only) in orientation of the slice
%var_plot='Wind direction (degrees)';
%var_plot='Density (kg m^{-3})';
%var_plot='Richardson Number (perpendicular wind)'; %dimensionles 
%var_plot='';

 end
 
% *** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recalc=1; %whether to recalculate the fields or just replot
% *** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recalc==0
    disp('----   Information :- recalc is set to zero  ------');
end

ncross_points=50;
ncross_points=150;
ncross_points=300;

y_cross_sec = 380; %position in y axis for cross section (km)
%y_cross_sec = 280; %position in y axis for cross section (km)
%y_cross_sec = 425; %position in y axis for cross section (km)     
%y_cross_sec = 295; %position in y axis for cross section (km) - location of jet (EC) at 12 UTC for NCEP
%y_cross_sec = 301; %position in y axis for cross section (km) - location of jet (EC) at 12 UTC for ECMWF
     

iadd_terrain=0;
iplot_latlon=0;
iadd_wind_quivers=0;
iadd_flight_path=0;
idraw_streamlines=0;
clab=0;
manclab=0;
clines=1;
icont=1;

%%% choose contour range %%%
iovride_conts=0;
conts_ovr=[270:2:310];
conts_ovr=[270:1:300];        

imaxovr=1;
iminovr=1;

switch var_plot
    case {'Component horizontal (UV) wind speed (m s^{-1})','Component (UVW) wind speed (m s^{-1})','Wind speed (m s^{-1})','Wind speed perpendicular component (m s^{-1})'}
        mincovOvr = -2;
        maxcovOvr = 20;

        mincovOvr = -2.5;
        maxcovOvr = 20.5;
        
        mincovOvr = -10;
        maxcovOvr = 5;

        iovride_conts=0;
        conts_ovr=[-2:1:20];
        cont_val=9.5;
%        conts_ovr=[cont_val cont_val]; %plotting a single contour

    case {'Relative humidity (%)'}
        mincovOvr = 40;
        maxcovOvr = 100;


        iovride_conts=0;
        conts_ovr=[40:5:100];
        cont_val=9.5;
%        conts_ovr=[cont_val cont_val]; %plotting a single contour


    case 'Potential temperature (K)'
        mincovOvr = 270;
        maxcovOvr = 300;
    case 'Richardson Number (perpendicular wind)'
         mincovOvr = 0;
        maxcovOvr = 30;
    case 'Vertical wind speed (m s^{-1})'
        mincovOvr = -3;
        maxcovOvr = 1;
end

%mincovOvr = 270;
%maxcovOvr = 307;



%mincovOvr = 0.8;
%maxcovOvr = 1.4;

%        conts_ovr=[0:0.5:19];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set time here for cross section
switch file_type
    case 'wrfout'
        if ~exist('man_choose_plotTimeHeight_graph') | man_choose_plotTimeHeight_graph == 0
        time=11; %6th Jan:- 11=6UTC, 13=12UTC, 15=18UTC, 17=00UTC
        time=12; %
%        time=idir;
        end
    case {'met_em','wrfinput'}
        time=1;
end

        
                                  
%        tit(1).tit=['Latent Heat Flux (W m^{-2})'];^M
%        tit(1).tit=['Sensible Heat Flux (W m^{-2})'];
%	tit(1).tit=['Radar reflectivity'];

%        if ilon_slice==1
%            tit(1).tit=['Radar reflectivity vertical slice at  ' num2str(lon_slice) '^{o} lon (dBZ)'];
%            tit(1).tit=['Total water vertical slice at  ' num2str(lon_slice) '^{o} lon (ppmv)'];




%	xlabelstr='Distance (km)';
%        else
%            tit(1).tit=['Radar reflectivity vertical slice at  ' num2str(lat_slice) '^{o} lat (dBZ)'];

switch x_axis_type
    case 'dist'
        xlabelstr='Distance (km)';
    case 'lat';
        xlabelstr='Latitude';
    case 'lon';
        xlabelstr='Longitude';    
end
                
            
%        end
        
        iylim=0;
        iylims=[0 25];
        
        ixlim=0;
        xlims=[-69 -67.15];

%         tit(1).tit=['Wind cross section at y=' num2str(y_cross_sec) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];              
%         tit(1).tit=['Potential temperature (K) cross section at y=' num2str(y_cross_sec) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
%         short_plot_name=['Equivalent potential temperature (K) cross section at y=' num2str(y_cross_sec) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];tit(1).tit=[short_plot_name ' for ' filestr];
%         short_plot_name=['Equivalent potential temperature (K) cross section for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];tit(1).tit=[short_plot_name ' for ' filestr];
short_plot_name=[var_plot ' cross section for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];tit(1).tit=[short_plot_name ' for ' filestr];         
%         short_plot_name=['Temperature (^{o}C) cross section at y=' num2str(y_cross_sec) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];tit(1).tit=[short_plot_name ' for ' filestr];
%         short_plot_name=['Wind speed (m s^{-1}) cross section at y=' num2str(y_cross_sec) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];tit(1).tit=[short_plot_name ' for ' filestr];    
%         short_plot_name=['Wind speed (m s^{-1}) cross section for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];tit(1).tit=[short_plot_name ' for ' filestr];    
        %   tit(1).tit=['ECMWF mass weighted mean temperature over lower 2.5 km at ' tlab ' UTC, 25th Feb (^{o}C)'];
        
        

        
        


        
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
        ylabelstr='Height (km)';
        
        sig=3;
        
        
	
      

   case 'wrf_radar_vert'
                
        tit(1).tit=['Latent Heat Flux (W m^{-2})'];
%        tit(1).tit=['Sensible Heat Flux (W m^{-2})'];
        tit(1).tit=['Radar reflectivity'];              
        tit(1).tit=['Total water'];              
         tit(1).tit=['Vertical velocity'];              
        if ilon_slice==1
            tit(1).tit=['Radar reflectivity vertical slice at  ' num2str(lon_slice) '^{o} lon (dBZ)'];
            tit(1).tit=['Total water vertical slice at  ' num2str(lon_slice) '^{o} lon (ppmv)'];
            xlabelstr='Distance (km)';

        else            
            tit(1).tit=['Radar reflectivity vertical slice at  ' num2str(lat_slice) '^{o} lat (dBZ)'];
            xlabelstr='Distance (km)';
        end
        
        iylim=1;
        iylims=[0 25];
        ixlim=0;
        xlims=[101 201];                
    
        %   tit(1).tit=['ECMWF mass weighted mean temperature over lower 2.5 km at ' tlab ' UTC, 25th Feb (^{o}C)'];
        
        imaxovr=1;
        iminovr=1;
        
        mincovOvr = 278;
        maxcovOvr = 298;
        clab=0;
        clines=0;
        icont=1;
        
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
        ylabelstr='Height (km)';
        
        sig=3;
        
        clab=0;
        
        iovride_conts=1;
        conts_ovr=[0:5:65];
		conts_ovr=[6:0.1:7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODIS plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'MODIS_plot'
        
        plot_type='pcolor';
%        plot_type='filled contour';
%        plot_type='3D surf'; %set iylim=0;        ixlim=0;
        
        CameraPosition= [2000 -690.8 14396];
        CameraTarget = [619.0630 219.0630 1.6905e+003];
        CameraViewAngle = 14.3;
        CameraViewAngle = 10.3;
        
        switch plot_type
            case '3D surf'
                colorbar_location = 'WestOutside';
                colorbar_location = 'EastOutside';
        end
        
        icolmap=1;
        lb_map=lbmap(256,'brownblue'); %nice colormap for colorblind people
        
        i_reverse_cmap=1;
        if i_reverse_cmap==1
            cmap=flipud(lb_map);
        else
            cmap=lb_map;
        end
        
%        cmap=jet;        
%        cmap=hsv;

        


        iylim=1;
        ixlim=1;
        
        xlims=[200 500];
        ylims=[250 550];                     

        
        
        iadd_terrain = 1;
        iadd_wind_quivers = 0;
        iplot_latlon=1;
        iadd_flight_path=1;
        
        isquare=1;
        
        iovride_conts=0;
        conts_ovr=[-70:10:40];
        
        clab=0;
        manclab=0;
        clines=0;
        icont=1;
        
        
        imaxovr=0;
        iminovr=0;
        
         mincovOvr = 0;
         maxcovOvr = 18;
                  
%         mincovOvr = 0;
%         maxcovOvr = 12;
        
        
    sig=4;       
	ncont=15; 
%	ncont=28;

            
      iset_locations_plotTimeHeight=0;
      
      if iset_locations_plotTimeHeight==1
          location_lab(1).l = 'A';
          location_lab(2).l = 'B';  %'M'
          location_lab(3).l = 'E';
          location_lab(4).l = 'AWS';
          location_lab(5).l = 'F';

          %times_flight_loc = [20.3755 21.9815 22.0205]; %20.3755 is the time where the max was seen for aircraft data (index=269899)
          %21.9815 is the time of second wind peak on the mini descent before the last ascent
          times_flight_loc = [20.3755 22.0205];

          for iflt=1:length(times_flight_loc)
              it_flt(iflt) = findheight(time_flt19,times_flight_loc(iflt));
          end
          LAT_extra = dat_flt19(it_flt,2)';
          LON_extra = dat_flt19(it_flt,3)';

          LAT_extra2 =[];  %the locations specified by extra_x and extra_y will be stored in LAT/LON_extra2
          LON_extra2 =[];

%          extra_x=[];
%          extra_y=[];

          %now calculate the lat/lon of the extra_x and extra_y points and put them in LAT_extra2
          for iflt=1:length(extra_x)
              i_extra = findheight_nearest(i_grid,extra_x(iflt));
              j_extra = findheight_nearest(j_grid,extra_y(iflt));
              LAT_extra2(iflt) = lat2d.var(j_extra,i_extra);
              LON_extra2(iflt) = lon2d.var(j_extra,i_extra);
          end

          LAT=[-67.01]; %Larsen AWS
          LON=[-61.55];

          iplot_aircraft_locs_model=1; %add LAT_extra points
          if iplot_aircraft_locs_model==1
              LAT = [LAT_extra LAT];
              LON = [LON_extra LON];
          end

          %now add the requested extra_x, extra_y locations for model plot
          LAT=[LAT LAT_extra2];
          LON=[LON LON_extra2];

          [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);





      end


     
         short_plot_name=['Cloud Top Pressure'];
%         short_plot_name=['First level height ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m)'];
         %         short_plot_name=['Wind speed (m s^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filestr='MOD06_L2.A2010039.1405.005.2010039225429.hdf';
tit(1).tit=[short_plot_name ' for ' filestr];
short_plot_name = tit(1).tit; %might want to comment this out if filename is too long
                                                  
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
%        ylabelstr='iLatitude';
%        xlabelstr='iLongitude';
        
        ylabelstr='y (km)';
        xlabelstr='x (km)';
        
        
%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 'wrf_plot'   %Choose xlims, clims, ih_wrf, etc below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear tit
        
        if ~exist('man_choose_plotTimeHeight_graph') | man_choose_plotTimeHeight_graph==0
            var_plot = 'Wind speed Nth model level';
%            var_plot = 'Temperature Nth model level';
%            var_plot = 'set below';  %this is set below for most vars, except 'Wind speed Nth model level' - need to add names for the others.
        end
        
        plot_type='pcolor';
        plot_type='filled contour';
%        plot_type='3D surf'; %set iylim=0;        ixlim=0;
        
        CameraPosition= [2000 -690.8 14396];
        CameraTarget = [619.0630 219.0630 1.6905e+003];
        CameraViewAngle = 14.3;
        CameraViewAngle = 10.3;
        
        switch plot_type
            case '3D surf'
                colorbar_location = 'WestOutside';
                colorbar_location = 'EastOutside';
        end
        
        icolmap=1;
        lb_map=lbmap(256,'brownblue'); %nice colormap for colorblind people
        
        i_reverse_cmap=1;
        if i_reverse_cmap==1
            cmap=flipud(lb_map);
        else
            cmap=lb_map;
        end
        
%        cmap=jet;        
%        cmap=hsv;

        colorbar_loc='vert'; %default
%        colorbar_loc='horiz';

 if ~exist('man_choose_plotTimeHeight_graph') | man_choose_plotTimeHeight_graph==0        
        ih_wrf=82; %z-level to plot for
        ih_wrf=0; %z-level to plot for    --- set to zero for the 10 m winds ---    
%        ih_wrf=10; %z-level to plot for
        ih_wrf=4; %z-level to plot for        
%        ih_wrf=idir; %z-level to plot for

 end

%% Set time here for the 'wrf_plot' plots
switch file_type
    case 'wrfout'
        if ~exist('man_choose_plotTimeHeight_graph') | man_choose_plotTimeHeight_graph==0
            time = 1; %11=06, 12=09, 13=12, 14=15, 15=18, 16=21 17=0, 18=03, 19=06 UTC
            time = 12; %11=06 UTC for NCEP run, 45=20 for ECMWF
            %      time = idir;

        end
     
    case {'met_em','wrfinput'}
     time=1;
    case 'ecmwf'
end

        iylim=1;
        ixlim=1;
        
        ylims=[900 1800];
        xlims=[1100 2000];     %for d02 (11th Aug, 2014).            
        
        ylims=[3200 4100];  %
        xlims=[3400 4300];     %for d01 (11th Aug, 2014). 
        
%        ylims=[125 525];
%        xlims=[300 700]; %used before 1st June, 2010 (d03)
                         % Still same for 11th Aug, 2014

%        ylims=[175 425];
%        xlims=[300 550];
                
%        ylims=[175 575];
%        xlims=[300 700];
               
       

%        ylims=[75 575];
%        xlims=[250 750];
        
%        ylims=[200 400];
%        xlims=[350 550];  
        
%        ylims=[390 460];
%        xlims=[330 400];  

%        ylims=[120 470];
%        xlims=[280 630];



%        xlims=[1200 2500];
%        ylims=[900 2200]; 
        
        h_wrf=5.57; %height (km) for which to interpolate for
        h_wrf=9.03; %height (km) for which to interpolate for
        h_wrf=2.3; %height (km) for which to interpolate for
%        h_wrf=7.5; %height (km) for which to interpolate for
%        h_wrf=1.3; %height (km) for which to interpolate for

        p_wrf=300; %pressure (mb) to interpolate for if i_p_int==1
        potemp_wrf = 310;   %6 UTC
%        potemp_wrf = 288.5; %6 UTC 
%        potemp_wrf = 282.1; %0 UTC 
%        potemp_wrf = 288; %0 UTC 
        
        iadd_terrain = 1;
        iadd_wind_quivers = 1;
        iplot_latlon=1;
        iadd_flight_path=1;
        
        isquare=1;
        

        
        clab=0;
        manclab=0;
        clines=0;
        icont=1;
        
%%%%%%%%  set the colour limits here  %%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imaxovr=1;  %on/off flags - set the limits below
        iminovr=1;
        
        iovride_conts=0;  %this overrides mincovOvr and maxcovOvr
        conts_ovr=[991:0.5:998];
        conts_ovr=[990.5:0.1:992.5];        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        mincovOvr = 970;
        maxcovOvr = 1000;
        
%        mincovOvr = 989;
%        maxcovOvr = 1000;
 
        
%        mincovOvr = 990.2; %for the Coriolis plots
%        maxcovOvr = 993;

%        mincovOvr = 990.5; %for the Coriolis plots, 15 UTC
%        maxcovOvr = 992.5;

%        mincovOvr = 954; %for level 4 pressure ecmwf_nudging case
%        maxcovOvr = 956.4;
        
         mincovOvr = -5;
         maxcovOvr = 10;
         
%         mincovOvr = 0;
%         maxcovOvr = 150;

%         
%         
%         mincovOvr = 725;
%         maxcovOvr = 750;
%         
% 
%         mincovOvr = 0;
%         maxcovOvr = 18;
         
%          mincovOvr = 0;
%         maxcovOvr = 15;
         
         %mincovOvr = 200;
         %maxcovOvr = 250;
                  

%         mincovOvr = 60;
%         maxcovOvr = 100;
         
%         mincovOvr = 990.5;
%         maxcovOvr = 992.5;
         
%         mincovOvr = 990.0;
%         maxcovOvr = 991.6;
         
%          mincovOvr = 954;
%         maxcovOvr = 956.4;
        
        
nx_quiver=25; %number of arrows to draw for x
ny_quiver=25; %number for y
scale_speed_quiver = [15 15]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
%scale_speed_quiver = [7 7]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
%scale_speed_quiver = [10 10]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed

nx_quiver=50; %number of arrows to draw for x
ny_quiver=50; %number for y

nx_quiver=25; %number of arrows to draw for x
ny_quiver=25; %number for y


    sig=4;       
	ncont=15; 
%	ncont=28;


      inversion_height1=1750;
      inversion_height2=2750;
      
%      inversion_height1=750;
%      inversion_height2=1500;

      inversion_height1=250;
      inversion_height2=1500; %lower inversion for max gradient calcualtion

      
      iadd_AWS=0; %add the AWS point to the end of the other locations as output from WaterVapour...
      if iadd_AWS==1
          LAT_orig=LAT;
          location_lab(length(LAT)+1).l = 'AWS'; 
         
          LAT=[LAT -67.01]; %Larsen AWS
          LON=[LON -61.55];
          
          [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);          
      end
      
      iset_locations_plotTimeHeight=0;
      
      if iset_locations_plotTimeHeight==1
          location_lab(1).l = 'A';
          location_lab(2).l = 'B';  %'M'
          location_lab(3).l = 'AWS';
          location_lab(4).l = 'E';
          location_lab(5).l = 'F';

          %times_flight_loc = [20.3755 21.9815 22.0205]; %20.3755 is the time where the max was seen for aircraft data (index=269899)
          %21.9815 is the time of second wind peak on the mini descent before the last ascent
          times_flight_loc = [20.3755 22.0205];

          for iflt=1:length(times_flight_loc)
              it_flt(iflt) = findheight_nearest(time_flt19,times_flight_loc(iflt));
          end
          LAT_extra = dat_flt19(it_flt,2)';
          LON_extra = dat_flt19(it_flt,3)';

          LAT_extra2 =[];  %the locations specified by extra_x and extra_y will be stored in LAT/LON_extra2
          LON_extra2 =[];

          extra_x=[];
          extra_y=[];

          %now calculate the lat/lon of the extra_x and extra_y points and put them in LAT_extra2
          for iflt=1:length(extra_x)
              i_extra = findheight_nearest(i_grid,extra_x(iflt));
              j_extra = findheight_nearest(j_grid,extra_y(iflt));
              LAT_extra2(iflt) = lat2d.var(j_extra,i_extra);
              LON_extra2(iflt) = lon2d.var(j_extra,i_extra);
          end

          LAT=[-67.01]; %Larsen AWS
          LON=[-61.55];
            
          %these are the descent and ascent locations (A and B)
          iplot_aircraft_locs_model=1; %add LAT_extra points
          if iplot_aircraft_locs_model==1
%              LAT = [LAT_extra LAT];
%              LON = [LON_extra LON];
               LAT = [-67.1984  -67.1557 LAT];
               LON = [-61.7080  -61.9106 LON];
          end

          %now add the requested extra_x, extra_y locations for model plot
          LAT=[LAT LAT_extra2];
          LON=[LON LON_extra2];

          [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);





      end





        
    
%	time = idir; %do this if want mulitple plots from multisave

        
       

% switch file_type
%     case {'wrfout','wrfinput'}
%         Zlev=WRFUserARW(nc,'Z',time,ih_wrf);
%         terrain = nc{'HGT'}(time,:,:);
%         medZ=num2str( mean(mean((squeeze(Zlev) - terrain))) ,3 );
%     case 'met_em'
%         try
%         Zlev=nc{'GHT'}(time,ih_wrf,:,:);
%         terrain = nc{'HGT_M'}(time,:,:);
%         medZ=num2str( mean(mean((squeeze(Zlev) - terrain))) ,3 );
%         catch
%             medZ='';
%         end
%     case 'ecmwf'
%         medZ='';
% end

switch file_type
            case 'wrfout'
                Zlev=WRFUserARW(nc,'Z',time,ih_wrf);
                terrain = nc{'HGT'}(time,:,:);
                medZ=['(median of ' num2str( median(median((squeeze(Zlev) - terrain))) ,3 ) ' m above terrain)'];
            case 'met_em'
                if prod(size(nc{'GHT'})) > 0
                    Zlev=nc{'GHT'}(time,ih_wrf,:,:);
                    %terrain = nc{'HGT_M'}(time,:,:);
                    medZ=['(mode of ' num2str( mode(mode((squeeze(Zlev)))) ,3 ) ' m above terrain)'];
                else
                    Plev=nc{'PRES'}(time,ih_wrf,:,:)/100;
                    medZ=['(median of ' num2str(median(median(Plev)),'%.0f') ' hPa)'];
                end
            otherwise
                medZ='';
end
        
        
        


year=Times(time,1:4);
day=Times(time,9:10);
month=Times(time,6:7);
hour=Times(time,12:16);
month_str=datestr(datenum([year '-' month '-' day]),'mmm');

%% What to plot for wrf_plot type
switch var_plot
    case 'Wind speed Nth model level'
          short_plot_name=['Wind speed (m s^{-1}) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC at level ' num2str(ih_wrf) ' ' medZ];  
          
    case 'Temperature Nth model level'
        short_plot_name=['Temperature (^{o}C) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC at level ' num2str(ih_wrf) ' ' medZ];
        mincovOvr = -5;
        maxcovOvr = 10;
        
    otherwise  %set these others to work with this system sometime.
        
%ones that have been modified so that just need to select the choice here and not in the other part too
         short_plot_name=['Terrain height ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m)'];
%         short_plot_name=['First level height ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m)']; %         


%         short_plot_name=['Sea Ice flag for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];
%         short_plot_name=['10 m wind speed for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];
%         short_plot_name=['Surface pressure (hPa) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC']; 
%         short_plot_name=['Pressure(hPa) at model level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC'];
%         short_plot_name=['Total water at level ' num2str(ih_wrf) ' (ppmv)'];
%          short_plot_name=['Skin temperature for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C)'];%          short_plot_name=['Air density at level ' num2str(ih_wrf) ' (kg m^{-3})'];
%        short_plot_name=['2m air temperature ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C)'];
%        short_plot_name=['2m vapour mixing ratio ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C)'];        
        % short_plot_name=['Ice supersaturation (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];
        % short_plot_name=['Hallet Mossop flag at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];
%          short_plot_name=['RH (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         %       short_plot_name=['Total water mixing ratio (g kg^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
         
%        short_plot_name=['Latent Heat Flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2})']; tit(1).tit=short_plot_name;
%        short_plot_name=['Sensible Heat Flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2})'];
%        short_plot_name=['Sensible + Latent Heat Flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2})'];
soil_level=4; %level one is closest to the surface, 4 is deep below
%        short_plot_name=['Soil temperature (^{o}C) at level ' num2str(soil_level) ' at ' Times(time,9:10) ' Jan ' Times(time,12:16) ];
%        short_plot_name=['2m air temperature difference (ECMWF-NCEP)' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C)'];        

%        short_plot_name=['2m vapour mixing ratio ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (g kg^{-1})'];

%        short_plot_name=['10m wind speed ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m s^{-1})'];

         
 %        short_plot_name=['Melt energy for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2}) for ' filestr];
%         short_plot_name=['SW downward flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2}) for ' filestr];

 %        short_plot_name=['Average melt rate (mm day^{-1} w.e.) for ' filestr];
 %        short_plot_name=['Total melting (mm w.e.) for ' filestr];
 %        short_plot_name=['Total melt flux (W m^{-2}) for ' filestr]; %don't forget to choose the time for the 
         %wind vectors (if requried)  %0 150 for clims
 %        short_plot_name=['SW net flux (W m^{-2}) for ' filestr];
%         short_plot_name=['SW downwards flux (W m^{-2}) for ' filestr]; %0 650 for clims
 %        short_plot_name=['LW net flux (W m^{-2}) for ' filestr];
%        short_plot_name=['LW down flux (W m^{-2}) for ' filestr]; %200 250 for clims
%         short_plot_name=['SH flux (W m^{-2}) for ' filestr];  % -30 10 for clims
%         short_plot_name=['LH flux (W m^{-2}) for ' filestr];   %-10 to +10 for clims
%         short_plot_name=['GRD flux (W m^{-2}) for ' filestr];  %-30 to +5 for clims             
         
%         short_plot_name=['ALBEDO for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' for ' filestr];

%          short_plot_name=['Mean condensate (g m^{-2}) up to level 30 for ' filestr];
%         short_plot_name=['Mean SW melt contribution (mm day{-1}) for ' filestr];
%          short_plot_name=['Mean 2m temperature (^{o}C) up to level 30 for ' filestr];
%         short_plot_name=['Total SH melt contribution (mm) for ' filestr];
%         short_plot_name=['Relative humidity at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain), for ' filestr];
         
%         short_plot_name=['LW net flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2}) for ' filestr];
%         short_plot_name=['LW downwards flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (mm day^{-1}) for ' filestr];

%         short_plot_name=['Ground heat flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2}) for ' filestr];
%         short_plot_name=['Sensible + latent heat flux for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (W m^{-2}) for ' filestr];

%         short_plot_name=['Max cloud mixing ratio for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (g kg^[-1}) for ' filestr];
         
%        tit(1).tit=['Surface tests for ' Times(time,9:10) ' Jan '        Times(time,12:16) ' (^{o}C) for ' filestr];
%        tit(1).tit=['SNOW for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]; 
%        tit(1).tit=['Snow cover flag for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];         


        


 

%         short_plot_name=['Temperature (^{o}C) at ' num2str(h_wrf) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];  

%         short_plot_name=['Potential temperature (K) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
%         short_plot_name=['Equivalent potential temperature (K) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
         
%       short_plot_name=['Vapour mixing ratio (g kg^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];

%       tit(1).tit=['Vapour mixing ratio (g kg^{-1}) at ' num2str(h_wrf) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];

%       short_plot_name=['Condensed water mixing ratio (g kg^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%       short_plot_name=['Ice water mixing ratio (g kg^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];       
%       short_plot_name=['Snow water mixing ratio (g kg^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];              
%       short_plot_name=['Rain water mixing ratio (g kg^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];              

       
%       short_plot_name=['Snow and ice precipitation (mm) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%       short_plot_name=['Rain precipitation (mm) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%       short_plot_name=['Precipitation (mm) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%      short_plot_name=['Precipitation tendency (mm) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
       
%       short_plot_name=['Total number concentration(L^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];       
%       short_plot_name=['Rain number concentration(L^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];       

%      short_plot_name=['Ice number concentration(L^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];       
%      short_plot_name=['Snow number concentration(L^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];       

%short_plot_name=['Total condensed water mass (kg) up to ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%        short_plot_name=['Max cloud mixing ratio(g kg^{-1}) up to ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%        short_plot_name=['Max ice number  (L^{-1}) up to level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];
%        short_plot_name=['Max snow number  (L^{-1}) up to level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC '];



%         tit(1).tit=['Vertical velocity at level ' num2str(ih_wrf) ' for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (m s^{-1})'];

%       short_plot_name=['Height (m) abv. terrain for level ' num2str(ih_wrf) ' for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];

%        short_plot_name=['Pressure at level ' num2str(ih_wrf) ' for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa)'];
%        short_plot_name=['Pressure at ' num2str(h_wrf) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa)'];
%        short_plot_name=['Wind speed at ' num2str(h_wrf) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa)'];
%        short_plot_name=['Wind speed at ' num2str(p_wrf) ' hPa for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa)'];
%        short_plot_name=['Wind speed at ' num2str(potemp_wrf) ' K for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa)'];

%        short_plot_name=['Pressure difference at ' num2str(h_wrf) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa)'];


%short_plot_name=['Surface pressure difference (hPa) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC']; 

%short_plot_name=['Potential temperature gradient between ' num2str(inversion_height1/1000) ' km and ' num2str(inversion_height2/1000) ' km for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (K km^{-1})'];


    tit(1).tit=short_plot_name;

%tit(1).tit=['Geopotential Height (m) '];
%short_plot_name=['Snow depth change (m) for ' filestr];
%short_plot_name=['Snow mass change (kg m^{-2})']; 

%short_plot_name=['ALBEDO12M']; 


%short_plot_name=['Downward SW flux (W m^{-2})' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]; tit(1).tit=short_plot_name;
%short_plot_name=['TOA OLR flux (W m^{-2})' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]; tit(1).tit=short_plot_name;
%short_plot_name=['LH+SW flux (W m^{-2})' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]; tit(1).tit=short_plot_name;



tit(1).tit=[short_plot_name ' for ' filestr];
short_plot_name = tit(1).tit; %might want to comment this out if filename is too long


    
        %   tit(1).tit=['ECMWF mass weighted mean temperature over lower 2.5 km at ' tlab ' UTC, 25th Feb (^{o}C)'];
        
        
        

    


        
        
        
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
%        ylabelstr='iLatitude';
%        xlabelstr='iLongitude';
        
        ylabelstr='y (km)';
        xlabelstr='x (km)';
     
    
    case 'wrf_wind2d'
        
    iadd_terrain=1;
	iplot_latlon=1;
	iadd_wind_quivers=1;
    iadd_flight_path=1;
    
    icolmap=1;
    lb_map=lbmap(256,'brownblue'); %nice colormap for colorblind people
    isquare=1;
    
        i_reverse_cmap=1;
        if i_reverse_cmap==1
            cmap=flipud(lb_map);
        else
            cmap=lb_map;
        end
        
        switch file_type
            case 'wrfout'
                time=6; %13=12UTC
                time=idir;
            case 'met_em'
                time=1;
        end
        
        ih_wrf=4; %z-level to plot for
        ih_wrf=11; %z-level to plot for
%        ih_wrf=18; %z-level to plot for
        
        iylim=0;
        ylims=[150 450];
        ixlim=0;
        xlims=[400 700];
        
%        xlims=[200 650];
%        ylims=[250 700];
      
        switch file_type
            case 'wrfout'
                Zlev=WRFUserARW(nc,'Z',time,ih_wrf);
                terrain = nc{'HGT'}(time,:,:);
                medZ=['(median of ' num2str( median(median((squeeze(Zlev) - terrain))) ,3 ) ' m above terrain)'];
            case 'met_em'
                if prod(size(nc{'GHT'})) > 0
                    Zlev=nc{'GHT'}(time,ih_wrf,:,:);
                    %terrain = nc{'HGT_M'}(time,:,:);
                    medZ=['(mode of ' num2str( mode(mode((squeeze(Zlev)))) ,3 ) ' m above terrain)'];
                else
                    Plev=nc{'PRES'}(time,ih_wrf,:,:)/100;
                    medZ=['(median of ' num2str(median(median(Plev)),'%.0f') ' hPa)'];
                end
            otherwise
                medZ='';
        end
        
        
        tit(1).tit=['WRF wind speed (m s^{-1}) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC at model level ' num2str(ih_wrf) ' ' medZ ' for ' filestr];
        short_plot_name=tit(1).tit;
        %   tit(1).tit=['ECMWF mass weighted mean temperature over lower 2.5 km at ' tlab ' UTC, 25th Feb (^{o}C)'];
        
        imaxovr=1;
        iminovr=1;
        
        mincovOvr = 0;
        maxcovOvr = 18;
        clab=0;
        clines=0;
        icont=1;
        
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
        ylabelstr='y (km)';
        xlabelstr='x (km)';
        
        sig=3;
        
        clab=0;

	
    

        
    case 'ecmwf_surf'
        
        switch it
        case 1
            tlab='21:00 UTC, 24th Feb';
        case 2
            tlab='0 UTC, 25th Feb';
        end
        
        tit(1).tit=['ECMWF mass weighted mean water vapour over lower 2.5 km at ' tlab ' (g kg^{-1})'];
        %   tit(1).tit=['ECMWF mass weighted mean temperature over lower 2.5 km at ' tlab ' UTC, 25th Feb (^{o}C)'];
        
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = -2;
        maxcovOvr = 2;
        clab=0;
        clines=0;
        icont=1;
        
        izovr=2;
        i2d=3;
        
        
        idirstamp=0;
        itimestamp=0;
        
        ylabelstr='Latitude';
        xlabelstr='Longitude';
        
        sig=3;
        
        clab=0;
        
    case 'cdensity'
        tit(1).tit='C value (density *pi/6) (kg m^{-3})';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = -2;
        maxcovOvr = 2;
        clab=0;
        clines=0;
        icont=1;
        
    case 'w_3d'
        tit(1).tit='Vertical Velocity (m s^{-1})';        
        imaxovr=1;
        iminovr=1;
        
        mincovOvr = -2;
        maxcovOvr = 2;
        clab=0;
        clines=0;
        icont=1;
        
    case 'vap_3d'
        tit(1).tit='Water Vapour (ppmv)';
         tit(1).tit='Vertical velocity (m s^{-1})';
        % tit(1).tit='Total ice mixing ratio (ppmv)';
%            tit(1).tit='Supersaturation (%)';  %run wrap_slice using 'ice supersat3' flag after loading in pressure_hslice etc.
            %from diags
        % tit(1).tit=['Tracer at ' num2str(GridDan(1).Z(ih)/1000+0.62,4) ...
        %   ' km (kg^{-1})'];
        % tit(1).tit=['Tot water at ' num2str(GridDan(1).Z(ih)/1000+0.62,4) ...
        %    ' km (ppmv)'];
        
        dlogflag=0;
        dlogmin=10;
        
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = -60;
        maxcovOvr = 2;
        
        %        mincovOvr = -0.5;
        %		maxcovOvr = 0.5;
        
        clab=0;
        clines=0;
        icont=1;  
        
        izovr=2;
        i2d=3;
        
        xlabelstr='X (km)';
        ylabelstr='Y (km)';
        
        isquare=1;
        i3d=1;
        idirstamp=0;
        
        
        
    case 'vap_3d_vert'
        i3d=1;
        tit(1).tit='Water vapour mixing ratio (ppmv)';
                tit(1).tit='Total ice mixing ratio (ppmv)';
                tit(1).tit='Potential temperature (K)';
        %        tit(1).tit='Vertical velocity (m s^{-1})';
        %        tit(1).tit='Tracer mixing ratio (g kg^{-1})';
        
        
        dlogflag=0;
        dlogmin=5;
        
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 25;
        clab=1;
        clines=1;
        icont=1;  
        
        ncont=25;
        
        sig=3;
        
        
        
    case 'vertvel'
        tit(1).tit='Vertical Velocity (m s^{-1})';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = 2;
        maxcovOvr = 20;
        clab=0;
        clines=0;
        icont=1;
        
    case 'lnb'
        tit(1).tit='Level of Neutral Buoyuancy (km)';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = 14.000000;
        maxcovOvr = 17.00000;
        clab=0';
        clines=0;
        icont=0;
        
    case 'si'
        tit(1).tit='Supersaturation wrt ice (%)';
        imaxovr=1;
        iminovr=0;
        
        maxcovOvr=40;
        mincovOvr=0;
        clab=0;
        clines=0;
        
        vectorf=1;
        
        
    case 'icesatMR'
        tit(1).tit='Ice saturation mixing ratio (ppmv)';
        imaxovr=0;
        iminovr=0;
        
        maxcovOvr=40;
        mincovOvr=0;
        clab=0;
        clines=0;
        
        vectorf=0;
        
        dlogflag=1;
        dlogmin=1;
        
    case 'potemp'
        normcbar=1;
        icolmap=1;
        cmap='hsv';
        tit(1).tit='Potential Temperature (K)';
        dlogflag=0;
        dlogmin=0.001;
        
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = 330.000000;
        maxcovOvr = 474.000000;
        
        %   maxcovOvr=750;
        
        maxcovOvr= 400;
        
        %        mincovOvr = dlog(330.000000,dlogmin);
        %		maxcovOvr = dlog(474.000000,dlogmin);
        
        clab=1;
        clines=1;
        ncont=75;
        
        
        ncont=25;
        clines=0;
        clab=0;
        normcbar=0;
        icolmap=0;
        
        sig=3;
        
        icont_extra=0;
        
        %ncont=100;
    case 'wind'
        tit(1).tit='Vertical Wind Speed (m/s)';
    case 'htracer'
        tit(1).tit='Height Dependent Tracer';
    case 'vapour'
        tit(1).tit='Vapour Mixing Ratio (ppmv)';
        %		tit(1).tit='Ice Mixing Ratio (ppmv)';
        imaxovr=1;
        iminovr=0;
        
        mincovOvr=1;
        maxcovOvr=8.5;
        
        %         maxcovOvr=25;
        % maxcovOvr=15;
        % 
        %         mincovOvr=3;        
        maxcovOvr=25;;
        
        dlogflag=0;
        dlogmin=0.1;
        
        %  mincovOvr=dlog(0.01,dlogmin);
        %  maxcovOvr=dlog(8.5,dlogmin);
        
        
        %         mincovOvr = dlog(3.600000,dlogmin);
        % 		maxcovOvr = dlog(540.000000,dlogmin);
        
        
        clines=0;
        clab=0;
        
        vectorf=0;
        
        icont_extra=0;
        
        
        
    case 'ozone'
        tit(1).tit='Ozone Mixing Ratio (ppmv)';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr=8;
        maxcovOvr=60;
        
    case 'lowtracer'
        tit(1).tit='Low level tracer (g kg^{-1})';
        
        dlogflag=1;
        dlogmin=1e-3;
        clines=0;
        
        imaxovr=0;
        iminovr=0;
        
        mincovOvr=0;
        maxcovOvr=0;
        
        
    case 'totwater'
        tit(1).tit='Total Water Mixing Ratio (ppmv)';
        imaxovr=1;
        iminovr=1;
        
        mincovOvr=1;
        maxcovOvr=8.5;
        %        mincovOvr=3;
        
        % mincovOvr = dlog(3.600000,dlogmin);
        %	maxcovOvr = dlog(540.000000,dlogmin);
        
        %      dlogflag=1;
        dlogmin=1e-2;
        
        %        maxcovOvr=15;
        
        ncont=25;
        
        clines=0;
        clab=1;
        
    case 'general'
        tit(1).tit='Total Condensate (g kg^{-1})';
%         tit(1).tit='Total ice mixing ratio (g kg^{-1})';
% %        tit(1).tit='Vapour mixing ratio (g kg^{-1})';
         tit(1).tit='Liquid mixing ratio (g kg^{-1})';
         tit(1).tit='Liquid mixing ratio (g m^{-3})';
                  
%                  tit(1).tit='Ice mixing ratio (g m^{-3})';
                  

%         tit(1).tit='Snow+graupel mixing ratio (g kg^{-1})';
%          tit(1).tit='Ice number concentration (L^{-1})';

         %      tit(1).tit='Graupel mixing ratio (g m^{-3})';

%         tit(1).tit='Liquid mixing ratio (g m^{-3})';

         %         tit(1).tit='Horizontal Wind (perpendicular to slice) (m s^{-1})';
%         tit(1).tit='Horizontal Wind (parallel to slice) (m s^{-1})';
%         tit(1).tit='Vertical Wind (m s^{-1})';
%         tit(1).tit='Temperature Perturbation (K)';
%         tit(1).tit='Pressure (hPa)';
%         
        dlogflag=0;
        dlogmin=200;
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(3.600000,dlogmin);
        maxcovOvr = dlog(540.000000,dlogmin);
        
        maxcovOvr = 25;
        
        vectorf=0;
        
        %        ncont=8;
        clab=0;
        clines=1;
        ncont=20;
        
    case 'inc'
        tit(1).tit='Ice Number Concentration (mg^{-1})';
        clines=0;
        dlogflag=0;
        dlogmin=1e-2;
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(3.600000,dlogmin);
        maxcovOvr = dlog(540.000000,dlogmin);
        
        vectorf=0;
        
    case 'temppert'
        tit(1).tit='Temperature Perturbation (K)';
        tit(1).tit='Temperature (K)';
        
        imaxovr=1;
        iminovr=0;
        
        mincovOvr=0;
        maxcovOvr=5;
        
        clines=0;
        
    case 'rhopert577'
        tit(1).tit='Density Perturbation (kg m^{-3})';
        clines=0;
        
    case 'hydbal'
        tit(1).tit='Upwards Non-hydrostatic Force Term (N m^{-3})';
        clines=0;
        dlogflag=0;
        dlogmin=2;
        
        clab=0;
        
    case 'dpdz'
        tit(1).tit='dpdz perturbation (N m^{-3})';
        clines=0;
        dlogflag=0;
        dlogmin=2;
        
        clab=0;
        
    case 'rhog'
        tit(1).tit='rho x G perturabtion (N m^{-3})';
        clines=0;
        dlogflag=0;
        dlogmin=2;
        
        clab=0;   
        
    case 'rad'
        tit(1).tit='Radiative Heating (kg day^{-1})';
        clines=1;
        clab=1;
        
    case 'radar'
        tit(1).tit='Radar reflectivity (dBZ)';
        clines=1;
        clab=0;   
        
        imaxovr=1;
        iminovr=1;
        
        mincovOvr=10;
        maxcovOvr=70;
        
        ncont=13;
        
        i3d=0;
        iabc=0;
        
        case 'ARM_radar'
            tit(1).tit='dbZ Echo Top (km)';
            clines=0;
            clab=0;   

            imaxovr=1;
            iminovr=1;

            mincovOvr=0;
            maxcovOvr=18;  

            izovr=2;

            idirstamp=0;
            ylabelstr='Distance (km)';
            timlab=fileNAME(26:29);

%--------------------------
end   %end of switch i577
%--------------------------    
    
    % tit(2).tit='Ice Number Concentration (mg^{-1})';
    nplots2d=1;
    
    
    
    
    
    
    if izovr==0 
        z=GridDan(idir).Z;
    end
    
    notsame=1;
    
    
    %    iminovr=[0 0];
    %    imaxovr=[0 0];
    
    
    figlab=[tit(1).tit ' 2d plot ' i577];

    
    if ~exist('savedir')
        savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savename=[savedir short_plot_name];   % ' ' dirname];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
case 56
    top_down_cumulative
case 55
    upflux_7thSep2005
    
case 54
    fallspeed_7thSep2005
    
case 53
    meanIce_7thSep2005
    
case 52
    microrate_7thSep2005; %microphysical source rate of vapour
    
case 51
    %idir=1;
    
    logflag=0;
    fact=1e6*28.97/18;
    
    %iminovr=1;
    mincovOvr=(3.5);
    
    imaxovr=1;
    maxcovOvr=(1);
    %tit(1).tit='Mean Vapour Change (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
    %tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    tit(1).t='Mean Tot Water time-height';
    figlab=tit(1).tit;
    
    minZ=14.7e3;
    maxZ=22e3;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    %i2d=2; %tells it to label x axis in km
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
    sig=1;
    
case 50
    prc=15; %percentile index
    
    manclab=0;
    
    idirstamp=1;
    icont_extra=0;
    
    dlogflag=0;
    logflag=0;
    fact=1e6*28.97/18;
    
    iminovr=0;
    %    mincovOvr=-0.1;
    %mincovOvr=3.8;
    
    imaxovr=0;
    %    maxcovOvr=0.1;
    
    iflux=0;
    
    ixlim=1;
    xlims=[0 1.75];
    
    
    
    %	tit(1).tit='Mean Vapour (ppmv)';
    %    tit(1).tit='Mean Total Water (ppmv)';
    %    tit(1).tit='Mean Lower Tracer Mixing Ratio (ppmv)';
    %     tit(1).tit='Microphysical Ice Number Source (kg^{-1} s^{-1})';
    %     tit(1).tit='Microphysical Ice Source (kg/kg s^{-1})';
    %    tit(1).tit='Min Vapour (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
    %tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    
    hrange=16;
    switch hrange
    case 1
        minZ=13e3;
        maxZ=22e3;
    case 2
        minZ=14.6e3;
        maxZ=17e3;
        %  maxZ=17e3;
    case 3
        minZ=0.2e3;
        maxZ=19e3;
    case 4
        minZ=0.2e3;
        maxZ=22e3;
    case 5
        minZ=0e3;
        maxZ=30.4e3;
    case 6
        minZ=15e3;
        maxZ=20e3;
    case 7
        minZ=14e3;
        maxZ=22e3;
    case 8
        minZ=14e3;
        maxZ=19e3;
    case 9
        minZ=0e3;
        maxZ=4e3;    
    case 10
        minZ=13e3;
        maxZ=18.6e3;
    case 11
        minZ=15e3;
        maxZ=30.4e3;    
    case 12
        minZ=15e3;
        maxZ=16.5e3;  
    case 13
        minZ=11e3;
        maxZ=17e3;  
    case 14
        minZ=0e3;
        maxZ=50.4e3; 
    case 15
        minZ=0e3;
        maxZ=15e3;    
    case 16
        minZ=3e3;
        maxZ=15e3;         
    end
    
    
    s50='adrate';
    s50='fallrate';
    s50='change'; %mean change in total water from initial
    %s50='change_from_dqtot'; %as above but calc'd from the dqtot and nntot & averaged over 1000km (to allow comparison with different domain sizes)
    %s50='topdowncum';
    %s50='fallflux';
    %s50='icedist';
    %s50='icendist';
    
    %s50='vapad';
    %s50='micronc';
    % %s50='fallnc';
    %s50='adnc';
    % %s50='changenc';
    %s50='icead';
    %s50='changevap';
    %s50='meanvap';
%    s50='changeice';
    %s50='icemass';
%    s50='minvap';
    %s50='mintot'; %also does the max
    % % % %s50='fall+ad';
    %s50='microice';
%    s50='rain';
%     s50='liq';
    s50='ice';
    s50='maxice';   %plus maximums for other HMs

    % % 
    %s50='allice';
    % % 
    % % % % % s50='snow';
    %s50='graupel';
    % % % % % 
%    s50='maxw';
    % s50='minw';
    % % %
    %s50='iceno';
    % % % % s50='snowno';
    % % % % s50='grano';
    % % % % 
    % % % % s50='minvap';
    % % % 
    %s50='mphys_process';
    %s50='PGMLT';
    % % % s50='racw';
    %s50='praut';
    % % % s50='pifrw';
    % % % s50='allpr';
    % % % %s50='piacw';
    %s50='pidep';
    %s50='pIdep';
    %s50='pIsub'; %subimation of all ice
    %s50='pisub';
    % % % s50='dqi';
    % % % s50='pgsub';
    %s50='prevp';
    %s50='dql';
    % % s50='pcond';
    % %  s50='dq_potemp';
    % %  s50='dq_non';
    %s50='dqtot';
    %s50='meanLT5tot'_;
    %s50='dqvap';
    %s50='dqvap_dist';
    %s50='dqvap_dist_abv';
    %s50='dqtot_dist_abv';
    
    
    %s50='nntot';
    %s50='nnvap';
    % %   s50='ratio_potemp';
    % %   s50='change_conv_potemp';
    % %   s50='combined_potemp';
    %s50='low_tracer';
    %s50='maxlowtracer';
    %s50='tracerflux';
    
    % 
    % s50='iceadcum';
    % s50='icefallcum';
    % s50='icemicrocum';
    % % s50='picesubcum';
    %s50='totadcum';
    %s50='vapadcum';
    
    % s50='si'; %supersat wrt ice
    % s50='si_diag'; %supersat wrt ice
    
    %s50='rhopert'; 
    % s50='drhodz';
    % 
    %s50='upflux';
%    s50='meanw';
    
    %s50='lnbbel';
    %s50='lnbdist';
    %s50='lnbdist_tot';
    %s50='lnbdist_vap';
    %s50='rad';
    %s50='lwrad';
    %s50='swrad';
    %s50='minTguess';
    %s50='vapdist';
    %s50='totdist';
    %s50='meanTemp';
    %s50='radar10dbz';
    %s50='mean_rho';
    %s50='radar_ndbzARM';
    
    %s50='mass flux';
    %s50='lwc_width';
    
    %s50='dT_conv';
    
    %pdat(i).p=length(GridDan(idir).Y1)*( dq_tot(idir).d(izmin:izmax,dumprange,1) ) *dy; %multiply br dy so is in ppmv*km since otherwise high res will mean there are more
    
    
    sig=2;
    clines=1; %makes black contour lines appear
    clab=0;
    
    %ncont=19;
    %points counted
    dgs={'PGDEP','PGMLT', ...  %1-2
            'PRAUT',   'PGSHD',   'PRACW',   'PSMLT', ... %3-6
            'PIMLT',   'PSAUT',   'PSDEP',   'PIACR_G', ... %7-10
            'PSACI',   'PGACI',   'PGACW',   'PGACS', ... %11-14
            'PGACR',   'PSACR',   'PRACS',   'PSACW', ... %15-18
            'PGAUT',   'PGFR',    'PRACI_G', 'PGWET', ... %19-22
            'PGDRY',   'PGSUB',   'PSSUB',   'PREVP', ... %23-26
            'PISUB',   'DQI  ',   'PIHAL',   'PIPRM', ... %27-30
            'PIDEP',   'PIACW',   'PICNT',   'PIFRW', ... %31-34
            'PIACR_S', 'PRACI_S', 'PRACI',   'PIACR', ... %35-38
            'RIACR_S', 'RSACR',   'RIACR_G', 'RGFR', ...  %39-42
            'RGACS',   'RGACR',   'RSAUT',   'RIACI', ... %43-46
            'RSACS',   'RSBRK'  ...                       %47-48
        };
    switch s50
    case 'dT_conv'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='dT conv';
        %        tit(1).tit='dT bubble';
        %        tit(1).tit='dT nonconv';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 20;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'radar10dbz'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='LWC cloud width (km)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 20;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'mass flux' 
        iflux=0;
        
        dlogflag=1;
        dlogmin=1e-3;
        % dlogmin=1e-9;
        
        tit(1).tit='Mass flux';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = -30;
        maxcovOvr = 20;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(1e-4,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
        
    case 'radar_ndbzARM' 
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='Number of points with radar echo above 10 dBZ';
        
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0;
        maxcovOvr = 1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 10;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
        izovr=1;
        
    case 'mean_rho'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        tit(1).tit='Mean Density in points with total water less than 5 ppmv (kg m^{-3})';
        
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 20;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'radar10dbz'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='Radar 20 dBZ length (km)';
        tit(1).tit='Radar 40 dBZ length (km)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 20;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'meanTemp'
        iflux=0;
        
        dlogflag=0;
        dlogmin=170;
        
        tit(1).tit='Mean potential temp in updraughts (K)';
        tit(1).tit='Mean temp change (K)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 190;
        maxcovOvr = 0;
        
        % mincovOvr = 0.000000;
        maxcovOvr = 0.460000;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
        
    case 'totdist'
        iflux=0;
        
        dlogflag=1;
        dlogmin=1e-4;
        
        tit(1).tit='Total water frequency distribution (ppmv^{-1})';
        tit(1).tit='Total water distribution (ppmv^{-1})';
        
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.460000;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
        ylabelstr='Total water mixing ratio (ppmv)';
        izovr=2;
        
        icont=0;
    case 'vapdist'
        iflux=0;
        
        dlogflag=1;
        dlogmin=1e-4;
        
        tit(1).tit='Vapour frequency distribution (ppmv^{-1})';
        tit(1).tit='Vapour frequency distribution (ppmv^{-1})';
        
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.460000;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.1,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
        ylabelstr='Vapour mixing ratio (ppmv)';
        izovr=2;
        
    case 'minTguess'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='Minimum Temperature (^{o}C)';
        tit(1).tit='Max Temperature Perturbation Estimate (^{o}C)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.460000;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'maxlowtracer'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='Max Low Tracer Value (g kg^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.460000;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'tracerflux'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        tit(1).tit='Low Tracer Flux';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.460000;
        
        %mincovOvr = dlog(0.000000,dlogmin);
        %maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;    
        
    case 'icedist'
        %run gamdistTimH first - makes ice dist from the icediagsALL averages rather than loading in 2d fields
        
        
        iflux=0;
        
        dlogflag=1;
        dlogmin=0.5e-10;
        dlogmin=1;
        
        H=15.25;
        
        iz=findheight(GridDan(idir).Z+620,H*1000);
        
        tit(1).tit=['Cloud Ice dq dD^{-1} (ppmv micron^{-1}) at ' num2str(GridDan(idir).Z(iz)/1000+add_ground_height,4) ' km'];
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0,dlogmin);
        maxcovOvr = dlog(2,dlogmin);
        
        clines=0;
        clab=1;
        
        manclab=0;
        
        ylabelstr='Diameter (microns)';
        
        izovr=2;
        
        clab=0;
        
    case 'icendist'
        iflux=0;
        
        dlogflag=1;
        dlogmin=0.5e-10;
        
        tit(1).tit=['Total Ice dN dD^{-1} (kg^{-1} micron^{-1}) at ' num2str(GridDan(idir).Z(iz)/1000+add_ground_height,4) ' km'];
        
        iminovr=0;
        imaxovr=0;
        
        %        mincovOvr = 999;
        %		maxcovOvr = 999;
        clines=0;
        clab=1;
        
        manclab=0;
        
        ylabelstr='Diameter (microns)';
        
        izovr=2;
        
        clab=0;
        
    case 'rad'
        iflux=0;
        
        dlogflag=0;
        dlogmin=0.5e-1;
        
        tit(1).tit='Radiative Forcing (K day^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        %        mincovOvr = 999;
        %		maxcovOvr = 999;
        clines=0;
        clab=1;
        
        manclab=0;
        
        ylabelstr='Radiative Forcing (K day^{-1})';
        
    case 'lwrad'
        iflux=0;
        
        dlogflag=0;
        dlogmin=0.5e-1;
        
        tit(1).tit='LW Radiative Forcing (K day^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        %        mincovOvr = 999;
        %		maxcovOvr = 999;
        clines=0;
        clab=1;
        
        manclab=0;
        
        ylabelstr='LW Radiative Forcing (K day^{-1})';
        
    case 'swrad'
        iflux=0;
        
        dlogflag=0;
        dlogmin=0.5e-1;
        
        tit(1).tit='SW Radiative Forcing (K day^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        %        mincovOvr = 999;
        %		maxcovOvr = 999;
        clines=0;
        clab=1;
        
        manclab=0;
        
        ylabelstr='SW Radiative Forcing (K day^{-1})';
        
    case 'lnbdist_tot'
        iflux=0;
        
        dlogflag=0;
        dlogmin=0.5e-1;
        
        tit(1).tit='Relative Frequency of LNBs for Negatively Bouyant Air with LT 5ppmv Total Water (km^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = GridDan(1).Z(100)/1000+add_ground_height;   %GridDan(1).Z(105)/1000+add_ground_height;
        maxcovOvr = GridDan(1).Z(160)/1000+add_ground_height;
        
        clines=0;
        clab=0;
        
        izovr=2;
        ylabelstr='LNB (km)';
        
        
    case 'lnbdist_vap'
        iflux=0;
        
        dlogflag=0;
        dlogmin=0.5e-1;
        
        tit(1).tit='Relative Frequency of LNBs for Negatively Bouyant Air with LT 5ppmv Vapour Content (km^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = GridDan(1).Z(100)/1000+add_ground_height;   %GridDan(1).Z(105)/1000+add_ground_height;
        maxcovOvr = GridDan(1).Z(160)/1000+add_ground_height;
        
        clines=0;
        clab=0;
        
        izovr=2;
        ylabelstr='LNB (km)';
        
        
        
        
    case 'lnbdist'
        iflux=0;
        
        dlogflag=0;
        dlogmin=0.5e-1;
        
        tit(1).tit='Distribution of LNBs for Negativaly Bouyant Air';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = GridDan(1).Z(100)/1000+add_ground_height;   %GridDan(1).Z(105)/1000+add_ground_height;
        maxcovOvr = GridDan(1).Z(160)/1000+add_ground_height;
        
        clines=0;
        clab=0;
        
        izovr=2;
        ylabelstr='LNB (km)';
        
        
    case 'lnbbel'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        % tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean LNB for negatively buoyant air with tot water below 5ppmv (km)';
        %tit(1).tit='Minimum LNB (km)';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = GridDan(1).Z(100)/1000+add_ground_height;   %GridDan(1).Z(105)/1000+add_ground_height;
        maxcovOvr = GridDan(1).Z(160)/1000+add_ground_height;
        
        clines=1;
        clab=1;
        
        sig=3;
        
    case 'dqtot'
        iflux=0;        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Sum of Total Water Points Deficit Below 5 ppmv x Grid Length (ppmv km)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        %	maxcovOvr = 150.000000;
        
        mincovOvr = 0.000000;
        %	maxcovOvr = 240.000000;
        
        clines=0;
        clab=0;
        
    case 'dqvap'
        iflux=0;
        
        dlogflag=1;
        dlogmin=1;
        
        tit(1).tit='Sum of Vapour Deficit Below 5 ppmv x Grid Resolution (ppmv km)';
        %tit(1).tit='Mean Vapour Value for Points Below 5 ppmv (ppmv)';
        %  tit(1).tit='Sum of vapour mass below 5 ppmv (kg m^{-1})';
        
        iminovr=1;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 500.000000;
        
        mincovOvr = dlog(10.000000,dlogmin);
        maxcovOvr = 450.000000;
        
        clines=0;
        clab=0;     
        
    case 'dqvap_dist'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Sum of Vapour Deficit Below 5 ppmv x Grid Length (ppmv km)';
        %tit(1).tit='Mean Vapour Value for Points Below 5 ppmv (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 500.000000;
        
        mincovOvr = 0.000000;
        maxcovOvr = 450.000000;
        
        clines=0;
        clab=0;  
        
        izovr=1;
        
    case 'dqvap_dist_abv'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Sum of Vapour Deficit Above 5 ppmv x Grid Length (ppmv km)';
        %tit(1).tit='Mean Vapour Value for Points Below 5 ppmv (ppmv)';
        
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 500.000000;
        
        mincovOvr = 0.000000;
        maxcovOvr = 3;
        
        clines=0;
        clab=0;  
        
        izovr=1;
        
    case 'dqtot_dist_abv'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Sum of Vapour Deficit Above 5 ppmv x Grid Length (ppmv km)';
        %tit(1).tit='Mean Vapour Value for Points Below 5 ppmv (ppmv)';
        
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 500.000000;
        
        mincovOvr = 0.000000;
        maxcovOvr = 3;
        
        clines=0;
        clab=0;  
        
        izovr=1;
        
        
    case 'meanLT5tot'
        iflux=0;        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Mean mixing ratio of tot water points LT 5 ppmv and covering more than 150 km (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 150.000000;
        
        mincovOvr = 0.000000;
        maxcovOvr = 240.000000;
        
        clines=0;
        clab=0;
        
        sig=3;
        
    case 'dqvap'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Sum of Vapour Deficit Below 5 ppmv x Grid Length (ppmv km)';
        %tit(1).tit='Mean Vapour Value for Points Below 5 ppmv (ppmv)';
        
        iminovr=1;
        imaxovr=0;
        
        mincovOvr = 20;
        maxcovOvr = 500.000000;
        
        mincovOvr = 0.000000;
        maxcovOvr = 400.000000;
        
        clines=0;
        clab=0;
        
        
        
    case 'nntot'
        iflux=0;        
        dlogflag=0;
        dlogmin=1e-2;
        
        iminovr=1;
        imaxovr=0;
        
        tit(1).tit='Total 2-D Length of Total Water Points Below 5 ppmv (km)';
        
        mincovOvr = 0.000000;
        maxcovOvr = 150.000000;
        
        clines=0;
        clab=0;
        
    case 'nnvap'
        iflux=0;
        dlogflag=1;
        dlogmin=1;
        
        tit(1).tit='Total 2-D Length of Vapour Points Below 5 ppmv (km)';
        
        iminovr=1;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 250.000000;
        
        clines=0;
        clab=0;
        
        
    case 'meanw'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-10;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean updraught (m s^{-1})';
        tit(1).tit='Mean cloudy updraught (m s^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 1e-20;  %0.000000;
        maxcovOvr = 0.870000;
        
        %         mincovOvr = dlog(0.000000,dlogmin);
        %         maxcovOvr = dlog(0.420000,dlogmin);
        
        %  ncont=25;
        
        clines=0;
        clab=0;
        
        ixlim=1;
        xlims=[0.7 1.55];
        
    case 'drhodz'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Density Gradient with Height (kg m^{-4})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=25;
        
        clines=0;
        clab=0;
        
    case 'rhopert'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-4;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean Density Change for Points with Water Vapour Below 5 ppmv (kg m^{-3})';
        % tit(1).tit='Mean Density Change (kg m^{-3})';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -0.0040;
        maxcovOvr = 0.01;
        
        %        mincovOvr = dlog(0.000000,dlogmin);
        %        maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=12;
        % ncont=25;
        
        clines=1;
        clab=0;
        
        sig=4;
        
    case 'upflux'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean Upwards Mass Flux (kg m^{-2} s^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.870000;
        
        %         mincovOvr = dlog(0.000000,dlogmin);
        %         maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=15;
        
        clines=1;
        clab=1;
        
    case 'si'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Maximum Supersaturation wrt ice';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        %  ncont=25;
        
        clines=1;
        clab=1;
        
    case 'si_diag'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Maximum Supersaturation wrt ice';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = -99.999;
        
        %    mincovOvr = dlog(-100.000000,dlogmin);
        %    maxcovOvr = dlog(-99,dlogmin);
        
        %  ncont=25;
        
        clines=1;
        clab=1;   
        
        sig=4;
        
    case 'totadcum'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Cumulative Advective Source of Total Water (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=25;
        
        clines=1;
        clab=1;
        
        sig=1;
        
    case 'vapadcum'
        
        manclab=0;
        
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Cumulative Advective Source of Vapour (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=15;
        
        clines=0;
        clab=1;
        
        sig=1;
        
        
        
    case 'icemicrocum'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Cumulative Ice Microphysical Source (ppmv)';
        
        iminovr=1;
        imaxovr=0;
        
        mincovOvr = 0;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=15;
        
        clines=1;
        clab=1;
        
    case 'icefallcum'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Cumulative Fall Speed Ice Loss (ppmv)';
        
        iminovr=1;
        imaxovr=0;
        
        mincovOvr = -0.7;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(-0.2e-3,dlogmin);
        maxcovOvr = dlog(1600.000000,dlogmin);
        
        ncont=25;
        
        clines=1;
        clab=1;
        
    case 'iceadcum'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Cumulative Advective Ice Source (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = -0.7;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(-0.003200,dlogmin);
        maxcovOvr = dlog(1300.000000,dlogmin);
        
        ncont=25;
        
        clines=1;
        clab=1;
        
    case 'low_tracer'
        dlogflag=0;
        dlogmin=1e-2;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean Low Tracer (kg kg^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = -0.7;
        maxcovOvr = 0.05;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.750000,dlogmin);
        
        %         mincovOvr = dlog(-0.000027,dlogmin);
        % 		maxcovOvr = dlog(0.000061,dlogmin);
        
        %maxcovOvr = 0.05;
        
        
        
        ncont=25;
        clines=1;
        clab=0;
        
        iflux=1;
        
    case 'combined_potemp'
        dlogflag=0;
        dlogmin=1e-6;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Reversible and Non-reversible effect combined (ppmv)';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -0.7;
        maxcovOvr = 0.05;
        
        %   mincovOvr = dlog(-3.000000,dlogmin);
        %	maxcovOvr = dlog(400.00000,dlogmin);
        
        ncont=25;
        clines=1;
        clab=1;
        
    case 'change_conv_potemp'
        dlogflag=0;
        dlogmin=1e-6;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Estimated Contribution of Convective Motion to Total Water Change (ppmv)';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -0.4;
        maxcovOvr = 0.1;
        
        %   mincovOvr = dlog(-3.000000,dlogmin);
        %	maxcovOvr = dlog(400.00000,dlogmin);
        
        ncont=15;
        clines=1;
        clab=1;
        
    case 'ratio_potemp'
        dlogflag=0;
        dlogmin=1e-6;
        
        %        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Contribution of Non-reversbile to changes from non and reversibile changes combined';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0;
        maxcovOvr = 1.00000;
        
        %   mincovOvr = dlog(-3.000000,dlogmin);
        %	maxcovOvr = dlog(400.00000,dlogmin);
        
        ncont=15;
        clines=0;
        clab=1;
        
    case 'dq_potemp'
        dlogflag=0;
        dlogmin=1e-6;
        
        tit(1).tit='Change in Initial Total Water Mixing Ratio due to Reversible Motions (ppmv km)';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -0.700000;
        maxcovOvr = 0.050000;
        
        %   mincovOvr = dlog(-3.000000,dlogmin);
        %	maxcovOvr = dlog(400.00000,dlogmin);
        
        ncont=15;
        clines=0;
        clab=1;
        
    case 'dq_non'
        dlogflag=0;
        dlogmin=1e-1;
        
        tit(1).tit='Change in Initial Total Water Mixing Ratio due to Non-Reversible Motions (ppmv km)';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -200.000000;
        maxcovOvr = 0.0000;
        
        %  mincovOvr = dlog(-5,dlogmin);
        %	maxcovOvr = dlog(0,dlogmin);
        
        ncont=15;
        clines=0;
        clab=1;
        
    case 'pcond'
        dlogflag=1;
        dlogmin=1e-2;
        
        tit(1).tit='Condensation Source of Liquid (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -2.000000;
        maxcovOvr = 1.900000;
        
        mincovOvr = dlog(-3.000000,dlogmin);
        maxcovOvr = dlog(1.900000,dlogmin);
        
        ncont=15;
        clines=0;
        
        
    case 'dql'
        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Microphysical Source of Liquid (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = -1.000000;
        maxcovOvr = 1.00000;
        
        %  mincovOvr = dlog(-2.000000,dlogmin);
        %  maxcovOvr = dlog(2.000000,dlogmin);
        
        ncont=15;
        clines=0;
        
    case 'prevp'
        
        
        dlogflag=0;
        dlogmin=1e-3;
        
        tit(1).tit='Evaporation of Rain (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.300000;
        
        %  mincovOvr = dlog(-0.079000,dlogmin);
        %	 maxcovOvr = dlog(0.550000,dlogmin);
        
        ncont=15;
        clines=0;
        
    case 'pgsub'
        
        
        dlogflag=0;
        dlogmin=1e-3;
        
        tit(1).tit='Sublimation of graupel (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.390000;
        
        %  mincovOvr = dlog(-0.079000,dlogmin);
        %	 maxcovOvr = dlog(0.550000,dlogmin);
        
        ncont=15;
        clines=0;
        
    case 'dqi'
        
        
        dlogflag=1;
        dlogmin=1e-6;
        
        tit(1).tit='Microphysical source of ice (ppmv s^{-1})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = -0.079000;
        maxcovOvr = 0.541000;
        
        mincovOvr = dlog(-0.120000,dlogmin);
        maxcovOvr = dlog(1.000000,dlogmin);
        
        ncont=15;
        clines=0;
        
    case 'pisub'
        
        
        dlogflag=1;
        dlogmin=1e-5;
        
        tit(1).tit='Sublimation of ice (ppmv s^{-1})';
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.100000;
        
        ncont=12;
        clines=0;
        
    case 'pIsub'
        
        
        dlogflag=1;
        dlogmin=1e-5;
        
        tit(1).tit='Sublimation of all ice (ppmv s^{-1})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(4e-4,dlogmin);
        
        ncont=12;
        clines=0;
        
    case 'picesubcum'
        
        
        dlogflag=1;
        dlogmin=1e-2;
        
        tit(1).tit='Cumulative Sublimation of Ice (ppmv s^{-1})'; %all ice (i+s+g)
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.100000;
        
        ncont=12;
        clines=1;
        clab=1;
        
        
    case 'pidep'
        
        
        dlogflag=1;
        dlogmin=1e-4;
        
        tit(1).tit='Deposition of vapour onto ice (ppmv s^{-1})';
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.600000;
        
        mincovOvr = dlog(0.010000,dlogmin);
        maxcovOvr = dlog(1e-2,dlogmin);
        
        ncont=19;
        clines=0;
        
    case 'pIdep'
        
        
        dlogflag=1;
        dlogmin=1e-4;
        
        tit(1).tit='Deposition of vapour onto all ice (ppmv s^{-1})';
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.600000;
        
        mincovOvr = dlog(0.010000,dlogmin);
        maxcovOvr = dlog(1e-2,dlogmin);
        
        ncont=19;
        clines=0;        
        
    case 'piacw'
        
        
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit='Accretion of water by ice (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.044000;
        
        ncont=12;
        clines=0;
        
    case 'allpr'
        
        
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit=dgs{iallpr};
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.036000;
        
        ncont=12;
        clines=0;
        
        
    case 'pifrw'
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit='Homogeneous freezing of liquid to ice (ppmv s^{-1})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.036000;
        
        ncont=15;
        clines=0;
        
    case 'praut'
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit='Autoconversion of liquid to rain (ppmv s^{-1})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.036000;
        
        ncont=12;
        
    case 'racw'
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit='Accretion of Liquid by Rain (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.700000;
        
        ncont=12;
        
    case 'mphys_process'
        
        pname='PRACW';
        pname='PRACW';
        pname='PRACW';
        pname='';
        %   pname='PGMLT';
        %    pname='PRAUT';
        
        dlogflag=0;
        dlogmin=0;
        
        for imp=1:length(dgs)
            if strcmp(dgs{imp},pname)==1
                imphys=imp;
                break
            end
        end
        
        
        tit(1).tit=[pname ' (ppmv s^{-1})'];
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.800000;
        
        ncont=12;    
        
    case 'PGMLT'
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit='Melting of Graupel (ppmv s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0.000000;
        maxcovOvr = 0.800000;
        
        ncont=12;
        
    case 'minvap'
        dlogflag=0;
        dlogmin=0;
        
        tit(1).tit='Minimum vapour mixing ratio (ppmv)';
        %     tit(1).tit='Maximum vapour mixing ratio (ppmv)';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = 0;
        maxcovOvr = 5.3;
        
        mincovOvr = 0;
        maxcovOvr = 15;
        
        clines=0;
        
        ncont=25;
        
    case 'iceno'
        dlogflag=1;
        dlogmin=1e-12;
        
        tit(1).tit='Ice Number Concentration (kg^{-1})';
        tit(1).tit='Ice Number Concentration in Cloudy Updraughts (m^{-3})';
        % tit(1).tit='Max Ice Number Concentration (m^{-3})';
        tit(1).tit='Mean ice number concentration x domain length (kg^{-1} km) x 10^9';
        
        iminovr=0;
        imaxovr=1;
        
        
        mincovOvr = dlog(3.800000,dlogmin);
        maxcovOvr = dlog(269999999.999999,dlogmin);
        
        mincovOvr = dlog(0,dlogmin);
        maxcovOvr = dlog(2e-3,dlogmin);
        
        clines=0;
        ncont=16;
        
    case 'grano'
        dlogflag=1;
        dlogmin=1;
        
        tit(1).tit='Grapuel Number Concentration (kg^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(1200.000000,dlogmin);
        
    case 'snowno'
        dlogflag=1;
        dlogmin=1;
        
        tit(1).tit='Snow Number Concentration (kg^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(13000.000000,dlogmin);
        
        % mincovOvr = 0.000000;
        %maxcovOvr = 15000.000000;
        
        
        
    case 'maxw'
        dlogflag=0;
        dlogmin=1;
        
        tit(1).tit='Maximum updraught (m s^{-1})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(480.000000,dlogmin);
        
        mincovOvr = 0.0;
        maxcovOvr = 8; %44.000000;
        
        
        clines=0;
        clab=1;
        
        ncont=12;
        ncont=25;
        
        %  minZ=0.2e3;
        %  maxZ=19e3;
        
%        ixlim=1;
%        xlims=[0 4];
        
        ixlim=1;
		xlims=[0.7 1.55];
		xlims=[0.7 2.4];
        xlims=[0.55 1.55];
        
    case 'minw'
        dlogflag=0;
        dlogmin=1;
        
        tit(1).tit='Maximum downdraught (ms^{-1})';
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(10.000000,dlogmin);
        
        %        mincovOvr = 0.0;
        %		maxcovOvr = 44.000000;
        
        
        clines=0;
        clab=1;
        
        ncont=12;
        ncont=20;
        

        
        
    case 'graupel'
        dlogflag=0;
        dlogmin=1;
        
        tit(1).tit='Graupel mixing ratio (ppmv)';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(480.000000,dlogmin);
        
        mincovOvr = 0.000000;
        maxcovOvr = 760.000000;
        
    case 'snow'
        dlogflag=0;
        dlogmin=1;
        
        tit(1).tit='Snow mixing ratio (ppmv)';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(480.000000,dlogmin);
        
        mincovOvr = 0.000000;
        maxcovOvr = 440.000000;
        
    case 'ice'
        dlogflag=0;
        dlogmin=1e-1;
        
%        tit(1).tit='Ice mixing ratio (g kg^{-1} km)';
        tit(1).tit='Mean ice mixing ratio (g kg^{-1})';        
         tit(1).tit='Mean total ice mixing ratio in cloudy air (g kg^{-1})';
         tit(1).tit='Mean snow mixing ratio in cloudy air (g kg^{-1})';
%         tit(1).tit='Mean ice mixing ratio in cloudy updraughts (g kg^{-1})';
		tit(1).tit='Mean tracer mixing ratio in cloudy air (g kg^{-1})';
      %  tit(1).tit='Mean tracer mixing ratio (g kg^{-1})';
        
        iminovr=0;
        imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(10.000000,dlogmin);
        
        %  mincovOvr = 0.000000;
        maxcovOvr = 250;
        
        ncont=15;
        clab=0;
        
         xlims=[0.55 1.55];
         
    case 'maxice'
        dlogflag=0;
        dlogmin=1e-1;
        
        tit(1).tit='Max ice mixing ratio (g kg^{-1})';        
        tit(1).tit='Max graupel mixing ratio (g kg^{-1})';        
%        tit(1).tit='Max snow mixing ratio (g kg^{-1})';
		tit(1).tit='Max total IWC mixing ratio (g kg^{-1})';
        tit(1).tit='Max tracer mixing ratio (kg^{-1})';
        tit(1).tit='Max rain mixing ratio (g kg^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(10.000000,dlogmin);
        
        %  mincovOvr = 0.000000;
        maxcovOvr = 1.500000;
        
        ncont=15;
        clab=0;
        
         xlims=[0.55 1.55];     

        
    case 'allice'
        dlogflag=0;
        dlogmin=1e-3;
        
        tit(1).tit='Total ice mixing ratio (ppmv km)';
        tit(1).tit='Mean ice mixing ratio (g m^{-3})';
        
        iminovr=0;
        imaxovr=0;
        
        %	mincovOvr = dlog(0.000000,dlogmin);
        %	maxcovOvr = dlog(670.000000,dlogmin);
        
        %  mincovOvr = 0.000000;
        %  maxcovOvr = 460.000000;
        clab=1;
        
    case 'rain'
        dlogflag=0;
        dlogmin=1e-5;
        
        tit(1).tit='Rain mixing ratio (g kg^{-1})';
        %  tit(1).tit='Mean rain water mixing ratio cloudy updraught regions (g m^{-3})';
        %  tit(1).tit='Mean rain water mixing ratio in W GT 1 m s^{-1} regions (g m^{-3})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.000120,dlogmin);
        
        mincovOvr = 0.000000;
        maxcovOvr = 3.500000;
        
        clines=0;
        
        ixlim=1;
        xlims=[0.7 1.55];
        xlims=[0.55 1.55];
        
    case 'liq'
        dlogflag=0;
        dlogmin=0;
        
        %tit(1).tit='Mean Liquid water mixing ratio in Cloudy Updraughts (g m^{-3})';
        %    tit(1).tit='Max Liquid water mixing ratio (g m^{-3})';
        %tit(1).tit='Mean Liquid water mixing ratio in W GT 1 m s^{-1} regions (g m^{-3})';
        tit(1).tit='Mean Liquid water mixing ratio in cloudy updraught regions (g m^{-3})';
        tit(1).tit='Mean Liquid water mixing ratio (g m^{-3})';
        
        %     tit(1).tit='Mean liquid water mixing ratio x domain size (ppmv km)';
         %    tit(1).tit='Mean liquid water mixing ratio x domain size (g kg^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        % 		mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = 60;
        %         
        %         
        %         
        % mincovOvr = 0.000000;
        %maxcovOvr = 1e5;
        %     
        ncont=15;
        
        clines=0;
        clab=0;
        
        ixlim=1;
 %       xlims=[0.75 1.75];
%         xlims=[0.7 1.7];
         xlims=[0.55 1.55];
         
    
    case 'minvap'
        dlogflag=0;
        dlogmin=0.5;
        
        tit(1).tit='Minimum vapour mixing ratio (ppmv)';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr=dlog(0.92,dlogmin);		
        maxcovOvr=dlog(34,dlogmin);
        
        mincovOvr = 0.770000;
        maxcovOvr = 5.200000;
        
    case 'mintot'
        dlogflag=0;
        dlogmin=2;
        
        tit(1).tit='Minimum total water mixing ratio (ppmv)';
        tit(1).tit='Maximum total water mixing ratio (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr=dlog(2.3,dlogmin);
        maxcovOvr=dlog(34,dlogmin);
        
        mincovOvr = 2.100000;
        maxcovOvr = 5.200000;
        
        clines=0;
        clab=0;
        
    case 'adrate' %rate of change of tot water due to advection as worked out from change and fall flux 
        iflux=1;
        tit(1).tit='Tot water source from advection (ppmv s^{-1})';
        
        iminovr=0;	
        imaxovr=0
        
        dlogflag=1;
        dlogmin=1e-3;
        
        mincovOvr = dlog(-0.700000,dlogmin);
        maxcovOvr = dlog(1.00000,dlogmin);
        
        
        switch hrange
        case 1
            dlogflag=1;
            dlogmin=1e-3;
            mincovOvr = dlog(-0.120000,dlogmin);
            maxcovOvr = dlog(0.450000,dlogmin);
        case 2
            dlogflag=1;
            dlogmin=1e-1;
            mincovOvr = dlog(-0.700000,dlogmin);
            maxcovOvr = dlog(3.200000,dlogmin);
        end
        
        dlogmin=1e-6;
        
        
        ncont=15;
        ncont=25;
        clines=1;
        
        clab=1;
        manclab=1;
        
        
    case 'fallrate'
        iflux=1;
        tit(1).tit='Fall Speed Flux Loss (ppmv s^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        dlogflag=1;
        dlogmin=1e-3;
        
        mincovOvr = dlog(-5.6e-6,dlogmin);
        maxcovOvr = dlog(0.21,dlogmin);
        
        switch hrange
        case 1
            dlogflag=1;
            dlogmin=1e-6;
            mincovOvr = dlog(-5.6e-6,dlogmin);
            maxcovOvr = dlog(0.21,dlogmin);
        case 2
            dlogflag=1;
            dlogmin=1e-7;
            mincovOvr = dlog(-0.000003,dlogmin);
            maxcovOvr = dlog(0.660000,dlogmin);
        end
        
        dlogmin=1e-7;
        mincovOvr = dlog(-0.000003,dlogmin);
        maxcovOvr = dlog(0.660000,dlogmin);
        
        clab=1;
        clines=1;
        
        manclab=0;
        
    case 'fallflux'
        iflux=1;
        tit(1).tit='Fall Speed Flux Loss Calculated from Ice Means (ppmv s^{-1})';
        
        iminovr=1;
        imaxovr=1;
        
        dlogflag=1;
        dlogmin=1e-3;
        
        mincovOvr = dlog(-5.6e-6,dlogmin);
        maxcovOvr = dlog(0.21,dlogmin);
        
        switch hrange
        case 1
            dlogflag=1;
            dlogmin=1e-6;
            mincovOvr = dlog(-5.6e-6,dlogmin);
            maxcovOvr = dlog(0.21,dlogmin);
        case 2
            dlogflag=1;
            dlogmin=1e-7;
            mincovOvr = dlog(-0.000003,dlogmin);
            maxcovOvr = dlog(0.660000,dlogmin);
        end
        
        mincovOvr = dlog(-0.000006,dlogmin);
        maxcovOvr = dlog(0.210000,dlogmin);
        
        clab=0;
        clines=0;    
        
    case 'change'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-2;
        dlogmin=1;
        
        dlogmin=0.1;
        
        tit(1).tit='Change in Total Water (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        %         mincovOvr = -1.200000;
        %         maxcovOvr = 10.000000;
        %         
        %         mincovOvr = dlog(-0.460000,dlogmin);
        % 		maxcovOvr = dlog(500.000000,dlogmin);
        %         
        %         mincovOvr = dlog(-0.4600,dlogmin);
        % 		maxcovOvr = dlog(203.000000,dlogmin);
        %         
        %       %  mincovOvr = dlog(-0.480400,dlogmin);
        % 	%	maxcovOvr = dlog(568.000000,dlogmin);
        %         
        %         
        %         %mincovOvr = 0.000000;
        %         %maxcovOvr = dlog(10,dlogmin);
        % 
        % 		maxcovOvr = dlog(1.1000000,dlogmin);
        %         
        % 		maxcovOvr = dlog(0.1000000,dlogmin);
        maxcovOvr = 0.1000000;
        maxcovOvr = 0.02000000;
        
        
        %		maxcovOvr = 20;
        
        %        maxcovOvr = dlog(0.02000000,dlogmin);
        
        mincovOvr = -0.45;
        %maxcovOvr = 10;
        
        maxcovOvr = dlog(2.000000,dlogmin);
        
        ncont=35;
                ncont=15;
        
        clines=0;
        clab=0;
        sig=2;
        
        manclab=0;
        
        icolmap=0;
        cmap=hsv;
        
    case 'change_from_dqtot'
        iflux=1;
        
        dlogflag=0;
        dlogmin=1e-1;
        
        dlogmin=0.1;
        
        tit(1).tit='Change in Total Water over 1000 km (ppmv)';
        
        iminovr=1;
        imaxovr=1;
        
        maxcovOvr = 0.1000000;
        mincovOvr = -0.6;
        
        
        ncont=35;
        
        clines=0;
        clab=1;
        sig=3;
        
        manclab=0;
        
        icolmap=1;
        cmap=hsv;    
        
    case 'topdowncum'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-1;
        
        tit(1).tit='Top down Cumulative Sum of Total Water Mass (kg)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = -1.200000;
        maxcovOvr = 10.000000;
        
        mincovOvr = dlog(-0.460000,dlogmin);
        maxcovOvr = dlog(500.000000,dlogmin);
        
        
        %mincovOvr = 0.000000;
        %maxcovOvr = dlog(10,dlogmin);
        
        ncont=35;
        
        clines=1;
        clab=1;
        sig=3;    
        
    case 'changevap'
        dlogflag=1;
        dlogmin=1e-3;
        dlogmin=1e-1;
        
        
        tit(1).tit='Change in Vapour (ppmv)';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(-0.183600,dlogmin);
        maxcovOvr = dlog(603.000000,dlogmin);
        
        mincovOvr = dlog(-1.00,dlogmin);
        maxcovOvr= dlog(2,dlogmin);
        
        
        clab=0;
        clines=0;
        
        manclab=0;
        
        iflux=1;
        
        icolmap=0;
        cmap=hsv;
        
    case 'meanvap'
        dlogflag=0;
        dlogmin=1e-3;
        
        tit(1).tit='Mean Vapour (ppmv)';
        iminovr=1;
        mincovOvr=dlog(-0.82,dlogmin);
        
        imaxovr=1;
        maxcovOvr=dlog(0.96,dlogmin);
        maxcovOvr=dlog(5.2,dlogmin);
        
        maxcovOvr=5.1;
        mincovOvr=4.99;
        
        clab=1;
        
        manclab=0;
        
        iflux=1;
        
        sig=4;
        
    case 'changeice'
        iflux=1;
        dlogflag=1;
        dlogmin=1e-1;
        dlogmin=1;
        
        tit(1).tit='Mean Total Ice (ppmv)';
        %tit(1).tit='Mean Total Ice (g kg^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(600.000000,dlogmin);
        
        mincovOvr = dlog(0.0,dlogmin);
        maxcovOvr = dlog(603.000000,dlogmin);
        
        maxcovOvr = dlog(0.300000,dlogmin);
        
        %  mincovOvr=0.05;
        %  maxcovOvr=0.3;
        
        clab=0;
        clines=1;
        
        manclab=0;   
        
        ncont=35;
        
        icolmap=1;
        cmap=hsv;
        
    case 'icemass'
        iflux=1;
        dlogflag=1;
        dlogmin=1e-3;
        %dlogmin=1e-1;
        
        tit(1).tit='Mean Cloud Ice Particle Mass (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.090000,dlogmin);
        
        clab=0;
        clines=0;
        
        manclab=0;   
        
        ncont=25;        
        
    case 'changerate'
        tit(1).tit='Rate of Change of Total Water (ppmv s^{-1})';
        iminovr=1;
        mincovOvr=-0.001;
        
        imaxovr=1;
        maxcovOvr=0.05;
        
    case 'microice'
        
        iflux=1;
        
        tit(1).tit='Microphysical Source of Ice (ppmv s^{-1})';
        
        dlogflag=0;
        dlogmin=10;        
        
        iminovr=1;
        imaxovr=0;
        
        mincovOvr = dlog(100,dlogmin);
        mincovOvr=100;
        maxcovOvr = dlog(50,dlogmin);
        
        
        
        switch hrange
        case 1
            dlogflag=1;
            dlogmin=1e-4;
            mincovOvr = dlog(-0.002200,dlogmin);
            maxcovOvr = dlog(0.038000,dlogmin);
            
            dlogmin=1e-5;
            maxcovOvr = dlog(1e-4,dlogmin);                
            mincovOvr = dlog(-1e-4,dlogmin);
            
            
        case 2
            dlogflag=1;
            dlogmin=1e-4;
            mincovOvr=dlog(-0.32,dlogmin);
            maxcovOvr=dlog(0.5,dlogmin);
            mincovOvr=dlog(-0.00038,dlogmin);
            maxcovOvr=dlog(0.00075,dlogmin);
            mincovOvr = dlog(-0.001400,dlogmin);
            maxcovOvr = dlog(0.004900,dlogmin);
        end
        
        %    dlogmin=1e-5;
        %    maxcovOvr = dlog(1e-4,dlogmin);                
        %    mincovOvr = dlog(-1e-4,dlogmin);
        
        
        
        
        clab=0;
        clines=0;
        
    case 'vapad'
        dlogflag=1;
        dlogmin=1e-8;
        
        tit(1).tit='Advective Source of Vapour (ppmv s^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr=-0.0006;
        maxcovOvr=0.0004;
        
        mincovOvr = -0.000480
        maxcovOvr = 0.003400
        
        mincovOvr = dlog(-0.000480,dlogmin);
        maxcovOvr = dlog(0.003400,dlogmin);
        
        mincovOvr = dlog(-0.240000,dlogmin);
        maxcovOvr = dlog(0.640000,dlogmin);
        
        mincovOvr = dlog(-0.002200,dlogmin);
        maxcovOvr = dlog(0.038000,dlogmin);
        
        clines=0;
        
    case 'icead'
        iflux=1;
        
        dlogflag=1;
        dlogmin=1e-5;
        
        tit(1).tit='Advective Source of Ice (ppmv s^{-1})';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(-0.240000,dlogmin);
        maxcovOvr = dlog(0.640000,dlogmin);
        
        clines=0;
        clab=0;
        
    case 'micronc'
        dlogflag=1;
        dlogmin=100;
        
        tit(1).tit='Microphysical Source of Ice Number (kg^{-1}s^{-1})';
        iminovr=0;
        imaxovr=0;
        
        mincovOvr=dlog(-7000,dlogmin);
        maxcovOvr=dlog(2500,dlogmin);
        
        clines=0;
        clab=0;
        
    case 'adnc'
        iflux=1;
        
        dlogflag=1;
        dlogmin=100;
        
        tit(1).tit='Advective Source of Ice Number (kg^{-1}s^{-1})';
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(-10000.000000,dlogmin);
        maxcovOvr = dlog(53000.000000,dlogmin);
        
        clines=0;
        
    case 'fallnc'
        iflux=1;
        
        dlogflag=1;
        dlogmin=100;
        
        tit(1).tit='Fall Speed Flux Source of Ice Number (kg^{-1}s^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr=dlog(-7e5,dlogmin);
        maxcovOvr=dlog(1.2e5,dlogmin);
        
        clab=0;
        
    case 'changenc'
        iflux=1;
        
        dlogflag=1;
        dlogmin=100;
        
        tit(1).tit='Mean Ice Number (g^{-1})';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr=dlog(-7e5,dlogmin);
        maxcovOvr=dlog(5.5e7,dlogmin);
        
        clab=0;
        clines=0;
        
    case 'fall+ad'
        dlogflag=1;
        dlogmin=1e-5;
        
        tit(1).tit='Fall speed + advective source of ice mixing ratio (ppmv s^{-1})';
        
        iminovr=1;   
        imaxovr=1;
        
        
        switch hrange
        case 1
            mincovOvr=dlog(-0.43,dlogmin);
            maxcovOvr=dlog(1.73,dlogmin);
        case 2
            dlogflag=1;
            dlogmin=1e-4;
            mincovOvr=dlog(-0.32,dlogmin);
            maxcovOvr=dlog(0.5,dlogmin);
            mincovOvr=dlog(-0.0031,dlogmin);
            maxcovOvr=dlog(0.0031,dlogmin);
            mincovOvr=dlog(-0.0075,dlogmin);
            maxcovOvr=dlog(0.0075,dlogmin);
        end
        
        
    end
    %     tit(1).tit='Max Radar Echo (dBZ)';
    %     mincovOvr=-20;
    %    maxcovOvr=70;
    
    %prcstr=num2str(100-(21-prc)*5);
    %tit(1).tit=[prcstr 'th Percentile Updraught (m/s)'];
    
    
    figlab=[tit(1).tit ' Time-height '];
    
%    savename=figlab;
    %    minZ=12e3;
    
    % minZ=620;
    %  minZ=14.7e3;
    
    
    %minZ=16.0e3;
    
    
    %i2d=2; %tells it to label x axis in km
    
    if izovr==0
        dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
        z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
        if iutc==1
            time=GridDan(idir).t;
        else
            time=GridDan(idir).t;
        end
        rho=repmat(GridDan(idir).RHON(izmin:izmax),[1 length(dumprange)]);
    end
    
    
    
    
case 49
    %fact=1e6*28.97/18;
    logflag=0;
    tit(1).tit='Temp (^oC)';
    tit(2).tit='Ice NC (mg^{-1}';
    nplots2d=2;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    i2d=2; %tells it to label x axis in km
    
    minZ=0e3;
    maxZ=30e3;
    ncont=40;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
case 48
    fact=1;
    logflag=0;
    clab=0;
    figlab='Max W TimH';
    clines=0;
    clab=1;
    idirstamp=1;
    
    
    iminovr=[0];
    mincovOvr=[0.5];
    
    imaxovr=[1];
    maxcovOvr=[36];
    
    
    z=GridDan(idir).Z;
    time=GridDan(idir).t+3;
    
    %minZ=0e3;
    %maxZ=max(z);
    
    %iylim=1;
    %iylims=[0.62 0.62+GridDan(2).Z(end)/1000];
    
    minZ=0e3;
    %minZ=12e3;
    maxZ=22e3;
    
    nplots2d=1;
    
    tit(1).tit=['Max Updraught (m/s)'];
    
    itimelab=1; 
    
case 488
    fact=1;
    logflag=0;
    clab=0;
    figlab='Max W TimH';
    clines=0;
    clab=1;
    idirstamp=0;
    
    
    iminovr=[0];
    mincovOvr=[0.5];
    
    imaxovr=[1];
    maxcovOvr=[21];
    
    
    z=hemm(:,1)*1000 - 620;
    time=temm(1,:);
    
    %minZ=0e3;
    %maxZ=max(z);
    

    
    minZ=z(1);
    %minZ=12e3;
    maxZ=z(end);
    
    nplots2d=1;
    
    tit(1).tit=['Max Updraught (m/s)'];
    
    itimelab=1; 
    
    i2d=3;
    timesTH(1).t=time;
    
    xlabelstr='Time UTC';
    
    
case 47
    fact=1;
    logflag=0;
    clab=1;
    figlab='temp t=46';
    
    
    iminovr=[1 1];
    mincovOvr=[-80];
    i2d=2;
    
    
    minZ=0e3;
    maxZ=25e3;
    
    nplots2d=1;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    
    tit(1).tit=['Temp (^oC)'];
    
case 46
    fact=1;
    logflag=0;
    clab=1;
    figlab='max sat ice total';
    
    
    %iminovr=[1 1];
    mincovOvr=[0.01];
    i2d=0;
    
    
    minZ=15e3;
    maxZ=22e3;
    
    ncont=40;
    
    tit(1).tit=['Max Total Water Sat Ratio wrt ice'];    
    %tit(1).tit=['Max Vapour in Updraughts (g/kg)'];
    %tit(1).tit=['Max W (m/s)'];
    %tit(1).tit=['Min Temp (^oC)'];
    %tit(2).tit=['Mean Ice Mixing Ratio Time Height for MPC (g/kg)'];
    
case 45
    fact=1;
    logflag=1;
    clab=1;
    figlab='Ice MR';
    ncont=15;
    
    
    iminovr=[1 1];
    mincovOvr=[-7];
    i2d=0;
    
    
    minZ=15e3;
    maxZ=22e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Mean Ice+Snow+Graupel Mixing Ratio Time Height for LEM (g/kg)'];
    tit(2).tit=['Mean Ice Mixing Ratio Time Height for MPC (g/kg)'];
    
    
case 44
    fact=1;
    logflag=0;
    clab=1;
    figlab='Ice MR';
    ncont=15;
    
    manclab=1;
    clines=1; %=1 for black lines
    
    %iminovr=[1 1];
    %mincovOvr=[-7];
    i2d=0;
    
    imaxovr=[1 1];
    maxcovOvr=6; %ppmv for max
    maxcovOvr=2; %ppmv for mean
    
    minZ=15e3;
    maxZ=22e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    %time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Max Ice+Snow+Graupel Mixing Ratio Time Height for LEM (ppmv)'];
    tit(2).tit=['Max Ice Mixing Ratio Time Height for MPC (ppmv)'];
    
    dumprange=[1:83];
    
case 43
    fact=1;
    logflag=1;
    clab=1;
    figlab='Ice MR';
    ncont=15;
    
    
    iminovr=[1 1];
    mincovOvr=[-2];
    i2d=0;
    
    
    minZ=0e3;
    maxZ=22e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    %    time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Mean Ice+Snow+Graupel Number Time Height for LEM (mg^{-1})'];
    tit(2).tit=['Mean Ice Number Time Height for MPC (mg^{-1})'];
    
    dumprange=[1:63];
    
case 42
    fact=1;
    logflag=1;
    clab=1;
    figlab='Ice MR';
    %ncont=15;
    
    iminovr=[1 1];
    mincovOvr=[-2];
    i2d=0;
    
    
    minZ=15e3;
    maxZ=22e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Mean Ice+Snow+Graupel Number Time Height for LEM (mg^{-1})'];
    tit(2).tit=['Mean Ice Number Time Height for MPC (mg^{-1})'];
    
case 41
    fact=1;
    logflag=1;
    clab=0;
    figlab='Ice MR';
    
    
    iminovr=[1 1];
    mincovOvr=[-6];
    i2d=2;
    
    
    minZ=15e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Ice+Snow+Graupel MR at time ' time2d ' for LEM (g/kg)'];
    tit(2).tit=['Ice+Snow+Graupel MR at time ' time2d ' for MPC (g/kg)'];
    
case 40
    fact=1;
    logflag=1;
    clab=0;
    figlab='Ice MR';
    
    
    iminovr=[1 1];
    mincovOvr=[-6];
    
    
    minZ=15e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Ice MR at time ' time2d ' for LEM (g/kg)'];
    tit(2).tit=['Ice MR at time ' time2d ' for MPC (g/kg)'];
    
case 39
    fact=1;
    logflag=1;
    clab=1;
    figlab='Ice NC';
    
    
    iminovr=[1 1];
    mincovOvr=[-1];
    
    imaxovr=[1 1];
    maxcovOvr=2.15;
    
    
    minZ=15e3;
    maxZ=22e3;
    
    nplots2d=2;
    %tt=17; %set time index
    izovr=1; %flag to say are setting own z axis
    %time2d=19.5;
    
    time2d=num2str(time(tt),'%.4f');
    
    tit(1).tit=['Ice NC at time ' time2d ' for LEM'];
    tit(2).tit=['Ice NC at time ' time2d ' for MPC'];
    
    ncont=15; 
    
case 38
    fact=1;
    logflag=1;
    clab=1;
    figlab='Mean Snow MR';
    tit(1).tit='Mean Snow MR';
    tit(2).tit='';
    
    iminovr=1;
    mincovOvr=-1;
    
    ncont=20;
    
case 37
    fact=1;
    logflag=0;
    clab=1;
    figlab='Ratio Snow MR';
    tit(1).tit='Ratio for Snow MR';
    tit(2).tit='';
    MI0=1e-15;
    
    ncont=20;
    
case 36
    fact=1;
    logflag=0;
    clab=1;
    figlab='Ratio Snow NC';
    tit(1).tit='Ratio for Snow NC';
    tit(2).tit='';
    MI0=1e-15;
    
    ncont=20;
    
    %imaxovr=1;
    maxcovOvr=1e-3;
    
case 35
    fact=1;
    logflag=0;
    clab=0;
    figlab='Ratio Ice NC';
    tit(1).tit='Ratio for Ice NC';
    tit(2).tit='';
    MI0=1e-15;
    
    %ncont=20;
    
    dumprange=[1:58];
    
    z=GridDan(1).Z;
    
    
case 34
    fact=1;
    logflag=0;
    clab=1;
    figlab='Ratio Ice MR';
    tit(1).tit='Ratio for Ice MR';
    tit(2).tit='';
    
    %imaxovr=1;
    %maxcovOvr=1e10;
    
    %ncont=20;
    
case 33
    fact=1;
    logflag=1;
    clab=1;
    figlab='PIDEP';
    tit(1).tit=titinput;
    tit(2).tit='';
    
    iminovr=1;
    mincovOvr=0.1;
    
    ncont=20;
    
case 32
    fact=1;
    logflag=1;
    clab=0;
    figlab='Entrainment due to ALu vertical flux';
    tit(1).tit='ALu vertical flux of Ice No. Conc.';
    tit(2).tit='';
    
    iminovr=0;
    mincovOvr=0.1;
    
    ncont=20;
    
case 31
    fact=1;
    logflag=1;
    clab=1;
    figlab='PIPRM / MI0';
    tit(1).tit='PIFRW';
    tit(2).tit='';
    
    iminovr=1;
    mincovOvr=0.1;
    
    ncont=20;
    
    MI0=1e-15;
    
case 30
    fact=1;
    logflag=1;
    clab=1;
    figlab='RSAUT';
    tit(1).tit='RSAUT';
    tit(2).tit='';
    
    iminovr=1;
    mincovOvr=0.1;
    
    ncont=20;
    
case 29
    fact=1;
    logflag=1;
    clab=1;
    figlab='RIACI';
    tit(1).tit='RIACI';
    tit(2).tit='';
    
    iminovr=1;
    mincovOvr=0.1;
    
    ncont=20;
    
case 1
    fact=1;
    %tit(1).tit='Max Ice NC (/mg) - No RHOMOG';
    %tit(2).tit='Max Ice NC (/mg) - Normal';
    tit(1).tit='Max Ice MR (g/kg) - No RHOMOG';
    tit(2).tit='Max Ice MR (g/kg) - Normal';
    logflag=1;
    clab=1;
    figlab='Max Ice MR No. Conc';
    
    iminovr=[1 1];
    mincovOvr=-5;
    
    %imaxovr=[1 1];
    maxcovOvr=0.5;
    
    minZ=16e3;
    
    nplots2d=2;
    
    
case 2
    tit(1).tit='Av Ice Process Rate';
    tit(2).tit='Max Ice Mixing Ratio (g/kg)';
    logflag=0;
    clab=0;
    imaxovr=0;
    maxcovOvr=1e-9;
    iproc=8;
    figlab=['Av Ice Process Rate No: ',num2str(iproc)]; 
    
case 3
    fact=1;
    tit(1).tit='Max Snow Number Concentration (/kg)';
    tit(2).tit='Max Snow Number Concentration (/kg)';
    logflag=0;
    clab=0;
    figlab='Max Snow No. Conc';
    
    minZ=15e3;
    
    
case 4
    fact=1;
    logflag=1;
    clab=1;
    figlab='Non-MPC processes * WQ07';
    tit(1).tit='Non MPC';
    tit(2).tit='Max Upper Tracer Mixing Ratio (g/kg)';
    %ncont=30;
    c
case 5
    fact=1;
    logflag=1;
    clab=1;
    figlab='(PSACI+PGACI+PRACI_G+PRACI_S) * q/v';
    tit(1).tit='(PGACI) * q/v';
    tit(2).tit='Max Upper Tracer Mixing Ratio (g/kg)';
    
    iminovr=1;
    mincovOvr=0.1;
    
    ncont=20;
    
case 6
    %fact=1e6*28.97/18;
    logflag=1;
    tit(1).tit='Low updraught 25th Percentile Ice Saturation Mixing Ratio (ppmv)';
    tit(2).tit='High updraught 25th Percentile Ice Saturation Mixing Ratio (ppmv)';
    
case 7
    %fact=1e6*28.97/18;
    logflag=1;
    tit(1).tit='10th Percentile Ice Saturation Mixing Ratio (ppmv)';
    tit(2).tit=tit(1).tit;
    
case 8
    %fact=1e6*28.97/18;
    logflag=0;
    %tit(1).tit='Potential temp (K)';
    tit(1).tit='Ice Sat MR (ppmv)';
    
    tit(2).tit='Ice NC (mg^{-1})';
    
    nplots2d=2;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    i2d=2; %tells it to label x axis in km
    
    minZ=0e3;
    maxZ=30e3;
    ncont=40;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
    notsame=1;
    
case 9
    i2d=1;
    %fact=1e6*28.97/18;
    logflag=1;
    tit(1).tit='Ice Sat MR (ppmv)';
    tit(2).tit=tit(1).tit;
    figlab='LEM dump 85 ice sat MR';
    
case 10
    logflag=1;
    tit(1).tit='Max Ice No. Conc (#/kg)';
    tit(2).tit=tit(1).tit;
    figlab='max ice NC time-height';
    
case 11
    logflag=1;
    iminovr=1;
    mincovOvr=-3;
    imaxovr=1;
    maxcovOvr=log10(15);
    tit(1).tit='Low Updraught Case Max Ice Mixing Ratio (ppmv vapour equivalent)';
    tit(2).tit='Low Updraught Case Max Ice Mixing Ratio (ppmv vapour equivalent)'
    figlab='max ice MR time-height';
    
case 12
    logflag=1;
    iminovr=1;
    mincovOvr=-8;
    tit(1).tit='Max Snow Mixing Ratio (g/kg)';
    tit(2).tit=tit(1).tit;
    figlab='max snow MR time-height';
    
case 13
    logflag=1;
    iminovr=1;
    mincovOvr=-8;
    tit(1).tit='Max Graupel Mixing Ratio (g/kg)';
    tit(2).tit=tit(1).tit;
    figlab='max graupel MR time-height';
    
case 14
    logflag=0;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-1.5;
    tit(1).tit='Max Vapour Deficit (ppmv)';
    tit(2).tit=tit(1).tit;
    figlab='Max Vapour Deficit time-height';
    
case 15
    logflag=1;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=log10(3.5);
    
    imaxovr=1;
    maxcovOvr=log10(8);
    tit(1).tit='Mean Water Vapour (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
    %tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='Min Vapour time-height';
    
    minZ=14e3;
    maxZ=25e3;
    
case 16
    logflag=0;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-1.5;
    tit(1).tit='Max Vapour Deficit (ppmv)';
    tit(2).tit=tit(1).tit;
    figlab='Max Vapour Deficit time-height';
    %imaxovr=1;
    %maxcovOvr=-6.1
    
case 17
    logflag=1;
    fact=1e6*28.97/18;
    %iminovr(1:2)=0;
    mincovOvr=-1.5;
    tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
    tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='Max Vapour time-height';
    %dumprange=[1:64];
    %nplots2d=1;
    %hrstartles=12.67;
    ncont=20;
    clab=1;
    
case 18
    logflag=1;
    %fact=1e6*28.97/18;
    %     iminovr=1;
    %     mincovOvr=-1.5;
    tit(1).tit='Min Sat Vap MR (ppmv)';
    tit(2).tit=tit(1).tit;
    figlab='Min Sat MR';
    %imaxovr=1;
    %maxcovOvr=-6.1
    
    
case 19
    logflag=0;
    fact=1e6*28.97/18;
    %     iminovr=1;
    %     mincovOvr=-1.5;
    tit(1).tit='Min Sat Vap MR (ppmv)';
    tit(2).tit='Min Vap MR (ppmv)';
    figlab='Min Sat + Vap MR';
    %imaxovr=1;
    %maxcovOvr=-6.1
    
    
case 20
    logflag=1;
    fact=1e6*28.97/18;
    %     iminovr=1;
    %     mincovOvr=-1.5;
    tit(1).tit='Max Ice MR (g/kg)';
    tit(2).tit='Max Snow MR (g/kg)';
    figlab='Min Sat + Vap MR';
    iminovr=1;
    mincovOvr=-6.1
    %nplots2d=1;
    a1=2; %so plots as though are two graphs
    notsame=1;
    
case 21
    logflag=1;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-5;
    tit(1).tit='Max Ice MR (g/kg)';
    tit(2).tit='Max Snow MR (g/kg)';
    figlab='Min Sat + Vap MR';
    %iminovr=1;
    %mincovOvr=-6.1
    %nplots2d=1;
    a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots2d=1;
    normcbar=0;
    
case 22
    logflag=1;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-5;
    tit(1).tit='Max Ice MR (g/kg)';
    tit(2).tit='Max Snow MR (g/kg)';
    figlab='Min Sat + Vap MR';
    %iminovr=1;
    %mincovOvr=-6.1
    %nplots2d=1;
    a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots2d=1;
    normcbar=0;
    
case 23
    logflag=0;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-1e-7;
    tit(1).tit='Mass Flux of Snow (kg/m^2/s)';
    tit(2).tit='Mass Flux of Snow (kg/m^2/s)';
    figlab='Snow Flux';
    imaxovr=1;
    maxcovOvr=0;
    %nplots2d=1;
    %a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots2d=1;
    %normcbar=0;
    
case 24
    logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
    tit(1).tit='Fall Speed Flux of Ice (kg/m^2/s)';
    tit(2).tit='Fall Speed Flux of Ice (kg/m^2/s)';
    figlab='Fall Speed Ice Flux';
    imaxovr=1;
    maxcovOvr=1e-9;
    %nplots2d=1;
    %a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots2d=1;
    %normcbar=0;
    
case 25
    logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
    tit(1).tit='Fall Speed Flux of Ice (kg/m^2/s)';
    tit(2).tit='Fall Speed Flux of Ice (kg/m^2/s)';
    figlab='Fall Speed Ice Flux';
    %imaxovr=1;
    maxcovOvr=30;
    %nplots2d=1;
    a1=2; %so plots as though are two graphs
    %notsame=1;
    
    nplots2d=1;
    %normcbar=0;
    
    dz=z(izmin+1:izmax+1)-z(izmin:izmax);
    
case 26
    logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
    tit(1).tit='Net Flux of All Ice Species (kg/m^2/s)';
    tit(2).tit='Net Flux of All Ice Species (kg/m^2/s)';
    figlab='Net Ice Flux';
    %imaxovr=[0 1];
    maxcovOvr=0.5e-3;
    %nplots2d=1;
    a1=2; %so plots as though are two graphs
    notsame=1;
    
    ncont=16;
    
    %nplots2d=1;
    %normcbar=0;
    
case 27
    logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
    tit(1).tit='Average Sublimation Rate of Ice (ppmv/s)';
    tit(2).tit='Average Sublimation Rate of Ice (ppmv/s)';
    figlab='Average Sublimation Rate of Ice';
    %imaxovr=1;
    maxcovOvr=6.5e-5;
    
case 28
    logflag=1;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-2;
    tit(1).tit='Average Mixing Ratio of All Ice Species (ppmv vapour equivalent)';
    tit(2).tit='Average Mixing Ratio of All Ice Species (ppmv vapour equivalent)';
    figlab='Average Ice MR';
    imaxovr=1;
    maxcovOvr=0.6;
    ncont=17;
    
                               
    
end  %switch(plotcase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of first switch command %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%timesTH=[hrstartles+((dumprange-1)*dumpint)/3600];




if noplot==0
    scrsz=get(0,'ScreenSize');
    %posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

    switch comp
        case 'uni'
            %    posit=[9 60 scrsz(3)/1.2 scrsz(4)/1.14];
            posit=[9 60 scrsz(3)/1.7 scrsz(4)/1.8];
        case 'lacieLap'
            posit=[9 50 scrsz(3)/1.46 scrsz(4)/2.07];
        case 'UWchallenger'
            posit=[9 50 scrsz(3)/1.46 scrsz(4)/1.46];
            posit=[9 50 scrsz(3)/1.46 scrsz(4)/1.26];            
    end


    if subplotting==1
        itemp_bypass=0;
        if itemp_bypass==0
        
        %    posit=[9 50 scrsz(3)/2.4 scrsz(4)/1.5];
        %posit=[9 50 scrsz(3)/2 scrsz(4)/1.1];
        posit = [5 30 1200 590]; 

        if nsub==1; hf=figure('position',posit,'name',figlab,'color','w'); end
        h(iplot).h=subplot(xsub,ysub,nsub);

        %     for ih=1:length(h)
        %         posh=get(h(ih).h,'position');
        %         set(h(ih).h,'position',[posh(1) posh(2) posh(3)-0.2 posh(4)]);
        % 	end
        
        [jsubplot,isubplot]=ind2sub([xsub ysub],nsub);
        yoff=0.12;
	xoff=0.08;
	xoff2=0.2;

	yoff_top=0.05;
        yoff_bot=0.1;
        ygap = 0.07; %y gap between plots
        Ly = 1-yoff_top-yoff_bot - ygap*(ysub-1); %total height available for plots
        dy = Ly/ysub + ygap;


        %posh=get(h(iplot).h,'position');
        posh=[xoff+(jsubplot-1)*1/xsub yoff+(ysub-isubplot)*(1-yoff)/ysub 0.9*(1-xoff-xoff2)/xsub 0.9*(1-yoff)/ysub];
	posh=[xoff+(jsubplot-1)*1/xsub yoff_bot+(ysub-isubplot)*dy 0.9*(1-xoff-xoff2)/xsub Ly/ysub];

        %set(h(iplot).h,'position',[posh(1) posh(2) posh(3)-0.2 posh(4)]);
        set(h(iplot).h,'position',posh);
       
	posh_save{nsub}=posh;
 
        end



    else
        hf=figure('position',posit,'name',figlab,'color','w');
        h(iplot).h=subplot(1,1,1);
    end
end


%set(hf,'renderer','painters');

%nplots2d=length(prof);


if nplots2d==1
    a=a1;
else
    a=a2;
end
b=ceil(min(nplots2d,jmax)/2);


phase=2;
for  i=1:nplots2d
    if subplotting==1
        isub=iplot;
    else
        isub=i;
    end
    try      
        if (exist('time') & i2d~=2 & i2d~=1 & i2d~=3); timesTH(i).t=time(dumprange); end
    catch
        disp(' **** WARNING "if (exist(''time'') & i2d~=2 & i2d~=1 & i2d~=3); timesTH(i).t=time(dumprange)" has failed *****');
    end
    
    %exname=strcat('c:/matlabr12/work/bauru/tracersjan2005/force+3_3th3qv/TracerTimH-',num2str(i));
    %xdat(i).x=time;
    %xdat(i).x=datenum(2004,2,24,hrstartles+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0);
    
    if noplot==0
        scrsz=get(0,'ScreenSize');
        posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
    end
    
    if izovr==0; [izmin izmax]=findheight(z,minZ,maxZ); end
    
    
    %     if length(iz)>=1
    % 		iz=iz(1);
    %     else
    %         iz=length(z);
    %     end
    
    
    %pcolor(9+time./3600,z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp
    
    if subplotting==1 & noplot==0   
%        h(isub).h=subplot(a,b,isub); %here
    end
    
switch plotcase
    case 65
        emmTimH
    case 64
        hms
    case 63
        vap_potemp
    case 62
        lem_min_temp
    case 61
        mpc_min_temp
        
    case 60
        mpc_tot_satmr
    case 59
        z=GridDan(idir).Z+620;   
        %[z0 zend]=findheight(z,15.8e3,17e3);
        zz(i).z=zzeq-620; %620 added later
        izmin=1;
        izmax=length(zzeq);
        
        pdat(1).p=qq;
        
        timesTH(1).t=Nevs;
        
        %     case 59
        %         zz(1).z=GridDan(idir).Z(izmin:izmax);
        %         timesTH(idir).TH=Nevs;
        %         pdat(1).p=
        
    case 58
        dq_dehyd_10thSep2005
        
    case 57
        switch i
            
        case 1
            [izmin izmax]=findheight(GridDan(jc).Z,minZ,maxZ);
            timesTH(1).t=GridDan(jc).Y1'/1000;
            zz(1).z=GridDan(jc).Z(izmin:izmax);
            timesTH(1).t=GridDan(jc).Y1(:)'/1000;
            %pdat(1).p=squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[13]),3)); %height dependent tracer
            switch i57
            case 1  %use jc so that works with allimpPldanSame2_3radar_new.m case 'vap MR & ice NC'
                pdat(1).p=fact*squeeze(sum(TwoDDan(jc).Q(izmin:izmax,:,[1:6]),3));  %total water
                %                pdat(1).p=fact*squeeze(sum(TwoDDan(jc).Q(izmin:izmax,:,[2:6]),3));  %total condensate 
                %pdat(1).p=squeeze(sum(TwoDDan(jc).Q(izmin:izmax,:,[13]),3));  %
            case 2
                pdat(1).p=fact*squeeze(sum(TwoDDan(jc).Q(izmin:izmax,:,[1]),3)); %vapour
            end
            
        case 2
            [izmin izmax]=findheight(GridDan(jc).Z,minZ,maxZ);
            zz(2).z=GridDan(jc).Z(izmin:izmax);
            timesTH(2).t=GridDan(jc).Y1(:)'/1000;
            pdat(2).p=1e-6*squeeze(sum(TwoDDan(jc).Q(izmin:izmax,:,7:9),3));
            %maxcovOvr=100;
        end 
        
        
    case 577
        
        
        
        if izovr==0      
            xinds=[1:150];
            %xinds=[160:230];
            
            %xinds=[130:260];
            
            %xinds=[300:430];
            
            xinds=[1:1000];
            
            ix=findheight(GridDan(idir).Y1,-477.5e3);
            ix2=findheight(GridDan(idir).Y1,-445.8e3);
            
            ix=findheight(GridDan(idir).Y1,-10e3);
            ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+100e3);
            %   ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+150e3);  %for radar plots
            %   ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+600e3);
            ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+50e3);  %for radar plots
            
            ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+150e3);  %for radar plots
            ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+300e3);
            
            
            
            
                ix=1;
            %   ix2=findheight(GridDan(idir).Y1,GridDan(1).Y1(1)+250e3);
            
            
            ix=2;
            %ix2=125;
            % ix2=floor((length(Grid.Y1))/2);
            
            xinds=[ix:ix2];
            
            
            LY1=length(GridDan(idir).Y1);
            
            
            %%%%  NOTE for 2d runs can get the halo output so should only use the data from TwoD.Q(:,2:end-1,iq) - BUT this was for when using
            %%%%       npes=1 using gcom so could be different (doubt it though)
            %       xinds=[length(GridDan(idir).Y1)-ix2:length(GridDan(idir).Y1)-1 ix:ix2 ]; %doing this so that can see vapour from the other side
            %         
            %         ix2=LY1-1;
            %        newpos=35;
            %        xinds=[newpos:LY1-1 2:ix2+newpos-(LY1-ix2) ]; %doing this so that can see vapour from the other side
            
            %  newpos=25;
            %  xinds=[newpos:2*ix2+newpos]; %doing this so that can see vapour from the other side
            
            %(appears in the middle of the domain)
            
%            i3d=1;
            if i3d==1
                X3d=GridDan(1).X1(2:end-1);
                ix2=findheight(X3d,X3d(1)+150e3);
                xinds=[length(X3d)-ix2:length(X3d) ix:ix2 ]; %doing this so that can see vapour from the other side
                xinds=[length(X3d)-length(X3d)/2:length(X3d) 1:length(X3d)/2 - 1]; %doing this so that can see vapour from the other side
                
            end
            
            wrap2d=1;
            if wrap2d==1
                % xmid=75e3/2;
                X3d=GridDan(idir).Y1;
                %ix2=findheight(X3d,X3d(1) + xmid);
                %ix=2;
                % xinds=[length(X3d)-ix2:length(X3d)-1 ix:ix2 ]; %doing this so that can see vapour from the other side

%%%% D is the half the total window size %%%%%                
                %D=62.5e3;
                D=75e3;
         %       D=32.5e3;
%                D=40e3;
                D=15e3;

                %DL=length(X3d); %if want all the domain
                DL=round(2*D/diff(GridDan(idir).Y1(1:2))); %if want 2*D m wide window
               % xinds=[length(X3d)-round(DL/2):length(X3d)-1 2:round(DL/2) - 1]; %doing this so that can see vapour from the other side    
                
                xinds=[round( length(X3d)/2 - DL/2 ):round( length(X3d)/2 + DL/2 ) ]; %doing this so that can see vapour from the other side    
                
            end
            
            
            %     xinds=1:length(GridDan(idir).Y1);
            %    xinds=2:length(GridDan(idir).Y1)-1;
            
            %        xmid=35e3; %x position for where want the new centre
            %        D=74e3;
            
            
            %    xinds=[37:75 1:36];
            
            [izmin izmax]=findheight(GridDan(idir).Z,minZ,maxZ);
            
            if length(izmin)==0 | length(izmax)==0
                disp('*********IZMIN or IZMAX not set properly***********');
            end
            
            %        timesTH(1).t=GridDan(1).Y1'/1000;
            zz(1).z=GridDan(idir).Z(izmin:izmax);
            
            %normal case
            if ~strcmp(i577,'vap_3d_vert')        
                %timesTH(1).t=GridDan(idir).Y1(xinds)'/1000;
                %timesTH(1).t=GridDan(idir).Y1(1:length(xinds))'/1000;
                clear diff
                dy=diff(GridDan(idir).Y1(1:2))/1000;
                L=(length(xinds)-1)*dy;
                timesTH(1).t = [ - L/2 : dy : L/2 ];
                if length(timesTH(1).t<length(xinds))
                    timesTH(1).t(end+1)=timesTH(1).t(end)+dy;
                end
                if length(timesTH(1).t>length(xinds))
                    timesTH(1).t(end)=[];
                end
                
            end
            %pdat(1).p=squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[13]),3)); %height dependent tracer
            %        pdat(1).p=squeeze(sum(TwoDDan(1).W(izmin:izmax,:),3));  %total water
            %                pdat(1).p=fact*squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[1]),3)); %vapour
            
        end %izovr
        

        
%%%%%%%%%%%%%%%%%%%%%        
switch i577
%%%%%%%%%%%%%%%%%%%%%   
    case 'GCM_longitude_transect_UW'
        
        switch lon_csection_case
            case 'CFAD longitude cross section'
%                 b=squeeze(meanNoNan(gcm_calipso_CFAD,1));
%                 sr = repmat(0.5*(CFAD_sr_edges(1:end-1)+CFAD_sr_edges(2:end)),[40 1 48 48]);
%                 sr = permute(sr,[1 3 4 2]);
% 
%                 ilat=13;
%                 lon_inds=[25:40];
%                 mean_sr = sum( sr.*b ,4 ) ./ sum(b,4);
%                 pdat(1).p = squeeze(mean_sr(:,ilat,lon_inds)); pdat(1).p(end+1,:)=NaN; pdat(1).p(:,end+1)=NaN;
% 
%                 lons = gcm_Plon2D_edges(ilat,lon_inds(1)-1:lon_inds(end));

pdat(1).p = cfad_csection'; pdat(1).p(end+1,:)=NaN; pdat(1).p(:,end+1)=NaN;

dlons = diff(cfad_csection_lons);
lons=[cfad_csection_lons(1)-dlons(1) cfad_csection_lons];

                timesTH(1).t=lons;
%                zz(1).z=eval(['csatindx_edges_' gcm_str '/1e3']);
                zz(1).z=eval(['cfad_alts_edges_' gcm_str '/1e3']);   
                
            case 'Model field longitude cross section'  %CPT

                %run make_z_lon_slice first
                
%                dpcolor(x_zslice,z_zslice,dat_zslice); set(gca,'ydir','reverse')
                
pdat(1).p = dat_zslice;
%cfad_csection'cfad_csection'; pdat(1).p(end+1,:)=NaN; pdat(1).p(:,end+1)=NaN;

% dlons = diff(cfad_csection_lons);
% lons=[cfad_csection_lons(1)-dlons(1) cfad_csection_lons];

                timesTH(1).t=x_zslice;
%                zz(1).z=eval(['csatindx_edges_' gcm_str '/1e3']);
                zz(1).z=z_zslice;     
                
                tit(1).tit = titlenam;
                
               
                

        end




    case 'MODIS_plot_UW'
        switch MODIS_plot
            case 'global flat map'

                timesTH(1).t=[1:360];
                zz(1).z=[1:180];

                pdat(1).p = modis_data;

            case '2D plot'
                
                thresh_str='xxx'; %default value
                
                switch modis2d_case
                    case 'general XY 2D PDF'
                        icolmap=1;
                        colormap_jet_no_white;
%                        cmap=jetNEW;
%                        cmap=colormap('bone');
%                        cmap=zeros([64 3]);
cmap=jet;  

savedir='~/modis_work/plots/';


bin2D_method = 'bin_data2D';
%bin2D_method = 'histcn';
bin2D_method = 'ndHistc';
%bin2D_method = '';
                        
%%%% this does the bulk of the data choosing for the plot    
MODIS_data_type = 'L3';
%MODIS_data_type = 'L2';
%MODIS_data_type ='L3 mock'; bin2D_method = ''; %already in the
%Tau_Reff_2D_PDF array

switch MODIS_data_type
    case 'L3'
        time_series_type = 'daily';
        pdf2D_plot_commands

        modis_year_str2='';
        if exist('modis_data_case')
            for ifiles=1:length(modis_data_case)
                modis_year_str2 = [modis_year_str2 char(modis_data_case{ifiles}) ' '];
            end
            modis_year_str2(end)='';
        end


        if length(modis_year_str2)>0
            %                short_plot_name=[modis2d_case ' for day ' modis_day_str ', Y' modis_year_str];
            short_plot_name=['Y' modis_year_str2];
        end


    case 'L2'
        pdf2D_plot_commands_L2   
        
    case 'L3 mock'
        X=0;
        Y=0;
        Xbins=Tau2D_bins;
        Ybins=Re2D_bins;
        
        idat=3;
              
        ilat=ilats(idat);
        ilon=ilons(idat);
        qh = squeeze(Tau_Re_2D_PDF(ilat,ilon,:,:));
        
        qh(:,end+1)=0;
        qh(end+1,:)=0;

        qh=qh';
        
        xlabelstr = 'Tau';
        ylabelstr = 'Reff (\mum)';
        
        short_plot_name=['Tau vs Reff for LAT=' num2str(LATS(ilat)) ', LON=' num2str(LONS(ilon))];
        tit(1).tit=[short_plot_name];
        savename=[savedir tit(1).tit];

end
%%%%
                                                                                                                     


if size(X,1)==1
    X=X';
end
if size(Y,1)==1
    Y=Y';
end

X_orig=X;
Y_orig=Y;

dX=(Xbins(end)-Xbins(1))/nXpdf;
Xbins(1)=Xbins(1)-dX/1e6;
Xbins(end)=Xbins(end)+dX/1e6;

dY=(Ybins(end)-Ybins(1))/nYpdf;
Ybins(1)=Ybins(1)-dY/1e6;
Ybins(end)=Ybins(end)+dY/1e6;

if ~exist('ndims_hist')
    ndims_hist = 2;
end

switch ndims_hist
    case 2
        [X,Y]=bin_data2D_remove_NaNs(X,Y);
    case 3
        [X,Y,Z]=bin_data3D_remove_NaNs(X,Y,Z);

        if size(Z,1)==1
            Z=Z';
        end
        dZ=(Zbins(end)-Zbins(1))/nZpdf;
        Zbins(1)=Zbins(1)-dZ/1e6;
        Zbins(end)=Zbins(end)+dZ/1e6;
end



%the ndHistc function seems to go wrong if any of the data is higher than
%the highest bin value - but makes sense to remove data outside of the bins
%to start with anyway
%   *** now this is done in the wrapper script ndhistc_run ***
% ixrem = find(X>Xbins(end));
% X(ixrem)=[];
% Y(ixrem)=[];
% 
% ixrem = find(X<Xbins(1));
% X(ixrem)=[];
% Y(ixrem)=[];
% 
% iyrem = find(Y>Ybins(end));
% X(iyrem)=[];
% Y(iyrem)=[];
% 
% iyrem = find(Y<Ybins(1));
% X(iyrem)=[];
% Y(iyrem)=[];



%                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);
                        
%qh=histcn([X Y],Xbins,Ybins);  %faster version from the central matab
%exchange

if ~exist('man_choose_plotTimeHeight_graph')
    logflag=1; %when setting logflag=1 make sure that iminovr and mincovOvr is set so that the min value is not zero
    dlogflag=0;
    if logflag==1
        iminovr=1;
        mincovOvr = 1e-4;
%        mincovOvr = 1e-2;
        imaxovr=1;
        maxcovOvr = 1;
    else
        iminovr=0;
        mincovOvr = 0;
        imaxovr=0;
        maxcovOvr = 1;
    end
    
    iarea_normalize=1; %Normalise each 2D bin by the area of X and Y bin
       % Probably important when using unevely spaced bins.    
end

if ~exist('man_choose_plotTimeHeight_graph2')
    iuseYnorm=0;
    iuseXnorm=0;
    iuse_overall_norm=1;
end







switch bin2D_method
    case 'bin_data2D'
        [qh,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);
    case 'histcn'
        [count edges mid loc] = histcn(X, varargin);
    case 'ndHistc'      
        switch ndims_hist
            case 2
                qh = ndhistc_run([X Y], Xbins, Ybins);  %even faster version from the central matab exchange
                %this is the mex version (i.e. C)
    
                pause(1)
                %qh=ones([length(Xbins)-1 length(Ybins)-1]);
                %add dummy values to the end so that pcolor doesn't chop off the last row
                %and column
                qh(:,end+1)=0;
                qh(end+1,:)=0;

                qh=qh';                
                
                %Also do a 1D histogram separately for the x and y bins in order to sample the full range of data and get an accurate mean for each x and y bin.

                [qh_X,bin_no_X] = histc([X], Xbins);  %matlab 1d histo function - returns the bin number that each value in X goes into (ignore last bin)
                clear Y_mean_accurate
                for ibin=1:length(qh_X)-1
                    Y_mean_accurate(ibin) = mean(Y(bin_no_X==ibin)); %mean of Y values for each X bin value.
                    Y_median_accurate(ibin) = prctile(Y(bin_no_X==ibin),50); %median of Y values for each X bin value.
                end
                
                [qh_Y,bin_no_Y] = histc([Y], Ybins);  %matlab 1d histo function - returns the bin number that each value in X goes into (ignore last bin)
                clear X_mean_accurate
                for ibin=1:length(qh_Y)-1
                    X_mean_accurate(ibin) = mean(X(bin_no_Y==ibin)); %mean of X values for each Y bin value.
                    X_median_accurate(ibin) = prctile(X(bin_no_Y==ibin),50); %mean of X values for each Y bin value.
                end
                              

            case 3               
                qh_out = ndhistc_run([X Y Z], Xbins, Ybins, Zbins);
                
                Xmid=0.5*(Xbins(1:end-1)+Xbins(2:end));
                Ymid=0.5*(Ybins(1:end-1)+Ybins(2:end));
                Zmid=0.5*(Zbins(1:end-1)+Zbins(2:end));
                
                Z3D = repmat(Zmid,[size(qh_out,1) 1 size(qh_out,2)]);
                Z3D = permute(Z3D,[1 3 2]);
                        
               
                                                               
                %this is what we want to plot - the mean of the value for
                %each 2D bin
                qh = sum(qh_out.*Z3D,3) ./ sum(qh_out,3);
                
                %calculate the mean in each bin weighted by the frequency
                %in each bin (i.e., contribution of each bin to overall
                %mean)
                freq_2D = sum(qh_out,3);
                qh_mean_contribution = qh.*freq_2D / sum(freq_2D(:));
                
                if ioverall_contribution_mean_2D==1
                    qh = qh_mean_contribution;
                end
                
                %the sum of the squares of the Z - if Z is the absolute
                %error then this is useful for working out the combined
                %error in X, Y, or overall means.
                sumsq2D = sum(qh_out.*Z3D.^2,3);
                N2D = sum(qh_out,3);  
                NX_3D = sum(N2D,2);
                NY_3D = sum(N2D,1);                

                sumsq1DY = sum(sumsq2D,1);
                PDFerror_Y = sqrt(sumsq1DY) ./ NY_3D;

                sumsq1DX = sum(sumsq2D,2);
                PDFerror_X = sqrt(sumsq1DX) ./ NX_3D;
                
                
                %pad the last dimensions with zeros so that all the array
                %is plotted
                qh(:,end+1)=0;
                qh(end+1,:)=0;

                qh=qh';
                
                logflag=0; %when setting logflag=1 make sure that iminovr and mincovOvr is set so thatt the min value is not zero
                dlogflag=0;

                mincovOvr = 0;
                maxcovOvr = 800;

                imaxovr=0;
                iminovr=0;

                if ~exist('ifold_over_ACSIS') 
                	ifold_over_ACSIS=0;
		end
                if ifold_over_ACSIS==1
                	lx = size(qh,1)-1;
			ly = size(qh,2)-1;
			for ix=1:lx
				for iy=1:ix-1
					qh(ix,iy) = qh(ix,iy) + qh(iy,ix); %sum the elements opposite from the diagnonal
					qh(iy,ix) = NaN;
				end
			end    
                end
		clear ifold_over_ACSIS
                

        end
                

       
end

if length(X)==0 | length(Y)==0
   fprintf(1,'\n*** WARNING - no non-Nan data in the fields for the PDF!!\n');
   return 
end
        
disp('*** NOTE - qh has dummy zero values added to the end of each row and column for pcolor plotting purposes. True data lies in qh(1:end-1,1:end-1) (for Xbins and Ybins edges)');

pdf_norm = qh; %default - number of counts


if iarea_normalize==1
    dX_bins = diff(Xbins);
    dY_bins = diff(Ybins);     
    ndX = length(dX_bins);
    ndY = length(dY_bins);    
    dX_bins = repmat(dX_bins,[ndY 1]);   
    dY_bins = repmat(dY_bins,[ndX 1]); 
    
    area_bins = dX_bins.*dY_bins';
    
    pdf_norm(1:end-1,1:end-1) = pdf_norm(1:end-1,1:end-1)./area_bins;
end

if ndims_hist==2
    pdf_norm = pdf_norm/sum(sum(pdf_norm));
end



%values if normalise along the x-axis (i.e. divide each row by the total N
%for that row)
NpX = sum(qh,2);
NpX2D = repmat(NpX,[1 size(qh,2)]);
pdf_normX = qh./NpX2D;

%or normalise along the y-axis
NpY = sum(qh,1);
NpY2D = repmat(NpY,[size(qh,1) 1]);
pdf_normY = qh./NpY2D;



if iuseXnorm==1
    pdat(1).p = pdf_normX;
    short_plot_name=[short_plot_name '-NORMALISED-ALONG-X'];
    tit(1).tit=[tit(1).tit '-NORMALISED-ALONG-X'];
    savename=[savename '-NORMALISED-ALONG-X'];
elseif iuseYnorm==1
    pdat(1).p = pdf_normY;
    short_plot_name=[short_plot_name '-NORMALISED-ALONG-Y'];
    tit(1).tit=[tit(1).tit '-NORMALISED-ALONG-Y'];
    savename=[savename '-NORMALISED-ALONG-Y'];
elseif iuse_overall_norm==1
    pdat(1).p = pdf_norm;
else
    pdat(1).p = qh;
end

if size(Xbins,2)==1
    Xbins = Xbins';
end
if size(Ybins,2)==1
    Ybins = Ybins';
end

mid_Xbins = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );
mid_Ybins = 0.5 * ( Ybins(1:end-1) + Ybins(2:end) );



%calculate the means in 1D for both X and Y as might be useful
X2d=repmat(mid_Xbins,[size(pdf_norm,1)-1 1]);
X_mean = sum(pdf_norm(1:end-1,1:end-1).*X2d,2)./sum(pdf_norm(1:end-1,1:end-1),2);
[maxval,imax]=max(pdf_norm,[],2);
X_mode = Xbins(imax);
NX_vals = sum(qh(1:end-1,1:end-1),2);
NX_vals(isnan(NX_vals))=0;
X_mean2=X_mean; X_mean2(isnan(X_mean))=0;
%X_mean_overall = sum(X_mean2.*NX_vals)./sum(NX_vals);
X_mean_overall = meanNoNan(X(:),1);

X_mean2D = repmat(X_mean2,[1 size(qh,2)-1]);
std_dev_X = sqrt (   sum( qh(1:end-1,1:end-1).*(X2d-X_mean2D).^2 , 2 ) ./ (NX_vals)   );

%calculate the sqrt of the sum of the squares, then divided by N - e.g. for
%combined error of a mean calculation (x here would be the absolute errors
%for the quantity
rms_X = sqrt (   sum( qh(1:end-1,1:end-1).*(X2d).^2 , 2 )  ) ./ NX_vals;
%when have zero counts will be doing 0/0 = NaN

Y2d=repmat(mid_Ybins,[size(pdf_norm,2)-1 1])';
Y_mean = sum(pdf_norm(1:end-1,1:end-1).*Y2d,1)./sum(pdf_norm(1:end-1,1:end-1),1);
[maxval,imax]=max(pdf_norm,[],1);
Y_mode = Ybins(imax);

%Takes a while, so might not want to run this all the time :-
%[Y_median,imed] = percentiles_from_PDF(mid_Ybins,pdf_norm(1:end-1,1:end-1),50,1);

NY_vals = sum(qh(1:end-1,1:end-1),1);
NY_vals(isnan(NY_vals))=0;
Y_mean2=Y_mean; Y_mean2(isnan(Y_mean))=0;
%Y_mean_overall = sum(Y_mean2.*NY_vals)./sum(NY_vals);
Y_mean_overall = meanNoNan(Y(:),1);  %The true mean of all the data (after screening)
%N.B. plot X_mean agains Ybins (since is the mean along the x direction for
%all the different y values)

Y_mean2D = repmat(Y_mean2,[size(qh,1)-1 1]);
std_dev_Y = sqrt (   sum( qh(1:end-1,1:end-1).*(Y2d-Y_mean2D).^2 , 1 ) ./ (NY_vals)   );
rms_Y = sqrt (   sum( qh(1:end-1,1:end-1).*(Y2d).^2 , 1 )  ) ./ (NY_vals);
%when have zero counts will be doing 0/0 = NaN

%Calcualte the correlation between X and Y.
if length(X)>0    
        corr_coeffXY = corr(X(:),Y(:));    
else
    fprintf(1,'\n*** WARNING - X has zero length...\n');
end

tit(1).tit = [tit(1).tit ', corr_coeff=' num2str(corr_coeffXY)];


timesTH(1).t=Xbins;
zz(1).z=Ybins;

if ikeep_X_above_zero==1
    timesTH(1).t(timesTH(1).t<=0) = 1e-6;
end
    
%                        pdat(1).p=nvals_norm; %normalised values

ione_to_one_line=0;

switch iplot_mean_XY
    case 'y'
        ipost_plotTime_commands=1;
        post_plotTime_commands{end+1} = 'plot(mid_Xbins,Y_mean_accurate,''wo'',''markerfacecolor'',''r'')';
    case 'x'
        ipost_plotTime_commands=1
        post_plotTime_commands{end+1} = 'plot(X_mean_accurate,mid_Ybins,''wo'',''markerfacecolor'',''r'')';
    case 'y median'
        ipost_plotTime_commands=1;
        post_plotTime_commands{end+1} = 'plot(mid_Xbins,Y_median_accurate,''wo'',''markerfacecolor'',''r'')';
    case 'x median'
        ipost_plotTime_commands=1
        post_plotTime_commands{end+1} = 'plot(X_median_accurate,mid_Ybins,''wo'',''markerfacecolor'',''r'')';
end


                        
%%%%%%                    
                      case 'vs. LAT'
                        
                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd
                        N_H_calc_histo
                        
                        thresh=0.8;
                        
                        [MLAT2d,MLON2d]=meshgrid(MLON,MLAT); %these are the mid points e.g. MLAT=89.5 means 89-90 I think

                        ihtot = find(abs(MLAT2d)<=30 & cf>=thresh);  thresh_str='LAT LTE 30 & CF GTE ';
                        ihtot = find(abs(Solar_Zenith_Mean.data)>=70 & cf>=thresh);  thresh_str='SZA GTE 70 & CF GTE ';
%                        ihtot = find(cf>=thresh);  thresh_str='CF GTE';                        
                        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
                        
%                        iInf=find(N~=Inf);
%                        X = N(iInf); %MOD06
                        
                        case_vs_lat = 'Solar Zenith Mean';
                        case_vs_lat = 'Sensor Zenith Mean';


                        switch case_vs_lat
                            case 'Solar Zenith Mean'
                                X = Solar_Zenith_Mean.data(ihtot);
                            case 'Sensor Zenith Mean'
                                X = Sensor_Zenith_Mean.data(ihtot);
                        end
                       
                       
                        xlabelstr = case_vs_lat;
                        

                        
                        nXpdf=100;
%                        dX=(maxALL(X)-minALL(X))/nXpdf;
                        max_Xval = maxALL(X);
%                        max_Xval=500;
                        dX=(max_Xval-0)/nXpdf;  

%                        dX=0.05;
                        Xend=(floor(max_Xval/dX) + 1)*dX; %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1                        
%                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[75:0.5:180];


                        Y = MLAT2d(ihtot);
                        ylabelstr = 'Latitude';
                        
                        nYpdf=100;
                        dY=(maxALL(Y)-minALL(Y))/nYpdf;
%                        dY=(maxALL(Y)-0)/nYpdf;    
                        
%                        dY=25;
 
%                        Yend=(floor(maxALL(Y)/dY) + 1)*dY;  %making sure the last bin is beyond the last data point
                        Ybins=[-90:dY:90]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                        
                        
                    case 'std(Nd) vs.'

                        thresh=0.8;

                        thresh2 = 50;

                        %this script does the calculations
                        N_H_calc_std_histo
                        MODIS_N_H_calc
                        total_npix=totN./cf;




                        %                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
                        %                ihtot = [1:prod(size(N_histo_mean))]; %all the data
                        ihtot = find(cf>=thresh);  thresh_str=['CF GTE ']; %only plot for data with more a threshold value (of CF in this case)
                        %ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                        ihtot = find(cf>=thresh & totN>thresh2);  thresh_str=['No. cloudy points GT ' num2str(thresh2) ' & CF GTE '];                                                
%                        ihtot = find(cf>=thresh & total_npix>thresh2);  thresh_str=['Total no. points GT ' num2str(thresh2) ' & CF GTE '];                        
%                        ihtot = find(cf>=thresh2 & total_npix>thresh);  thresh_str='CF GTE ' num2str(thresh2) ' & Total no. points GT ';
                        %ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data

                        
                        X = Nstd(ihtot);
%                        X = Nstd_norm(ihtot);
                        xlabelstr= 'Std. dev of N_d divided by mean N_d';
%                        xlabelstr= 'Std. dev of N_d';
                        
                        X = Wstd(ihtot);
                        X = Wstd_norm(ihtot);
                        xlabelstr= 'Std. dev of W divided by mean W';
%                        xlabelstr= 'Std. dev of W';



                        nXpdf=100;
                        %dX=(maxALL(X)-minALL(X))/nXpdf;
                        dX=(maxALL(X)-0)/nXpdf;

%                        dX=25;
                        Xend=(floor(maxALL(X)/dX) + 1)*dX;  %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
                        %                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        %                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                        Nvs_case='WMOD';
                        Nvs_case='CF';
                        %Nvs_case='Scattering angle';
                        Nvs_case='Nd';

                        plot_type=[plot_type ' ' Nvs_case];

                        switch Nvs_case
                            case 'WMOD'
                                ylabelstr = 'W times five sixths (kg m^{-2})';
                                Y = WMOD(ihtot); %MOD06
                            case 'CF'
                                ylabelstr = 'Cloud Fraction';
                                Y = cf(ihtot); %MOD06 CF
                            case 'Scattering angle'
                                ylabelstr = 'Scattering angle';
                                Y = sangle(ihtot);
                            case 'Nd'
                                ylabelstr = 'N_d (cm^{-3})';
                                Y = N(ihtot);
                        end


                        nYpdf=100;
                        %                        dY=(maxALL(Y)-minALL(Y))/nYpdf;
                        dY=(maxALL(Y)-0)/nYpdf;
                        %                        dY=0.05;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY; %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[75:0.5:180];

                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                        
%                        short_plot_name=[xlabelstr ' for day ' modis_day_str ', Y' modis_year_str];
                        savename=[xlabelstr ' for day ' modis_day_str ', Y' modis_year_str];
                        tit(1).tit=[xlabelstr ' for day ' modis_day_str ', Y' modis_year_str];
                        
                    case 'No. pixels vs CF'
                        N_H_calc_histo
                        MODIS_N_H_calc %runs the routine to calculate Nd

                        ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                        ihtot = find(totN>25);  thresh_str='Optical depth';
                        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
                        
%                        X=totN(ihtot);
                        X=totN(ihtot)./cf(ihtot);
                        xlabelstr = 'Total No. pixels';

                        nXpdf=100;
                        %dX=(maxALL(X)-minALL(X))/nXpdf;
                        dX=(maxALL(X)-0)/nXpdf;

%                        dX=25;
                        Xend=(floor(maxALL(X)/dX) + 1)*dX;  %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
                        %                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        %                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)


                        Y = cf(ihtot);
%                        Y = Cloud_Fraction_Liquid_Pixel_Counts.data(ihtot); 
                       
                        ylabelstr = 'Cloud fraction';

                        nYpdf=100;
                        %                        dY=(maxALL(Y)-minALL(Y))/nYpdf;
                        dY=(maxALL(Y)-0)/nYpdf;
                        %                        dY=0.05;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY; %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[75:0.5:180];

                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;

                        
                    case 'Total No. pixels vs tau'

                        N_H_calc_histo
                        MODIS_N_H_calc %runs the routine to calculate Nd
                        
                        thresh=0.8;

                        ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                        ihtot = find(totN>thresh);  thresh_str='Optical depth';
                        ihtot = find(cf>thresh);  thresh_str='Cloud fration GT ';                        
%                        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data

%                        X=totN(ihtot);
                        X=totN(ihtot)./cf(ihtot);
                        xlabelstr = 'Total No. pixels';

                        nXpdf=100;
                        %dX=(maxALL(X)-minALL(X))/nXpdf;
                        dX=(maxALL(X)-0)/nXpdf;
                        dX=(500-0)/nXpdf;
                        
%                        dX=25;
                        Xend=(floor(maxALL(X)/dX) + 1)*dX;  %making sure the last bin is beyond the last data point
                        Xend=(floor(500/dX) + 1)*dX;  %making sure the last bin is beyond the last data point                        
                        Xbins=[0:dX:Xend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
                        %                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        %                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)


                        Y = tau(ihtot);
                        ylabelstr = 'Optical depth';

                        nYpdf=100;
                        %                        dY=(maxALL(Y)-minALL(Y))/nYpdf;
                        dY=(maxALL(Y)-0)/nYpdf;
                        %                        dY=0.05;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY; %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        %                        Ybins=[75:0.5:180];

                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                     
                        
                        
                    case 'Nd vs tau'
                        
                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd
                        
                        ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                        ihtot = find(totN>25);  thresh_str='Optical depth';                        
%                        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data


%                        iInf=find(N~=Inf);
%                        X = N(iInf); %MOD06

                        X=N(ihtot);
                        
                        nXpdf=100;
                        %dX=(maxALL(X)-minALL(X))/nXpdf;
                        dX=(maxALL(X)-0)/nXpdf;    
                        
                        dX=25;
                        Xend=(floor(maxALL(X)/dX) + 1)*dX;  %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                      
                        Y = tau(ihtot);
                        
                        nYpdf=100;
%                        dY=(maxALL(Y)-minALL(Y))/nYpdf;
                        dY=(maxALL(Y)-0)/nYpdf;                        
%                        dY=0.05;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY; %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1                        
%                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[75:0.5:180];
                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                        
                    case 'Nd vs total no. pixels'
                        
                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd
                        N_H_calc_histo
                        
                        thresh=0.8;
                        
                        [MLAT2d,MLON2d]=meshgrid(MLON,MLAT); %these are the mid points e.g. MLAT=89.5 means 89-90 I think

                        ihtot = find(abs(MLAT2d)<=30 & cf>=thresh);  thresh_str='LAT LTE 30 & CF GTE ';
                        ihtot = find(abs(Solar_Zenith_Mean.data)>=70 & cf>=thresh);  thresh_str='SZA GTE 70 & CF GTE ';
%                        ihtot = find(cf>=thresh);  thresh_str='CF GTE';                        
%                        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
                        
%                        iInf=find(N~=Inf);
%                        X = N(iInf); %MOD06

                       
                     
                        X = totN(ihtot)./cf(ihtot);
                        xlabelstr = 'Total number of pixels';
                        
                        nXpdf=100;
%                        dX=(maxALL(X)-minALL(X))/nXpdf;
                        max_Xval = maxALL(X);
                        max_Xval=500;
                        dX=(max_Xval-0)/nXpdf;  

%                        dX=0.05;
                        Xend=(floor(max_Xval/dX) + 1)*dX; %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1                        
%                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[75:0.5:180];


                        Y = N(ihtot);
                        ylabelstr = 'Number of droplets (cm^{-3})';
                        
                        nYpdf=100;
                        %dX=(maxALL(X)-minALL(X))/nXpdf;
                        dY=(maxALL(Y)-0)/nYpdf;    
                        
                        dY=25;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY;  %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                        
                    case 'Nd vs cloud fraction'
                        
                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd

                        
                        
                        iInf=find(N~=Inf);
                        X = N(iInf); %MOD06
                        
                        nXpdf=100;
                        %dX=(maxALL(X)-minALL(X))/nXpdf;
                        dX=(maxALL(X)-0)/nXpdf;    
                        
                        dX=25;
                        Xend=(floor(maxALL(X)/dX) + 1)*dX;  %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                                               
                        Y = Cloud_Fraction_Liquid.data;
                        
                        nYpdf=100;
                        dY=(maxALL(Y)-minALL(Y))/nYpdf;
                        dY=(maxALL(Y)-0)/nYpdf;                        
%                        dY=0.05;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY; %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1                        
%                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[75:0.5:180];
                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                        
                     case 'Nd vs scattering angle'
                        
                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd

                        iInf=find(N~=Inf);
                        X = N(iInf); %MOD06
                        
                        dX=25;
                        Xend=(floor(maxALL(X)/dX) + 1)*dX;  %making sure the last bin is beyond the last data point
                        Xbins=[0:dX:Xend]; %Xbins(end)=Xbins(end)*(1+1e-10); %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[0:25:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
%                        Xbins=[75:0.5:180]; %Xbins(end)=Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                        
                        
                        sangle = Scattering_Angle_Mean.data;                                                
                        Y = sangle; 
                        
                        dY=25;
                        Yend=(floor(maxALL(Y)/dY) + 1)*dY; %making sure the last bin is beyond the last data point
                        Ybins=[0:dY:Yend]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1                        
%                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
%                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        Ybins=[75:0.5:180];
                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=0;
                        
                    case 'W as funciton of tau&reff'
                        N_H_calc_histo
                        
                      %use the bin edges for the plot  
                        tau = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Ybins;
                        reff =Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Xbins*1e-6; %convert to metres

%add strips to the end of the 2D array so that the last column and row is
%plotted by pcolor
                        W2d(:,end+1)=NaN;
                        W2d(end+1,:)=NaN;
                        
                        timesTH(1).t=tau;
                        zz(1).z=reff*1e6;
                        pdat(1).p=W2d;

                    case 'CF35 vs CF06'                                           
                        Xbins=[0:0.01:1]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        Ybins=[0:0.01:1]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1

                        X=Cloud_Fraction_Day_Mean.data(:); %MOD35
                        Y=Cloud_Fraction_Combined.data(:); %MOD06
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm;

                        ione_to_one_line=1;
                
                    case 'Nd with W calculated vs with MODIS W'                                               
                        
                        set_MODIS_NH_flags=1;
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd
                        
                        
%                        xlabel_name = 'N_d, W calculated';
                        iInf=find(N~=Inf);
                        X = N(iInf); %MOD35
                        Xbins=[0:10:maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        Xbins=[0:1:300 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='MODIS'; %use the MODIS LWP
                        MODIS_N_H_calc %runs the routine to calculate Nd

 %                       ylabel_name = 'N_d, MODIS W';
                        iInf=find(N~=Inf);
                        Y = N(iInf); %MOD06
                        Ybins=[0:10:maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        Ybins=[0:1:400 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=1;
                        
                   
                        
                     
                    case 'Nd from histogram vs using mean tau&reff'
                        set_MODIS_NH_flags=1;
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        N_H_calc_histo %runs the routine to calculate Nd
                        N=N_histo_mean;
                                                
%                        xlabel_name = 'N_d, W calculated';
                        iInf=find(N~=Inf);
                        X = N(iInf); %MOD35
                        Xbins=[0:25:maxALL(X)]; %Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        Xbins=[0:5:800 maxALL(X)]; Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)
                        Xbins=[0:1:100]; %Xbins(end)=+Xbins(end)+1e-10; %adding this to last bin as finds X<Xbins(end)

                        set_MODIS_NH_flags=1; %gets reset every time
                        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                        MODIS_N_H_calc %runs the routine to calculate Nd

%                        ylabel_name = 'N_d, MODIS W';
                        iInf=find(N~=Inf);
                        Y = N(iInf); %MOD06
                        Ybins=[0:25:maxALL(Y)]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1                        
                        Ybins=[0:5:800 maxALL(Y)]; Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        Ybins=[0:1:100]; %Ybins(end)=+Ybins(end)+1e-10; %and want to include X==1 and Y==1
                        
                        [nvals,nvals_norm,nvals_area,nvals_norm_area]=bin_data2D(X,Y,Xbins,Ybins);

                        timesTH(1).t=Xbins;
                        zz(1).z=Ybins;
                        pdat(1).p=nvals_norm; %normalised values

                        ione_to_one_line=1;
                        
                  
                        
                        
                    case 'Zonal max SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        pdat(1).p=squeeze(max(Solar_Zenith_Maximum2.timeseries3,[],2)); %
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
                        short_plot_name=['Maximum Solar Zenith Angle'];
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        cont_val=[65 65];
                        
                    case 'Daily max SZA for one location vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
                        pdat(1).p=squeeze(Solar_Zenith_Maximum.timeseries3(:,1,:)); %
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Daily max Solar Zenith Angle for one longitude'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        cont_val=[65 65];     
                        
              case 'Zonal mean of daily max SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
                        pdat(1).p=squeeze(meanNoNan(Solar_Zenith_Maximum.timeseries3,2)); %
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal mean of daily max Solar Zenith Angle'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        cont_val=[65 65];     
                        
                        
                   case 'Zonal std dev of daily max SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
                        [cont_dat_alt,tmp,pdat(1).p]=meanNoNan(Solar_Zenith_Maximum2.timeseries3,2); %                         
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal std dev of daily max Solar Zenith Angle'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=0;
                        mincovOvr=0;

                        imaxovr=0;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=1;
                        xcont2 = timesTH(1).t;
                        ycont2 = zz(1).z;
                        cont_val=[65 65];     
                        
                        
                   case 'Zonal mean of daily min SZA minus one std dev vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
                        [pdat(1).p,tmp,std_dev]=meanNoNan(Solar_Zenith_Minimum.timeseries3,2); %                        
                        pdat(1).p = pdat(1).p - 1.*std_dev;
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal mean of daily min Solar Zenith Angle minus one std dev'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=0;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];   
                        
                        
                    case 'Zonal mean of daily max SZA minus one std dev vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
                        [pdat(1).p,tmp,std_dev]=meanNoNan(Solar_Zenith_Maximum2.timeseries3,2); %                        
                        pdat(1).p = pdat(1).p - 1.*std_dev;
%                        pdat(1).p = prctile(Solar_Zenith_Maximum2.timeseries,10,2);
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal mean of daily max Solar Zenith Angle minus one std dev'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=0;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];            
                        
                        
                  case '10^{th} percentile of daily max SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
%                        [pdat(1).p,tmp,std_dev]=meanNoNan(Solar_Zenith_Maximum2.timeseries3,2); %                        
%                        pdat(1).p = pdat(1).p - 1.*std_dev;
                        pdat(1).p = squeeze(prctile(Solar_Zenith_Maximum2.timeseries3,10,2));
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['10^{th} percentile (zonal) of daily max Solar Zenith Angle'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];                                    
                        
                  case '10^{th} percentile of daily min SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
%                        [pdat(1).p,tmp,std_dev]=meanNoNan(Solar_Zenith_Maximum2.timeseries3,2); %                        
%                        pdat(1).p = pdat(1).p - 1.*std_dev;
                        pdat(1).p = squeeze(prctile(Solar_Zenith_Minimum.timeseries3,10,2));
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['10^{th} percentile (zonal) of daily min Solar Zenith Angle'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];   
                        
                    case 'Percent of longitudes with minSZA lte 65 and maxSZA gt 65 vs time'
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        max_thresh=75;
                        
                        icount=find(Solar_Zenith_Maximum2.timeseries3>max_thresh & Solar_Zenith_Minimum.timeseries3<=65);
                        a=zeros(size(Solar_Zenith_Maximum2.timeseries3));
                        a(icount)=1;
                        pdat(1).p=squeeze(sum(a,2))/360*100;
                                               
%                        pdat(1).p = squeeze(prctile(Solar_Zenith_Minimum.timeseries3,90,2));
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Percentage of longitudes with daily min SZA<=65 and max SZA>' num2str(max_thresh)];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=0;
                        mincovOvr=0;

                        imaxovr=0;
                        maxcovOvr=100;
                        
                        iadd_cont=0;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];   

                  case '90^{th} percentile of daily min SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
%                        [pdat(1).p,tmp,std_dev]=meanNoNan(Solar_Zenith_Maximum2.timeseries3,2); %                        
%                        pdat(1).p = pdat(1).p - 1.*std_dev;
                        pdat(1).p = squeeze(prctile(Solar_Zenith_Minimum.timeseries3,90,2));
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['90^{th} percentile (zonal) of daily min Solar Zenith Angle'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];   
                        
                        
                  case '90^{th} percentile of daily max SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                                                %previously was plotting the zonal (i.e across all lons) min of the
       %daily max SZA. However, because of the orbits at the end of the day
       %that get put into the next day there is always a small region where
       %there is a lower SZA. This will be accessible to the L3 product,
       %but it will only be over a very small region and so perhaps not
       %much use.
%       pdat(1).p=squeeze(min(Solar_Zenith_Maximum.timeseries3,[],2)); %
             %so instead will plot the max SZA for one location and show how this
             %varies with time. Will give some idea of the range of the lowest SZA
             %available across all lons (excluding the small regions just mentioned).
%                        [pdat(1).p,tmp,std_dev]=meanNoNan(Solar_Zenith_Maximum2.timeseries3,2); %                        
%                        pdat(1).p = pdat(1).p - 1.*std_dev;
                        pdat(1).p = squeeze(prctile(Solar_Zenith_Maximum.timeseries3,90,2));
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['90^{th} percentile (zonal) of daily max Solar Zenith Angle'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  

                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];                           

                        
                        
                   case 'Zonal min SZA vs time'
%                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        pdat(1).p=squeeze(min(Solar_Zenith_Minimum.timeseries3,[],2)); %
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
                        short_plot_name=['Minimum Solar Zenith Angle'];
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  
                        
                        iminovr=1;
                        mincovOvr=0;

                        imaxovr=1;
                        maxcovOvr=81.4;
                        
                        iadd_cont=1;
                        cont_val=[65 65];
                        
                    case 'Percentage Nd difference of allSZA vs lowSZA vs time and lat'
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        prc=90;
                        prc_str=num2str(prc);
                        band_str='21'; band_str2='2.1 mum';
                        band_str='16'; band_str2='1.6 \mum';
                        band_str='37'; band_str2='3.7 \mum';

                        
                        pdat(1).p = eval(['squeeze(prctile(100*(Nd_' band_str '_allSZA-Nd_' band_str '_lowSZA)./Nd_' band_str '_allSZA,prc,2));']);
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=[prc_str '^{th} percentile (zonal) of percentage N_{d} difference for ' band_str2];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  
                        
                        imaxovr=1;
                        iminovr=1;
                        switch prc
                            case 10
                                mincovOvr=-30;
                                maxcovOvr=5;
                            case 50
                                mincovOvr=-10;
                                maxcovOvr=15;
                            otherwise
                                imaxovr=0;
                                iminovr=0;
                        end
                        
                        iadd_cont=0;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];
                        
                    case 'Mean percentage Nd difference of allSZA vs lowSZA vs time and lat'
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
%                        prc=90;
%                        prc_str=num2str(prc);
                        band_str='21'; band_str2='2.1 \mum';
                        band_str='16'; band_str2='1.6 \mum';
%                        band_str='37'; band_str2='3.7 \mum';
                        
                        pdat(1).p = eval(['squeeze(meanNoNan(100*(Nd_' band_str '_allSZA-Nd_' band_str '_lowSZA)./Nd_' band_str '_allSZA,2));']);
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal mean of percentage N_{d} difference for ' band_str2];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  
                        
                        imaxovr=1;
                        iminovr=1;

                        mincovOvr=-5;
                        maxcovOvr=15;
 
                        
                        iadd_cont=0;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];         
                        
                        
                    case 'Percentage Re difference of allSZA vs lowSZA vs time and lat'
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
                        prc=10;
                        prc_str=num2str(prc);
                        band_str='21'; band_str2='2.1 mum';
                        band_str='16'; band_str2='1.6 \mum';
                        band_str='37'; band_str2='3.7 \mum';

                        
                        pdat(1).p = eval(['squeeze(prctile(100*(Re_' band_str '_allSZA-Re_' band_str '_lowSZA)./Re_' band_str '_allSZA,prc,2));']);
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=[prc_str '^{th} percentile (zonal) of percentage r_{e} difference for ' band_str2];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  
                        
                        imaxovr=1;
                        iminovr=1;
                        switch prc
                            case 100
                                mincovOvr=-30;
                                maxcovOvr=5;
                            case 50
                                mincovOvr=-10;
                                maxcovOvr=15;
                            case 90
                                mincovOvr=-5;
                                maxcovOvr=50;
                            otherwise
                                imaxovr=0;
                                iminovr=0;
                        end
                        
                        iadd_cont=0;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];
                        
                    case 'Mean percentage Re difference of allSZA vs lowSZA vs time and lat'
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
%                        prc=90;
%                        prc_str=num2str(prc);
                        band_str='21'; band_str2='2.1 \mum';
                        band_str='16'; band_str2='1.6 \mum';
%                        band_str='37'; band_str2='3.7 \mum';
                        
                        pdat(1).p = eval(['squeeze(meanNoNan(100*(Re_' band_str '_allSZA-Re_' band_str '_lowSZA)./Re_' band_str '_allSZA,2));']);
                        
                        pdat(1).p(abs(pdat(1).p)<1e-5)=NaN;
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal mean of percentage r_{e} difference for ' band_str2];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  
                        
                        imaxovr=1;
                        iminovr=1;

                        mincovOvr=-5;
                        maxcovOvr=5;
 
                        
                        iadd_cont=0;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];      
                        
                        
                        
                    case 'Mean percentage Tau difference of allSZA vs lowSZA vs time and lat'
                        zz(1).z=[MLAT-0.5 MLAT(end)+0.5];
                        timesTH(1).t=[daynum_timeseries3-0.5 daynum_timeseries3(end)+0.5];
                        
%                        prc=90;
%                        prc_str=num2str(prc);
%                        band_str='21'; band_str2='2.1 \mum';
%                        band_str='16'; band_str2='1.6 \mum';
%                        band_str='37'; band_str2='3.7 \mum';
                        
                        pdat(1).p = eval(['squeeze(meanNoNan(100*(Cloud_Optical_Thickness_Liquid_Mean_allSZA.timeseries3 - Cloud_Optical_Thickness_Liquid_Mean_lowSZA.timeseries3)./Cloud_Optical_Thickness_Liquid_Mean_allSZA.timeseries3,2));']);
                        
                        idpcolor=1;
                        
                        xlabelstr = 'Day of year';
                        ylabelstr = 'Latitude';   
                        
%                        short_plot_name=['Zonal minimum of daily max Solar Zenith Angle'];
                        short_plot_name=['Zonal mean of percentage optical depth difference'];                        
                        tit(1).tit = short_plot_name;                        
                        savename = [savename short_plot_name ' time-lat'];  
                        
                        imaxovr=1;
                        iminovr=1;

                        mincovOvr=-5;
                        maxcovOvr=25;
 
                        
                        iadd_cont=0;
                        idiff_cont_dat=0;
%                        xcont2 = timesTH(1).t;
%                        ycont2 = zz(1).z;
                        cont_val=[65 65];                                  
                        
                        
                        
                        
                        
                end %-------------    end modis2d_case  ------------------------------------
                
                if strcmp(thresh_str,'xxx')~=1
                    thresh_str2=[' for ' thresh_str];
                else
                    thresh_str2='';
                end
                
                
                short_plot_name=[short_plot_name thresh_str2];
                savename=[savename thresh_str2];
                tit(1).tit=[tit(1).tit thresh_str2];
                        

        end
        
    case 'Particle size dist vs time'
            
%            if ~exist('icas_count')
                icas_count=1;
%            end

            cut_off_size=1; %size (microns) below which to ignore counts for particle concentration
            CAS_LWC_cut_off_sizes = [0 50];
            
            switch instrument
                case {'CAS','CIP'}
                    air_speed_type='CIP probe';
                    airspeed_constant=60; %not needed if are using the CIP_probe
                case {'CAS MAN','CAS no CIP'}
                    air_speed_type='constant';
                    airspeed_constant=10; %for Manchester 19th June tests
                    airspeed_constant=17; %for Manchester 29th July tests                    
            end
            
            switch instrument
                case {'CAS','CIP','CAS no CIP'}                                                
            %get the sample volume and total concentrations, plus air speed if required.
            [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number...
                ,CAS_total_number_cutoff ...                                    
                ,CIP_total_number,LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                        ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                =cas_sample_volume_and_stats2(dat_flt,CAS_time_all,...
                CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,airspeed_constant);   
            %note air speed is in cm/s           
            

            
            
                case {'CAS MAN'} %is using Manchester data from PACS
                       data_particle = data_CAS_PACS(41:41+29,:)';  
                       time_CAS_PACS=data_CAS_PACS(1,:); 
                       
                        [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number...
                        ,CAS_total_number_cutoff ...    
                        ,CIP_total_number,LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                        ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                        =cas_sample_volume_and_stats2(0,time_CAS_PACS,...
                         bins_PACS',data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                         ,CAS_LWC_cut_off_sizes,airspeed_constant); 
            %note air speed is in cm/s
            %using constant airspeed here because this was a lab study with presumed constant flow

            
            end
            

%            timesTH(1).t = (CAS_time_all - CAS_start_time)/3600;
                    timesTH(1).t = time_plot2D;
            
            switch instrument
                case {'CAS','CAS no CIP'}
                    data_particle = CAS_counts_all;
                    log_bins = log10(CAS_bins);
                    dlogD = repmat(diff(log_bins),[1 size(CAS_counts_all,1)])';
                    orig_bins = CAS_bins;
                    sample_volume = sample_volume_CAS;

                case 'CIP'                    
                    data_particle = CIP_counts_all;
%                    data_particle = data;
                    
                    orig_bins = ([0; CIP_bins(1:end-1)] + CIP_bins) / 2 ; %CIP_bins are the mid-points so calculate boundaries
                    log_bins = log10( orig_bins );
                    dlogD = repmat(diff(log_bins),[1 size(CIP_counts_all,1)])';
                    sample_volume = sample_volume_CIP;   
                    
                case 'CDP'
%                     clear data_particle
%                     data_particle(:,1)=zeros(size(time_plot2D))';
%                     
%                     for ichannel=1:30    
%                         if ichannel<10
%                             channel_str=['0' num2str(ichannel)];
%                         else
%                             channel_str=[num2str(ichannel)];
%                         end
%                         
%                         eval_str = ['data_particle(:,ichannel+1) = nc{''CDP_' channel_str '''}(:);'];
%                         eval(eval_str);
%                     end
%                     
%                     CDP_bins = [2.0000000, 3.0000000, 4.0000000, 5.0000000, 6.0000000,...
%                         7.0000000, 8.0000000, 9.0000000, 10.000000,...
%                         11.000000, 12.000000, 13.000000, 14.000000, 16.000000, 18.000000,...
%                         20.000000, 22.000000, 24.000000, ...
%                         26.000000, 28.000000, 30.000000, 32.000000, 34.000000, 36.000000, ...
%                         38.000000, 40.000000, 42.000000, ...
%                         44.000000, 46.000000, 48.000000, 50.000000];

                    data_particle = CDP_all;

                        
                    log_bins = log10(CDP_bins);
                    dlogD = repmat(diff(log_bins)',[1 size(data_particle,1)])';
                    orig_bins = CDP_bins;
                    sample_volume = ones(size(data_particle));    
                    
                    
              case 'CAS MAN'                  
                                   
                    log_bins = log10(bins_PACS);
                    dlogD = repmat(diff(log_bins)',[1 size(data_particle,1)])';
                    orig_bins = bins_PACS;
                                                                                             
                    sample_volume=sample_volume_CAS;
                    
               case 'Welas'
                   
                    data_particle = 1e-6*welas_dat.conc1'; %*1e-6 to convert to cm^-3

                        
                    log_bins = log10(welas_dat.size(:,1));
                    dlogD = repmat(diff(log_bins)',[size(data_particle,1) 1]);
                    orig_bins = welas_dat.size(:,1);
                    sample_volume = ones(size(data_particle));        
                    
               case 'FSSP'
                   data_particle = data_FSSP(23:42,:)'; %
                        
                    log_bins = log10(bins_FSSP)';
                    dlogD = repmat(diff(log_bins)',[size(data_particle,1) 1]);
                    orig_bins = bins_FSSP;
                    sample_volume = ones(size(data_particle));        
     
                    
            end
            
            
            switch var_plot
                case 'dN/dlogD (cm^{-3} \mum^{-1})'
%                    pdat(1).p = (data_particle(:,1:end-1)./dlogD./sample_volume(:,1:end-1))';
%according to Tom LC and DMT the sizes are the upper limit of the bins - so for first bin we don't know the lower
%bin size - so will ignore the first count
                    pdat(1).p = (data_particle(:,2:end)./dlogD./sample_volume(:,2:end))';
                    pdat(1).p(end+1,:) = pdat(1).p(end,:); %add dummy data (repeated from previous bin)
                    
                    switch y_axis_type
                        case 'log10'
                            %        ylabelstr='log_{10} D (microns)';
                            zz(1).z = log_bins(1:end);
                        case {'','log10_matlab'};
%                            zz(1).z = 10.^log_bins; %back to the normal bins (diameter)
                            zz(1).z = orig_bins;
                    end
                 
                 case 'LWC size dist (g cm^{-3})'
%                    pdat(1).p = (data_particle(:,1:end-1)./dlogD./sample_volume(:,1:end-1))';
%according to Tom LC and DMT the sizes are the upper limit of the bins - so for first bin we don't know the lower
%bin size - so will ignore the first count
                    pdat(1).p = LWC_size_dist'; %LWC for each bin and time. pdat needs to be in array order [size time]
                    pdat(1).p(end+1,:) = pdat(1).p(end,:); %add dummy data (repeated from previous bin)
                    
                    switch y_axis_type
                        case 'log10'
                            %        ylabelstr='log_{10} D (microns)';
                            zz(1).z = log_bins(1:end);
                        case {'','log10_matlab'};
%                            zz(1).z = 10.^log_bins; %back to the normal bins (diameter)
                            zz(1).z = orig_bins;
                    end
                    
                 case 'dM/dlogD (\mug m^{-3} \mum^{-1})'
%                    pdat(1).p = (data_particle(:,1:end-1)./dlogD./sample_volume(:,1:end-1))';
%according to Tom LC and DMT the sizes are the upper limit of the bins - so for first bin we don't know the lower
%bin size - so will ignore the first count
                    mid_points_D = 1e-4*(CDP_bins(1:end-1) + CDP_bins(2:end) )/2; %convert from microns to cm
                    D_CDP = repmat(mid_points_D,[size(data_particle,1) 1]);
                    mass = 2.65*pi.*D_CDP.^3/6 .* data_particle(:,2:end);
%                    pdat(1).p = 1e12*mass'; %convert to microgrammes per m^3
                    pdat(1).p = 1e12*(mass./dlogD)'; %convert to microgrammes per m^3
                    pdat(1).p(end+1,:) = pdat(1).p(end,:); %add dummy data (repeated from previous bin)
                    
                    switch y_axis_type
                        case 'log10'
                            %        ylabelstr='log_{10} D (microns)';
                            zz(1).z = log_bins(1:end);
                        case {'','log10_matlab'};
%                            zz(1).z = 10.^log_bins; %back to the normal bins (diameter)
                            zz(1).z = orig_bins;
                    end   
                    
                    case 'M (\mug m^{-3})'
%                    pdat(1).p = (data_particle(:,1:end-1)./dlogD./sample_volume(:,1:end-1))';
%according to Tom LC and DMT the sizes are the upper limit of the bins - so for first bin we don't know the lower
%bin size - so will ignore the first count
                    mid_points_D = 1e-4*(CDP_bins(1:end-1) + CDP_bins(2:end) )/2; %convert from microns to cm
                    D_CDP = repmat(mid_points_D,[size(data_particle,1) 1]);
                    mass = 2.65*pi.*D_CDP.^3/6 .* data_particle(:,2:end);
                    pdat(1).p = 1e12*mass'; %convert to microgrammes per m^3
%                    pdat(1).p = (mass./dlogD)';
                    pdat(1).p(end+1,:) = pdat(1).p(end,:); %add dummy data (repeated from previous bin)
                    
                    switch y_axis_type
                        case 'log10'
                            %        ylabelstr='log_{10} D (microns)';
                            zz(1).z = log_bins(1:end);
                        case {'','log10_matlab'};
%                            zz(1).z = 10.^log_bins; %back to the normal bins (diameter)
                            zz(1).z = orig_bins;
                    end   
                    
                    
                    
                 case 'Backscatter dN/dlogD (cm^{-3} \mum^{-1})'
%                    pdat(1).p = (data_particle(:,1:end-1)./dlogD./sample_volume(:,1:end-1))';
%according to Tom LC and DMT the sizes are the upper limit of the bins - so for first bin we don't know the lower
%bin size - so will ignore the first count

                    data_particle = CAS_back_all;
                    log_bins = log10(bins_CAS_back{2});
                    dlogD = repmat(diff(log_bins),[1 size(CAS_back_all,1)])';
                    orig_bins = bins_CAS_back{2};
                    sample_volume = sample_volume_CAS;
                    
                    pdat(1).p = (data_particle(:,2:end)./dlogD./sample_volume(:,2:end))';
                    pdat(1).p(end+1,:) = pdat(1).p(end,:); %add dummy data (repeated from previous bin)
                    
                    switch y_axis_type
                        case 'log10'
                            %        ylabelstr='log_{10} D (microns)';
                            zz(1).z = log_bins(1:end);
                        case {'','log10_matlab'};
%                            zz(1).z = 10.^log_bins; %back to the normal bins (diameter)
                            zz(1).z = orig_bins;
                    end    
                    
                case 'N (cm^{-3})'
                    pdat(1).p = (data_particle(:,2:end)./sample_volume(:,2:end))'; %note the transpose ' here
                    pdat(1).p(end+1,:) = pdat(1).p(end,:); %this is dummy data- will not be plotted by
                    %pcolor since it chops off the top of the data and plots the blocks between
                    %the two adjacent axis points in zz(1).z.                 
                    %Really for pcolor plots the data should be one row less than the y axis length
                    %but that's not how Matlab does it!
                    
                    switch y_axis_type
                        case 'log10'
                            %        ylabelstr='log_{10} D (microns)';
                            zz(1).z = log_bins(1:end);
                        case {'','log10_matlab'};
                            zz(1).z = 10.^log_bins; %back to the normal bins (diameter)
                    end
                    
                case 'log10 of particle separation'
                        pdat(1).p = CAS_psep_all';
                        zz(1).z = [1:size(CAS_psep_all,2)];
                    
%                    zz(1).z=([1:length(log_bins)]);
            end
            

            
            %can also calculate the LWC by assuming the mid-point is the bin average diameter and 
            %assuming spherical droplets of liquid water
            
            
        case 'wrf_radar_vert'
            
            zz(1).z = 0:0.1:maxALL(z_slice)/1e3; %set up regular vertical grid
            
            for ixdim=1:size(z_slice,1)
                pdat(1).p(:,ixdim) = interp1 (z_slice(ixdim,:)/1e3,Z_slice(ixdim,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before
            end
            
            clear timesTH         
            timesTH(1).t(1)=0;
            
            if ilon_slice==1
                lats=lat2d.var(:,ilon2);

                for ilats=2:length(lats)
                    timesTH(1).t(ilats)= distlatlon(lats(ilats),lon2d.var(ilats,ilon2),lats(1),lon2d.var(1,ilon2));        
                end
                
            else
                lons=lon2d.var(ilat2,:);
                
                for ilons=2:length(lons)
                    timesTH(1).t(ilons)= distlatlon(lat2d.var(ilat2,ilons),lons(ilons),lat2d.var(ilat2,1),lons(1));        
                end
            end
            
        %    itime=idir;   
        
        pdat(1).p=pdat(1).p*f; %convert to ppmv if total water
        
        for ipot=1:size(pot_slice,1)
            i380=min( find(pot_slice(ipot,:)>380) );
            z380(ipot)=z_slice(ipot,i380)/1e3;  %height of the 380K contour
        end

%%  Set the position of the cross section line
%*******************************************************
case 'wrf_vert_cross_section'
 
	    x_grid = ([1:size(lat2d(1).var,2)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));     
	    y_grid = ([1:size(lat2d(1).var,1)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));   
        
        DY=nc.DY(1)/1000;
        DX=nc.DX(1)/1000;
        
        %new x_grid formulation - might cause problems?
        x_grid = ([1:size(lat2d(1).var,2)]-1) * DX;
        y_grid = ([1:size(lat2d(1).var,1)]-1) * DY;
        
                    
%        zz(1).z = [0:0.05:4];  %maxALL(Z.var)/1e3; %set up regular vertical grid
        zz(1).z = [0:0.05:1];  %maxALL(Z.var)/1e3; %set up regular vertical grid
        
        zz(1).z = [0:0.01:0.05 0.1:0.1:12];  %maxALL(Z.var)/1e3; %set up regular vertical grid
        
        zz(1).z = [0:0.01:0.05 0.1:0.1:8];  %maxALL(Z.var)/1e3; %set up regular vertical grid
        
        zz(1).z = [0:0.01:0.05 0.1:0.1:5];  %maxALL(Z.var)/1e3; %set up regular vertical grid
        
%        zz(1).z = [0:0.01:0.05 0.1:0.1:2];  %maxALL(Z.var)/1e3; %set up regular vertical grid
        
        switch hor_vert
            case 0  %cross sections along any line selected
                %	x_line = [350 400];
                %	y_line = [550 365];

                %    x_line = [592 625 ]; %
                %    y_line = [222 120 ];

                x_line = [570 690]; %the x and y co-ords of the line (left side of slice then right)
                y_line = [275 10];

                x_line = [400 400]; %the x and y co-ords of the line (left side of slice then right)
                y_line = [50 3000];

                x_line = [200 1200]; %the x and y co-ords of the line (left side of slice then right)
                y_line = [150 1900];

                x_line = [277.1 544.4]; %the x and y co-ords of the line (left side of slice then right)
                y_line = [379.4 214.3];

                x_line = [302.6745 662.3461];
                y_line = [426.6242 192.6690];

                x_line=[250 700];  %horizontal line that follows the Fohn flow on east side at 9UTC on 6th Jan.
                y_line=[320 320];

                %    x_line=[250 700];
                %    y_line=[355 355];

                %    x_line=[250 800];
                %    y_line=[355 255];

                    x_line=[250 800];
                    y_line=[375 235];  %latest line as of pre 31st July, 2009

                x_line=[10 800];
                y_line=[435 235]; %latest line as of pre 31st July, 2009 but bit longer *** main one used in write-up as of 26th Oct, 2009 ***
                

                
%                x_line=[10 800];
%                y_line=[434.3978 264.3978];                
%                y_line =[423.3063 273.3063];
%                y_line =[412.2148 282.2148];
%                y_line =[401.1232 291.1232];
%                y_line =[390.0317 300.0317];

%                y_line = [478.7639 228.7639];

%                y_line =[506.4928 206.4928]; %rotated clockwise                                
                
               
%                y_line = y_line+10;       
%                y_line = y_line+5;     %when applied to the 31st July line this gives the diving streamline height of 3.1 km (292 K)     
%                y_line = y_line-10;
%                y_line = y_line-4;


                x_line=[10 800];
                y_line=[510 170]; %modified from 31st July (moved north a little to get the 3.1 km streamline height)
                
                x_line=[10 800];
                y_line=[525 185]; %modified from 31st July (moved north a little to get the 3.1 km streamline height) + 15

                x_line=[235.5 800];
                y_line=[428 185]; %modified from 31st July (moved north a little to get the 3.1 km streamline height) + 15 
                %and shortened it to only include ~200 km upwind of the mountain - code 08-12-09_15_x235-5_800_y428_185
                
                x_line=[235.5 800];
                y_line=[413 170]; %modified from 31st July (moved north a little to get the 3.1 km streamline height) 
                %and shortened it to only include ~200 km upwind of the mountain - code 08_12_09_16_x235_5_800_y413_170
                
                %need to convert from old use of DX and DY
                DX_old = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));     
                DY_old = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));
                
        
        
                x_line=[235.5 602.9]*DX/DX_old;
                y_line=[413 254.8]*DY/DY_old; %directly through the jet and at the right angle for the 2-3 km winds. Also shortened
                                  %on both sides to represent only the bits of line parallel to those winds 
                                  %- code 09_12_09_12_x235_5_602_9_y413_254_8
                                  % *** this is the one used in the report
                                  % as of Oct 2013 ***
                                  
                x_line=[224.6659 600]; %as above, but extended to cover the flight track
                y_line=[393.8815 232.3137]; %N.B. - have applied the DX/DX_old, etc. factors here.
                
                x_line = [3400 4200];  %for d01 at 9 UTC 6th Jan, through upper jet
                y_line = [3850 3500];  %For response to Editor's comments RE effect of resolution, Orr (2008), etc.

                x_line = [3400 4200];  %for d01 at 9 UTC 6th Jan, through lower jet
                y_line = [3750 3400];  %For response to Editor's comments RE effect of resolution, Orr (2008), etc.
                
%                x_line=[224.6659 600]; % Oct 2013 - to go through jet 1 (furthest north one)
%                y_line=[393.8815+25+22+12 232.3137+25+22-12]; %
                
%                x_line=[535 650]; % May 2014 - in response to ACPD ref#1 - perpendicular to jet
%                y_line=[347 265]; % at AWS location for 15 UTC to look at 10m and 300m winds                
                                  
%                 x_line=[235.5 602.9];                 
%                 y_line=[458.0000 299.8000];
                
                                  
%               x_line=[235.5 602.9];
%                y_line=[453 304.8]; %As for 09_12_09_12_x235_5_602_9_y413_254_8 (same angle)
                %   but for the northerlymost jet. Code 14_12_09_13_x235_5_602_9_453_304_8
                
%                x_line=[255.5 622.9]; %Same angle as for 09_12_09_12_x235_5_602_9_y413_254_8 
                %   but for the southerlymost jet. Code 16_12_09_15_x255_5_622_9_317_168_8
%                y_line=[317.0 168.8];
                
%                  x_line=[255.5 622.9]; %Same as above but 5 km down. Code 16_12_09_16_x255_5_622_9_312_163_8
%                 y_line=[312.0 163.8];
%                 
%                 x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                 y_line=[382.0 233.8]; 
%                 % Code 12_01_10_12_x255_5_622_9_382_233_8
%                 
%                 x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                 y_line=[362.0 213.8]; 
                % Code 12_01_10_12_x255_5_622_9_362_213_8
                
%                x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                y_line=[342.0 193.8]; 
                % Code 12_01_10_12_x255_5_622_9_342_193_8, lat ~ -68.3664
                
%                 x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                 y_line=[332.0 183.8]; 
%                 % Code 12_01_10_12_x255_5_622_9_332_183_8, lat ~-68.4132
%                 
%                 x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                 y_line=[330.0 181.8]; 
%                 % Code 12_01_10_12_x255_5_622_9_330_181_8, lat ~-68.4772 
%                 
%                 x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                 y_line=[327.0 178.8]; 
%                 % Code 12_01_10_12_x255_5_622_9_327_178_8, lat ~-68.5070                                               
%                 
%                 x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                 y_line=[322.0 173.8]; 
%                 % Code 12_01_10_12_x255_5_622_9_322_173_8, lat ~-68.5424 
                
%                x_line=[255.5 622.9]; %Trying to get a non-jet area 
%                y_line=[317.0 168.8]; 
                % Code 12_01_10_12_x255_5_622_9_317_168_8, lat ~-68.5813
                

                                  
                
                                  
%                y_line = [518.2216 168.2216];
                
%                y_line = [487.4928 187.4928]; %28th Oct, 2009 - one that is a bit steeper than 31st July line but goes 
%                               through the windstorm at 6UTC


                
                %    x_line=[350 500];  %horizontal line that follows the Fohn flow on east side at 9UTC on 6th Jan.
                %    y_line=[320 320];

                %    x_line=[250 600];
                %    y_line=[380 150];
                
%                x_line=[200 350]; 
%                y_line=[250 400];      %across the bay parallel to mountain

                
if ilon_line==1                
                 %LON=-66; %do cross section along this longitude line (set earlier)  
                 [ilat,ilon] = getind_latlon_quick(lat2d_fine,lon2d_fine,-69,LON,0.1);
                 x1=x_fine(ilon); y1=y_fine(ilat);
                 [ilat2,ilon2] = getind_latlon_quick(lat2d_fine,lon2d_fine,-67,LON,0.1);
                 x2=x_fine(ilon2); y2=y_fine(ilat2);                 
                
                 x_line = [x1 x2];
                 y_line = [y1 y2];

             %now get the height of the ridge (max along each constant latitude line)
                hgt=nc{'HGT'}(1,:);
                np2=80; %number of points along the line of longitude to make averages for
                np2=320;
                %         latA=-70;
                %         latB=-64.8;
                latA=-69;
                latB=-67;
                LATS=latA:(latB-latA)/np2:latB;
                
                for imean=1:length(LATS)
                    icons_inds = get_inds_constant_lat(LATS(imean),lat2d,lon2d); %get indices for a constant latitude slice
                    ieast=find(lon2d.var(icons_inds)>-67.5); %to remove the mountains to the west of peninsula from being included
                    [peak_height ipeak]=max(hgt(icons_inds(ieast))); %find the position of the peninsula mountain
                    lon_peak(imean)=lon2d.var(icons_inds(ieast(ipeak)));
                    [ilat_peak(imean),ilon_peak(imean)] = ind2sub(size(lat2d.var),icons_inds(ieast(ipeak))); %get lat and lon indices for peak
                    peak_vs_lat(imean)=peak_height;

                end
                
                
                 
                 iadd_overlay=iadd_overlay+1;
                 dat_overlay(iadd_overlay).x=LATS;
                 dat_overlay(iadd_overlay).y=peak_vs_lat/1000;
                 dat_overlay(iadd_overlay).overlay_style=['''w--'',''linewidth'',4'];
            

end
                
                if x_line(1)>x_line(2)
                    fprintf(1,'ERROR - need to put x_line values in ascending order - changing...');
                    xtemp=x_line(2);
                    x_line(2)=x_line(1);                    
                    x_line(1)=xtemp;
                    
                    ytemp=y_line(2);
                    y_line(2)=y_line(1);
                    y_line(1)=ytemp;
                end


                n_line = [1:50]; %number of points that want along the cross-section line
                n_line = [1:ncross_points]; %number of points that want along the cross-section line
                

                L_line = sqrt( (x_line(2) - x_line(1))^2 + (y_line(2) - y_line(1))^2 );
                d_line = [0:L_line/(length(n_line)-1):L_line];

                angle = atan( (y_line(2) - y_line(1)) / (x_line(2) - x_line(1)) );
                xd = x_line(1) + d_line*(cos(angle));
                yd = y_line(1) + d_line*(sin(angle)); % x and y positions of the equally spaced points along the cross-section line

                if xd(end)>xd(1)
                    dx=1;
                else
                    dx=-1;
                end
                if yd(end)>yd(1)
                    dy=1;
                else
                    dy=-1;
                end

                x_inds = [ findheight_nearest(x_grid,xd(1)) - dx : dx : findheight_nearest(x_grid,xd(end)) + dx ];
                y_inds = [ findheight_nearest(y_grid,yd(1)) - dy : dy : findheight_nearest(y_grid,yd(end)) + dy ];

                x_inds = sort(x_inds);
                y_inds = sort(y_inds);

                grid_x = x_grid(x_inds);
                grid_y = y_grid(y_inds); %smaller grid over area concerned

                timesTH(1).t = d_line;
                
                switch var_plot
                    case {'Wind speed perpendicular component (m s^{-1})','Richardson Number (perpendicular wind)'}                        
                        angle = angle-pi/2; %convert to be perpendicular to the line
                end
                
                angle = angle * ones(size(n_line)); %make a vector
%                angle = angle * ones(size(zz(1).z));


                %length of cross section
            case 1
                n_line = [findheight(x_grid,275):findheight(x_grid,625)];
                %n_line = [findheight(x_grid,250):findheight(x_grid,675)];  %for y=295 km

                angle = 0 * ones(size(n_line)); %make zeros to ignore the new wind component metric



            case 2
                n_line = [findheight(x_grid,275):findheight(x_grid,625)];  %previous sections at e.g. y=380 km
                n_line = [findheight(x_grid,250):findheight(x_grid,675)];  %for y=295 km
                timesTH(1).t = x_grid(n_line);

                angle = 0 * ones(size(n_line)); %make zeros to ignore the new wind component metric

            case 3  %%%%%  cross sections along a streamline

                n_line = [1:ncross_points]; %number of points that want along the cross-section line - note might be changed slightly later
                
                

                ih_wrf=22;
                %    ih_wrf=32;
                i_limit_stream=1; %flag to cut short the cross section at d_max km

                d_max=5500;
                d_max=300;

                u_wind_dat=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
                v_wind_dat=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
                [X,Y] = MESHGRID(x_grid,y_grid);

                %    STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,0,500);  %get streamline through the wind field
                %    STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,1500,1600);  %get streamline through the wind field
%                STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,1100,900);  %get streamline through the wind field
%                STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,250,378);  %get streamline through the wind field
                %STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,1,417);  %get streamline through the wind field starting at x,y
             %% if want to use a streamline that has already been computed by for e.g. streamlines_threeD_draw.m :-
%                STREAM{1}=STREAM_z0{str_succ(1)}/1e3; %convert from m to km for this plotting routine
%                STREAM{1}=STREAM_z0{str_succ(3)}/1e3; %convert from m to km for this plotting routine                
                    %111 for z0=3100m for streamline 2
%                STREAM{1}=STREAM_z0{str_succ(1)}/1e3; %convert from m to km for this plotting routine                

                STREAM{1}=STREAM_z0{171}/1e3; %convert from m to km for this plotting routine                    
                %171 for jet streamline at z0=1200m for jet
                
                STREAM{1}=STREAM_z0{150}/1e3; %convert from m to km for this plotting routine                    
                %150 and 151 for jet streamline at z0=1250m for non-jet
                STREAM{1}=STREAM_z0{134}/1e3; %convert from m to km for this plotting routine                    
                %129,133,134,135 for jet streamline at z0=650m for near jet - ones that go beyond the peak
                STREAM{1}=STREAM_z0{118}/1e3; %convert from m to km for this plotting routine                    
                %118 for streamline at z0=550m for near jet - first height ones that go beyond the peak in that region
                STREAM2{1}=STREAM_z0{127}/1e3; %convert from m to km for this plotting routine                    
                %118 for streamline at z0=600m for near jet - first height ones that go beyond the peak in that region
                %using one cell resolution for streamline calculation
                STREAM2{1}=STREAM_z0{152}/1e3; %convert from m to km for this plotting routine                    
                %152 for streamline at z0=800m for near jet - first height ones that go beyond the peak in that region
                %using one cell resolution for streamline calculation
                
                STREAM2{1}=STREAM_z0{129}/1e3; %z0=650 for write up as of 14th June 2010
                
%                STREAM2{1}=STREAM_z0{148}/1e3; %z0=1050 for write up as of 14th June, 2010
                
%                STREAM2{1}=STREAM_z0{180}/1e3; %z0=1050 for write up as of 14th June, 2010. For 68.1 lat jet
                
%                STREAM2{1}=STREAM_z0{200}/1e3; %z0=1050 for write up as of 14th June, 2010. For 68.6 lat jet
%                STREAM2{1}=STREAM_z0{169}/1e3; %z0=1200 for write up as of 14th June, 2010. For 68.6 lat jet
                %*** above one is only for d_perp = 400 ****
%                STREAM2{1}=STREAM_z0{170}/1e3; %z0=950 for write up as of 14th June, 2010. For 68.6 lat jet
                %*** above one is only for d_perp = 400 ****
%                STREAM2{1}=STREAM_z0{145}/1e3; %z0=1200 for write up as of 14th June, 2010. For 68.6 lat jet
                %*** above one is only for d_perp = 400 ****
                
%                STREAM2{1}=STREAM_z0{129}/1e3; %z0=950 for write up as of 14th June, 2010. For 67.5 lat jet
                %*** above one is only for d_perp = 400 ****
                
%                STREAM2{1}=STREAM_z0{149}/1e3; %z0=950 for write up as of 14th June, 2010. For 68.1 lat jet
                %*** above one is only for d_perp = 400 ****
                
%N.B. the str_succ numbers above will only apply if use the same
%nstream=200 value in streamlines_threeD_draw.m
                
                inan=isnan(STREAM2{1}(:,1));
                inan=find(inan==0);
                STREAM{1}=STREAM2{1}(1:inan(end),:);

                %here we allow to shorten the line if don't want it all - runs to d_max in terms of distance along it
                if i_limit_stream==1
                    xd2 = STREAM{1}(:,1);
                    yd2 = STREAM{1}(:,2);
                    d_line2 = cumsum([0; sqrt(diff(xd2).^2 + diff(yd2).^2)]); %calculate the distance along the line by assuming straight lines between
                    i_limit = findheight(d_line2,d_max); %find the indices of the end of the line required
                else
                    i_limit = length(STREAM{1}(:,1));
                end




                stream_inds = [1:max([1 round(i_limit/length(n_line))]):i_limit]; %equally spaced indices along the streamline
                n_line=[1:length(stream_inds)];
                xd = STREAM{1}(stream_inds,1);
                yd = STREAM{1}(stream_inds,2);
                if size(STREAM{1},2)==3
                    zd = STREAM{1}(stream_inds,3);
                end
                d_line = cumsum([0; sqrt(diff(xd).^2 + diff(yd).^2)]); %calculate the distance along the line by assuming straight lines between


                angle = atan( (diff(yd)) ./ (diff(xd)) );
                angle = [angle(1); angle]; %but this contains discontinuities that will affect
                %the wind speed at those points
                angle=mean(angle)* ones(size(n_line)); %might be better to just take the mean and use the same for all?

                
%                angle = 0 * ones(size(n_line)); %make zeros to ignore the new wind component metric
                
                %     if xd(end)>xd(1)
                %         dx=1;
                %     else
                %         dx=-1;
                %     end
                %     if yd(end)>yd(1)
                %         dy=1;
                %     else
                %         dy=-1;
                %     end

                dx=1;
                dy=1;

                x_inds = [ findheight(x_grid,min(xd)) : dx : findheight(x_grid,max(xd))+dx ]; %find indices for a smaller grid covering the streamline
                y_inds = [ findheight(y_grid,min(yd)) : dy : findheight(y_grid,max(yd))+dy ];

                x_inds = sort(x_inds);
                y_inds = sort(y_inds);

                grid_x = x_grid(x_inds);
                grid_y = y_grid(y_inds); %smaller grid over area concerned

                timesTH(1).t = d_line';


                %length of cross section
                
            case 4  %%%%%  cross section along the mountain crest (at max mountain height)

                

                %                ih_wrf=15;
                %    ih_wrf=32;
                i_limit_stream=1; %flag to cut short the cross section at d_max km

                d_max=5500;

                %                u_wind_dat=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
                %                v_wind_dat=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
                [X,Y] = MESHGRID(x_grid,y_grid);

                %%%%%%%%%%%% find the crest of the mountain
                hgt=nc{'HGT'}(1,:);

                np=80; %number of points along the line of longitude to make averages for
                np=320;
                %         latA=-70;
                %         latB=-64.8;
                latA=-69;
                latB=-67.1;
                LATS=latA:(latB-latA)/np:latB;
                


                for imean=1:length(LATS)
                    icons_inds = get_inds_constant_lat(LATS(imean),lat2d,lon2d); %get indices for a constant latitude slice
                    ieast=find(lon2d.var(icons_inds)>-67.5); %to remove the mountains to the west of peninsula from being included
                    [peak_height ipeak]=max(hgt(icons_inds(ieast))); %find the position of the peninsula mountain
                    lon_peak(imean)=lon2d.var(icons_inds(ieast(ipeak)));
                    [ilat_peak(imean),ilon_peak(imean)] = ind2sub(size(lat2d.var),icons_inds(ieast(ipeak))); %get lat and lon indices for peak
                    peak_vs_lat(imean)=peak_height;

                end
            
                
                
                
%%%%%%%%%                

                %    STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,0,500);  %get streamline through the wind field
                %    STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,1500,1600);  %get streamline through the wind field
%                STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,1100,900);  %get streamline through the wind field
%                STREAM=stream2(X,Y,u_wind_dat,v_wind_dat,250,378);  %get streamline through the wind field
%                n_line = [1:50]; %number of points that want along the cross-section line

                %here we allow to shorten the line if don't want it all - runs to d_max in terms of distance along it
%                 if i_limit_stream==1
%                     xd2 = STREAM{1}(:,1);
%                     yd2 = STREAM{1}(:,2);
%                     d_line2 = cumsum([0; sqrt(diff(xd2).^2 + diff(yd2).^2)]); %calculate the distance along the line by assuming straight lines between
%                     i_limit = findheight(d_line2,d_max); %find the indices of the end of the line required
%                 else
%                     i_limit = length(STREAM{1}(:,1));
%                 end




%                stream_inds = [1:round(i_limit/length(n_line)):i_limit]; %equally spaced indices along the streamline
                stream_inds = [1:length(LATS)];
%                n_line=[1:length(stream_inds)];
                n_line=stream_inds;
%                xd = STREAM{1}(stream_inds,1);
%                yd = STREAM{1}(stream_inds,2);
                xd = x_grid(ilon_peak)';
                yd = y_grid(ilat_peak)';
                
                angle = 0 * ones(size(n_line)); %make zeros to ignore the new wind component metric
                    
                d_line = cumsum([0; sqrt(diff(xd).^2 + diff(yd).^2)]); %calculate the distance along the line by assuming straight lines between


                %     if xd(end)>xd(1)
                %         dx=1;
                %     else
                %         dx=-1;
                %     end
                %     if yd(end)>yd(1)
                %         dy=1;
                %     else
                %         dy=-1;
                %     end

                dx=1;
                dy=1;

                x_inds = [ findheight(x_grid,min(xd)) : dx : findheight(x_grid,max(xd))+dx ]; %find indices for a smaller grid covering the streamline
                y_inds = [ findheight(y_grid,min(yd)) : dy : findheight(y_grid,max(yd))+dy ];

                x_inds = sort(x_inds);
                y_inds = sort(y_inds);

                grid_x = x_grid(x_inds);
                grid_y = y_grid(y_inds); %smaller grid over area concerned

                switch x_axis_type
                    case 'dist'
                        timesTH(1).t = d_line';
                    case 'lat';
                        timesTH(1).t = LATS;    
        
                end

                %length of cross section

        end

%    u = 0.5* (nc{'U'}(time,:,:,1:end-1) + nc{'U'}(time,:,:,2:end) ); %2d wind at one height                
%	 w = 0.5 * ( nc{'W'}(time,2:end,:,:) + nc{'W'}(time,1:end-1,:,:) );
%	Z = squeeze(WRFUserARW(nc,'Z',time));

                       






if recalc==1
    clear wrf_dat u_slice v_slice w_slice z_slice u_quiver v_quiver cont_slice terr_slice lon_slice lat_slice      

    switch hor_vert
        case {0,3,4}   %will have to take components of u and v along the direction of the cross section 
            
            u = nc{'U'}(time,:,y_inds,x_inds);
            v = nc{'V'}(time,:,y_inds,x_inds);
            w = nc{'W'}(time,:,y_inds,x_inds);
            z_temp=(nc{'PH'}(time,:,y_inds,x_inds) + nc{'PHB'}(time,:,y_inds,x_inds) )./9.81;
            z = 0.5.*(z_temp(1:end-1,:,:)+z_temp(2:end,:,:));
            terr_all = nc{'HGT'}(time,:);
            terr = terr_all(y_inds,x_inds);
            
            lat_restrict = lat2d.var(y_inds,x_inds);
            lon_restrict = lon2d.var(y_inds,x_inds);
            
%                     v_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(v(iz,:,:)),yd(i_line),xd(i_line));
%                     w_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(w(iz,:,:)),yd(i_line),xd(i_line));
%                     z_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(Z.var(iz,:,:)),yd(i_line),xd(i_line));
%                     terr_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(nc{'HGT'}(time,:,:)),yd(i_line),xd(i_line));



            switch var_plot
                case 'Relative humidity (%)'
                    potemp = nc{'T'}(time,:,y_inds,x_inds) + 300;
                    P = nc{'P'}(time,:,y_inds,x_inds) + nc{'PB'}(time,:,y_inds,x_inds);
                    T = potemp ./ ( (1e5./P).^0.286 );
                    qv = nc{'QVAPOR'}(time,:,y_inds,x_inds);
                    qs = SatVapPress(T,'goff','liq',P,1)/f;
                    cont_data = 100*( qv./qs );

                case 'Potential temperature (K)'
                    cont_data = nc{'T'}(time,:,y_inds,x_inds) + 300;    
                case 'Equivalent potential temperature (K)'
                    switch file_type
                        case 'met_em'
                            potemp = nc{'T'}(time,:,y_inds,x_inds) + 300;
                            P=eta_get_p(nc);
                            P = P(:,y_inds,x_inds);
                            T = potemp ./ ( (1e5./P).^0.286 );
                            qv = nc{'QVAPOR'}(time,:,y_inds,x_inds);
                            cont_data = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 );
                        otherwise
                            potemp = nc{'T'}(time,:,y_inds,x_inds) + 300;
                            P = nc{'P'}(time,:,y_inds,x_inds) + nc{'PB'}(time,:,y_inds,x_inds);
                            T = potemp ./ ( (1e5./P).^0.286 );
                            qv = nc{'QVAPOR'}(time,:,y_inds,x_inds);
                            cont_data = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 );
                    end
                case 'Temperature (^{o}C)'
                    potemp = nc{'T'}(time,:,y_inds,x_inds) + 300;
                    P = nc{'P'}(time,:,y_inds,x_inds) + nc{'PB'}(time,:,y_inds,x_inds);
                    cont_data = potemp ./ ( (1e5./P).^0.286 ) -273.15;   
                case 'Wind speed (m s^{-1})'
                    cont_data = sqrt( (0.5 * ( nc{'U'}(time,:,y_inds,x_inds(1):x_inds(end)) + nc{'U'}(time,:,y_inds,x_inds(2):x_inds(end)+1) )).^2 ...
                    + (0.5 * ( nc{'V'}(time,:,y_inds(1):y_inds(end),x_inds) + nc{'V'}(time,:,y_inds(2):y_inds(end)+1,x_inds) ) ).^2 );
                case 'Vertical wind speed (m s^{-1})'
                    cont_data = 0.5 * ( nc{'W'}(time,1:end-1,y_inds,x_inds) + nc{'W'}(time,2:end,y_inds,x_inds) );
                case 'Horizontal wind speed (m s^{-1})'
                    cont_data = 0.5 * ( nc{'U'}(time,:,y_inds,x_inds(1):x_inds(end)) + nc{'U'}(time,:,y_inds,x_inds(2):x_inds(end)+1) );
                case 'Wind direction (degrees)'
                    Udir = 0.5 * ( nc{'U'}(time,:,y_inds,x_inds(1):x_inds(end)) + nc{'U'}(time,:,y_inds,x_inds(2):x_inds(end)+1) );
                    Vdir = 0.5 * ( nc{'V'}(time,:,y_inds,x_inds(1):x_inds(end)) + nc{'V'}(time,:,y_inds,x_inds(2):x_inds(end)+1) );
                    theta2 = 180/pi * atan ( Udir ./ Vdir );
                    cont_data = theta2;
                        
                            iuv=find(Udir==0 & Vdir==0);
                            cont_data(iuv) = 0;
                            
                            iuv=find(Udir>=0 & Vdir>=0);
                            cont_data(iuv) = theta2(iuv);
                            
                            iuv = find(Udir>0 & Vdir<0);  %theta2 is negative
                            cont_data(iuv) = 180 + theta2(iuv);
                            
                            iuv = find(Udir<=0 & Vdir<=0);
                            cont_data(iuv) = 180 + theta2(iuv);
                            
                            iuv=find(Udir<0 & Vdir>0);
                            cont_data(iuv) = 360 + theta2(iuv); %theta2 is negative
                            
                            cont_data=cont_data+180; %convert to direction from
                            iuv=find(cont_data>=360);
                            cont_data(iuv) = cont_data(iuv)-360; 
                            
                case 'Density (kg m^{-3})'
                    potemp = nc{'T'}(time,:,y_inds,x_inds) + 300;
                    P = nc{'P'}(time,:,y_inds,x_inds) + nc{'PB'}(time,:,y_inds,x_inds);
                    T = potemp ./ ( (1e5./P).^0.286 );
                    cont_data = density(P,T);   
                    
                case 'Richardson Number (perpendicular wind)'   
                    cont_data = nc{'T'}(time,:,y_inds,x_inds) + 300;

                    
                    
                otherwise %inc 'Component wind speed...' - calculated later
                    cont_data =  ones(size(nc{'T'}(time,:,y_inds,x_inds)));
                    
                    

                    
                                        

            end
            
            
            
            
            for i_line=n_line
        %N.B. ordering the 2D grids as (lat,lon) or (y,x) is correct as in Matlab the first index is
        %row number and so (1,:) should be the first row of the data, i.e. all the x values for y=1
                terr_slice(i_line) = interp2(grid_x,grid_y,squeeze(terr),xd(i_line),yd(i_line));
                lat_slice(i_line) = interp2(grid_x,grid_y,squeeze(lat_restrict),xd(i_line),yd(i_line));
                lon_slice(i_line) = interp2(grid_x,grid_y,squeeze(lon_restrict),xd(i_line),yd(i_line));
                
                
                for iz=1:size(z,1)
%                    u_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(u(iz,:,:)),yd(i_line),xd(i_line));
                    u_slice(i_line,iz) = interp2(grid_x,grid_y,squeeze(u(iz,:,:)),xd(i_line),yd(i_line));
                    v_slice(i_line,iz) = interp2(grid_x,grid_y,squeeze(v(iz,:,:)),xd(i_line),yd(i_line));
                    w_slice(i_line,iz) = interp2(grid_x,grid_y,squeeze(w(iz,:,:)),xd(i_line),yd(i_line));
                    z_slice(i_line,iz) = interp2(grid_x,grid_y,squeeze(z(iz,:,:)),xd(i_line),yd(i_line));                                        
                    
%                     v_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(v(iz,:,:)),yd(i_line),xd(i_line));
%                     w_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(w(iz,:,:)),yd(i_line),xd(i_line));
%                     z_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(Z.var(iz,:,:)),yd(i_line),xd(i_line));
%                     terr_slice(i_line,iz) = interp2(x_grid,y_grid,squeeze(nc{'HGT'}(time,:,:)),yd(i_line),xd(i_line));
                    
                    cont_slice(i_line,iz) = interp2(grid_x,grid_y,squeeze(cont_data(iz,:,:)),xd(i_line),yd(i_line));                    

                end %now should have profiles at this position along the cross-section line

            end
            
            


        case 1
            ipos_line = findheight(y_grid,y_cross_sec);
            for i_line=n_line

                switch var_plot
                    case 'v'
                        cont_slice(i_line-n_line(1)+1,:) = 0.5* (nc{'V'}(time,:,ipos_line,i_line) + nc{'V'}(time,:,ipos_line+1,i_line) ); %2d wind at one height
                    case 'potemp'
                        cont_slice(i_line-n_line(1)+1,:) = nc{'T'}(time,:,ipos_line,i_line) + 300;  %potential temp (300 is base potemp for WRF and T the pert)
                    
                        

                end

                %           w_slice(i_line,:) = 0.5 * ( nc{'W'}(time,2:end,ipos_line,i_line) + nc{'W'}(time,1:end-1,ipos_line,i_line) );
                %           u_slice(i_line,:) = 0.5* (nc{'U'}(time,:,ipos_line,i_line) + nc{'U'}(time,:,ipos_line,i_line+1) ); %2d wind at one height
                z_slice(i_line-n_line(1)+1,:) = WRFUserARW(nc,'Z',time,ipos_line,i_line);

                %w_slice(i_line,:) = w(:,ipos_line,i_line);
                %u_slice(i_line,:) = u(:,ipos_line,i_line);
                %    cont_slice(i_line,:) = v(:,ipos_line,i_line);
                %cont_slice(i_line,:) = potemp(:,ipos_line,i_line);
                %z_slice(i_line,:) = Z.var(:,ipos_line,i_line);

                %       terr_slice(i_line-n_line(1)+1,:) = nc{'HGT'}(time,ipos_line,i_line);
            end

        case 2
            ipos_line = findheight(y_grid,y_cross_sec);
            
            switch file_type
                case {'wrfout','wrfinput'}
                    z_temp=(nc{'PH'}(time,:,ipos_line,n_line) + nc{'PHB'}(time,:,ipos_line,n_line) )./9.81;
                    z_slice = 0.5.*(z_temp(1:end-1,:)+z_temp(2:end,:))';
                case 'met_em'
                    isigma_levs=0;
                    if prod(size(nc{'GHT'})) > 0
                        %NOTE - GHT height is the height of the given pressure level above height zero for the model - i.e. abv mean sea level height.
                        %The first pressure level given in the met_em files is the surface pressure (pressure at the terrain level)
                        %For analysis on pressure levels the heights of the pressures can be below the surface height
                        
                        %pres=nc{'PRES'}(time,2:end,ipos_line,n_line)/100;
                        %z_slice = (find_height_from_p_ant_d03(pres))';
                        
                        z_slice=(nc{'GHT'}(time,2:end,ipos_line,n_line)');                        
                    else
                        pres=nc{'PRES'}(time,:,ipos_line,n_line);
                        psfc=nc{'PSFC'}(:); psfc=psfc(ipos_line,n_line);
                        pmsl=nc{'PMSL'}(:); pmsl=pmsl(ipos_line,n_line);
                        temp=nc{'TT'}(1,:,ipos_line,n_line);
                        soilhgt=nc{'SOILHGT'}(:); soilhgt=soilhgt(ipos_line,n_line);  %height of PSFC level
                        %See watervapourMay2005 for explanation of this height determination
                        
                        clear z_slice
                        for islice2=1:size(pres,2)
                            PSPAN=[pres(1,islice2) pres(end,islice2)]; %pressure range for integration
                            [P,H] = ODE45(@hydrostatic2,PSPAN,soilhgt(islice2),[],pres(:,islice2),temp(:,islice2)); %enter the initial height for the given the first value of PSPAN
                            %note that some of the pressure levels in the PRES array will be below the
                            z_slice(islice2,:) = interp1(P,H,pres(:,islice2));
                        end
                        
                        isigma_levs=1;

                    end
                    
                
                       
                    
                    
            end
                    
                    
            
            
            switch var_plot                    
                case 'potemp'
                    switch file_type
                        case {'wrfout','wrfinput'}
                            cont_slice = nc{'T'}(time,:,ipos_line,n_line)' + 300;  %potential temp (300 is base potemp for WRF and T the pert)
                        case 'met_em'
                            if isigma_levs==0
                                cont_slice = ( nc{'TT'}(time,2:end,ipos_line,n_line) .* (1e5/nc{'PRES'}(time,2:end,ipos_line,n_line)).^0.286 )';  %
                            else
                                cont_slice = ( nc{'TT'}(time,1:end,ipos_line,n_line) .* (1e5/nc{'PRES'}(time,1:end,ipos_line,n_line)).^0.286 )';  %
                            end
                    end

                    
                    
                case 'equiv_potemp'
                    switch file_type
                        case 'wrfout'
                            potemp = nc{'T'}(time,:,ipos_line,n_line) + 300;
                            P = nc{'P'}(time,:,ipos_line,n_line) + nc{'PB'}(time,:,ipos_line,n_line);
                            T = potemp ./ ( (1e5./P).^0.286 );
                            qv = nc{'QVAPOR'}(time,:,ipos_line,n_line);
                            cont_slice = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                            lowest_height = z_slice(:,1);
                        case 'met_em'
                            if isigma_levs==0
                                P = nc{'PRES'}(time,2:end,ipos_line,n_line);
                                T = nc{'TT'}(time,2:end,ipos_line,n_line);
                                rh = nc{'RH'}(time,2:end,ipos_line,n_line);
                            else
                                P = nc{'PRES'}(time,1:end,ipos_line,n_line);
                                T = nc{'TT'}(time,1:end,ipos_line,n_line);
                                rh = nc{'RH'}(time,1:end,ipos_line,n_line);
                            end
                            
                            qsat = satvappress(T,'goff','liq',P,1)/f;
                            qv = rh/100 .* qsat;

                            cont_slice = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                            lowest_height = z_slice(:,1);
                    end
                    
                    case 'Relative humidity (%)'
                            potemp = nc{'T'}(time,:,ipos_line,n_line) + 300;
                            P = nc{'P'}(time,:,ipos_line,n_line) + nc{'PB'}(time,:,ipos_line,n_line);
                            T = potemp ./ ( (1e5./P).^0.286 );
                            qv = nc{'QVAPOR'}(time,:,ipos_line,n_line);
                            qs = SatVapPress(T,'goff','liq',P,1)/f;
                            cont_slice = 100*( qv./qs )';
                            lowest_height = z_slice(:,1); 
                    
             end
                    
                    
                    switch file_type
                        case {'wrfout','wrfinput'}
                            u_slice=0.5 * (nc{'U'}(time,:,ipos_line,n_line) + nc{'U'}(time,:,ipos_line+1,n_line) )';
                            w_slice=0.5 * (nc{'W'}(time,1:end-1,ipos_line,n_line) + nc{'W'}(time,2:end,ipos_line,n_line) )';
                        case 'met_em'
                            if isigma_levs==0
                                u_slice=0.5 * (nc{'UU'}(time,2:end,ipos_line,n_line) + nc{'UU'}(time,2:end,ipos_line+1,n_line) )';
                                w_slice=zeros(size(u_slice)); %no vertical velocity from analysis
                            else
                                u_slice=0.5 * (nc{'UU'}(time,1:end,ipos_line,n_line) + nc{'UU'}(time,1:end,ipos_line+1,n_line) )';
                                w_slice=zeros(size(u_slice)); %no vertical velocity from analysis
                            end
                    end
                    
                    
            end
                        

end
    
%allow the axis type to be changed
switch hor_vert
    case 0
        switch x_axis_type
            case 'dist'
                %timesTH(1).t = d_line'; %already done earlier
            case 'lon';
                timesTH(1).t = lon_slice;
            case 'lat';
                timesTH(1).t = lat_slice;    
        end
        
    case {3}        
                
        switch x_axis_type
            case 'dist'
                %timesTH(1).t = d_line'; %already done earlier
            case 'lon';
                timesTH(1).t = lon_slice;
        end
        
        if iadd_streamline_z0==1
            iadd_overlay=iadd_overlay+1;
            dat_overlay(iadd_overlay).x=timesTH(1).t;
            dat_overlay(iadd_overlay).y=zd;
            dat_overlay(iadd_overlay).overlay_style=['''k--'',''linewidth'',2'];
        end
        if iadd_effective_surface==1
            iadd_overlay=iadd_overlay+1;
            dat_overlay(iadd_overlay).x=X2_pot;
            dat_overlay(iadd_overlay).y=Y2_pot;
            dat_overlay(iadd_overlay).overlay_style=['''w--'',''linewidth'',2'];
        end
        
end

    for i_line=n_line
        %pdat(1).p(:,i_line-n_line(1)+1) = interp1 (z_slice(i_line,:)/1e3,v_slice(i_line,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before^M
        
        wrf_dat(:,i_line-n_line(1)+1) = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,cont_slice(i_line-n_line(1)+1,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before
        u = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,u_slice(i_line-n_line(1)+1,:),zz(1).z);
        v = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,v_slice(i_line-n_line(1)+1,:),zz(1).z);
%%%%%%%%%%%%%%%%%
                    sp = sqrt( u.^2 + v.^2 );
                                                                        
                    clear dir                    
%                    for iuv=1:length(sp)
                        theta2 = 180/pi * atan ( u ./ v );
                        
                            iuv=find(u==0 & v==0);
%                        if u(iuv)==0 & v(iuv)==0
                            dir(iuv) = 0;
                            iuv=find(u>=0 & v>=0);
%                        elseif u(iuv)>=0 & v(iuv)>=0
                            dir(iuv) = theta2(iuv);
                            iuv = find(u>0 & v<0);  %theta2 is negative
%                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            dir(iuv) = 180 + theta2(iuv);
                            iuv = find(u<=0 & v<=0);
%                        elseif u(iuv)<=0 & v(iuv)<=0
                            dir(iuv) = 180 + theta2(iuv);
                            iuv=find(u<0 & v>0);
%                        elseif u(iuv)<0 & v(iuv)>0
                            dir(iuv) = 360 + theta2(iuv); %theta2 is negative
%                        end
%                    end
                                       
                    dir = dir.*pi/180; %convert to radians
                        %%%%% easterly component
%                    spE = sp.*sin(dir); %get the wind direcion in the east component so that can then convert to the horizontal component for the cross section
                        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
        u_quiver(:,i_line-n_line(1)+1) = sp.*cos(pi/2-angle(i_line)-dir); %find component for horizontal axis of the cross section
        %sp is the wind speed magnitude from u and v
        v_quiver(:,i_line-n_line(1)+1) = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,w_slice(i_line-n_line(1)+1,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before
        %or rather w_quiver in this case
        
        switch var_plot
            case 'Component horizontal (UV) wind speed (m s^{-1})' %rather wind speed (vertical and horizontal combined) in direction of cross section
                wrf_dat(:,i_line-n_line(1)+1) = sqrt( u_quiver(:,i_line-n_line(1)+1).^2 + v_quiver(:,i_line-n_line(1)+1).^2 );
            case {'Component (UVW) wind speed (m s^{-1})','Wind speed perpendicular component (m s^{-1})'} %the horizontal wind component (of sqrt(u^2+v^2) in direction of cross section
                wrf_dat(:,i_line-n_line(1)+1) = u_quiver(:,i_line-n_line(1)+1);
            case 'Richardson Number (perpendicular wind)'
                N_col = sqrt( 9.81./wrf_dat(2:end,i_line-n_line(1)+1) .* diff(wrf_dat(:,i_line-n_line(1)+1))./diff(1e3*zz(1).z') );
                dUdz = diff(u_quiver(:,i_line-n_line(1)+1))./diff(1e3*zz(1).z');
                %wrf_dat contains the potential temperature
                wrf_dat(2:end,i_line-n_line(1)+1) = N_col.^2./dUdz.^2;
                wrf_dat(1,i_line-n_line(1)+1) = NaN;

                
        end
    end

%end

pdat(1).p = wrf_dat;

switch hor_vert
    case {0,3,4}   %wi
        terr_height = terr_slice;
    case 2
        switch file_type
            case {'wrfout','wrfinput'}
                terr_height = nc{'HGT'}(time,ipos_line,n_line);
            case 'met_em'
                terr_height = nc{'HGT_M'}(time,ipos_line,n_line);
        end

end
	
	

 
%zz(1).z=zz(1).z*10;


%u_quiver=u_quiver;
%v_quiver=20*v_quiver;
x_quiver=timesTH(1).t;
y_quiver=zz(1).z;
nx_quiver=25; %number of arrows to draw for x
ny_quiver=25; %number for y
scale_speed_quiver = [15 0.2]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed

iadd_overlay=iadd_overlay+1;
dat_overlay(iadd_overlay).x=x_quiver;
dat_overlay(iadd_overlay).y=terr_height/1000;
dat_overlay(iadd_overlay).overlay_style=['''k-'''];

%%%% MODIS plots %%%%%%
%********************************************************************
case 'MODIS_plot'

x_grid = ([1:size(lat2d(1).var,2)]-1) * DX;
y_grid = ([1:size(lat2d(1).var,1)]-1) * DY;

zz(1).z = YLAT;
timesTH(1).t = XLAT;
       

%N.B. when setting short_plot_name earlier, take off the file_str at the end but keep on here (below)
switch short_plot_name
    case  ['Cloud Top Pressure for ' filestr]        
        pdat(1).p =CTP; %hPa   
end









                  
                  
                  


       
%%%%%%%%%%%%%%        

                
           
    


            

    







if(ixlim==1 & iylim==1)
    [xinds(1),xinds(2)] = findheight(timesTH(1).t,xlims(1),xlims(2));
    [yinds(1),yinds(2)] = findheight(timesTH(1).t,ylims(1),ylims(2));
    
    xy_inds = find( timesTH(1).t>=xlims(1) & timesTH(1).t<=xlims(2) &zz(1).z>=ylims(1) & zz(1).z<=ylims(2) );

%    xinds = [xinds(1):xinds(2)];
%    yinds = [yinds(1):yinds(2)];
    timesTH(1).t = timesTH(1).t(xy_inds);
    zz(1).z = zz(1).z(xy_inds);
    ixlim=0;
    iylim=0;
    pdat(1).p=pdat(1).p(xy_inds);
else
%    xinds=1:length(timesTH(1).t(1:pend));
    xinds=1:size(pdat(1).p,1);    
    yinds=1:length(zz(1).z);    
end





%********************************************************************
case 'wrf_plot'
%             zz(1).z = 1:size(lat2d(1).var,1);
%             timesTH(1).t = 1:size(lat2d(1).var,2);                       
%             
% 	    dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));     
% 	    dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));     
% 
% 	zz(1).z = (zz(1).z - 1)*dy_grid;
% 	timesTH(1).t = (timesTH(1).t - 1)*dx_grid;
    
%xinds = 1:size(lat2d(1).var,2); %lat2d is arranged in LAT,LON order so x is second index (LON)
%yinds = 1:size(lat2d(1).var,1);


%zz(1).z = yinds;
%timesTH(1).t = xinds;

%dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
%dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

%zz(1).z = (zz(1).z - 1)*dy_grid;
%timesTH(1).t = (timesTH(1).t - 1)*dx_grid;

%x_grid = ([1:size(lat2d(1).var,2)]-1) * DX;
%y_grid = ([1:size(lat2d(1).var,1)]-1) * DY;

zz(1).z = y_grid;
timesTH(1).t = x_grid;

switch plot_type
    case '3D surf'
        zsurf = squeeze(Zlev);
end

            
    

%             savemem=0;
%             switch savemem
%             case 0
%              pdat(1).p = squeeze(WRFUserARW(nc,'Z',time,ih_wrf)); %height
%		pdat(1).p = WRFUserARW(nc,'p',time,ih_wrf);
%		pdat(1).p = WRFUserARW(nc,'tc',time,ih_wrf);   %temperature
%        pdat(1).p = nc{'TT'}(time,ih_wrf,:,:)-273.15;   %temperature
%        pdat(1).p = nc{'RH'}(time,ih_wrf,:,:);   %RH
%        pdat(1).p = nc{'T'}(time,ih_wrf,:) + 300;   %potential temperature

%         pdat(1).p = ( nc{'PRES'}(time,ih_wrf,:) );
%        pdat(1).p = squeeze(Zlev) - terrain; %height above terrain

%    pdat(1).p = sp_nudging_z8 - sp_no_nudging_z8;

%        pdat(1).p=squeeze(Zlev);
%		pdat(1).p = nc{'W'}(time,ih_wrf,:,:);
%                pdat(1).p = -nc{'LH'}(time,:,:); %latent heat (W/m2)  -minus so that positive gives energy to surface
%                 pdat(1).p = -nc{'HFX'}(time,:,:); %sensible heat (W/m2)

%                pdat(1).p = -nc{'LH'}(time,:,:) - nc{'HFX'}(time,:,:); %latent heat (W/m2)  minus so that positive gives energy to surface
%               pdat(1).p = nc{'TSLB'}(time,soil_level,:,:)-273.15; %soil temperature (degC) - level one is closest to the surface, 4 is deep below
%                pdat(1).p = nc{'Q2'}(time,:,:)*1000; %vapour MR air temperature (degC)

%%%%%%%%%%%%%%%%%%%%%%%%%% difference plots %%%%%%%%%%%%%%%%%%
%                 pdat(1).p=temp_ecmwf-temp_ncep
%                 pdat(1).p=(p270_ecmwf_d02-p270_ncep_d02);
%                 pdat(1).p=(psfc_ecmwf_d01-psfc_ncep_d01)/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   





                
%                pdat(1).p = nc{'GRDFLX'}(time,:,:);

%                pdat(1).p = nc{'QFX'}(time,:,:);

%                pdat(1).p = nc{'SWDOWN'}(time,:,:);
%                pdat(1).p = nc{'ALBEDO'}(time,:,:);
                
%                 pdat(1).p = nc{'SEAICE'}(time,:,:);
%                  pdat(1).p = nc{'LU_INDEX'}(time,:,:);
%                 pdat(1).p = nc{'SNOWC'}(time,:,:);

%                 pdat(1).p = nc{'LANDMASK'}(time,:,:); %includes the ice shelves - same in ncep and ecmwf
%                 pdat(1).p = nc{'LANDSEA'}(time,:,:); %large diffs - ncep looks more blocky and doesn't cover the ice shelf like ecmwf does
                                                       %ecmwf looks better - but perhaps LANDMASK takes precedent - this looks good and the same in both
%                 pdat(1).p = nc{'SEAICE'}(time,:,:);   %bit different but perhaps not enough to explain diffs - same over the ice shelf
%                 pdat(1).p = nc{'SLOPECAT'}(time,:,:);   %large difference especially over the ice shelf
%                 pdat(1).p = nc{'SNOALB'}(time,:,:);     %large changes similar to above               
%                 pdat(1).p = nc{'GREENFRAC'}(time,1,:,:); %large changes similar to above               

%                 pdat(1).p = nc{'ALBEDO12M'}(time,12,:,:); %%large changes similar to above               
%                 pdat(1).p = nc{'SOILCBOT'}(time,16,:,:); % %same
%                 pdat(1).p = nc{'SOILCTOP'}(time,14,:,:); % %same  
%                  pdat(1).p = nc{'SPECHUMD'}(time,:,:); % %not present for ECMWF
%                  pdat(1).p = nc{'SLPY'}(time,:,:); %SLPX and SLPY very similar for both ncep and ecmwf
%                  pdat(1).p = nc{'HGT_M'}(time,:,:); %same
%                  pdat(1).p = nc{'SNOW'}(time,:,:); %very different - units wrong for ecmwf - 10 kg/m2 instead of 10m!
%                  pdat(1).p = nc{'SM'}(time,1,:,:); %think is soil moisture - only two levels for ncep, 4 for ecmwf
%                  pdat(1).p = nc{'ST'}(time,1,:,:)-273.15; %soil temps - since have 4 depth layers in ecmwf, 2 in ncep                  
                  
%                pdat(1).p = nc{'ST010200'}(time,:,:);
%                pdat(1).p = nc{'ST000010'}(time,:,:);                
%                pdat(1).p = nc{'SM000010'}(time,:,:);
%                pdat(1).p = nc{'SM010200'}(time,:,:);

%                pdat(1).p = nc{'ST010200'}(time,:,:);
%                pdat(1).p = nc{'ST000010'}(time,:,:);                
%                pdat(1).p = nc{'SM000010'}(time,:,:);
%                pdat(1).p = nc{'ST100255'}(time,:,:);
%                 pdat(1).p = nc{'SOILTEMP'}(time,:,:)-273.15;
%                 pdat(1).p = nc{'LANDUSEF'}(time,16,:,:);
%                 pdat(1).p = nc{'SOILHGT'}(time,:,:); %large diffs between ncep and ecmwf


%                pdat(1).p = nc{'SKINTEMP'}(time,:,:) - 273.15; %skin temperature for met_em files
%                pdat(1).p = nc{'SST'}(time,:,:) - 273.15; %skin temperature for met_em files - only for ecmwf files - values are very strange
                                                           %range from 1.6e29 to -1.2e30


                
%                pdat(1).p = 0.01* nc{'PMSL'}(time,:,:); %pressure

%                pdat(1).p = nc{'TSLB'}(time,1,:,:) - 273.15; %skin temperature for wrfout files
%                 pdat(1).p = nc{'XLAND'}(time,:,:);

%                pdat(1).p = nc{'ALBEDO'}(time,:,:); %
%                pdat(1).p = nc{'SOILTB'}(time,:,:); % 
%                pdat(1).p = nc{'VEGFRA'}(time,:,:);
%                pdat(1).p = nc{'ISLTYP'}(time,:,:);
%                pdat(1).p = f* (nc{'QICE'}(itime,ih_wrf,:,:)+nc{'QSNOW'}(itime,ih_wrf,:,:)+nc{'QGRAUP'}(itime,ih_wrf,:,:)+nc{'QVAPOR'}(itime,ih_wrf,:,:) );

if ih_wrf~=0
        switch file_type
            case 'wrfinput'
                P=eta_get_p(nc);
                P=squeeze(P(ih_wrf,:,:));
                T = nc{'T'}(1,ih_wrf,:,:) + 300; %this fits the WRFout potential temp - so don't seem to need T_INIT
                T = T./(1e5./P).^0.286;
            otherwise
                P = 100*WRFUserARW(nc,'p',time,ih_wrf);    %Pa
                T = WRFUserARW(nc,'tc',time,ih_wrf)+273.15;   %temperature K
        end
        
        rho=density(P,T);
        
end
%        pdat(1).p = (nc{'QNRAIN'}(time,ih_wrf,:,:)+nc{'QNICE'}(time,ih_wrf,:,:)+nc{'QNSNOW'}(time,ih_wrf,:,:)+nc{'QNGRAUPEL'}(time,ih_wrf,:,:)).*rho/1000; %total number
%        pdat(1).p = (nc{'QNICE'}(time,ih_wrf,:,:)).*rho/1000; %total number
%        pdat(1).p = (nc{'QNSNOW'}(time,ih_wrf,:,:)).*rho/1000; %total number
%        pdat(1).p = (nc{'QNRAIN'}(time,ih_wrf,:,:)).*rho/1000; %total number
%        pdat(1).p = (nc{'QNGRAUPEL'}(time,ih_wrf,:,:)).*rho/1000; %total number
    
%        pdat(1).p = 1000*(nc{'QCLOUD'}(time,ih_wrf,:,:));
%        pdat(1).p = 1000*(nc{'QICE'}(time,ih_wrf,:,:));
%        pdat(1).p = 1000*(nc{'QSNOW'}(time,ih_wrf,:,:));
%        pdat(1).p = 1000*(nc{'QVAPOR'}(time,ih_wrf,:,:));
%        pdat(1).p = 1000*(nc{'QGRAUP'}(time,ih_wrf,:,:));
%        pdat(1).p = 1000*(nc{'QRAIN'}(time,ih_wrf,:,:));
        
%        pdat(1).p = (nc{'SNOWNC'}(time,:,:));
%        pdat(1).p = (nc{'RAINNC'}(time,:,:));   
        
%        pdat(1).p = (nc{'RAINNC'}(time,:,:)) - (nc{'RAINNC'}(time-1,:,:));   %precip tendency (since the last output time)

%        pdat(1).p = 1000*(nc{'QCLOUD'}(time,ih_wrf,:,:)+nc{'QICE'}(time,ih_wrf,:,:)+nc{'QSNOW'}(time,ih_wrf,:,:)+nc{'QGRAUP'}(time,ih_wrf,:,:)); %total condensate



%N.B. when setting short_plot_name earlier, take off the file_str at the end but keep on here (below)
switch short_plot_name
    case ['RH (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
        if ih_wrf==0
            qv=nc{'Q2'}(time,:,:);
            T2=nc{'T2'}(time,:,:);
            P2=nc{'PSFC'}(time,:,:);
            qvs=satvappress(T2,'goff','liq',P2,1)/f;
            pdat(1).p =100*(qv./qvs);
        else
            qv=nc{'QVAPOR'}(time,ih_wrf,:,:);
            qvs=satvappress(T,'goff','liq',P,1)/f;
            pdat(1).p =100*(qv./qvs);
        end
    case ['Ice supersaturation (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
        qv=nc{'QVAPOR'}(time,ih_wrf,:,:);
        qvs=satvappress(T,'goff','ice',P,1)/f;
        pdat(1).p=100*(qv./qvs-1);        
%    case ['Temperature (^{o}C) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
    case ['Temperature (^{o}C) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC at level ' num2str(ih_wrf) ' ' medZ ' for ' filestr];                 
        switch file_type
            case 'met_em'
                pdat(1).p = nc{'TT'}(time,ih_wrf,:,:)-273.15;   %temperature
            otherwise
                pdat(1).p=WRFUserARW(nc,'tc',time,ih_wrf);
        end

    case ['Hallet Mossop flag at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]    
        qv=nc{'QVAPOR'}(time,ih_wrf,:,:);
        qvs=satvappress(T,'goff','ice',P,1)/f;
        si=100*(qv./qvs-1);  
                    
        qvsw=satvappress(T,'goff','liq',P,1)/f;
        sw =100*(qv./qvsw-1);    
        
        ql=1000*nc{'QCLOUD'}(time,ih_wrf,:,:); %g/kg
        qi=1000*nc{'QICE'}(time,ih_wrf,:,:);        
        
        tc=WRFUserARW(nc,'tc',time,ih_wrf);
        ihm=find(tc<=-3 & tc>=-8 & ql>0.01 & qi>2e-4);
        pdat(1).p=zeros(size(tc));
        pdat(1).p(ihm)=1;
        
    case ['Sea Ice flag for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
        pdat(1).p=nc{'SEAICE'}(time,:,:);
    case ['Skin temperature for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C) for ' filestr];    
                switch file_type
            case 'met_em'
                    pdat(1).p = nc{'SKINTEMP'}(time,:,:) - 273.15; %skin temperature for met_em files
            otherwise                        
                    pdat(1).p = nc{'TSK'}(time,:,:) - 273.15; %skin temperature for wrfout files
            end
    case ['Wind speed (m s^{-1}) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC at level ' num2str(ih_wrf) ' ' medZ ' for ' filestr];                 
        switch file_type
            case 'met_em'
                u = nc{'UU'}(1,ih_wrf,:,:);
                v = nc{'VV'}(1,ih_wrf,:,:);
                u=0.5 * ( u(:,1:end-1) + u(:,2:end) );
                v=0.5 * ( v(1:end-1,:) + v(2:end,:) );
            otherwise
                u=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
                v=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
        end
        pdat(1).p = sqrt(u.^2+v.^2);
                        
    case ['Terrain height ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m) for ' filestr]   
        switch file_type
            case 'met_em'
                pdat(1).p = nc{'HGT_M'}(time,:,:); %terrain height
            otherwise
                pdat(1).p = nc{'HGT'}(time,:,:); %terrain height
        end
    case ['First level height ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m) for ' filestr]   
                pdat(1).p = terr_fine; %terrain height of first available data for streamlines 
    case ['10 m wind speed for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
                 u10 = nc{'U10'}(time,:,:); %10 m winds
                 v10 = nc{'U10'}(time,:,:); %
                 pdat(1).p = sqrt(u10.^2+v10.^2);
    case ['Surface pressure (hPa) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]             
                 pdat(1).p = 0.01* nc{'PSFC'}(time,:,:); %surface pressure
                 
    case ['Pressure(hPa) at model level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' ' month_str ' ' Times(time,12:16) ' UTC for ' filestr];                 
                 pdat(1).p = ( nc{'P'}(time,ih_wrf,:) + nc{'PB'}(time,ih_wrf,:) ) / 100; %pressure
    case ['Total water at level ' num2str(ih_wrf) ' (ppmv) for ' filestr]
                 pdat(1).p = f*(nc{'QVAPOR'}(time,ih_wrf,:,:)+nc{'QCLOUD'}(time,ih_wrf,:,:)+nc{'QICE'}(time,ih_wrf,:,:)+nc{'QSNOW'}(time,ih_wrf,:,:)+nc{'QGRAUP'}(time,ih_wrf,:,:)); %total water (ppmv)
    case ['Air density at level ' num2str(ih_wrf) ' (kg m^{-3}) for ' filestr]
                 pdat(1).p = rho;
    case ['2m air temperature ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C) for ' filestr];             
                 pdat(1).p = nc{'T2'}(time,:,:)-273.15; %2m air temperature (degC)
    case ['2m vapour mixing ratio ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C) for ' filestr];             
                 pdat(1).p = nc{'Q2'}(time,:,:); %2m vapour MR (kg/kg)
                 
                 
end

    if length(strfind(short_plot_name,'Total melting (mm w.e.)'))>0
        pdat(1).p = melt_tot; %run calc_melt_tot to get this
    end    
    if length(strfind(short_plot_name,'Total melt flux (W m^{-2})'))>0
        pdat(1).p = melt_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'SW net flux (W m^{-2})'))>0
        pdat(1).p = sw_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'SW downwards flux (W m^{-2})'))>0
        pdat(1).p = sw_tot./(1-ALBEDO) * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'LW down flux (W m^{-2})'))>0
        pdat(1).p = lwdown_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'LW net flux (W m^{-2})'))>0
        pdat(1).p = lw_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'SH flux (W m^{-2})'))>0
        pdat(1).p = sh_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'LH flux (W m^{-2})'))>0
        pdat(1).p = lh_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    if length(strfind(short_plot_name,'GRD flux (W m^{-2})'))>0
        pdat(1).p = grd_tot * 3.34e5/(time_step*3600); %run calc_melt_tot to get this
    end
    


%pdat(1).p = nc{'GLW'}(time,:,:); %Downward LW flux
%pdat(1).p = lwdown_tot; %Net LW flux - mean over whole 3 days
%pdat(1).p = temp_tot; %Net LW flux
%pdat(1).p = nc{'SWDOWN'}(time,:,:); %Downward SW flux
%pdat(1).p = nc{'OLR'}(time,:,:); %TOA outgoing LW

%pdat(1).p = sw_tot; %
%pdat(1).p = lw_tot; %
%pdat(1).p = 1000*cond_tot; %average total condensate (column total g/m2)
% pdat(1).p =  temp_tot-273.15; %
%pdat(1).p = sh_tot; %average melt in mm for 12,15,18 and 21 UTC
%pdat(1).p = rh_iz_tot; %average relative humidity


% iadd_overlay=0;
% i2_save2=i2_save(:);
% i2_save2(i2_save2==0)=[];
% i2_save2=unique(i2_save2);
% [ix_overlay, iy_overlay] = ind2sub(size(melt_tot),i2_save2);
% x_overlay=timesTH(1).t(iy_overlay);
% y_overlay=zz(1).z(ix_overlay);

% [ix_overlay, iy_overlay] = find(n_tot2_restrict>=1);
% x_overlay=timesTH(1).t(ix_overlay);
% y_overlay=zz(1).z(iy_overlay);

iwind_field=0;
if iwind_field==1
    u=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
    v=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
    pdat(1).p = sqrt(u.^2+v.^2);
end

i_cloud_layer=0;
if i_cloud_layer==1;
    pdat(1).p = 1000*nc{'QCLOUD'}(time,ih_wrf,:,:) + 1000*nc{'QRAIN'}(time,ih_wrf,:,:)...
        + 1000*nc{'QICE'}(time,ih_wrf,:,:) + 1000*nc{'QSNOW'}(time,ih_wrf,:,:)...
        + 1000*nc{'QGRAUP'}(time,ih_wrf,:,:);
end

i_condensate=0;
recalc=1;
if i_condensate==1
    if recalc==1
        ih_inds=1:10;

        M = 28.97*1.67E-27;
        k = 1.38E-23;
        P = WRFUserARW(nc,'p',time); P.var=P.var*100;
        T = WRFUserARW(nc,'tc',time); T.var=T.var+273.15;
        rho=P.var(ih_inds,:,:).*M./k./T.var(ih_inds,:,:);
        Z = squeeze(WRFUserARW(nc,'Z',time));
        dz_grid = Z.var(ih_inds+1,:,:)-Z.var(ih_inds,:,:);



        cond_dat = squeeze (sum( (nc{'QCLOUD'}(time,ih_inds,:,:) + nc{'QRAIN'}(time,ih_inds,:,:)...
            + nc{'QICE'}(time,ih_inds,:,:) + nc{'QSNOW'}(time,ih_inds,:,:)...
            + nc{'QGRAUP'}(time,ih_inds,:,:)) .* rho*1000*dx_grid*1000*dy_grid.*dz_grid ,1) );
    else
        fprintf(1,'\n****** WARNING - not recalculating melt energy values *******');
    end


    pdat(1).p = cond_dat;

end

i_max_cloud=0;
recalc=1;
if i_max_cloud==1
    if recalc==1  

%         M = 28.97*1.67E-27;
%         k = 1.38E-23;
%         
%         P = nc{'P'}(time,1:ih_inds,:) + nc{'PB'}(time,1:ih_inds,:);       
%         potemp = nc{'T'}(time,1:ih_inds,:) + 300;
%         T = potemp ./ ( (1e5./P).^0.286 );
%         
%         rho=P.*M./k./T;
%         z_temp=(nc{'PH'}(time,1:ih_inds+1,:) + nc{'PHB'}(time,1:ih_inds+1,:) )./9.81;
%         Z = 0.5.*(z_temp(1:end-1,:)+z_temp(2:end,:));
%         
%         dz_grid = Z(2:end,:,:)-Z.var(1:end-1,:,:);
%         cond_dat = squeeze ( sum( nc{'QCLOUD'}(time,1:ih_inds,:,:) .* rho*1000*dx_grid*1000*dy_grid.*dz_grid , 2) );
        
        max_var='cloud';
%        max_var='ice no';
%        max_var='snow no';
        
        %whether to calculate the density or not
        switch max_var
            case {'ice no','snow no'}
                potemp = nc{'T'}(time,1:ih_wrf,:,:) + 300;
                P = nc{'P'}(time,1:ih_wrf,:,:) + nc{'PB'}(time,1:ih_wrf,:,:);
                T = potemp ./ ( (1e5./P).^0.286 );
                rho=density(P,T);
                
        end
        
        %no calculate the data
        switch max_var
            case 'cloud'
                var_dat = 1000*squeeze(nc{'QCLOUD'}(time,1:ih_wrf,:,:));
            case 'ice no'                
                var_dat = squeeze(nc{'QNICE'}(time,1:ih_wrf,:,:)).*rho/1000;
            case 'snow no'                
                var_dat = squeeze(nc{'QNSNOW'}(time,1:ih_wrf,:,:)).*rho/1000;    
        end
        
        %calculate the max of the data
        if ih_wrf>1
            [max_var_dat, i_max_cloud_height] = (max(var_dat));
            max_var_dat=squeeze(max_var_dat);
        else
            fprintf(1,'ih_wrf=1 SO MAX OVER HEIGHT DOESN''T MAKE SENSE!');
            break
        end

        
    else
        fprintf(1,'\n****** WARNING - NOT recalculating MAX values *******');
    end


    pdat(1).p = max_var_dat;

end

iequiv=0;
if iequiv==1 %equivalent potential temperature
    potemp = nc{'T'}(time,ih_wrf,:) + 300;
    P = nc{'P'}(time,ih_wrf,:) + nc{'PB'}(time,ih_wrf,:);
    T = potemp ./ ( (1e5./P).^0.286 );
    qv = nc{'QVAPOR'}(time,ih_wrf,:);
    pdat(1).p = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
end
      

isnow_depth=0;
if isnow_depth==1;
%    snowH=nc{'SNOWH'}(:); %snow height is wrong for polar WRF Noah set up (density applied to whole snow layer)
%    pdat(1).p = squeeze(snowH(end,:,:)-snowH(1,:,:));
    
    snow=nc{'SNOW'}(:,:);
    pdat(1).p = squeeze(snow(end,:,:)-snow(1,:,:));
end

imelt_energy=0;
recalc=1;
if imelt_energy==1
    if recalc==1 %to avoid unneccessary loading and recalculation
        LW=nc{'GLW'}(time,:); %downwelling LW at surface
        SW=nc{'SWDOWN'}(time,:); %downwelling SW at surface
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
        LW_NET=EMISS.*(LW-LW_UP); %Should be LW_NEW = LW - EMISS.*LW_UP; ??Since emissivity just affects the surface upwelling?
          %and not the downwelling (LW)?

        MNET=LW_NET + SW_NET + LH+SH + GRDFLX; %total melt rate W/m2
%        MNET=SW_NET;
%        MNET=SH+LH;
    else
        fprintf(1,'\n****** WARNING - not recalculating melt energy values *******');
    end
    
    if length(strfind(short_plot_name,'Melt energy'))>0
        pdat(1).p=MNET;
    end

    if length(strfind(short_plot_name,'Average melt rate (mm day^{-1} w.e.'))>0
        pdat(1).p = 24*3600*MNET / 3.34e5;
    end
            
    
end
                  
                  
                  


              %  v = 0.5* (nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) ); %2d wind at one height
              
%              pdat(1).p = nc{'GHT'}(time,ih_wrf,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%
inversion_find=0; %

inversion_method = 'quick';
inversion_method = 'slow - actual_heights';

inversion_type='max over height';
%inversion_type='average';
%%%%%%%%%%%%%%%%%%%%%%%%    
    recalc=1;
    
    if inversion_find==1

        if recalc==1

            clear Z_inv potemp wrf_dat


            switch file_type
                case {'wrfout','wrfinput'}
                    Zlev=WRFUserARW(nc,'Z',time,1,1);
                    terrain = nc{'HGT'}(time,:);
                case 'met_em'
                    Zlev=nc{'GHT'}(time,:);
                    terrain = nc{'HGT_M'}(time,:);
            end

            switch inversion_method
                case 'quick'

                    [temp i0]=minALL(terrain);

                    switch file_type
                        case {'wrfout','wrfinput'}
                            heights=WRFUserARW(nc,'Z',time,i0(1),i0(2)); %get a profile over low terrain
                        case 'met_em'
                            heights=get_ecmwf_ml_met_em_heights(nc,i0(1),i0(2));
                    end


                    [ih_wrf_1,ih_wrf_2]=findheight_nearest(heights,inversion_height1,inversion_height2);
                    %find the approx model levels for the required heights
                    %things may vary over terrain though -best to adopt a wide layer

                    %now get the potemp and heights so can calculate the potemp gradient


                    switch file_type
                        case {'wrfout','wrfinput'}
                            Z_inv =(nc{'PH'}(time,ih_wrf_1:ih_wrf_2,:) + nc{'PHB'}(time,ih_wrf_1:ih_wrf_2,:) )./9.81;
                            potemp = nc{'T'}(time,ih_wrf_1:ih_wrf_2,:) + 300;
                        case 'met_em'
                            for ilat=1:size(terrain,1)
                                for ilon=1:size(terrain,2)
                                    profile=get_ecmwf_ml_met_em_heights(nc,ilat,ilon);
                                    Z_inv(1,ilat,ilon)=profile(ih_wrf_1);
                                    Z_inv(2,ilat,ilon)=profile(ih_wrf_2);
                                end
                            end

                            Tinv = nc{'TT'}(1,i0(1),:,:);
                            Pinv = nc{'PRES'}(1,i0(1),:,:);
                            potemp(1,:,:) = Tinv.*(1000./Pinv).^0.286;

                            Tinv = nc{'TT'}(1,i0(2),:,:);
                            Pinv = nc{'PRES'}(1,i0(2),:,:);
                            potemp(2,:,:) = Tinv.*(1000./Pinv).^0.286;

                    end

                    wrf_dat = squeeze( 1000*( potemp(end,:,:) - potemp(1,:,:) ) ./ ( Z_inv(end,:,:) - Z_inv(1,:,:) ) ); %(K km^{-1})

                case 'slow - actual_heights'
                    [temp i0]=minALL(terrain);



                    switch file_type
                        case {'wrfout','wrfinput'}
                            height_low=WRFUserARW(nc,'Z',time,i0(1),i0(2)); %get a profile over low terrain
                              %to estimate the height levels we might require so don't have to load them all in
                            [ih_wrf_1,ih_wrf_2]=findheight_nearest(height_low,inversion_height1,inversion_height2);
                            ih_wrf_1_2 = max(ih_wrf_1-2,1);
                            ih_wrf_2_2 = min(ih_wrf_2+2,length(height_low)); %go 2 indices either side

                            height =(nc{'PH'}(time,ih_wrf_1_2:ih_wrf_2_2,:) + nc{'PHB'}(time,ih_wrf_1_2:ih_wrf_2_2,:) )./9.81;
                            potemp = nc{'T'}(time,ih_wrf_1_2:ih_wrf_2_2,:) + 300;
                        case 'met_em'
                            %%
                    end
                    
                    
                switch inversion_type
                    case 'max over height'
                        dz_inv=25;
                        sample_heights = [inversion_height1:50:inversion_height2];
                        for ilat=1:size(terrain,1)
                            for ilon=1:size(terrain,2)
                                pots = interp1(height(:,ilat,ilon),potemp(:,ilat,ilon),sample_heights);
                                wrf_dat(ilat,ilon) = max(1000* diff(pots)./diff(sample_heights));
                            end
                        end

                    case 'average'
                        sample_heights = [inversion_height1:25:inversion_height2];
                        for ilat=1:size(terrain,1)
                            for ilon=1:size(terrain,2)
                                pots = interp1(height(:,ilat,ilon),potemp(:,ilat,ilon),[inversion_height1 inversion_height2]);                                
                                wrf_dat(ilat,ilon) = 1000* diff(pots)./diff([inversion_height1 inversion_height2]);
                            end
                        end

                end

            end %switch inversion_method




        end

        pdat(1).p=wrf_dat;



    end
%%%%%%%%%%%%%%        

                
%%%%%%%%%%%%%%%%%%%%%%%%
constant_height=0; %for plotting a constant height (rather than model level) by in interpolation
%%%%%%%%%%%%%%%%%%%%%%%%    
    recalc=0;
    reread=1;
    con_height_winds=1; %set to one to also interpolate the winds - otherwise it uses the closest model level to h_wrf
    
    %choose the surface to interpolate onto
    surface='height';
%    surface='pressure';
%    surface='potemp';            
    
    %choose the data to interpolate
%    con_height_case = 'wind';
    con_height_case = 'pressure';
    
	if constant_height==1
        
        if con_height_winds==0
            switch file_type
                case {'wrfout','wrfinput'}
                    Zlev=WRFUserARW(nc,'Z',time,1,1);
                    terrain = nc{'HGT'}(time,1,:);
                    terrain = terrain(1);                    
                case 'met_em'
                        Zlev=nc{'GHT'}(time,:,1,1);
                        terrain = nc{'HGT_M'}(time,1,:);
                        terrain = terrain(1);                                            
            end

            ih_wrf_quiver = findheight_nearest(Zlev-terrain,h_wrf*1000);
        end
        
		if recalc==1
            clear wrf_zint wrf_zint_u wrf_zint_v
            if reread==1
                %	wrf_dat = 1000*nc{'QVAPOR'}(time,:,:,:);
                %                wrf_dat = WRFUserARW(nc,'tc',time); %load in all data for interpolation (set h_wrf in first set of settings)
                %				wrf_dat =  squeeze(WRFUserARW(nc,'p',time));
                %                wrf_dat=eta_get_p(nc)/100; %for wrfinput files

                switch con_height_case
                    case 'wind'
                        wrf_dat = sqrt( (0.5 * ( nc{'U'}(time,:,:,1:end-1) + nc{'U'}(time,:,:,2:end) )).^2 ...
                            + (0.5 * ( nc{'V'}(time,:,1:end-1,:) + nc{'V'}(time,:,2:end,:) ) ).^2 );

                    case 'pressure'
                        switch file_type
                            case 'wrfinput'
                                wrf_dat = eta_get_p(nc)/100; %get pressure from the eta levels (quite accurate but can be slightly off? - 2mb or so)
                                % and convert Pa output from eta_get_p to mb
                            otherwise
                                wrf_dat =  squeeze(WRFUserARW(nc,'p',time));
                        end
                end

                if con_height_winds==1
                    u_wind_dat=0.5 * ( nc{'U'}(time,:,:,1:end-1) + nc{'U'}(time,:,:,2:end) );
                    v_wind_dat=0.5 * ( nc{'V'}(time,:,1:end-1,:) + nc{'V'}(time,:,2:end,:) );
                end


                if prod(size(wrf_dat))==1  %if WRFUserARW returns the structure as wrf_dat.var
                    wrf_dat = wrf_dat.var;
                end

                switch surface
                    case 'pressure'
                        switch file_type
                            case 'wrfinput'
                                Z_wrf = eta_get_p(nc); %if want plot at a constant pressure
                                h_wrf = p_wrf; %
                            otherwise
                                Z_wrf = squeeze(WRFUserARW(nc,'p',time)); %if want plot at a constant pressure
                                h_wrf = p_wrf; %
                        end
                    case 'height'
                        Z_wrf.var = squeeze(WRFUserARW(nc,'Z',time));
                        Z_wrf.var = Z_wrf.var/1000;
                    case 'potemp'
                        Z_wrf.var = nc{'T'}(time,:) + 300;   %potential temperature;
                        h_wrf = potemp_wrf;
                end
                
            end

            
%h_wrf set in first settings bit

igriddata=0;  %don't use gridddata as it takes too long and/or crashes
if igriddata==1  %if want to use the griddata3 method of interpolation instead
    
    x_grid = ([1:size(lat2d(1).var,2)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
    y_grid = ([1:size(lat2d(1).var,1)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

    [X2D,Y2D] = MESHGRID(x_grid,y_grid);
    X=permute(repmat(X2D,[1 1 size(Z_wrf.var,1)]) , [3 1 2]);
    Y=permute(repmat(Y2D,[1 1 size(Z_wrf.var,1)]) , [3 1 2]);


%    I=[1:size(Z_wrf.var,1)];

        %                         wrf_zint(ilat,ilon) = interp1(Z_wrf.var(:,ilat,ilon),wrf_dat(:,ilat,ilon),h_wrf*1000,[],'extrap');
        wrf_zint = griddata3(X,Y,Z_wrf.var,wrf_dat,X2D,Y2D,h_wrf*ones(size(X2D)));
        if con_height_winds==1
            wrf_zint_u(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),u_wind_dat(1:length(I),ilat,ilon),h_wrf);
            wrf_zint_v(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),v_wind_dat(1:length(I),ilat,ilon),h_wrf);
        end

                         
else

                 for ilat=1:size(wrf_dat,2);
                     for ilon=1:size(wrf_dat,3);
    
                         switch surface
                             case 'potemp'
                                 [temp,I]=unique(Z_wrf.var(:,ilat,ilon));
                                 Z_wrf.var(1:length(I),ilat,ilon)=temp;
                                 wrf_dat(1:length(I),ilat,ilon)=wrf_dat(I,ilat,ilon);
                                 if con_height_winds==1
                                     u_wind_dat(1:length(I),ilat,ilon)=u_wind_dat(I,ilat,ilon);
                                     v_wind_dat(1:length(I),ilat,ilon)=v_wind_dat(I,ilat,ilon);
                                 end
                             otherwise
                                 I=[1:size(Z_wrf.var,1)];
                         end

%                         wrf_zint(ilat,ilon) = interp1(Z_wrf.var(:,ilat,ilon),wrf_dat(:,ilat,ilon),h_wrf*1000,[],'extrap');
                         wrf_zint(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),wrf_dat(1:length(I),ilat,ilon),h_wrf);
                         if con_height_winds==1
                             wrf_zint_u(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),u_wind_dat(1:length(I),ilat,ilon),h_wrf);
                             wrf_zint_v(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),v_wind_dat(1:length(I),ilat,ilon),h_wrf);
                         end
                     end
                 end
                 
end

		end

	    pdat(1).p = wrf_zint;
        if con_height_winds==1
            u_quiver = wrf_zint_u;
            v_quiver = wrf_zint_v;
        end
	end
            %    pdat(1).p=squeeze(WRFUserARW(nc,'Z',time,0,0,ih_wrf));  %horiz_slice at first level
            
        if constant_height==0
            ih_wrf_quiver = ih_wrf;
        end
    

if constant_height==0 | con_height_winds==0

    %ih_wrf_quiver = 4;

    switch file_type
        case {'wrfout','wrfinput'}
            if ih_wrf==0 %if have requested 10m wind speeds
                u_quiver=nc{'U10'}(time,:);
                v_quiver=nc{'V10'}(time,:);
            else
                u_quiver=0.5 * ( nc{'U'}(time,ih_wrf_quiver,:,1:end-1) + nc{'U'}(time,ih_wrf_quiver,:,2:end) );
                v_quiver=0.5 * ( nc{'V'}(time,ih_wrf_quiver,1:end-1,:) + nc{'V'}(time,ih_wrf_quiver,2:end,:) );
            end
        case 'met_em'
            UU = nc{'UU'}(1,ih_wrf_quiver,:,:);
            VV = nc{'VV'}(1,ih_wrf_quiver,:,:);
            u_quiver=0.5 * ( UU(:,1:end-1) + UU(:,2:end) );
            v_quiver=0.5 * ( VV(1:end-1,:) + VV(2:end,:) );
    end

end





pdat_orig = pdat(1).p; %save the full field as could be useful
if(ixlim==1 & iylim==1)
    [xinds(1),xinds(2)] = findheight_nearest(timesTH(1).t,xlims(1),xlims(2));
    [yinds(1),yinds(2)] = findheight_nearest(timesTH(1).t,ylims(1),ylims(2));
        u_quiver = u_quiver(yinds(1):yinds(2),xinds(1):xinds(2));
        v_quiver = v_quiver(yinds(1):yinds(2),xinds(1):xinds(2));
            
    xinds = [xinds(1):xinds(2)];
    yinds = [yinds(1):yinds(2)];
    timesTH(1).t = timesTH(1).t(xinds);
    zz(1).z = zz(1).z(yinds);
    ixlim=0;
    iylim=0;
    pdat(1).p=pdat(1).p(yinds,xinds);
    

else
%    xinds=1:length(timesTH(1).t(1:pend));
    xinds=1:size(pdat(1).p,1);    
    yinds=1:length(zz(1).z);    
end

x_quiver=timesTH(1).t;
y_quiver=zz(1).z;



if is_met_em==1
    short_plot_name = [short_plot_name ' (analysis)'];
    tit(1).tit = [tit(1).tit ' (analysis)'];
end

if strcmp(file_type,'wrfinput')==1
    short_plot_name = [short_plot_name ' (analysis)'];
    tit(1).tit = [tit(1).tit ' (analysis)'];
end



    case 'wrf_wind2d'

%		xinds = 83:166;
%		yinds = 83:166;

%		xinds = 1:156;
%		yinds = 20:176;

%		xinds = 204:357;
%		yinds = 77:230;

        
        
            xinds = 1:size(lat2d(1).var,2);
            yinds = 1:size(lat2d(1).var,1);

        
            zz(1).z = yinds;
            timesTH(1).t = xinds;

	    dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));     
	    dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));     

	zz(1).z = (zz(1).z - 1)*dy_grid;
	timesTH(1).t = (timesTH(1).t - 1)*dx_grid;
    
    if(ixlim==1 & iylim==1)
            [xinds(1),xinds(2)] = findheight(timesTH(1).t,xlims(1),xlims(2));
            [yinds(1),yinds(2)] = findheight(timesTH(1).t,ylims(1),ylims(2));
            xinds = [xinds(1):xinds(2)];
            yinds = [yinds(1):yinds(2)];
            timesTH(1).t = timesTH(1).t(xinds);
            zz(1).z = zz(1).z(yinds);      
            ixlim=0;
            iylim=0;
    end
            
 

            
    
            
            %time=46;
            savemem=0;
            switch savemem
            case 0
%                Z = WRFUserARW(nc,'Z',time);


switch file_type
    case 'wrfout'
        u = 0.5* (nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) ); %2d wind at one height
        v = 0.5* (nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) ); %2d wind at one height
    case 'met_em'
        u = 0.5* (nc{'UU'}(time,ih_wrf,:,1:end-1) + nc{'UU'}(time,ih_wrf,:,2:end) ); %2d wind at one height
        v = 0.5* (nc{'VV'}(time,ih_wrf,1:end-1,:) + nc{'VV'}(time,ih_wrf,2:end,:) ); %2d wind at one height
end
                
%                 for ilat=1:length(zz(1).z)
%                     for ilon=1:length(timesTH(1).t)
%                         ih_wrf=findheight(Z(:,ilat,ilon),h_wrf);
%                         if length(ih_wrf)==0
%                             ih_wrf=1;
%                         end
%                         pdat(1).p(ilat,ilon) = sqrt( u(ih_wrf,ilat,ilon)^2 + u(ih_wrf,ilat,ilon)^2 );
%                     end
%                 end

                pdat(1).p=sqrt(u(yinds,xinds).^2 + v(yinds,xinds).^2);
             %   pdat(1).p=squeeze(WRFUserARW(nc,'Z',time,0,0,ih_wrf));  %horiz_slice at first level

                
            case 1
                
                for ilat=1:length(zz(1).z)
                    for ilon=1:length(timesTH(1).t)
                        Z = WRFUserARW(nc,'Z',time,ilat,ilon);
                        ih_wrf=findheight(Z,h_wrf);
                        if length(ih_wrf)==0
                            ih_wrf=1;
                        end
                        pdat(1).p(ilat,ilon) = sqrt( (nc{'U'}(time,ih_wrf,ilat,ilon))^2 + (nc{'V'}(time,ih_wrf,ilat,ilon))^2 );
                    end
                end
                
            end
            

           
ih_wrf_quiver = ih_wrf;
%ih_wrf_quiver = 4;

%u_quiver=nc{'U'}(time,ih_wrf_quiver,:,:);
%v_quiver=nc{'V'}(time,ih_wrf_quiver,:,:);
u_quiver = u(yinds,xinds);
v_quiver = v(yinds,xinds); %is yinds,xinds as is lat,lon
x_quiver=timesTH(1).t;
y_quiver=zz(1).z;

%nx_quiver=75; %25; %number of arrows to draw for x
%ny_quiver=75; %25; %number for y
nx_quiver=25; %25; %number of arrows to draw for x
ny_quiver=25; %25; %number for y

scale_speed_quiver = [15 15]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed


            
%            ih=9;
%            ih2=17;
            

            
            
        case 'vap_3d_vert'
            pdat(1).p=slice(izmin:izmax,:);    
            timesTH(1).t=GridDan(idir).Y1(2:end-1)'/1000;
            %  pdat(1).p=slice(izmin:izmax,:);    
            
            
        case 'ecmwf_surf'
            zz(1).z=(ecmwf(1).lat + 0.62)*1000; %fix to get round fact that conversion done to km and ground height added
            timesTH(1).t = ecmwf(1).lon'-360;
            
            pdat(1).p=squeeze(ecmwf(1).q(1,:,:))*1000;
            pdat(1).p=squeeze(ecmwf(1).t(1,:,:))-273.15;
            
            
            ih=9;
            ih2=17;
            
            sq=size(ecmwf(1).q);
            
            p=repmat(flipud(ecmwf(1).p),[1 25 27]);
            t=squeeze(ecmwf(1).t(it,:,:,:));
            rho=p*100.*28.97e-3/8.3144./t;
            dz=repmat(diff(z)',[1 sq(3) sq(4)]);
            meecm=sum(dz(ih:ih2,:,:).*rho(ih:ih2,:,:).*squeeze(ecmwf(1).q(it,ih:ih2,:,:)));
            air=sum(dz(ih:ih2,:,:).*rho(ih:ih2,:,:) );
            
            meecm2=meecm./air;
            
            
            meecmT=sum( dz(ih:ih2,:,:).*rho(ih:ih2,:,:).*squeeze(ecmwf(1).t(it,ih:ih2,:,:)) ) ./ air;
            
            pdat(1).p=squeeze(meecm2)*1000;
            
            %  pdat(1).p=squeeze(meecmT)-273.15;
            
        case 'cdensity'
            m=TwoD.Q(izmin:izmax,xinds,5);
            vol=TwoD.Q(izmin:izmax,xinds,11); %total volume of graupel / kg air
            m(vol<6e-10)=NaN; % to avoid dividing by very small volumes
            
            pdat(1).p=pi/6 * m./vol;
            
        case 'w_3d'
            pdat(1).p=squeeze(ThreeD.W(2,:,izmin:izmax))';
            zz(1).z=GridDan(1).Y1(2:end-1) / 1000;
            timesTH(1).t = GridDan(1).X1 / 1000;
            
            
        case 'vap_3d'
            
            wrap_slice
            
            [sx sy]=size(slice);
            
            
            
            D=40e3; %for vapour cross section at dump 8
            D=150e3;
%            D=75e3;
            X=D/diff(GridDan(idir).Y1(1:2));
            inds=[sx/2-X:sx/2+X]; 
            
            if min(inds)==0; inds(1)=[]; end
            
            inds=1:size(slice,1);            
            pdat(1).p=slice(inds,:);
           
            
            %     zz(1).z=( (GridDan(1).Y1(1:length(inds)) )/1000 - 0.62) *1000;
            %     timesTH(1).t = GridDan(1).X1(1:length(inds))' / 1000;
            
            
            zz(1).z=( (GridDan(idir).Y1(inds) )/1000 - 0.62) *1000;
            timesTH(1).t = GridDan(idir).X1(2:end-1)' / 1000;
            
        case 'ARM_radar'
            echo_tops;
            d=[0.1:0.1:20];
            zz(1).z=(xar - add_ground_height)*1000; %fix to get round fact that conversion done to km and ground height added
            timesTH(1).t = yar
            pdat(1).p=tops;   
            
        case 'radar'
            ifilter=0;
            
            afilter=1;
            Lfilter=17; %length of filter in km
            clear diff
            nfilter=Lfilter*1000/diff(Grid.Y1(1:2));
            bfilter=ones([1 nfilter])*1/nfilter;
            
            clear raddat
            %            raddat(1)=TwoD;
            %             for iq=1:11
            %                     Lx=size(TwoD.Q,2);
            %                     haloQ=[TwoD.Q(:,:,iq) TwoD.Q(:,:,iq) TwoD.Q(:,:,iq)];
            %                     haloQ=filter(bfilter,afilter,haloQ,[],2); %data with halo either side for filter averaging of edges
            %                     raddat(1).Q(:,:,iq)=haloQ(:,Lx+1:Lx*2); %average using averaging 
            %                     %window. Filters over 2nd dimension (i.e. horizontal). Trying to simulate the 
            %                     %radar dispersion effect (2 deg beam width). Points at left of domain not filtered as though
            %             end
            
            irecalc=1;
            if irecalc==1
                
                RHO=repmat(GridDan(idir).RHON,[1 size(TwoDDan(idir).Q,2)]);  %using this approximation to rho (instead of 2d pressure and temp
                % makes little difference (< 0.5 dbz when tested for typical field)
                %Using this makes life easier for 3d runs (less mem needed)
                
                %             [r,c,p]=size(TwoD.P);
                %			 PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
                %	 Tempera=TwoD.TH2.*((TwoD.PP)./100000).^0.286;
                % 
                % Calculate density
                %	 RHO=(TwoD.PP./287)./Tempera;
                
                ztot=Radar_new(GridDan(idir),TwoDDan(idir),izmin,izmax,RHO);
                
                % raddat(1).Q(:,:,[5 6 7 8])=0;
                
                ifilter=0;
                if ifilter==1
                    Lx=size(ztot,2);
                    haloQ=[ztot ztot ztot];
                    %                 ztot=filter(bfilter,afilter,haloQ,[],2);
                    filter2d=ones([6 17]);filter2d=filter2d/prod(size(filter2d));
                    ztot=filter2(filter2d,haloQ);
                    ztot=ztot(:,Lx+1:Lx*2);
                    
                    
                end
                
                ztot=10*log10(ztot); 
                
                
                %pdat(1).p=filter(b,a,pdat(1).p);
                
            end
            
            %xinds=[150-37:150 1:37];
            xinds=xinds-1;
            
            
            pdat(1).p=ztot(:,xinds);
            
            
            rads=[30 20 15 10 40];  %was [30 20 15 10] until 8th Jan, 2007
            
            istats=0;       %%%%%%%%%%% set istats here for statistics on n10dbz %%%%%%%%%%%
            if istats==1    
                for irad=1:length(rads)
                    [arad brad]=find(ztot>=rads(irad)); %all points with >=n dBz
                    ihs=unique(arad); %all the different height indices with points > 10 dBz
                    
                    n10dbz(1).n(1:izmax,irad,jj)=0;
                    for iradar=1:length(ihs)
                        n10dbz(1).n(izmin-1+ihs(iradar),irad,jj)=length(find(arad==ihs(iradar))); %number of points > n dBz at each height index with > 10 dBz
                    end
                end
            end
            
            
            
            
            
            
        case 'lnb'       
            %			pdat(1).p=repmat(zz(1).z/1000+add_ground_height,[1 size(lnb2d,2)])-lnb2d(izmin:izmax,:);    
            %LNB_2d; %do new LNB calculations            
            
            %pdat(1).p=lnb2d_tot(izmin:izmax,:);
            
            bins=[GridDan(1).Z(izmin)/1000+0.62:0.05:GridDan(1).Z(izmax)/1000+0.62];
            
            for iLNB=izmin:izmax
                pdat(1).p(iLNB-izmin+1,:)=binner(lnb2d_tot(iLNB,:),bins);
            end
            
            timesTH(1).t=bins(2:end);
            %             
            %             minlnb(j).m(:,jj)=min(lnb2d,[],2);
            %             maxlnb(j).m(:,jj)=max(lnb2d,[],2);
            %             
            %             zref=repmat(GridDan(1).Z(1:size(lnb2d,1))/1000+add_ground_height,[1 size(lnb2d,2)]);
            %             lnbdiff=lnb2d-zref;
            %             meanlnb_abv(j).m(:,jj)=zref(:,1)+meanselect(lnbdiff,'dat>0'); %calculate the mean only for points where lnb is lower than where air at
            %             meanlnb_bel(j).m(:,jj)=zref(:,1)+meanselect(lnbdiff,'dat<0'); %calculate the mean only for points where lnb is higher than where air at            
            %             
            %             bins(j).b=GridDan(1).Z/1000+add_ground_height;
            %             ipos=find(lnbdiff>=0);
            %             ineg=find(lnbdiff<0);
            %             
            %             lnbtemp=lnb2d;
            %             lnbtemp(ineg)=0;
            %             lnbbins_pos(j).l(:,jj)=binner(lnbtemp,bins(j).b); %put lnbs into bins - positiviely buoyant only
            %             
            %             lnbtemp=lnb2d;
            %             lnbtemp(ipos)=0;
            %             lnbbins_neg(j).l(:,jj)=binner(lnbtemp,bins(j).b); %put lnbs into bins - negatively buoyant only
            
            %[minlnb_vap,maxlnb_vap,meanlnb_abv_vap,meanlnb_bel_vap,lnbbins_pos_vap,lnbbins_neg_vap,bins_vap]=lnb_calcs(lnb2d_vap,GridDan,add_ground_height,jj,j);
            %[minlnb_tot,maxlnb_tot,meanlnb_abv_tot,meanlnb_bel_tot,lnbbins_pos_tot,lnbbins_neg_tot,bins_tot]=lnb_calcs(lnb2d_tot,GridDan,add_ground_height,jj,j);
            
            % if  jj==fnmin %only search for column numbers on first pass 
            %     clear dgcol
            % 	id=1;
            % 	dgcol(j).d(id)=getDGcol('ALL_LW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ALL_SW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ACC_LW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ACC_SW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ALu_LW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ALu_SW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ALd_LW',dgstrDan(jc).dg);
            % 	id=id+1;
            % 	dgcol(j).d(id)=getDGcol('ALd_SW',dgstrDan(jc).dg);
            % 	id=id+1;
            % end
            % 
            % for idg=1:length(dgcol(j).d)
            %         icediagsRAD(j).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol(j).d(idg)); %saved as icediagsALL_a-b from Nov 2005
            % end
        case 'vertvel'
            %            lenY1=length(GridDan(idir).Y1);
            %            xinds_wrap=[lenY1-floor(xinds(end)/2):lenY1 xinds(1):floor(xinds(end)/2) ];
            
            pdat(1).p=TwoDDan(idir).W(izmin:izmax,xinds); %vertical velocity
            
            %vertvel=TwoD(idir).W;
            %exdirB=[direcDan(jc).dir 'results/diags/vertvel/'];
            %exname=strcat(exdirB,'vertvel-',int2str(jj),'.mat');
            %save(exname,'vertvel');
            
            %            getLEM_up_width;
            
            
            
            
            
            
            
            
        case 'ozone'
            pdat(1).p=TwoDDan(idir).Q(izmin:izmax,xinds,15); 
        case 'potemp'
           % tref=repmat(GridDan(idir).THREF(izmin:izmax),[1 length(GridDan(idir).Y1(xinds))]);
           % pdat(1).p=squeeze(sum(TwoD.TH1(izmin:izmax,xinds),3))+tref ...
           %     - (squeeze(sum(TwoD_init.TH1(izmin:izmax,xinds),3)) + tref) ; %potemp
            
            pdat(1).p=TwoDDan(idir).TH2(izmin:izmax,xinds,[1]);
            
        case 'wind'
            pdat(1).p=squeeze(sum(TwoDDan(idir).W(izmin:izmax,:),3)); %vertical velocity
        case 'lowtracer'                
            pdat(1).p=squeeze(sum(TwoD.Q(izmin:izmax,xinds,[10]),3)); %low tracer
        case 'totwater'                
            pdat(1).p=f*squeeze(sum(TwoD.Q(izmin:izmax,xinds,[1:6]),3)); %total water
        case 'vapour'                
            pdat(1).p=f*squeeze(sum(TwoDDan(idir).Q(izmin:izmax,xinds,[1]),3)); %ice MR
        case 'general'                
            %            rho=repmat(GridDan(idir).RHON(izmin:izmax),[1 ix2-ix+1]);
            %          pdat(1).p=1000*squeeze(sum(TwoDDan(idir).Q(izmin:izmax,xinds,[2:6]),3)); %tot condensate
            % pdat(1).p=f*squeeze(sum(TwoDDan(idir).Q(izmin:izmax,xinds,[2:6]),3)); %tot condensate            
            %          pdat(1).p=1000*squeeze(sum(TwoD.Q(izmin:izmax,xinds,[2]),3)); %tot condensate
            
%             tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
%             pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
%             rhoref=pref.*28.97e-3/8.3144./tref;     
%    ****** N.B. the above method using THREF and PREFN is inaccurate - much better to use Grid.RHON as compares better with actual 2D rho values ****

               rhoref = lemrho(TwoDDan(idir),GridDan(idir)); %this calculates the actual 2D density from the temperature and pressure within cloud
              % rhoref=ones(length(GridDan(idir).Z),length(GridDan(idir).Y1) ); % set this to remove the conversion from kg/kg to g/m3
            
            
            
          %  pdat(1).p=1000*squeeze(sum(TwoDDan(idir).Q(izmin:izmax,xinds,:),3)).*rhoref(izmin:izmax,xinds,:); %tot condensate
%            pdat(1).p=1000*squeeze(TwoDDan(idir).Q(izmin:izmax,xinds)); %
            pdat(1).p=1000*squeeze(TwoDDan(idir).Q(izmin:izmax,xinds)).*rhoref(izmin:izmax,xinds,:); %tot condensate per m3
            
%            pdat(1).p=squeeze(TwoDDan(idir).Q(izmin:izmax,xinds)); %number per kg            
%            aone=find(pdat(1).p<=1.1); pdat(1).p(aone)=0; %removes all the 1 per kg values (default value for numbers in LEM)    
%            pdat(1).p=pdat(1).p.*rhoref(izmin:izmax,xinds,:)*1e-3; %convert to per Litre values

%            pdat(1).p=TwoDDan(idir).W(izmin:izmax,xinds,:); 
           
%pdat(1).p=TwoDDan(idir).T(izmin:izmax,xinds,:); 
%pdat(1).p=TwoDDan(idir).P(izmin:izmax,xinds,:); 
            
        case 'inc'                
            %            pdat(1).p=squeeze(sum(TwoDDan(idir).Q(izmin:izmax,:,[7:9]),3))/1e8; %tot condensate
            pdat(1).p=1e-6*squeeze(sum(TwoD.Q(izmin:izmax,xinds,[7:9]),3)); %tot condensate    
            
            
        case 'si'
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            T=squeeze(sum(TwoD(idir).TH1,3))+tref;
            P=TwoD(idir).PP;
            T=T./(1e5./P).^0.286;
            qsi=satvapPress(T,'lem','ice',P,1)/fact; %satvappress gives in ppmv if 5th argument=1
            si=100*(TwoD(idir).Q(:,:,1)-qsi)./qsi;
            pdat(1).p=si(izmin:izmax,xinds);
            
            simaxTimH(1).s(:,jj)=max(si,[],2);
            siminTimH(1).s(:,jj)=min(si,[],2);
            simean(1).s(:,jj)=mean(si,2);
            
        case 'icesatMR'
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            T=squeeze(sum(TwoD(idir).TH1,3))+tref;
            P=TwoD(idir).PP;
            T=T./(1e5./P).^0.286;
            qsi=satvapPress(T,'lem','ice',P,1); %satvappress gives in ppmv if 5th argument=1
            %si=100*(TwoD(idir).Q(:,:,1)-qsi)./qsi;
            pdat(1).p=qsi(izmin:izmax,xinds);
            
            
        case 'temppert'
            
            TwoD=TwoDDan(idir);
            
            %tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]); %ref potemp
            tref=repmat(GridDan(idir).THREF + GridDan(idir).OLTHBAR,[1 length(GridDan(idir).Y1)]); %ref potemp
            
            
            T=squeeze(sum(TwoD.TH1,3))+tref; %potemp
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            P=TwoD.PP; %actual p
            T=T./(1e5./P).^0.286; %actual T
            tref=tref./(1e5./pref).^0.286; %ref temp
            
            Tp=T-tref; %perturbation of temperature
            
            pdat(1).p=Tp(izmin:izmax,xinds);
            
            tpertTimH(1).t(:,jj)=mean(Tp,2); %mean temp pert
            
        case 'rhopert577'
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            
            
            T=TwoDDan(idir).TH1+tref; %tot potemp
            Tav=repmat(mean(T,2),[1 length(GridDan(idir).Y1)]); %mean T at this point in time
            
            P=TwoDDan(idir).PP; %tot P
            Pav=repmat(mean(P,2),[1 length(GridDan(idir).Y1)]); %mean P at this time
            
            T=T./(1e5./P).^0.286; %tot temp
            tref=tref./(1e5./pref).^0.286; %tot temp
            
            Tav=Tav./(1e5./Pav).^0.286; %tot temp
            
            rho=P.*28.97e-3/8.3144./T;
            rhoref=pref.*28.97e-3/8.3144./tref;
            rhoref=Pav.*28.97e-3/8.3144./Tav;
            
            rhoref=rhoref.*(1+TwoDDan(idir).Q(:,:,1))./(1+1.608*TwoDDan(idir).Q(:,:,1));
            rhomoist=rho.*(1+TwoDDan(idir).Q(:,:,1))./(1+1.608*TwoDDan(idir).Q(:,:,1));
            
            rhopert=rhomoist-rhoref;
            
            pdat(1).p=rhopert(izmin:izmax,xinds); %so is the density change as calculated from the mean temp and pressure profiles
            
            %rhopertTimH(1).t(:,jj)=mean(rhopert,2); %mean temp pert
            %rhopertTimHmax(1).t(:,jj)=max(rhopert,[],2); %mean temp pert
            %rhopertTimHmin(1).t(:,jj)=min(rhopert,[],2); %mean temp pert
            
        case 'hydbal'
            %pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref presssure
            % pdat(1).p=TwoDDan(idir).PP(izmin:izmax,:)-pref(izmin:izmax,:);
            %zrefdiff=repmat(diff(GridDan(idir).Z(izmin-1:izmax)),[1 length(GridDan(idir).Y1)]); %height diff
            zref=repmat(GridDan(idir).Z(izmin-1:izmax+1),[1 length(GridDan(idir).Y1)]);
            
            %dpdz=diff(TwoDDan(idir).PP(izmin-1:izmax,:),1,1)./zrefdiff;
            
            dpdz=diffdan(TwoD(idir).PP(izmin-1:izmax+1,:),zref,1);
            
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            T=TwoD(idir).TH1+tref; %tot potemp
            P=TwoD(idir).PP; %tot P
            T=T./(1e5./P).^0.286; %tot temp
            rho=P.*28.97e-3/8.3144./T;
            
            
            rhomoist=rho.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1));
            
            pdat(1).p=-rhomoist(izmin:izmax,:)*9.81 - dpdz; %hydrostatic balance : remaining upwards force
            
        case 'rhog'
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            
            
            T=TwoD(idir).TH1+tref; %tot potemp
            Tav=repmat(mean(T,2),[1 length(GridDan(idir).Y1)]); %mean T at this point in time
            
            P=TwoD(idir).PP; %tot P
            Pav=repmat(mean(P,2),[1 length(GridDan(idir).Y1)]); %mean P at this time
            
            T=T./(1e5./P).^0.286; %tot temp
            tref=tref./(1e5./pref).^0.286; %tot temp
            
            Tav=Tav./(1e5./Pav).^0.286; %tot temp
            
            rho=P.*28.97e-3/8.3144./T;
            rhoref=pref.*28.97e-3/8.3144./tref;
            %rhoref=Pav.*28.97e-3/8.3144./Tav;
            
            
            rhomoist=rho.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1));
            
            rhopert=rhomoist-rhoref;
            
            pdat(1).p=9.81*rhopert(izmin:izmax,:);
            
            %            rhopertTimH(1).t(:,jj)=mean(rhopert,2); %mean temp pert
            
            
        case 'dpdz'
            %pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref presssure
            % pdat(1).p=TwoDDan(idir).PP(izmin:izmax,:)-pref(izmin:izmax,:);
            %zrefdiff=repmat(diff(GridDan(idir).Z(izmin-1:izmax)),[1 length(GridDan(idir).Y1)]); %height diff
            zref=repmat(GridDan(idir).Z(izmin-1:izmax+1),[1 length(GridDan(idir).Y1)]);
            
            %dpdz=diff(TwoDDan(idir).PP(izmin-1:izmax,:),1,1)./zrefdiff;
            
            dpdz=diffdan(TwoD(idir).PP(izmin-1:izmax+1,:),zref,1);
            
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            T=TwoD(idir).TH1+tref; %tot potemp
            P=TwoD(idir).PP; %tot P
            T=T./(1e5./P).^0.286; %tot temp
            rho=P.*28.97e-3/8.3144./T;
            
            dpdz=diffdan(TwoD(idir).PP(izmin-1:izmax+1,:),zref,1);
            
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            refdpdz=diffdan(pref(izmin-1:izmax+1,:),zref,1);
            
            pdat(1).p= - dpdz + refdpdz; %pressure grad perturbation
            
            
            
        case 'rad'
            pdat(1).p=TwoDDan(idir).twdrad(izmin:izmax,:);
%%%            
end  %switch i577  
%%%        
    case 56
        top_down_cumulative
    case 55
        upflux_7thSep2005
    case 54
        fallspeed_7thSep2005
        
    case 53
        meanIce_7thSep2005
        
    case 52
        microrate_7thsep2005
    case 51
        start=fact*repmat(sum(icediag4(idir).i(izmin:izmax,3,35:36),3),[1 length(dumprange)]);
        pdat(i).p=fact*squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,35:36),3)) - start;
        
    case 50
        %pdat(i).p=fact*squeeze(tot_prctiles(idir).t(izmin:izmax,dumprange,1));
        %pdat(i).p=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,37),3); %vapour
        %       pdat(i).p=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,37:42),3); %mean total water - need to divide diags by NPES for 2-D multi-processor run
        
        %pdat(i).p=fact*squeeze(vap_prctiles(idir).t(izmin:izmax,dumprange,1));
        
        %       pdat(i).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,34:36),3); %microphysical number rate
        %pdat(i).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,31:33),3); %microphysical mass rate
        
        %        pdat(1).p=squeeze(w_prctiles(idir).w(prc,izmin:izmax,dumprange));
        
        %		pdat(i).p=zmax(idir).z(izmin:izmax,dumprange); %mean total water - need to divide diags by NPES for 2-D multi-processor run
        
        
        %  pdat(1).p=f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[22 23 24]),3),GridDan(idir).t,izmin-1,izmax); %fall speed flux
        
        %  pdat(1).p=f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[1:6]),3),GridDan(idir).t,izmin-1,izmax); %dq/dt from wq flux
        
        
        if iflux==1
            ad_calcs4timeseries;
        end
        
        
        
        switch s50
        case 'dT_conv' 
            pdat(1).p= dT_conv(1).dat(izmin:izmax,dumprange);
            %    pdat(1).p= dT_nonconv(1).dat(izmin:izmax,dumprange);
            %	pdat(1).p= dT_bubble(1).dat(izmin:izmax,dumprange);
            
            %distance covered by contour (km)
            
        case 'lwc_width' 
            pdat(1).p=squeeze( lwc_width(idir).dat(izmin:izmax,1,dumprange) ) * diff(GridDan(idir).Y1(1:2))/1000; %radar dbzs of 30 20 15 10 40
            %distance covered by contour (km)
            
        case 'mass flux'
            
            flux_case='ALL_WTH';
            flux_case='vap_neg_up';
            flux_case='vap_neg_tot';
            flux_case='ALu_W';
            
            switch flux_case
                
            case 'ALL_WQtot'
                
                area=icediagsALL(idir).i(izmin:izmax,dumprange,[283])/npess2(idir); %W>1_A
                pdat(1).p=area*diff(GridDan(idir).Y1([1 end])) /1000;
                area(find(area==0))=1;
                
                % pdat(1).p=icediagsALL(idir).i(izmin:izmax,dumprange,353)./area /npess2(idir) ; %W>1_W
                
            case 'ALL_WTH'
                
                pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[353 354]),3) /npess2(idir) ; %WTHAD + WTHSG
                pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[353 354]),3) /npess2(idir) ; %WTHAD + WTHSG
                
                tit(1).tit='WTH';
                
            case 'ALL_WQtot'
                
                %        pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[355 356]),3) /npess2(idir) ; %VWAD + VWSG
                tit(1).tit='VW';
                %        pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[357 358]),3) /npess2(idir) ; %WKE + WKESG
                tit(1).tit='WKE';
                
            case 'ALu_W'
                
                %        pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[359:362]),3) /npess2(idir) ; %VV + WW
                %       tit(1).tit='V''V'' + W''W''';
                
                %        pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[359:360]),3) /npess2(idir) ; %VV + WW
                %        tit(1).tit=' V''V'' ';
                
                %  pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[363 373]),3) /npess2(idir) ; %W'Qv' 363(norm) 373(sub)
                %  tit(1).tit=' W''Qv'' ';
                
                %       area=icediagsALL(idir).i(izmin:izmax,dumprange,[280]) / npess2(idir); %ALu_A
                %	   area(find(area==0))=1;
                
                pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[302]),3) /npess2(idir) ; %ALu_W
                pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[285]),3) /npess2(idir) ; %ACu_W
                tit(1).tit='ACu_W ';
                
            case 'ALL_WQtot'
                
                %         pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[46]),3) /npess2(idir) ; %ALu_WQ01
                %         tit(1).tit=' ALu_WQ01 ';
            case 'ALL_WQtot'
                
                %         pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[49:51]),3) /npess2(idir) ; %ALu_WQtot
                %         tit(1).tit=' ALu_WQice ';
            case 'ALL_WQtot'       
                %        pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[302]),3) /npess2(idir) ; %ACu_W
                tit(1).tit=' ACu_W ';
                
            case 'vap_neg_up'
                pdat(1).p=vapneg_up(idir).dat(izmin:izmax,dumprange) *diff(GridDan(idir).Y1(1:2))/1000 ; %ACu_W
                tit(1).tit=' Vap neg up ';
                
            case 'vap_neg_d'          
                pdat(1).p=vapneg_d(idir).dat(izmin:izmax,dumprange) *diff(GridDan(idir).Y1(1:2))/1000 ; %ACu_W
                tit(1).tit=' Vap neg down ';
                
            case 'vap_neg_tot'          
                pdat(1).p= ( vapneg_up(idir).dat(izmin:izmax,dumprange) + ...
                    vapneg_d(idir).dat(izmin:izmax,dumprange)  )       *diff(GridDan(idir).Y1(1:2))/1000 ; %ACu_W
                tit(1).tit=' Vap neg total ';   
                
            case 'ALL_WQtot'
                
                %       pdat(1).p=vapneg_d(idir).dat(izmin:izmax,dumprange)  + vapneg_up(idir).dat(izmin:izmax,dumprange) ; %ACu_W
                %       tit(1).tit=' Vap neg total ';
                
                area=icediagsALL(idir).i(izmin:izmax,dumprange,[284]) / npess2(idir); %ALu_A
                area(find(area==0))=1;                    
                
            case 'ALL_WQtot'
                
                %   pdat(1).p=( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[247]),3)./area - icediagsALL(idir).i(izmin:izmax,dumprange,246) ) /npess2(idir); %ACu_TH
                tit(1).tit=' ACC_TH ';       
                
            case 'ALL_WQtot'
                pdat(1).p=( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[1:6]),3) ) /npess2(idir); %ALL_WQtot
                tit(1).tit=' ALL_WQtot ';
                
                
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%            
%savename=tit(1).tit;
%%%%%%%%%%%%%%%%%%%%%%%%%            
            
        case 'mean_rho'
            pdat(1).p=rho_prof(idir).tot(izmin:izmax,dumprange)...
                - repmat(rho_prof(idir).mean(izmin:izmax,1) , [1 length(dumprange)]);    
        case 'radar10dbz' 
            pdat(1).p=squeeze( n10dbz(idir).n(izmin:izmax,1,dumprange) ) * diff(GridDan(idir).Y1(1:2))/1000; %radar dbzs of 30 20 15 10 40
            %distance covered by contour (km)
            
        case 'radar_ndbzARM' 
            pdat(1).p=squeeze( ntop(idir).n(:,2,1:end) ); %radar dbzs of 5:5:70    
            zz(1).z=(zar - add_ground_height)*1000; %fix to get round fact that conversion done to km and ground height added
            timesTH(1).t = [0:10/60:10/60*(size(ntop(1).n,3)-1)];
            
        case 'meanTemp' 
            area=icediagsALL(idir).i(izmin:izmax,dumprange,[284]); %280
            area(find(area==0))=1;
            
            P=repmat(GridDan(idir).PREFN(izmin:izmax),[1 length(dumprange)]);
            T=icediagsALL(idir).i(izmin:izmax,dumprange,247)./area./(1e5./P).^0.286; %248
            
            area=1;
            T0=repmat(icediagsALL(idir).i(izmin:izmax,1,246)./area./(1e5./P(:,1)).^0.286 ,...
                [1 length(dumprange)]);
            T=icediagsALL(idir).i(izmin:izmax,dumprange,246)./area./(1e5./P).^0.286 - T0; %246=ALL_TH
            
            % T=icediagsALL(idir).i(izmin:izmax,dumprange,246)./(1e5./P).^0.286; %248
            pdat(1).p=T;
            
        case 'totdist'
            d=[0.1:0.1:20];
            zz(1).z=(d - add_ground_height)*1000; %fix to get round fact that conversion done to km and ground height added
            
            pdat(1).p=totdist(idir).v(1:200,dumprange); 
            
        case 'vapdist'
            d=[0.1:0.1:20];
            zz(1).z=(d - add_ground_height)*1000; %fix to get round fact that conversion done to km and ground height added
            
            pdat(1).p=vapdist(idir).v(:,dumprange); 
            
        case 'minTguess'
            pref=repmat(GridDan(idir).PREFN(izmin:izmax),[1 length(dumprange)]); %ref p            
            rhoref=repmat(GridDan(idir).RHON(izmin:izmax),[1 length(dumprange)]); %ref p            
            
            rhopert=rhopertTimHmax(idir).t(izmin:izmax,dumprange);
            rhopert=rhopertTimHmin(idir).t(izmin:izmax,dumprange);
            
            rho=rhopert+rhoref;
            tref=pref.*28.97e-3/8.3144./rhoref;
            
            pdat(1).p=pref.*28.97e-3/8.3144./rho - 273.15; 
            
            
        case 'maxlowtracer'
            pdat(1).p=low_prctiles(idir).t(izmin:izmax,dumprange,21); 
            masspr_i=10;
            tracer_i=10;
            nums=tra(idir).num(izmin:izmax,dumprange,tracer_i,masspr_i);
            izero=find(nums==0);
            nums(izero)=1;
            pdat(1).p=tra(idir).mass(izmin:izmax,dumprange,tracer_i,masspr_i)./nums; 
            
        case 'tracerflux'
            pdat(1).p=icediagsALL(idir).i(izmin:izmax,dumprange,151)/npess2(idir);
            
            
        case 'icedist'     
            
            %run gamdistTimH first - makes ice dist from the icediagsALL averages rather than loading in 2d fields
            gamdistTimH 
            
            
            iend=2800;
            iend=2500;
            %iend=3500;
            %iend=3400;
            d=[D(1):D(iend)/500:D(iend)];
            sum_dm=0;
            for it=1:1
                sum_dm=sum_dm+distIce(it).dm(:,dumprange);
            end
            pdat(i).p=1e-6*f*interp1(D,sum_dm,d);
            
            zz(1).z=(d'*1e6 - add_ground_height)*1000; %fix to get round fact that conversion done to km and ground height added
        case 'icendist'
            iend=2800;
            iend=3500;
            iend=3400;
            d=[D(1):D(iend)/500:D(iend)];
            sum_dn=0;
            for it=1:3
                sum_dn=sum_dn+distIce(it).dn(:,dumprange);
            end
            pdat(i).p=1e-6*interp1(D,sum_dn,d);
            
            zz(1).z=(d'*1e6 - add_ground_height)*1000; %fix to get round fact that conversion done to km and ground height added
            
        case 'fallflux'       
            pdat(i).p=-fall_from_mean;
            %		pdat(i).p=tot_fallflux(izmin:izmax,dumprange);
            
        case 'lnbdist_tot' 
            ilnbmin=findheight(bins_tot(idir).b,minZ/1000+0.62);
            ilnbmax=findheight(bins_tot(idir).b,maxZ/1000+0.62);
            
            zz(1).z=1000*(bins_tot(idir).b(ilnbmin+1:ilnbmax)+bins_tot(idir).b(ilnbmin:ilnbmax-1))/2 - 1000*add_ground_height; %since for plotting adds ground height and /1000
            pdat(i).p=lnbbins_neg_tot(idir).l(ilnbmin:ilnbmax-1,dumprange);
            
        case 'lnbdist_vap' 
            ilnbmin=findheight(bins_tot(idir).b,minZ/1000+0.62);
            ilnbmax=findheight(bins_tot(idir).b,maxZ/1000+0.62);
            
            zz(1).z=1000*(bins_vap(idir).b(ilnbmin+1:ilnbmax)+bins_vap(idir).b(ilnbmin:ilnbmax-1))/2 - 1000*add_ground_height; %since for plotting adds ground height and /1000
            pdat(i).p=lnbbins_neg_vap(idir).l(ilnbmin:ilnbmax-1,dumprange);    
            
        case 'rad'       
            pdat(i).p=sum(icediagsRAD(idir).i(izmin:izmax,dumprange,[1 2]),3); 
        case 'swrad'       
            pdat(i).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[2]),3); 
            
            clear diff
            
            pdat(i).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[239]),3) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; 
            
            
        case 'lwrad'       
            %pdat(i).p=sum(icediagsRAD(idir).i(izmin:izmax,dumprange,[1]),3); 
            clear diff
            
            pdat(i).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[238]),3) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; 
            
        case 'lnbdist' 
            ilnbmin=findheight(bins_tot(idir).b,minZ/1000+0.62);
            ilnbmax=findheight(bins_tot(idir).b,maxZ/1000+0.62);
            
            zz(1).z=1000*(bins_tot(idir).b(ilnbmin+1:ilnbmax)+bins_tot(idir).b(ilnbmin:ilnbmax-1))/2 - 1000*add_ground_height; %since for plotting adds ground height and /1000
            pdat(i).p=lnbbins_neg_tot(idir).l(ilnbmin:ilnbmax-1,dumprange);
            
        case 'lnbbel' 
            zref=repmat(GridDan(1).Z(1:size(meanlnb_bel_tot(idir).m,1))/1000+add_ground_height,[1 size(meanlnb_bel_tot(idir).m,2)]);
            %pdat(i).p=meanlnb_abv(idir).m(izmin:izmax,dumprange)-zref(izmin:izmax,dumprange);
            pdat(i).p=meanlnb_bel_tot(idir).m(izmin:izmax,dumprange);
            % pdat(i).p=minlnb(idir).m(izmin:izmax,dumprange);
            % pdat(i).p=minlnb_tot(idir).m(izmin:izmax,dumprange);
        case 'dqtot'       
            pdat(i).p=length(GridDan(idir).Y1)*( dq_tot(idir).d(izmin:izmax,dumprange,2) ) *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
            %was divided by length of .Y1 in Allimp...
        case 'meanLT5tot'
            nntot=nn(idir).n(izmin:izmax,dumprange,2);
            pdat(i).p= 5 - length(GridDan(idir).Y1)*dq_tot(idir).d(izmin:izmax,dumprange,2)./nntot ; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
            %was divided by length of .Y1 in Allimp... 
            ilt100=find(nntot<150);
            pdat(i).p(ilt100)=NaN;
            %        pdat(i).p=dq_tot(idir).d(izmin:izmax,dumprange,2);
        case 'nntot'       
            pdat(i).p=nn(idir).n(izmin:izmax,dumprange,1)*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
        case 'nnvap'       
            pdat(i).p=nn2(idir).n(izmin:izmax,dumprange,2)*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
            
            if idir==1
                pdat(i).p=sqrt(4/pi*nn2(idir).n(izmin:izmax,dumprange,2)*dy^2); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
            else
                pdat(i).p=nn2(idir).n(izmin:izmax,dumprange,2)*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
            end
            
        case 'dqvap'       
            pdat(i).p=length(GridDan(idir).Y1) * ( dq_vaps(idir).d(izmin:izmax,dumprange,2) ) *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
            %dq_vaps originally in ppmv * length(Grid.Y1)
            %      dm=repmat( diff(GridDan(idir).Z(izmin-1:izmax)) .* GridDan(idir).RHON(izmin:izmax) , [1 length(dumprange)]);
            %     dm=repmat( GridDan(idir).RHON(izmin:izmax) , [1 length(dumprange)]);
            
            %	pdat(i).p=length(GridDan(idir).Y1) * ( dq_vaps(idir).d(izmin:izmax,dumprange,2) ) .*dm*dy*1000 / f; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
            
            %     ipps=1;
            %	pdat(i).p= ( length(GridDan(idir).Y1)*( dq_vaps(idir).d(izmin:izmax,dumprange,ipps) )...
            %        ./nn2(idir).n(izmin:izmax,dumprange,ipps) ); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more        
            
            
            if idir==1
                %			pdat(i).p=sqrt(length(GridDan(idir).Y1) *4/pi*dq_vaps(idir).d(izmin:izmax,dumprange,2)*dy^2); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
                pdat(i).p=length(GridDan(idir).Y1)*dq_vaps(idir).d(izmin:izmax,dumprange,2)*dy^2; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
                
            else
                pdat(i).p=length(GridDan(idir).Y1) *(dq_vaps(idir).d(izmin:izmax,dumprange,2))*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
                pdat(i).p= pi*(length(GridDan(idir).Y1)*dq_vaps(idir).d(izmin:izmax,dumprange,2)).^2 * dy^2; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
            end
            
            
        case 'dqvap_dist' 
            bins=[0:0.1:20];
            bins=(bins(2:end)+bins(1:end-1))/2;
            
            ih1=findheight(GridDan(idir).Z/1000+0.62,16);
            ih2=findheight(GridDan(idir).Z/1000+0.62,19);
            
            
            for it=dumprange
                
                for k=1:size(vapdist(1).v,3)
                    
                    [avap bvap]=max(vapdist(1).v(:,it,k),[],1);
                    if bvap==51
                        pdat(i).p(k,it-dumprange(1)+1)=-sum((bins(1:bvap-1)-bins(bvap)).*squeeze(vapdist(1).v(1:bvap-1,it,k))' );
                    else
                        pdat(i).p(k,it-dumprange(1)+1)=0;
                    end
                    
                end        
                
            end
            
            ih1=findheight(GridDan(idir).Z/1000+0.62,12);
            ih2=findheight(GridDan(idir).Z/1000+0.62,19);
            zz(1).z=( GridDan(1).Z(ih1:ih2) ) ;
            
        case 'dqvap_dist_abv' 
            bins=[0:0.1:20];
            bins=(bins(2:end)+bins(1:end-1))/2;
            
            for it=dumprange
                
                for k=1:size(vapdist(idir).v,3)
                    
                    [avap bvap]=max(vapdist(idir).v(:,it,k),[],1);
                    if bvap==51
                        bvap=51;
                        pdat(i).p(k,it-dumprange(1)+1)=sum((bins(bvap+1:end)-bins(bvap)).*squeeze(vapdist(idir).v(bvap+1:end,it,k))' );
                    else
                        pdat(i).p(k,it-dumprange(1)+1)=0;
                    end
                    
                end
                
            end
            ih1=findheight(GridDan(idir).Z/1000+0.62,12);
            ih2=findheight(GridDan(idir).Z/1000+0.62,19);
            zz(1).z=( GridDan(1).Z(ih1:ih2) ) ;
            
            
        case 'dqtot_dist_abv' 
            bins=[0:0.1:200];
            bins=(bins(2:end)+bins(1:end-1))/2;
            
            for it=dumprange
                
                for k=1:size(totdist(idir).v,3)
                    
                    [avap bvap]=max(totdist(idir).v(:,it,k),[],1);
                    if bvap==51
                        bvap=51;
                        pdat(i).p(k,it-dumprange(1)+1)=sum((bins(bvap+1:end)-bins(bvap)).*squeeze(totdist(idir).v(bvap+1:end,it,k))' );
                    else
                        pdat(i).p(k,it-dumprange(1)+1)=0;
                    end
                    
                end
                
            end
            ih1=findheight(GridDan(idir).Z/1000+0.62,12);
            ih2=findheight(GridDan(idir).Z/1000+0.62,19);
            zz(1).z=( GridDan(1).Z(ih1:ih2) ) ;
            
            
            
        case 'rhopert'
            pdat(1).p=rhopertTimH(idir).t(izmin:izmax,dumprange); 
            pdat(1).p=rho5ppmv_tot(idir).r(izmin:izmax,dumprange); 
            %		pdat(1).p=rho5ppmv_totpos(idir).r(izmin:izmax,dumprange); 
            %pdat(1).p=low_prctiles(idir).t(izmin:izmax,dumprange,21); 
            %    pdat(1).p=rho5ppmv_vap(idir).r(izmin:izmax,dumprange); 
            
        case 'drhodz'
            
            tref=repmat(GridDan(idir).THREF,[1 length(dumprange)]);
            pref=repmat(GridDan(idir).PREFN,[1 length(dumprange)]); %ref p
            tref=tref./(1e5./pref).^0.286; %tot temp
            rhoref=pref.*28.97e-3/8.3144./tref;
            rho=rhoref+rhopertTimH(idir).t;
            
            dz=GridDan(1).Z(izmin:izmax)-GridDan(1).Z(izmin-1:izmax-1);
            dz=repmat(dz,[1 length(dumprange)]);
            %             
            pdat(1).p=( rho(izmin:izmax,dumprange)-rho(izmin-1:izmax-1,dumprange) ) ./ dz;
            
        case 'upflux'
            pdat(1).p=rho.*icediagsALL(idir).i(izmin:izmax,dumprange,137);
        case 'meanw'
            pdat(1).p=icediagsALL(idir).i(izmin:izmax,dumprange,137);
            
            ihm2=302; %137=ALu_W 302=Acu_W  %303=W>1_W  
            aind=285; %280=ALu_A 285=ACu_A  %283=W>1_A
            
            area=icediagsALL(idir).i(izmin:izmax,dumprange,aind)/npess2(idir);
            ilow=find(area<1e-5);
            %area(area==0)=1;
            area(area==0)=1e99;
            
            pdat(1).p=icediagsALL(idir).i(izmin:izmax,dumprange,ihm2)./area;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors
            
            
            
            
        case 'si'
            %       pdat(1).p=simaxTimH(idir).s(izmin:izmax,dumprange);
            pdat(1).p=simean(idir).s(izmin:izmax,dumprange);
            
        case 'si_diag'
            %       pdat(1).p=simaxTimH(idir).s(izmin:izmax,dumprange);
            pot = icediagsALL(1).i(izmin:izmax,dumprange,246);
            pbar = repmat(GridDan(idir).PREFN(izmin:izmax),[1 size(pot,2)]);
            
            T=pot./(1e5./pbar).^0.286;
            qsi=satvapPress(T,'lem','ice',pbar,1)/f; %satvappress gives in ppmv if 5th argument=1
            
            f=1e6*28.97/18;        
            vap=icediagsALL(idir).i(izmin:izmax,dumprange,37);        
            
            pdat(1).p=100*(vap-qsi)./qsi;
            
        case 'vapadcum'    
            pdat(1).p=-vapadcum;  %vapadcum is the advective loss 
        case 'totadcum'    
            pdat(1).p=ad;    
        case 'icemicrocum'    
            pdat(1).p=cumsum(microicerate,2)*300;    
        case 'icefallcum'    
            pdat(1).p=-cumsum(fallrate,2)*300;
        case 'iceadcum'
            pdat(1).p=iceadcum;     
        case 'low_tracer'
            
            pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[151]),3)/npess2(idir);
            
            Y=diff(GridDan(idir).Y1([1 end]));
            if idir==1
                fact=Y^2; %mulitplying by area covered in 3d and then by area in 2d if assume distance of 3rd dim.
            else
                fact=Y*1000; %1 km in 3rd dimension for 2d
            end
            
            fact=1;
            
            fact=fact/1e6;
            
            pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[151]),3)/npess2(idir);
            
            % pdat(1).p=adlow(izmin:izmax,dumprange);
            
        case 'combined_potemp'
            
            pdat(1).p=-(dqpo(idir).d(izmin:izmax,dumprange)/f) + -(dqnon(idir).d(izmin:izmax,dumprange)/f);
            
            
        case 'change_conv_potemp'
            dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
            apot=abs(-length(GridDan(idir).Y1)*dqnon(idir).d(izmin:izmax,dumprange)*dy/f);
            bpot=abs(-length(GridDan(idir).Y1)*dqpo(idir).d(izmin:izmax,dumprange)*dy/f);
            ratio=apot./(apot+bpot);
            pdat(1).p=change.*ratio;
            
        case 'ratio_potemp'
            dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
            apot=abs(-length(GridDan(idir).Y1)*dqnon(idir).d(izmin:izmax,dumprange)*dy/f);
            bpot=abs(-length(GridDan(idir).Y1)*dqpo(idir).d(izmin:izmax,dumprange)*dy/f);
            pdat(1).p=apot./(apot+bpot);
        case 'dq_non'
            dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
            pdat(1).p=-length(GridDan(idir).Y1)*dqnon(idir).d(izmin:izmax,dumprange)*dy/f;
        case 'dq_potemp'
            dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
            pdat(1).p=-length(GridDan(idir).Y1)*dqpo(idir).d(izmin:izmax,dumprange)*dy/f;
        case 'pcond'
            pdat(1).p=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[29]),3) - f*sum(icediag(idir).i(izmin:izmax,dumprange,liqsink),3);
        case 'dql'
            pdat(1).p=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[29]),3);
        case 'prevp'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[26]),3); 
        case 'pgsub'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[24]),3); 
        case 'dqi'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[28]),3); 
        case 'pisub'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[27]),3); 
        case 'pIsub'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[24 25 27]),3);       
        case 'picesubcum'
            pdat(1).p=cumsum(f*sum(icediag(idir).i(izmin:izmax,dumprange,[24 25 27]),3),2);   
            %  pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[24 25 27]),3);   
        case 'pidep'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[31]),3); 
        case 'pIdep'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[31 9 1]),3);       
        case 'piacw'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[32]),3);  
        case 'allpr'
            'check using correct icediag'
            %pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[iallpr]),3);
            pdat(1).p=f*sum(icediag_nums(idir).i(izmin:izmax,dumprange,[iallpr]),3);
        case 'pifrw'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[34]),3);   
        case 'praut'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[3]),3);   
        case 'racw'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[5]),3);  
        case 'mphys_process'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,imphys),3);       
        case 'PGMLT'
            pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[2]),3);       
        case 'minvap'
            pdat(1).p=f*vap_prctiles(idir).t(izmin:izmax,dumprange,1); %change last index to 1 for min, end for max
        case 'grano'
            pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[44]),3);
        case 'iceno'
            pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43]),3);
            
            rho=repmat(GridDan(1).RHON(izmin:izmax),[1 length(dumprange)]);
            
            area=icediagsALL(idir).i(izmin:izmax,dumprange,[282]);
            area(find(area==0))=1;
            %       pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[38]),3).*rho;
            %       pdat(1).p=icediagsALL(idir).i(izmin:izmax,dumprange,[256]).*rho./area;
            
            %        pdat(1).p=squeeze(q_prctiles.q(end,izmin:izmax,dumprange,7)).*rho;
            
            
            pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43]),3)./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000/1e9;
            
            
        case 'snowno'
            pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[45]),3);
        case 'maxw'
            pdat(1).p=MaxW(idir).w(izmin:izmax,dumprange);
        case 'minw'
            pdat(1).p=MinW(idir).w(izmin:izmax,dumprange);    
        case 'graupel'
            pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[41]),3);
        case 'snow'
            pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40]),3);
        case 'ice'
            % area=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[284]),3)./npess2(idir);
            % area(area==0)=1;
            % pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[226:228]),3)./area./npess2(idir);
            
            % pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[226:228]),3)./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000;
%            pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[42]),3)./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000;


			rho=repmat(GridDan(idir).RHON(izmin:izmax),[1 length(dumprange)]);
            
            ihm2=[232]; %83=ALu_Q02 302=Acu_W   266=W>1_Q02  228=ACC_Q06   291=ACu_Q06  42=ALL_Q06
            aind=[284]; %280=ALu_A 285=ACu_A   283=W>1_A    284=ACC_A      282=ACu_A    []=ALL_A      
            
            if length(aind)==0
                area=1;  %if are using the ALL partition
            else
                area=icediagsALL(idir).i(izmin:izmax,dumprange,[aind])./npess2(idir);
                area(find(area==0))=1;

            end
            
            
            pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[ihm2]),3)./npess2(idir)./area;


            %       pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[42]),3);
        case 'maxice'    
%            pdat(1).p=1000*squeeze(sum(q_prctiles(idir).q(end,izmin:izmax,dumprange,[4 5 6]),4) );
             pdat(1).p=1000*squeeze(sum(q_prctiles(idir).q(end,izmin:izmax,dumprange,[3]),4) );
            
        case 'allice'
            %   pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40:42]),3);
            
            %       pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40:42]),3)./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000;
            
            
            pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40:42]),3)./npess2(idir);
            
            rho=repmat(GridDan(1).RHON(izmin:izmax),[1 length(dumprange)]);
            area=icediagsALL(idir).i(izmin:izmax,dumprange,[285])./npess2(idir);
            area(find(area==0))=1;
            %    pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[289:291]),3)./npess2(idir) .*rho ./area;
            
        case 'liq'
            rho=repmat(GridDan(idir).RHON(izmin:izmax),[1 length(dumprange)]);
            
            ihm2=287; %83=ALu_Q02 302=Acu_W   266=W>1_Q02  224=ACC_Q02   287=ACu_Q02  38=ALL_Q02
            aind=282; %280=ALu_A 285=ACu_A   283=W>1_A    284=ACC_A      282=ACu_A    []=ALL_A      
            
            if length(aind)==0
                area=1;  %if are using the ALL partition
            else
                area=icediagsALL(idir).i(izmin:izmax,dumprange,[aind])./npess2(idir);
                area(find(area==0))=1;

            end
            
%            pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[ihm2]),3)./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000;
            pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[ihm2]),3)./npess2(idir) ./ area;
            
            
            izrem=findheight(GridDan(idir).Z,1e3);
            ixrem=findheight(GridDan(idir).t+3,23);     
            %pdat(1).p(1+izmin-1:izrem+izmin-1,1:ixrem)=0;
            
            
            
            
            
        case 'rain'
            
            ihm2=267; %82=ALu_Q01 302=Acu_W   265=W>1_Q01  223=ACC_Q01  286=ACu_Q01
            aind=283; %280=ALu_A 285=ACu_A   283=W>1_A    284=ACC_A                
            
            area=icediagsALL(idir).i(izmin:izmax,dumprange,[aind])./npess2(idir);
            area(find(area==0))=1;                     
            
            clear diff
            
            %pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[39]),3) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; %39 = ALL_Q03
            pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[ihm2]),3);
            
        case 'minvap' %min vapour
            pdat(i).p=fact*squeeze(vap_prctiles(idir).t(izmin:izmax,dumprange,1)); 
        case 'mintot' %min tot water
            pdat(i).p=fact*squeeze(tot_prctiles(idir).t(izmin:izmax,dumprange,end)); %total water prcs - min =1 max=end
        case 'adrate' %rate of change of tot water due to advection as worked out from change and fall flux 
            pdat(1).p=fluxrate2;
        case 'fallrate'
            pdat(1).p=-fallrate; %source of ice MR from fall speed flux - made it minus so can be put on the same colour scale as advection
        case 'changevap'            
            pdat(1).p=changevap;
        case 'meanvap'            
            pdat(1).p=f*icediagsALL(idir).i(izmin:izmax,dumprange,[37]) /npes;  
        case 'changeice'                      
            pdat(1).p=changeice;
        case 'change'
            pdat(1).p=change; 
        case 'change_from_dqtot'                    
            ylen=length(GridDan(idir).Y1);  
            init=f*repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[37:42]),3),[1 length(dumprange)])/npes; %inital total water MR (ppmv)
            mean_tot = 5 - ylen/1000 * dq_tot(idir).d(izmin:izmax,dumprange,2); %ppmv
            pdat(1).p = mean_tot - init;
            
            %dq_tot was divided by length of .Y1 in Allimp...
            
        case 'change'
            pdat(1).p=change;   
            
        case 'topdowncum'
            pdat(1).p=topdown;           
        case 'changerate'
            pdat(1).p=changerate; 
        case 'icemass'
            ice=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[42]),3) )/npes;
            icenc=(sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43]),3) )/npes;
            
            %ice=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40:42]),3) )/npes;
            %icenc=(sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43:45]),3) )/npes;
            pdat(1).p=changeice./changenc*1e3;      
            pdat(1).p=ice./icenc*1e3;      
            
        case 'microice'
            pdat(1).p=microicerate * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000;
        case 'vapad'
            % pdat(1).p=-vapad;
            pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[1 10]),3);
        case 'micronc'
            pdat(1).p=micronc;
        case 'fallnc'
            pdat(1).p=fallnc;
        case 'adnc'
            pdat(1).p=adnc;
        case 'changenc'
            pdat(1).p=1e-3*changenc;
        case 'icead'
            pdat(1).p=icead; %advective source of ice mixing ratio
        case 'icead_num'
            pdat(1).p=icead; %advective source of ice mixing ratio    
        case 'fall+ad'; %advective source + fall speed source
            pdat(1).p=icead+fallrate;
        end
        
    case 49
        
        T=TwoD.TH2./(1e5./TwoD.PP).^0.286;
        ei=SatVapPress(T,'goff','ice'); %Pa
        P=GridDan(2).PREFN; %Pa
        
        xdat(6).x=f*0.622*ei./(P-ei);
        
    case 48
        pdat(1).p=abs( squeeze( maxw1km(idir).m(izmin:izmax,dumprange) ) ); 
        
    case 488
        pdat(1).p=wemm; 
        
    case 47
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            zz(1).z=z(izmin:izmax);
            timesTH(1).t=x(x1:x2)'/1000;
            pdat(1).p=squeeze(TT(zi:zi2,x1:x2,46))-273.15; 
            
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            % timesTH(2).t=xmpc(x1:x2)/1000;
            [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
            zz(2).z=zmpc(izminmpc:izmaxmpc);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=squeeze(max(icempc((x1:x2),1,izminmpc:izmaxmpc,:),[],1));
        end 
        
    case 46
        zz(1).z=z(izmin:izmax);
        %tt=findheight(time,time2d)
        pdat(1).p=maxSatTot;
        
    case 45    
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=timHlemMR;
            
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            % timesTH(2).t=xmpc(x1:x2)/1000;
            [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
            zz(2).z=zmpc(izminmpc:izmaxmpc);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=timHmpcMR;
        end
        
    case 44 %mean time height mr plots for LEM and MPC
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=f*1e-3*squeeze(max(icemr(1).i(izmin:izmax,(x1:x2),:)+snowmr(1).i(izmin:izmax,(x1:x2),:)+graupelmr(1).i(izmin:izmax,(x1:x2),:),[],2)); 
            %lem values in g/kg so *1e-3 - *f= conversion for kg/kg to ppmv 
            
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            % timesTH(2).t=xmpc(x1:x2)/1000;
            [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
            zz(2).z=zmpc(izminmpc:izmaxmpc);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=f*1e-3*squeeze(  max ( icemrmpc(x1:x2,1,izminmpc:izmaxmpc,:) ./ rhompc(x1:x2,1,izminmpc:izmaxmpc,:) ,[],1)  )*1e3;
            %think mpc values in g/cm^3 so need to divide by the density - defo in kg/m^3 to give g/kg. then *1e-3 to kg/kg
            %so then multiply by 1e-3 to get g/kg 
        end
        
    case 43
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            x1=1;
            x2=length(x);
            
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=1e-6*squeeze(max(icenc(1).i(izmin:izmax,(x1:x2),dumprange)+snownc(1).i(izmin:izmax,(x1:x2),dumprange)+graupelnc(1).i(izmin:izmax,(x1:x2),dumprange),[],2)); 
            %LEM numbers in kg^-1 so times 1e-6 to give (mg)^-1
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            x1=1;
            x2=length(xmpc);
            % timesTH(2).t=xmpc(x1:x2)/1000;
            [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
            zz(2).z=zmpc(izminmpc:izmaxmpc);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=squeeze(max(icempc((x1:x2),1,izminmpc:izmaxmpc,dumprange) ./ rhompc(x1:x2,1,izminmpc:izmaxmpc,dumprange) ,[],1));
            pdat(2).p=squeeze(max(icempc((x1:x2),1,izminmpc:izmaxmpc,dumprange)  ,[],1));
            
            % pretty sure mpc number conc in cm^-3 so divide by density (kg/m^3) to give (mg)^-1
            % would compare much better if mpc valules were in 10^3 cm^-3.
        end 
        
    case 42
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=squeeze(mean(icenc(1).i(izmin:izmax,(x1:x2),:)+snownc(1).i(izmin:izmax,(x1:x2),:)+graupelnc(1).i(izmin:izmax,(x1:x2),:),2))*1e-6; 
            
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            % timesTH(2).t=xmpc(x1:x2)/1000;
            [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
            zz(2).z=zmpc(izminmpc:izmaxmpc);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=squeeze(mean(icempc((x1:x2),1,izminmpc:izmaxmpc,:),1));
        end 
        
    case 41 %2d plots of total ice mr at time tt for LEM and MPC
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            %if ~exist('rhoLEM'); rhoLEM=pressure(1).p.^(1-0.286).*28.97e-3.*1e5^0.286 ./8.31 ./potemp(1).p; end
            timesTH(1).t=x(x1:x2)'/1000;
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=icemr(1).i(izmin:izmax,(x1:x2),tt)+snowmr(1).i(izmin:izmax,(x1:x2),tt)+graupelmr(1).i(izmin:izmax,(x1:x2),tt); 
            
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            timesTH(2).t=xmpc(x1:x2)/1000;
            izminmpc=findheight(zmpc,minZ);
            zz(2).z=zmpc(izminmpc:end);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=squeeze(icemrmpc((x1:x2),1,izminmpc:end,tt)*1e3)';
        end 
        
    case 40
        switch i
        case 1
            %if ~exist('rhoLEM'); rhoLEM=pressure(1).p.^(1-0.286).*28.97e-3.*1e5^0.286 ./8.31 ./potemp(1).p; end
            timesTH(1).t=x'/1000;
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=icemr(1).i(izmin:izmax,:,tt); 
            
        case 2
            timesTH(2).t=xmpc/1000;
            izminmpc=findheight(zmpc,minZ);
            zz(2).z=zmpc(izminmpc:end);
            %tt=findheight(TimeMPC,time2d);
            pdat(2).p=squeeze(icemrmpc(:,1,izminmpc:end,tt)*1e3)';
        end 
        
    case 39
        switch i
        case 1
            [x1 x2]=findheight(x,200e3,500e3);
            if ~exist('rhoLEM'); rhoLEM=pressure(1).p.^(1-0.286).*28.97e-3.*1e5^0.286 ./8.31 ./potemp(1).p; end
            timesTH(1).t=x(x1:x2)'/1000;
            zz(1).z=z(izmin:izmax);
            %tt=findheight(time,time2d)
            pdat(1).p=icenc(1).i(izmin:izmax,x1:x2,tt)/1e6.*rhoLEM(izmin:izmax,x1:x2,tt); 
            
        case 2
            [x1 x2]=findheight(xmpc,200e3,500e3);
            timesTH(2).t=xmpc(x1:x2)/1000;
            [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
            zz(2).z=zmpc(izminmpc:izmaxmpc);
            %tt=findheight(TimeMPC,time2d);
            
            pdat(2).p=squeeze(icempc(x1:x2,1,izminmpc:izmaxmpc,tt).*rhompc(x1:x2,1,izminmpc:izmaxmpc,tt))';
        end    
    case 38
        pdat(i).p=icediag2(1).i(izmin:izmax,dumprange,12); %see note for description of icediag2
        
    case 37
        %ratio of process rates included in MPC model to those not included for ice MR calc in LEM
        %ptot= sum(icediag(1).i(izmin:izmax,dumprange,[2 16 3 4 17:23]),3);
        
        
        %mpc = sum(icediag(1).i(izmin:izmax,dumprange,[2 16 19 22 23]),3);
        
        smrsour=[8 9 18 11 35 36]; %sources of snow mixing ratio 
        smrsink=[25 14 17 19 6];
        smrmpc=[8 9 25 19 6]; %processes included in the MPC
        
        %brat(brat==0)=1e-10; %ensures that when b=0 a value shows up in the ratio
        
        ptot= sum(icediag(idir).i(izmin:izmax,dumprange,[smrsour smrsink]),3);
        mpc = sum(icediag(idir).i(izmin:izmax,dumprange,[smrmpc]),3);
        
        pdat(i).p=mpc./ptot;
        
    case 36
        %ratio of process rates included in MPC model to those not included for snow NC calc in LEM
        lem= sum( icediag(1).i(izmin:izmax,dumprange,[snowncsour snowncsink] ) ,3);
        
        
        mpc = sum( icediag(1).i(izmin:izmax,dumprange,[sncmpc] ) ,3);
        
        %brat(brat==0)=1e-10; %ensures that when b=0 a value shows up in the ratio
        mpc(mpc<1e-3)=0;
        lem(lem<1e-3)=0;
        
        pdat(i).p=mpc./lem;  %./(lem);
        
    case 35
        %ratio of process rates included in MPC model to those not included for ice NC calc in LEM
        
        icencsour=[59:62]; incmpc=[59 61 62 45]; %RIPRM - heterogeneous nucleation not in MPC, but homog would replace it
        icencsink=[52 53 63 64 45 46];
        
        snowncsour=[39 45 48]; sncmpc=[45 48 58 50 56]; %not sure whether RSBRK(48) is in MPC
        snowncsink=[58 50 56 43 40 47];
        
        
        
        lem=sum( icediag(1).i(izmin:izmax,dumprange,[icencsour icencsink] ) ,3);
        mpc=sum( icediag(1).i(izmin:izmax,dumprange,[incmpc]) ,3);
        
        mpc(mpc<1e-3)=0;
        lem(lem<1e-3)=0;
        
        %brat(brat==0)=1e-10; %ensures that when b=0 a value shows up in the ratio
        
        pdat(i).p=mpc./lem;
        
        
        
        
    case 34
        %ratio of process rates included in MPC model to those not included for ice MR calc in LEM
        %mpc=( sum(icediag(1).i(izmin:izmax,dumprange,[10 8]),3) + sum(icediag(1).i(izmin:izmax,dumprange,[2 15 1]),3) );
        %non=( sum(icediag(1).i(izmin:izmax,dumprange,[7 9 11]),3) + sum(icediag(1).i(izmin:izmax,dumprange,[5 3 12 6 4]),3) ) ;
        
        imrsour=[29:34]; %sources of ice mixing ratio as in 'oldnew' process numbering
        imrsink=[8 11 12 27 7 21 36];
        imrnonmpc=[11 12 21 36 29 32]; %accretion type processes not included in the MPC
        imrmpc=[8 27 7 30 31 33 34]; %non-accretion type processes prob accounted for in the MPC
        
        mpc=sum(icediag(idir).i(izmin:izmax,dumprange,[imrmpc]),3);
        %non=sum(icediag(idir).i(izmin:izmax,dumprange,[imrnonmpc]),3);
        non=sum(icediag(idir).i(izmin:izmax,dumprange,[imrsour imrsink]),3);
        
        pdat(i).p=mpc./(non);
        
    case 33
        pdat(i).p=inputted(izmin:izmax,dumprange); %icediag(1).i(izmin:izmax,dumprange,[11]); 
        
    case 32
        pdat(i).p=(wq(izmin:izmax,dumprange));
        offset=abs(minALL(pdat(i).p)); 
        pdat(i).p=pdat(i).p + offset;
        pdat(i).p(pdat(i).p<=0)=NaN;
        %pdat(i).p=-pdat(i).p;
    case 31
        pdat(i).p=sum(icediag(1).i(izmin:izmax,dumprange,[12]),3)/MI0; 
    case 30
        pdat(i).p=sum(icediag(1).i(izmin:izmax,dumprange,[13]),3); 
    case 29
        pdat(i).p=sum(icediag(1).i(izmin:izmax,dumprange,[8]),3)/MI0;  %*ndivqav(izmin:izmax,dumprange);  
        
        
    case 1
        pdat(i).p=squeeze(max(icemr(i).i(izmin:izmax,:,dumprange),[],2));
    case 2
        pdat(i).p=icediag(1).i(izmin:izmax,dumprange,iproc);
    case 3
        pdat(i).p=squeeze(max(snownc(i).i(izmin:izmax,:,dumprange),[],2));
    case 4
        %         pdat(i).p=sumLiq(izmin:izmax,dumprange);
        pdat(i).p=sum(icediag(1).i(izmin:izmax,dumprange,[3]),3)/MI0;
        
    case 5
        pdat(i).p=sum(icediag(1).i(izmin:izmax,dumprange,[5]),3).*ndivqav(izmin:izmax,dumprange);
    case 6
        pdat(i).p=pcents_icemr(i).p(dumprange,izmin:izmax,4)';
    case 7
        pdat(i).p=pcents_icemr(i).p(dumprange,izmin:izmax,3)';
    case 8
        switch i
        case 1
            %pdat(i).p=permute(pcents_potemp(i).p(izmin:izmax,dumprange,8),[2 1 3])';
            [izmin izmax]=findheight(Grid.Z,14e3,22e3);
            timesTH(1).t=Grid.Y1'/1000;
            zz(1).z=Grid.Z(izmin:izmax);
            %pdat(i).p=TwoDDan(1).TH2(izmin:izmax,:);
            T=potemp(i).p(izmin:izmax,:,tt)./(1e5./pressure(i).p(izmin:izmax,:,tt)).^0.286;
            pdat(1).p=SatVapPress(T,'goff','ice',pressure(i).p(izmin:izmax,:,tt),1); %ppmv
            
            izovr=1;
            %pdat(1).p=T(izmin:izmax,:)-273.15;
            logflag=1;
            
            %imaxovr=[1 0];
            maxcovOvr=2;
            
            ncont=20;
            clines=1;
            
            
        case 2
            timesTH(2).t=Grid.Y1'/1000;
            [izmin izmax]=findheight(Grid.Z,14e3,22e3);
            zz(2).z=Grid.Z(izmin:izmax);
            pdat(2).p=icenc(1).i(izmin:izmax,:,tt)/1e6;
            logflag=1;
            
            iminovr=[0 1];
            mincovOvr=0;
            clines=0;
            
        end
        
    case 9
        pdat(i).p=satmr(i).s(izmin:izmax,:,85);
    case 10
        pdat(i).p=squeeze(max(icenc(i).i(izmin:izmax,:,dumprange),[],2));
    case 11
        pdat(i).p=squeeze(max(icemr(i).i(izmin:izmax,:,dumprange),[],2)) * fact/1000 ...
            + squeeze(max(snowmr(i).i(izmin:izmax,:,dumprange),[],2)) * fact/1000 ;
    case 12
        pdat(i).p=squeeze(max(snowmr(i).i(izmin:izmax,:,dumprange),[],2));
    case 13
        pdat(i).p=squeeze(max(graupelmr(i).i(izmin:izmax,:,dumprange),[],2));
    case 14
        vv=fact*repmat(mean(vap(i).v(izmin:izmax,:,dumprange(1)),2),[1 length(dumprange)]);
        pdat(i).p=squeeze(min(fact*vap(i).v(izmin:izmax,:,dumprange),[],2))-vv;
    case 15
        pdat(i).p=fact*squeeze(mean(vap(i).v(izmin:izmax,:,dumprange),2));
    case 16
        vv=fact*repmat(mean(vap(i).v(izmin:izmax,:,dumprange(1)),2),[1 length(dumprange)]);
        pdat(i).p=squeeze(pcents(i).p(dumprange,izmin:izmax,2))'-vv;
    case 17
        pdat(i).p=fact*squeeze(min(vap(i).v(izmin:izmax,:,dumprange),[],2));
    case 18
        pdat(i).p=squeeze(min(satmr(i).s(izmin:izmax,:,dumprange),[],2));
    case 19
        switch i
        case 1
            pdat(i).p=squeeze(min(satmr(1).s(izmin:izmax,:,dumprange),[],2));
        case 2
            pdat(i).p=fact*squeeze(min(vap(1).v(izmin:izmax,:,dumprange),[],2));
        end
        
    case 20
        switch i
        case 1
            pdat(i).p=squeeze(max(icemr(i).i(izmin:izmax,:,dumprange),[],2));
        case 2
            pdat(i).p=fact*squeeze(max(snowmr(1).i(izmin:izmax,:,dumprange),[],2));
        end
        
    case 21
        pdat(i).p=squeeze(max(V(i).v(izmin:izmax,:,dumprange),[],2));
        
    case 22
        %pdat(i).p=squeeze(max(Vsnow(i).v(izmin:izmax,:,dumprange),[],2));
        pdat(i).p=pcents_vsnow(i).p(izmin:izmax,dumprange,6);
        
    case 23
        %pdat(i).p=squeeze(max(Vsnow(i).v(izmin:izmax,:,dumprange),[],2));
        dgfind=findhead('ALL_WQ04',dgstrDan(1).dg)
        pdat(i).p=squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange));
        
    case 24
        pdat(i).p=squeeze(Falldiag(i).dg(izmin:izmax,6,dumprange));
        
    case 25
        dzz=repmat(dz,[1 length(dumprange)]);
        rho=repmat(GridDan(i).RHON(izmin:izmax),[1 length(dumprange)]);
        dgfind=findhead('ACC_A',dgstrDan(1).dg);
        A=squeeze(diag(i).dg(izmin:izmax,dgfind(1),1:length(dumprange)));
        az=find(A<0.001);
        A(az)=1;
        
        
        pdat(i).p=fact*squeeze(Falldiag(i).dg(izmin:izmax,6,dumprange))./dzz./rho./A*300;
        
    case 26
        pdat(i).p=squeeze(Fluxdiag(i).dg(izmin:izmax,6,dumprange) +Fluxdiag(i).dg(izmin:izmax,6+14,dumprange) - Falldiag(i).dg(izmin:izmax,6,dumprange))...
            + squeeze(Fluxdiag(i).dg(izmin:izmax,4,dumprange) +Fluxdiag(i).dg(izmin:izmax,4+14,dumprange) - Falldiag(i).dg(izmin:izmax,4,dumprange))...
            + squeeze(Fluxdiag(i).dg(izmin:izmax,5,dumprange) +Fluxdiag(i).dg(izmin:izmax,5+14,dumprange) - Falldiag(i).dg(izmin:izmax,5,dumprange));
        
    case 27
        pdat(i).p=fact*squeeze(pimlt(i).dg(izmin:izmax,1,dumprange) + psmlt(i).dg(izmin:izmax,1,dumprange) - pgmlt(i).dg(izmin:izmax,1,dumprange));
        
    case 28
        dgfind=findhead('ACC_A',dgstrDan(1).dg);
        A=squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange));
        az=find(A<0.001);
        A(az)=1;
        
        dgfind=findhead('ALL_Q06',dgstrDan(1).dg);
        pdat(i).p=fact*squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange) )./A; 
        
        dgfind=findhead('ALL_Q04',dgstrDan(1).dg);
        pdat(i).p=pdat(i).p + fact*squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange) )./A; 
        
        dgfind=findhead('ALL_Q05',dgstrDan(1).dg);
        pdat(i).p=pdat(i).p + fact*squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange) )./A; 
        
%%%%%%        
    end %switch plotcase
%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of switch plotcase %%%%%%%%%%%%%%%%%%%%%%%%%*******************    
    
    maxC(i)=max(max(pdat(i).p));
    minC(i)=min(min(pdat(i).p));
    
    if maxC(i)==minC(i)  %if have a constant field then use pcolor otherwise contourf crashes
        icont=0;

        if minC(i)==0
            minC(i)=-1;
            maxC(i)=1;
        else
            minC(i)=minC(i)-minC(i)/2;
            maxC(i)=maxC(i)+maxC(i)/2;
        end

    end
    
    if logflag==1
        
        mincovOvr=log10(mincovOvr);
        maxcovOvr=log10(maxcovOvr);

        if iminovr(i)==1
            minC(i)=(mincovOvr);
        else
            minC(i)=log10(min(min(pdat(i).p)));
        end
        
        if imaxovr(i)==1
            maxC(i)=(maxcovOvr);
        else
            maxC(i)=log10(max(max(pdat(i).p)));
        end
        
        minC(minC==-Inf)=NaN;
        maxC(maxC==-Inf)=NaN;
        
        pdat(i).p=log10(pdat(i).p);
        pdat(i).p(pdat(i).p==-Inf)=NaN;
    end
    
    if dlogflag==1
        maxC(i)=dlog(max(max(pdat(i).p)),dlogmin);
        minC(i)=dlog(min(min(pdat(i).p)),dlogmin);
        pdat(i).p=dlog(pdat(i).p,dlogmin);
    end
    
    
    
    %to put extrememes of contours on right side but put mostly a continuation of the right hand value to avoid labels at right edge
%    half(i)=abs((izmax-izmin)/2); %changed Nov 2008 -Dan
    half(i) = round(size(pdat(i).p,1)/2);
    pend(i)=size(pdat(i).p,2);
    npend(i)=max([round(pend*1.25) pend+20]);
    nend(i)=max([round((npend-pend)*0.25) 10]);
    
    dumpint=(timesTH(i).t(end)-timesTH(i).t(1))/10;
    if dumpint==0
        fprintf(1,'***********  PROBLEM - dumpint==0. Check the x-axis values for non-monotonic increase ***********');
        break
    end
    x_grid_orig = timesTH(i).t;
    
    switch plot_type
        case '3D surf'
            %don't change pdat or timesTH(1).t
        otherwise

            switch right_side_extra_bits
                case 1
                    if size(timesTH(1).t,2)==1
                        timesTH(i).t=[timesTH(i).t; timesTH(i).t(end)+(dumpint.*(1:npend(i)-pend(i)))'];
                    else
                        timesTH(i).t=[timesTH(i).t timesTH(i).t(end)+(dumpint.*(1:npend(i)-pend(i)))];
                    end
                case 2
                if size(timesTH(1).t,2)==1
                    timesTH_dummy=[timesTH(i).t; timesTH(i).t(end)+(dumpint.*(1:npend(i)-pend(i)))'];
                else
                    timesTH_dummy=[timesTH(i).t timesTH(i).t(end)+(dumpint.*(1:npend(i)-pend(i)))];
                end
            end

    end
    
    
end





maxCov=max(maxC);
minCov=min(minC);

if isnan(maxCov); maxCov=0; end 
if isnan(minCov); minCov=0; end





%need to make sure highest contour is higher than highest data value for colorbarf

if i2d==1
    timesTH(1).t=Grid.Y1/1000;
    xlabelstr='Horizontal Distance (km)';
elseif i2d==2
    xlabelstr='Horizontal Distance (km)';
elseif i2d==0
    if iutc==1
        xlabelstr='UTC Time (hrs)';
        xlabelstr='Time (hrs)';
        
    else
        xlabelstr='Local Time (hrs)';
    end        
end


if isamescale==1
    iminovr=1;
    imaxovr=1;
    if dlogflag==1
        mincovOvr = dlog(minVal,dlogmin);
        maxcovOvr = dlog(maxVal,dlogmin);
    elseif logflag==1
        mincovOvr = log10(minVal);
        maxcovOvr = log10(maxVal);
    else
        mincovOvr = minVal;
        maxcovOvr = maxVal;
    end
end

for i=1:nplots2d
       
    
    if notsame==1 
        maxCov=maxC(i);
        minCov=minC(i);
    end
    
    nsigfig=1+ceil(abs(log10(maxC-minC))); %attempt to guess a sensible number of sigfigs
    
    if abs(minC(i))>=100 | abs(maxC(i))>=100 %if length(num2str(round(maxC(i)))) - length(num2str(round(minC(i)))) == 0
        nsigfig = max(nsigfig,3);
    elseif abs(minC(i))>=10 | abs(maxC(i))>=10
        nsigfig = max(nsigfig,2);
    else
        nsigfig = max(nsigfig,1);
    end
    
    if iminovr(i)==1    %imincovovr==
        minCov=mincovOvr(i);
    else        
        mincovOvr(i) = sigfig(minC(i),nsigfig);   %Dan - added to make nicer numbers on contours and colorbar
        minCov=mincovOvr(i);
        iminovr(i)=1;
    end
    
    if imaxovr(i)==1    %imaxcovovr==
        maxCov=maxcovOvr(i);
    else	                
        maxcovOvr(i) = sigfig(maxC(i),nsigfig);
        maxCov=maxcovOvr(i);
        imaxovr(i)=1;
    end
    
    minc=minCov-0.0*abs(minCov);
    maxc=maxCov+abs(maxCov*0.0);
    
    if minc>maxc; m=minc; minc=maxc; maxc=m; end
    
    %     fixmin=fix2(minc,abs(round(log10(min(abs(minc))))));
    %     dd=(maxc-minc)/ncont;
    %     dfix=fix2(dd,ceil(abs(log10(min(abs(dd))))));
    %     conts=[fixmin:dfix:maxc];
    %conts=[minc:(maxc-minc)/ncont:maxc];
    %conts=round2(conts,abs(round(log10(min(abs(conts)))))+1);
    
    if logflag==0 & dlogflag==0
        dd=(maxc-minc)/ncont;
        %dfix=fix2(dd,ceil(abs(log10(min(abs(dd))))));
        
        dfix=sigfig(dd,0);
        if dfix==0;
            dfix=sigfig(dd,1);
        end
        if minc>=100
            fixmin=sigfig(minc,sig);
        else
            fixmin=sigfig(minc,sig);     
        end
        if abs(fixmin)<dfix/100; fixmin=0; end
        conts=[fixmin:dfix:maxc];
        iszero=zeros([1 length(conts)+1]); %flag to say whether value is exactly zero
        
        if length(conts)==0;
            conts(1)=0;
            conts(2)=0;
        end
        
        if sign(conts(1))~=sign(conts(end))
            [zeromin izero]=min(abs(conts));
            if abs(conts(izero))<dfix/2.2
                conts(izero)=0;
                iszero(izero)=1; %flag to say that value is exactly zero
            else
                if conts(izero)<0
                    conts(izero+2:end+1)=conts(izero+1:end);
                    conts(izero+1)=0;
                    iszero(izero+1)=1;
                else
                    conts(izero+1:end+1)=conts(izero:end);
                    conts(izero)=0;
                    iszero(izero)=1;
                end
            end
        end
    elseif dlogflag==1
        conts=[minc:(maxc-minc)/ncont:maxc];
        iszero=zeros([1 length(conts)+1]);
        
        unlog=idlog(conts,dlogmin);
        if sign(unlog(1))~=sign(unlog(end))
            [zeromin izero]=min(abs(unlog));
            if abs(unlog(izero))<((maxc-minc)/ncont)/2.2
                unlog(izero)=0;
                iszero(izero)=1;
            else
                if unlog(izero)<0
                    unlog(izero+2:end+1)=unlog(izero+1:end);
                    unlog(izero+1)=0;
                    iszero(izero+1)=1;
                else
                    unlog(izero+1:end+1)=unlog(izero:end);
                    unlog(izero)=0;
                    iszero(izero)=1;
                end
            end
        end
        
        unlog=sigfig(unlog,sig);
        conts=dlog(unlog,dlogmin);
        
    else
        conts=[minc:(maxc-minc)/ncont:maxc]; %logflag==1
        iszero=zeros([1 length(conts)+1]);
    end
    
    
    if length(conts)==0
        conts=[0 1];
    end
    
    %     ac=find(conts>0);
    %     ac2=find(conts<0);
    %     if length(ac)>0 & ac(1)>1
    %         ac=ac(1);
    %         conts=[conts(1:ac-1) 0 conts(ac:end)];
    %     end
    
    if iovride_conts==1
        conts=conts_ovr;
    end
    
    
    
% put in the extra bits at the side to get the filled contour/filled colorbar right    
    psame=repmat(pdat(i).p(:,pend(i)),[1 npend(i)-pend(i)-nend(i)-1]);
%    psame(isnan(psame))=0;

switch plot_type
    case '3D surf'
        %do nothing
    otherwise
        pdat(i).p(:,pend(i)+1:npend(i)-nend(i)-1)=psame;
end
    
    
    mindat=minCov*(1-0.01*sign(minCov));   
    maxdat=maxCov*(1+0.01*sign(maxCov));

    switch plot_type
        case '3D surf'
            %do nothing
        otherwise
            switch right_side_extra_bits
                case 1
                    pdat(i).p(1:half(i),npend(i)-nend(i):npend(i))=mindat;
                    pdat(i).p(half(i)+1:end,npend(i)-nend(i):npend(i))=maxdat;
                case 2
                    pdat_dummy = pdat(i).p;
                    pdat_dummy(1:half(i),npend(i)-nend(i):npend(i))=mindat;
                    pdat_dummy(half(i)+1:end,npend(i)-nend(i):npend(i))=maxdat;
            end

    end



    
    if subplotting==1 & noplot==0
        %h(isub).h=subplot(a,b,isub);
    end
    
    
    if izovr==0;
        zz(i).z=z(izmin:izmax);
    end
    
    %pcolor(timesTH,0.62+z(izmin:izmax)./1000,pdat);
    if ilem==1;
        height=add_ground_height+zz(i).z./1000;
    else
        height=zz(i).z;
    end        
    
    
    
    if noplot==1
        return %exit if don't want to plot
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%    draw the plot    %%%%%%%%%%%%%%%%%%%%%%%%%
%    set(gcf,'renderer','painters');
    
%     if iplot_latlon==1
%         plot_latlon_lines2;
%     end
    
    switch plot_type
        case 'filled contour'
            if clines==1
                %                [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts); %old version
                %                 [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
                %                 [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
                %                 hold on

                [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
                hold on
                %[cbfA(i).cpos cbfB(i).cpos]=contourf(timesTH(i).t,height,pdat(i).p,conts(conts<0),'k','linestyle',':');
                %[cbfA(i).cneg cbfB(i).cneg]=contourf(timesTH(i).t,height,pdat(i).p,conts(conts>=0),'k','linestyle','-');
                %changed the above two contourf calls to contour calls 26/03/09
                [cbfA(i).cpos cbfB(i).cpos]=contour(timesTH(i).t,height,pdat(i).p,conts(conts<0),'k','linestyle',':');
                [cbfA(i).cneg cbfB(i).cneg]=contour(timesTH(i).t,height,pdat(i).p,conts(conts>=0),'k','linestyle','-');
                icont_neg=1;
                %hold on
                % [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');



            else
                [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
                cbfA(i).cpos = cbfA(i).c;
                cbfB(i).cpos = cbfB(i).c;
            end
            
        case '3D surf'
%            surf(x_grid,y_grid,squeeze(terrain),'facecolor','red','edgecolor','none');
%            camlight headlight; lighting phong
%            camlight headlight; lighting phong
            
            surf(x_grid,y_grid,squeeze(Zlev),pdat(1).p(:,1:pend),'edgecolor','none');
            
%            alpha(0.7);
%            shading interp
            set(gca,'CameraPosition',CameraPosition);
            set(gca,'CameraTarget',CameraTarget);
            set(gca,'CameraViewAngle',CameraViewAngle);
            
            camlight headlight; lighting phong
            camlight headlight; lighting phong
            alpha(0.95);
%            camlight headlight; lighting phong

            iadd_terrain = 1;
            iadd_wind_quivers = 0;
            iplot_latlon=0;
            iadd_flight_path=0;
            isquare=0;
        otherwise
            if idpcolor==1
                [cbfA(i).c hback]=dpcolor(timesTH(i).t,height,pdat(i).p);
            else               
                 hback2=pcolor(timesTH(i).t,height,ones(size(pdat(i).p))); %if add two backgrounds on then it reduces the aliasing hatching
                 %in the final PDF
                 hold on
                 hback=pcolor(timesTH(i).t,height,ones(size(pdat(i).p)));
                 hold on
                 [cbfA(i).c]=pcolor(timesTH(i).t,height,pdat(i).p);
            end
            
            switch pcolor_shading
                case 'interp'
                    shading interp;
                otherwise
                    shading flat;
            end
    end
    
    
end
    

   
    if icont_extra==1
        hold on
        [c2,h2]=contour(timesTH(i).t(1:length(xinds)),add_ground_height+zz(i).z/1000 ...
        ,f/1000*sum(TwoDDan(idir).Q(izmin:izmax,xinds,4:6),3),[0:5:20],'w');
        set(h2,'linewidth',1.5);
        clabel(c2,h2,'color','w','fontsize',max([fsize-8 6]));
%        clabel(c2,h2,'labelspacing',72,'color','w');
    end
    

    if iadd_cont==1
        if idiff_cont_dat~=1
            cont_dat_alt=pdat(1).p;
            xcont2 = timesTH(1).t;
            ycont2 = zz(1).z;
        end
        
        if idpcolor==1
            xcont = 0.5*(xcont2(1:end-1)+xcont2(2:end));
            ycont = 0.5*(ycont2(1:end-1)+ycont2(2:end));
        else
           xcont =  timesTH(1).t;
           ycont = zz(1).z;            
        end
        
        
        [ccont,hcont] = contour(xcont,ycont,cont_dat_alt,cont_val,'k--','linewidth',4);
    end
    
    
    
    if isquare==1
        axis square
    end
    

    switch plot_type
        case 'filled contour'
            %nothing
        otherwise
            cax=get(gca,'clim');
            if iminovr(i)==1
                if logflag==1
                    cax(1)=mincovOvr;
                    if isnan(cax(1))
                        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                        disp('****   Set a mincovOvr for the minimum value when logflag==1 ****');
                        disp('***************************************');
                    end
                else
                    cax(1)=mincovOvr;
                end
            end
            
            if imaxovr(i)==1
                if logflag==1
                    cax(2)=maxcovOvr;
                else
                    cax(2)=maxcovOvr;
                end
            end
            caxis(cax);
            %        hc=colorbar;
            %        set(hc,'fontsize',fsize-3);


            normcbar=1;


    end
    
    %hc=colorbarf(cbfA,cbfB);
    
    
    
    %xti=set(h(i).h,'xticklabels',timestxt);
    
    set(h(isub).h,'fontsize',fsize);
    
    if onexlabel==0 | (isub==nplots2d & onexlabel==1 & subplotting==0) | (isub==length(idirs) & onexlabel==1 & subplotting==1)
        xlabel(xlabelstr);
    end
    
    if izovr~=2
        ylabel('Height (km)');
    else
        ylabel(ylabelstr);
    end
    
    tit(i).tit = remove_character(tit(i).tit,'_','-');
    
    if subplotting==1 & isub==1
        title(tit(i).tit,'fontsize',fsize_tit);
    elseif subplotting==0
        titw=textwrap({tit(i).tit},100); %wrap the title onto two lines if >100 chars
        title(titw,'fontsize',fsize_tit);
    end
    
    
    if exist('cbfA')
        if length(cbfA(i).c)==0
            nocbar=1;
            normcbar=1;
        end
    end
    

    %caxis(h(i).h,[minCov maxCov*1.05]);
    
    
    
    
    
    if exist('idirs');
        if isub==length(idirs)
            %    if isamescale==1 &  isub==length(idirs)

            pos1=get(h(1).h,'position'); %[left bottom width height]
            posend=get(h(end).h,'position');
            
            height_pos=pos1(2)-posend(2)+pos1(4);
            pos=[posend(1) posend(2) posend(3)+0.15 height_pos];
        else
            pos=[0 0 1 1];
        end
    end
    

    if isamescale==1 &  isub==length(idirs)  
        if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            hc=colorbar( 'peer' , h(isub).h , 'location',colorbar_location);
        elseif normcbar==1 & bigcbar==0 | (bigcbar==1 & isub==length(idirs)) | lememm==1
            axdan=axes('position',pos,'visible','off');
            %colbar=colorbar; %if colorbar already in place then colorbarf will replace it
            %set(colbar,'tag','Colorbar');
            hc=colorbarf(cbfA(i).c,cbfB(i).c,colorbar_loc);
            
            
        end
    else
        if nocbar~=1 & strcmp(plot_type,'filled contour')
            if normcbar==0 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
          %%%%%   usual one %%%%%%%%      
                hc=colorbarf_30Sep09_2(cbfA(i).c,cbfB(i).c,colorbar_loc); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!  
%                hc=colorbarf(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!                  
%    hc=cbarf(conts,conts,'vertical','linear');                
            elseif normcbar==0 &(bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
                axdan=axes('position',pos,'visible','off')
                %   colbar=colorbar;
                %   set(colbar,'tag','Colorbar');
                hc=colorbarf(cbfA(i).c,cbfB(i).c,colorbar_loc);
%             
                
            end 
        end
        
        if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            hc=colorbar( 'peer' , h(isub).h , 'location',colorbar_location);
        end
        
    end
    
    %   if normcbar==0
    
    if (logflag==1 | dlogflag==1 ) & nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
        %re-label colorbar ticks if have log settings  
        clear ctickstr
        
        ctick=get(hc,'yticklabel');
        if normcbar==1 %since for normal cbar ctick returned as matrix of characters
            ctick2=ctick;
            clear ctick
            for j=1:size(ctick2,1)
                ctick{j}=ctick2(j,:);
            end
            jinds=1:length(ctick);
        else
            jinds=2:length(ctick)-1;
        end
        
        ctickstr(1,1)=' ';
        for j=jinds
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            nu=str2num(ctick{j});
            
            
            
            %te=num2str(10^nu,'%2.2f');       %'%2.2e');
            if dlogflag==0
                if 10^nu>99
                    te=num2str(sigfig(10^nu - offset ,3));
                else
                    te=num2str(sigfig(10^nu - offset,2));
                end
            else
                %                     if idlog(nu,dlogmin)>99
                %                         te=num2str(sigfig(idlog(nu,dlogmin) - offset ,4));
                %                     else
                %                         te=num2str(sigfig(idlog(nu,dlogmin) - offset,4));
                %                     end
                te=num2str(sigfig(idlog(nu,dlogmin) - offset,2),4);
                if normcbar==1 
                    if iszero(j-1)==1  %if is suppossed to be exactly zero then set to zero
                        te='0.0';
                    end
                else
                    if iszero(j-1)==1  %if is suppossed to be exactly zero then set to zero
                        te='0.0';
                    end
                end
                
            end
            
            ctickstr(j,1:length(te))=te;
        end 
        
        set(hc,'yticklabel',ctickstr);
        set(hc,'Fontsize',fsize);
        
        %         add=str2num(ctick{end-1})/50;
        %         
        %         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        %         
        %         for i=2:length(ctick)-1
        %             cticknums(i)=str2num(ctick{i});
        %         end
        %         text(  ones( length(cticknums),1 )*1.05,cticknums+add,ctickstr, 'fontsize',fsize-6  );
        
    elseif (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) ) %if (logflag==1 | dlogflag==1 ) & nocbar=...
        %re-label colorbar if not log plot
        clear ctickstr
        switch colorbar_loc
            case 'vert'
                ctick=get(hc,'yticklabel');
            case 'horiz'
                ctick=get(hc,'xticklabel');
        end
        
        if normcbar==1 %since for normal cbar ctick returned as matrix of characters
            ctick2=ctick;
            clear ctick
            for j=2:size(ctick2,1)+1
                ctick{j}=ctick2(j-1,:);
            end
            jinds=2:length(ctick);
        else
            jinds=2:length(ctick)-1;
        end
        
        
        ctick2=ctick;
        clear ctick
        for j=jinds
            ctick(j-1)=str2num(ctick2{j});
        end
        ctickstr(1,:)=' ';
        
        for j=2:length(ctick)+1
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            %te=num2str(ctick(j),'%2.2e');
            te=num2str(sigfig(ctick(j-1),sig));
            if  normcbar==0
                ctickstr(j,1:length(te))=te;
            else
                ctickstr(j-1,1:length(te))=te;
            end                
        end
        
        if no_set_ctick_str==0 %the changing of the tick strings was causing problems when the font
            %size was being changed - caused the colorbar values to change!

            fprintf(1,'\n******* WARNING - colorbar tick strings are being changed. Causes problems with values changing when fontsize is changed!!! ******* \n\n\n');
            switch colorbar_loc
                case 'vert'
                    set(hc,'yticklabel',ctickstr);
                case 'horiz'
                    set(hc,'xticklabel',ctickstr);
            end

        end
        
        
        
        
        %         add=ctick(end)/50;
        %         
        %         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        %         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
        
        
        
        
        
        %set(hc,'fontsize',fsize-6);
        
        
        
        %else
        
        
        
        %     if logflag==1 & nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
        %         clear ctickstr
        %         ctick=get(hc,'ytick');
        %         for j=1:length(ctick)
        %             %te=strcat('10^','{',num2str(ctick(j)),'}' );
        %             %te=num2str(10^ctick(j),'%2.2g');
        %             if nu>99
        %                 te=num2str(sigfig(10^ctick(j) - offset,2));
        %             else
        %                 te=num2str(sigfig(10^ctick(j) - offset,2));
        %             end
        %             
        %             ctickstr(j,1:length(te))=te;
        %         end
        %         
        %         set(hc,'yticklabel','');
        %         
        %         add=ctick(end)/50;
        %         
        %         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        %         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize);
        %         
        %     
        %     
        %     elseif nocbar==0 & icont==1 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
        %         clear ctickstr
        %         ctick=get(hc,'ytick');
        %         for j=1:length(ctick)
        %             nu=ctick(j);
        %             %te=strcat('10^','{',num2str(ctick(j)),'}' );
        %             if nu>99
        %                 %te=num2str(ctick(j),'%2.2e');
        %                 te=num2str(sigfig(ctick(j),2));
        % 
        %             else
        %                % te=num2str(ctick(j),'%2.2e');
        %                 te=num2str(sigfig(ctick(j),2));
        %             end
        %             ctickstr(j,1:length(te))=te;
        %         end
        %      end
        
        
        %end
        
        if nocbar==0 & subplotting==1 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            set(hc,'fontsize',fsize);
        elseif nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            set(hc,'fontsize',fsize);
        end
        
        
        
        
    end
    


        ylims_temp=get(h(isub).h,'ylim');   %re-scale to hide extra column put in to get the colorbars the same in both plots
        
        switch plot_type
            case '3D surf'
                %do nothing
            otherwise
                axis(h(isub).h,[timesTH(i).t(1) timesTH(i).t(pend(i)) ylims_temp]);
        end
        

    if dan_test==1 %%% dan test %%%
        cax=get(gca,'clim');
        if iminovr(i)==1
            cax(1)=mincovOvr;
        end
        if imaxovr(i)==1
            cax(2)=maxcovOvr;
        end
        caxis(cax);
%        hc=colorbar;
%        set(hc,'fontsize',fsize-3);
        
    end
 





    

if clab==1 & strcmp(plot_type,'filled contour') %round up contour values so labels match those on colorbar

    for iclabel_count=1:icont_neg+1
        
        if iclabel_count==1
            ch=cbfA(i).cpos;
            cb_dat = cbfA(i).cpos;
            c_handle = cbfB(i).cpos;
        else
            ch=cbfA(i).cneg;
            cb_dat = cbfA(i).cneg;
            c_handle = cbfB(i).cneg;
        end
        jc=1;
        while jc<size(cb_dat,2)
            %ch(1,jc)=str2num(num2str(10^(cbfA(i).c(1,jc)),'%2.2e'));

            if dlogflag==0 & logflag==1
                if 10^cb_dat(1,jc)>99
                    ch(1,jc)=sigfig(10^(cb_dat(1,jc) ),3);
                else
                    ch(1,jc)=sigfig(10^(cb_dat(1,jc) ),2);
                end
            elseif dlogflag==1
                %                    if idlog(cbfA(i).c(1,jc),dlogmin)>99
                %                         ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin ),2);
                %                     else
                %                         ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin ),2);
                %                     end
                ch(1,jc)=sigfig(idlog(cb_dat(1,jc),dlogmin),3);
            else
                ch(1,jc)=sigfig(cb_dat(1,jc),sig);
            end

            jc=jc+cb_dat(2,jc)+1;
        end

        if length(ch)>0
            if manclab==0
                clabel(ch,c_handle,'labelspacing',144); %default 144
            else
                clabel(ch,c_handle,'manual');
            end
        end

    end  %for iclabel_count=1:icont_neg+1
    icont_neg=0;

end %clab==1 & icont==1
    
    
    %if clines==0
     %   shading flat; %shading stops black contour lines from appearing
    %end  %NOTE - needs to be after the contour labelling as this puts the
    %lines back! 
    %This method produced incorrect results when have small NaN areas
    %filled white - using contourf(.....,'linestyle','none') is better



if vectorf==1
    spy=25;
    spz=15;
    
    %    spy=20;
    %    spz=10;
    
    sqy=size(GridDan(idir).Y1(xinds),1)/spy;
    sqz=round((izmax-izmin)/spz);
    
    zinds=[izmin:sqz:izmax+2*sqz];
    yinds=[xinds(1):sqy:xinds(end)+sqy];
    
    sf=max(max(TwoD.V))/max(max(TwoD.W));
    hold on;
    quiver(GridDan(idir).Y1(yinds)./1000,GridDan(idir).Z(zinds)./1000,TwoD.V(zinds,yinds),TwoD.W(zinds,yinds),'w');
end


if itimestamp==1
    %    text(timesTH(i).t(1)-(timesTH(i).t(end)-timesTH(i).t(1))*0.12,((zz(1).z(end)-zz(1).z(1))*0.18+zz(1).z(end))/1000,['Time = ' timlab ' UTC'],'fontsize',18);
    %    text(timesTH(i).t(1)-(timesTH(i).t(end)-timesTH(i).t(1))*0.12,((zz(1).z(end)-zz(1).z(1))*0.08+zz(1).z(end))/1000,['Time = ' timlab ' UTC'],'fontsize',18);
    
    if subplotting==1
 %       text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-1.5 16],'fontsize',fsize);
  %     text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-1.5 13.2],'fontsize',fsize);
 %       text(0,0,['Time = ' timlab],'units','centimeters','position',[-1.5 13.2],'fontsize',fsize);
      text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-1.5 18],'fontsize',fsize);
        
    else     
       %  text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-2.5 11.5],'fontsize',fsize);
        text(0,0,['Time = ' timlab],'units','centimeters','position',[-2.5 13],'fontsize',fsize);
        
    end   
    
    
    % text(timesTH(i).t(1)*1.2,23.0,['Time = ' f ' UTC'],'fontsize',18);
end




% 
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');

if lememm==1   %rescale if doing lem/emm comparison so as to get same axes
    set(gca,'xlim',timelims);
    set(gca,'ylim',zlims);
    if length(clims)==2
        %        set(gca,'clim',clims);             
    end
end  


if idirstamp==1
    %     if subplotting==1
    %         ylims=get(h(iplot).h,'ylim'); 
    %         xlims=get(h(iplot).h,'xlim'); 
    %     else
    %         ylims=get(h.h,'ylim'); 
    %         xlims=get(h.h,'xlim'); 
    %     end        
    
    ylims=get(h(iplot).h,'ylim'); 
    xlims=get(h(iplot).h,'xlim'); 
    
    %    text(timesTH(i).t(1)-0.5,ylims(2)*1.002,[direcDan(idir).dir],'fontsize',12);
    if subplotting==1
        dist=(ylims(2)-ylims(1))/7.8;
	else
        dist=(ylims(2)-ylims(1))/10.1;
    end
    
    %dist2=(timesTH(i).t(pend(i))-timesTH(i).t(1))/20;
    dist2=(xlims(end)-xlims(1))/8;
%      dist2=(xlims(end)-xlims(1))/20; %last one used
      
    %    text(timesTH(i).t(1)-0.5,ylims(1)*0.988,[direcDan(idir).dir],'fontsize',12);    
    axes(h(end).h);
    
    if plotcase==65
        dirname=run_name_emm_select;
    else
        dirname=runName(idir).nam
    end
    
    if iabc==1
        abc={'(a) ','(b) ','(c) ','(d) '};
        dirstr=[abc{iplot} dirname];
    else
        dirstr=dirname;
    end
    
%    text(xlims(1)-dist2,ylims(1)-dist,[dirstr],'fontsize',fsize-2);    
    text(xlims(1)-dist2,ylims(1)-dist,[dirstr],'fontsize',fsize-2);    
    %	text(xlims(1)-4*dist2,ylims(2)+dist,[dirstr],'fontsize',fsize+4);
end


if (i2d~=1 & i2d~=2 & i2d~=3)
    xx=get(h(isub).h,'xticklabels');
    xx=str2num(xx);
    xx=num2str(mod(xx,24));
    set(h(isub).h,'xticklabels',xx);
end

%set(gcf,'paperpositionmode','auto');
set (gcf, 'Papersize',[3 3],'Paperposition',[0 0 3 3 ]); 

%set (gca, 'Fontsize', 10); 

if isave==1
    set(gcf,'paperpositionmode','auto');
    print(gcf,'-djpeg','-r350',exname);
    %print(gcf,'-dmeta',exname);
    %close(gcf);
end



clims_terr = get(gca,'clim'); %need to remember and reset the color lims as goes weird when do contours

try

    if iadd_terrain==1
        add_terrain;
    end

    if icolmap==1
        colormap(cmap);
    end

     if iplot_latlon==1
         plot_latlon_lines2;
     end

    if iadd_overlay>=1
        hold on
        for iover=1:iadd_overlay
            eval_str=['plot(dat_overlay(iover).x,dat_overlay(iover).y,' dat_overlay(iover).overlay_style ');'];
            eval(eval_str);
        end
    end

    if iadd_flight_path==1 & exist('time_flt19')
        plot_flight_path;
        if iadd_AWS==1
            LAT=LAT_orig;
        end
    end

catch error_flight_path

    disp('*********** ERROR - error occurred in adding terrain or flight path **********');

end
    
caxis(clims_terr); %resest colourscale as sometimes contour causes it to go wrong




if iylim==1
    set(gca,'ylim',ylims);
end
if ixlim==1
    set(gca,'xlim',xlims);
end

switch y_axis_type        
    case {'log10_matlab'} %log10_matlab is when the 'yscale' is set to 'log' using matlab built in feature
        set(gca,'yscale','log');
        iytick_override=0;
end
    
    
if iytick_override==1
    set(gca,'YTick',tick_locs);
    set(gca,'Yminortick','off')
end

if iytick_relabel==1
    ticklabels = get(gca,'yticklabel');
    tick_nums = str2num(ticklabels);
    
    switch y_axis_type        
        case {'log10','log10_matlab'} %log10_matlab is when the 'yscale' is set to 'log' using matlab built in feature
            new_tick_nums = 10.^tick_nums;            
    end
    
    clear new_ytick_text
    
    for itick=1:length(new_tick_nums)
        %round to one sig fig
        if i_only_major_log_ticks==1
            if round(log10(new_tick_nums(itick)))-log10(new_tick_nums(itick)) == 0
                ytick_text_i = num2str(sigfig(new_tick_nums(itick),1));
                new_ytick_text(itick,1:length(ytick_text_i)) = ytick_text_i;
            else
                new_ytick_text(itick,1)=' ';
            end
        else
            ytick_text_i = num2str(sigfig(new_tick_nums(itick),1));
            new_ytick_text(itick,1:length(ytick_text_i)) = ytick_text_i;
        end
    end
    
    set(gca,'yticklabel',new_ytick_text);
    
end
        
        
if iadd_wind_quivers==1

    [ax_orig]=plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');
    %for some reason it only plots the reference arrow head if run twice!
%    plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1);
    axes(ax_orig);  %revert back to the original axis for zooming in etc
    set(gca,'ylim',[zz(1).z(1) zz(1).z(end)]); %need to rescale the vertical axis
    
    set(gca,'outerposition',[0 0 0.92 1]); %this moves the outer box around the plot sideways so that the wind quiver scale
    %arrow doesn't get cut off.
end

if idatetick==1; %flag to say the want the xaxis in proper time format rather than decimal time
    %specify the type with datetick_type (see help datetick)
    datetickzoom('x',datetick_type,'keeplimits');
    %note need to plot time in days after day 0 of year 0
    
    if ixtickmode_auto==1
        set(gca,'xtickmode','auto');
    end
end

       
        
        if idraw_streamlines==1
            draw_streamlines;
        end

%pause(60);

set(gca,'TickDir',tick_direction);

if ione_to_one_line==1
    ylim_scatter=get(gca,'ylim');
    xlim_scatter=get(gca,'xlim');
    if(strcmp(get(gca,'ydir'),'reverse'))
        ytemp=ylim_scatter(1);
        ylim_scatter(1)=ylim_scatter(2);
        ylim_scatter(2)=ytemp;
    end
    if(strcmp(get(gca,'xdir'),'reverse'))
        xtemp=xlim_scatter(1);
        xlim_scatter(1)=xlim_scatter(2);
        xlim_scatter(2)=xtemp;
    end
    
    L1=min(xlim_scatter(1),ylim_scatter(1));
    L2=min(xlim_scatter(2),ylim_scatter(2));    
    
    hone=line([L1 L2],[L1 L2]);  %    line([0 10],[0 10]);
    set(hone,'color','k');
    set(hone,'linewidth',2);
    
end


mincovOvr=mincovOvr_save; %reset these to avoid confusion for outside scripts
maxcovOvr=maxcovOvr_save;

titlenam = tit(1).tit;
picname=[savename];

%set(gca,'Color',[0 0 0]);  %perversely, this actually seems to undo the
%good work of the hback color below... maybe it makes [0 0 0] transparent
%or something??
%set(gcf,'inverthardcopy','off'); %this stops the black form being turned
%white, but it gives you a gray figure background...

if exist('hback')
    set(hback,'facecolor',[0 0 0]);
end

title_nice = short_plot_name;

if ipost_plotTime_commands==1
   for ipost=1:length( post_plotTime_commands )
      eval( post_plotTime_commands{ipost} );
   end
end

%clear flags in case of no errors
clear man_choose_plotTimeHeight_graph man_choose_plotTimeHeight_graph2

%and in case of errors too
catch plotTimeHeight_ERROR
    %clear the flags that might cause trouble in case of an error
    ikeep_X_above_zero=0; %reset this flag
    iexec_post_commands=0;
    clear man_choose_plotTimeHeight_graph man_choose_plotTimeHeight_graph2
    rethrow(plotTimeHeight_ERROR); %"re-issue" the error (also creates the links to the error etc)
    
end
    

   
