
plotcase=50;
plotcase=65; %emm
plotcase=577;


iabc=1; %flag to tell it to add (a), (b), (c), etc. to runName(i).nam
comp='uni';
iovride_conts=0;

isquare=0; %flag to make axes square

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

clear zz timesTH diff

if ~exist('idir'); idir=1; end
if ~exist('ieps'); ieps=0; end
if ~exist('isamescale'); isamescale=0; end
if ~exist('subplotting'); subplotting=0; end
if ~exist('onexlabel'); onexlabel=0; end      %flag to make it so that only the bottom plot has the xlabel on it
if ~exist('iplot'); iplot=1; end 

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
    fsize=14; %26 best for paper figures when more than one plot is being used
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

fprintf('\n**** add_ground_height = %f ********\n',add_ground_height);

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

%iminovr=zeros([1 10]);
%imaxovr=zeros([1 10]);
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
logflag=0;
dlogflag=0;




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
    logflag=0;
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
    i577 = 'wrf_wind2d';
%    i577 = 'wrf_plot';
%    i577 = 'wrf_radar_vert';
    
%    i3d=0;
    
    
    imaxovr=0;
    iminovr=0;
    
    clines=1; %makes black contour lines appear
    clab=0;
    
    ilem=0;
    
    
    i2d=2; %tells it to label x axis in km
    
    
    switch i577
    case 'wrf_radar_vert'
                
        tit(1).tit=['Latent Heat Flux (W m^{-2})'];
%        tit(1).tit=['Sensible Heat Flux (W m^{-2})'];
        tit(1).tit=['Radar reflectivity'];              
        tit(1).tit=['Total water'];              
                        
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
        
        ylabelstr='Height (km)';
        
        sig=3;
        
        clab=0;
        
        iovride_conts=1;
        conts_ovr=[0:5:65];
		conts_ovr=[6:0.1:7];
        
    case 'wrf_plot'
        
        
        tit(1).tit=['Latent Heat Flux (W m^{-2})'];
%        tit(1).tit=['Sensible Heat Flux (W m^{-2})'];
        tit(1).tit=['Sea Ice Flag'];
        
        ih_wrf=90; %z-level to plot for

                        
        tit(1).tit=['Total water at level ' num2str(ih_wrf) ' (ppmv)'];
        
        iylim=0;
        iylims=[1 100];
        ixlim=0;
        xlims=[101 201];
    
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
        
        ylabelstr='iLatitude';
        xlabelstr='iLongitude';
        
        sig=3;
        
        clab=0;
        
    imaxovr=1;
    maxcovOvr=995;
    
    iminovr=1;
    mincovOvr=985;
    


    case 'wrf_wind2d'
        
        time=47;
        ih_wrf=5; %z-level to plot for
	Zlev=WRFUserARW(nc,'Z',time,ih_wrf);

        tit(1).tit=['WRF wind speed at ' num2str(Zlev(1,1,400)) 'm (m s^{-1}) at ' Times(time,12:16) ' UTC'];
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
        
        ylabelstr='iLatitude';
        xlabelstr='iLongitude';
        
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
        
        
    end
    
    
    % tit(2).tit='Ice Number Concentration (mg^{-1})';
    nplots2d=1;
    
    
    
    
    
    
    if izovr==0 
        z=GridDan(idir).Z;
    end
    
    notsame=1;
    
    
    %    iminovr=[0 0];
    %    imaxovr=[0 0];
    
    
    figlab=[tit(1).tit ' 2d plot ' i577];
    savename=[figlab];   % ' ' dirname];
    
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
    
    savename=figlab;
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
    imaxovr=1;t
    maxcovOvr=0.6;
    ncont=17;
    
    
    
    
    
    
    
    
    
    
end
%timesTH=[hrstartles+((dumprange-1)*dumpint)/3600];

scrsz=get(0,'ScreenSize');
%posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

switch comp
case 'uni'
%    posit=[9 60 scrsz(3)/1.2 scrsz(4)/1.14];
    posit=[9 60 scrsz(3)/1.7 scrsz(4)/1.8];
case 'lacieLap'
    posit=[9 50 scrsz(3)/1.46 scrsz(4)/2.07];
end

if subplotting==1
%    posit=[9 50 scrsz(3)/2.4 scrsz(4)/1.5];
    posit=[9 50 scrsz(3)/2 scrsz(4)/1.1];
    
    if nsub==1; hf=figure('position',posit,'name',figlab); end
    h(iplot).h=subplot(xsub,ysub,nsub); 
    
    %     for ih=1:length(h)
    %         posh=get(h(ih).h,'position');
    %         set(h(ih).h,'position',[posh(1) posh(2) posh(3)-0.2 posh(4)]);
    % 	end
    
    posh=get(h(iplot).h,'position');
    set(h(iplot).h,'position',[posh(1) posh(2) posh(3)-0.2 posh(4)]);
    
    
    
else
    hf=figure('position',posit,'name',figlab);
    h(iplot).h=subplot(1,1,1); 
end

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
    
    scrsz=get(0,'ScreenSize');
    posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
    
    if izovr==0; [izmin izmax]=findheight(z,minZ,maxZ); end
    
    
    %     if length(iz)>=1
    % 		iz=iz(1);
    %     else
    %         iz=length(z);
    %     end
    
    
    %pcolor(9+time./3600,z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp
    
    if subplotting==0    
        h(isub).h=subplot(a,b,isub);
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
        
        switch i577
            
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
        
            
            
        case 'wrf_plot'
            zz(1).z = 1:size(lat2d(1).var,1);
            timesTH(1).t = 1:size(lat2d(1).var,2);                       
            
            time=15;
%             savemem=0;
%             switch savemem
%             case 0
%                Z = WRFUserARW(nc,'Z',time);


%                pdat(1).p = nc{'LH'}(time,:,:);
%                pdat(1).p = nc{'GRDFLX'}(time,:,:);
%                pdat(1).p = nc{'HFX'}(time,:,:);
%                pdat(1).p = nc{'QFX'}(time,:,:);
               % pdat(1).p = nc{'SEAICE'}(time,:,:);

		 pdat(1).p = WRFUserARW(nc,'p',time,1);

                itime=idir;
%                pdat(1).p = f* (nc{'QICE'}(itime,ih_wrf,:,:)+nc{'QSNOW'}(itime,ih_wrf,:,:)+nc{'QGRAUP'}(itime,ih_wrf,:,:)+nc{'QVAPOR'}(itime,ih_wrf,:,:) );

                
              %  v = 0.5* (nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) ); %2d wind at one height
                
%                 for ilat=1:length(zz(1).z)
%                     for ilon=1:length(timesTH(1).t)
%                         ih_wrf=findheight(Z(:,ilat,ilon),h_wrf);
%                         if length(ih_wrf)==0
%                             ih_wrf=1;
%                         end
%                         pdat(1).p(ilat,ilon) = sqrt( u(ih_wrf,ilat,ilon)^2 + u(ih_wrf,ilat,ilon)^2 );
%                     end
%                 end

            %    pdat(1).p=squeeze(WRFUserARW(nc,'Z',time,0,0,ih_wrf));  %horiz_slice at first level
    
            
        case 'wrf_wind2d'
             zz(1).z = 1:size(lat2d(1).var,1);
            timesTH(1).t = 1:size(lat2d(1).var,2);
            
            h_wrf=250; %(m)
            
            %time=46;
            savemem=0;
            switch savemem
            case 0
%                Z = WRFUserARW(nc,'Z',time);



                u = 0.5* (nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) ); %2d wind at one height            
                v = 0.5* (nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) ); %2d wind at one height
                
%                 for ilat=1:length(zz(1).z)
%                     for ilon=1:length(timesTH(1).t)
%                         ih_wrf=findheight(Z(:,ilat,ilon),h_wrf);
%                         if length(ih_wrf)==0
%                             ih_wrf=1;
%                         end
%                         pdat(1).p(ilat,ilon) = sqrt( u(ih_wrf,ilat,ilon)^2 + u(ih_wrf,ilat,ilon)^2 );
%                     end
%                 end

                pdat(1).p=sqrt(u.^2 + v.^2);
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
            
        end    
        
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
            
            savename=tit(1).tit;
            
            
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
        
    end
    
    maxC(i)=max(max(pdat(i).p));
    minC(i)=min(min(pdat(i).p));
    
    if logflag==1
        maxC(i)=log10(max(max(pdat(i).p)));
        minC(i)=log10(min(min(pdat(i).p)));
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
    half(i)=abs((izmax-izmin)/2);
    pend(i)=size(pdat(i).p,2);
    npend(i)=max([round(pend*1.25) pend+20]);
    nend(i)=max([round((npend-pend)*0.25) 10]);
    
    dumpint=timesTH(i).t(end)-timesTH(i).t(end-1);
    timesTH(i).t=[timesTH(i).t timesTH(i).t(end)+(dumpint.*(1:npend(i)-pend(i)))];    
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
    
    if iminovr(i)==1
        minCov=mincovOvr(i);
    end
    
    if imaxovr(i)==1
        maxCov=maxcovOvr(i);
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
    
    
    
    
    psame=repmat(pdat(i).p(:,pend(i)),[1 npend(i)-pend(i)-nend(i)-1]);
    
 %   pdat(i).p(:,pend(i)+1:npend(i)-nend(i)-1)=psame;
    
    mindat=minCov*(1-0.01*sign(minCov));   
    maxdat=maxCov*(1+0.01*sign(maxCov));
    
     pdat(i).p(1:half(i),npend(i)-nend(i):npend(i))=mindat;
     pdat(i).p(half(i)+1:end,npend(i)-nend(i):npend(i))=maxdat;



    
    if subplotting==0
        h(isub).h=subplot(a,b,isub);
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
    
    if icont==1
        [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts);
    else
        [cbfA(i).c]=pcolor(timesTH(i).t,height,pdat(i).p);shading interp;
    end
    
    if icont_extra==1
        hold on
        [c2,h2]=contour(timesTH(i).t(1:length(xinds)),add_ground_height+zz(i).z/1000 ...
        ,f/1000*sum(TwoDDan(idir).Q(izmin:izmax,xinds,4:6),3),[0:5:20],'w');
        set(h2,'linewidth',1.5);
        clabel(c2,h2,'color','w','fontsize',max([fsize-8 6]));
%        clabel(c2,h2,'labelspacing',72,'color','w');
    end
    
    
    
    if isquare==1
        axis square
    end
    
    if clab==1 & icont==1 %round up contour values so labels match those on colorbar
        
        ch=cbfA(i).c;
        jc=1;
        while jc<size(cbfA(i).c,2)
            %ch(1,jc)=str2num(num2str(10^(cbfA(i).c(1,jc)),'%2.2e'));
            
            if dlogflag==0 & logflag==1 
                if 10^cbfA(i).c(1,jc)>99
                    ch(1,jc)=sigfig(10^(cbfA(i).c(1,jc) ),3);
                else
                    ch(1,jc)=sigfig(10^(cbfA(i).c(1,jc) ),2);
                end
            elseif dlogflag==1
                %                    if idlog(cbfA(i).c(1,jc),dlogmin)>99
                %                         ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin ),2);
                %                     else
                %                         ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin ),2);
                %                     end
                ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin),3);
            else
                ch(1,jc)=sigfig(cbfA(i).c(1,jc),sig);
            end
            
            jc=jc+cbfA(i).c(2,jc)+1;
        end
        
        if length(ch)>0
            if manclab==0
                clabel(ch,cbfB(i).c,'labelspacing',144); %default 144
            else
                clabel(ch,cbfB(i).c,'manual');
            end
        end
    end %clab==1 & icont==1
    
    if icont==0
        cax=get(gca,'clim');
        if iminovr==1
            cax(1)=mincovOvr;
        end
        if imaxovr==1
            cax(2)=maxcovOvr;
        end
        caxis(cax);
        hc=colorbar;
        set(hc,'fontsize',fsize-3);
        
    end
    
    %hc=colorbarf(cbfA,cbfB);
    
    if clines==0
        shading flat; %shading stops black contour lines from appearing
    end
    
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
    
    if subplotting==1 & isub==1
        title(tit(i).tit,'fontsize',fsize-4);
    elseif subplotting==0
        title(tit(i).tit,'fontsize',fsize-4);
    end
    
    
    
    if length(cbfA(i).c)==0
        nocbar=1;
        normcbar=1;
    end
    
    if icont==0
        normcbar=1;
    end
    %caxis(h(i).h,[minCov maxCov*1.05]);
    
    if exist('idirs');
        if isub==length(idirs)
            %    if isamescale==1 &  isub==length(idirs)
            pos1=get(h(1).h,'position'); %[left bottom width height]
            posend=get(h(end).h,'position');
            
            height=pos1(2)-posend(2)+pos1(4);
            pos=[posend(1) posend(2) posend(3)+0.15 height];
        else
            pos=[0 0 1 1];
        end
    end
    
    if isamescale==1 &  isub==length(idirs)  
        if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            hc=colorbar( 'peer' , h(isub).h );
        elseif normcbar==1 & bigcbar==0 | (bigcbar==1 & isub==length(idirs)) | lememm==1
            axdan=axes('position',pos,'visible','off');
            %colbar=colorbar; %if colorbar already in place then colorbarf will replace it
            %set(colbar,'tag','Colorbar');
            hc=colorbarf(cbfA(i).c,cbfB(i).c);
            
            
        end
    else
        if nocbar~=1 & icont==1
            if normcbar==0 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
                hc=colorbarf(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!  
            elseif normcbar==0 &(bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
                axdan=axes('position',pos,'visible','off')
                %   colbar=colorbar;
                %   set(colbar,'tag','Colorbar');
                hc=colorbarf(cbfA(i).c,cbfB(i).c);
                
            end 
        end
        
        if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            hc=colorbar( 'peer' , h(isub).h );
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
        ctick=get(hc,'yticklabel');
        
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
        
        set(hc,'yticklabel',ctickstr);
        
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
            set(hc,'fontsize',fsize-2);
        end
        
        
        
        
    end
    
    ylims=get(h(isub).h,'ylim');   %re-scale to hide extra column put in to get the colorbars the same in both plots
    if iylim==1
        ylims=iylims;
    end
    axis(h(isub).h,[timesTH(i).t(1) timesTH(i).t(pend(i)) ylims]);
    
    if ixlim==1
        set(gca,'xlim',xlims);
    end
    
end



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

set(gcf,'paperpositionmode','auto');


if isave==1
    set(gcf,'paperpositionmode','auto');
    print(gcf,'-djpeg','-r350',exname);
    %print(gcf,'-dmeta',exname);
    %close(gcf);
end

if icolmap==1
    colormap(cmap);
end



