plotcase=50;

clear zz timesTH

if ~exist('idir'); idir=1; end
if ~exist('ieps'); ieps=0; end
if ~exist('isamescale'); isamescale=0; end
if ~exist('subplotting'); subplotting=0; end
if ~exist('onexlabel'); onexlabel=0; end      %flag to make it so that only the bottom plot has the xlabel on it


npes=npess2(idir);

%plots time height plots

icont=1;

if ieps==1
    fsize=12;
elseif subplotting==1
    fsize=9;
else
	fsize=15;
end

isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

icolmap=0; %flag to set to another colormap (defined by cmap)

iutc=1; %flag for time axis to be labelled as UTC (otherwise is labelled Local Time)
add_ground_height=0.62; %height to add to the vertical axis to account for level of ground abv msl.
%add_ground_height=0.0; %height to add to the vertical axis to account for level of ground abv msl.

minZ=0e3;
maxZ=25e3;  %19000;
maxtr=1.0;
%timesTH=Time;
hrstartles=18.67;
dumprange=[1:31];
dumprange=[1:49];
%dumprange=[1:39];
dumprange=[1:72];
%dumprange=[1:62]; %62

vectorf=0;

timesTH(1).t=(dumprange-1)*300/3600 + hrstartles;
jmax=5; %max no. plots on one screen

a1=1;
a2=2; %values for subplot(ai,b,i)
   
izmin=2;
izmax=2;
f=1e6*28.97/18;

iminovr=zeros([1 10]);
imaxovr=zeros([1 10]);
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
	tit(1).tit='Equilibrium Vapour Mixing Ratio (ppmv)';

    nplots2d=1;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    i2d=3; %tells it are labelling own x axis
    xlabelstr='dN/dt - Number of overshoots per month';
    
    minZ=15.8e3;
    maxZ=17e3;
    ncont=15;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
    
    notsame=1;
    
    izovr=1; %flag to set own z axis
    
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
    

    
    z=GridDan(1).Z;
    
    notsame=1;
    
    itimestamp=1;
    
        imaxovr=[1 0];
        maxcovOvr=8;
        iminovr=[1 0];
        %mincovOvr=2.5;
        mincovOvr=1.0;
        
        
case 577
    hrange=1;
    switch hrange
    case 1
        minZ=14e3;
        maxZ=23e3;
    case 2
        minZ=10e3;
        maxZ=27e3;
    case 3
        minZ=0.2e3;
        maxZ=23e3;
    case 4
        minZ=14e3;
        maxZ=19e3;    
	case 5
        minZ=14e3;
        maxZ=17e3; 
    case 6
        minZ=14e3;
        maxZ=30e3;  
    case 7
        minZ=0e3;
        maxZ=8e3;      
	end
    %ncont=15;
    
    
	fact=1e6*28.97/18;
	logflag=0;
    itimestamp=1;
    
    timlab=num2str(SER(end,1)/3600 + 19.75,3);
    
%      i577='potemp';
%     i577='lowtracer';
% %      i577='totwater';
%        i577='vapour';
% %  i577='tot_condensate';
% %    i577='si'; %supersat wrt ice
i577='temppert';
% %    i577='rhopert';
%   %  i577='ozone';
%  % i577='hydbal';
%   i577='icesatMR';
  
 
% i577='dpdz';
 % i577='rhog';

%  i577='rad';
%i577='lnb';
i577='vertvel';

idirstamp=1;
    
       imaxovr=0;
       iminovr=0;

  clines=1; %makes black contour lines appear
    clab=0;
    
    
    switch i577
    case 'vertvel'
        tit(1).tit='Vertical Velocity (m s^{-1})';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = 14.000000;
		maxcovOvr = 20.300000;
        clab=0;
        clines=0;
        icont=0;
        
    case 'lnb'
        tit(1).tit='Level of Neutral Buoyuancy (km)';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr = 14.000000;
		maxcovOvr = 17.00000;
        clab=0;
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
        
        imaxovr=1;
        iminovr=1;
        
        mincovOvr = 330.000000;
		maxcovOvr = 474.000000;
        
        maxcovOvr=750;
        
      %  maxcovOvr= 400;
        
%        mincovOvr = dlog(330.000000,dlogmin);
%		maxcovOvr = dlog(474.000000,dlogmin);
        
        clab=0;
        clines=1;
        ncont=75;
        
        ncont=100;
    case 'wind'
		tit(1).tit='Vertical Wind Speed (m/s)';
    case 'htracer'
        tit(1).tit='Height Dependent Tracer';
    case 'vapour'
		tit(1).tit='Vapour Mixing Ratio (ppmv)';
        imaxovr=1;
        iminovr=1;
        
        mincovOvr=1;
        maxcovOvr=8.5;
        
%         maxcovOvr=25;
         maxcovOvr=15;
% 
%         mincovOvr=3;        
%        % maxcovOvr=5.5;
       
       dlogflag=0;

        
%         mincovOvr = dlog(3.600000,dlogmin);
% 		maxcovOvr = dlog(540.000000,dlogmin);

        
        clines=0;
        clab=0;
        
        vectorf=0;
        
    case 'ozone'
		tit(1).tit='Ozone Mixing Ratio (ppmv)';
        imaxovr=0;
        iminovr=0;
        
        mincovOvr=8;
        maxcovOvr=60;
        
    case 'lowtracer'
		tit(1).tit='Low level tracer (g kg^{-1})';
        clines=0;
        
    case 'totwater'
		tit(1).tit='Total Water Mixing Ratio (ppmv)';
        imaxovr=1;
        iminovr=1;
        
        mincovOvr=1;
        maxcovOvr=8.5;
        
       % mincovOvr = dlog(3.600000,dlogmin);
	%	maxcovOvr = dlog(540.000000,dlogmin);
        
  %      dlogflag=1;
        dlogmin=1e-2;
        
%        maxcovOvr=15;
        
        ncont=25;
        
        clines=0;
        
    case 'tot_condensate'
		tit(1).tit='Total Condensate (ppmv)';
        clines=0;
        dlogflag=1;
        dlogmin=1e-2;
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(3.600000,dlogmin);
		maxcovOvr = dlog(540.000000,dlogmin);
        
     case 'temppert'
		tit(1).tit='Temperature Perturbation (K)';
		tit(1).tit='Temperature (K)';
        
        imaxovr=0;
        iminovr=0;
        
        mincovOvr=0;
        maxcovOvr=20;
        
        clines=0;
        
    case 'rhopert'
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
        
     end
    
    
   % tit(2).tit='Ice Number Concentration (mg^{-1})';
    nplots2d=1;
    
  
    
    i2d=2; %tells it to label x axis in km
    
    
 
    
    z=GridDan(1).Z;
    
    notsame=1;
    

%    iminovr=[0 0];
%    imaxovr=[0 0];
        
    
    figlab=[tit(1).tit ' 2d plot '];
    savename=figlab;
    
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
    
    dlogflag=0;
	logflag=0;
    fact=1e6*28.97/18;
    
    iminovr=1;
    mincovOvr=-0.1;
    %mincovOvr=3.8;
    
    imaxovr=1;
    maxcovOvr=0.1;
    
    iflux=0;
    
    
    
%	tit(1).tit='Mean Vapour (ppmv)';
%    tit(1).tit='Mean Total Water (ppmv)';
%    tit(1).tit='Mean Lower Tracer Mixing Ratio (ppmv)';
%     tit(1).tit='Microphysical Ice Number Source (kg^{-1} s^{-1})';
%     tit(1).tit='Microphysical Ice Source (kg/kg s^{-1})';
%    tit(1).tit='Min Vapour (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	%tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';

hrange=11;
switch hrange
	case 1
        minZ=13e3;
        maxZ=22e3;
	case 2
        minZ=15.3e3;
        maxZ=19e3;
      %  maxZ=17e3;
    case 3
        minZ=15.3e3;
        maxZ=18e3;
    case 4
        minZ=0.2e3;
        maxZ=22e3;
    case 5
        minZ=0e3;
        maxZ=30.5e3;
    case 6
        minZ=15.3e3;
        maxZ=30e3;
    case 7
        minZ=10e3;
        maxZ=22e3;
    case 8
        minZ=14e3;
        maxZ=19e3;
    case 9
        minZ=0e3;
        maxZ=20e3;    
	case 10
        minZ=13e3;
        maxZ=18.6e3;
    case 11
        minZ=15e3;
        maxZ=30.4e3;    
end

    
s50='adrate';
s50='fallrate';
s50='change';
%s50='topdowncum';
%s50='fallflux';
s50='icedist';
%s50='icendist';

% % s50='vapad';
% %s50='micronc';
% %s50='fallnc';
%s50='adnc';
%s50='changenc';
%s50='icead';
%s50='changevap';
%s50='meanvap';
%s50='changeice';
%s50='icemass';
%s50='minvap';
%s50='mintot';
% % %s50='fall+ad';
%s50='microice';
% % s50='microvap';
s50='rain';
s50='liq';
% % s50='ice';
% 
% s50='allice';
% 
% % % % s50='snow';
%s50='graupel';
% % % % 
 s50='maxw';
% s50='minw';
% % % 
% % % % s50='snowno';
%s50='iceno'
% % % % s50='grano';
% % % % 
% % % % s50='minvap';
% % % 
% % % s50='PGMLT';
% % % s50='racw';
% % % s50='praut';
% % % s50='pifrw';
% % % s50='allpr';
% % % %s50='piacw';
% %  s50='pidep';
%s50='pIsub'; %sublimation of all ice
% %s50='pisub';
% % % s50='dqi';
% % % s50='pgsub';
% % % %s50='prevp';
%s50='dql';
% % s50='pcond';
% %  s50='dq_potemp';
% %  s50='dq_non';
%s50='dqtot';
%s50='dqvap';
%s50='nntot';
%s50='nnvap';
% %   s50='ratio_potemp';
% %   s50='change_conv_potemp';
% %   s50='combined_potemp';
%s50='low_tracer';
%s50='maxlowtracer';
% 
% s50='iceadcum';
% s50='icefallcum';
% s50='icemicrocum';
% % s50='picesubcum';
%s50='totadcum';
%s50='vapadcum';

% s50='si'; %supersat wrt ice
% s50='rhopert'; 
% s50='drhodz';
% 
% s50='upflux';
% s50='meanw';

%s50='lnbbel';
%s50='lnbdist';
%s50='lnbdist_tot';
%s50='lnbdist_vap';
%s50='rad';
%s50='lwrad';
%s50='swrad';
%s50='minTguess';
s50='vapdist';
s50='totdist';


    dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
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
        
    case 'icedist'
        iflux=0;
        
        dlogflag=1;
        dlogmin=0.5e-10;
        dlogmin=1e-3;
        
        tit(1).tit=['Total Ice dq dD^{-1} (ppmv micron^{-1}) at ' num2str(GridDan(idir).Z(iz)/1000+add_ground_height,4) ' km'];

		iminovr=1;
		imaxovr=1;
        
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
        tit(1).tit='Mean LNB for negatively buoyant air (km)';
        tit(1).tit='Minimum LNB (km)';

		iminovr=1;
		imaxovr=1;
        
        mincovOvr = GridDan(1).Z(100)/1000+add_ground_height;   %GridDan(1).Z(105)/1000+add_ground_height;
		maxcovOvr = GridDan(1).Z(160)/1000+add_ground_height;
      
        clines=1;
        clab=1;
        
        case 'dqtot'
        iflux=0;        
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Sum of Total Water Points Deficit Below 5 ppmv x Grid Length (ppmv km)';

		iminovr=0;
		imaxovr=0;
        
        mincovOvr = 0.000000;
		maxcovOvr = 150.000000;
        
        mincovOvr = 0.000000;
		maxcovOvr = 240.000000;
        
        clines=0;
        clab=0;
        
        case 'dqvap'
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
		maxcovOvr = 400.000000;
         
        clines=0;
        clab=0;
        
        case 'nntot'
        iflux=0;        
        dlogflag=0;
        dlogmin=1e-2;

		iminovr=1;
		imaxovr=1;
        
        tit(1).tit='Total 2-D Length of Total Water Points Below 5 ppmv (km)';
        
        mincovOvr = 0.000000;
		maxcovOvr = 150.000000;
        
        clines=0;
        clab=0;
        
        case 'nnvap'
        iflux=0;
        dlogflag=0;
        dlogmin=1e-2;
        
        tit(1).tit='Total 2-D Length of Vapour Points Below 5 ppmv (km)';

		iminovr=1;
		imaxovr=1;
        
		mincovOvr = 0.000000;
		maxcovOvr = 250.000000;
      
        clines=0;
        clab=0;
        
        
        case 'meanw'
        iflux=0;
        
        dlogflag=0;
        dlogmin=1e-2;
        
%        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean Updraught (ms^{-1})';

		iminovr=1;
		imaxovr=0;
        
        mincovOvr = 0.000000;
		maxcovOvr = 0.870000;
         
%         mincovOvr = dlog(0.000000,dlogmin);
%         maxcovOvr = dlog(0.420000,dlogmin);
        
      %  ncont=25;
      
        clines=1;
        clab=1;
        
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
        tit(1).tit='Mean Density Change for Points with Total Water Below 5 ppmv (kg m^{-3})';
        tit(1).tit='Mean Density Change (kg m^{-3})';

		iminovr=0;
		imaxovr=0;
        
        mincovOvr = 0;
		maxcovOvr = 0.05;
         
%        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        ncont=12;
       % ncont=25;
      
        clines=1;
        clab=0;
        
        sig=2;
        
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
        dlogflag=1;
        dlogmin=1e-2;
        
%        tit(1).tit='Ratio of Non-reversbile to reversbile changes';
        tit(1).tit='Mean Low Tracer (kg kg^{-1})';

		iminovr=1;
		imaxovr=1;
        
        mincovOvr = -0.7;
		maxcovOvr = 0.05;
         
        mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.420000,dlogmin);
        
        mincovOvr = dlog(0.000000,dlogmin);
		maxcovOvr = dlog(0.690000,dlogmin);
        
        mincovOvr = dlog(-0.000027,dlogmin);
		maxcovOvr = dlog(0.000061,dlogmin);
        
        
        
        ncont=25;
        clines=1;
        clab=1;
        
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
		imaxovr=0;
        
         mincovOvr = 0.000000;
		maxcovOvr = 0.100000;
        
        ncont=12;
        clines=0;
        
        case 'pIsub'
        
       
        dlogflag=1;
        dlogmin=1e-5;
        
        tit(1).tit='Sublimation of all ice (ppmv s^{-1})';
		iminovr=1;
		imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
		maxcovOvr = dlog(0.051000,dlogmin);
        
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
        dlogmin=1e-2;
        
        tit(1).tit='Deposition of vapour onto ice (ppmv s^{-1})';
		iminovr=1;
		imaxovr=1;
        
         mincovOvr = 0.000000;
		 maxcovOvr = 0.600000;
         
         mincovOvr = dlog(0.010000,dlogmin);
		maxcovOvr = dlog(0.610000,dlogmin);
        
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
		iminovr=1;
		imaxovr=1;
        
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
        
        tit(1).tit='Minimum Vapour Mixing Ratio (ppmv)';
        tit(1).tit='Maximum Vapour Mixing Ratio (ppmv)';
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
        dlogmin=1;
        
        tit(1).tit='Ice Number Concentration (kg^{-1})';
		iminovr=1;
		imaxovr=1;
        
        mincovOvr = dlog(3.800000,dlogmin);
maxcovOvr = dlog(269999999.999999,dlogmin);
        

        
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
        
        tit(1).tit='Maximum updraught (ms^{-1})';
		iminovr=0;
		imaxovr=0;

		mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(480.000000,dlogmin);
        
        mincovOvr = 0.0;
		maxcovOvr = 44.000000;
		
		
		clines=0;
		clab=0;
		
		ncont=12;
		ncont=25;

case 'minw'
        dlogflag=0;
        dlogmin=1;
        
        tit(1).tit='Maximum downdraught (ms^{-1})';
		iminovr=0;
		imaxovr=0;

		mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(480.000000,dlogmin);
        
        mincovOvr = 0.0;
		maxcovOvr = 44.000000;
		
		
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
        dlogflag=1;
        dlogmin=10;
        
        tit(1).tit='Ice mixing ratio (ppmv)';
		iminovr=1;
		imaxovr=1;

		mincovOvr = dlog(0.000000,dlogmin);
		maxcovOvr = dlog(670.000000,dlogmin);
        
 %  mincovOvr = 0.000000;
 %  maxcovOvr = 460.000000;
 clab=1;

     case 'allice'
        dlogflag=1;
        dlogmin=1e-3;
        
        tit(1).tit='Ice mixing ratio (ppmv)';
		iminovr=0;
		imaxovr=0;

		mincovOvr = dlog(0.000000,dlogmin);
		maxcovOvr = dlog(670.000000,dlogmin);
        
 %  mincovOvr = 0.000000;
 %  maxcovOvr = 460.000000;
 clab=1;
 
    case 'rain'
        dlogflag=0;
        dlogmin=1e-5;
        
        tit(1).tit='Rain mixing ratio (ppmv)';
		iminovr=1;
		imaxovr=1;

		mincovOvr = dlog(0.000000,dlogmin);
        maxcovOvr = dlog(0.000120,dlogmin);
        
   mincovOvr = 0.000000;
maxcovOvr = 250.000000;
    
    case 'liq'
        dlogflag=0;
        dlogmin=1;
        
        tit(1).tit='Mean Liquid water mixing ratio (g kg^{-1})';
		iminovr=1;
		imaxovr=1;

		mincovOvr = dlog(0.000000,dlogmin);
		maxcovOvr = dlog(0.1000000,dlogmin);
        
        
        
mincovOvr = 0.000000;
maxcovOvr = 0.100000;
    
    ncont=15;
    
    clines=0;
    clab=0;
        
    case 'minvap'
        dlogflag=0;
        dlogmin=0.5;
        
        tit(1).tit='Minimum vapour mixing ratio (ppmv)';
		iminovr=1;
		imaxovr=1;

		mincovOvr=dlog(0.92,dlogmin);		
		maxcovOvr=dlog(34,dlogmin);
        
        mincovOvr = 0.770000;
		maxcovOvr = 5.200000;
        
    case 'mintot'
        dlogflag=0;
        dlogmin=2;
        
        tit(1).tit='Minimum total water mixing ratio (ppmv)';
		
        iminovr=1;
        imaxovr=1;
        
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
        dlogmin=1e-1;
        
        dlogmin=0.1;
                
        tit(1).tit='Change in Total Water (ppmv)';
        
        iminovr=0;
        imaxovr=0;
        
        mincovOvr = -1.200000;
        maxcovOvr = 10.000000;
        
        mincovOvr = dlog(-0.460000,dlogmin);
		maxcovOvr = dlog(500.000000,dlogmin);
        
        mincovOvr = dlog(-0.183600,dlogmin);
		maxcovOvr = dlog(203.000000,dlogmin);
        
      %  mincovOvr = dlog(-0.480400,dlogmin);
	%	maxcovOvr = dlog(568.000000,dlogmin);
        
        
        %mincovOvr = 0.000000;
        %maxcovOvr = dlog(10,dlogmin);

        ncont=35;
        
        clines=1;
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
        
        clab=0;
        clines=0;
        
        manclab=0;
        
        iflux=1;
        
        icolmap=1;
        cmap=hsv;
        
    case 'meanvap'
        dlogflag=1;
        dlogmin=1e-3;
        
        tit(1).tit='Mean Vapour (ppmv)';
        iminovr=0;
        mincovOvr=dlog(-0.82,dlogmin);
    
        imaxovr=0;
        maxcovOvr=dlog(0.96,dlogmin);
        
        clab=1;
        
        manclab=0;
        
        iflux=1;
        
    case 'changeice'
        iflux=1;
        dlogflag=1;
        dlogmin=1e-1;
         %dlogmin=1e-4;
        
        tit(1).tit='Mean Total Ice (ppmv)';
        %tit(1).tit='Mean Total Ice (g kg^{-1})';
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(0.000000,dlogmin);
		maxcovOvr = dlog(600.000000,dlogmin);
        
        mincovOvr = dlog(-0.183600,dlogmin);
		maxcovOvr = dlog(603.000000,dlogmin);
        
        clab=1;
        clines=0;
        
        manclab=0;   
        
        ncont=35;
        
        icolmap=1;
        cmap=hsv;
        
    case 'icemass'
        iflux=1;
        dlogflag=1;
        dlogmin=1e-3;
        %dlogmin=1e-1;
        
        tit(1).tit='Mean Ice Particle Mass (ppmv)';
        
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
        
        dlogflag=1;
        dlogmin=1e-2;        
        
        iminovr=1;
        imaxovr=1;
        
        mincovOvr = dlog(-0.02200,dlogmin);
		maxcovOvr = dlog(0.38000,dlogmin);
        

        
        switch hrange
			case 1
                dlogflag=1;
                dlogmin=1e-4;
                mincovOvr = dlog(-0.002200,dlogmin);
				maxcovOvr = dlog(0.038000,dlogmin);

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
        
        clab=0;
        clines=0;
        
    case 'vapad'
        dlogflag=1;
        dlogmin=1e-4;
        
        tit(1).tit='Advective Source of Vapour (ppmv s^{-1})';
        
        iminovr=1;
        imaxovr=1;
        
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
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    if iutc==1
        time=GridDan(idir).t+3;
    else
        time=GridDan(idir).t;
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
    
    iylim=1;
    iylims=[0.62 22.69];
    
    ixlim=1;
    xlims=[GridDan(1).t(1)+3 GridDan(1).t(18)+3];
    
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
	posit=[9 50 scrsz(3)/1.2 scrsz(4)/1.7];
case 'lacieLap'
	posit=[9 50 scrsz(3)/1.46 scrsz(4)/2.07];
end

if subplotting==1
    posit=[9 50 scrsz(3)/2.4 scrsz(4)/1.5];
    
    if nsub==1; hf=figure('position',posit,'name',figlab); end
    h(idir).h=subplot(xsub,ysub,nsub); 
    
%     for ih=1:length(h)
%         posh=get(h(ih).h,'position');
%         set(h(ih).h,'position',[posh(1) posh(2) posh(3)-0.2 posh(4)]);
% 	end
    
    posh=get(h(idir).h,'position');
    set(h(idir).h,'position',[posh(1) posh(2) posh(3)-0.2 posh(4)]);
            
    
    
else
	hf=figure('position',posit,'name',figlab);
    h(idir).h=subplot(1,1,1); 
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
        isub=idir;
    else
        isub=i;
    end
    
    if (exist('time') & i2d~=2 & i2d~=1 & i2d~=3); timesTH(i).t=time(dumprange); end
    
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
        zz(i).z=z(z0:zend)-620; %620 added later
        izmin=1;
        izmax=length(zz(1).z);
        
        pdat(1).p=qq;
        imaxovr=[1];
        maxcovOvr=5;
        iminovr=[1];
        mincovOvr=3;
        
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
        
        xinds=[1:200];
        %xinds=[160:230];
        
        %xinds=[130:260];
        
        %xinds=[300:430];

       %xinds=[1:1000];
       
       ix=findheight(GridDan(idir).Y1,-477.5e3);
       ix2=findheight(GridDan(idir).Y1,-445.8e3);

		xinds=[ix:ix2];
       
       xinds=[3900:4000 1:500];
       
       
       
       %xinds=1:length(GridDan(idir).Y1);
        
        [izmin izmax]=findheight(GridDan(1).Z,minZ,maxZ);
%        timesTH(1).t=GridDan(1).Y1'/1000;
        zz(1).z=GridDan(idir).Z(izmin:izmax);
        timesTH(1).t=GridDan(idir).Y1(1:length(xinds))'/1000;
        %pdat(1).p=squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[13]),3)); %height dependent tracer
%        pdat(1).p=squeeze(sum(TwoDDan(1).W(izmin:izmax,:),3));  %total water
%                pdat(1).p=fact*squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[1]),3)); %vapour
        
        switch i577
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
            pdat(1).p=TwoD.W(izmin:izmax,xinds); %vertical velocity
            
            %vertvel=TwoD(idir).W;
            %exdirB=[direcDan(jc).dir 'results/diags/vertvel/'];
            %exname=strcat(exdirB,'vertvel-',int2str(jj),'.mat');
            %save(exname,'vertvel');
        case 'ozone'
            pdat(1).p=TwoDDan(idir).Q(izmin:izmax,xinds,15); 
        case 'potemp'
            tref=repmat(GridDan(idir).THREF(izmin:izmax),[1 length(GridDan(idir).Y1(xinds))]);
            pdat(1).p=squeeze(sum(TwoD.TH1(izmin:izmax,xinds),3))+tref; %potemp
        case 'wind'
            pdat(1).p=squeeze(sum(TwoDDan(idir).W(izmin:izmax,:),3)); %vertical velocity
        case 'lowtracer'                
            pdat(1).p=squeeze(sum(TwoDDan(idir).Q(izmin:izmax,xinds,[13]),3)); %low tracer
        case 'totwater'                
            pdat(1).p=f*squeeze(sum(TwoD.Q(izmin:izmax,xinds,[1:6]),3)); %total water
        case 'vapour'                
            pdat(1).p=f*squeeze(sum(TwoD.Q(izmin:izmax,xinds,[1]),3)); %vapour
        case 'tot_condensate'                
%            pdat(1).p=1000*squeeze(sum(TwoDDan(idir).Q(izmin:izmax,:,[2:6]),3)); %tot condensate
            pdat(1).p=f*squeeze(sum(TwoD.Q(izmin:izmax,xinds,[2:6]),3)); %tot condensate
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
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]); %ref potemp
            T=squeeze(sum(TwoD.TH1,3))+tref; %potemp
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            P=TwoD.PP; %actual p
            T=T./(1e5./P).^0.286; %actual T
            tref=tref./(1e5./pref).^0.286; %ref temp
            
            Tp=T-tref; %perturbation of temperature
            
            pdat(1).p=T(izmin:izmax,xinds);
            
            tpertTimH(1).t(:,jj)=mean(Tp,2); %mean temp pert
            
        case 'rhopert'
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

rho=repmat(GridDan(idir).RHON(izmin:izmax),[1 length(dumprange)]);

   switch s50
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
       
   case 'icedist'       
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
		pdat(i).p=sum(icediagsRAD(idir).i(izmin:izmax,dumprange,[1 2]),3); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
   case 'swrad'       
		pdat(i).p=sum(icediagsRAD(idir).i(izmin:izmax,dumprange,[2]),3); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
   case 'lwrad'       
		pdat(i).p=sum(icediagsRAD(idir).i(izmin:izmax,dumprange,[1]),3); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more     
   case 'lnbdist' 
       ilnbmin=findheight(bins_tot(idir).b,minZ/1000+0.62);
       ilnbmax=findheight(bins_tot(idir).b,maxZ/1000+0.62);

       zz(1).z=1000*(bins_tot(idir).b(ilnbmin+1:ilnbmax)+bins_tot(idir).b(ilnbmin:ilnbmax-1))/2 - 1000*add_ground_height; %since for plotting adds ground height and /1000
       pdat(i).p=lnbbins_neg_tot(idir).l(ilnbmin:ilnbmax-1,dumprange);
  
   case 'lnbbel' 
        zref=repmat(GridDan(1).Z(1:size(meanlnb_abv(idir).m,1))/1000+add_ground_height,[1 size(meanlnb_abv(idir).m,2)]);
		%pdat(i).p=meanlnb_abv(idir).m(izmin:izmax,dumprange)-zref(izmin:izmax,dumprange);
        pdat(i).p=meanlnb_bel(idir).m(izmin:izmax,dumprange);
        pdat(i).p=minlnb(idir).m(izmin:izmax,dumprange);
    case 'dqtot'       
		pdat(i).p=length(GridDan(idir).Y1)*( dq_tot(idir).d(izmin:izmax,dumprange,2) ) *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
        %was mulitplied by length of .Y1 in Allimp...
    case 'nntot'       
		pdat(i).p=nn(idir).n(izmin:izmax,dumprange)*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
    case 'nnvap'       
		pdat(i).p=nn2(idir).n(izmin:izmax,dumprange)*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more    
   case 'dqvap'       
		pdat(i).p=length(GridDan(idir).Y1)*( dq_vaps(idir).d(izmin:izmax,dumprange,2) ) *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
        
   %     ipps=1;
	%	pdat(i).p= ( length(GridDan(idir).Y1)*( dq_vaps(idir).d(izmin:izmax,dumprange,ipps) )...
    %        ./nn2(idir).n(izmin:izmax,dumprange,ipps) ); %multiply by dy so is in ppmv*km since otherwise high res will mean there are more        
   case 'rhopert'
		pdat(1).p=rhopertTimH(idir).t(izmin:izmax,dumprange); 
		%pdat(1).p=rho5ppmv_vappos(idir).r(izmin:izmax,dumprange); 
		pdat(1).p=rho5ppmv_totpos(idir).r(izmin:izmax,dumprange); 
		%pdat(1).p=low_prctiles(idir).t(izmin:izmax,dumprange,21); 
                    
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
   case 'si'
%       pdat(1).p=simaxTimH(idir).s(izmin:izmax,dumprange);
         pdat(1).p=simean(idir).s(izmin:izmax,dumprange);
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
       
       pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[151]),3);

      pdat(1).p=adlow(izmin:izmax,dumprange);
       
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
   case 'PGMLT'
        pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[2]),3);       
   case 'minvap'
       pdat(1).p=f*vap_prctiles(idir).t(izmin:izmax,dumprange,end);
   case 'grano'
       pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[44]),3);
   case 'iceno'
       pdat(1).p=sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43]),3);
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
       pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[42]),3);
       
   case 'allice'
       pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40:42]),3);
       
   case 'liq'
       pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[38]),3);
       %pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[224]),3);
      %pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[83]),3);

   case 'rain'
       pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[39]),3);
       %pdat(1).p=fact*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[84]),3);
    case 'minvap' %min vapour
        pdat(i).p=fact*squeeze(vap_prctiles(idir).t(izmin:izmax,dumprange,1)); 
    case 'mintot' %min tot water
        pdat(i).p=fact*squeeze(tot_prctiles(idir).t(izmin:izmax,dumprange,1)); %min total water
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
	case 'topdowncum'
        pdat(1).p=topdown;           
    case 'changerate'
        pdat(1).p=changerate; 
   case 'icemass'
        ice=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[42]),3) )/npes;
        icenc=(sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43]),3) )/npes;
        pdat(1).p=changeice./changenc*1e3;      
        pdat(1).p=ice./icenc*1e3;      

	case 'microice'
        pdat(1).p=microicerate;
    case 'vapad'
        pdat(1).p=-vapad;
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

    
    
    
    psame=repmat(pdat(i).p(:,pend(i)),[1 npend(i)-pend(i)-nend(i)-1]);
    
    pdat(i).p(:,pend(i)+1:npend(i)-nend(i)-1)=psame;
    
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
    if icont==1
        [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,add_ground_height+zz(i).z./1000,pdat(i).p,conts);
    else
        [cbfA(i).c]=pcolor(timesTH(i).t,add_ground_height+zz(i).z./1000,pdat(i).p);shading interp;
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
                    ch(1,jc)=sigfig(cbfA(i).c(1,jc),2);
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

    %hc=colorbarf(cbfA,cbfB);
    
    if clines==0
        shading flat; %shading stops black contour lines from appearing
    end
	
    %xti=set(h(i).h,'xticklabels',timestxt);
    
    set(h(isub).h,'fontsize',fsize);
    
    if onexlabel==0 | (isub==nplots2d & onexlabel==1 & subplotting==0) | (isub==idirs(end) & onexlabel==1 & subplotting==1)
        xlabel(xlabelstr);
    end
        
    if izovr~=2
		ylabel('Height (km)');
    else
        ylabel(ylabelstr);
    end
    
    if subplotting==1 & isub==1
        title(tit(i).tit);
    elseif subplotting==0
        title(tit(i).tit);
    end
    
    

    if length(cbfA(i).c)==0
        nocbar=1;
        normcbar=1;
    end
    
    if icont==0
        normcbar=1;
    end
    %caxis(h(i).h,[minCov maxCov*1.05]);
    
    if isamescale==1 &  isub==idirs(end)
        pos1=get(h(1).h,'position'); %[left bottom width height]
        posend=get(h(end).h,'position');
        
        height=pos1(2)-posend(2)+pos1(4);
        pos=[posend(1) posend(2) posend(3)+0.15 height];
    else
        pos=[0 0 1 1];
    end
 
  if isamescale==1 &  isub==idirs(end)  
    if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
            hc=colorbar( 'peer' , h(isub).h );
    elseif (normcbar==1 & bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
            axdan=axes('position',pos,'visible','off');
            %colbar=colorbar; %if colorbar already in place then colorbarf will replace it
            %set(colbar,'tag','Colorbar');
            hc=colorbarf(cbfA(i).c,cbfB(i).c);
            
           
    end
  else
        if nocbar~=1 & icont==1
            if bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
                hc=colorbarf(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!  
            elseif (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
                axdan=axes('position',pos,'visible','off')
             %   colbar=colorbar;
             %   set(colbar,'tag','Colorbar');
                hc=colorbarf(cbfA(i).c,cbfB(i).c);
            
            end 
        end
  end
    
    if normcbar==0
    
      if (logflag==1 | dlogflag==1 ) & nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
      %re-label colorbar ticks if have log settings  
        clear ctickstr
        
        ctick=get(hc,'yticklabel');
        ctickstr(1,1)=' ';
        for j=2:length(ctick)-1
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
                    if iszero(j-1)==1  %if is suppossed to be exactly zero then set to zero
                        te='0.0';
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

      elseif (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) ) %if (logflag==1 | dlogflag==1 ) & nocbar=...
          %re-label colorbar if not log plot
        clear ctickstr
        ctick=get(hc,'yticklabel');
        ctick2=ctick;
        clear ctick
        for j=2:length(ctick2)-1
            ctick(j-1)=str2num(ctick2{j});
        end
        ctickstr(1,:)=' ';
        for j=2:length(ctick)+1
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            %te=num2str(ctick(j),'%2.2e');
            te=num2str(sigfig(ctick(j-1),sig));
            ctickstr(j,1:length(te))=te;
        end
        
        set(hc,'yticklabel',ctickstr);
        
%         add=ctick(end)/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
    
    
        
   
        
        %set(hc,'fontsize',fsize-6);
        
     
    
      %else
    
    
    
    if logflag==1 & nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
        clear ctickstr
        ctick=get(hc,'ytick');
        for j=1:length(ctick)
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            %te=num2str(10^ctick(j),'%2.2g');
            if nu>99
                te=num2str(sigfig(10^ctick(j) - offset,2));
            else
                te=num2str(sigfig(10^ctick(j) - offset,2));
            end
            
            ctickstr(j,1:length(te))=te;
        end
        
        set(hc,'yticklabel','');
        
        add=ctick(end)/50;
        
        set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize);
        
    
    
    elseif nocbar==0 & icont==1 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
        clear ctickstr
        ctick=get(hc,'ytick');
        for j=1:length(ctick)
            nu=ctick(j);
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            if nu>99
                %te=num2str(ctick(j),'%2.2e');
                te=num2str(sigfig(ctick(j),2));

            else
               % te=num2str(ctick(j),'%2.2e');
                te=num2str(sigfig(ctick(j),2));
            end
            ctickstr(j,1:length(te))=te;
        end
     end
     
    
 end
    
    if nocbar==0 & subplotting==1 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
        set(hc,'fontsize',fsize);
    elseif nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==idirs(end)) )
        set(hc,'fontsize',fsize-6);
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
     text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-2.5 11.5],'fontsize',fsize);
    
    

   % text(timesTH(i).t(1)*1.2,23.0,['Time = ' f ' UTC'],'fontsize',18);
end

if idirstamp==1
%    text(timesTH(i).t(1)-0.5,ylims(2)*1.002,[direcDan(idir).dir],'fontsize',12);
dist=(ylims(2)-ylims(1))/10;
dist2=(timesTH(i).t(pend(i))-timesTH(i).t(1))/10;
%    text(timesTH(i).t(1)-0.5,ylims(1)*0.988,[direcDan(idir).dir],'fontsize',12);    
    axes(h(end).h);
    text(timesTH(i).t(1)-dist2,ylims(1)-dist,[runName(idir).nam],'fontsize',fsize-2);    

end


% 
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');

if (i2d~=1 & i2d~=2 & i2d~=3)
    xx=get(h(idir).h,'xticklabels');
    xx=str2num(xx);
    xx=num2str(mod(xx,24));
    set(h(idir).h,'xticklabels',xx);
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



