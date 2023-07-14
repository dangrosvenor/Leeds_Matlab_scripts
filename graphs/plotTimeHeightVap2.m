clear zz timesTH

%plots time height plots
fsize=18;
isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

minZ=0e3
maxZ=25e3;  %19000;
maxtr=1.0;
%timesTH=Time;
hrstartles=18.67;
dumprange=[3:86];
timesTH(1).t=(dumprange-1)*300/3600 + hrstartles;
jmax=5; %max no. plots on one screen
a1=1;
a2=2; %values for subplot(ai,b,i)
izmin=2;
izmax=2;

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


normcbar=0;

dumpint=300; %dump interval in seconds

i2d=0;
izovr=0; %flag to say are setting own z axis
itimelab=0; %flag to say that x axis is time and so should be restricted to <24
figlab='2d contour plot';
iylim=0; %flag to override the setting of the y part of axis (e.g. so can have axis larger than that plotted)

nplots2d2d=1;

ncont=25;
clab=1; %flag to label contours
i=1;

phase=1;


clear h max timestxt pdat minC maxC



% for i=1:length(times)
%     te=num2str(times(i),3);
% 	timestxt(i,1:length(te))=te;
% end
% tit(1).tit='Max of Low Level Tracer';
% tit(2).tit='Max of Low Level Tracer'
logflag=0;


plotcase=48;

switch(plotcase)
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
    
    case 57
	fact=1e6*28.97/18;
	logflag=0;
	tit(1).tit='Tot Water Mixing Ratio (ppmv)';
	tit(2).tit='Ice Number Concentration (mg^{-1})';
    nplots2d=2;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    i2d=2; %tells it to label x axis in km
    
    minZ=0e3;
    maxZ=30e3;
    ncont=15;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
    z=GridDan(1).Z;
    
    notsame=1;
    
    itimestamp=1;
    
    
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
	tit(1).tit='Mean Vapour Change (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	%tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='Mean Tot Water time-height';
    
    minZ=14.7e3;
    maxZ=22e3;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    %i2d=2; %tells it to label x axis in km
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
    sig=1;
    
    case 50
        
    dumprange=[3:86];    
	logflag=0;
    fact=1e6*28.97/18;
    
    iminovr=1;
    mincovOvr=2;
    
    imaxovr=1;
    maxcovOvr=6;
	tit(1).tit='Mean Vapour (ppmv)';
  %  tit(1).tit='Min of Total Water (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	%tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='5th per tot  Vapour time-height';
    
%    minZ=12e3;
    minZ=14.7e3;
    maxZ=22e3;
    
    sig=2;
    clines=1; %makes black contour lines appear
    clab=1;
    
    ncont=15;
    
    %i2d=2; %tells it to label x axis in km
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
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
    clab=1;
    figlab='Min W TimH';
	
    
    %iminovr=[1 1];
    mincovOvr=[-80];
    
    z=GridDan(idir).Z;
    time=GridDan(idir).t+3;
    
    minZ=0e3;
    maxZ=max(z);
    
    iylim=1;
    iylims=[0.62 0.62+GridDan(2).Z(end)/1000];
    
    nplots2d=1;
    
    tit(1).tit=['Max Updraught (m/s)'];
    
    itimelab=1; 

    
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
    
    case 43
	fact=1;
	logflag=1;
    clab=1;
    figlab='Ice MR';
    ncont=15;
	
    
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
    
    case 35
	fact=1;
	logflag=0;
    clab=1;
    figlab='Ratio Ice NC';
	tit(1).tit='Ratio for Ice NC';
	tit(2).tit='';
    MI0=1e-15;
    
    ncont=20;
    
    dumprange=[1:64];
    
    case 34
	fact=1;
	logflag=0;
    clab=1;
    figlab='Ratio Ice MR';
	tit(1).tit='Ratio for Ice MR';
	tit(2).tit='';
    
    %imaxovr=1;
    %maxcovOvr=1e10;
    
    ncont=20;
    
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
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
hf=figure('position',posit,'name',figlab);

%nplots2d=length(prof);


    if nplots2d==1
        a=a1;
    else
        a=a2;
    end
    b=ceil(min(nplots2d,jmax)/2);
    
phase=2;
for  i=1:nplots2d
    
    if (exist('time') & i2d~=2 & i2d~=1); timesTH(i).t=time(dumprange); end
    
	%exname=strcat('c:/matlabr12/work/bauru/tracersjan2005/force+3_3th3qv/TracerTimH-',num2str(i));
	%xdat(i).x=time;
	%xdat(i).x=datenum(2004,2,24,hrstartles+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0);
	
	scrsz=get(0,'ScreenSize');
	posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
	
	[izmin izmax]=findheight(z,minZ,maxZ);
  
  
%     if length(iz)>=1
% 		iz=iz(1);
%     else
%         iz=length(z);
%     end
	
	
	%pcolor(9+time./3600,z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp
	
    h(i).h=subplot(a,b,i);
    
    switch plotcase
    case 62
        lem_min_temp
    case 61
        mpc_min_temp
    
    case 60
        mpc_tot_satmr
        close all
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
        
    case 58
        dq_dehyd_10thSep2005
        
    case 57
        switch i
            
        case 1
        [izmin izmax]=findheight(Grid.Z,14e3,22e3);
        timesTH(1).t=Grid.Y1'/1000;
        zz(1).z=GridDan(1).Z(izmin:izmax);
        timesTH(1).t=GridDan(1).Y1(:)'/1000;
        pdat(1).p=fact*squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[1:6]),3));
        imaxovr=[1 0];
        maxcovOvr=8;
        iminovr=[1 0];
        mincovOvr=1;
        
        case 2
        [izmin izmax]=findheight(Grid.Z,14e3,22e3);
        zz(2).z=GridDan(1).Z(izmin:izmax);
        timesTH(2).t=GridDan(1).Y1(:)'/1000;
        pdat(2).p=1e-6*squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,7:9),3));
        %maxcovOvr=100;
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
       pdat(i).p=f*sum(icediag4(idir).i(izmin:izmax,dumprange,35:36),3);
    case 49
        
        T=TwoD.TH2./(1e5./TwoD.PP).^0.286;
        ei=SatVapPress(T,'goff','ice'); %Pa
        P=GridDan(2).PREFN; %Pa
    
    xdat(6).x=f*0.622*ei./(P-ei);
    
        case 48
        pdat(1).p=abs( squeeze( MaxW(idir).w(izmin:izmax,dumprange) ) ); 
        
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
        
        case 44
        switch i
        case 1
        [x1 x2]=findheight(x,200e3,500e3);
        zz(1).z=z(izmin:izmax);
        %tt=findheight(time,time2d)
        pdat(1).p=squeeze(mean(icemr(1).i(izmin:izmax,(x1:x2),:)+snowmr(1).i(izmin:izmax,(x1:x2),:)+graupelmr(1).i(izmin:izmax,(x1:x2),:),2)); 
        
        case 2
        [x1 x2]=findheight(xmpc,200e3,500e3);
       % timesTH(2).t=xmpc(x1:x2)/1000;
        [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
        zz(2).z=zmpc(izminmpc:izmaxmpc);
        %tt=findheight(TimeMPC,time2d);
        pdat(2).p=squeeze(mean(icemrmpc((x1:x2),1,izminmpc:izmaxmpc,:),1))*1e3;
        end
        
        case 43
        switch i
        case 1
        [x1 x2]=findheight(x,200e3,500e3);
        zz(1).z=z(izmin:izmax);
        %tt=findheight(time,time2d)
        pdat(1).p=squeeze(max(icenc(1).i(izmin:izmax,(x1:x2),:)+snownc(1).i(izmin:izmax,(x1:x2),:)+graupelnc(1).i(izmin:izmax,(x1:x2),:),[],2))*1e-6; 
        
        case 2
        [x1 x2]=findheight(xmpc,200e3,500e3);
       % timesTH(2).t=xmpc(x1:x2)/1000;
        [izminmpc izmaxmpc]=findheight(zmpc,minZ,maxZ);
        zz(2).z=zmpc(izminmpc:izmaxmpc);
        %tt=findheight(TimeMPC,time2d);
        pdat(2).p=squeeze(max(icempc((x1:x2),1,izminmpc:izmaxmpc,:),[],1));
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
        
        case 41
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
        ptot= sum(icediag(1).i(izmin:izmax,dumprange,[2 16 3 4 17:23]),3);

        
        mpc = sum(icediag(1).i(izmin:izmax,dumprange,[2 16 19 22 23]),3);
        
        %brat(brat==0)=1e-10; %ensures that when b=0 a value shows up in the ratio
        
        pdat(i).p=mpc./ptot;
        
    case 36
         %ratio of process rates included in MPC model to those not included for snow NC calc in LEM
        ptot= sum(icediag(1).i(izmin:izmax,dumprange,[25 13 24 26 27 28]),3) ...
                + sum(icediag(1).i(izmin:izmax,dumprange,[19 23 22]),3).*ndivqavSnow(izmin:izmax,dumprange);

        
        mpc = sum(icediag(1).i(izmin:izmax,dumprange,[13 24]),3) ...
              + sum(icediag(1).i(izmin:izmax,dumprange,[19 23 22]),3).*ndivqavSnow(izmin:izmax,dumprange);
        
        %brat(brat==0)=1e-10; %ensures that when b=0 a value shows up in the ratio
        
        pdat(i).p=mpc./ptot;
        
     case 35
         %ratio of process rates included in MPC model to those not included for ice NC calc in LEM
        arat=( icediag(1).i(izmin:izmax,dumprange,8)/MI0 + icediag(1).i(izmin:izmax,dumprange,13) );

        
        brat= ( sum(icediag(1).i(izmin:izmax,dumprange,[7 9 12]),3)/MI0 + sum(icediag(1).i(izmin:izmax,dumprange,[3:6]),3).*ndivqav(izmin:izmax,dumprange)+ icediag(1).i(izmin:izmax,dumprange,14) );
        
        %brat(brat==0)=1e-10; %ensures that when b=0 a value shows up in the ratio
        
        pdat(i).p=arat./(arat+brat);
        
      
              
     case 34
         %ratio of process rates included in MPC model to those not included for ice MR calc in LEM
        mpc=( sum(icediag(1).i(izmin:izmax,dumprange,[10 8]),3) + sum(icediag(1).i(izmin:izmax,dumprange,[2 15 1]),3) );
        non=( sum(icediag(1).i(izmin:izmax,dumprange,[7 9 11]),3) + sum(icediag(1).i(izmin:izmax,dumprange,[5 3 12 6 4]),3) ) ;
        
        pdat(i).p=mpc./(mpc+non);
        
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
    
        
    if logflag==1
        pdat(i).p=log10(pdat(i).p);
        pdat(i).p(pdat(i).p==-Inf)=NaN;
    end
    
    maxC(i)=max(max(pdat(i).p));
    minC(i)=min(min(pdat(i).p));
    
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
    xlabelstr='UTC Time (hrs)';
end




for i=1:nplots2d
    
    if notsame==1 
        maxCov=maxC(i);
	    minCov=minC(i);
    end
    
    if iminovr(i)==1
        minCov=mincovOvr;
    end

    if imaxovr(i)==1
        maxCov=maxcovOvr;
    end

    minc=minCov-0.0*abs(minCov);
    maxc=maxCov+abs(maxCov*0.0);
    
%     fixmin=fix2(minc,abs(round(log10(min(abs(minc))))));
%     dd=(maxc-minc)/ncont;
%     dfix=fix2(dd,ceil(abs(log10(min(abs(dd))))));
%     conts=[fixmin:dfix:maxc];
    %conts=[minc:(maxc-minc)/ncont:maxc];
    %conts=round2(conts,abs(round(log10(min(abs(conts)))))+1);
    
    if logflag==0
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
        if sign(conts(1))~=sign(conts(end))
            [zeromin izero]=min(abs(conts));
            if abs(conts(izero))<dfix/2.2
                conts(izero)=0;
            else
                if conts(izero)<0
                    conts(izero+n2:end+1)=conts(izero+1:end);
                    conts(izero+1)=0;
                else
                    conts(izero+1:end+1)=conts(izero:end);
                    conts(izero)=0;
                end
            end
        end
                
    else
        conts=[minc:(maxc-minc)/ncont:maxc];
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
    
    pdat(i).p(1:half(i),npend(i)-nend(i):npend(i))=minCov;
    pdat(i).p(half(i)+1:end,npend(i)-nend(i):npend(i))=maxCov;
    
    h(i).h=subplot(a,b,i);
    
    if izovr==0;
        zz(i).z=z(izmin:izmax);
    end
    
	%pcolor(timesTH,0.62+z(izmin:izmax)./1000,pdat);
    
    [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,0.62+zz(i).z./1000,pdat(i).p,conts);
    if clab==1
          
        
        
       if logflag==1
            ch=cbfA(i).c;
            jc=1;
            while jc<size(cbfA(i).c,2)
                %ch(1,jc)=str2num(num2str(10^(cbfA(i).c(1,jc)),'%2.2e'));
                
                    
                    if 10^cbfA(i).c(1,jc)>99
                        ch(1,jc)=sigfig(10^(cbfA(i).c(1,jc) ),3);
                    else
                        ch(1,jc)=sigfig(10^(cbfA(i).c(1,jc) ),2);
                    end
           
              
                
                jc=jc+cbfA(i).c(2,jc)+1;
            end
        else %logflag==0
            ch=cbfA(i).c;
        end
        if length(ch)>0
            if manclab==0
                clabel(ch,cbfB(i).c);
            else
                clabel(ch,cbfB(i).c,'manual');
            end
        end
    end
    %hc=colorbarf(cbfA,cbfB);
    
    if clines==0
        shading flat; %shading stops black contour lines from appearing
    end
	
    %xti=set(h(i).h,'xticklabels',timestxt);
    
    set(h(i).h,'fontsize',fsize);
    xlabel(xlabelstr);
	ylabel('Height (km)');
    title(tit(i).tit);
    
    

    if length(cbfA(i).c)==0
        nocbar=1;
        normcbar=1;
    end
    %caxis(h(i).h,[minCov maxCov*1.05]);
    if normcbar==1
        hc=colorbar( 'peer' , h(i).h );
	else
        if nocbar~=1
            hc=colorbarf(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!  
        end    
    end
    
    if normcbar==0
    
      if logflag==1 & nocbar==0
        
        clear ctickstr
        
        ctick=get(hc,'yticklabel');
        ctickstr(1,1)=' ';
        for j=2:length(ctick)-1
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            nu=str2num(ctick{j});
            
                %te=num2str(10^nu,'%2.2f');       %'%2.2e');
                if 10^nu>99
                    te=num2str(sigfig(10^nu - offset ,3));
                else
                    te=num2str(sigfig(10^nu - offset,2));
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

      else
%         clear ctickstr
%         ctick=get(hc,'ytick');
%         for j=1:length(ctick)
%             %te=strcat('10^','{',num2str(ctick(j)),'}' );
%             te=num2str(ctick(j),'%2.2e');
%             ctickstr(j,1:length(te))=te;
%         end
%         
%         set(hc,'yticklabel','');
%         
%         add=ctick(end)/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
%     
%     
%         
%         set(hc,'yticklabel','');
%         
%         add=ctick(end)/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
%         
%         set(hc,'fontsize',fsize-6);
        
      end
    
    else
    
    if nocbar==0
        set(hc,'fontsize',fsize-6);
    end
    
    if logflag==1 & nocbar==0
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
        text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
    
    
        
        set(hc,'yticklabel','');
        
        add=ctick(end)/50;
        
        set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
    
    elseif nocbar==0
        clear ctickstr
        ctick=get(hc,'ytick');
        for j=1:length(ctick)
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
    
    
    
    set(hc,'fontsize',fsize-6);
    
    end
    
ylims=get(h(i).h,'ylim');   %re-scale to hide extra column put in to get the colorbars the same in both plots
if iylim==1
    ylims=iylims;
end
axis(h(i).h,[timesTH(i).t(1) timesTH(i).t(pend(i)) ylims]);
    
end

if itimestamp==1
    text(-300,34.5,['Time = ' timlab ' UTC'],'fontsize',18);
end


% 
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');

if (i2d~=1 & i2d~=2 & i2d~=3)
    xx=get(gca,'xticklabels');
    xx=str2num(xx);
    xx=num2str(mod(xx,24));
    set(gca,'xticklabels',xx);
end

set(gcf,'paperpositionmode','auto');


if isave==1
     set(gcf,'paperpositionmode','auto');
	print(gcf,'-djpeg','-r350',exname);
    %print(gcf,'-dmeta',exname);
	%close(gcf);
end

