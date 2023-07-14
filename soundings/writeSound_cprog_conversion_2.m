function writeSound_cprog_conversion_2(comp,HOUR,local_diff,XDMP,HDMP,HGD,KGD,...
    DXX,DYY,RLAT,RLON,heights,pots,press,qvap,vels,velsU,fluxfactorS,fluxfactorL,DAY,fluxTIM,i3d)
%converted c program for writing soundings converted to matlab 16/10/06
%function writeSound_cprog_conversion_2(HOUR,heights,pots,press,rh,vels,velsU)
%HOUR=time of day in UTC for sens and latent heat fluxes
%H 31225.; %top height in sounding - heights are as they appear in sounding (not as in namelist)
H=54565.; % 13/02 at 12UTC
%HOUR=20.25; %15.0 for 13/02  19.75=24thFeb(used most recently 03/03/06) %20.25=24thFeb  %UTC sounding time (float so need decimal point)
%DAY=55; %44=13th Feb 55=24thFeb%day no.(integer)
NSPINUP=0; %no hrs to delay changes in flux by in order to allow for spin up
%automatically adjusts for summertime depending on day no.

ETA=0; %flag to read eta files
every=0; %flag to read in every point in sounding
Uvel=0; % flag to decide between u and v components of wind for 2D LES 0=v 1=u

QnotRH=1; %flag to use q fields directly (=1) instead of converting RH to q (!=1) (for +n deg runs)
LES=0; %flag to read in LES style sounding with winds already as components

if ~exist('i3d'); i3d=1; end
%i3d=1; %flag for 3d namelist

%DXX=1000.0; %500000.0
%DYY=1000.0; %101.8016 %resolution (float)

%H1=0  %changing H1 will give you the height of second point
%H1=335.4



NOP=260; %49 %no. of points
NOP=49;
NPmax=260; %max no points for arrays
NPmax=49; %max no points for arrays

PI=3.14159;

%KGD1+KGD2=NOP-2=16

SHFLX_SEN=500.; %sensible surface heat flux (W/m^2)
SHFLX_LAT=800.; %sens surf latent heat flux (W/m^2)
ZTOP=0; %domain top - put as zero if want it as top of sounding data
RHADD=0; %amount to add to RH to get super sat by cheating - set to zero normally
RMLMAX=23; %taken out so LES calculates automatically

%HDMP=2000.; %damping height wave thing - Be sure to include decimal or won't output properly
%XDMP=3000.; %height below domain top for damping layer  orig=2500

DTM=1.0; %time step
DTMMAX=50.0; %max time step
TINIT=0.; %temp perturbation of bubble - need to set this to zero to swtich off bubble
ZTINIT=1300.; %height of bubble
WTINIT=1500.; %width of bubble
NTINIT=1.0; %number of bubbles
Z0=0.14; %surface roughness for momentum
Z0TH=0.014; %and for scalars




%fluxfactorS=0.813; %factor to multiply sens flux by
%fluxfactorL=1.2672; %factor to multiply lat flux by

finebound=0; %flag to give lowest 2km a high resoultion profile

NTMHF=36; %no. times for surface flux variation
midSEN=450.; %sensible heat flux at time defined by fluxTIM sens=328.6 lat=347.232
midLAT=600.; %latent heat flux at time defined by fluxTIM (W/m2) 300 eta-bauru 450

%midSEN 450.; %sensible heat flux at 4pm 155 eta-bauru 225
%midLAT 600.; %latent heat flux at 4pm (W/m2) 300 eta-bauru 450

%FLUXTIM=1;  %flag for time varying surface fluxes
%fluxTIM=4;  %(4=multiply by factors below)
%fluxTIM=0; %0=keep as in input curve 2=max value required(at 1pm)
noCurvePs=24; %no points in the flux curves



WINDf=1; %reduce wind by a factor 1=no reduction
WIND=0; %set to one for atificial wind, otherwise 0 for sounding wind
windlapse=1; %rate of wind increase per km
dircomp=37.5; %135.0 %135=SE for 24th Feb degrees clockwise for which to take component from - put as -1 for wind as before
                         %13th =37.5

j=0;
hnext=0;

switch comp
case 'laptop'
	output=fopen('c:/documents and settings/g/My Documents/Nml.txt','w')
case 'pc'
	output=fopen('c:/documents and settings/Login/My Documents/cprogs/Nml.txt','w')
end    


times=[0:1:24]; %times of day for flux curve given as UTC-3


senflux(1)=-54.040; %sensible heat flux for flux curve
senflux(2)=-51.640;
senflux(3)=-49.920;
senflux(4)=-48.072;
senflux(5)=-45.256;
senflux(6)=-42.792;
senflux(7)=-40.224;
senflux(8)=-25.840;
senflux(9)=52.504;
senflux(10)=137.52;
senflux(11)=218.192;
senflux(12)=274.296;
senflux(13)=315.720;
senflux(14)=328.600;
senflux(15)=318.472;
senflux(16)=281.720;
senflux(17)=219.232;
senflux(18)=140.976;
senflux(19)=48.968;
senflux(20)=-22.320;
senflux(21)=-45.488;
senflux(22)=-51.312;
senflux(23)=-54.080;
senflux(24)=-53.040;
senflux(25)=-54.040;
 
latflux(1)=21.160; %latent heat flux for flux curve
latflux(2)=17.496;
latflux(3)=13.344;
latflux(4)=9.904;
latflux(5)=6.112;
latflux(6)=1.928;
latflux(7)=-0.336;
latflux(8)=13.712;
latflux(9)=80.464;
latflux(10)=151.768;
latflux(11)=229.600;
latflux(12)=290.640;
latflux(13)=330.784;
latflux(14)=347.232;
latflux(15)=342.304;
latflux(16)=317.464;
latflux(17)=279.160;
latflux(18)=230.184;
latflux(19)=160.432;
latflux(20)=78.016;
latflux(21)=43.072;
latflux(22)=40.632;
latflux(23)=38.376;
latflux(24)=34.304;
latflux(25)=21.160;

%%%%%%%%%% define times to be interpolated (timesIN)
local=HOUR+local_diff; %time curve set up for local time

for (i=1:NTMHF)
    timesIN(i)=local+i-1; %input times of day for fluxes in namelist for fluxes every hour
end

timesOUT(1)=0.;
for (i=2:NTMHF)
    timesOUT(i)=(timesIN(i)-timesIN(1))*3600;   %timesOUT is converted to seconds from the first time as LEM namelist requires
    timesIN(i)=timesIN(i)-floor(timesIN(i)/24.)*24; %convert ones that go over 24:00 to wrap around
end
    
    
    
    for (i=1:noCurvePs) 
    
        if (strcmp(fluxTIM,'timvar_fluxes')==1) 
                senflux(i)=senflux(i)*fluxfactorS;
                latflux(i)=latflux(i)*fluxfactorL;
        end
                
        if (fluxTIM==1) 
                senflux(i)=senflux(i)*midSEN/161;
                latflux(i)=latflux(i)*midLAT/250;
        end
        if (fluxTIM==2) 
                senflux(i)=senflux(i)*midSEN/328.6;
                latflux(i)=latflux(i)*midLAT/347.232;
        end
        if (fluxTIM==3) 
                senflux(i)=senflux(i)*midSEN/179;
                latflux(i)=senflux(i)*midLAT/420;
        end
        %printf("sen(%d)=%f lat(%d)=%f\n",i,senflux(i),i,latflux(i));
    end


    
%spline(times,senflux,noCurvePs,1E30,1E30,y);
%spline(times,latflux,noCurvePs,1E30,1E30,y2);
%splint(times,senflux,y,noCurvePs,timesIN(0),&ansy);
%splint(times,latflux,y2,noCurvePs,timesIN(0),&ansy2);

senOUT(1)=interp1(times,senflux,timesIN(1));
latOUT(1)=interp1(times,latflux,timesIN(1));

for (i=2:NSPINUP+1) %if want to delay changes in fluxes
    senOUT(i)=senOUT(1);
    latOUT(i)=latOUT(1);
    %printf("\nsenOUT(%d)=%lf latOUT(%d)=%lf timesIN=%lf",i,senOUT(i),i,latOUT(i),timesIN(i));
    %printf("CHECK");
end

for (i=NSPINUP+2:NTMHF)
%     splint(times,senflux,y,noCurvePs,timesIN(i-NSPINUP),&ansy);
%     splint(times,latflux,y2,noCurvePs,timesIN(i-NSPINUP),&ansy2);
%     senOUT(i)=ansy;
%     latOUT(i)=ansy2;
    
    senOUT(i)=interp1(times,senflux,timesIN(i-NSPINUP));
	latOUT(i)=interp1(times,latflux,timesIN(i-NSPINUP));
    
    %printf('\time=%f nsen=%f\n',timesIN(i-NSPINUP),senOUT(i));
    %printf("\nsenOUT(%d)=%lf latOUT(%d)=%lf timesIN=%lf",i,senOUT(i),i,latOUT(i),timesIN(i));
    %printf("CHECK");
end
  
    



%if ((input=fopen("c:/program files/pcgrads/win32e/30.01-00/profs/prof-18","r"))==NULL)
%if ((input=fopen("c:/documents and settings/Login/My Documents/hibiscus/baurufield/soundings/ascii/3/04_02_13_12AED.txt","r"))==NULL)

%if ((input=fopen("c:/documents and settings/Login/My Documents/Leeds_MMOCCA/MileCity.AED","r"))==NULL)

%if ((input=fopen("c:/documents and settings/Login/my documents/HIBISCUS/troccibras/radiosondesnetwork/LES_dump45_SDLA","r"))==NULL)
%if ((input=fopen("c:/documents and settings/Login/my documents/HIBISCUS/troccibras/radiosondesnetwork/CampoGrande_LES_dump45","r"))==NULL)
% if ((input=fopen("c:/documents and settings/Login/my documents/HIBISCUS/troccibras/radiosondesnetwork/DMI_1715_5ppmv_16.6-19km","r"))==NULL)
%if ((input=fopen("c:/documents and settings/Login/my documents/HIBISCUS/troccibras/radiosondesnetwork/sdla_24thFeb_5ppmv_2nd","r"))==NULL)

%input=fopen('c:/documents and settings/Login/my documents/HIBISCUS/troccibras/radiosondesnetwork/lem_ecmwf_hightop','r');


fprintf(output,'&CNTRL NN=4,NNDIAG=4,NNDUMP=3,ISTART=1 / &END\n');
fprintf(output,'&TIMENML NTMHALT=1,TIMHALT(1)=1000.,IPRTDG=1,NSTEPMAX=10000 / &END\n');
fprintf(output,'&JOBINFO FILEA=''filea'',FILEB=''fileb'',FILEZ=''filez'',MSGFILE=''mesg'' / &END\n');

press0=press(1);
fprintf(output,'&INPUT NRUN=1,Z0=%.4f,Z0TH=%.4f,PSF=%.1f,',Z0,Z0TH,press0);

if strcmp(fluxTIM,'timvar_fluxes')==1 
    fprintf(output,'NTMHF=%d',NTMHF);
	for (i=1:NTMHF) 
        fprintf(output,',\nTIMHF(%d)=%.1f',i,timesOUT(i));
	end
	for (i=1:NTMHF) 
        fprintf(output,',\nFSHFLX_SEN(%d)=%.1f',i,senOUT(i));
	end
	for (i=1:NTMHF) 
        fprintf(output,',\nFSHFLX_LAT(1,%d)=%.1f',i,latOUT(i));
	end    
	
	fprintf(output,' / &END\n');
	
elseif strcmp(fluxTIM,'const_fluxes')==1 
    fprintf(output,'SHFLX_SEN=%.1f,\nSHFLX_LAT(1)=%.1f  / &END\n',SHFLX_SEN,SHFLX_LAT);
    
elseif strcmp(fluxTIM,'const_surface')==1     
	fprintf(output,'THSURF=%.1f,\nQSURF=%.1f  / &END\n',fluxfactorS,fluxfactorL);
        
elseif strcmp(fluxTIM,'const_surface_sat')==1     
	fprintf(output,'THSURF=%.1f  / &END\n',fluxfactorS); %no need for moisture, as saturated
    
end

%%%%%%%%%%% GRID
zztop=heights(end);
fprintf(output,'&GRID NSMTH=5,\n');
for i=1:length(KGD)-1
	fprintf(output,'HGD(%d)=%f,KGD(%d)=%d\n,',i,HGD(i),i,KGD(i));
end
fprintf(output,'HGD(%d)=%.1f,KGD(%d)=%d\n,ZZTOP=%.1f, DYY=%.4f, ',i+1,zztop,i+1,KGD(i+1),zztop,DYY);

if (i3d==1) 
    fprintf(output,'DXX=%.4f ',DXX);
end


%%%%%%%%%%%% THPROF
fprintf(output,'/ &END\n&THPROF\n');

nps=length(heights);
for (j=1:nps)
    fprintf(output,'ZNREF_READ(%d) = %.1f, THREF_READ(%d) = %.1f,\n',j,heights(j),j,pots(j));
end

fprintf(output,'THREF0=%.1f / &END\n&INITPROF\n',pots(1));
if i3d==1
	fprintf(output,'L_INITU=.TRUE.\nL_INITV=.TRUE.\nL_INITTH=.TRUE.\nL_INITQ(1)=.TRUE.\n');
else
	fprintf(output,'L_INITU=.FALSE.\nL_INITV=.TRUE.\nL_INITTH=.TRUE.\nL_INITQ(1)=.TRUE.\n');
end
    

for (j=1:nps)
    fprintf(output,'ZNTHI_READ(%d) = %.1f, THINIT_READ(%d) = %.1f,\n',j,heights(j),j,pots(j));
end

    
%%%%%%%%%% U winds
if (i3d==1)
    for (j=1:nps)
        fprintf(output,'ZNUI_READ(%d) = %.1f, UINIT_READ(%d) = %.2f,\n',j,heights(j),j,velsU(j));
    end
end

%%%%%%%%%% V winds
for (j=1:nps)
    fprintf(output,'ZNVI_READ(%d) = %.1f, VINIT_READ(%d) = %.2f,\n',j,heights(j),j,vels(j));
end

%%%%%%%%%% vapour
for (j=1:nps)    
    fprintf(output,'ZNQI_READ(%d,1) = %.1f, QINIT_READ(%d,1) = %.3e,\n',j,heights(j),j,qvap(j));
end

fprintf(output,'NT_INIT=%.1f,T_INIT=%.1f,W_INIT=%.1f,ZT_INIT=%.1f\n/ &END\n&SUBMODEL / &END\n',NTINIT,TINIT,WTINIT,ZTINIT);
if (i3d==1) 
    fprintf(output,'&DIAGNOST NITEST=1,XTEST(1)=0.,IDGP=3,IDGU=3,IDGV=3,IDGW=3,IDGTH=3,\nIDGQMR_ALL=3,IDGQNC_ALL=3,IDGQ(1)=3,IDGTHETA=3,IDGTHV=3,\nIDGAV=1,IDGPV=2,');
else
    fprintf(output,'&DIAGNOST NITEST=1,XTEST(1)=0.,IDGP=2,IDGU=2,IDGV=2,IDGW=2,IDGTH=2,\nIDGQMR_ALL=2,IDGQNC_ALL=2,IDGQ(1)=2,IDGTHETA=2,IDGTHV=2,\nIDGAV=1,IDGPV=2,');
end

fprintf(output,'IDGPD(1)=2,IEXDG=1,');
fprintf(output,'nupdown=27,\ncupdown(1)=''ALL'',cupdown(2)=''ALu'',cupdown(3)=''ALd'',cupdown(4)=''ACC'',');
fprintf(output,'cupdown(5)=''W>1'',\ncupdown(6)=''W<1'',cupdown(7)=''CLu'',cupdown(8)=''CLd'',cupdown(9)=''ACu'',');
fprintf(output,'cupdown(10)=''ACd'',\ncupdown(11)=''BYu'',cupdown(12)=''PPd'',cupdown(13)=''VPd'',cupdown(14)=''PVd'',');
fprintf(output,'cupdown(15)=''AAd'',\ncupdown(16)=''AVd'',cupdown(17)=''AHM'',cupdown(18)=''MOu'',cupdown(19)=''BMu'',');
fprintf(output,'cupdown(20)=''M_1'',\ncupdown(21)=''M1u'',cupdown(22)=''M1d'',cupdown(23)=''M_2'',cupdown(24)=''M2u'',');
fprintf(output,'cupdown(25)=''M2d'',cupdown(26)=''M_3'',cupdown(27)=''BCu'' / &END\n');
fprintf(output,'&DYNAMICS UG0=0.,VG0=0.,DUGDZ=0.0,DVGDZ=0.0,FCORIOL=0.0 / &END\n');
fprintf(output,'&NUMERICS DTM=%.2f,DTMMAX=%.2f,DTMMIN=0.01 / &END\n',DTM,DTMMAX);
fprintf(output,'&DAMPNML DMPTIM=0.001,ZDMP=%.1f,HDMP=%.1f / &END\n&OVRIDE1 / &END\n',XDMP,HDMP);
fprintf(output,'&RADCNTL RLAT=%f, RLONG=%f, HOUR=%.2f,IYEAR=2004, IDAY=%d,\nRAD_INT_TIME=300., VAR_ALB=.true. / &END\n&MOVDATA / &END ',RLAT,RLON,HOUR,DAY);

fclose(output);

fprintf(1,'\nFinished');
