 H=22755.;   %top height in sounding - heights are as they appear in sounding (not as in namelist)
 HOUR=12.0; %20.25=24thFeb  %UTC sounding time (float so need decimal point)
 DAY=55;  %55=24thFeb%day no.(integer)
 NSPINUP=0; %no hrs to delay changes in flux by in order to allow for spin up
%automatically adjusts for summertime depending on day no.

 ETA=0; %flag to read eta files
 every=1; %flag to read in every point in sounding
 Uvel=0; % flag to decide between u and v components of wind for 2D LES 0=v 1=u

 i3d=0; %flag for 3d namelist

 DXX=500000.0;
 DYY=1000.0; %101.8016 %resolution (float)

% H1=0;  %changing H1 will give you the height of second point
% H1=335.4;



 NOP=19; %50 %no. of points
 NPmax=60; %max no points for arrays
 PI=3.14159;

%KGD1+KGD2=NOP-2=16;

 SHFLX_SEN=500.;  %sensible surface heat flux (W/m^2)
 SHFLX_LAT=800.;   %sens surf latent heat flux (W/m^2)
 ZTOP=0;  %domain top - put as zero if want it as top of sounding data
 RHADD=0; %amount to add to RH to get super sat by cheating - set to zero normally
 RMLMAX=23; %taken out so LES calculates automatically

 HDMP=2000.; %damping height wave thing - Be sure to include decimal or won't output properly
 XDMP=2500.; %height below domain top for damping layer

 DTM=1.0; %time step
 DTMMAX=50.0; %max time step
 TINIT=0.; %temp perturbation of bubble - need to set this to zero to swtich off bubble
 ZTINIT=1300.; %height of bubble
 WTINIT=1500.; %width of bubble
 NTINIT=1.0; %number of bubbles
 Z0=0.14; %surface roughness for momentum
 Z0TH=0.014; %and for scalars

 FLUXTIM=4; %flag for time varying surface fluxes (4=multiply by factors below)
 fluxfactorS=0.813; %factor to multiply sens flux by
 fluxfactorL=1.2672; %factor to multiply lat flux by

 finebound=0; %flag to give lowest 2km a high resoultion profile

 NTMHF=36; %no. times for surface flux variation
 midSEN=450.; %sensible heat flux at time defined by fluxTIM sens=328.6 lat=347.232
 midLAT=600.; %latent heat flux at time defined by fluxTIM (W/m2) 300 eta-bauru 450

% midSEN=450.; %sensible heat flux at 4pm 155 eta-bauru 225
% midLAT=600.; %latent heat flux at 4pm (W/m2) 300 eta-bauru 450

 fluxTIM=0; %0=keep as in input curve 2=max value required(at 1pm)
 noCurvePs=24; %no points in the flux curves



 WINDf=1; %reduce wind by a factor 1=no reduction
 WIND=0; %set to one for atificial wind, otherwise 0 for sounding wind
 windlapse=1; %rate of wind increase per km
 dircomp=135.0; %degrees for which to take component from - put as -1 for wind as before




input2=fopen("heights","r");
output=fopen("Nml.txt","w");




if (ETA==1)
    nps=38;
else 
    nps=NOP;
end



%if (DAY<46) {  %adjustment for clocks going back on 15th Feb, 2004.
%    local=HOUR-2;
%    }
%else {
%    local=HOUR-3;
%    }

local=HOUR-3; %time curve set up as UTC-3 HOUR=UTC time

for (i=0;i<=NTMHF-1;i++) 
    timesIN(i)=local+i; %input times of day for fluxes in namelist for fluxes every hour
end


    


%times=(0 2 4 6 10 12.3 13.5 16 18.5 20.5 22 24);
%sens=(-20 -23.3 -26.7 -30 161 250 250 179 -10 -10 -13.3 -16.7);
%lat=(0 0 0 5 250 450 500 420 180 0 0 0);

times(1)=0.; %times of day for flux curve given as UTC-3
times(2)=1.;
times(3)=2.;
times(4)=3.;
times(5)=4.;
times(6)=5.;
times(7)=6.;
times(8)=7.;
times(9)=8.;
times(10)=9.;
times(11)=10.;
times(12)=11.;
times(13)=12.;
times(14)=13.;
times(15)=14.;
times(16)=15.;
times(17)=16.;
times(18)=17.;
times(19)=18.;
times(20)=19.;
times(21)=20.;
times(22)=21.;
times(23)=22.;
times(24)=23.;

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


%senflux(0)=238.3913;
%senflux(1)=250.;
%senflux(2)=235.8;
%senflux(3)=207.4;
%senflux(4)=179.;
%senflux(5)=103.4;
%senflux(6)=27.8;
%senflux(7)=-10.;
%
%latflux(0)=423.913;
%latflux(1)=479.1667;
%latflux(2)=484.0;
%latflux(3)=452.;
%latflux(4)=420.;
%latflux(5)=324.;
%latflux(6)=228.;
%latflux(7)=135.;

timesOUT(0)=0.;
for (i=1:NTMHF-1)
    timesOUT(i)=(timesIN(i)-timesIN(0))*3600;
    timesIN(i)=timesIN(i)-floor(timesIN(i)/24.)*24;
end
    
    
    
    for (i=1:noCurvePs)
    
        if (fluxTIM==4) 
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
                latflux(i)=latflux(i)*midLAT/420;
            end
        %printf("sen(%d)=%f lat(%d)=%f\n",i,senflux(i),i,latflux(i));
    end


    
spline(times,senflux,noCurvePs,1E30,1E30,y);
spline(times,latflux,noCurvePs,1E30,1E30,y2);


%for (i=0;i<=NTMHF-1;i++) {
%printf("y(%d)=%f y2(%d)=%f\n",i,y(i),i,y2(i));
%}

splint(times,senflux,y,noCurvePs,timesIN(0),&ansy);
splint(times,latflux,y2,noCurvePs,timesIN(0),&ansy2);
senOUT(0)=ansy;
latOUT(0)=ansy2;

splint(times,senflux,y,noCurvePs,15.0,&ansy);
printf("\nsen=%f\n",ansy);

for (i=1;i<=NSPINUP;i++) { %if want to delay changes in fluxes
    senOUT(i)=senOUT(0);
    latOUT(i)=latOUT(0);
    %printf("\nsenOUT(%d)=%lf latOUT(%d)=%lf timesIN=%lf",i,senOUT(i),i,latOUT(i),timesIN(i));
    %printf("CHECK");
    }

for (i=NSPINUP+1;i<=NTMHF-1;i++) {
    splint(times,senflux,y,noCurvePs,timesIN(i-NSPINUP),&ansy);
    splint(times,latflux,y2,noCurvePs,timesIN(i-NSPINUP),&ansy2);
    senOUT(i)=ansy;
    latOUT(i)=ansy2;
    
    printf("\time=%f nsen=%f\n",timesIN(i-NSPINUP),ansy);
    %printf("\nsenOUT(%d)=%lf latOUT(%d)=%lf timesIN=%lf",i,senOUT(i),i,latOUT(i),timesIN(i));
    %printf("CHECK");
    }

  
    



%if ((input=fopen("c:/program files/pcgrads/win32e/30.01-00/profs/prof-18","r"))==NULL)
%if ((input=fopen("c:/documents and settings/Login/My Documents/hibiscus/baurufield/soundings/ascii/5/04_02_24_2015_dmi.AED.txt","r"))==NULL)

if ((input=fopen("c:/documents and settings/Login/my documents/HIBISCUS/troccibras/radiosondesnetwork/CampoGrande_040224_12","r"))==NULL)
	{printf("Error reading input file1! \n");
	printf("Press any key to continue\n");
	system("PAUSE");
    exit(0);}

%for (i=1;i<=55;i++) fscanf(input,"%s",&s);
%printf("\n%s\n",s);


%fscanf(input,"%f %f %f %f %f %f\n",&h0,&press,&temp,&rh,&dir,&speed);
if (ETA==1) {
    fscanf(input,"%f %f\n",&h0,&press0);    %"%f %f\n",&press,&h0);
    do {
        fscanf(input,"%f %f %f %f %f %f\n",&height,&press,&temp,&rh,&uv,&vv); 
        }
    while (height<h0);
    if (Uvel==1) {speed=uv;}
    else {speed=vv;}
}
else {
    fscanf(input,"%f %f %f %f %f %f %f %f %f\n",&skip,&skip,&h0,&press,&temp,&rh,&skip,&dir,&speed);
    dir*=PI/180;  %convert direction to radians
    if (dircomp>-0.5) {
       speedU=speed*-sin(dir-dircomp*PI/180);
       speed*=-cos(dir-dircomp*PI/180); }%take component of speed at dircomp degrees
    else {
       speedU=speed*sin(dir);
       speed*=cos(dir);  } %take component of speed at 0 degrees
    press0=press;
}

speed*=WINDf; %multiply wind speed by a factor
speedU*=WINDf;
rh/=100; %convert percentage to fraction - CHECK if RH right for LEM!


H1=h0+200; %define second point at 200m above surface

/*
fscanf(input2,"%d",&totnlev);
for (i=1;i<=totnlev;i++) {
fscanf(input2,"%f %d\n",&hlev(i),&nlev(i));
printf("%f %d\n",hlev(i),nlev(i));
}
%fclose(input2);
*/

if (finebound==1) {
    hlev(1)=2000.+h0;
    nlev(1)=10;
    hlev(2)=H;
    nlev(2)=nps-2-10; 
    }
else {
    hlev(1)=H;
    nlev(1)=nps-2; 
    }

zztop=ZTOP;
if (zztop<1) zztop=H-h0;

%lambda=log(H/H1)/(nps-2);
potemp=(temp+273.15)*pow(1e3/press,0.286);

pots(0)=potemp;
        Y=373.16/(temp+273.15);
        ew=100*1013.246*pow(     10.0 , (    -7.90298*(Y-1.0)  +  5.02808*log10(Y)
     		    -1.3816E-07*(   pow( 10.0,(11.344*(1.0-(1.0/Y))) ) - 1.0   )
     		    +8.1328E-03*(   pow( 10.0,(-3.49149*(Y-1.0)) ) - 1.0   )    )     ); %*100 to convert to Pascals
hums(0)=rh*0.622*ew/(press*100-rh*ew);
vels(0)=speed;
velsU(0)=speedU;

if (ETA==1) heights(0)=height-h0;
else heights(0)=0.0;




%fprintf(output,"ZNREF_READ(1)=1.,THREF_READ(1)=%.1f,\nZNINIT_READ(1)=1.,THINIT_READ(1)=%.1f,\n",potemp,potemp);
%fprintf(output,"ZVINIT(1)=0.0,VINIT(1)=0.0,\nZHINIT(1)=0.0,HINIT(1)=%.2f,\n",rh);



i=1;
hlev(0)=H1;
hnext=H1;
do {
    %hnext=H1*exp(lambda*j);  %takes heights at exponentially increasing positions
    %if (fabs(hnext-hlev(i))>0.001 && hnext<HGD1-0.001) hnext=H1+(j*(HGD1-H1)/KGD1);
    %else hnext=HGD1+((j-KGD1)*(H-HGD1)/KGD2);
    
    
    
    
    
    %printf(" hnext=%lf j=%d %lf %lf %lf\n",hnext,j,lambda*j,lambda,A);
    if (ETA==1) {fscanf(input,"%f %f %f %f %f %f\n",&height,&press,&temp,&rh,&uv,&vv);
    /*printf("\nh=%f h0=%f",height,h0);*/}
    else { fscanf(input,"%f %f %f %f %f %f %f %f %f\n",&skip,&skip,&height,&press,&temp,&rh,&skip,&dir,&speed);}
    
    
    %printf(" %f ",height,hlev(1));
    %getch();
    if ( height>hnext || fabs(height-hnext)<0.01 || ETA==1 || every==1) { %read in every time for ETA soundings
        
        if (hnext>hlev(i) || (fabs(hlev(i)-hnext)<0.01 & fabs(H-hnext)>1.0) || ETA==1 || every==1) i++;
        
        hnext+=(hlev(i)-hlev(i-1))/nlev(i);
        
        %if (hnext<hlev(i)) {
        %hnext+=(hlev(i)-hlev(i-1))/nlev(i);
        %}
        %else {
        %i++;
        %hnext+=(hlev(i)-hlev(i-1))/nlev(i);
        %}
        
        if (ETA != 1) {
                dir*=PI/180;  %convert direction to radians
                if (dircomp>-0.5) {
                speedU=speed*-sin(dir-dircomp*PI/180); %take component of speed at dircomp degrees
                speed*=-cos(dir-dircomp*PI/180); }%take component of speed at dircomp degrees
                else {
                speedU=speed*sin(dir);
                speed*=cos(dir);
                } %take component of speed at 0 degrees
    
        }
        else {
                if (Uvel==1) {speed=uv;}
                else {speed=vv;}
        }
        
        speed*=WINDf; %multiply wind speed by a factor
        speedU*=WINDf;
        
        if (WIND==1) {  % if using atrificial wind
        speed=(height-h0)*windlapse/1000;
        speedU=(height-h0)*windlapse/1000;
        }
        
        
        
        rh/=100; %convert percentage to fraction 
        rh+=RHADD; %adds some to RH if are cheating & artificially increasing RHADD!
        potemp=(temp+273.15)*pow(1e3/press,0.286);  %calc potential temperature
        
        %fprintf(output,"ZNREF_READ(%d)=%.1f,THREF_READ(%d)=%.1f,\nZNINIT_READ(%d)=%.1f,THINIT_READ(%d)=%.1f,\n",j+2,height-h0,j+2,potemp,j+2,height-h0,j+2,potemp);
        %fprintf(output,"ZVINIT(%d)=%.1f,VINIT(%d)=%.2f,\nZHINIT(%d)=%.1f,HINIT(%d)=%.3f,\n",j+2,height-h0,j+2,speed,j+2,height-h0,j+2,rh);
        j++;
        heights(j)=height-h0;
        
        printf("\n%d %d %f %f %f dH=%f i=%d %f",nps,j,heights(j),h0,height,heights(j)-heights(j-1),i,hnext); 
        
        vels(j)=speed;
        velsU(j)=speedU;
        pots(j)=potemp;
        Y=373.16/(temp+273.15);
        ew=100*1013.246*pow(     10.0 , (    -7.90298*(Y-1.0)  +  5.02808*log10(Y)
     		    -1.3816E-07*(   pow( 10.0,(11.344*(1.0-(1.0/Y))) ) - 1.0   )
     		    +8.1328E-03*(   pow( 10.0,(-3.49149*(Y-1.0)) ) - 1.0   )    )     ); %*100 to convert to Pascals
        hums(j)=rh*0.622*ew/(press*100-rh*ew);
        %hums(j)=rh;    
    }
}
while (j<nps-1 && !feof(input));

printf("finished height hunt!");

nps=j+1;
zztop=heights(j);

if (ETA==1) press=press0;

fprintf(output,"&CNTRL NN=4,NNDIAG=4,NNDUMP=3,ISTART=1 &END\n");
fprintf(output,"&TIMENML NTMHALT=1,TIMHALT(1)=1000.,IPRTDG=1,NSTEPMAX=10000 &END\n");
fprintf(output,"&JOBINFO FILEA='filea',FILEB='fileb',FILEZ='filez',MSGFILE='mesg' &END\n");
fprintf(output,"&INPUT NRUN=1,Z0=%.4f,Z0TH=%.4f,PSF=%.1f,",Z0,Z0TH,press0*100);

if (FLUXTIM>0.99) {
fprintf(output,"NTMHF=%d",NTMHF);
for (i=0;i<=NTMHF-1;i++) {
    fprintf(output,",\nTIMHF(%d)=%.1f",i+1,timesOUT(i));
    }
for (i=0;i<=NTMHF-1;i++) {
    fprintf(output,",\nFSHFLX_SEN(%d)=%.1f",i+1,senOUT(i));
    }
for (i=0;i<=NTMHF-1;i++) {
    fprintf(output,",\nFSHFLX_LAT(1,%d)=%.1f",i+1,latOUT(i));
    }
fprintf(output," &END\n");
}
else {
fprintf(output,"SHFLX_SEN=%.1f,\nSHFLX_LAT(1)=%.1f  &END\n",SHFLX_SEN,SHFLX_LAT);
}


    
fprintf(output,"&GRID NSMTH=5, HGD(1)=1500,KGD(1)=21,HGD(2)=18000,KGD(2)=121,\nHGD(3)=%.1f,KGD(3)=150,ZZTOP=%.1f, DYY=%.4f, ",zztop,zztop,DYY);

if (i3d==1) {
    fprintf(output,"DXX=%.4f ",DXX);
    }
    
fprintf(output,"&END\n&THPROF\n");

for (j=0;j<nps;j++) {
    fprintf(output,"ZNREF_READ(%d) = %.1f, THREF_READ(%d) = %.1f,\n",j+1,heights(j),j+1,pots(j));
    }
    
fprintf(output,"THREF0=%.1f &END\n&INITPROF\nL_INITU=.FALSE.\nL_INITV=.TRUE.\nL_INITTH=.TRUE.\nL_INITQ(1)=.TRUE.\n",pots(0));


for (j=0;j<nps;j++) {
    fprintf(output,"ZNTHI_READ(%d) = %.1f, THINIT_READ(%d) = %.1f,\n",j+1,heights(j),j+1,pots(j));
    }
    
%fprintf(output,"UINIT(1)=0.0,ZUINIT(1)=0.0,\nUINIT(2)=0.0,ZUINIT(2)=%.1f,\n",zztop);

if (i3d==1) {
    for (j=0;j<NOP;j++) {
        fprintf(output,"ZNVI_READ(%d) = %.1f, UINIT_READ(%d) = %.2f,\n",j+1,heights(j),j+1,velsU(j));
    }
}
    
for (j=0;j<nps;j++) {
    fprintf(output,"ZNVI_READ(%d) = %.1f, VINIT_READ(%d) = %.2f,\n",j+1,heights(j),j+1,vels(j));
    }
for (j=0;j<nps;j++) {
    
    fprintf(output,"ZNQI_READ(%d,1) = %.1f, QINIT_READ(%d,1) = %.3e,\n",j+1,heights(j),j+1,hums(j));
    }





fprintf(output,"NT_INIT=%.1f,T_INIT=%.1f,W_INIT=%.1f,ZT_INIT=%.1f\n&END\n&SUBMODEL &END\n",NTINIT,TINIT,WTINIT,ZTINIT);
if (i3d==1) {
    fprintf(output,"&DIAGNOST NITEST=1,XTEST(1)=0.,IDGP=3,IDGU=3,IDGV=3,IDGW=3,IDGTH=3,\nIDGQMR_ALL=3,IDGQNC_ALL=3,IDGQ(1)=3,IDGTHETA=3,IDGTHV=3,\nIDGAV=1,IDGP=1,IDGPV=0,");
}
else {
    fprintf(output,"&DIAGNOST NITEST=1,XTEST(1)=0.,IDGP=2,IDGU=2,IDGV=2,IDGW=2,IDGTH=2,\nIDGQMR_ALL=2,IDGQNC_ALL=2,IDGQ(1)=2,IDGTHETA=2,IDGTHV=2,\nIDGAV=1,IDGP=1,IDGPV=0,");
}
fprintf(output,"IDGPD(1)=2,IEXDG=1,");
fprintf(output,"nupdown=27,\ncupdown(1)='ALL',cupdown(2)='ALu',cupdown(3)='ALd',cupdown(4)='ACC',");
fprintf(output,"cupdown(5)='W>1',\ncupdown(6)='',cupdown(7)='CLu',cupdown(8)='CLd',cupdown(9)='ACu',");
fprintf(output,"cupdown(10)='ACd',\ncupdown(11)='BYu',cupdown(12)='PPd',cupdown(13)='VPd',cupdown(14)='PVd',");
fprintf(output,"cupdown(15)='AAd',\ncupdown(16)='AVd',cupdown(17)='AHM',cupdown(18)='MOu',cupdown(19)='BMu',");
fprintf(output,"cupdown(20)='M_1',\ncupdown(21)='M1u',cupdown(22)='M1d',cupdown(23)='M_2',cupdown(24)='M2u',");
fprintf(output,"cupdown(25)='M2d',cupdown(26)='M_3',cupdown(27)='BCu' &END\n");
fprintf(output,"&DYNAMICS UG0=0.,VG0=0.,DUGDZ=0.0,DVGDZ=0.0,FCORIOL=0.0 &END\n");
fprintf(output,"&NUMERICS DTM=%.2f,DTMMAX=%.2f,DTMMIN=0.01 &END\n",DTM,DTMMAX);
fprintf(output,"&DAMPNML DMPTIM=0.001,ZDMP=%.1f,HDMP=%.1f &END\n&OVRIDE1 &END\n",zztop-XDMP,HDMP);
fprintf(output,"&RADCNTL RLAT=-22.36, RLONG=-49.03, HOUR=%.2f,IYEAR=2004, IDAY=%d,\nRAD_INT_TIME=300., VAR_ALB=.true. &END\n&MOVDATA &END ",HOUR,DAY);

fclose(output);
fclose(input);
printf("Finished, press a key");
getch();    
}