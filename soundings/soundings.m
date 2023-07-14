%prepare pr from ecmwf output combined with LES for writing sounding
%load in first LEM dump first
f=1e6*28.97/18;

clear diff

comp='pc';
%comp='laptop';

clear press temp qvap heights vels velsU HGD KGD

velsU=0;

isounding=1; % flag to say whether want to prepare a souding from the data

sound='ecmwf high top 2';
%sound='Miles City';
sound='Miles City_inversion';
%sound='24Feb';
sound='24Feb_orig'; %using the original sounding as from TROCCIBRAS

% sound='24Feb - hot BL';
% sound='TOGA-COARE';
% sound='newmexico';
sound='texas';
%sound='TOGA-COARE_less_CIN';
%sound='New_Mexico_EMM';
sound='Ant_6thJan06';
sound='mac3_output';
%sound='radar_titan';


DXX=0; %default x res for 2-d cases

switch sound
    
case 'radar_titan'

    isounding=0; %to stop it from writing a sounding  
    
   
        
        exname=['c:\documents and settings\Login\my documents\Gerhard\titan_tracks\titan_27thFeb.txt'];
        
        fid=fopen(exname,'rt');
        %temp_var=fscanf(fid,'%g',[1]);
        
        
        headings=fscanf(fid,'%s',[73 1]);
        dat=fscanf(fid,'%g',[73 inf]);
        %	dat=dlmread(exname,'\t');
        %    M=length(dat)/15;
        %    dat=reshape(dat,15,M);


        
        fclose(fid);
        

    
case 'mac3_output'

    isounding=0; %to stop it from writing a sounding
    dt=5;
    nz=126;  %(dz=120 m)
    nr=101;  %(dx=60 m)    
    
    for it=1:12
        t=it*dt;
        if t<10
            tstr=['0' num2str(t)];
        else
            tstr=[num2str(t)];
        end   
        
        %[b9 b5 b2 b4]= aerosol * [0.1 0.5 2.0 4.0]
        
        exname=['c:\documents and settings\Login\my documents\Leeds_MMOCCA\MAC3_output\wu_' tstr '.dat'];
        
        fid=fopen(exname,'rt');
        %temp_var=fscanf(fid,'%g',[1]);
        
        
        
        dat=fscanf(fid,'%g',[nr inf]);
        %	dat=dlmread(exname,'\t');
        %    M=length(dat)/15;
        %    dat=reshape(dat,15,M);
        
        sidat=(size(dat,2));
        
        
        i=1;
        mac3(it).w=dat(:,(i-1)*nz+1:i*nz); i=i+1;  %updraught
        mac3(it).u=dat(:,(i-1)*nz+1:i*nz); i=i+1;  %horiz wind
        mac3(it).ql=dat(:,(i-1)*nz+1:i*nz); i=i+1; %drop mass mixing ratio (g/kg)
        mac3(it).qi=dat(:,(i-1)*nz+1:i*nz); i=i+1; %ice mass
        mac3(it).qg=dat(:,(i-1)*nz+1:i*nz); i=i+1; %graupel mass
        mac3(it).qs=dat(:,(i-1)*nz+1:i*nz); i=i+1; %snow mass
        mac3(it).nd=dat(:,(i-1)*nz+1:i*nz); i=i+1; %drop number concentration (#/m3)
        mac3(it).ni=dat(:,(i-1)*nz+1:i*nz); i=i+1; %ice number
        mac3(it).ng=dat(:,(i-1)*nz+1:i*nz); i=i+1; %graupel number
        mac3(it).ns=dat(:,(i-1)*nz+1:i*nz); i=i+1; %snow number
        
        fclose(fid);
        
    end


	
%    time=dat(1,:)'/3.6e6; %time is in milliseconds

    
case 'Ant_6thJan06'
    
    isounding=0; %to stop it from writing a sounding
    exname='c:\documents and settings\Login\my documents\logbook\BAS\6th_Jan_flight_data';
    
    fid=fopen(exname,'rt');
    %temp_var=fscanf(fid,'%g',[1]);
    
%	dat=fscanf(fid,'%g',[15 inf]);
	dat=dlmread(exname,'\t');
    M=length(dat)/15;
    dat=reshape(dat,15,M);

	sidat=(size(dat,2));
	
    time=dat(1,:)'/3.6e6; %time is in milliseconds
	lat=dat(2,:)'; %
	lon=dat(3,:)'; %
	press=dat(4,:)'; %
	surf_temp=dat(5,:)'; %from infra-red thermometer              
	air_temp=dat(6,:)'; %
    tdew_hyg=dat(7,:)'; %from frost point hygrometer
    tdew_humicap=dat(8,:)'; %from humicap
    wind_speed=dat(9,:)'; %
    wind_direction=dat(10,:)'+180; %180 added to convert to conventional direction
    alt=dat(11,:)'; %from gps (not very good)
    sw_up=dat(12,:)'; %upwelling short wave
    sw_down=dat(13,:)'; %downwelling short wave
    lw_up=dat(14,:)'; %upwelling long wave
    lw_down=dat(15,:)'; %downwelling long wave
    

     
     
case 'EMM_storm_case2'
    EMM_storm_case2
	
	temp=temp_K;
	press=dat(1,:)'*100; %pressure
	tdew=dat(3,:)'+273.15; %temp
    
    qvap=satvappress(tdew,'goff','liq',press,1)/f;                
    
    Y0=press(1); %pressure at ground for Mile City
	TSPAN=temp; %temperature range
	PSPAN=[press(1) press(end)]; %height range for temperature interpolation above ground
	[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],press,temp); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
   
    heights=interp1(P,h,press);    
	
        
    qsw=satvapPress(temp,'buck2','liq',press,1)/f;        
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=20.25;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=heights(end)-h0; KGD(2)=250;
    
    vels=zeros(length(heights));
    velsU=zeros(length(heights));
 
    
    DYY=2000;
    DXX=2000;
        
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb

    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by

    
     fluxtim='timvar_fluxes';  
     
case 'New_Mexico_EMM'
	clear p
	h0=620;    
    %think this is New Mexcio but it doesn't match the one I got from Stewart
    exname='c:\documents and settings\Login\my documents\emm_0013a\text\STORM_CASE_1_ENV.txt';
    
    fid=fopen(exname,'rt');
    temp_var=fscanf(fid,'%g',[1]);
	dat=fscanf(fid,'%g',[3 inf]);
	sidat=(size(dat,2));
	
	temp=dat(2,:)'+273.15; %temp
	press=dat(1,:)'*100; %pressure
	tdew=dat(3,:)'+273.15; %temp
    
    qvap=satvappress(tdew,'goff','liq',press,1)/f;                
    
    Y0=press(1); %pressure at ground for Mile City
	TSPAN=temp; %temperature range
	PSPAN=[press(1) press(end)]; %height range for temperature interpolation above ground
	[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],press,temp); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
   
    heights=interp1(P,h,press);    
	
        
    qsw=satvapPress(temp,'buck2','liq',press,1)/f;        
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=20.25;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=heights(end)-h0; KGD(2)=250;
    
    vels=zeros(length(heights));
    velsU=zeros(length(heights));
 
    
    DYY=2000;
    DXX=2000;
        
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb

    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by

    
     fluxtim='timvar_fluxes';    
        
case '24Feb_orig'
    %    load in first dump of LEM run
	clear p
	h0=620;
    
    exname='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\5\04_02_24_2015_dmi.AED.txt';
    

    sound='combination';
    sound='dmi';
%    sound='sdla';

    %%%flag to use the sdla sounding and to put in the constant 5 ppmv in the TTL and stratosphere here
    % N.B. orginal LEM runs used the DMI sounding for everything including water vapour, except in the stratosphere where a constant 5 ppmv
    % was imposed based on hte SDLA sounding. However, there seems to have been a mistake in the ground height offset as the LEM and dat data only
    % match when the 620 m is left in.
    switch sound
    case 'sdla'
        readsdla; %reads the sdla sounding into sdla array
        temp=sdla(6,:)'+273.15; %K
        qvap=sdla(7,:)'/f;
        press=sdla(5,:)'*100;
        heights=sdla(2,:)';
        
    case 'combination'
        fid=fopen(exname,'rt');
		dat=fscanf(fid,'%g',[9 inf]);	       
		temp=dat(5,:)'+273.15; %temp
		press=dat(4,:)'*100; %pressure        
        
		heights=dat(3,:)'-h0; %height
        
        readsdla; %reads the sdla sounding into sdla array
        temp=sdla(6,:)+273.15; %K
        qvap=sdla(6,:)'/f;
        
    case 'dmi'
    
        fid=fopen(exname,'rt');
		dat=fscanf(fid,'%g',[9 inf]);
		sidat=(size(dat,2));
	
        
		temp=dat(5,:)'+273.15; %temp
		press=dat(4,:)'*100; %pressure
		rh=dat(6,:)'/100; %RH
        dir=dat(8,:)'*pi/180; %wind dir (rads)
        
        
		heights=dat(3,:)'-h0; %height
        
        qsw=satvapPress(temp,'buck2','liq',press,1)/f; 
        qvap=rh.*qsw;        
        
        ih=findheight(heights,14.72e3);
        qvap(ih+1:length(heights))=5/f; %have a constant 5 ppmv from 14.7 km upwards (where dmi sounding is equal to 5 ppmv)
    
    end

    %%%% sin(dir) = wind from east (no matter where on the CAST diagram) 
    %%%% cos(dir) = wind from north
    % dir=90-180 --- cos(180-dir)=-cos(dir) = wind from north since dirs are for where wind is from
    %            sin(180-dir)= sin(dir) = wind from east  etc.
    % so vels is then wind from N and velsU(=-dat(9,:)*sin(dir)) is then wind from west
    % then pos winds in j-dir are going from bottom of domain to top, which is wind travelling south
    % so vertical is arranged SN (top to bottom)
    % and horizontal is arranged WE (left to right)
    % BUT then the other direction is applied to get  values with dircomp degrees as zero (see below)
    
    
	dircomp=135 %135 %37.5
    vels=-dat(9,:)'.*sin(dir-dircomp*pi/180); %take component of speed at dircomp degrees - towards SW if dircomp=135
    velsU=-dat(9,:)'.*cos(dir-dircomp*pi/180); %take component of speed at dircomp degrees - towards SE if dircomp=135   
    
    %positive V speeds are therefore towards SW (up domain) so domain is SW at top and NE at bottom
    %positive U speeds are towards SE and to the right of the domain so left is NW and right is SE
    
    
    %%%%%%%%%%%%%%%%%%%%%% 3D flag - set properly!    
    i3d=1;
	%%%%%%%%%%%%%%%%%%%%%% 3D flag - set properly! 
	
        
    
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=20.25;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=heights(end)-h0; KGD(2)=250;
 
    
    DYY=2000;
    DXX=2000;
        
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb

    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
    
     
                
    
    
    %%%aside - code here to give the boundary layer 90% RH for a set of heights and potemps (as in namelist)
    hs=[0 225 900 1532 2134.8 2896.5 3531.3 4166.1 4800.9 5435.6 6070.4];
    potemps=[300.7 301.3 305.2 308.0 310.4 313.4 316.7 319.8 323.1 326.1 329.7]+1.5;
    
    f=1e6*28.97/18;
    
    
    
    
		extend=1;
		switch extend
		case 1
            
            h=15.8e3-620; %height at which to start the constant MR
            ih=findheight(heights,h)
        
        
            endZ=50e3; %final height that want to extend to
            Lz=length(heights);
            dz=heights(Lz)-heights(Lz-1);
            nz = (endZ-heights(Lz))/dz;
            endZi= Lz + nz;
            heights(Lz+1:endZi)=[heights(Lz)+dz:dz:heights(Lz)+nz*dz];
            
            
            temp(Lz+1:endZi) = temp(Lz); %for tropopause give a constant temperature            
            qvap(Lz+1:endZi) = qvap(Lz);
            vels(Lz+1:endZi) = vels(Lz);
            velsU(Lz+1:endZi) = velsU(Lz);
            
            Y0=press(Lz); %
			HSPAN=[heights(Lz+1) heights(endZi)]; %height range for temperature interpolation above ground
			[h,pdat] = ODE45( @hydrostatic,HSPAN,Y0,[],heights(Lz+1:endZi),temp(Lz+1:endZi) ); %solve hydrostatic equation - uses TSPAN to interpolate temp            
            press(Lz+1:endZi)=interp1(h,pdat,heights(Lz+1:endZi));
            
            
            
            imax=length(press);
            
		end



    
    
    
% this is stuff for modifying the sounding
% 
%             
%     for k=1:length(hs)
%         ih=findheight(heights,hs(k));
%         p=press(ih); %take pressure from LEM field
%         t=potemps(k)./(1000e2/p)^0.286; %use potemp supplied with hs vector
%         qsw=satvapPress(t,'buck2','liq',p,1)/f; %satvappress gives in ppmv if final flag=1
%         qnew(k)=0.9*qsw; %set to 90% RH
%     end
%     
% %values for removing dry layer from ~6.2 - 10.6 km (above ground)
%     hs2=[6705.2 7340.0 7974.8 8609.5 9244.3 9879.1 10513.9];
%     potemps2=[333.7 337.8 339.6 341.0 342.1 342.6 344.5];
%     
%     ih1=findheight(heights,6200);
% 	ih2=findheight(heights,10600);
%     rh2=rh;
%     rh2(ih1:ih2)=[];
%     
%     heights2=heights;
%     heights2(ih1:ih2)=[];
%     
%     rhnew=interp1(heights2,rh2,hs2); %interpolate to find values in between the removed heights
%     
%     for k=1:length(hs2)
%         ih=findheight(heights,hs2(k));
%         p=press(ih); %take pressure from LEM field
%         t=potemps2(k)./(1000e2/p)^0.286;
%         qsw=satvapPress(t,'buck2','liq',p,1)/f;
%         
%         qnew2(k)=rhnew(k)*qsw; %new mxing ratio from interpolated RH
%     end
    
     fluxtim='timvar_fluxes';
     
     
case 'TOGA-COARE_less_CIN'
    
    exname='c:\documents and settings\Login\my documents\Leeds_MMOCCA\gerry_tphidata_new.dat';
    fid=fopen(exname,'rt');
    fscanf(fid,'%s',[3 1]);
    dat=fscanf(fid,'%g',[3 inf]);
    
    tdat=dat(2,:)+273.15;
    pdat=dat(1,:)*100;
    tdew=dat(3,:)+273.15;
    
    qvap=satvappress(tdew,'goff','liq',pdat,1)/f;
    
    qisat=satvappress(tdat,'goff','ice',pdat,1)/f;
    
    
    rhi=qvap./qisat;
    
    
    Y0=pdat(1); %pressure at ground for Mile City
	TSPAN=[tdat]+273.15; %temperature range
	PSPAN=[pdat(1) pdat(end)]; %height range for temperature interpolation above ground
	[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],pdat,tdat); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
   
    heights=interp1(P,h,pdat);

    
	h0=0;

	press=pdat;
	temp=tdat;   
	
	vels=zeros([1 length(pdat)]);
    

    HDMP=2e3;  %height constant in damping layer equatiom
    XDMP=12e3; %height of the damping layer %%%%%%%%%%%%%%%%%%%%%check this!!!!!!!!!!!!!!!!!

    HGD(1)=1500;
    KGD(1)=21;
    KGD(2)=120;
    %HGD(end)=max height of sounding
    
    DYY=1000; %y res in m    

    
    RLAT=-2.0;
    RLON=156.0;
    IDAY=365;
    
    HOUR=0;
    local_diff=-7;  %??
    
 %   switch miles_case
 %   case 'more_stable_bl' %case where boundary layer was made less unstable (added inversion)
 
 fluxfactorS=302.8; %surface temp taken from the LEM test case 4 toga coare namlist
	fluxfactorL=1; %factor to multiply lat flux by (not need for saturated case)
 fluxtim='const_surface_sat';    
    
    %change BL - start with warmer surface layer, go slightly warmer than dry adiabat to p1,t1 and then dry ad to orig sounding
    Ts=tdat(1); %add 1 degC to lowest temp
    p1=980e2; t1=25+273;
    ih1=findheight(press,p1);
    ih1=findheight(heights,1000); %do dry adiabat up to 1km
%    temp(1:ih1)=interp1([press(1) press(ih1)],[Ts t1],press(1:ih1));
    thstart0=temp(1)*(1000e2/press(1))^0.286;
    temp(1:ih1)=thstart0./(1000e2./press(1:ih1)).^0.260; %reduced ad constant from 0.286 to make warmer than ad
        
    thstart=temp(ih1)*(1000e2/press(ih1))^0.286;
    
    
    i=ih1;
    temp_old=temp;
    
    %now do dry adiabat to meet oringal sounding
    temp_old(i)=temp_old(i)-1; %reduce old temp for the first pass through while so condition is true
    while temp(i)>temp_old(i)
        i=i+1;
        temp(i)=thstart/(1000e2/press(i))^0.286;
    end
    
       

        
    
    
case 'texas'
    exname='c:\documents and settings\Login\my documents\Leeds_MMOCCA\texas.dat';
    fid=fopen(exname,'rt');
%    fscanf(fid,'%s',[3 1]);
    dat=fscanf(fid,'%g',[4 inf]);
    
    tdat=dat(2,:);
    rh=dat(1,:); %relative humidity
    zdat=dat(3,:);
    
add_trop=1;
switch add_trop
case 1
    endZ=18e3; %height that want to extend to
    tropZ=12e3; %height where want tropopause to start   - NOTE orginally the tropopause was added to end of the sounding (e.g. tx_supp8 run) - 
    %it was lowered to 12 km on 29th April, 2008 so since LEM cloud looked to be going much higher than MAC3 cloud.
    Lz=findheight(zdat,tropZ);
    %Lz=length(zdat);
    dz=zdat(Lz)-zdat(Lz-1);
    nz = (endZ-zdat(Lz))/dz;
    endZi= Lz + nz;
    zdat(Lz+1:endZi)=[zdat(Lz)+dz:dz:zdat(Lz)+nz*dz];
    
    
    tdat(Lz+1:endZi) = tdat(Lz); %for tropopause give a constant temperature
    
end
            
    Y0=920e2; %pressure at ground 
	HSPAN=[zdat(1) zdat(end)]; %height range for temperature interpolation above ground
	[h,pdat] = ODE45(@hydrostatic,HSPAN,Y0,[],zdat,tdat); %solve hydrostatic equation - uses TSPAN to interpolate tempe
    
    pdat=interp1(h,pdat,zdat);
    press=pdat;
    
    qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
    
    
switch add_trop
case 1   
    qdat=rh(1:Lz).*qsat(1:Lz); 
    qdat(Lz+1:endZi) = qdat(Lz); %for tropopause give a constant q MR   
otherwise
    qdat=rh.*qsat;
end
    
	groundheight=667; 

%%%%%%%%%%%%%%%%%%%%%% 3D flag - set properly!    
    i3d=1;
%%%%%%%%%%%%%%%%%%%%%% 3D flag - set properly!    

	h0=groundheight;
    
    imin=1;
    imax=length(pdat);
    inds=[1:imax]; 
    
	heights=zdat;
	press=pdat;
	temp=tdat;
   
	qvap=qdat;  
   	
	vels=zeros([1 imax-imin+1]);
    velsU=vels;
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=22;
    local_diff=-6; %Houston, Texas standard time = -6hrs from GMT
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=heights(end)-HDMP; %height of the start of the damping layer

    HGD(1)=1500;
    KGD(1)=21;
    KGD(2)=150;
    %HGD(end)=max height of sounding
    
    DYY=1000; %y res in m
    DXX=DYY; %x res in m
    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
    
    RLAT=29.45;
    RLON=-95.23; %for Houston - check exact location
    IDAY=183; %???? - check exact date        
 
    fluxtim='timvar_fluxes';  
    
        
    %change BL - start with warmer surface layer, go slightly warmer than dry adiabat to p1,t1 and then dry ad to orig sounding
    Ts=tdat(1); %
%     p1=980e2; t1=25+273;
%     ih1=findheight(press,p1);
    ih1=findheight(heights,1000); %do dry adiabat up to 1km
%    temp(1:ih1)=interp1([press(1) press(ih1)],[Ts t1],press(1:ih1));

inv_type = 'supressed start temp';
%inv_type = 'same start temp';
switch inv_type
case 'same start temp'
    thstart0=temp(1)*(1000e2/press(1))^0.286;
    temp(1:ih1)=thstart0./(1000e2./press(1:ih1)).^0.260; %reduced ad constant from 0.286 to make warmer than ad

case 'supressed start temp'
    Tsupress=8; %amount to reduce surface temp by - then interpolate in potemp back to sounding (to create stability)
    thstart0=(temp(1)-Tsupress)*(1000e2/press(1))^0.286;
%    thend=thstart0*1.01;
	thend=(temp(ih1))*(1000e2/press(ih1))^0.286;
%comment the following line for normal BL
	temp(1:ih1)=(thstart0+[thend-thstart0].*([1:ih1]-1)/ih1)./(1000e2./press(1:ih1)).^0.286; %reduced ad constant from 0.286 to 0.26 make warmer than ad
end



    
    
    
case 'newmexico'
    exname='c:\documents and settings\Login\my documents\Leeds_MMOCCA\newmexico.dat';
    fid=fopen(exname,'rt');
%    fscanf(fid,'%s',[3 1]);
    dat=fscanf(fid,'%g',[4 inf]);
    
    i3d=0;  %%%%%%% MAKE SURE to set this properly - controls the 2d/3d field output
    
    tdat=dat(2,:);
    rh=dat(1,:); %relative humidity
    zdat=dat(3,:);
            
    Y0=821e2; %pressure at ground for New Mexico sounding
	HSPAN=[zdat(1) zdat(end)]; %height range for temperature interpolation above ground
	[h,pdat] = ODE45(@hydrostatic,HSPAN,Y0,[],zdat,tdat); %solve hydrostatic equation - uses TSPAN to interpolate tempe
    
    pdat=interp1(h,pdat,zdat);
    
    qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
    qdat=rh.*qsat; 
    
	groundheight=1660; 

	h0=groundheight;
    
    imin=1;
    imax=length(pdat);
    inds=[1:imax]; 
    
	heights=zdat;
	press=pdat;
	temp=tdat;
   
	qvap=qdat;  %take from the first dump of previous LEM run
   	
	vels=zeros([1 imax-imin+1]);
    pots=temp.*(1e5./press).^0.286;
    
    velsU=zeros(length(vels));
    
    %%%%%%%%%%%%%%%%%%% extrapolate the top of the sounding to allow a higher domain top
    itop=25e3; %new domain top to be extrapolated to
    L=length(heights); % end index so far
    dz=diff(heights(end-1:end)); %last height step
    extra=(itop-heights(end) )/dz + 1;
    
    heights(L+1:L+extra)=heights(end) + [1:extra]*ones(length(extra))*dz;    
    
    temp(L+1:L+extra)=interp1(heights(1:L),temp(1:L),heights(L+1:L+extra),'linear','extrap');
%    press(L+1:L+extra)=interp1(heights(1:L),press(1:L),heights(L+1:L+extra),'linear','extrap');
    pots(L+1:L+extra)=interp1(heights(1:L),pots(1:L),heights(L+1:L+extra),'linear','extrap');
 %   pots(L+1:L+extra)=pots(L);
    qvap(L+1:L+extra)=qvap(L);
    vels(L+1:L+extra)=vels(L);
    velsU(L+1:L+extra)=velsU(L);
    
    [h,p] = ODE45(@hydrostatic,[heights(L) heights(end)],press(L),[],heights,temp); %solve hydrostatic equation - solves pressure from temp array    
    press(L+1:L+extra)=interp1(h,p,heights(L+1:L+extra)); %interpolate for required values in heights array
    
    HOUR=22;
    local_diff=-7; %New Mexico standard time = -7hrs from GMT
    HDMP=7e3;  %height constant in damping layer equatio
    XDMP=heights(end)-HDMP; %height of the damping layer

    HGD(1)=1500;
    KGD(1)=21;
    KGD(2)=120;
    %HGD(end)=max height of sounding
    
    DYY=1000; %y res in m
    DXX=DYY;
    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
    
    RLAT=35.07;
    RLON=-106.4; %for Albuquerque
    IDAY=183; %???? - check exact date        
 
    fluxtim='timvar_fluxes';   
    
%%%change BL - start with warmer surface layer, go slightly warmer than dry adiabat to p1,t1 and then dry ad to orig sounding
    Ts=tdat(1); 
%     p1=980e2; t1=25+273;
%     ih1=findheight(press,p1);
    ih1=findheight(heights,1000); %do dry adiabat up to 1km
%    temp(1:ih1)=interp1([press(1) press(ih1)],[Ts t1],press(1:ih1));
    Tsupress=4; %amount to reduce surface temp by - then interpolate in potemp back to sounding (to create stability)
    thstart0=(temp(1)-Tsupress)*(1000e2/press(1))^0.286;
%    thend=thstart0*1.01;
	thend=(temp(ih1))*(1000e2/press(ih1))^0.286;
%comment the following line for normal BL
	temp(1:ih1)=(thstart0+[thend-thstart0].*([1:ih1]-1)/ih1)./(1000e2./press(1:ih1)).^0.286; %reduced ad constant from 0.286 to 0.26 make warmer than ad
    

%	temp(2:ih1)=thstart0.*(1:ih1)*./(1000e2./press(2:ih1)).^0.286; %reduced ad constant from 0.286 to 0.26 make warmer than ad


    
    
case 'TOGA-COARE'
    
    exname='c:\documents and settings\Login\my documents\Leeds_MMOCCA\gerry_tphidata_new.dat';
    fid=fopen(exname,'rt');
    fscanf(fid,'%s',[3 1]);
    dat=fscanf(fid,'%g',[3 inf]);
    
    tdat=dat(2,:)+273.15;
    pdat=dat(1,:)*100;
    tdew=dat(3,:)+273.15;
    
    qvap=satvappress(tdew,'goff','liq',pdat,1)/f;
    
    qisat=satvappress(tdat,'goff','ice',pdat,1)/f;
    
    rhi=qvap./qisat;
    
    
    Y0=pdat(1); %pressure at ground for Mile City
	TSPAN=[tdat]+273.15; %temperature range
	PSPAN=[pdat(1) pdat(end)]; %height range for temperature interpolation above ground
	[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],pdat,tdat); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
   
    heights=interp1(P,h,pdat);
    
	h0=0;

	press=pdat;
	temp=tdat;   
	
	vels=zeros([1 length(pdat)]);
    pots=temp.*(1e5./press).^0.286;        
    

    HDMP=2e3;  %height constant in damping layer equatiom
    XDMP=12e3; %height of the damping layer %%%%%%%%%%%%%%%%%%%%%check this!!!!!!!!!!!!!!!!!

    HGD(1)=1500;
    KGD(1)=21;
    KGD(2)=120;
    %HGD(end)=max height of sounding
    
    DYY=1000; %y res in m    

    
    RLAT=-2.0;
    RLON=156.0;
    IDAY=365;
    
    HOUR=0;
    local_diff=-7;  %??
    
 %   switch miles_case
 %   case 'more_stable_bl' %case where boundary layer was made less unstable (added inversion)
 
 fluxfactorS=302.8; %surface temp taken from the LEM test case 4 toga coare namlist
	fluxfactorL=1; %factor to multiply lat flux by (not need for saturated case)
 fluxtim='const_surface_sat';
 
case '24Feb - hot BL'
    %    load in first dump of LEM run
	clear p
	h0=620;
    
	heights=GridDan(1).Z;
	press=GridDan(1).PREFN;
	temp=TempLES(GridDan(1));
    
   
	qvap=GridDan(1).OLQBAR(:,1);  %take from the first dump of previous LEM run
	vels=GridDan(1).VBAR;
    
    qsw=satvapPress(temp,'buck2','liq',press,1)/f; 
    rh=qvap./qsw;
    
    
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=1;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=heights(end)-h0; KGD(2)=250;
 
    
    DYY=1000;
        
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb
    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=2; %factor to multiply lat flux by
    
    
    %%%aside - code here to give the boundary layer 90% RH for a set of heights and potemps (as in namelist)
    hs=[1000];
    ih=findheight(heights,hs);

    
    f=1e6*28.97/18;
    
    Tpert=1;
    for k=1:ih
        p=press(k); %take pressure from LEM field
        temp(k)=temp(k) + Tpert; %use potemp supplied with hs vector
        pots(k)=temp(k) * (1000e2/p)^0.286;
        qsw=satvapPress(temp(k),'buck2','liq',p,1)/f; %satvappress gives in ppmv if final flag=1
        qvap(k)=0.9*qsw; %set to 90% RH
    end

    %create inversion
    hs=[1500];
    ih2=findheight(heights,hs);
    Tinv=temp(ih);
    temp(ih:ih2)=Tinv; %use constant temp
    

 fluxtim='timvar_fluxes';

    
case '24Feb'
    %    load in first dump of LEM run
	clear p
	h0=620;
    
	heights=GridDan(1).Z;
	press=GridDan(1).PREFN;
	temp=TempLES(GridDan(1));
    
   
	qvap=GridDan(1).OLQBAR(:,1);  %take from the first dump of previous LEM run
	vels=GridDan(1).VBAR;
    
    qsw=satvapPress(temp,'buck2','liq',press,1)/f; 
    rh=qvap./qsw;
    
    
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=20.25;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=heights(end)-h0; KGD(2)=250;
 
    
    DYY=1000;
        
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb

    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
    
    
    %%%aside - code here to give the boundary layer 90% RH for a set of heights and potemps (as in namelist)
    hs=[0 225 900 1532 2134.8 2896.5 3531.3 4166.1 4800.9 5435.6 6070.4];
    potemps=[300.7 301.3 305.2 308.0 310.4 313.4 316.7 319.8 323.1 326.1 329.7]+1.5;
    
    f=1e6*28.97/18;
    


            
    for k=1:length(hs)
        ih=findheight(heights,hs(k));
        p=press(ih); %take pressure from LEM field
        t=potemps(k)./(1000e2/p)^0.286; %use potemp supplied with hs vector
        qsw=satvapPress(t,'buck2','liq',p,1)/f; %satvappress gives in ppmv if final flag=1
        qnew(k)=0.9*qsw; %set to 90% RH
    end
    
%values for removing dry layer from ~6.2 - 10.6 km (above ground)
    hs2=[6705.2 7340.0 7974.8 8609.5 9244.3 9879.1 10513.9];
    potemps2=[333.7 337.8 339.6 341.0 342.1 342.6 344.5];
    
    ih1=findheight(heights,6200);
	ih2=findheight(heights,10600);
    rh2=rh;
    rh2(ih1:ih2)=[];
    
    heights2=heights;
    heights2(ih1:ih2)=[];
    
    rhnew=interp1(heights2,rh2,hs2); %interpolate to find values in between the removed heights
    
    for k=1:length(hs2)
        ih=findheight(heights,hs2(k));
        p=press(ih); %take pressure from LEM field
        t=potemps2(k)./(1000e2/p)^0.286;
        qsw=satvapPress(t,'buck2','liq',p,1)/f;
        
        qnew2(k)=rhnew(k)*qsw; %new mxing ratio from interpolated RH
    end
    
     fluxtim='timvar_fluxes';

    
case 'ecmwf high top 2'
    %run read_profiles_netcdf and load in first dump from LEM
	clear p
    parr=flipud(ecmwf(1).p);
    
    hinterp=interp1(GridDan(1).PREFN,GridDan(1).Z,parr*100); 
    %find heights of pressure points within LEM data by interpolation
    inan=isnan(hinterp); %some will be NaN due to LEM pressure data being out of range
    inan_zero=find(inan==0); %find first proper point (lower levels < 620m not in LEM output)
    inan_one=find(inan==1);
    inan_one=inan_one(inan_one>inan_zero(1));
    ips=inan_one(1);
        

    
    
    Tarr=ecmwf(1).t(:,1,2);
  
    z = hinterp(1:ips-1);
	z(ips-1:length(hinterp)) = hyd_pressure_height( parr(ips-1:end),Tarr(ips-1:end),hinterp(ips-1) );  %get pressure height profile for ecmwf data

        
	%iz1=49; %height index from which to use the ECMWF data
	iztop=60;
	%iend=length(temp); %end point of LEM data
    
    hmatch=27e3;
    hmatch=25.4e3;
    
    iz1=findheight(z,hmatch); %height index from which to use the ECMWF data
    iend=findheight(GridDan(1).Z,hmatch); %end point of LEM data
    if z(iz1)<GridDan(1).Z(iend);   
       iz1=iz1+1;
    end
    
    heights=GridDan(1).Z(1:iend);
	press=GridDan(1).PREFN(1:iend);
	tLES=TempLES(GridDan(1));
    temp=tLES(1:iend);
	qvap=GridDan(1).OLQBAR(1:iend,1);  %take from the first dump of previous LEM run
	vels=GridDan(1).VBAR(1:iend);
    
	%%% now for ecmwf bit at the top
        
    ih=findheight(heights+620,17e3);
    qvap(ih:iend+iztop-iz1+1)=qvap(ih); %make the vapour in the strat equal to the constant TTL value
	
	heights(iend+1:iend+iztop-iz1+1)=z(iz1:iztop)';
	press(iend+1:iend+iztop-iz1+1)=parr(iz1:iztop)*100;
    
    irelax=1;
    
    if irelax==1
        ii=iend;
        tdiff_new=9e9;
        tdiff_old=9e10;
        while tdiff_new < tdiff_old
            ii=ii+1;
            dtdz=diff(temp(end-1:end))./diff(heights(ii-2:ii-1));
            temp(end+1)=temp(end)+dtdz*1*(heights(ii)-heights(ii-1));
            tdiff_old=tdiff_new;
            tdiff_new=abs(temp(end) - Tarr(ii-iend+iz1-1))       
        end
    
        temp(ii:iend+iztop-iz1+1)=Tarr(ii-iend+iz1-1:iztop);
    else
 	    temp(iend+1:iend+iztop-iz1+1)=Tarr(iz1:iztop);
    end

    	

	%qvap(iend+1:iend+iztop-iz1+1)=ecmwf(1).q(iz1:iztop,1,2); %not using qvap from ecmwf, just assuming constant
	%vapour above TTL
        
	u=ecmwf(1).u(iz1:iztop,1,2);
	v=ecmwf(1).v(iz1:iztop,1,2);
	bear=bearing(u,v);
	mag=sqrt(u.^2+v.^2);
    
    [velsU vels_ecmwf]=windcomp(mag,bear,135); %gets components for domain orientation of 135 degrees. vels is the output required for 2-d runs
    vels(iend+1:iend+iztop-iz1+1)=-vels_ecmwf; % need to make it negative to match the LEM winds (not sure exactly why but get very good match like this)
    pots=temp.*(1e5./press).^0.286;
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=30556; KGD(2)=250;
    KGD(3)=330;
    
    HOUR=20.25;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    DYY=1000;
    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
        
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb
    
     fluxtim='timvar_fluxes';

    
case 'ecmwf high top'
    %changed this with 'ecmwf high top 2' since was producing a large spike of dth(potemp)/dz at the cross over
    %point that may have been causing a stable layer
    %run soundings/read_profiles_netcdf and load in first dump from LEM
	clear p
	iz1=49;
	iztop=60;
	%izlem=findheight(GridDan(1).Z,z(iz1)); %z=ecmwf height
	h0=620;
    
	heights=GridDan(1).Z;
	press=GridDan(1).PREFN;
	temp=TempLES(GridDan(1));
   
	qvap=GridDan(1).OLQBAR(:,1);  %take from the first dump of previous LEM run
   
	
	vels=GridDan(1).VBAR;
    
	%now for ecmwf bit at the top
	iend=length(temp);
    	
	hyd_pressure_height;  %get pressure height profile for ecmwf data
    
    ih=findheight(heights+620,17e3);
    qvap(ih:iend+iztop-iz1+1)=qvap(ih); %make the vapour in the strat equal to the constant TTL value
	
	heights(iend+1:iend+iztop-iz1+1)=z(iz1:iztop)'-620;
	press(iend+1:iend+iztop-iz1+1)=parr(iz1:iztop)*100;
	temp(iend+1:iend+iztop-iz1+1)=Tarr(iz1:iztop);
	%qvap(iend+1:iend+iztop-iz1+1)=ecmwf(1).q(iz1:iztop,1,2); %not using qvap from ecmwf, just assuming constant
	%vapour above TTL
    
	u=ecmwf(1).u(iz1:iztop,1,2);
	v=ecmwf(1).v(iz1:iztop,1,2);
	bear=bearing(u,v);
	mag=sqrt(u.^2+v.^2);
    
    [velsU vels_ecmwf]=windcomp(mag,bear,135); %gets components for domain orientation of 135 degrees. vels is the output required for 2-d runs
    vels(iend+1:iend+iztop-iz1+1)=-vels_ecmwf; % need to make it negative to match the LEM winds (not sure exactly why but get very good match like this)
    pots=temp.*(1e5./press).^0.286;
    
    HGD(1)=1500; KGD(1)=21;
    HGD(2)=30556; KGD(2)=250;
    KGD(3)=330;
    
    HOUR=20.25;
    local_diff=-3;
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=24e3; %height of the damping layer
    
    DYY=1000;

    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
    
    RLAT=-22.36;
    RLON=-49.03;
    IDAY=55; %55=24th feb
    
     fluxtim='timvar_fluxes';

    
case 'Miles City'
%    run Sound_MileCity_MMOCCA first
	clear p
	h0=1000;
    imax=length(pr(3).p(:,1));
    imin=2;          %%%%%%%%%%%%  NOTE - =2 means ignoring first point as a test. Change back to 1 if want full sounding
    inds=[imin:imax]; 
    
	heights=pr(3).p(inds,1);
	press=pr(3).p(inds,2)*100;
	temp=pr(3).p(inds,3)+273.15;
   
	qvap=pr(3).p(inds,10);  %take from the first dump of previous LEM run
   
	
	vels=zeros([1 imax-imin+1]);
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=22;
    local_diff=-7; 
    HDMP=2e3;  %height constant in damping layer equatio
    XDMP=12e3; %height of the damping layer

    HGD(1)=1500;
    KGD(1)=21;
    KGD(2)=120;
    %HGD(end)=max height of sounding
    
    DYY=500; %y res in m
    
    fluxfactorS=1; %factor to multiply sens flux by
	fluxfactorL=1; %factor to multiply lat flux by
    
    RLAT=46.4;
    RLON=-105.8;
    IDAY=183; %approx July ???? - check exact date
    
    
 %   switch miles_case
 %   case 'more_stable_bl' %case where boundary layer was made less unstable (added inversion)
 
  fluxtim='timvar_fluxes';

 case 'Miles City_inversion'
%    run Sound_MileCity_MMOCCA first
    Sound_MileCity_MMOCCA
    i3d=1;
    
	clear p
	h0=1000;
    imax=length(pr(3).p(:,1));
    imin=1;          %%%%%%%%%%%%  NOTE - =2 means ignoring first point as a test. Change back to 1 if want full sounding
    inds=[imin:imax]; 
    
	heights=pr(3).p(inds,1);
	press=pr(3).p(inds,2)*100;
	temp=pr(3).p(inds,3)+273.15;
   
	qvap=pr(3).p(inds,10); 
    
add_trop=1;
switch add_trop
case 1
    endZ=18e3; %final height that want to extend to
    Lz=length(heights);
    dz=heights(Lz)-heights(Lz-1);
    nz = (endZ-heights(Lz))/dz;
    endZi= Lz + nz;
    heights(Lz+1:endZi)=[heights(Lz)+dz:dz:heights(Lz)+nz*dz];
    
    
    temp(Lz+1:endZi) = temp(Lz); %for tropopause give a constant temperature
    
    qvap(Lz+1:endZi) = qvap(Lz);
    
    Y0=press(Lz); %pressure at ground for New Mexico sounding
	HSPAN=[heights(Lz+1) heights(endZi)]; %height range for temperature interpolation above ground
	[h,pdat] = ODE45( @hydrostatic,HSPAN,Y0,[],heights(Lz+1:endZi),temp(Lz+1:endZi) ); %solve hydrostatic equation - uses TSPAN to interpolate tempe
    
    press(Lz+1:endZi)=interp1(h,pdat,heights(Lz+1:endZi));
    
    imax=length(press);
    
end
            
    
    
    

    
    %this next line was included for the NM_tests and NM_tests_noheat runs
   % qvap(1)=qvap(2); %qv values seem very high compared to those of sounding from paper (Farley, 1992)
                     %so are using the value from second point, which seems more comparable (=11 g/kg in that
                     %paper)
    
   
% 	ip830=findheight(press,830e2); %index for 830 mb where want to start inversion
%     
%     %temperature that want to push ground to
%     groundTemp=28+273.15; %as was the case for NM_2e-3_240ccn
%     groundTemp=30+273.15; %changed to try and get stronger convection
%     %remove this line if don't want inversion
%     temp(1:ip830)=interp1([heights(1) heights(ip830)],[groundTemp temp(ip830)],heights(1:ip830));
%     
	vels=zeros([1 imax-imin+1]);
    velsU=vels;
    pots=temp.*(1e5./press).^0.286;
    
    HOUR=11+7; %start at 11am for plenty of spin up opportunity from surface fluxes
    HOUR=14.67+7; %actual time of sounding - Farley paper says cloud developed 1.5 hours later
%    HOUR=16+7; %actual time of sounding - Farley paper says cloud developed 1.5 hours later    
    
    local_diff=-7; 
    HDMP=2e3;  %height constant in damping layer equation
    XDMP=12e3; %height of the start of the damping layer

    HGD(1)=1500;
    KGD(1)=21;
    KGD(2)=150;
    %HGD(end)=max height of sounding
    
    DYY=1000; %y res in m
    DXX=DYY;
    
%    fluxfactorS=1; %factor to multiply sens flux by
%	 fluxfactorL=1; %factor to multiply lat flux by
    
    fluxfactorS=0.4; %factor to multiply sens flux by - reduced them relative to 24th Feb case so that convection doesn't start everywhere
                     %really just want one cloud
	fluxfactorL=1; %factor to multiply lat flux by
    
    RLAT=46.4;
    RLON=-105.8;   
    IDAY=183; %approx July ???? - check exact date
    
 
     fluxtim='timvar_fluxes';
     
     
ih1=findheight(heights,1000); %do dry adiabat up to 1km     
     
     inv_type = 'supressed start temp';
%inv_type = 'same start temp';
switch inv_type
case 'same start temp'
    thstart0=temp(1)*(1000e2/press(1))^0.286;
    temp(1:ih1)=thstart0./(1000e2./press(1:ih1)).^0.260; %reduced ad constant from 0.286 to make warmer than ad

case 'supressed start temp'
    Tsupress=6; %amount to reduce surface temp by - then interpolate in potemp back to sounding (to create stability)
    thstart0=(temp(1)-Tsupress)*(1000e2/press(1))^0.286;
%    thend=thstart0*1.01;
	thend=(temp(ih1))*(1000e2/press(ih1))^0.286;
%comment the following line for normal BL
	temp(1:ih1)=(thstart0+(thend-thstart0).*([1:ih1]'-1)/ih1)./(1000e2./press(1:ih1)).^0.286; %reduced ad constant from 0.286 to 0.26 make warmer than ad
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
end


if isounding==1
%for sounding program
    zdat=heights;
    pdat=press;
    tdat=temp;
    qdat=qvap;
    qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
    
    
    pots=temp.*(1e5./press).^0.286; 
    rh=qdat./qsat;
    
    


	%now run writesound
    writeSound_cprog_conversion_2(comp,HOUR,local_diff,XDMP,HDMP,HGD,KGD,DXX,DYY,RLAT,RLON,...
        heights,pots,press,qvap,vels,velsU,fluxfactorS,fluxfactorL,IDAY,fluxtim,i3d);
    %note for velocities need to specify the u and v components in velsU and vels
    
end

disp('*****  Finished soundings  *****');
    