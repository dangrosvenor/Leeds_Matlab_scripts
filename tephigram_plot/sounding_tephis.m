%example plot of tephigram

    f=1e6*28.97/18;
    
plot_case='MilesCity';
plot_case='lem';
plot_case='24thFebDMI_orig';
% plot_case='24thFebDMI_heated_3d_dump4';
% plot_case='MilesCity_ecmwf';
% plot_case='MilesCity_altered_sounding';
%plot_case='24thFebDMI_maxpert';
%plot_case='TOGA';
%plot_case='TOGA_testcase4';
plot_case='newmexico';
plot_case='texas';
%plot_case='EMM';
%plot_case='emm_storm_case2';

switch plot_case
case 'emm_storm_case2'
    emm_storm_case2; %writes values needed
	
	nz=length(ps_Pa);
    inds=[1:nz];
	
    zdat=zRead(inds)*1000;
	pdat=ps_Pa(inds)*100;    
	tdat=temp_K(inds)+273.15;
    
    qdat=satvappress(td_K+273.15,'goff','liq',pdat,1)/f;   
    qsat=satvappress(tdat,'goff','liq',pdat,1)/f; 
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1); %1 flag at end says to plot on top of what's already there
    
%    ps_Pa(1) = 642.86; temp_K(1) = 5.77; td_K(1) = 1.44;			zRead(1) = 3.;
    
    %check orientation of wind
  %  plot_wind_barbs(pdat,pr(1).p(:,5),pr(1).p(:,6)); %plot_wind_barbs(p,u,v)
  
case 'EMM'
    soundings; %run this to read in sounding (set flag in there to the case required)
	
	nz=length(press);
    inds=[1:nz];
	
    zdat=heights(inds);
	pdat=press(inds);    
	qdat=qvap(inds);
	tdat=temp(inds);
    qsat=qsw(inds);
	
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1); %1 flag at end says to plot on top of what's already there
    
    %check orientation of wind
  %  plot_wind_barbs(pdat,pr(1).p(:,5),pr(1).p(:,6)); %plot_wind_barbs(p,u,v)
    
case 'texas'
    soundings; %run this to read in sounding (set flag in there to the case required)
	
	nz=length(pdat);
    inds=[1:nz];
	
    zdat=zdat(inds);
	pdat=pdat(inds);    
	qdat=qdat(inds);
	tdat=tdat(inds);
	
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1); %1 flag at end says to plot on top of what's already there
    
    %check orientation of wind
    plot_wind_barbs(pdat,pr(1).p(:,5),pr(1).p(:,6)); %plot_wind_barbs(p,u,v)
        
case 'newmexico'
    soundings; %run this to read in sounding (set flag to 'newmexico')
	
	nz=length(pdat);
    inds=[1:nz];
	
    zdat=zdat(inds);
	pdat=pdat(inds);    
	qdat=qdat(inds);
	tdat=tdat(inds);
	qsat=qsat(inds);
	
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,0);
    
    

case 'TOGA_testcase4'
    TOGA_COARE_DATA_testcase4; %puts toga coare data into arrays (as given in LEM test case 4 - runTEST4_2.3.f)
    
    tdat=TITOG;
    pdat=PSTOG;    
    qvap=QITOG;
        
    Y0=pdat(1); %pressure at ground 
	TSPAN=[tdat]+273.15; %temperature range
	PSPAN=[pdat(1) pdat(end)]; %height range for temperature interpolation above ground
	[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],pdat,tdat); %solve hydrostatic equation (dz/dp as func p)
    
    zdat=interp1(P,h,pdat); %interpolate onto required pressure grid as in pdat
    
	groundheight=0; %oceanic convection
	
	nz=length(pdat);
    inds=[1:nz];
	
    zdat=h(inds);
	pdat=pdat(inds);    
	qdat=qvap(inds);
	tdat=tdat(inds);


	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1);
    
 

case 'TOGA'
    exname='c:\documents and settings\Login\my documents\Leeds_MMOCCA\gerry_tphidata_new.dat';
    fid=fopen(exname,'rt');
    fscanf(fid,'%s',[3 1]);
    dat=fscanf(fid,'%g',[3 inf]);
    
    tdat=dat(2,:)+273.15;
    pdat=dat(1,:)*100;
    tdew=dat(3,:)+273.15;
    
    qvap=satvappress(tdew,'goff','liq',pdat,1)/f;
    
    
    Y0=pdat(1); %pressure at ground for Mile City
	TSPAN=[tdat]+273.15; %temperature range
	PSPAN=[pdat(1) pdat(end)]; %height range for temperature interpolation above ground
	[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],pdat,tdat); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
        
	groundheight=0; %oceanic convection
	
	nz=length(pdat);
    inds=[1:nz];
	
    zdat=h(inds);

	pdat=pdat(inds);    
	qdat=qvap(inds);
	tdat=tdat(inds);


	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1);
                
    
case '24thFebDMI_maxpert'
	groundheight=620;	
	
	nz=length(GridDan(idir).Z);
	zdat=GridDan(idir).Z(1:nz);
	pdat=GridDan(idir).PREFN(1:end);
	qdat=GridDan(idir).OLQBAR(1:nz,1);
%	qdat(1)=qdat(2);

clear diff
% dy=diff(GridDan(1).Y1(1:2));
% D=14e3;
% %D=5e3;
% ny=D/dy;
ny=7;

    plotline='maxT';
    plotline='background';
   
   idir=3;
    
    switch plotline
    case 'maxT'
        
        [th thind]= max( TwoDDan(idir).TH2( :, 2:end-1 ),[],2 )  ;
        %    th=th+GridDan(idir).THREF;
        %    qdat=squeeze(   TwoDDan(idir)Q01( [end-ny:end 2:2+ny] , [end-ny:end 1:ny+1] , : )  ))  );
        
        %    th=GridDan(idir).THREF + squeeze ( TH1(151,1,:) ) ;
        %    tdat=th ./ (1000e2./pdat).^0.286 ;    
        %    qdat=squeeze ( Q01(151,1,:) ) ;
        
%         D=1e3;
%         X=round(D/diff(GridDan(idir).Y1(1:2))/2);
%         X=3;
        
        k=1;
        
        D=7e3;
%        D=1e3;
        
        [yinds,xinds]=find_inds_feature(GridDan(idir),GridDan(idir).Y1(thind(k)),GridDan(idir).Y1(thind(k)),D);

        
        for k=1:length(th)
            
                                    
            th(k)=mean(TwoDDan(idir).TH2(k,yinds));
%            qdat(k)=squeeze(  TwoDDan(idir).Q(k,thind(k),1)  );
            qdat(k)=mean(TwoDDan(idir).Q(k,yinds,1));
            
        end
        
        tdat=th ./ (1000e2./pdat).^0.286 ;
        
        add_to_previous=0;


        
  
    case 'background'
        th=squeeze(   median(TwoDDan(idir).TH2 , 2 )  );
        tdat=th ./ (1000e2./pdat).^0.286 ;
        qdat=squeeze(  median(TwoDDan(idir).Q(:,thind) ,2) );
        
        add_to_previous=1;
        
    end
    



    
	f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;

    inds=[3:length(pdat)];

	tdat=tdat(inds);
	zdat=zdat(inds);
	qdat=qdat(inds);
	qsat=qsat(inds);
	pdat=pdat(inds);
    
    if add_to_previous==0        
        [tp,h]=plot_tephi(-40+273,50+273,153,5000);
        plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,0);
    else
        plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,2);
    end

    
    
  
    
    
case 'MilesCity_altered_sounding'
	groundheight=1000;
	%[tp,h]=plot_tephi(-40+273,50+273,153,5000);
	
	nz=length(press);
    inds=[1:nz];
	zdat=heights;

	pdat=press(inds);    
	qdat=qvap(inds);
	tdat=temp(inds);

    f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1);
    
    %Miles City is 7 hours behind UTC.
    
case 'MilesCity_ecmwf'
	groundheight=1000;
	%[tp,h]=plot_tephi(-40+273,50+273,153,5000);
	
    iground=findheight(ecmwf(1).z,groundheight);
	nz=length(ecmwf(1).z);
    inds=[iground:nz];
	zdat=ecmwf(1).z(inds);

	pdat=flipud(ecmwf(1).p)*100; %convert to Pa
    pdat=pdat(inds);
	qdat=ecmwf(1).q(2,inds,1,1)';
	tdat=ecmwf(1).t(2,inds,1,1)';

    f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat);
    
    %Miles City is 7 hours behind UTC. Have times of 18 (19th July) and 00 UTC (20th) = 11am LT and 17:00 LT. So miss the midday sun 
    %heating and perhaps explains why temperatures are quite low in comparison to 14:40 sounding.
    
case '24thFebDMI_heated_3d_dump4'
	groundheight=620;	
	
	nz=length(GridDan(idir).Z);
	zdat=GridDan(idir).Z(1:nz);
	pdat=GridDan(idir).PREFN(1:end);
	qdat=GridDan(idir).OLQBAR(1:nz,1);
%	qdat(1)=qdat(2);

clear diff
% dy=diff(GridDan(1).Y1(1:2));
% D=14e3;
% %D=5e3;
% ny=D/dy;
ny=7;

pert_back='pert2';
switch pert_back
case 'pert'

	th=GridDan(idir).THREF + squeeze(   mean(mean(TH1( [end-ny:end 2:2+ny] , [end-ny:end 1:ny+1] , : )  ))  );
    tdat=th ./ (1000e2./pdat).^0.286 ;
    qdat=squeeze(   mean(mean(Q01( [end-ny:end 2:2+ny] , [end-ny:end 1:ny+1] , : )  ))  );
    
%    th=GridDan(idir).THREF + squeeze ( TH1(151,1,:) ) ;
%    tdat=th ./ (1000e2./pdat).^0.286 ;    
%    qdat=squeeze ( Q01(151,1,:) ) ;
    
case 'back'
	th=GridDan(idir).THREF + squeeze(   median(median(TH1(:,:,:)  ))  );
    tdat=th ./ (1000e2./pdat).^0.286 ;
    qdat=squeeze(   median(median(Q01(:,:,:)  ))  );
    
case 'pert2'

	th=squeeze( mean(TwoD.TH2( :, [end-ny:end-1 2:2+ny]  ) ,2 )  );
    tdat=th ./ (1000e2./pdat).^0.286 ;
    qdat=squeeze(  mean(TwoD.Q( : , [end-ny:end-1 2:2+ny] ,1 ),2  )  );
  
    
end



    
	f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;

    ind_ignore=2;
    ignore_lowest=1;
    if ignore_lowest==1
     	tdat=tdat(ind_ignore:end);
     	zdat=zdat(ind_ignore:end);
     	qdat=qdat(ind_ignore:end);
     	qsat=qsat(ind_ignore:end);
     	pdat=pdat(ind_ignore:end);
    end

add_to_previous=1;    
if add_to_previous==0        
    [tp,h]=plot_tephi(-40+273,50+273,153,5000);
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,0);
else
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1);
end



case 'lem'
    idir=1;
    
	groundheight=1000;
	groundheight=620;
	
	%%[tp,h]=plot_tephi(-40+273,50+273,153,5000);

	%[tp,h]=plot_tephi(-70+273,50+273,273-100,20e2);

	nz=length(GridDan(idir).Z);
    
    inds=[2:nz];
    
	zdat=GridDan(idir).Z(inds);
	pdat=GridDan(idir).PREFN(inds);
	qdat=GridDan(idir).OLQBAR(inds,1);
%	qdat(1)=qdat(2);
	tdat=tempLES(GridDan(idir));
    tdat=tdat(inds);
    
	%tdat=tdat(2:end);
	f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
% 	isup=find(qdat./qsat > 1);
%     if length(isup)>=1	
% 		qdat(isup)=qsat(isup);
% 		qdat(isup(1))=qsat(isup(1))*0.98;
%     end
	
	%plot_tephi_data(tdat,pdat,qdat,qsat);
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat+groundheight,1);

case 'MilesCity'
	groundheight=1000;	
	[tp,h]=plot_tephi(-40+273,50+273,153,5000);
	
	nz=length(HSPAN);
    inds=[2:nz];
	zdat=HSPAN(inds)+groundheight;
	pdat=pr(3).p(inds,4);
	qdat=pr(3).p(inds,10);
	tdat=pr(3).p(inds,3)+273.15;

    f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
	isup=find(qdat./qsat > 1);
	
	qdat(isup)=qsat(isup);	
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat);

case '24thFebDMI_orig'
	groundheight=620;	
	
    zinds=1:600;
	
	zdat=pr(1).p(zinds,1);
	pdat=pr(1).p(zinds,2)*100;
	qdat=pr(1).p(zinds,10);
	tdat=pr(1).p(zinds,3)+273.15;

    f=1e6*28.97/18;
	qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
	
    
%    [tp,h]=plot_tephi(-40+273,50+273,153,5000);
	plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1);
    
    %check orientation of wind
    u=pr(1).p(zinds,6)/0.514; %convert to knots (one knot = 0.514 m/s)
    v=pr(1).p(zinds,5)/0.514; %convert to knots (one knot = 0.514 m/s)
    
    plot_wind_barbs(pdat,u,v); %plot_wind_barbs(p,u,v)
    
    
end