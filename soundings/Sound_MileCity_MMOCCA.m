%read in AED.txt style sounding
clear dire

dire(3).d='c:\documents and settings\Login\my documents\Leeds_MMOCCA\Milecityo.txt';


%for i=1:size(dire,2)

for i=3:3
	pr(i).p=[]; %clear arrays
        
	fid=fopen(dire(i).d,'rb');
	
	dat=fscanf(fid,'%g',[5 inf]);
	sidat=(size(dat,2));
	
	pr(i).p(1:sidat,3)=dat(3,:)' - 273.15; %temp (degC)
	%pr(i).p(1:sidat,2)=dat(4,:)'; %pressure
	pr(i).p(1:sidat,4)=100*dat(1,:)'; %RH in percent
	%pr(i).p(1:sidat,5)=dat(9,:)'.*cos( dat(8,:)'.*pi./180 ); %vwind from N I think
	%pr(i).p(1:sidat,6)=dat(9,:)'.*sin( dat(8,:)'.*pi./180 );%uwind from E
	pr(i).p(1:sidat,5)=0; %vwind from N I think
	pr(i).p(1:sidat,6)=0;%uwind from E
	
	pr(i).p(1:sidat,1)=dat(4,:)'; %height - add ground height
	
	pr(i).p(1:sidat,7)=log10(pr(i).p(1:sidat,3).*(1000./pr(i).p(1:sidat,2)).^0.286);%log theta
	%pr(i).p(1:sidat,8)=dat(9,1:sidat)'; %wind mag
	%pr(i).p(1:sidat,9)=dat(8,1:sidat)'; %wind dir
	
    %HabvMSL=838.12; %worked out from hydrostatic for 921 hPa
    
	%Y0=1013.25*100; %mean sea level pressure
    Y0=921*100; %pressure at ground for Mile City
	TSPAN=pr(i).p(:,3)+273.15; %temperature range
	HSPAN=pr(i).p(:,1); %height range for temperature interpolation above ground
	[h,P] = ODE45(@hydrostatic,HSPAN,Y0,[],HSPAN,TSPAN); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H
	
	pr(i).p(1:sidat,2)=P/100; %needs to be in hPa for AED stlye sounding that c program reads in
	
	es=GGEW(pr(i).p(:,3)+273.15);
	rh=pr(i).p(:,4)/100;
	pr(i).p(1:sidat,10)=(0.622*(es.*rh)./(pr(i).p(:,2)*100-es.*rh) )'; %mixing ratio
	
	
	
	
	fclose(fid);

end



'mwng'
%rhstart=dat(6,1); %relative humidity at ground (%)

