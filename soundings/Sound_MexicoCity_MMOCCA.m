%read in AED.txt style sounding
clear dire

dire(1).d='c:\documents and settings\Login\my documents\Leeds_MMOCCA\Milecityo.txt';


for i=1:size(dire,2)
pr(i).p=[]; %clear arrays
    
fid=fopen(dire(i).d,'rb');

dat=fscanf(fid,'%g',[5 inf]);
sidat=(size(dat,2));

pr(i).p(1:sidat,3)=dat(3,:)'; %temp (K)
%pr(i).p(1:sidat,2)=dat(4,:)'; %pressure
pr(i).p(1:sidat,4)=dat(1,:)'; %RH
%pr(i).p(1:sidat,5)=dat(9,:)'.*cos( dat(8,:)'.*pi./180 ); %vwind from N I think
%pr(i).p(1:sidat,6)=dat(9,:)'.*sin( dat(8,:)'.*pi./180 );%uwind from E
pr(i).p(1:sidat,5)=0; %vwind from N I think
pr(i).p(1:sidat,6)=0;%uwind from E

pr(i).p(1:sidat,1)=dat(4,:)'; %height - add ground height

pr(i).p(1:sidat,7)=log10(pr(i).p(1:sidat,3).*(1000./pr(i).p(1:sidat,2)).^0.286);%log theta
%pr(i).p(1:sidat,8)=dat(9,1:sidat)'; %wind mag
%pr(i).p(1:sidat,9)=dat(8,1:sidat)'; %wind dir

Y0=1013.25*100; %mean sea level pressure
TSPAN=pr(i).p(:,3); %temperature range
HSPAN=pr(i).p(:,1); %height range
[h,P] = ODE45(@hydrostatic,HSPAN,Y0,[],HSPAN,TSPAN); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H

pr(i).p(1:sidat,2)=P;

es=GGEW(pr(i).p(:,3)+273.15);
rh=pr(i).p(:,4)/100;
pr(i).p(1:sidat,10)=(0.622*(es.*rh)./(pr(i).p(:,2)*100-es.*rh) )'; %mixing ratio




fclose(fid);

end



'mwng'
%rhstart=dat(6,1); %relative humidity at ground (%)

