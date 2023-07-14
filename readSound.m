%read in AED.txt style sounding


clear dire



% dire(1).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_01_30_12AED.txt';
% dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_01_31_12AED.txt';
% dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
% dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';


%dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_06_18AED.txt';
%dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
%dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';

dire(1).d='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\5\04_02_24_2015_dmi.AED.txt';
%dire(2).d='c:\documents and settings\Login\my documents\HIBISCUS\TrocciBras\radiosondesNetwork\CampoGrande_040224_12';
%dire(2).d='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\5\04_02_26_12AED.txt';

% dire(1).d='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\3\04_02_13_12AED.txt';
% dire(2).d='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\3\04_02_13_18AED.txt';

for i=1:size(dire,2)
pr(i).p=[]; %clear arrays
    
fid=fopen(dire(i).d,'rb');

dat=fscanf(fid,'%g',[9 inf]);
sidat=(size(dat,2));

pr(i).p(1:sidat,3)=dat(5,:)'; %temp
pr(i).p(1:sidat,2)=dat(4,:)'; %pressure
pr(i).p(1:sidat,4)=dat(6,:)'; %RH
pr(i).p(1:sidat,5)=dat(9,:)'.*cos( dat(8,:)'.*pi./180 ); %vwind from N I think
pr(i).p(1:sidat,6)=dat(9,:)'.*sin( dat(8,:)'.*pi./180 );%uwind from E
pr(i).p(1:sidat,1)=dat(3,:)'; %height

pr(i).p(1:sidat,7)=log10(pr(i).p(1:sidat,3).*(1000./pr(i).p(1:sidat,2)).^0.286);%log theta
pr(i).p(1:sidat,8)=dat(9,1:sidat)'; %wind mag
pr(i).p(1:sidat,9)=dat(8,1:sidat)'; %wind dir

es=GGEW(dat(5,:)+273.15);
rh=dat(6,:)/100;
pr(i).p(1:sidat,10)=(0.622*(es.*rh)./(dat(4,:)*100-es.*rh) )';




fclose(fid);

end



'finito'
%rhstart=dat(6,1); %relative humidity at ground (%)

