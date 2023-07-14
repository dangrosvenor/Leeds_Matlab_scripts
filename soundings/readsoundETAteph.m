%read in AED.txt style sounding

clear dire
pr(1).p=[];

% dire(1).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_01_30_12AED.txt';
% dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_01_31_12AED.txt';
% dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
% dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';

fid=fopen('c:/documents and settings/user/my documents/hibiscus/baurufield/soundings/ascii/1/30.01-18etamod.txt','rb');

%dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_06_18AED.txt';
%dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
%dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';


%for i=1:size(dire,2)

%fid=fopen(dire(i).d,'rb');

a=fscanf(fid,'%g',[1 2]);

dat=fscanf(fid,'%g',[6 inf]);
sidat=(size(dat,2));

dat=dat';

%pr(i).p(1:sidat,3)=dat(3,:)'; %temp
%pr(i).p(1:sidat,2)=dat(2,:)'; %pressure
%pr(i).p(1:sidat,4)=dat(4,:)'; %RH
%pr(i).p(1:sidat,5)=dat(9,:)'.*cos( dat(8,:)'.*pi./180 ); %uwind
%pr(i).p(1:sidat,6)=dat(9,:)'.*sin( dat(8,:)'.*pi./180 );%vwind
%pr(i).p(1:sidat,1)=dat(3,:)'; %height

%pr(i).p(1:sidat,7)=log10(pr(i).p(1:sidat,3).*(1000./pr(i).p(1:sidat,2)).^0.286);%log theta
%pr(i).p(1:sidat,8)=dat(9,1:sidat)'; %wind mag



fclose(fid);

fla=4;
admulti;

%rhstart=dat(6,1); %relative humidity at ground (%)

