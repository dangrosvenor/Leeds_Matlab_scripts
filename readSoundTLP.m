%read in .tlp style sounding (from Brazilian radiosonde network)

clear dire temp dat pro
pr(1).p=[];

readsound;
datB=dat;
clear dat;

rhbauru=0;   %flag for adding rh from Bauru sounding where =100% RH
windbauru=0; %flag for adding wind from bauru sounding where missing
isave=0; %flag to output sounding or not

% dire(1).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_01_30_12AED.txt';
% dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_01_31_12AED.txt';
% dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
% dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';

%dire(1).d='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\5\04_02_24_2015_dmi.AED.txt';
%dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_06_18AED.txt';
%dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
%dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';

dire(1).d='c:\documents and settings\Login\my documents\HIBISCUS\troccibras\radiosondesnetwork\040224_12.tlp';
dire(2).d='c:\documents and settings\Login\my documents\HIBISCUS\troccibras\radiosondesnetwork\040225_00.tlp';

% dire(1).d='c:\documents and settings\Login\my documents\HIBISCUS\troccibras\radiosondesnetwork\040213_12.tlp';
% dire(2).d='c:\documents and settings\Login\my documents\HIBISCUS\troccibras\radiosondesnetwork\040214_00.tlp';

savedir='c:\documents and settings\Login\my documents\HIBISCUS\troccibras\radiosondesnetwork\CampoGrande_040213_12';

station=[83554 83612 83779 83840 83827 83899]; %Corumba, CG, SP, Curitiba, Foz do Igacu, Florianopolis

%station=83554;

for i=1:size(dire,2)





for jj=1:length(station)
    
    fid=fopen(dire(i).d,'rb'); %re-open file so can scan from beginning again

cont=1;
stop=0;
temp='ok';
finish=0;
while stop==0
    clear dat
    clear dat2
	stat=fscanf(fid,'%g',[1 1]);
    
    if length(stat)==0; %if is a string read in temp string
       temp=fscanf(fid,'%s',[1 1]);
    end
	if stat==station(jj); stop=1;end %stop if correct station found
    if length(temp)==0; 
        stop=1;
        finish=1; %station not found
    end
end

cont=1;
while (cont==1 & finish~=1) 
temp=fscanf(fid,'%s',[1 1]);
test2=abs(str2num(temp)-1111); %1111 is marker for start of sounding data

    if length(test2)==0 
        cont=1;
	else
        if test2(1)>0.01;cont=1;
        else cont=0; end
	end

end



ii=1;
cont=1;
cont2=1;
while (cont==1 & finish~=1) 
	temp=fscanf(fid,'%s',[1 1]);
	test2=abs(str2num(temp)-2222); %2222 is marker for end of sounding data
	if abs(str2num(temp)-9999)<0.01; %some don't have a 2222, just a 9999
        test2(1)=0;
        cont2=0;
	end

        if test2(1)>0.01
            cont=1;
            dat(1,ii)=str2num(temp);
            dat(2:4,ii)=fscanf(fid,'%g',[3 1]);
            ii=ii+1;
        else 
            cont=0; 
        end
        
end

ii=1;
cont=1;
while (cont==1 & cont2==1 & finish~=1) 
temp=fscanf(fid,'%s',[1 1]);
test2=abs(str2num(temp)-9999); %9999 is marker for end of wind sounding data
	
        if test2(1)>0.01;
            cont=1;
            dat2(1,ii)=str2num(temp);
            dat2(2:3,ii)=fscanf(fid,'%g',[2 1]);
            ii=ii+1;
        else 
            cont=0; 
        end
        
end

%dat=fscanf(fid,'%g',[4 inf]);
%sidat=(size(dat,2));


%end %main while for station number

if finish~=1

	%fill in wind data with zeroes
	%zrows=find( dat(1,:)>dat2(1,1) | dat(1,:)<dat2(1,end) );
	dat(9,:)=0;
	dat(8,:)=0;
	
	%add wind in dat2 to dat
	for iw=1:size(dat2,2)
		datrows(iw)=find( dat(1,:)==dat2(1,iw) );
        dat(8,datrows(iw))=dat2(2,iw);
        dat(9,datrows(iw))=dat2(3,iw)*0.514791; 
	end
	
        
	% datrows=find( dat(1,:)<=dat2(1,1) & dat(1,:)>=dat2(1,end) );
	% dat(8,datrows)=dat2(2,:);
	% dat(9,datrows)=dat2(3,:)*0.514791;    
	
	%work out RH from dew point temp
	ev=GGEW(dat(4,:)+273.15); %vap pressure
	es=GGEW(dat(3,:)+273.15); %sat pressure
	dat(6,:)=100*ev./es; %rh

end

if rhbauru==1
	sat=find(abs(ev-es)<0.01); %make sure that any points that are saturated are not to be kept
	for j=1:length(sat)
        ii=sat(j);
        [a b]=min( abs(dat(1,ii)-datB(4,:)) ); %find nearest pressure level to that missing
        dat(6,ii)=datB(6,b); %rh
	end
end

%ei=SatVapPress(dat(3,:)+273.15,'buck2','ice');
%icesat=find(ev>=ei);
%dat(6,icesat)=100*ei(icesat)./es(icesat);

if windbauru==1
	%add wind fields from other sounding for missing data
	for j=1:length(zrows)
        ii=zrows(j);
        [a b]=min( abs(dat(1,ii)-datB(4,:)) ); %find nearest pressure level to that missing
        dat(8,ii)=datB(8,b); %dir
        dat(9,ii)=datB(9,b); %mag
	end
end

if finish~=1
	clear pr;
	sidat=size(dat,2);
	pr(i).p(1:sidat,3)=dat(3,:)'; %temp
	pr(i).p(1:sidat,2)=dat(1,:)'; %pressure
	pr(i).p(1:sidat,4)=dat(6,:)'; %RH
	pr(i).p(1:sidat,5)=dat(9,:)'.*cos( dat(8,:)'.*pi./180 ); %vwind from N I think
	pr(i).p(1:sidat,6)=dat(9,:)'.*sin( dat(8,:)'.*pi./180 );%uwind from E
	pr(i).p(1:sidat,1)=dat(2,:)'; %height
	
	pr(i).p(1:sidat,7)=log10(pr(i).p(1:sidat,3).*(1000./pr(i).p(1:sidat,2)).^0.286);%log theta
	pr(i).p(1:sidat,8)=dat(9,1:sidat)'; %wind mag
	pr(i).p(1:sidat,9)=dat(8,1:sidat)'; %wind dir
	
	es=GGEW(dat(3,:)+273.15);
	rh=dat(6,:)/100;
	pr(i).p(1:sidat,10)=(0.622*(es.*rh)./(dat(1,:)*100-es.*rh) )';
	
	if i>1
        si1=size(pr(i).p,1);
        si2=size(pro(jj).p,1);
        if si1>si2
            
            for i2=1:size(pro(jj).p,2)
                pro(jj).p(si2:si1,i2,i-1)=pro(jj).p(si2,i2,i-1);
            end
	
        end
	end
	
	pro(jj).p(1:size(pr(i).p,1),1:size(pr(i).p,2),i)=pr(i).p;
end
    
fclose(fid);
end %end for repeating over more than one station (jj)


end %for for dirs (i)

if isave==1
	fid=fopen(savedir,'wt');
	for i=1:size(dat,2)
        fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',999,999,dat(2,i),dat(1,i),dat(3,i),dat(6,i),dat(4,i),dat(8,i),dat(9,i));
	end
end

lpro=length(pro);
readsound;
for i=1:length(pr)
    pro(lpro+1).p(1:size(pr(i).p,1),1:size(pr(i).p,2),i)=pr(i).p;
end








