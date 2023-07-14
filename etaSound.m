% takes eta values and makes a sounding and plots it on a tephigram
clear dat;
fla=2;

dire='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ETA-53-43-020512+03-2.txt';
fid=fopen(dire,'w');

%xs=63; %Bauru sounding site
%ys=38;

xs=53; %100km west and 50 km north of Bauru
ys=43;

ps=[1000 980 975 960 950 925 900 875 850 825 800 775 750 725 700 675 650 625 600 575 550 525 500 475 450 425 400 375 350 325 300 275 250 225 200 175 150 125 100];
ord=[5 5 3 5 6 7 5 9 8]; %order of dat column nos. for AED.txt format output 
%pot=thd(2).x(xs,ys,:).*(1000./ps)^k;


    dat(:,3)=squeeze(thd(6).x(xs,ys,:)); %height
    dat(:,5)=ps'; %pressure
    dat(:,6)=squeeze(thd(2).x(xs,ys,:))-273; %temp (C)
    dat(:,7)=squeeze(thd(3).x(xs,ys,:)); %relative humidity
    dat(:,8)=sqrt( squeeze(thd(4).x(xs,ys,:)).^2 + squeeze(thd(5).x(xs,ys,:)).^2 ); %wind speed
    dat(:,9)=180.*atan( squeeze(thd(4).x(xs,ys,:)) ./ squeeze(thd(5).x(xs,ys,:)) ) ./ pi; %dir degrees
    %written so as to be able to insert into multitep with fla=2
    
sp=xdan(16).x(xs,ys)/100; %surface pressure in hPa
pi=find(ps<sp); %find first pressure above surface pressure


fprintf(fid,'%g %g\n',sp,xdan(26).x(xs,ys)); %print surface pressure and surface height 
for i=pi(1):39
    
for j=1:size(ord,2); 
    fprintf(fid,'%g ',dat(i,ord(j)) ); 
end

fprintf(fid,'\n');

end

fclose(fid);

fid=fopen(dire,'rb');
fscanf(fid,'%g',[1 2]);
admulti;
fclose(fid);

    

