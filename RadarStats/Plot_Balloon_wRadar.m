%readALLsss; %read in ozone sss sounde data - puts into data(i) where i=1:4 for SF numbers
%time =data(i).sss(1,:) in secs after 20:03 UT = 17:03 LT. Lat in degs in col 5 lon in col 6. Lat/Lons are positive so are W & S.

dirout='c:/documents and settings/login/my documents/hibiscus/troccibras/radar/bauru/echotops/24.02/balloonTraj/';



xc=253; %centre of radar image in im (origin at top left)
yc=305;

y400=504-106; %400km circumference of radar pic in x and y dirs
x400=453-53;

Blat=22.36;
Blon=49.03;

npoints=15;
	ntot=length(data(4).sss(1,:));
	
	i=1:ntot/npoints:ntot; %array for places to plot balloon position
	
	lats=data(4).sss(5,i);
	lons=data(4).sss(6,i);
	
	
	
	xd=sign(Blon-lons).*distlatlon(Blat,Blon,Blat,data(4).sss(6,i)); %horizontal distances from Bauru
	yd=sign(Blat-lats).*distlatlon(Blat,Blon,data(4).sss(5,i),Blon); %vertical distances from Bauru (km)
	
	x=x400*xd/400;
	y=y400*yd/400;
	
	x0=xc+x; %distance from top left
	y0=yc-y; %reason for diff signs is that x and y were postive for above centre point and to right
    
    times=17+(0.03/60) + data(4).sss(1,i)/3600; %times in LT decimal hours

dire='c:/documents and settings/login/my documents/hibiscus/troccibras/radar/bauru/echotops/24.02/';
d=dir(dire);


tstart=17+(0.03/60);
times=tstart  + data(4).sss(1,:)/3600; %times in LT decimal hours

lats=data(4).sss(5,:);
lons=data(4).sss(6,:);
		
		
		
		xd=sign(Blon-lons).*distlatlon(Blat,Blon,Blat,data(4).sss(6,:)); %horizontal distances from Bauru
		yd=sign(Blat-lats).*distlatlon(Blat,Blon,data(4).sss(5,:),Blon); %vertical distances from Bauru (km)
		
		x=x400*xd/400;
		y=y400*yd/400;
		
		x0=xc+x; %distance from top left
		y0=yc-y; %reason for diff signs is that x and y were postive for above centre point and to right

        
        
for ii=1:length(d)   %119
    if length(d(ii).name)>4
    if strcmp('.gif',d(ii).name(end-3:end))==1
    %if strcmp(d(ii).name,'ppi_bru_240204_2100h.gif')==1
        date(ii).d=(d(ii).name(11:16));
        day(ii)=str2num(date(ii).d(1:2));
        month(ii)=str2num(date(ii).d(3:4));
        year(ii)=2000+str2num(date(ii).d(5:6));
        hrs(ii)=str2num(d(ii).name(18:19));
        %hrs=mod(hrs+3,24); %add 3 hrs to convert to UTC
        mins(ii)=str2num(d(ii).name(20:21));
    end
    end
    dh=hrs(ii)+(mins(ii)/60); %radar image times in decimal hours
    
    if dh>=tstart;
    
        %read in and draw image
        fi=strcat(dire,d(ii).name);
        [im map]=imread(fi,'gif');
        
        scrsz=get(0,'ScreenSize');
        posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
        figure('position',posit);
        image(im);
        colormap(map);
        set(gca,'position',[0 0 1 1]); %set so axis doesn't appear on screen 1=full screen
        
        [hmin i]=min(abs(dh-times)); %index for closest matching time in sss flight

        h=data(4).sss(4,i)/1000; 
        h=num2str(h,'%.1f');
        h=strcat(h,'km');
        
        hold on;
        plot(x0,y0,'ko','markersize',1);
        
        plot(x0(i),y0(i),'kx','markersize',20,'linewidth',2);
        plot(x0(i),y0(i),'k+','markersize',20,'linewidth',2);
        text(x0(i)+sign(x0(i))*0.025*x0(i),y0(i),h,'fontweight','bold');
        
        exname=strcat(dirout,'BallTraj_',d(ii).name(1:end-4));
        set(gcf,'paperpositionmode','auto');
        print(gcf,'-djpeg',exname);
        
        close(gcf);
        
    end
    
    
end

% dh=hrs+(mins/60); %radar image times in decimal hours
% 
% for j=1:1; %length(times)
% 	[dhmin imin]=min(abs(dh-times(j)));
%     imin=imin+2;
%     fi=strcat(dire,lis(imin).name);
%     [im map]=imread(fi,'gif');
%     
%     scrsz=get(0,'ScreenSize');
%     posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
%     figure('position',posit);
%     image(im);colormap(map);
%     set(gca,'position',[0 0 1 1]); %set so axis doesn't appear on screen 1=full screen
%     h=data(4).sss(4,i(j))/1000; 
%     h=num2str(h,1);
%     
%     hold on;
%     plot(x(j),y(j),'kx');
%     text(x(j),y(j)-sign(y(j)),h);
%     
%     exname=strcat(dirout,'BallTraj_',num2str(hrs(j)),num2str(mins(j)));
%     set(gcf,'paperpositionmode','auto');
%     print(gcf,'-djpeg',exname);
%         
% %    fi=strcat(dire,lis(3).name);
% %    [im map]=imread(fi,'gif');
% 
% 
% end

	
























