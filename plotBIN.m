%plots divergence binned
clear bav xdat;

nplots=input('Enter number of plots: ');
nb=input('Enter number bins: ');
sb=(size(Grid.Y1)-2)/nb; %size of bin for divergences

tend=30;



xdat=Grid.Y1(2:sb:sb*nb+1)./1000;
figure;
    for jj=1:nplots
        
%         dav(j).d(1:sb*nb-2)=dav(j).d(2:sb*nb-1);
        
        for ih=1:nb
            tmean=mean(divs(jj).d(1:15,(ih-1)*sb+1:ih*sb,1:tend),3);
            htmean=mean(tmean,1);
            bav(jj).b(ih)=mean(htmean); %average over bin
        end

        hs(jj).h=subplot(3,1,jj); %subplotting
        plot(xdat,bav(jj).b);
        hold on;
        plot(xdat,bav(jj).b,'kx')
        title(direcDan(jj).dir);
        
        maxy(jj)=max(bav(jj).b)*1.1;
        miny(jj)=min(bav(jj).b)*1.1;
        
        
	end
    
    for jj=1:nplots
        axis(hs(jj).h,[min(xdat) max(xdat) min(miny) max(maxy)]);   
    end