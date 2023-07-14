t=size(time,2);   %35
dt=1;

exdir='field\cons\divergence\profilesSAME\';

clear dmax

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
ovmax=0; %overall max
ovmin=0;

for tt=1:dt:t
    

    for jj=1:nplots
        

        
        

        
        dmax(jj,tt)=max(max(divs(jj).d(1:15,:,tt)));
        dmin(jj,tt)=min(min(divs(jj).d(1:15,:,tt)));
        %me=divp(jj).d(tt,:);
        
        
        
        
        
        
	end
%     if (max(maxx)>ovmax)
%         ovmax=max(maxx)
%     end
%     if (min(minx)<ovmin)
%         ovmin=min(minx)
%     end
end


    hfdp=figure('Position',posit);
    set(gcf,'paperpositionmode','auto');

    
        maxx=max(max(dmax))*1.1;
        minx=min(min(dmax))*1.1;
        maxx2=max(max(dmin))*1.1;
        minx2=min(min(dmin))*1.1;
       
    
    for jj=1:nplots
        
        hs(jj).h=subplot(3,1,jj); %subplotting
         %me=mean(divs(jj).d(:,:,tt),2);
         plot(time./3600,dmax(jj,:));
         hold on;
         plot(time./3600,dmax(jj,:),'kx');
         axis(hs(jj).h,[0 max(time)/3600 minx maxx]); 
         title(direcDan(jj).dir);
    end
   
        
        
        
%         maxy(jj)=max(me)*1.1;
%         miny(jj)=min(me)*1.1;

hfdp=figure('Position',posit);
    set(gcf,'paperpositionmode','auto');

    
       
    
    for jj=1:nplots
        
        hs(jj).h=subplot(3,1,jj); %subplotting
         %me=mean(divs(jj).d(:,:,tt),2);
         plot(time./3600,dmin(jj,:));
         hold on;
         plot(time./3600,dmin(jj,:),'kx');
         axis(hs(jj).h,[0 max(time)/3600 minx2 maxx2]);
         title(direcDan(jj).dir);
    end
    

%     text(min(minx),-max(Grid.Z)/10,strcat('TIME=',num2str(time(tt)/3600,2),' hrs'));
%     
%     exname=strcat('c:\matlabR12\work\',exdir,'divp','-',int2str(tt+1),'.jpg');
%     
% 	%gcf=hfdp;
% 	
% 	print(gcf,'-djpeg','-r150',exname);
% 	close(gcf);
    
   