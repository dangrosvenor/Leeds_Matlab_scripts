t=59;   %35
dt=1;

exdir='field\cons\rainMRprofs\';


scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
ovmax=0; %overall max
ovmin=0;

for tt=1:dt:t
    

    for jj=1:nplots
        

        
        

        
        %me=mean(divs(jj).d(:,:,tt),2);
        
        %me=divp(jj).d(tt,:);
        %me=imrprof(jj).prof(:,tt)+smrprof(jj).prof(:,tt)+gmrprof(jj).prof(:,tt);
        
        me=vlrprof(jj).prof(3,:,tt);
        
        maxx(jj)=max(me)*1.1;
        minx(jj)=min(me)*1.1;
        
        
        
        
	end
    if (max(maxx)>ovmax)
        ovmax=max(maxx)
    end
    if (min(minx)<ovmin)
        ovmin=min(minx)
    end
end

for tt=1:dt:t
    hfdp=figure('Position',posit);
    set(gcf,'paperpositionmode','auto');

    for jj=1:nplots
        hs(jj).h=subplot(1,3,jj); %subplotting
         %me=mean(divs(jj).d(:,:,tt),2);
         %me=imrprof(jj).prof(:,tt)+smrprof(jj).prof(:,tt)+gmrprof(jj).prof(:,tt);
         
         me=vlrprof(jj).prof(3,:,tt);
         
         plot(me,Grid.Z);
        axis(hs(jj).h,[ovmin ovmax min(Grid.Z) max(Grid.Z)]);   
    end
   
        
        title(direcDan(jj).dir);
        
%         maxy(jj)=max(me)*1.1;
%         miny(jj)=min(me)*1.1;

    text(min(minx),-max(Grid.Z)/10,strcat('TIME=',num2str(time(tt)/3600,2),' hrs'));
    
    exname=strcat('c:\matlabR12\work\',exdir,'rainMR','-',int2str(tt+1),'.jpg');
    
	%gcf=hfdp;
	
	print(gcf,'-djpeg','-r150',exname);
	close(gcf);
    
end