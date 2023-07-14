%imports several eta files and makes flux curves for Bauru

outdir='c:\documents and settings\user\my documents\HIBISCUS\baurufield\mesoeta\bin\05.02.04\prec\';
patt='c:\documents and settings\user\my documents\HIBISCUS\baurufield\mesoeta\bin\05.02.04\';
%etabru10.2004020512';

xx=90;   %62;
yy=5;    %37;


lis=dir(patt);
[aeta beta]=size(lis);

for ieta=4:aeta
    pat=strcat(patt,lis(ieta).name);
    [siza sizb]=size(lis(ieta).name);
    if sizb==26
        impetashort;
   
        sen(ieta)=-xdan(5).x(yy,xx); %point where Bauru is
        lat(ieta)=-xdan(10).x(yy,xx);
        gro(ieta)=-xdan(9).x(yy,xx);
        
        senmax(ieta)=max(max(-xdan(5).x(:,1:95) )); 
        latmax(ieta)=max(max(-xdan(10).x));
        gromax(ieta)=max(max(-xdan(9).x));
        
       
        
	%     figure;
	%     ba=-xdan(5).x(27:47,52:72);
	%     pcolor(ba);shading interp;
	%     caxis([0 450]);
	%     colorbar;
	%     hold on;plot(11,11,'kx');
	%     
	%     set(gcf,'paperpositionmode','auto');
	%     exname=strcat(outdir,'SEN-',lis(ieta).name(1:end-4),'.jpg');
	%     print(gcf,'-djpeg','-r150',exname);
	%     close(gcf);
	%     
	%     figure;
	%     ba=-xdan(10).x(27:47,52:72);
	%     pcolor(ba);shading interp;
	%     caxis([0 660]);
	%     colorbar;
	%     hold on;plot(11,11,'kx');
        
	
        figure;
        ba=xdan(14).x(27:47,52:72);
        pcolor(ba);shading interp;
        %caxis([0 660]);
        colorbar;
        hold on;plot(11,11,'kx');
        
	
	
        set(gcf,'paperpositionmode','auto');
        exname=strcat(outdir,'PW-',lis(ieta).name(1:end-4),'.jpg');
        print(gcf,'-djpeg','-r150',exname);
        close(gcf);
        
    end %if sizb==26
    
    
    
    
end