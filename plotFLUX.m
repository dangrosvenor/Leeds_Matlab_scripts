exdir='field\onesine\';

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

h=figure('Position',posit);

% plot(Grid.Y1(2:end-1)./1000,lflux(3).l);
% xlabel('X (km)');
% ylabel('Latent Heat Flux (Wm^-2)');
% title('Surface variation of Lat. Heat Flux at 13:00');
% 
% 
% gcf=h;
% set(gcf,'paperpositionmode','auto');
% exname=strcat('c:\matlabR12\work\',exdir,'latFlux.jpg');	
% print(gcf,'-djpeg','-r150',exname);


% plot(Grid.Y1(2:end-1)./1000,sflux(3).s);
% xlabel('X (km)');
% ylabel('Sensible Heat Flux (Wm^-2)');
% title('Surface variation of Sens. Heat Flux at 13:00');
% 
% 
% gcf=h;
% set(gcf,'paperpositionmode','auto');
% exname=strcat('c:\matlabR12\work\',exdir,'sensFlux.jpg');	
% print(gcf,'-djpeg','-r150',exname);

plot(SerDan(1).SER(1:end-1000,1)./3600,SerDan(1).SER(1:end-1000,2));
hold on;
%plot(SerDan(1).SER(:,1)./3600,SerDan(1).SER(:,3),'kx');
xlabel('Time (hrs)');
ylabel('Average Latent Heat Flux (Wm^-2)');
title('Time variation of Av. Lat. Heat Flux');


gcf=h;
set(gcf,'paperpositionmode','auto');
exname=strcat('c:\matlabR12\work\',exdir,'latsFluxTvar.jpg');	
print(gcf,'-djpeg','-r150',exname);