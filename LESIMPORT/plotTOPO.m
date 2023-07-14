%draws topography map around Brazil area
%x is 1200x1200 and each cell is approx 1km

diff1=100; %number of grid spaces either side of Bauru - each ~1km - max 283
diff2=100; %for x direction - max 115


outdir='c:\documents and settings\user\my documents\HIBISCUS\baurufield\mesoeta\topo\';

b1=917; %co-ordinates for Bauru -  Y
b2=116; % - X

figure;
pcolor(x(b1-diff1:b1+diff1,b2-diff2:b2+diff2));
%pcolor(x(1:end,1:end)); %for some reason wouldn't accept pcolor(x) plot as a valid figure handle - had to do x(1:end,1:end)
shading interp;
colorbar;
hold on;
%plot(diff2+1,diff1+1,'kx');
plot(b2,b1,'kx');

set(gcf,'paperpositionmode','auto');
    exname=strcat(outdir,'pic-',num2str(diff1),'-',num2str(diff2),'.jpg');
    %exname=strcat(outdir,'pic-whole.jpg');
    print(gcf,'-djpeg','-r150',exname);
    close(gcf);