figure
%ecmwf_dat = nc{'al'}(:);
%ecmwf_dat = nc{'istl1'}(:);
ecmwf_dat = nc{'stl1'}(:)-273.15;
%ecmwf_dat = nc{'sd'}(:);
%ecmwf_dat = nc{'asn'}(:);
%ecmwf_dat = nc{'skt'}(:)-273.15;

pcolor(nc{'longitude'}(:)-360,nc{'latitude'}(:),squeeze(ecmwf_dat));shading flat;
set(gca,'clim',[0 5]);
colorbar

hold on

hline=line([285-360;315-360],[-67;-67]); set(hline,'color','w');
hline=line([285-360;315-360],[-68;-68]); set(hline,'color','w');
hline=line([-60;-60],[-75;-60]); set(hline,'color','w');
hline=line([-65;-65],[-75;-60]); set(hline,'color','w');

%set(gca,'xlim',[-65 -58]);
%set(gca,'ylim',[-70 -65]);