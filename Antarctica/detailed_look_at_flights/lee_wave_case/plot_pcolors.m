x_start = 146;

DT_acpim=0.5; %timestep in seconds        
time_acpim=[1:size(acpim.Z,1)]*DT_acpim;   
dist_acpim = 20*time_acpim/1000 + x_start;

fsize=12;

tag='dz_30';

picname=['Dist-Height LWC ACPIM' tag];
figure('name',picname);
pcolor(dist_acpim,acpim.Z(1,:),rho'.*acpim.qc'*1e3);shading flat;colorbar
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
title('LWC (g m^{-3})','fontsize',fsize);
set(gca,'fontsize',fsize);
saveas_ps_fig_emf(gcf,[filedir picname]);

picname=['Dist-Height Ice No. ACPIM' tag];
figure('name',picname);
pcolor(dist_acpim,acpim.Z(1,:),rho'.*acpim.conci'*1e-3);shading flat;colorbar
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Ice number (L^{-1})','fontsize',fsize);
saveas_ps_fig_emf(gcf,[filedir picname]);

picname=['Dist-Height Droplet No. ACPIM per cc' tag];
figure('Name',picname);
pcolor(dist_acpim,acpim.Z(1,:),rho'.*acpim.conc2'*1e-6);shading flat;colorbar
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Droplet number (cm^{-3})','fontsize',fsize);
%saveas_ps_fig_emf(gcf,[filedir picname]);


picname=['Dist-Height Droplet No. ACPIM per mg' tag];
figure('Name',picname);
pcolor(dist_acpim,acpim.Z(1,:),acpim.conc2'*1e-6);shading flat;colorbar
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Drople number (mg^{-1})','fontsize',fsize);
%saveas_ps_fig_emf(gcf,[filedir picname])