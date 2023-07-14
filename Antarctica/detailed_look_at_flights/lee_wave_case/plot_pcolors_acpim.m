x_start = 146;
x_end=156.9;
%get col_alt, etc.
set_column_numbers_for_flight_data

%calculate the distances flown by the plane
eval(['X_flt = X_flt' flight_no ';']);
eval(['Y_flt = Y_flt' flight_no ';']);
dist = [0; cumsum(sqrt( (diff(X_flt)).^2 + (diff(Y_flt)).^2 ))];

[inds_1 inds_2] = findheight_nearest(dist,x_start,x_end);
inds2=inds_1:inds_2;

X_segment=dist(inds2);
Z_segment = dat_flt(inds2,col_alt); %altitude (m)

DT_acpim=0.5; %timestep in seconds
time_acpim=[1:size(acpim.Z,1)]*DT_acpim;
dist_acpim = 20*time_acpim/1000 + x_start;
rho=density(acpim.P,acpim.T);

Zseg_acpim = interp1(X_segment,Z_segment,dist_acpim);


fsize=12;

tag='dz_30';

picname=['Dist-Height LWC ACPIM' tag];
figure('name',picname);
pcolor(dist_acpim,acpim.Z(1,:),rho'.*acpim.qc'*1e3);shading flat;colorbar
hold on
plot(dist_acpim,Zseg_acpim,'color','k','linewidth',2);
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
title('LWC (g m^{-3})','fontsize',fsize);
set(gca,'fontsize',fsize);
saveas_ps_fig_emf(gcf,[filedir picname]);

picname=['Dist-Height Ice No. ACPIM' tag];
figure('name',picname);
pcolor(dist_acpim,acpim.Z(1,:),rho'.*acpim.conci'*1e-3);shading flat;colorbar
hold on
plot(dist_acpim,Zseg_acpim,'color','k','linewidth',2);
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Ice number (L^{-1})','fontsize',fsize);
saveas_ps_fig_emf(gcf,[filedir picname]);

picname=['Dist-Height Droplet No. ACPIM per cc' tag];
figure('Name',picname);
pcolor(dist_acpim,acpim.Z(1,:),rho'.*acpim.conc2'*1e-6);shading flat;colorbar
hold on
plot(dist_acpim,Zseg_acpim,'color','k','linewidth',2);
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Droplet number (cm^{-3})','fontsize',fsize);
%saveas_ps_fig_emf(gcf,[filedir picname]);


picname=['Dist-Height Droplet No. ACPIM per mg' tag];
figure('Name',picname);
pcolor(dist_acpim,acpim.Z(1,:),acpim.conc2'*1e-6);shading flat;colorbar
hold on
plot(dist_acpim,Zseg_acpim,'color','k','linewidth',2);
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Droplet number (mg^{-1})','fontsize',fsize);
%saveas_ps_fig_emf(gcf,[filedir picname])

picname=['Dist-Height Mean Ice Diameter ACPIM' tag];
figure('Name',picname);
pcolor(dist_acpim,acpim.Z(1,:),acpim.diami');shading flat;colorbar
hold on
plot(dist_acpim,Zseg_acpim,'color','k','linewidth',2);
xlabel('Distance along flight track (km)','fontsize',fsize);
ylabel('Height (m)','fontsize',fsize);
set(gca,'fontsize',fsize);
title('Mean Ice Diameter (\mum)','fontsize',fsize);
%saveas_ps_fig_emf(gcf,[filedir picname])