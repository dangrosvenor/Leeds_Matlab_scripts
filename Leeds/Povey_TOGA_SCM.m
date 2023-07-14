%TOGA-COARE SCM quick plots

output_dir = '/home/disk/eos15/d.grosvenor/UM/Povey/';

file_01 = 'u-ci899_TOGAhydP_L85_ga8.nc';
file_02 = 'u-ci900_TOGAhydP_L85_ga8.nc';

nc_01 = netcdf([output_dir file_01]);
nc_02 = netcdf([output_dir file_02]);

z=[1:85];

%%
%N.B. - cc_qcf_rad values are all zero
qcf_01 = nc_01{'qcf'}(:);
qcf_02 = nc_02{'qcf'}(:);

mean_01 = meanNoNan(qcf_01,1);
mean_02 = meanNoNan(qcf_02,1);

dq = mean_01 - mean_02;
dq_pct = 100*(mean_01 - mean_02 ) ./ mean_02;

figure('color','w');
plot(dq_pct,z,'bo-','linewidth',3);
xlabel('\Deltaq_{ice} (%)');    
ylabel(['Model Level']);
fontsize_figure(gcf,gca,18);
title('100*(u-ci899 minus u-ci900)/u-ci900');
set(gca,'xlim',[-20 20]);
set(gca,'ylim',[1 60]);
grid on

figure('color','w');
plot(dq,z,'bo-','linewidth',3);
xlabel('\Deltaq_{ice} (kg kg^{-1})');    
ylabel(['Model Level']);
fontsize_figure(gcf,gca,18);
title('u-ci899 minus u-ci900')
%set(gca,'xlim',[-20 20]);
set(gca,'ylim',[1 60]);
grid on

clear leg_str
figure('color','w');
plot(mean_01,z,'bo-','linewidth',3); leg_str{1}='u-ci899';
hold on
plot(mean_02,z,'ro-','linewidth',3); 
xlabel('q_{ice} (kg kg^{-1})'); leg_str{2}='u-ci900';    
ylabel(['Model Level']);
fontsize_figure(gcf,gca,18);
%set(gca,'xlim',[-20 20]);
set(gca,'ylim',[1 60]);
grid on
legend(leg_str);



%%
%N.B. - cc_qcl_rad values are all zero
qcl_01 = nc_01{'qcl'}(:);
qcl_02 = nc_02{'qcl'}(:);

mean_01 = meanNoNan(qcl_01,1);
mean_02 = meanNoNan(qcl_02,1);

dq = mean_01 - mean_02;
dq_pct = 100*(mean_01 - mean_02 ) ./ mean_02;

figure('color','w');
plot(dq_pct,z,'bo-','linewidth',3);
xlabel('\Deltaq_{liq} (%)');    
ylabel(['Model Level']);
fontsize_figure(gcf,gca,18);
title('100*(u-ci899 minus u-ci900)/u-ci900');
%set(gca,'xlim',[-20 20]);
set(gca,'ylim',[1 60]);
grid on

figure('color','w');
plot(dq,z,'bo-','linewidth',3);
xlabel('\Deltaq_{liq} (kg kg^{-1})');    
ylabel(['Model Level']);
fontsize_figure(gcf,gca,18);
title('u-ci899 minus u-ci900')
%set(gca,'xlim',[-20 20]);
set(gca,'ylim',[1 60]);
grid on

clear leg_str
figure('color','w');
plot(mean_01,z,'bo-','linewidth',3); leg_str{1}='u-ci899';
hold on
plot(mean_02,z,'ro-','linewidth',3); 
xlabel('q_{liq} (kg kg^{-1})'); leg_str{2}='u-ci900';    
ylabel(['Model Level']);
fontsize_figure(gcf,gca,18);
%set(gca,'xlim',[-20 20]);
set(gca,'ylim',[1 60]);
grid on
legend(leg_str);








%%
nc_01 = close(nc_01);
nc_02 = close(nc_02);