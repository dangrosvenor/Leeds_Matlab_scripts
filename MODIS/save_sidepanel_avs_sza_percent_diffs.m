%See save_vars_for_lowSZA_allSZA_comparison.m used to save the variables.
%and load_vars_for_lowSZA_allSZA_comparison.m used to load the variables.
%for saving the variables
%create sidebar averages for the seasonal plots of %diff of <65 SZA data vs
%all SZA

ioverride_time_selection=1;
days_required_for_mean = [354-45:366 1:mod(354+45,366)]; time_mean_str = 'solar DJF';
plot_global_maps;
close(gcf);
waterVapourMay2005 %case 139 I think
set(gca,'position',[0.3280    0.1100    0.170    0.68141]);
set(gca,'xlim',[0 20]);
set(gca,'ytick',yticks_miller);
set(gca,'yticklabels','');
saveas_ps_fig_emf(gcf,[savename],'',0,0);

ioverride_time_selection=1;
days_required_for_mean = [79-45:79+45]; time_mean_str = 'solar MAM';
plot_global_maps;
close(gcf);
waterVapourMay2005
set(gca,'position',[0.3280    0.1100    0.170    0.68141]);
set(gca,'xlim',[0 20]);
set(gca,'ytick',yticks_miller);
set(gca,'yticklabels','');
saveas_ps_fig_emf(gcf,[savename],'',0,0);


ioverride_time_selection=1;
days_required_for_mean = [171-46:171+45]; time_mean_str = 'solar JJA';
plot_global_maps;
close(gcf);
waterVapourMay2005
set(gca,'position',[0.3280    0.1100    0.170    0.68141]);
set(gca,'xlim',[0 20]);
set(gca,'ytick',yticks_miller);
set(gca,'yticklabels','');
saveas_ps_fig_emf(gcf,[savename],'',0,0);


ioverride_time_selection=1;
days_required_for_mean = [263-46:263+45]; time_mean_str = 'solar SON';
plot_global_maps;
close(gcf);
waterVapourMay2005
set(gca,'position',[0.3280    0.1100    0.170    0.68141]);
set(gca,'xlim',[0 20]);
set(gca,'ytick',yticks_miller);
set(gca,'yticklabels','');
saveas_ps_fig_emf(gcf,[savename],'',0,0);


