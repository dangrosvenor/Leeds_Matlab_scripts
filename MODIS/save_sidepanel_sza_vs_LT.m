%See save_vars_for_lowSZA_allSZA_comparison.m used to save the variables.
%and load_vars_for_lowSZA_allSZA_comparison.m used to load the variables.
%for saving the variables
%create sidebar averages for the seasonal plots of %diff of <65 SZA data vs
%all SZA


man_choose_water_graph=1;
graph=145;

waterVapourMay2005  %case 145

set(gca,'position',[0.3280    0.1100    0.170    0.68141]);
set(gca,'ylim',[0 25]);
set(gca,'xlim',[45 90]);
%set(gca,'ytick',yticks_miller);
set(gca,'yticklabels','');
hylab=get(gca,'ylabel');
set(hylab,'string','');
%saveas_ps_fig_emf(gcf,[savename],'',0,0);




