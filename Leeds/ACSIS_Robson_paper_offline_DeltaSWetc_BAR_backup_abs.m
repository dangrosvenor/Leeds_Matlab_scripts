
%% Bar plot of delta SW (rather than trend contribution - i.e., taking into account the time over which the trend has occurred for).
%% Could consider doing this analysis for whole downward SW trend part after 1971 rather than just post 1985.
%% Could have a first plot with just the ensemble mean on and no obs, AMIP, etc. showing the longer trend fits for period 2.

%xlabs_str={'Period 1','Period 2'};

%full_calc = 4.64;
dT_P1 = (yr_end_trend_box2(1) - yr_start_trend_box2(1)) * 1e-2; %divided by 100 to account for the scale factor used for the trends
dT_P2 = (yr_end_trend_box2(2) - yr_start_trend_box2(2)) * 1e-2;

%Period 1
full_calc_PA = table_vals.FullCalcPA * dT_P1;
%full_calc_uncer = sw_trends.trend_dat_box_obs2{ibox,it_trend}.uncer_max;

clear bar_dat xlabs_str

ibar=1;
xlabs_str{ibar}='';
bar_dat(ibar,1) = full_calc_PA - table_vals.NdconstPA * dT_P1;
bar_dat(ibar,2) = full_calc_PA - table_vals.CFconstPA * dT_P1;
bar_dat(ibar,3) = full_calc_PA - table_vals.LWPconstPA * dT_P1;
bar_dat(ibar,4) = 0; %to represent the calculated trend
bar_dat(ibar,5) = 0; %to represent the actual model trend

%Next bar just for calcaulted trend
ibar=ibar+1;
xlabs_str{ibar}='Period 1';
bar_dat(ibar,1) = 0;
bar_dat(ibar,2) = 0;
bar_dat(ibar,3) = 0;
bar_dat(ibar,4) = full_calc_PA;
bar_dat(ibar,5) = 0;

%Next bar just for model (true) trend
ibar=ibar+1;
xlabs_str{ibar}='';
bar_dat(ibar,1) = 0;
bar_dat(ibar,2) = 0;
bar_dat(ibar,3) = 0;
bar_dat(ibar,4) = 0;
bar_dat(ibar,5) = table_vals.FullModelPA * dT_P1;


% -- Period 2

ibar=ibar+1;
full_calc_PB = table_vals.FullCalcPB * dT_P2;
bar_dat(ibar,1) = (full_calc_PB - table_vals.NdconstPB * dT_P2);
bar_dat(ibar,2) = (full_calc_PB - table_vals.CFconstPB * dT_P2);
bar_dat(ibar,3) = (full_calc_PB - table_vals.LWPconstPB * dT_P2);
bar_dat(ibar,4) = 0; %to represent the calculated trend
bar_dat(ibar,5) = 0; %to represent the actual model trend

%Next bar just for calcalated trend
ibar=ibar+1;
xlabs_str{ibar}='Period 2';
bar_dat(ibar,1) = 0;
bar_dat(ibar,2) = 0;
bar_dat(ibar,3) = 0;
bar_dat(ibar,4) = full_calc_PB;
bar_dat(ibar,5) = 0;

%Next bar just for model (true) trend
ibar=ibar+1;
xlabs_str{ibar}='';
bar_dat(ibar,1) = 0;
bar_dat(ibar,2) = 0;
bar_dat(ibar,3) = 0;
bar_dat(ibar,4) = 0;
bar_dat(ibar,5) = table_vals.FullModelPB * dT_P2;



figure
set(gcf,'color','w');
bar(bar_dat,'stacked'); %will give N bars each with M components for bar_dat[N,M]
increase_font_size_map_figures
%title(['SW_{in}=' num2str(SW_in) ', A_{clear}=' num2str(A_clear)]);
ylabel('\Delta SW_{up TOA} (W m^{-2})');
set(gca,'xticklabel',xlabs_str);
leg_str = {'N_{d}','f_{c}','LWP_{ic}','Full Calculated','Full Actual'};
pad_list=[1 2];
for i=1:length(pad_list)
    iL=pad_list(i);
    leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
end
pad_list=[3 4 5];
for i=1:length(pad_list)
    iL=pad_list(i);
    leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
end
L=legend(leg_str,'location','NorthEast','fontsize',14);
%set(L,'string',leg_str)


set(gca,'xlim',[0 ibar+1]);
set(gca,'ylim',[-6.5 10.5]);


savename=[savedir_date titlenam ' ' ens_str ' ' land_ocean '_BAR_PLOT_dSW'];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

