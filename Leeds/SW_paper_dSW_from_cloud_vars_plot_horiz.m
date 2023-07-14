ylabs={'\DeltaF_{SW\uparrow} (W m^{-2})'};
period_lab = period_labs{itr};
msize = 14;
fsize_text = 18;
end_cap_fraction_size = 0.0005; %size of end caps for error bars 
% - note, the ylim scale is quite large before the tidy routine, so need to make this quite small

bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]};

clear ytick_labs




iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by



%% Contributions

ngaps = 5; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;


i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = eval(['table_vals.altClearSkycont' period_str_char ';']); uncer = eval(['table_vals.altClearSkycont' period_str_char 'un;']);
plot(val,x,'o','color',bar_cols{1},'markerfacecolor',bar_cols{1},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{1},'o',3,end_cap_fraction_size);
ytick_labs{i}='F_{SW\uparrow}^{clear-sky}';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = eval(['table_vals.altLWPcont' period_str_char ';']); uncer = eval(['table_vals.altLWPcont' period_str_char 'un;']);
plot(val,x,'o','color',bar_cols{1},'markerfacecolor',bar_cols{1},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{1},'o',3,end_cap_fraction_size);
ytick_labs{i}='L';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = eval(['table_vals.altCFcont' period_str_char ';']); uncer=eval(['table_vals.altCFcont' period_str_char 'un;']);
plot(val,x,'o','color',bar_cols{1},'markerfacecolor',bar_cols{1},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{1},'o',3,end_cap_fraction_size);
ytick_labs{i}='f_{c}';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = eval(['table_vals.altNdcont' period_str_char ';']); uncer=eval(['table_vals.altNdcont' period_str_char 'un;']);
plot(val,x,'o','color',bar_cols{1},'markerfacecolor',bar_cols{1},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{1},'o',3,end_cap_fraction_size);
ytick_labs{i}='N_d';


%% Tidy plot
%Hawaii_ERF_bar_plot_horiz_TIDY
ivar=1;
SW_paper_dSW_from_cloud_vars_plot_horiz_TIDY


savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_' DAMIP_runs{1} '_Region4_' var_str_tab ' ' land_ocean '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_jpg=0;
opts.iplot_eps=1;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% END FIG