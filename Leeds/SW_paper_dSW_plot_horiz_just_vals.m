% Called from ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script.m

%ylabs={'\DeltaF_{SW} (W m^{-2})'};
%period_lab = period_labs{itr};
msize = 14;
fsize_text = 18;
end_cap_fraction_size = 0.0005; %size of end caps for error bars 
% - note, the ylim scale is quite large before the tidy routine, so need to make this quite small

bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75],[0.5 0.5 0.5]};

clear ytick_labs



%% New figure - Emission contributions
iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
% N.B, - Changing the gca position in the TIDY script is the key to making the plots more
    %compact.
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break ytick_labs
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by
%fac=0.5;


% --- Data ----

ngaps = 4; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = AerChemGHG; uncer=AerChemGHG_un;
icol=4;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='GHG-only';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = AerChemAerosol; uncer = AerChemAerosol_un;
icol=11;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='Aerosol-only';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = netAerChemMIP; uncer = netAerChemMIP_un;
icol=3;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='All emissions';




%% Tidy plot
%Hawaii_ERF_bar_plot_horiz_TIDY
%ivar=1;
SW_paper_dSW_from_cloud_vars_plot_horiz_TIDY


savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_emissions_Region4_' var_str_tab ' ' land_ocean_str '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
opts.iplot_jpg = 0;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

% --- END OF FIG ---

%% New figure - Feedback terms
iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
% N.B, - Changing the gca position in the TIDY script is the key to making the plots more
    %compact.
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break ytick_labs
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by
%fac=0.5;


% --- Data ---
switch vars{ivar}
    case {'SW'}
        ngaps = 7; %No. bars + 1
    otherwise
        ngaps = 5;
end
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = AerChemAerosol; uncer = AerChemAerosol_un;
icol=11;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='Aerosol-only';

switch vars{ivar}
        case {'SW','LWPic','CF'} %Can't estimate a forcing term on Nd since Nd is used to estimate the forcing.
            
            switch vars{ivar}
                case {'SW'}
            
                    i=i+1;
                    x=x+gap{j}; xx(i)=x;
                    val = ukesm_local_indirect; uncer=0;
                    icol=9;
                    plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
                    %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
                    ytick_labs{i}='ACI Forcing';
                    
                    
                    i=i+1;
                    x=x+gap{j}; xx(i)=x;
                    val = ukesm_local_direct; uncer=0;
                    icol=9;
                    plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
                    %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
                    ytick_labs{i}='ARI Forcing';
            end
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_local_indirect + ukesm_local_direct; uncer=0;
            icol=9;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ACI+ARI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_non_local; uncer=0;
            icol=9;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = feedback_estimate_AerChemAerosol; uncer=0;
            icol=9;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback from \DeltaT';
            
%             i=i+1;
%             x=x+gap{j}; xx(i)=x;
%             val = ukesm_local_indirect + ukesm_local_direct + feedback_estimate_AerChemAerosol; uncer=0;
%             plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
%             %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
%             ytick_labs{i}='Aerosol Forcing + Feedback';
            
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
            i=i+1;
            x=x+gap{j-1}+gap{j}; xx(i)=x;
            %val = AerChemGHG + feedback_estimate_AerChemAerosol; uncer=AerChemGHG_un; %using feedback from deltaT
            val = AerChemGHG + ukesm_non_local; uncer=AerChemGHG_un; %using feedback from total change minus forcing.
            icol=1;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            %ytick_labs{i}='Aerosol Feedback from \DeltaT + GHG Feedbacks';
            ytick_labs{i}='Total (Aerosol+GHG) Feedback';
            
            % -- end of the aerosol section for UKESM/AerChemMIP                       
            
%             ngaps = 2; %No. bars + 1
%             j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
%             %Bars within section
%             i=i+1;
%             x=x+gap{j-1}+gap{j}; xx(i)=x;
%             %Feedback using estimate from dT
%             %val = AerChemGHG + feedback_estimate_AerChemAerosol + ukesm_local_indirect + ukesm_local_direct; uncer=AerChemGHG_un;            
%             %Feedback estimated from total minus forcing - not really much
%             %point in showing this since it will equal the total.
%             val = AerChemGHG + ukesm_non_local + ukesm_local_indirect + ukesm_local_direct; uncer=AerChemGHG_un;
%             icol=2;
%             plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
%             herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
%             %ytick_labs{i}='Aerosol Forcing+Feedback + GHG Feedback';
%             ytick_labs{i}='Aerosol Forcing + Total Feedback';
            
end

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = UKESMAerChemMIP; uncer=UKESMAerChemMIP_un;
icol=3;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
%ytick_labs{i}='All emissions (UKESM model subset)';
ytick_labs{i}='All emissions';
    






%% Tidy plot
%Hawaii_ERF_bar_plot_horiz_TIDY
%ivar=1;
SW_paper_dSW_from_cloud_vars_plot_horiz_TIDY_AerChemMIP


savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_feedbacks_Region4_' var_str_tab ' ' land_ocean_str '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
opts.iplot_jpg = 0;


savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% END FIG

%% New figure - Emission contributions for HADGEM and DAMIP
iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
% N.B, - Changing the gca position in the TIDY script is the key to making the plots more
    %compact.
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break ytick_labs
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by
%fac=0.5;


% --- Data ----

ngaps = 5; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = HistNat; uncer=HistNat_un;
icol=1;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='Natural aerosol-only';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = HistGhg; uncer=HistGhg_un;
icol=4;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='GHG-only';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = HistAer; uncer = HistAer_un;
icol=11;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='Aerosol-only';

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = HAD; uncer = HAD_un;
icol=3;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='All emissions';




%% Tidy plot
%Hawaii_ERF_bar_plot_horiz_TIDY
%ivar=1;
SW_paper_dSW_from_cloud_vars_plot_horiz_TIDY


savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_emissions_DAMIP_Region4_' var_str_tab ' ' land_ocean_str '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
opts.iplot_jpg = 0;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

% --- END OF FIG ---


%% New figure - Feedback terms for DAMIP/HADGEM
iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
% N.B, - Changing the gca position in the TIDY script is the key to making the plots more
    %compact.
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break ytick_labs
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by
%fac=0.5;


% --- Data ---
switch vars{ivar}
    case {'SW'}
        ngaps = 7; %No. bars + 1
    otherwise
        ngaps = 5;
end
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x = x + gap{j-1} + gap{j}; xx(i)=x;
val = HistAer; uncer = HistAer_un;
icol=11;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='Aerosol-only';

switch vars{ivar}
        case {'SW','LWPic','CF'} %Can't estimate a forcing term on Nd since Nd is used to estimate the forcing.
            
            switch vars{ivar}
                case {'SW'}
            
                    i=i+1;
                    x=x+gap{j}; xx(i)=x;
                    val = local_indirect; uncer=0;
                    icol=9;
                    plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
                    %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
                    ytick_labs{i}='ACI Forcing';
                    
                    i=i+1;
                    x=x+gap{j}; xx(i)=x;
                    val = local_direct; uncer=0;
                    icol=9;
                    plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
                    %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
                    ytick_labs{i}='ARI Forcing';
            end
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = local_indirect + local_direct; uncer=0;
            icol=9;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ACI+ARI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = non_local; uncer=0;
            icol=9;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = feedback_estimate_aer; uncer=0;
            icol=9;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback from \DeltaT';
            
%             i=i+1;
%             x=x+gap{j}; xx(i)=x;
%             val = ukesm_local_indirect + ukesm_local_direct + feedback_estimate_AerChemAerosol; uncer=0;
%             plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
%             %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
%             ytick_labs{i}='Aerosol Forcing + Feedback';
            
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
            i=i+1;
            x=x+gap{j-1}+gap{j}; xx(i)=x;
            %val = AerChemGHG + feedback_estimate_AerChemAerosol; uncer=AerChemGHG_un; %using feedback from deltaT
            val = HistGhg + non_local; uncer=HistGhg_un; %using feedback from total change minus forcing.
            icol=1;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            %ytick_labs{i}='Aerosol Feedback from \DeltaT + GHG Feedbacks';
            ytick_labs{i}='Total (Aerosol+GHG) Feedback';
            
            % -- end of the aerosol section for UKESM/AerChemMIP                       
            
%             ngaps = 2; %No. bars + 1
%             j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
%             %Bars within section
%             i=i+1;
%             x=x+gap{j-1}+gap{j}; xx(i)=x;
%             %Feedback using estimate from dT
%             %val = AerChemGHG + feedback_estimate_AerChemAerosol + ukesm_local_indirect + ukesm_local_direct; uncer=AerChemGHG_un;            
%             %Feedback estimated from total minus forcing - not really much
%             %point in showing this since it will equal the total.
%             val = AerChemGHG + ukesm_non_local + ukesm_local_indirect + ukesm_local_direct; uncer=AerChemGHG_un;
%             icol=2;
%             plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
%             herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
%             %ytick_labs{i}='Aerosol Forcing+Feedback + GHG Feedback';
%             ytick_labs{i}='Aerosol Forcing + Total Feedback';
            
end

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = HAD; uncer=HAD_un;
icol=3;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
%ytick_labs{i}='All emissions (HADGEM)';
ytick_labs{i}='All emissions';
    






%% Tidy plot
%Hawaii_ERF_bar_plot_horiz_TIDY
%ivar=1;
SW_paper_dSW_from_cloud_vars_plot_horiz_TIDY_AerChemMIP


savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_DAMIP_feedbacks_Region4_' var_str_tab ' ' land_ocean_str '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
opts.iplot_jpg = 0;


savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% END FIG

%% ACSIS synthesis paper - Emission contribution and feedback terms combined (and simplified).
iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
% N.B, - Changing the gca position in the TIDY script is the key to making the plots more
    %compact.
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break ytick_labs
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by
%fac=0.5;


% --- Data ---
switch vars{ivar}
    case {'SW'}
        ngaps = 4; %No. bars + 1
    otherwise
        ngaps = 5;
end
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;



i=i+1;
%x = x + gap{j-1} + gap{j}; xx(i)=x;
x = x + gap{j}; xx(i)=x;
val = AerChemGHG; uncer=AerChemGHG_un;
icol=4;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='GHG-only';

i=i+1;
%x = x + gap{j-1} + gap{j}; xx(i)=x;
x = x + gap{j}; xx(i)=x;
val = AerChemAerosol; uncer = AerChemAerosol_un;
icol=11;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='Aerosol-only';

i=i+1;
%x = x + gap{j-1} + gap{j}; xx(i)=x;
x = x + gap{j}; xx(i)=x;
val = netAerChemMIP; uncer = netAerChemMIP_un;
icol=3;
plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,end_cap_fraction_size);
ytick_labs{i}='All emissions';


switch vars{ivar}
        case {'SW'} %Can't estimate a forcing term on Nd since Nd is used to estimate the forcing.                        
            ngaps = 3; %No. bars + 1
            %j=j+1; section_width{j} = fac * 0.125; gap{j} = section_width{j}/ngaps;
            j=j+1; section_width{j} = 0.3; gap{j} = section_width{j}/ngaps;
            
            i=i+1;
            %x=x+gap{j}; xx(i)=x;
            %When starting a new section need to also add on the gap from
            %the last section to get to the horizontal line.
            x = x + gap{j-1} + gap{j}; xx(i)=x;
            val = ukesm_local_indirect + ukesm_local_direct; uncer=0;
            icol=1;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol ERF';
                                    
            
            %ngaps = 2; %No. bars + 1
            %j=j+1; section_width{j} = fac * 0.1; gap{j} = section_width{j}/ngaps;
            i=i+1;
            %x=x+gap{j-1}+gap{j}; xx(i)=x;           
            x=x+gap{j}; xx(i)=x;
            %val = AerChemGHG + feedback_estimate_AerChemAerosol; uncer=AerChemGHG_un; %using feedback from deltaT
            val = AerChemGHG + ukesm_non_local; uncer=AerChemGHG_un; %using feedback from total change minus forcing.
            icol=8;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol},'markersize',msize);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            %ytick_labs{i}='Aerosol Feedback from \DeltaT + GHG Feedbacks';
            ytick_labs{i}='Total (Aerosol+GHG) Feedback';
            

            
end


    






%% Tidy plot
%Hawaii_ERF_bar_plot_horiz_TIDY
%ivar=1;
SW_paper_dSW_from_cloud_vars_plot_horiz_TIDY_AerChemMIP


savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_feedbacks_Region4_SYNTHESIS_paper' var_str_tab ' ' land_ocean_str '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
opts.iplot_jpg = 0;


savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% END FIG

