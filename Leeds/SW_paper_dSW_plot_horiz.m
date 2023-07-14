%Is run from :- ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script.m

bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]};



clear ytick_labs
fsize_text = 14;

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



%% DAMIP sections

% GHG section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = ghg; uncer=ghg_un;
plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{6});
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{6},'o',3,0.005);
ytick_labs{i}='DAMIP Greenhouse';

%if isimplified==0 %Don't plot this if want a streamlined plot
% Nat section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = nat; uncer=nat_un;
plot(val,x,'o','color',bar_cols{8},'markerfacecolor',bar_cols{8});
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{8},'o',3,0.005);
ytick_labs{i}='DAMIP Natural';
%end

%% Aerosol section - N bars, N+1 gaps
switch vars{ivar}
    case 'SW'
        ngaps = 8; %No. bars + 1
        j=j+1; section_width{j} = fac * 1.0; gap{j} = section_width{j}/ngaps;
    otherwise
        ngaps = 2;
        j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps
end

%Bars within section
i=i+1;
x=x+gap{j-1} + gap{j}; xx(i)=x;
val = aer; uncer=aer_un;
plot(val,x,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4},'linewidth',1);
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{4},'o',3,0.005);
ytick_labs{i}='DAMIP Aerosol';

switch vars{ivar}
    case 'SW'
        i=i+1;
        x=x+gap{j}; xx(i)=x;
        val = local_indirect; uncer=0;
        plot(val,x,'s','color',bar_cols{6},'markerfacecolor',bar_cols{2});
        %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
        ytick_labs{i}='ACI Forcing';
        
        i=i+1;
        x=x+gap{j}; xx(i)=x;
        val = local_direct; uncer=0;
        plot(val,x,'^','color',bar_cols{6},'markerfacecolor',bar_cols{2});
        %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
        ytick_labs{i}='ARI Forcing';
        
        i=i+1;
        x=x+gap{j}; xx(i)=x;
        val = local_indirect+local_direct; uncer=0;
        plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{2},'markersize',12);
        %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
        ytick_labs{i}='ACI+ARI Forcing';
        
        i=i+1;
        x=x+gap{j}; xx(i)=x;
        val = non_local; uncer=0;
        plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{5});
        %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
        ytick_labs{i}='Aerosol Feedback';
        
        i=i+1;
        x=x+gap{j}; xx(i)=x;
        val = feedback_estimate_aer; uncer=0;
        plot(val,x,'s','color',bar_cols{6},'markerfacecolor',bar_cols{5});
        %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
        %ytick_labs{i}='DAMIP Aerosol Feedback Est';
        ytick_labs{i}='Aerosol Feedback from \DeltaT';
        
        i=i+1;
        x=x+gap{j}; xx(i)=x;
        val = local_indirect+local_direct + feedback_estimate_aer; uncer=0;
        plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
        %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
        ytick_labs{i}='Aerosol Forcing + Feedback';
        
        % -- end of aerosol section --
        
        
        ngaps = 2; %No. bars + 1
        j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
        %Bars within section
        i=i+1;
        x=x+gap{j-1}+gap{j}; xx(i)=x;
        val = ghg + feedback_estimate_aer; uncer=ghg_un;
        plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{9});
        %herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{6},'o',3,0.005);
        ytick_labs{i}='DAMIP Aerosol Feedback + GHG';
        
        
        
        
        ngaps = 2; %No. bars + 1
        j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
        %Bars within section
        i=i+1;
        x=x+gap{j-1}+gap{j}; xx(i)=x;
        val = ghg + local_indirect+local_direct + feedback_estimate_aer; uncer=ghg_un;
        icol=10;
        plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
        herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
        ytick_labs{i}='DAMIP Aerosol Forcing+Feedback + GHG';

        
        ngaps = 2; %No. bars + 1
        j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
        %Bars within section
        i=i+1;
        x=x+gap{j-1}+gap{j}; xx(i)=x;
        val = local_indirect+local_direct + feedback_estimate_all; uncer=0;
        icol=10;
        plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
        %herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
        ytick_labs{i}='Aerosol Forcing+All-emissions Feedback from \DeltaT';
        
end

%Net section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = net; uncer=net_un;
plot(val,x,'o','color',bar_cols{7},'markerfacecolor',bar_cols{7});
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{7},'o',3,0.005);
ytick_labs{i}='DAMIP Net (GHG+Natural+Aerosols)';

% HADGEM section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = had; uncer=had_un;
plot(val,x,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{3},'o',3,0.005);
ytick_labs{i}='All emissions (HADGEM model)';


%% Tidy plot
SW_paper_dSW_plot_horiz_TIDY

savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_DAMIP_Region4_' var_str_tab ' ' land_ocean '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% END FIG

%% New figure


clear ytick_labs
fsize_text = 14;

iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

clear xx section_width gap big_section_break
x=0;
i=0;
j=2; gap{1}=0;
k=1;
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by


if iagu==0
    
    %% Full minus DAMIP-hist-aer section (alternative method for calculating aerosol effects).
    %Aerosol section - 5 bar, 6 gaps
    switch vars{ivar}
        case 'SW'
            ngaps = 8; %No. bars + 1
            j=j+1; section_width{j} = fac * 1.0; gap{j} = section_width{j}/ngaps;
        otherwise
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    end
    
    
    %Bars within section
    i=i+1;
    x=x+gap{j-1} + gap{j}; xx(i)=x;
    val = aer2; uncer=aer2_un;
    plot(val,x,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4},'linewidth',1);
    herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{4},'o',3,0.005);
    ytick_labs{i}='Proxy Aerosol';
    
    
    switch vars{ivar}
        case 'SW'
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = local_indirect2; uncer=0;
            plot(val,x,'s','color',bar_cols{6},'markerfacecolor',bar_cols{2});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ACI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = local_direct2; uncer=0;
            plot(val,x,'^','color',bar_cols{6},'markerfacecolor',bar_cols{2});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ARI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = local_indirect2+local_direct2; uncer=0;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{2},'markersize',12);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ACI+ARI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = non_local2; uncer=0;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{5});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = feedback_estimate_aer2; uncer=0;
            plot(val,x,'s','color',bar_cols{6},'markerfacecolor',bar_cols{5});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback from \DeltaT';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = local_indirect2+local_direct2 + feedback_estimate_aer2; uncer=0;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Forcing + Feedback';
            
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
            i=i+1;
            x=x + gap{j-1} + gap{j}; xx(i)=x;
            val = ghg + feedback_estimate_aer2; uncer=ghg_un;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{9});
            %herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{6},'o',3,0.005);
            ytick_labs{i}='Proxy Aerosol Feedback + GHG';
            
            
            if isimplified==0 %Don't plot this if want a streamlined plot
                % Aerosol trend where do full trend minus GHG trend (aer3)
                ngaps = 2; %No. bars + 1
                j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
                i=i+1;
                x=x + gap{j-1} + gap{j}; xx(i)=x;
                val = aer3; uncer=aer3_un;
                plot(val,x,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4});
                herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{4},'o',3,0.005);
                ytick_labs{i}='DAMIP Aerosol Proxy';
            end
            
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
            %Bars within section
            i=i+1;
            x=x+gap{j-1}+gap{j}; xx(i)=x;
            val = ghg + local_indirect2+local_direct2 + feedback_estimate_aer2; uncer=ghg_un;
            icol=10;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
            herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
            ytick_labs{i}='Proxy Aerosol Forcing+Feedback+GHG';
                        
            
    end
    
end

% No point plotting this one since it is the same as the previous since the
% dT from the aerosol proxy is the difference between the full run and the
% GHG run.
% ngaps = 2; %No. bars + 1
% j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
% %Bars within section
% i=i+1;
% x=x+gap{j-1}+gap{j}; xx(i)=x;
% val = local_indirect2+local_direct2 + feedback_estimate_all; uncer=0;
% icol=10;
% plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
% %herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
% ytick_labs{i}='Proxy Aerosol Forcing+All-emissions Feedback from \DeltaT';


% Net using aer2
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = net2; uncer=net2_un;
plot(val,x,'o','color',bar_cols{7},'markerfacecolor',bar_cols{7});
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{7},'o',3,0.005);
ytick_labs{i}='Proxy Net (GHG+Natural+Proxy Aerosols)';

% HADGEM section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = had; uncer=had_un;
plot(val,x,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{3},'o',3,0.005);
ytick_labs{i}='All emissions (HADGEM model)';

 

%% Tidy and save plot
SW_paper_dSW_plot_horiz_TIDY

savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_Proxy_Region4_' var_str_tab ' ' land_ocean '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% New figure - UKESM / AerCheMIP


clear ytick_labs
fsize_text = 14;

iagu=0; %Simplified plot showing only the HADGEM and DAMIP values.
isimplified=1;

hf=figure;
set(gcf,'position',[1          26        1200         574]);
set(gcf,'color','w');
set(gca,'ylim',[0 15]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

clear xx section_width gap big_section_break
x=0;
i=0;
j=2; gap{1}=0;
k=1;
gap{j}=0; section_width{j}=0;
big_section_break=[];

fac=1.6; %factor to scale the gaps by


if iagu==0
    
   
    
    % % Test for a gap
    % ngaps = 2; %No. bars + 1
    % j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    % i=i+1;
    % x=x + gap{j-1} + gap{j}; xx(i)=x;
    % %val = net2; uncer=net2_un;
    % %plot(val,x,'o','color',bar_cols{7},'markerfacecolor',bar_cols{7});
    % %herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{7},'o',3,0.005);
    % ytick_labs{i}='';
    
    
    
    
    %% UKESM (AerChemMIP) section
    big_section_break(k) = x+gap{j-1}; k=k+1;
    % Aerosol (AerChemMIP) section
    ngaps = 2; %No. bars + 1
    j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    %Bars within section
    i=i+1;
    x=x+gap{j-1}+gap{j}; xx(i)=x;
    val = AerChemGHG; uncer=AerChemGHG_un;
    plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{6});
    herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{6},'o',3,0.005);
    ytick_labs{i}='AerChemMIP Greenhouse';
    
    % % Aerosol (AerChemMIP) section
    % ngaps = 2; %No. bars + 1
    % j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    % %Bars within section
    % i=i+1;
    % x=x+gap{j-1}+gap{j}; xx(i)=x;
    % val = AerChemAerosol; uncer=AerChemAerosol_un;
    % plot(val,x,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4});
    % herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{4},'o',3,0.005);
    % ytick_labs{i}='Aerosol AerChemMIP';
    
    
    
    
    %% Aerosol section - N bars, N+1 gaps
    switch vars{ivar}
        case 'SW'
            ngaps = 8; %No. bars + 1
            j=j+1; section_width{j} = fac * 1; gap{j} = section_width{j}/ngaps;
        otherwise
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    end
    
    %Bars within section
    i=i+1;
    x=x+gap{j-1} + gap{j}; xx(i)=x;
    val = AerChemAerosol; uncer=AerChemAerosol_un;
    plot(val,x,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4},'linewidth',1);
    herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{4},'o',3,0.005);
    ytick_labs{i}='AerChemMIP Aerosol';
    
    switch vars{ivar}
        case 'SW'
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_local_indirect; uncer=0;
            plot(val,x,'s','color',bar_cols{6},'markerfacecolor',bar_cols{2});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ACI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_local_direct; uncer=0;
            plot(val,x,'^','color',bar_cols{6},'markerfacecolor',bar_cols{2});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ARI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_local_indirect + ukesm_local_direct; uncer=0;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{2},'markersize',12);
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='ACI+ARI Forcing';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_non_local; uncer=0;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{5});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = feedback_estimate_AerChemAerosol; uncer=0;
            plot(val,x,'s','color',bar_cols{6},'markerfacecolor',bar_cols{5});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Feedback from \DeltaT';
            
            i=i+1;
            x=x+gap{j}; xx(i)=x;
            val = ukesm_local_indirect + ukesm_local_direct + feedback_estimate_AerChemAerosol; uncer=0;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='Aerosol Forcing + Feedback';
            
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
            i=i+1;
            x=x+gap{j-1}+gap{j}; xx(i)=x;
            val = AerChemGHG + feedback_estimate_AerChemAerosol; uncer=AerChemGHG_un;
            plot(val,x,'o','color',bar_cols{6},'markerfacecolor',bar_cols{9});
            %herr=errorbarYY('horiz',val,x,uncer,gca,'b','o',3,0.005);
            ytick_labs{i}='AerChemMIP Aerosol Feedback + GHG';
            
            % -- end of the aerosol section for UKESM/AerChemMIP
            
            if isimplified==0 %Don't plot this if want a streamlined plot
                % Aerosol (AerChemMIP) diff in trends section
                ngaps = 2; %No. bars + 1
                j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
                %Bars within section
                i=i+1;
                x=x+gap{j-1}+gap{j}; xx(i)=x;
                val = AerChemAerosol2; uncer=AerChemAerosol2_un;
                plot(val,x,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4});
                herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{4},'o',3,0.005);
                ytick_labs{i}='Aerosol: UKESM minus piAer trend diff';
            end
            
            ngaps = 2; %No. bars + 1
            j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
            %Bars within section
            i=i+1;
            x=x+gap{j-1}+gap{j}; xx(i)=x;
            val = AerChemGHG + feedback_estimate_AerChemAerosol + ukesm_local_indirect + ukesm_local_direct; uncer=AerChemGHG_un;
            icol=10;
            plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
            herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
            ytick_labs{i}='AerChemMIP Aerosol Forcing+Feedback + GHG';
            
    end
    
% No point plotting this one since it is the same as the previous since the
% dT from the aerosol proxy is the difference between the full run and the
% GHG run.
%     ngaps = 2; %No. bars + 1
%     j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
%     %Bars within section
%     i=i+1;
%     x=x+gap{j-1}+gap{j}; xx(i)=x;
%     val = ukesm_local_indirect + ukesm_local_direct + feedback_estimate_AerChem_all; uncer=0;
%     icol=10;
%     plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
%     %herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
%     ytick_labs{i}='Aerosol Forcing+All-emissions Feedback from \DeltaT';

    
    ngaps = 2; %No. bars + 1
    j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    %Bars within section
    i=i+1;
    x=x+gap{j-1}+gap{j}; xx(i)=x;
    val = netAerChemMIP; uncer=netAerChemMIP_un;
    icol=7;
    plot(val,x,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
    herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{icol},'o',3,0.005);
    ytick_labs{i}='AerChemMIP Net';
    
    
    ngaps = 2; %No. bars + 1
    j=j+1; section_width{j} = fac * 0.25; gap{j} = section_width{j}/ngaps;
    %Bars within section
    i=i+1;
    x=x+gap{j-1}+gap{j}; xx(i)=x;
    val = UKESMAerChemMIP; uncer=UKESMAerChemMIP_un;
    plot(val,x,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
    herr=errorbarYY('horiz',val,x,uncer,gca,bar_cols{3},'o',3,0.005);
    ytick_labs{i}='All emissions (UKESM model subset)';
    

    
end

%% Tidy plot
SW_paper_dSW_plot_horiz_TIDY

savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_UKESM_Region4_' var_str_tab ' ' land_ocean '_Period' period_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
