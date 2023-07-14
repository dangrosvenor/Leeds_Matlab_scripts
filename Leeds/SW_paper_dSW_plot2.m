 bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]};



clear xtick_labs


    
hf=figure;
set(gcf,'color','w');
set(gca,'xlim',[0 5]); %Need to do this since errobarYY plots as a fraction of the xlim at the time of plotting...
hold on
%set(gcf,'position',[3         297        1256         422]);

x=0;
i=0;
j=1;
k=1;
clear xx section_width gap big_section_break
gap{j}=0; section_width{j}=0;



% GHG section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = ghg; uncer=ghg_un;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{6});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{6},'o',3,0.005);
xtick_labs{i}='DAMIP Greenhouse';

% Nat section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = nat; uncer=nat_un;
plot(x,val,'o','color',bar_cols{8},'markerfacecolor',bar_cols{8});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{8},'o',3,0.005);
xtick_labs{i}='DAMIP Natural';

%% Aerosol section - N bars, N+1 gaps
ngaps = 9; %No. bars + 1
j=j+1; section_width{j} = 1.0; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1} + gap{j}; xx(i)=x;
val = aer; uncer=aer_un;
plot(x,val,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4},'linewidth',1);
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{4},'o',3,0.005);
xtick_labs{i}='DAMIP Aerosol';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_indirect; uncer=0;
plot(x,val,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ACI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_direct; uncer=0;
plot(x,val,'o','color',bar_cols{2},'markerfacecolor',bar_cols{2});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ARI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_indirect+local_direct; uncer=0;
plot(x,val,'o','color',bar_cols{9},'markerfacecolor',bar_cols{9});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ACI+ARI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = non_local; uncer=0;
plot(x,val,'o','color',bar_cols{2},'markerfacecolor',bar_cols{5});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='DAMIP Aerosol Feedback';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = feedback_estimate_aer; uncer=0;
plot(x,val,'o','color',bar_cols{2},'markerfacecolor',bar_cols{5});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='DAMIP Aerosol Feedback Est';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_indirect+local_direct + feedback_estimate_aer; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='Aerosol Forcing + Feedback';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ghg + feedback_estimate_aer; uncer=ghg_un;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
%herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{6},'o',3,0.005);
xtick_labs{i}='Aerosol Feedback + GHG';

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = ghg + local_indirect+local_direct + feedback_estimate_aer; uncer=ghg_un;
icol=7;
plot(x,val,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{icol},'o',3,0.005);
xtick_labs{i}='Aerosol Forcing+Feedback + GHG';



%% Full minus DAMIP-hist-aer section (alternative method for calculating aerosol effects).
%Aerosol section - 5 bar, 6 gaps
ngaps = 9; %No. bars + 1
j=j+1; section_width{j} = 1.0; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1} + gap{j}; xx(i)=x;
val = aer2; uncer=aer2_un;
plot(x,val,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4},'linewidth',1);
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{4},'o',3,0.005);
xtick_labs{i}='HADGEM minus DAMIP GHG';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_indirect2; uncer=0;
plot(x,val,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ACI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_direct2; uncer=0;
plot(x,val,'o','color',bar_cols{2},'markerfacecolor',bar_cols{2});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ARI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_indirect2+local_direct2; uncer=0;
plot(x,val,'o','color',bar_cols{9},'markerfacecolor',bar_cols{9});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ACI+ARI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = non_local2; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{5});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='HADGEM minus DAMIP GHG Aerosol Feedback';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = feedback_estimate_aer2; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{5});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='Aerosol Feedback Est';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = local_indirect2+local_direct2 + feedback_estimate_aer2; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='Aerosol Forcing + Feedback';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ghg + feedback_estimate_aer2; uncer=ghg_un;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
%herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{6},'o',3,0.005);
xtick_labs{i}='Aerosol Feedback + GHG';

% Aerosol trend where do full trend minus GHG trend (aer3)
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = aer3; uncer=aer3_un;
plot(x,val,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{4},'o',3,0.005);
xtick_labs{i}='HADGEM trend minus DAMIP-GHG trend';

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = ghg + local_indirect2+local_direct2 + feedback_estimate_aer2; uncer=ghg_un;
icol=7;
plot(x,val,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{icol},'o',3,0.005);
xtick_labs{i}='Aerosol Forcing+Feedback + GHG';

% HADGEM section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = had; uncer=had_un;
plot(x,val,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{3},'o',3,0.005);
xtick_labs{i}='HADGEM';

%Net section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = net; uncer=net_un;
plot(x,val,'o','color',bar_cols{7},'markerfacecolor',bar_cols{7});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{7},'o',3,0.005);
xtick_labs{i}='DAMIP Net';

% Net using aer2
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;

i=i+1;
x=x + gap{j-1} + gap{j}; xx(i)=x;
val = net2; uncer=net2_un;
plot(x,val,'o','color',bar_cols{7},'markerfacecolor',bar_cols{7});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{7},'o',3,0.005);
xtick_labs{i}='Net2';





%% UKESM (AerChemMIP) section
big_section_break(k) = x+gap{j-1}; k=k+1;
% Aerosol (AerChemMIP) section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = AerChemGHG; uncer=AerChemGHG_un;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{6});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{6},'o',3,0.005);
xtick_labs{i}='GHG AerChemMIP';

% % Aerosol (AerChemMIP) section
% ngaps = 2; %No. bars + 1
% j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
% %Bars within section
% i=i+1;
% x=x+gap{j-1}+gap{j}; xx(i)=x;
% val = AerChemAerosol; uncer=AerChemAerosol_un;
% plot(x,val,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4});
% herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{4},'o',3,0.005);
% xtick_labs{i}='Aerosol AerChemMIP';




%% Aerosol section - N bars, N+1 gaps
ngaps = 9; %No. bars + 1
j=j+1; section_width{j} = 1; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1} + gap{j}; xx(i)=x;
val = AerChemAerosol; uncer=AerChemAerosol_un;
plot(x,val,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4},'linewidth',1);
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{4},'o',3,0.005);
xtick_labs{i}='Aerosol AerChemMIP';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ukesm_local_indirect; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{3});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ACI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ukesm_local_direct; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{2});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ARI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ukesm_local_indirect + ukesm_local_direct; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{2});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='ACI+ARI Forcing';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ukesm_non_local; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{2});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='AerChemMIP Aerosol Feedback';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = feedback_estimate_AerChemAerosol; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{5});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='AerChemMIP Aerosol Feedback Est';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = ukesm_local_indirect + ukesm_local_direct + feedback_estimate_AerChemAerosol; uncer=0;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='Aerosol Forcing + Feedback';

i=i+1;
x=x+gap{j}; xx(i)=x;
val = AerChemGHG + feedback_estimate_AerChemAerosol; uncer=AerChemGHG_un;
plot(x,val,'o','color',bar_cols{6},'markerfacecolor',bar_cols{10});
%herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.005);
xtick_labs{i}='Aerosol Feedback + GHG';

% -- end of the aerosol section for UKESM/AerChemMIP


% Aerosol (AerChemMIP) diff in trends section
ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = AerChemAerosol2; uncer=AerChemAerosol2_un;
plot(x,val,'o','color',bar_cols{4},'markerfacecolor',bar_cols{4});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{4},'o',3,0.005);
xtick_labs{i}='Aerosol: UKESM minus piAer trend diff';

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = AerChemGHG + feedback_estimate_AerChemAerosol + ukesm_local_indirect + ukesm_local_direct; uncer=AerChemGHG_un;
icol=7;
plot(x,val,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{icol},'o',3,0.005);
xtick_labs{i}='Aerosol Forcing+Feedback + GHG';

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = UKESMAerChemMIP; uncer=UKESMAerChemMIP_un;
plot(x,val,'o','color',bar_cols{3},'markerfacecolor',bar_cols{3});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{3},'o',3,0.005);
xtick_labs{i}='Full UKESM AerChemMIP';

ngaps = 2; %No. bars + 1
j=j+1; section_width{j} = 0.25; gap{j} = section_width{j}/ngaps;
%Bars within section
i=i+1;
x=x+gap{j-1}+gap{j}; xx(i)=x;
val = netAerChemMIP; uncer=netAerChemMIP_un;
icol=7;
plot(x,val,'o','color',bar_cols{icol},'markerfacecolor',bar_cols{icol});
herr=errorbarYY('vert',x,val,uncer,gca,bar_cols{icol},'o',3,0.005);
xtick_labs{i}='Net AerChemMIP';



%% Tidy plot
yoff=0.28;
set(gca,'position',[0.1367    0.1100+yoff    0.7683    0.8150-(yoff*1.2)])
ylabel(ylabs{ivar});
%xlabel('Latitude');
fontsize_figure(gcf,gca,18);
title(period_lab); %set in ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot.m
grid off

set(gca,'xlim',[0 xx(end)+0.2]);

%set(gca,'xtick',[1:length(xtick_labs)]);
set(gca,'xtick',xx);
set(gca,'xticklabel',[]);

% Add tick labels at 45 degree angles.
ylims=get(gca,'ylim');
xlims=get(gca,'xlim');
dylim = ylims(2)-ylims(1);
ypos=ylims(1)-dylim*0.04;
for i=1:length(xtick_labs)
    text(xx(i),ypos,xtick_labs{i},'rotation',-45,'fontsize',14);
end

%Plot the vertical line breaks
x=0;
plot([x x],[ylims(1) ylims(2)],'k-');
for i=2:length(section_width)
    x = x + section_width{i};
    plot([x x],[ylims(1) ylims(2)],'k-');
    
    %hL=line([i-0.5 i-0.5],[ylims(1) ylims(2)]);
    %set(hL,'color','k');
end

%Plot the zero line at y=0
plot([xlims(1) xlims(2) ],[0 0],'k--');
set(gca,'XminorTick','off');
%hL=line([xlims(1) xlims(2)],[0 0]);
%set(hL,'color','k');
% 

for i=1:length(big_section_break)
    x = big_section_break(i);
    plot([x x],[ylims(1) ylims(2)],'k-','linewidth',6);
    
    %hL=line([i-0.5 i-0.5],[ylims(1) ylims(2)]);
    %set(hL,'color','k');
end


