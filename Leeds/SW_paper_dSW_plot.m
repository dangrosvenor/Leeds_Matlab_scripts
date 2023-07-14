

dT1 = (1970-1850);
dT2 = (2014-1971);

clear xtick_labs

%% Period 1
hf=figure;
set(gcf,'color','w');
%set(gcf,'position',[3         297        1256         422]);

dT = dT1; p=1; ibox=1;

x=1;
val = dT*trend_dat_box{ibox,p}.coeffs(2); uncer = dT*trend_dat_box{ibox,p}.uncer_max;
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='UKESM1';

hold on
yoff=0.28;
set(gca,'position',[0.1367    0.1100+yoff    0.7683    0.8150-(yoff*1.2)])
ylabel('\DeltaF_{SW} (W m^{-2})');
%xlabel('Latitude');
fontsize_figure(gcf,gca,18);
title('1850 to 1970');
grid on


x=2;
val = 3.6; uncer=0.99;
plot(x,3.6,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='HADGEM';

x=3;
val = 5.3; uncer = 0.81;
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='DAMIP hist-aer';

x=4;
aero_piaero = dT*(trend_dat_box{ibox,p}.coeffs(2) - trend_dat_box_hist_piaer{ibox,p}.coeffs(2) ); 
val = aero_piaero;
uncer = dT* (trend_dat_box{ibox,p}.uncer_max + trend_dat_box_hist_piaer{ibox,p}.uncer_max ); 
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='Aerosol (UKESM1 - histPIaer)';

x=5;
val = -1.9; uncer = 0.46;
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='DAMIP hist-GHG';

x=6;
ghg_piaero = dT*(trend_dat_box_hist_piaer{ibox,p}.coeffs(2) ); 
val = ghg_piaero;
uncer = dT* (trend_dat_box_hist_piaer{ibox,p}.uncer_max ); 
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='GHGs (histPIaer)';

set(gca,'xtick',[1:length(xtick_labs)]);
set(gca,'xticklabel',[]);

ylims=get(gca,'ylim');
dylim = ylims(2)-ylims(1);
ypos=ylims(1)-dylim*0.04;
for i=1:length(xtick_labs)
    text(i,ypos,xtick_labs{i},'rotation',-45,'fontsize',14);
end


%% Period 2
hf=figure;
set(gcf,'color','w');
%set(gcf,'position',[3         297        1256         422]);

dT = dT2; p = 2; ibox=1;

x=1;
plot(x,dT*trend_dat_box{ibox,p}.coeffs(2),'o','markerfacecolor','b');
herr=errorbarYY('vert',x,dT*trend_dat_box{ibox,p}.coeffs(2),dT*trend_dat_box{ibox,p}.uncer_max,gca,'b','o',3,0.01);
xtick_labs{x}='UKESM1';

hold on
yoff=0.28;
set(gca,'position',[0.1367    0.1100+yoff    0.7683    0.8150-(yoff*1.2)])
ylabel('\DeltaF_{SW} (W m^{-2})');
%xlabel('Latitude');
fontsize_figure(gcf,gca,18);
title('1971 to 2014');
grid on


x=2;
val = -4.6; uncer = 1.8;
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='HADGEM';

x=3;
val = 0.09; uncer = 1.1;
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='DAMIP hist-aer';

x=4;
aero_piaero = dT*(trend_dat_box{ibox,p}.coeffs(2) - trend_dat_box_hist_piaer{ibox,p}.coeffs(2) ); 
val = aero_piaero;
uncer = dT* (trend_dat_box{ibox,p}.uncer_max + trend_dat_box_hist_piaer{ibox,p}.uncer_max ); 
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='Aerosol (UKESM1 - histPIaer)';

x=5;
val = -2.3; uncer = 0.69;
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='DAMIP hist-GHG';

x=6;
ghg_piaero = dT*(trend_dat_box_hist_piaer{ibox,p}.coeffs(2) ); 
val = ghg_piaero;
uncer = dT* (trend_dat_box_hist_piaer{ibox,p}.uncer_max ); 
plot(x,val,'o','markerfacecolor','b');
herr=errorbarYY('vert',x,val,uncer,gca,'b','o',3,0.01);
xtick_labs{x}='GHGs (histPIaer)';

set(gca,'xtick',[1:length(xtick_labs)]);
set(gca,'xticklabel',[]);

ylims=get(gca,'ylim');
dylim = ylims(2)-ylims(1);
ypos=ylims(1)-dylim*0.04;
for i=1:length(xtick_labs)
    text(i,ypos,xtick_labs{i},'rotation',-45,'fontsize',14);
end














