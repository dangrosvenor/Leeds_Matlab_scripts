% Plot PDF of trends from the ensemble.
Nens = size(trend_dat_box_ens,3);

clear dat
for iens=1:Nens
   
    dat(iens) = trend_dat_box_ens{1,1,iens}.coeffs(2); %trend values
    
end

figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.

%Box plot of the data
boxplot(dat,'orientation','horizontal');
set(gca,'fontsize',18);
me = mean(dat);
hold on
plot(me,1,'ko');



xlabel('SW TOA trend (W m^{-2} yr^{-1})');

%% Plot on observed data and range of 95% CI
obs_trend = trend_dat_box_obs2{ibox,1}.coeffs(2);
obs_CI = trend_dat_box_obs2{ibox,1}.uncer_max;
yloc=1.3;
herr=errorbarYY('horiz',obs_trend,yloc,obs_CI,gca,'k','o',2,0.05);
plot(obs_trend,yloc,'ko','markerfacecolor','k');



%% Plot smallest trend member
[maxval,imax] = max(dat);

obs_trend = dat(imax);
obs_CI = trend_dat_box_ens{1,1,imax}.uncer_max;
yloc2=1.15;
herr=errorbarYY('horiz',obs_trend,yloc2,obs_CI,gca,'k','o',2,0.05);
plot(obs_trend,yloc2,'ko','markerfacecolor','k');



%% Tidy plot
minx = min(dat); minx = min(minx,obs_trend-obs_CI);
maxx = max(dat); maxx = max(maxx,obs_trend+obs_CI);
xrange = maxx-minx;
set(gca,'xlim',[minx-xrange*0.1 maxx+xrange*0.1]);
set(gca,'ylim',[0.85 yloc+0.1]);

set(gca,'ytick',[1 yloc2 yloc]);
set(gca,'yticklabel',{'UKESM ensemble','Extreme member','Observed'});