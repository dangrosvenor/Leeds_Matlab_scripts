% Plot bar plot, error bars, etc. of trends from the ensemble.

hbar_size=0.02;


Nens = size(trend_dat_box_ens,3);

clear dat
for iens=1:Nens
   
    dat(iens) = trend_dat_box_ens{1,1,iens}.coeffs(2); %trend values
    
end

figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
set(gcf,'position',[374   236   838   384]);
%Box plot of the data
boxplot(dat,'orientation','horizontal');
set(gca,'fontsize',18);
me = mean(dat);
hold on
%Plot ensemble mean and CI range
plot(me,1,'ko','markerfacecolor','k','markersize',12);

%% Plot trend from pre-averaged ensemble
yloc_mean = 0.85;
plot(trend_dat_box{1,1}.coeffs(2),yloc_mean,'ko','markerfacecolor','k','markersize',12);
herr=errorbarYY('horiz',trend_dat_box{1,1}.coeffs(2),yloc_mean,trend_dat_box{1,1}.uncer_max,gca,'k','o',2,hbar_size);


xlabel('Total cloud fraction trend (yr^{-1})');

%% Plot on observed data and range of 95% CI
obs_trend = trend_dat_box_obs_cci{ibox,1}.coeffs(2);
obs_CI = trend_dat_box_obs_cci{ibox,1}.uncer_max;
yloc=1.3;
herr=errorbarYY('horiz',obs_trend,yloc,obs_CI,gca,'k','o',2,hbar_size);
plot(obs_trend,yloc,'ko','markerfacecolor','k','markersize',12);



%% Plot largest trend member
[maxval,imax] = max(dat);

obs_trend = dat(imax);
obs_CI = trend_dat_box_ens{1,1,imax}.uncer_max;
yloc_max=1.17;
herr=errorbarYY('horiz',obs_trend,yloc_max,obs_CI,gca,'k','o',2,hbar_size);
plot(obs_trend,yloc_max,'ko','markerfacecolor','k','markersize',12);


%% Plot smallest trend member
[maxval,imax] = min(dat);

obs_trend = dat(imax);
obs_CI = trend_dat_box_ens{1,1,imax}.uncer_max;
yloc_min=1.13;
herr=errorbarYY('horiz',obs_trend,yloc_min,obs_CI,gca,'r','o',2,hbar_size);
plot(obs_trend,yloc_min,'ro','markerfacecolor','r','markersize',12);

yloc2=mean([yloc_min yloc_max]);

%% Tidy plot
minx = min(dat); minx = min(minx,obs_trend-obs_CI);
maxx = max(dat); maxx = max(maxx,obs_trend+obs_CI);
xrange = maxx-minx;
set(gca,'xlim',[minx-xrange*0.1 maxx+xrange*0.1]);
set(gca,'ylim',[0.75 yloc+0.1]);

set(gca,'ytick',[yloc_mean 1 yloc2 yloc]);
set(gca,'yticklabel',{'Trend from ensemble mean' 'UKESM ensemble','Extreme members','Observed'});
set(gca,'xlim',[-20e-4 0]);