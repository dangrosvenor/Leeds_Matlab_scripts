% Plot bar plot, error bars, etc. of trends from the ensemble.

iadd_obs=1;
iadd_amip=1;

switch corr_vars
    case {'ts,sw'}
        hbar_size=0.02;
        xlab = 'Spatial r for SWTOA vs ts trends'
        %        obs_trend = trend_dat_box_obs2{ibox,1}.coeffs(2);
        %        obs_CI = trend_dat_box_obs2{ibox,1}.uncer_max;
        switch land_ocean
            case 'land only'
                xlims = [-0.8 0];
            case 'ocean only'
                xlims = [-0.5 0.1];
        end
        
    case 'ts,ts_obs'
        hbar_size=0.02;
        xlab = 'Spatial r for ts trends model vs obs'
        %        obs_trend = trend_dat_box_obs2{ibox,1}.coeffs(2);
        %        obs_CI = trend_dat_box_obs2{ibox,1}.uncer_max;
        switch land_ocean
            case 'land only'
                xlims = [-0.8 0];
            case 'ocean only'
                xlims = [0 0.7];
        end
        
                
    otherwise
        error('Need to select a variable here');
end



%% Model box and whisker plot of indvidual ensemble means
Nens = length(corr_vals_ens);
dat = corr_vals_ens;

figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
set(gcf,'position',[374   236   838   384]);
%Box plot of the data
boxplot(dat,'orientation','horizontal','color','b');
set(gca,'fontsize',18);
me = mean(dat);
hold on
%Plot ensemble mean and CI range
plot(me,1,'bo','markerfacecolor','b','markersize',12);

%% Plot trend from pre-averaged ensemble
yloc_mean = 0.85;
%plot(trend_dat_box{ibox,1}.coeffs(2),yloc_mean,'bo','markerfacecolor','b','markersize',12);
%herr=errorbarYY('horiz',trend_dat_box{ibox,1}.coeffs(2),yloc_mean,trend_dat_box{ibox,1}.uncer_max,gca,'b','o',2,hbar_size);



%% Plot on observed data and range of 95% CI
yloc=1.3;
if iadd_obs==1
    %herr=errorbarYY('horiz',corr_val_obs,yloc,obs_CI,gca,'k','o',2,hbar_size);
    plot(corr_val_obs,yloc,'ko','markerfacecolor','k','markersize',12);
    %obs_min = obs_trend - obs_CI;
    %obs_max = obs_trend + obs_CI;
    obs_str = 'Observed';
else
    obs_str = '';
end


%% Plot largest model trend member
%[maxval,imax] = max(dat);

%obs_trend = dat(imax);
%obs_CI = trend_dat_box_ens{ibox,1,imax}.uncer_max;
yloc_max=1.17;
%herr=errorbarYY('horiz',obs_trend,yloc_max,obs_CI,gca,'b','o',2,hbar_size);
%plot(obs_trend,yloc_max,'bo','markerfacecolor','b','markersize',12);

%ens_max = obs_trend + obs_CI;


%% Plot smallest trend member
% [maxval,imax] = min(dat);
% 
% obs_trend = dat(imax);
% obs_CI = trend_dat_box_ens{ibox,1,imax}.uncer_max;
 yloc_min=1.13;
% herr=errorbarYY('horiz',obs_trend,yloc_min,obs_CI,gca,'b','o',2,hbar_size);
% plot(obs_trend,yloc_min,'bo','markerfacecolor','b','markersize',12);
% 
 yloc2=mean([yloc_min yloc_max]);
% 
% ens_min = obs_trend - obs_CI;

%% AMIP
 if iadd_amip==1
     yloc_amip = 0.7;
     plot(corr_val_amip,yloc_amip,'ro','markerfacecolor','r','markersize',12);
%     herr=errorbarYY('horiz',trend_dat_box_amip{ibox,1}.coeffs(2),yloc_amip,trend_dat_box_amip{ibox,1}.uncer_max,gca,'r','o',2,hbar_size);   
 end

%% Tidy plot
%minx = min(dat); minx = min(minx,obs_trend-obs_CI);
%maxx = max(dat); maxx = max(maxx,obs_trend+obs_CI);
%minx = min([min(dat) obs_min ens_min]);
%maxx = max([min(dat) obs_max ens_max]);

%xrange = maxx-minx;
%set(gca,'xlim',[minx-xrange*0.1 maxx+xrange*0.1]);
set(gca,'xlim',xlims);

if iadd_amip==1
    set(gca,'ytick',[yloc_amip yloc_mean 1 yloc2 yloc]);
    set(gca,'yticklabel',{'UKESM-AMIP' '' 'UKESM ensemble','',obs_str});
    set(gca,'ylim',[0.65 yloc+0.1]);
else
    set(gca,'ytick',[yloc_mean 1 yloc2 yloc]);
    %set(gca,'yticklabel',{'Trend from ensemble mean''UKESM ensemble','Extreme members','Observed'});
    set(gca,'yticklabel',{'','UKESM ensemble','',obs_str});
    set(gca,'ylim',[0.75 yloc+0.1]);
end
yrange = get(gca,'ylim');
xlabel({' ' xlab}); %Make a gap between x-axis and xlabel
title({['Region=' box_region ' ' land_ocean],'',[box_region_str ' ' season ' ' land_ocean]});


savename = [savedir_date xlab ' ' box_region_str ' ' season ' ' land_ocean];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
