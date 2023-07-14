iplot_extreme=0;
spacing = 0.15;
iplot_delta = 0; %Whether to plot the delta or the trend.
itrends_per_decade = 1; %Whether to plot per decade (=1) or per year (=0)
msize = 8; %Originally was 10, but it covered the box and whisker lines in the sub-plot version.
iadd_Tcorr_points=1;

if iplot_delta==1
    % *** Assuming here that the time periods are the same for each bar plotted
    %here. ***
    dT = trend_dat_box{ibox,itrend_box_whisker}.x(end) - trend_dat_box{ibox,itrend_box_whisker}.x(1);
end


minx=1e99;
maxx=-1e99;
Nbox_plots = 0;
itick = 0;
clear ypos_box yticks


xscale = 1; %factor to multiply the x values by
% Plot bar plot, error bars, etc. of trends from the ensemble.
switch var_ukesm
    case {'calipso_total_cloud_amount','clt'}
        hbar_size=0.02;
           
        %xlims = [-20e-4 0];
        if itrends_per_decade==1 & iplot_delta==0
            xscale = 1e3;
            %xlab_per_decade = 'Total cloud fraction trend (10^{-3} decade^{-1})';
            xlab_per_decade = [var_str ' trend (10^{-3} decade^{-1})'];
        else
            xscale = 1e4; %factor to multiply the x values by
            xlab = [var_str ' trend (10^{-4} yr^{-1})'];
            xlab_per_decade = [var_str ' trend (10^{-4} decade^{-1})'];
            xlab_delta = ['\Delta' var_str ' (10^{-4})'];
        end
        
    case {'SW_up_TOA','rsut'}
        hbar_size=0.02;       
        %xlims = [-20e-4 0];
        xlab = [var_str ' trend (W m^{-2} yr^{-1})'];
        xlab_per_decade = [var_str ' trend (W m^{-2} decade^{-1})'];
        xlab_delta = ['\Delta' var_str ' (W m^{-2})'];
        
    case {'ts'}
        hbar_size=0.02;       
        %xlims = [-20e-4 0];
        xlab =  [var_str ' trend (K yr^{-1})'];
        xlab_per_decade = [var_str ' trend (K decade^{-1})'];
        xlab_delta = ['\Delta' var_str ' (K)'];        
        
    case {'Nd_cf_weighted_UKESM','scldncl','Nd_clw_weighted_ESGF'}
        hbar_size=0.02;
        xlab = 'N_d trend (cm^{-3} yr^{-1})';
        xlab_per_decade = 'N_d trend (cm^{-3} decade^{-1})';
        xlab_delta ='\DeltaN_d (cm^{-3})';
        
    case 'dust_od550'
        hbar_size=0.02;
        xlab = 'Dust 550nm optical depth trend (yr^{-1})';
        xlab_per_decade = 'Dust 550nm optical depth trend (decade^{-1})';
        xlab_delta ='\DeltaDust 550nm optical depth'
        
    case 'od550tot'
        hbar_size=0.02;
        xlab = [var_str ' trend (yr^{-1})'];
        xlab_per_decade = [var_str ' trend (decade^{-1})'];
        xlab_delta = ['\Delta' var_str];
        
    case {'ts'}
        hbar_size=0.02;
        %xlims = [-20e-4 0];
        xlab = [var_str ' trend (K yr^{-1})'];   
        xlab_per_decade = [var_str ' trend (K decade^{-1})'];           
        xlab_delta = ['\Delta' var_str ' (K)']; 
        
    case 'lwp'
        hbar_size=0.02;
        %xlims = [-20e-4 0];
        xlab = [var_str ' trend (g m^{-2} yr^{-1})'];
        xlab_per_decade = [var_str ' trend (g m^{-2} decade^{-1})'];
        xlab_delta = ['\Delta' var_str ' (g m^{-2})'];
        
    case 'lwpic'
        hbar_size=0.02;        
        %xlims = [-20e-4 0];
        xlab = 'LWPic trend (g m^{-2} yr^{-1})'; 
        xlab_per_decade = 'LWPic trend (g m^{-2} decade^{-1})'; 
        xlab_delta = '\DeltaLWPic (g m^{-2})';
        
    case 'prw'
        hbar_size=0.02;        
        %xlims = [-20e-4 0];
        xlab = 'Water Vapour Path trend (kg m^{-2} yr^{-1})';     
        xlab_per_decade = 'Water Vapour Path trend (kg m^{-2} decade^{-1})';     
        xlab_delta = '\DeltaWater Vapour Path (kg m^{-2})';
        
    otherwise
        %error('Need to select a variable here');
        hbar_size=0.02;
        xlab = [var_str 'trend'];
        xlab_per_decade = '';
        xlab_delta = ['\Delta' var_str];
end

if iplot_delta==1
   xscale = xscale * dT; %multiply by the time delta to get the change rather than the trend. 
   xlab = xlab_delta;
elseif itrends_per_decade
    xscale = xscale * 10;
    xlab = xlab_per_decade;
end

%% Set up figure
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
set(gcf,'position',[374   236   838   384]);
set(gca,'ylim',[0.6 1.6]); %temporary setting of ylim to get the errorbars the right size

if iplot_extreme==1
    yloc=1.3 + 0.15;
else
    yloc=1.15 + 0.15;
end

hold on

%% Plot on observed data no. 2 and range of 95% CI
%% Plot on observed data no.1 and range of 95% CI
if inc_obs_01==1
    yloc = yloc - spacing;
    
    hbar_size=0.02;
    obs_trend_01 = xscale*trend_dat_box_obs1{ibox,1}.coeffs(2);
    obs_CI_01 = xscale*trend_dat_box_obs1{ibox,1}.uncer_max;
    
    [herr,minx2,maxx2]=errorbarYY('horiz',obs_trend_01,yloc,obs_CI_01,gca,col_01,'o',2,hbar_size);
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    plot(obs_trend_01,yloc,'o','color',col_01,'markerfacecolor',col_01,'markersize',msize);
    itick = itick + 1;
    ytick_strs{itick}=[obs_str1];
    
    obs_min = obs_trend_01 - obs_CI_01;
    obs_max = obs_trend_01 + obs_CI_01;
    
    yticks(itick)=[yloc];
    %ytick_strs=[ytick_strs 'Observed'];    
end


if inc_obs_02==1
    yloc = yloc - spacing;
    
    hbar_size=0.02;
    obs_trend_02 = xscale*trend_dat_box_obs2{ibox,1}.coeffs(2);
    obs_CI_02 = xscale*trend_dat_box_obs2{ibox,1}.uncer_max;
    
    [herr,minx2,maxx2]=errorbarYY('horiz',obs_trend_02,yloc,obs_CI_02,gca,'k','o',2,hbar_size);
    hold on
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    plot(obs_trend_02,yloc,'ko','markerfacecolor','k','markersize',msize);
    hold on
    itick = itick + 1;
    ytick_strs{itick}=[obs_str2];
    
    obs_min = obs_trend_02 - obs_CI_02;
    obs_max = obs_trend_02 + obs_CI_02;
    
    yticks(itick)=[yloc];
    %ytick_strs=[ytick_strs 'Observed'];    
end




%% Model box and whisker plot of indvidual ensemble means
Nens = size(trend_dat_box_ens,3);

yloc = yloc - spacing*1.2;

Nbox_plots = Nbox_plots + 1;
clear dat
for iens=1:Nens   
    dat{Nbox_plots}(iens) = xscale * trend_dat_box_ens{ibox,itrend_box_whisker,iens}.coeffs(2); %trend values    
end
ypos_box(Nbox_plots) = yloc;
cols_box{Nbox_plots} = [0 0 1];

minx = min([minx; dat{Nbox_plots}(:)]);
maxx = max([maxx; dat{Nbox_plots}(:)]);


%Box plot of the data
% boxplot(dat,'orientation','horizontal','color','b');
% set(gca,'fontsize',18);
 me = mean(dat{Nbox_plots});
% hold on
% %Plot ensemble mean and CI range
plot(me,yloc,'bo','markerfacecolor','b','markersize',msize);

itick = itick + 1;
yticks(itick)=[yloc];
%ytick_strs={'UKESM ensemble'};
ytick_strs{itick}=[model_str ' ensemble'];

%% Plot trend from pre-averaged ensemble
yloc = yloc - spacing; 
plot(xscale*trend_dat_box{ibox,itrend_box_whisker}.coeffs(2),yloc,'bo','markerfacecolor','b','markersize',msize);
[herr,minx2,maxx2]=errorbarYY('horiz',xscale*trend_dat_box{ibox,itrend_box_whisker}.coeffs(2),yloc,xscale*trend_dat_box{ibox,itrend_box_whisker}.uncer_max,gca,'b','o',2,hbar_size);

minx=min(minx,minx2);
maxx=max(maxx,maxx2);


itick = itick + 1;
ytick_strs{itick} = [model_str  ' ensemble mean'];    
yticks(itick)=[yloc];



%% Plot largest model trend member
if iplot_extreme==1
    
yloc = yloc - spacing*1.2;
yloc_max=yloc+0.05;
[maxval,imax] = max(dat{Nbox_plots});
ens_trend_max = dat{Nbox_plots}(imax);
ens_max_CI = xscale*trend_dat_box_ens{ibox,1,imax}.uncer_max;

%% Plot smallest trend member
yloc_min=yloc-0.05;
[minval,imin] = min(dat{Nbox_plots});
ens_trend_min = dat{Nbox_plots}(imin);
ens_min_CI = xscale*trend_dat_box_ens{ibox,1,imin}.uncer_max;
%yloc2=mean([yloc_min yloc_max]);


    
    [herr,minx2,maxx2]=errorbarYY('horiz',ens_trend_max,yloc_max,ens_max_CI,gca,'b','o',2,hbar_size);
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    plot(ens_trend_max,yloc_max,'bo','markerfacecolor','b','markersize',msize);
    
    [herr,minx2,maxx2]=errorbarYY('horiz',ens_trend_min,yloc_min,ens_min_CI,gca,'b','o',2,hbar_size);
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    plot(ens_trend_min,yloc_min,'bo','markerfacecolor','b','markersize',msize);
    
   % ex_str='Extreme members';
        
    itick = itick + 1;
    ytick_strs{itick}=['Extreme members'];
    yticks(itick)=[yloc];
    
    
else    
    %ex_str = '';    
end


%% Nudged
if iadd_nudged==1    
    yloc = yloc - spacing;
    plot(xscale*trend_dat_box_nudged{ibox,1}.coeffs(2),yloc_nudged,'go','markerfacecolor','g','markersize',msize);
    [herr,minx2,maxx2]=errorbarYY('horiz',xscale*trend_dat_box_nudged{ibox,1}.coeffs(2),yloc_nudged,xscale*trend_dat_box_nudged{ibox,1}.uncer_max,gca,'g','o',2,hbar_size);   
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    
    
    itick = itick + 1;
    ytick_strs{itick}=['UKESM-nudged'];  
    yticks(itick)=[yloc];

end

%% HADGEM
if iadd_HADGEM==1                      
    
    % Box and whisker plot of indvidual ensemble means
    yloc = yloc - spacing;
    Nens = eval(['size(trend_dat_box_ens_' had_str ',3);']);
    
    Nbox_plots = Nbox_plots + 1;
    %clear dat2
    for iens=1:Nens
        dat{Nbox_plots}(iens) = xscale * eval(['trend_dat_box_ens_' had_str '{ibox,itrend_box_whisker,iens}.coeffs(2);']); %trend values
    end
    ypos_box(Nbox_plots) = yloc;
    cols_box{Nbox_plots} = cols{1};
    
    minx = min([minx; dat{Nbox_plots}(:)]);
    maxx = max([maxx; dat{Nbox_plots}(:)]);
    
    %plot the mean on top of the box and whisker
    me = mean(dat{Nbox_plots});
    hold on
    %Plot ensemble mean and CI range
    plot(me,yloc,'o','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
    
    itick = itick + 1;
    ytick_strs{itick}=['HadGEM3 ensemble'];
    yticks(itick)=[yloc];
    
    %Plot the bar for the ensemble mean    
    yloc = yloc - spacing;                        
    
    plot(xscale* eval(['trend_dat_box_' had_str '{ibox,itrend_box_whisker}.coeffs(2)']),yloc,'o','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
    [herr,minx2,maxx2]=errorbarYY('horiz',xscale*eval(['trend_dat_box_' had_str '{ibox,itrend_box_whisker}.coeffs(2)']),yloc,xscale*eval(['trend_dat_box_' had_str '{ibox,itrend_box_whisker}.uncer_max']),gca,cols{1},'o',2,hbar_size);
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    
    itick = itick + 1;
    ytick_strs{itick}=['HadGEM3 ensemble mean'];
    yticks(itick)=[yloc];  
    
end
%% Plot SW trends when have replaced dT with the observed dT using an estimate of dSW/dT from GHG-only run
% Using uncertainty from trend from pre-averaged ensemble

switch var_ukesm
    case 'rsut'
        if iadd_Tcorr_points==1
                        
            yloc = yloc - spacing;
            val = -3.6 /(2014-1985); %Value estimated in text of the paper
            %plot(xscale*trend_dat_box{ibox,itrend_box_whisker}.coeffs(2),yloc,'bo','markerfacecolor','b','markersize',msize);
            plot(xscale*val,yloc,'bs','markerfacecolor','b','markersize',msize);
            [herr,minx2,maxx2]=errorbarYY('horiz',xscale*val,yloc,xscale*trend_dat_box{ibox,itrend_box_whisker}.uncer_max,gca,'b','o',2,hbar_size);
            
            minx=min(minx,minx2);
            maxx=max(maxx,maxx2);
            
            
            itick = itick + 1;
            ytick_strs{itick} = [model_str  ' using \DeltaT_{obs}'];
            yticks(itick)=[yloc];
            
            if iadd_HADGEM==1
                yloc = yloc - spacing;
                val = -3.1 /(2014-1985); %Value estimated in text of the paper
                %plot(xscale*trend_dat_box{ibox,itrend_box_whisker}.coeffs(2),yloc,'bo','markerfacecolor','b','markersize',msize);
                plot(xscale*val,yloc,'s','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
                [herr,minx2,maxx2]=errorbarYY('horiz',xscale*val,yloc,xscale*trend_dat_box{ibox,itrend_box_whisker}.uncer_max,gca,cols{1},'o',2,hbar_size);
                
                minx=min(minx,minx2);
                maxx=max(maxx,maxx2);
                
                
                itick = itick + 1;
                ytick_strs{itick} = ['HADGEM using \DeltaT_{obs}'];
                yticks(itick)=[yloc];
            
            end
            
            
        end
        
end




%% AMIP
if iadd_amip==1    
    yloc = yloc - spacing;    
    plot(xscale*trend_dat_box_amip{ibox,itrend_box_whisker}.coeffs(2),yloc,'ro','markerfacecolor','r','markersize',msize);
    [herr,minx2,maxx2]=errorbarYY('horiz',xscale*trend_dat_box_amip{ibox,itrend_box_whisker}.coeffs(2),yloc,xscale*trend_dat_box_amip{ibox,itrend_box_whisker}.uncer_max,gca,'r','o',2,hbar_size);   
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    
    
    itick = itick + 1;
    ytick_strs{itick}=['UKESM-AMIP'];  
    yticks(itick)=[yloc];

end

%% DAMIP hist-GHG
if iadd_DAMIP==1                      
    
    % Box and whisker plot of indvidual ensemble means
    yloc = yloc - spacing;
    run_str = 'hist_GHG';
    Nens = eval(['size(trend_dat_box_ens_' run_str ',3);']);
    
    Nbox_plots = Nbox_plots + 1;
    %clear dat2
    for iens=1:Nens
        dat{Nbox_plots}(iens) = xscale * eval(['trend_dat_box_ens_' run_str '{ibox,itrend_box_whisker,iens}.coeffs(2);']); %trend values
    end
    ypos_box(Nbox_plots) = yloc;
    cols_box{Nbox_plots} = cols{1};
    
    minx = min([minx; dat{Nbox_plots}(:)]);
    maxx = max([maxx; dat{Nbox_plots}(:)]);
    
    %plot the mean on top of the box and whisker
    me = mean(dat{Nbox_plots});
    hold on
    %Plot ensemble mean and CI range
    plot(me,yloc,'o','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
    
    itick = itick + 1;
    ytick_strs{itick}=['DAMIP Hist-GHG ens'];
    yticks(itick)=[yloc];
    
    %Plot the bar for the ensemble mean    
    yloc = yloc - spacing;                        
    
    plot(xscale* eval(['trend_dat_box_' run_str '{ibox,itrend_box_whisker}.coeffs(2)']),yloc,'o','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
    [herr,minx2,maxx2]=errorbarYY('horiz',xscale*eval(['trend_dat_box_' run_str '{ibox,itrend_box_whisker}.coeffs(2)']),yloc,xscale*eval(['trend_dat_box_' run_str '{ibox,itrend_box_whisker}.uncer_max']),gca,cols{1},'o',2,hbar_size);
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    
    itick = itick + 1;
    ytick_strs{itick}=['DAMIP Hist-GHG ens mean'];
    yticks(itick)=[yloc];  
    
    
    % hist-aer too
    
     % Box and whisker plot of indvidual ensemble means
    yloc = yloc - spacing;
    run_str = 'hist_aer';
    Nens = eval(['size(trend_dat_box_ens_' run_str ',3);']);
    
    Nbox_plots = Nbox_plots + 1;
    %clear dat2
    for iens=1:Nens
        dat{Nbox_plots}(iens) = xscale * eval(['trend_dat_box_ens_' run_str '{ibox,itrend_box_whisker,iens}.coeffs(2);']); %trend values
    end
    ypos_box(Nbox_plots) = yloc;
    cols_box{Nbox_plots} = cols{1};
    
    minx = min([minx; dat{Nbox_plots}(:)]);
    maxx = max([maxx; dat{Nbox_plots}(:)]);
    
    %plot the mean on top of the box and whisker
    me = mean(dat{Nbox_plots});
    hold on
    %Plot ensemble mean and CI range
    plot(me,yloc,'o','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
    
    itick = itick + 1;
    ytick_strs{itick}=['DAMIP Hist-Aer ens'];
    yticks(itick)=[yloc];
    
    %Plot the bar for the ensemble mean    
    yloc = yloc - spacing;                        
    
    plot(xscale* eval(['trend_dat_box_' run_str '{ibox,itrend_box_whisker}.coeffs(2)']),yloc,'o','color',cols{1},'markerfacecolor',cols{1},'markersize',msize);
    [herr,minx2,maxx2]=errorbarYY('horiz',xscale*eval(['trend_dat_box_' run_str '{ibox,itrend_box_whisker}.coeffs(2)']),yloc,xscale*eval(['trend_dat_box_' run_str '{ibox,itrend_box_whisker}.uncer_max']),gca,cols{1},'o',2,hbar_size);
    minx = min(minx,minx2); maxx = max(maxx,maxx2);
    
    itick = itick + 1;
    ytick_strs{itick}=['DAMIP Hist-Aer ens mean'];
    yticks(itick)=[yloc];  
    
end


%% Plot the mulitple box plots at once since can't seem to do them individually
L_max=-1;
for i=1:length(dat)
   L_max = max(L_max,length(dat{i}));         
end

dat2 = NaN * ones([L_max,Nbox_plots]);
clear col_matrix
for i=1:length(dat)
    dat2(1:length(dat{i}),i) = dat{i};
    col_matrix(i,:) = cols_box{i};
end
col_matrix = flipdim(col_matrix,1); %for some reason the order of the colours is reversed??
%Probably because the positions in ypos_box run from high to low
boxplot(dat2,'colors',col_matrix,'Positions',ypos_box,'orientation','horizontal');

%% Tidy plot

set(gca,'fontsize',18);
xrange = maxx-minx;
set(gca,'xlim',[minx-xrange*0.1 maxx+xrange*0.1]);

if inc_obs_02==1
    %obs_str='Observed';
    obs_str=obs_str2;
else
    obs_str = '';
end

[yticks2,b]=sort(yticks);
ytick_strs2 = ytick_strs(b);



% if iadd_amip==1
    set(gca,'ytick',yticks2);
    %set(gca,'yticklabel',{'UKESM-AMIP' 'Trend from ensemble mean' 'UKESM ensemble',ex_str,obs_str});
    %set(gca,'yticklabel',ytick_strs2,'FontName','symbol');
    
%    set(gca,'yticklabel',ytick_strs2); %Doesn't work with latex so do this :-
    yts=get(gca,'ytick');
    xlims = get(gca,'xlim');
    xrange2 = diff(xlims);
    xval = xlims(1) - xrange2*0.01;
    
    set(gca,'yticklabel',{[]});
    for i=1:length(yts); 
        text(xval,yts(i),ytick_strs2{i},'HorizontalAlignment','right','fontsize',18); 
    end

    %set(gca,'ylim',[0.65 yloc+0.1]);
    %set(gca,'ylim',[0.45 yloc+0.1]);
% else
%     set(gca,'ytick',[yloc_mean 1 yloc2 yloc]);
%     set(gca,'yticklabel',{'Trend from ensemble mean' 'UKESM ensemble',ex_str,obs_str});
%     set(gca,'ylim',[0.75 yloc+0.1]);
% end

%yrange = range(yticks);
set(gca,'ylim',[min(yticks)-spacing max(yticks)+spacing]);
yrange = get(gca,'ylim');
xlabel({' ' xlab}); %Make a gap between x-axis and xlabel
title({['Region=' box_region ' ' land_ocean],'',[box_region_str ' ' season]});


set(gca,'position',[0.3    0.2057    0.5    0.5833]);


