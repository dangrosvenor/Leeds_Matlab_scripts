%% Minghuai Nd trends

iscreen_sig=1; %Whether to screen for signficance
marker_size=1; %For non-signficant points
iplot_mgrid_lines_DRIVER=0; %whether to plot the grid lines for maps using m_grid

p_conf = 95; % Confidence limit (%) for the trend significance
nthresh_days = 3;
%nthresh_days = 0;

LAT_val = [30 42]; LON_val =[-68 -52]; %
LAT_val = [30 40]; LON_val =[-73 -63]; %box in plots sent to Ken
%LAT_val = [30 45]; LON_val =[-73 -40]; %wider box region in NW Atlantic
LAT_val = [30 40]; LON_val =[-20 -10]; %Region SW of Spain down African coast

UKESM_Nd_case = 'to ztop';
UKESM_Nd_case = 'to 3.2km';

savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
eval(['!mkdir ' savedir_date]); 

icoarse_grain=0;

time_round='';
time_format_str='';
icontour_DRIVER=0;
isave_plot=0;
iplot_wind_arrows=0;

cont_col_str_DRIVER='k';

region_choice = 'Southern NA';
%region_choice = 'Northern NA';
region_choice = 'Rosenfeld VOCALS';
region_choice = 'Rosenfeld ALL';
%region_choice = 'VOCALS CPT';
%region_choice = 'VOCALS coastal';
region_choice = 'Northern Hemisphere Yawen';
region_choice = 'Northern Hemisphere Yawen, further north';
region_choice = 'North Atlantic further north';
%region_choice = 'global';

irestrict_domain_DRIVER=1;
switch region_choice
    case 'global'
        irestrict_domain_DRIVER=0;
end

% Run the script to pick the region
[LAT_val_DRIVER2, LON_val_DRIVER2, region_shortname] = UM_ACSIS_choose_region(region_choice);

LAT_val_DRIVER_override = LAT_val_DRIVER2;
LON_val_DRIVER_override = LON_val_DRIVER2;

% var_UM = 'Nd_cf_weighted_UKESM_ztop';
% um_case='UKESM/r1i1p1f2_u-bc179/'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
% dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
% dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);


%% Run the script to load the 9 member ensemble
%ACSIS_Nd_trends_load_ensemble(); %switch to loading the .mat file

switch UKESM_Nd_case
    case 'to ztop'
        load('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat');
    case 'to 3.2km'
        load('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm_Nd_cf_weighted_UKESM.mat');
end




%% plots
ioverride_LAT_plots=1;


istart=1; iend=1;
    dat_modis = squeeze (Nd_annual(istart,:,:) );     
    var_UM = ['Nd for y' num2str(years_ukesm_1d(istart)) '; (cm^{-3})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-6 6]);   
    %caxis([-2 2]);

    
    
 %y2015 minus y2000
    istart=1; iend=15;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart); 
    var_UM = ['Nd trend between y' num2str(years_ukesm_1d(istart)) ' and y' num2str(years_ukesm_1d(iend)) '; (cm^{-3} yr^{-1})'];    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);
    %caxis([-2 2]);
        
    istart=2; iend=15;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart);     
    var_UM = ['Nd trend between y' num2str(years_ukesm_1d(istart)) ' and y' num2str(years_ukesm_1d(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);
    caxis([-2 2]);
            
    
    
    istart=1; iend=14;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart);     
    var_UM = ['Nd trend between y' num2str(years_ukesm_1d(istart)) ' and y' num2str(years_ukesm_1d(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);   
    caxis([-2 2]);

%% Trend analysis (linear least squares fit) - model
    yr_start=2003; yr_end=2014;
    istart=find(years_ukesm_1d==yr_start);
    iend=find(years_ukesm_1d==yr_end);    
    x = [yr_start:yr_end]';
%     y = Nd_annual(istart:iend,100,100);
%     [Nd_trend,bint,residuals,rint,stats] = regress(y,x);
%     
%     [coeffs,bint,residuals,rint,stats] = trend_linear_fit_nd_data(Nd_annual(istart:iend,:,:),1);
    
   [coeffs,t_trend] = trend_linear_fit_nd_data(x,Nd_annual(istart:iend,:,:),1); 
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    %t_thresh = tinv(p_conf/100,n_dof); %find the t value needed for 95% confidence using a one-tailed t-test
        %N.B. - here is the function for a 2-tailed test :-
            % E.g. 
            % t=4; v=10;
            % tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));   % 2-tailed t-distribution function
            % tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2; %1-tailed t-distribution function
            % OR just :-  
            % tail2P = 2*tcdf(-abs(t),v);
            % tail1P = tcdf(-abs(t),v);
            % Test :-            
            % T2 = [1-tdist2T(t,v)  tail2P]; %give the same answer
            % T1 = [1-tdist1T(t,v)  tail1P]; %give the same answer
    %itrend_not_sig=find(abs(t_trend)<t_thresh); 
    
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend),n_dof);
    T1 = T2/2;    
    itrend_not_sig = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig = find((1-T1)<=p_conf/100); %1-tailed
   
%% Map of 2003 to 2014 linear MODEL trend    
    dat_modis = squeeze(coeffs(2,:,:)); 
    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis(itrend_not_sig)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    var_UM = ['UKESM Nd trend of ensemble mean between y' num2str(years_ukesm_1d(istart)) ' and y' num2str(years_ukesm_1d(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);   
    %caxis([-2 2]);
   
    
    % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    plot_box_on_map
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_UM(itrend_not_sig),gcm_Plat2D_UM(itrend_not_sig),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end
    
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
    
%% Map of 2003 to 2014 difference
iplot_diff=0;
if iplot_diff==1
    
    istart=4; iend=15;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart);     
    var_UM = ['UKESM Nd trend of ensemble mean between y' num2str(years_ukesm_1d(istart)) ' and y' num2str(years_ukesm_1d(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);   
    %caxis([-2 2]);
    
    
% Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    plot_box_on_map
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
end

    
%% calculate some timeseries of Nd in the box
    ilat = find(gcm_Plat2D_UM(:,1)>LAT_val(1) & gcm_Plat2D_UM(:,1)<LAT_val(2));
    ilon = find(gcm_Plon2D_UM(1,:)>LON_val(1) & gcm_Plon2D_UM(1,:)<LON_val(2));
    clear Nd_time_mean
    for iy=1:size(Nd_ukesm,1)
        Nd_time_mean(iy,:) = meanNoNan(meanNoNan(Nd_ukesm(iy,:,ilat,ilon),4),2);        
    end
    
  
    clear Nd_annual_box_ukesm Nd_annual_box_ukesm_ens
    for it=1:size(Nd_annual,1)
        Nd_annual_box_ukesm(it) = meanNoNan(meanNoNan(Nd_annual(it,ilat,ilon),3),2);
        dat = meanNoNan(meanNoNan(Nd_annual_ens(:,it,ilat,ilon),4),2);
        Nd_annual_box_ukesm_ens_std(it) = std(dat); %std dev across the ensemble
    end
    

    
    
%% MODIS data
% Using the data given to Jane - screened for sea-ice etc.
cf_screen_str = 'CF>80';
cf_screen_str = 'CF>0';

str_2137='21';
str_2137='37';

res_str='1deg';
res_str='1km';


switch cf_screen_str
    case 'CF>80'
         file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/';
         dataset_str = 'SZA_LT_65_CF_GT_80_CTH_LT_3.2km_screened_for_seaice__2week_max';
         dataset_str = 'SZA_LT_65_CF_GT_80_CTH_LT_3.2km';
        
    case 'CF>0'        
        file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.0_meanCTT_173_meanCTH_3.2km_SZA_65/';
        dataset_str = 'SZA_LT_65_CF_GT_0_CTH_LT_3.2km';
end




switch str_2137
    case '21'
        str_label_2137='2.1 \mum';
    case '37'
        str_label_2137='3.7 \mum';
end

years_MODIS2=[2003:2014];
clear mon_me_MODIS2 mon_me_MODIS2_Ndatap
for iy=1:length(years_MODIS2)
    year_str = num2str(years_MODIS2(iy));
    filename = [file_dir 'Nd_monthly_' str_2137 '_' res_str '_' year_str '_' dataset_str '.mat.nc'];
    nc=netcdf(filename);
    if iy==1
        lat=nc{'lat'}(:);
        lon=nc{'lon'}(:);
        [gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE]=meshgrid(lon,lat);
    end

    mon_me_MODIS2{iy} = nc{['Nd_' res_str '_mean']}(:);
    mon_me_MODIS2_Ndatap{iy} = nc{['Nd_' res_str '_Ndatap']}(:);
    inan = find(mon_me_MODIS2_Ndatap{iy} < nthresh_days);
    mon_me_MODIS2{iy}(inan)=NaN;
    %mon_me_filter{iy}(inan)=NaN;
end

[mon_me_region_MODIS2] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,LAT_val,LON_val);
[mon_me_region_MODIS2_Ndatap] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2_Ndatap,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,LAT_val,LON_val);
%[mon_me_filter_SO] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_filter,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,thresh_LAT,thresh_LON);

for iy=1:length(years_MODIS2)
    Nd_MODIS(iy,:) = mon_me_region_MODIS2{iy}(:); %order is [year month]
end
 

clear Nd_annual_box_MODIS Nd_annual_MODIS
for iy=1:size(Nd_MODIS,1)
    Nd_annual_box_MODIS(iy) = meanNoNan(Nd_MODIS(iy,:),2);   
    Nd_annual_MODIS(iy,:,:) = meanNoNan(mon_me_MODIS2{iy},3);
end

%% Trend analysis (linear least squares fit) - MODIS
    yr_start=2003; yr_end=2014;
    istart=find(years_MODIS2==yr_start);
    iend=find(years_MODIS2==yr_end);                
    x = [yr_start:yr_end]';
    %The following takes a litte while...
   [coeffs_MODIS,t_trend_MODIS] = trend_linear_fit_nd_data(x,Nd_annual_MODIS(istart:iend,:,:),1); 
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    %t_thresh = tinv(p_conf/100,n_dof); %find the t value needed for 95% confidence using a one-tailed t-test
        %N.B. - here is the function for a 2-tailed test :-
            % E.g. 
            % t=4; v=10;
            % tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));   % 2-tailed t-distribution function
            % tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2; %1-tailed t-distribution function
            % OR just :-  
            % tail2P = 2*tcdf(-abs(t),v);
            % tail1P = tcdf(-abs(t),v);
            % Test :-            
            % T2 = [1-tdist2T(t,v)  tail2P]; %give the same answer
            % T1 = [1-tdist1T(t,v)  tail1P]; %give the same answer
    %itrend_not_sig_MODIS = find(abs(t_trend_MODIS)<=t_thresh); %make NaN for now, but can put a dot on, etc.
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend_MODIS),n_dof);
    T1 = T2/2;    
    itrend_not_sig_MODIS = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig_MODIS = find((1-T1)<=p_conf/100); %1-tailed
    
    %% Plot MODIS fitted trend
    dat_modis2 = squeeze(coeffs_MODIS(2,:,:));
    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd trend between y' num2str(years_MODIS2(istart)) ' and y' num2str(years_MODIS2(iend)) '; (cm^{-3} yr^{-1})' add_str];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);
    
        % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    plot_box_on_map
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_AMSRE(itrend_not_sig_MODIS),gcm_Plat2D_AMSRE(itrend_not_sig_MODIS),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end
    
    
 %% Plot MODIS Nd for 2006
    istart=find(years_MODIS2==2006);
    Nd_annual_MODIS(istart,:,:);
    dat_modis2 = squeeze(Nd_annual_MODIS(istart,:,:));
    %dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
    
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd for y' num2str(years_MODIS2(istart)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 300]);



%% plot MODIS difference map 2003 to 2014
%2003 to 2014    
    yr_start=2003; yr_end=2014;
    istart=find(years_MODIS2==yr_start);
    iend=find(years_MODIS2==yr_end);
    dat = squeeze (Nd_annual_MODIS(iend,:,:) - Nd_annual_MODIS(istart,:,:) ) / (iend-istart);   
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);
    var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd trend between y' num2str(years_MODIS2(istart)) ' and y' num2str(years_MODIS2(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);   
    %caxis([-2 2]);
    
        
% Plot the box    
    ioverride_box_colour=1;    
    col_str='k-';
    box_lwidth = 3;
    plot_box_on_map
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%% Linear trend calc for the annual means over the box region.
yr_start=2003; yr_end=2014;
istart=find(years_ukesm_1d==yr_start);
iend=find(years_ukesm==yr_end);
y = Nd_annual_box_ukesm(istart:iend);
x1 = [yr_start:yr_end];
x = [ones(size(y)); x1];
[Nd_trend,bint,residuals,rint,stats] = regress(y',x');
ylin=Nd_trend(1)+Nd_trend(2).*x1; %this is the linear trend fit line

    
%% plot annual mean timeseries
titlenam = 'Annual mean N_d';
figure
clear leg_str
ileg=1;
h1=errorbar(years_ukesm_1d,Nd_annual_box_ukesm,Nd_annual_box_ukesm_ens_std*2,'bo-');
leg_str{ileg}='UKESM'; ileg=ileg+1;
set(h1,'linewidth',3);
set(h1,'markerfacecolor','b');
hold on


h2=plot(years_MODIS2,Nd_annual_box_MODIS,'rs-');
leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
set(h2,'linewidth',3);
set(h2,'markerfacecolor','r');
loc='SouthWest';
loc='NorthEast';

set(gca,'xlim',[1999 2015]);
set(gca,'XMinorTick','off');
fontsize_figure(gcf,gca,18);

legend(leg_str,'location',loc);
xlabel('Year');
ylabel('N_d (cm^{-3})');
title(titlenam);
grid on

savename=[savedir_date titlenam ' ' str_label_2137 ' ' cf_screen_str];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);



%% Find a correction for the MODIS data to adjust for the fact that we are using y2000 emissions, but comparing to y2009.
% Will apply the correction to tbe obs - so need to take the obs back to y2000
% So, if the trend is positive then need to reduce the obs and so dNd will
% be negative
MODIS_trends = squeeze(coeffs_MODIS(2,:,:));
MODIS_trends(itrend_not_sig_MODIS)=0; %make the non-significant trends zero
ny_corr = 2009-2000; %number of years that are correcting for.
dNd_corr = -MODIS_trends*ny_corr;
dNd_corr_UM_grid = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dNd_corr,gcm_Plat2D_UM,gcm_Plon2D_UM); %ADD theses to MODIS values

save_file_dNd_trend = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/dNd_trend_correction.mat'];
save(save_file_dNd_trend,'dNd_corr','dNd_corr_UM_grid','-V7.3');
%ref year for percentage diff plot
yr=2009;
%yr=2010;
iyr=find(years_MODIS2==yr);
Nd_ref_year = squeeze (Nd_annual_MODIS(iyr,:,:));

prc_diff = 100*dNd_corr./Nd_ref_year;


dat_modis = dNd_corr_UM_grid;
%dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,prc_diff,gcm_Plat2D_UM,gcm_Plon2D_UM);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-60 60]);
    
%%
dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,Nd_ref_year,gcm_Plat2D_UM,gcm_Plon2D_UM);
var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd for y' num2str(years_MODIS2(iyr)) ' (cm^{-3})'];
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 300]);
