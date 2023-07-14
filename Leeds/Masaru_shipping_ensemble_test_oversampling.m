
var_name = 'rsut'; var_str = 'SW TOA';
var_name = 'rlut'; var_str = 'LW TOA';

% switch var_name
%     case 'rsut'
%         data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/30yr_avs/';
%         data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/rsut/';
%     case 'rlut'
%         data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/rsut/';
% end

data_dir = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/' var_name '/'];

%PI data runs from year 1960 to 3839
%(3839-1960+1)/30 = 62.667

%N = 62;
Nyears = 30;
N = 3839-1960+1 - Nyears + 1; %Need to allow a 30-year average at the end
dN = 1/N; %interval of 62 bins between 0 and 1
start_year = 1960;

%choose 2 sets of 11 members randomly from the 62 * 30-year samples

Nmembers=11;

lat_file ='/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';
load(lat_file,'gcm_Plat2D_UM','gcm_Plon2D_UM');

LATS = [-50 50]; LONS = [140 270];
gcm_Plon2D_UM_360 = gcm_Plon2D_UM;
ilon180 = find(gcm_Plon2D_UM<0);
gcm_Plon2D_UM_360(ilon180) = gcm_Plon2D_UM_360(ilon180) + 360;
ilat = find(gcm_Plat2D_UM(:,1)>LATS(1) & gcm_Plat2D_UM(:,1)<LATS(2));
ilon = find(gcm_Plon2D_UM_360(1,:)>LONS(1) & gcm_Plon2D_UM_360(1,:)<LONS(2));


clear pi_ens_diff_std
for isample=1:1000 %run random ens selection 1000 times.

    fprintf(1,'%d ',isample);

for i=1:2
   %Does it matter if we choose the same one twice? 
   %rand chooses a random number between 0 and 1   
   clear pi_ens_dat_temp
   for j=1:Nmembers
       n = min( floor(rand/dN) + 1 , N);
       n = max(1,n); %just in case it gets set to zero
       %year = start_year + (n-1)*30; %year for the file
       year = start_year + (n-1); %year for the file
       sum_dat = 0;
       for iyear=year:year+Nyears-1
           %filename = [data_dir 'rsut_UKESM1-0-LL_piControl_r1i1p1f2_30yr_av_' num2str(year) '.nc.nc3'];
           filename = [data_dir var_name '_UKESM1-0-LL_piControl_r1i1p1f2_ANNUAL_av_' num2str(iyear) '.nc.nc3'];
           nc=netcdf(filename);
           sum_dat = sum_dat + nc{var_name}(:);
           nc=close(nc);
       end
       pi_ens_dat(i).ens_dat(:,:,j) = sum_dat / Nyears; %Average over Nyears
       pi_ens_dat(i).start_year(j) = iyear; %Average over Nyears
   end   
   [pi_ens_dat(i).ens_mean, nums, pi_ens_dat(i).ens_std] = meanNoNan(pi_ens_dat(i).ens_dat,3); %ensemble mean
end

pi_ens_diff = pi_ens_dat(1).ens_mean - pi_ens_dat(2).ens_mean;

dat = pi_ens_diff(ilat,ilon);

[me,nnum,pi_ens_diff_std(isample)] = meanNoNan(dat(:),1); %
pi_ens_diff_range(isample) = range(dat(:)); %

end

save_name = [data_dir var_name '_random_1000_sample_stds.mat'];
save(save_name,'-V7.3','pi_ens_diff_std','pi_ens_diff_range','LATS','LONS');

min_std = min(pi_ens_diff_std,[],2);
max_std = max(pi_ens_diff_std,[],2);
[me_std,nnum,std_dev] = meanNoNan(pi_ens_diff_std,2);
median_std = prctile(pi_ens_diff_std,50);
fprintf(1,'\n min=%f max=%f mean=%f median=%f\n',min_std,max_std,me_std,median_std);


min_range = min(pi_ens_diff_range,[],2);
max_range = max(pi_ens_diff_range,[],2);
[me_range,nnum,std_dev] = meanNoNan(pi_ens_diff_range,2);
median_range = prctile(pi_ens_diff_range,50);
fprintf(1,'\n min=%f max=%f mean=%f median=%f\n',min_range,max_range,me_range,median_range);

%% t-test

dat01 = permute(pi_ens_dat(1).ens_dat,[3 1 2]); %pi_ens_dat(1).ens_dat is of size [144x192x11 double]
dat02 = permute(pi_ens_dat(2).ens_dat,[3 1 2]); %are permuting here to put the 11 sized dim first
[h,p] = ttest(dat01(:,:),dat02(:,:)); %then this function works to give the h-value (1 if is signficantly different at 95% and 0 if not)
% and the p-value (lower values are more significantly different).
h = reshape(h,size(pi_ens_dat(1).ens_dat(:,:,1))); %re-shape to global map
p = reshape(p,size(pi_ens_dat(1).ens_dat(:,:,1)));

%% Plot on map

dat_modis = pi_ens_diff;
var_UM = var_name;
%tit_str_clean='% f_c change';
subtitle_str=['Ensemble ' var_str ' difference (W m^{-2})'];
%dat_modis(dat_modis>999)=999; %to deal with gridpoints with Inf values from where f0_mean=0

irestrict_domain_DRIVER=0;
ioverride_proj_type=0;
proj_type_DRIVER='ortho';
ioverride_LAT_plots=0;
iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
ioverride_ticks_DRIVER=1;
icoarse_grain=0;
icontour_DRIVER=0; cont_col_str_DRIVER='';
time_round=0; time_format_str='';
isave_plot=0;
iplot_wind_arrows=0;
 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-0.1 0.1]);
caxis([-10 10]);
caxis([-1.5 1.5]);

idiff_sig = find(p<0.05); %need ot search for the low p values for significance
idiff_sig = find(p<0.1); %need ot search for the low p values for significance
%idiff_sig = find(h==1);

m_plot(gcm_Plon2D_UM_360(idiff_sig),gcm_Plat2D_UM(idiff_sig),'ko','markersize',1,'markerfacecolor','k'); %m_plot works using lon,lat



%% Plot on map
dat_modis = p;
var_UM = var_str;
%tit_str_clean='% f_c change';
subtitle_str='p-value';
%dat_modis(dat_modis>999)=999; %to deal with gridpoints with Inf values from where f0_mean=0

irestrict_domain_DRIVER=0;
 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-0.1 0.1]);
caxis([0 1]);

%% PDF of SW diffs
Ybins_DRIVER = [-15:0.1:15];
Y_driver = pi_ens_diff;
%Y_driver = dat;
ylabelstr=[var_str ' difference between ensembles (W m^{-2})'];
%run script
Masaru_shipping_PDF

%% PDF of std dev of SW diffs for Pacific
Ybins_DRIVER = [-1:0.01:1];
Y_driver = pi_ens_diff_std;
%Y_driver = dat;
ylabelstr=['Std. dev. of ' var_str ' difference in Pacific for samples (W m^{-2})'];
%run script
Masaru_shipping_PDF

