
var_name = 'rsut'; var_str = 'SW TOA';
%var_name = 'rlut'; var_str = 'LW TOA';

% switch var_name
%     case 'rsut'
%         data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/30yr_avs/';
%         data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/rsut/';
%     case 'rlut'
%         data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/rsut/';
% end

%data_dir = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/' var_name '/'];
data_dir = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/' var_name '/'];

load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsut.mat';
mat_obj = matfile(load_file);

load(load_file,'gcm_Plat2D_UM','gcm_Plon2D_UM','gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM','years_ukesm_1d');

start_year = 2014-30+1;
Nyears=30;
istart = find(years_ukesm_1d==start_year);
iend = find(years_ukesm_1d==start_year+Nyears-1);

LATS = [-50 50]; LONS = [140 270];
gcm_Plon2D_UM_360 = gcm_Plon2D_UM;
ilon180 = find(gcm_Plon2D_UM<0);
gcm_Plon2D_UM_360(ilon180) = gcm_Plon2D_UM_360(ilon180) + 360;
ilat = find(gcm_Plat2D_UM(:,1)>LATS(1) & gcm_Plat2D_UM(:,1)<LATS(2));
ilon = find(gcm_Plon2D_UM_360(1,:)>LONS(1) & gcm_Plon2D_UM_360(1,:)<LONS(2));

Nens = size(mat_obj,'dat_ens',1);

dat_in = mat_obj.dat_annual_ens(:,istart:iend,:,:);

clear pi_ens_diff_std pi_ens_diff_range
for isample=1:100 %run random ens selection N times.
    
    fprintf(1,'%d ',isample);
    
    
    %split the ensemble into 2 halves, assigned randomly
    ileft=[1:Nens];
    Nleft = length(ileft);
    i=0;
    clear groups
    while Nleft>0
        i=i+1;
        iens = max(ceil(rand*Nleft),1); %randomly assign between 1 and Nleft
        groups{1}(i) = ileft(iens);
        ileft(iens)=[];
        Nleft = length(ileft);
        if Nleft>0
            iens = max(ceil(rand*Nleft),1); %randomly assign between 1 and Nleft
            groups{2}(i) = ileft(iens);
            ileft(iens)=[];
            Nleft = length(ileft);
        end
    end
    
    
    
    
    
    
    for i=1:2
        Nmembers = length(groups{i});
        clear pi_ens_dat_temp
        for j=1:Nmembers
            iens = groups{i}(j);
            pi_ens_dat(i).ens_dat(:,:,j) = meanNoNan(dat_in(iens,:,:,:),2); %av over Nyears for each ensemble member
        end
        [pi_ens_dat(i).ens_mean, nums, pi_ens_dat(i).ens_std] = meanNoNan(pi_ens_dat(i).ens_dat,3); %ensemble mean
    end
    
    pi_ens_diff = pi_ens_dat(1).ens_mean - pi_ens_dat(2).ens_mean;
    
    dat = pi_ens_diff(ilat,ilon); %restrict to Pacific region
    
    [me,nnum,pi_ens_diff_std(isample)] = meanNoNan(dat(:),1); %
    pi_ens_diff_range(isample) = range(dat(:)); %
    
end

[me,nnum,std_dev] = meanNoNan(pi_ens_diff_range,2);
min_val = min(pi_ens_diff_range,[],2)
max_val = max(pi_ens_diff_range,[],2)

save_name = [data_dir var_name '_Masaru_random_samples_UKESM_historical_stds.mat'];
save(save_name,'-V7.3','pi_ens_diff_std','pi_ens_diff_range');

%% t-test

dat01 = permute(pi_ens_dat(1).ens_dat,[3 1 2]); %pi_ens_dat(1).ens_dat is of size [144x192x11 double]
dat02 = permute(pi_ens_dat(2).ens_dat,[3 1 2]); %are permuting here to put the 11 sized dim first
[h,p] = ttest(dat01(:,:),dat02(:,:)); %then this function works to give the h-value (1 if is signficantly different at 95% and 0 if not)
% and the p-value (lower values are more significantly different).
h = reshape(h,size(pi_ens_dat(1).ens_dat(:,:,1))); %re-shape to global map
p = reshape(p,size(pi_ens_dat(1).ens_dat(:,:,1)));

%% Plot on map

dat_modis = pi_ens_diff;
var_UM = 'rsut';
%tit_str_clean='% f_c change';
subtitle_str=['Ensemble ' var_str ' outgoing difference (W m^{-2})'];
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



%% Plot on p-value map
dat_modis = p;
var_UM = 'rsut';
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

