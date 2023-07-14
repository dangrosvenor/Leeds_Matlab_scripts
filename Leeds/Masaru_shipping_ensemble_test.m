

data_dir = '/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/r1i1p1f2/30yr_avs/';
var_name = 'rsut';
%PI data runs from year 1960 to 3839
%(3839-1960+1)/30 = 62.667

N = 62;
dN = 1/N; %interval of 62 bins between 0 and 1
start_year = 1960;

%choose 2 sets of 11 members randomly from the 62 * 30-year samples

Nmembers=11;


for i=1:2
   %Does it matter if we choose the same one twice? 
   %rand chooses a random number between 0 and 1
   sum_dat = 0;
   for j=1:Nmembers
       n = min( floor(rand/dN) + 1 , N)
       year = start_year + (n-1)*30; %year for the file
       filename = [data_dir 'rsut_UKESM1-0-LL_piControl_r1i1p1f2_30yr_av_' num2str(year) '.nc.nc3']; 
       nc=netcdf(filename); 
       sum_dat = sum_dat + nc{var_name}(:);
       nc=close(nc);
       
   end
   
   ens_av{i} = sum_dat / Nmembers;         
end

%% Plot on map

dat_modis = ens_av{1} - ens_av{2};
var_UM = 'rsut';
%tit_str_clean='% f_c change';
subtitle_str='Ensemble SWTOA outgoing difference (W m^{-2})';
%dat_modis(dat_modis>999)=999; %to deal with gridpoints with Inf values from where f0_mean=0

irestrict_domain_DRIVER=0;
 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-0.1 0.1]);
caxis([-10 10]);

