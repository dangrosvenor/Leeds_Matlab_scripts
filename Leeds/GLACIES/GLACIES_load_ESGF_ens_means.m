%function []=ACSIS_dat_trends_load_ensemble_esgf()
%Run from ACSIS_dat_trends_load_ESGF_ensemble_multi_vars.m

%output_period = 'all';
%output_period = 'recent';
%output_period = 'SW_down_TOA_partial';

%i_annual_data=0; %(setting default). Flag to say that the data is already annual averages (not monthly).

%Loading and scaling options
opts.isort_time=1;
iuse_mat_file=0;
fscale_fac = 1;
opts.lat_var = 'lat';
opts.lon_var = 'lon';
    
i_single_ens=1;

%Directory to process :-
dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/'];

%Directory to save to
%savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' MIP '_' expt '_' output_period ens_str];

%load_type = 'mat';
%load_type = 'merged netCDF';
%load_type = 'individual netCDF files'; opts.cat_dim=1;
load_type =  'named netCDF';       

var_UM = 'cli';        


%savefile = [savefile_pre_str '_' var_UM '.mat'];

dirUM = [dir_data var_UM '_all_ens/'];
%error('Need to set something here!')


opts.named_file = [var_UM '_ens_mean.nc3'];
opts.named_dir = [dirUM];

var_UM_load = var_UM;

opts.time_var='time';
opts.time_ref = datenum('01-Jan-1850');
opts.time_fconv = 1; %conversion multiplier to get to days
pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
dat_global = UM_load_merged_netCDF(dirUM,var_UM_load,run_type,load_type,[],[],var_UM_load,opts); %data is ordered [time lat lon]. 180 times (monthly over 15 years)

%Load the height data for ESGF data
nc = netcdf([dirUM var_UM '_ens_mean.nc3']);
lev = nc{'lev'}(:);
b = nc{'b'}(:);
orog = nc{'orog'}(:);
lev_bnds = nc{'lev_bnds'}(:);
b_bnds = nc{'b_bnds'}(:);

nc=close(nc);

b3d = repmat(b,[1 size(orog,1) size(orog,2)]);
b3d = permute(b3d,[2 3 1]);
lev3d = repmat(lev,[1 size(orog,1) size(orog,2)]);
lev3d = permute(lev3d,[2 3 1]);
orog3d = repmat(orog,[1 1 size(b,1)]);
altitude = lev3d + b3d.*orog3d; %according to the NetCDF "z = lev + b*orog" 

b_edges = [b_bnds(:,1); b_bnds(end,2)];
lev_edges = [lev_bnds(:,1); lev_bnds(end,2)];

b3d_edges = repmat(b_edges,[1 size(orog,1) size(orog,2)]);
b3d_edges = permute(b3d_edges,[2 3 1]);
lev3d_edges = repmat(lev_edges,[1 size(orog,1) size(orog,2)]);
lev3d_edges = permute(lev3d_edges,[2 3 1]);
orog_edges = griddata(
orog3d_edges = repmat(orog,[1 1 size(b_edges,1)]);
altitude_edges= lev3d_edges + b3d_edges.*orog3d_edges; %according to the NetCDF "z = lev + b*orog" 


% nc = netcdf([dir_data 'ta_all_ens/ta_ens_mean.nc3']);
% ta = nc{'ta'}(:);
% nc=close(nc);
% ta(ta>1e19)=NaN;
% 
% nc = netcdf([dir_data 'pfull_all_ens/pfull_ens_mean.nc3']);
% pfull = nc{'pfull'}(:);
% nc=close(nc);
% pfull(pfull>1e19)=NaN;
%  
% theta = potemp(ta,pfull); %Can't do this at the moment because the air temp is only available over lower numbers of levels.
    %Could get the airmass variable, which is the air mass for each model
    %level. Then can calculate the mean air density given the layer depths.
    %Then with the pressure and rho can calculate the temperature.
 
    
%% Calculate lats and times, etc.
gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
%Had this the wrong way around :-
%[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
[gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

start_year=2000; end_year=2014;
years_ukesm_1d = [start_year:end_year];
years_ukesm = repmat(years_ukesm_1d,[12 1]);
years_ukesm = years_ukesm(:);
    
%% Calculate zonal mean ice MMRs
years_to_plot = [2006 2014];
iy=find(years_ukesm>=years_to_plot(1) & years_ukesm<=years_to_plot(2));

cli_time_av = meanNoNan(dat_global.dat(iy,:,:,:),1); %[85 144 192]
cli_time_av = permute(cli_time_av,[2 3 1]);
lat_3d = repmat(gcm_Plat2D_UM,[1 1 size(altitude,3)]);
lon_3d = repmat(gcm_Plon2D_UM,[1 1 size(altitude,3)]);
%interpolate to a regular height abv msl grid.
dz=200;
z_new=[0:dz:20e3];
lat_3d_new = repmat(gcm_Plat2D_UM,[1 1 length(z_new)]); %[144 192 z]
lon_3d_new = repmat(gcm_Plon2D_UM,[1 1 length(z_new)]);
z_3d_new = repmat(z_new',[1 size(gcm_Plat2D_UM)]); %[z 144 192]
z_3d_new = permute(z_3d_new,[2 3 1]);

cli = griddata(altitude,lat_3d,lon_3d,cli_time_av,z_3d_new,lat_3d_new,lon_3d_new);

cli_zonal = meanNoNan(dat_global.dat,4);
alt_zonal = meanNoNan(altitude_edges,2);   alt_zonal = alt_zonal';  
%add an extra profile of alt_zonal to make it the same size as lat_edges
alt_zonal(:,end+1) = alt_zonal(:,end);
lat_zonal = repmat(gcm_Plat2D_edges_UM(:,1),[1 size(alt_zonal,1)]); lat_zonal = lat_zonal';
cli_zonal_timeav = meanNoNan(cli_zonal(iy,:,:),1);


cli_zonal_timeav2 = meanNoNan(cli,2); %[144 z]
alt_zonal2 = squeeze(z_3d_new(:,1,:)); alt_zonal2(:,end+1) = alt_zonal2(:,end)+dz; alt_zonal2(end+1,:) = alt_zonal2(end,:);
lat_zonal2 = repmat(gcm_Plat2D_edges_UM(:,1),[1 size(alt_zonal2,2)]); 

%altidue = [144 192 85] cli_time_av=[144   192    85]
cli_int = NaN*ones([size(cli_time_av,1) size(cli_time_av,2) size(z_new,2)]);
for ilat=1:size(cli_time_av,1)
    for ilon=1:size(cli_time_av,2)
        cli_int(ilat,ilon,:) = interp1(squeeze(altitude(ilat,ilon,:)),squeeze(cli_time_av(ilat,ilon,:)),z_new);
    end
end

cli_zonal_timeav3 = meanNoNan(cli_int,2); %[144 z]
alt_zonal3 = alt_zonal2;
lat_zonal3 = lat_zonal2;


figure; set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%set(gca,'position',[0.1300    0.100    0.7750    0.9150]);
set(gcf,'position',[3         297        1256         600]);
%dpcolor(lat_zonal,alt_zonal/1e3,cli_zonal_timeav); shading flat;
%dpcolor(lat_zonal2,alt_zonal2/1e3,cli_zonal_timeav2); shading flat;
dpcolor(lat_zonal3,alt_zonal3/1e3,cli_zonal_timeav3); shading flat;
%increase_font_size_map_figures
colorbar
set(gca,'ylim',[0 20]);
set(gca,'fontsize',18)
title({'UKESM1 ensemble mean ice MMR','2006-2014 time-average (kg/kg)'},'fontsize',18);
ylabel('Altitude abv. msl (km)','fontsize',18);
xlabel('Latitude (degrees)','fontsize',18);
%title('UKESM1 ensemble mean ice MMR 2006-2014 time-mean (kg/kg)');


%Do some maps for different height ranges
iz=find(z_new>17e3);
cli_17 = meanNoNan(cli_int(:,:,iz),3);

dat_modis = cli_17;
   
%run plotting script
figure; set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
UM_ACSIS_SW_vs_cloud_properties_global_DEFAULTS
icoarse_grain=0; M_coarse_grain=3; N_coarse_grain=3;
time_round=''; time_format_str='';
iplot_mgrid_lines_DRIVER=0;
ioverride_ticks_DRIVER=0;
irestrict_domain_DRIVER=0;
isave_plot=0;
iplot_wind_arrows=0;
subtitle_str = '17-20km mean ice MMR (kg/kg)';
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global


iz=find(z_new>14e3 & z_new<=17e3);
cli_14_17 = meanNoNan(cli_int(:,:,iz),3);
subtitle_str = '14-17km mean ice MMR (kg/kg)';
dat_modis = cli_14_17;
figure; set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global


iz=find(z_new>11e3 & z_new<=14e3);
cli_11_14 = meanNoNan(cli_int(:,:,iz),3);
subtitle_str = '11-14km mean ice MMR (kg/kg)';
dat_modis = cli_11_14;
figure; set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global

iz=find(z_new>8e3 & z_new<=11e3);
cli_8_11 = meanNoNan(cli_int(:,:,iz),3);
subtitle_str = '8-11km mean ice MMR (kg/kg)';
dat_modis = cli_8_11;
figure; set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global

iz=find(z_new>0.2e3 & z_new<=5e3);
cli_0_5 = meanNoNan(cli_int(:,:,iz),3);
subtitle_str = '0.2-5km mean ice MMR (kg/kg)';
dat_modis = cli_0_5;
figure; set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global

return


%% Calculate annual means, etc.
gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
%Had this the wrong way around :-
%[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
[gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

%manually set the time since matlab time doesn't work for the 360 day
%calendar (I think this is why the times are weird).
years_ukesm = repmat([start_year:end_year],[12 1]);
years_ukesm = years_ukesm';
years_ukesm_1d = years_ukesm(:,1);
months_ukesm = repmat([1:12],[length(years_ukesm) 1]);
days_ukesm = ones(size(months_ukesm));
time_ukesm = datenum(years_ukesm,months_ukesm,days_ukesm);

%For ISCCP pressure vs tau cloud fraction data
%size(dat_ens) = [nens nmonths ntau plev nlat nlon]

sdat=size(dat_ens);

%dat_std_ukesm = NaN*ones([length(years_ukesm_1d) size(dat_ens_mean,2) size(dat_ens_mean,3)]); %not used

if i_annual_data==1
    dat_annual = dat_ens_mean;
    dat_annual_ens = dat_ens;
    dat_ukesm = dat_ens_mean;
else
    dat_annual = NaN*ones([length(years_ukesm_1d) sdat(3:end)]);
    dat_ukesm = NaN*ones([length(years_ukesm_1d) 12 sdat(3:end)]);
    dat_annual_ens = NaN*ones([sdat(1) length(years_ukesm_1d) sdat(3:end)]);
    
    for iy=1:length(years_ukesm_1d) %Loop over years and use indices for the monthly data to create annual averages
        %for im=1:12
        tind_01 = ind_start_offset + (iy-1)*12 + 1;
        tind_02 = tind_01+11;
        %dat_annual(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %%calculate the annual means keeping each ensemble separate
        %dat_annual_ens(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        %Replaced the above with just one : at the end to allow for cases where
        %have more than just lat and lon dimensions (e.g., pressure, optical
        %depth for ISCCP). Should still work the same.
        dat_annual(iy,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:),1);
        a = dat_ens_mean(tind_01:tind_02,:);
        %dat_ukesm(iy,:) = dat_ens_mean(tind_01:tind_02,:);
        dat_ukesm(iy,:) = a(:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        %end
        
        %DJF means
        if iy==1 %skip the first year since don't have December
            dat_annual_DJF(iy,:,:) = NaN*ones([size(dat_ens_mean,2) size(dat_ens_mean,3)]);
            %dat_ukesm_DJF(iy,:,:,:) = NaN;
            dat_annual_ens_DJF(:,iy,:,:) = NaN*ones([size(dat_ens,1) size(dat_ens,3) size(dat_ens,4)]);
        else
            tind_01 = ind_start_offset + (iy-1)*12 + 0; %+1 is Jan, so +0 is Dec
            tind_02 = tind_01+2; %Feb
            dat_annual_DJF(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
            %dat_ukesm_DJF(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
            %calculate the annual means keeping each ensemble separate
            dat_annual_ens_DJF(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        end
        
        % MAM
        tind_01 = ind_start_offset + (iy-1)*12 + 3; %+1 is Jan
        tind_02 = tind_01+2; %May
        dat_annual_MAM(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_MAM(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_MAM(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        % JJA
        tind_01 = ind_start_offset + (iy-1)*12 + 6; %+1 is Jan
        tind_02 = tind_01+2; %Aug
        dat_annual_JJA(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_JJA(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_JJA(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        % SON
        tind_01 = ind_start_offset + (iy-1)*12 + 9; %+1 is Jan
        tind_02 = tind_01+2; %Nov
        dat_annual_SON(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_SON(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_SON(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
    end
    
end


%save(savefile,'/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat','-V7.3'...
%    ,'Nd_annual','Nd_ukesm','Nd_annual_ens','Nd_ens','Nd_ens_mean','Nd_ens_Ndatap','Nd_ens_std','gcm_Plat2D_UM','gcm_Plon2D_UM'...
%    ,'gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM','years_ukesm_1d',);

save(savefile,'-V7.3'); %save all variables


