try

savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';



%NOTE - AM3 fields say they are (grid_yt, grid_xt) orientated.
%However, it looks like (at least for tile 5) that grid_yt is actually
%latitude (or at least this is how Matlab reads it).

comp='UWchallenger';
%all have same lat range +44.2 to -44.2
%except two polar sterographic plots
%1= 305.6-34.2 (centred on lon=0)
%2= 35.8 - 124.2
%3 = north pole spherical 90 to 36 in lat approx 
%4 =  125.8-214.2 in lon
%5 is the best for the CAPT region, 215.8-304.2 (or -144.2 to -55.8)
%6 = South pole -36 to -90 in lat

if ~exist('ioverride_read_am3') | ioverride_read_am3==0
    iread_normal_vars=1;
    iread_cosp=1;
    tile_no = '5';
    gcm_idays=[1:324];
    gcm_idays=[1:10];
%    gcm_idays=[1];    
%    gcm_idays=[-1]; %set to <0 for all available times

% -----------------------------------
% select the required file in this script
       am3_choose_load_file
% -----------------------------------

end

gcm_str='AM3';
gcm_str_select='AM3';



%have output every 27 hours (324 in total).
if gcm_idays(1)<0
    gcm_time_read = nc_inst{'time'}(:); %days since 1st Jan 1980
    gcm_idays = [1:length(gcm_time_read)];
else
    gcm_time_read = nc_inst{'time'}(gcm_idays); %days since 1st Jan 1980
end


gcm_time_matlab = gcm_time_read + datenum('01-Jan-1980'); %convert to Matlab time
[Y,MO,D,H,MI,S] = datevec(gcm_time_matlab); %this outputs numbers for the date components


gcm_years_loaded_str=num2str(Y(1));

gcm_time_days = D;

%will be consistent with the actual calendar day rather than nearest day

%the hour of the day
gcm_time_UTC = H;
gcm_time_UTC_AM3 = H;
%the month
gcm_month = MO;

gcm_decimal_days = gcm_time_read - gcm_time_read(1);
%days of year since beginning of the year
daynum_timeseries3 = floor(gcm_time_matlab - datenum(Y(1),1,1) + 1);
daynum_timeseries3_AM3 = daynum_timeseries3;
%zero difference means day 1



% -----------------------------------
% lat lon grids
% -----------------------------------
gcm_phalf_ref=nc_grid{'phalf'}(:);

gcm_lat = load_am3('grid_latt','(:,:)',nc_grid,gcm_idays);
gcm_lat = squeeze(gcm_lat);
gcm_lon = load_am3('grid_lont','(:,:)',nc_grid,gcm_idays);
gcm_lon = squeeze(gcm_lon);

gcm_lat2 = load_am3('grid_lat','(:,:)',nc_grid,gcm_idays);
gcm_lat2 = squeeze(gcm_lat2);
gcm_lon2 = load_am3('grid_lon','(:,:)',nc_grid,gcm_idays);
gcm_lon2 = squeeze(gcm_lon2);



% -----------------------------------
%  now for the pressure & density
% -----------------------------------
%pfull(i,j,k) = (phalf(i,j,k+1)-phalf(i,j,k))) / ln(phalf(i,j,k+1)/phalf(i,j,k))  ;
gcm_pfull_ref = (gcm_phalf_ref(2:end)-gcm_phalf_ref(1:end-1)) ./ log(gcm_phalf_ref(2:end)./gcm_phalf_ref(1:end-1));

gcm_ps = load_am3('ps','(gcm_idays,:,:)',nc_inst,gcm_idays); %surface pressure (time,lat,lon)
gcm_bk = nc_inst{'bk'}(:); % bk and pk are constant with time and location
gcm_pk = nc_inst{'pk'}(:); % just vary vertically

gcm_phalf = am3_make_pressure(gcm_ps,gcm_pk,gcm_bk);
gcm_pfull = (gcm_phalf(:,2:end,:,:)-gcm_phalf(:,1:end-1,:,:)) ./ log(gcm_phalf(:,2:end,:,:)./gcm_phalf(:,1:end-1,:,:));
%this is the proper pressure in Pa

%gcm_temp = nc_inst{'temp'}(gcm_idays,:,:,:); %in K
gcm_temp = load_am3('temp','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
gcm_rho = density(gcm_pfull,gcm_temp); %now load this for speed
gcm_tsurf = load_am3('t_surf','(gcm_idays,:,:)',nc_inst,gcm_idays); %surface temperature in K
gcm_zsurf = load_am3('zsurf','(:,:)',nc_inst,gcm_idays); %surface height in m
gcm_theta = gcm_temp.*( (1000e2./gcm_pfull).^0.286 ); %calculate the potential temp
gcm_theta700 = squeeze(lininterp1f_multidim_RUN(gcm_pfull,gcm_theta,700e2,2)); %interpolate to find theta at 700mb
gcm_theta0 = squeeze(gcm_theta(:,end,:,:)); %theta at first model level
gcm_LTS = gcm_theta700 - gcm_theta0;

gcm_theta1000 = squeeze(lininterp1f_multidim_RUN(gcm_pfull,gcm_theta,1000e2,2)); %interpolate to find theta at 700mb
gcm_LTS1000 = gcm_theta700 - gcm_theta1000;

%qv
gcm_RH = load_am3('rh','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
f=1e6*28.97/18;
gcm_qsat = SatVapPress(gcm_temp,'goff','liq',gcm_pfull,1)/f; %kg/kg
gcm_qv = gcm_RH./100 .* gcm_qsat;
gcm_qv700 = squeeze(lininterp1f_multidim_RUN(gcm_pfull,gcm_qv,700e2,2)); %interpolate to find qv at 700mb


gcm_landmask = load_am3('land_mask','(:,:)',nc_inst,gcm_idays); %land fraction




% ---------------------------------------------
%  Calculate the height, assuming hydrostatic
% ---------------------------------------------
%the calculation below takes a long time - have saved the resulting .mat
%variables in the file below.

%[h_full,h_half]=am3_calc_height(gcm_pfull,gcm_temp,gcm_phalf);

%mat_file='AM3_height_grid_19900101.mat';
mat_file='AM3_height_grid_19900101_transpose.mat';
%save([nc_dir mat_file],'h_full','h_half','gcm_rho','-V7.3');

%%load([nc_dir mat_file],'h_full','h_half','gcm_rho');
%%h_half(:,end,:,:)=0; %first level of hhalf is the surface

% -----------------------------------
%  COSP only variables - don't forget to set the COSP flag to =1 in
%  read_am3
% -----------------------------------
if iread_cosp==1

cllcalipso = load_am3('cllcalipso','(gcm_idays,:,:,:)',nc_inst,gcm_idays,1); %CALIPSO simulated cloud fraction between 0 and 1
clmcalipso = load_am3('clmcalipso','(gcm_idays,:,:,:)',nc_inst,gcm_idays,1); %CALIPSO simulated cloud fraction between 0 and 1
clhcalipso = load_am3('clhcalipso','(gcm_idays,:,:,:)',nc_inst,gcm_idays,1); %CALIPSO simulated cloud fraction between 0 and 1

liqCF_modis = load_am3('lclmodis','(gcm_idays,:,:,:)',nc_inst,gcm_idays,1); %MODIS total (all heights) liquid cloud fraction
liqRe_modis = load_am3('lremodis','(gcm_idays,:,:,:)',nc_inst,gcm_idays,1)/100; %MODIS liquid water Re * CPCT
liqTau_modis = load_am3('ltaumodis','(gcm_idays,:,:,:)',nc_inst,gcm_idays,1)/100; %MODIS liquid water Tau * CPCT
%CPCT must be the cloud percentage (up to 100). This is the only way that
%the values make sense - so they are grid box means, not in-cloud means

%CFADs - each value is the fractional area with that dbz for the given
%height. So they should sum to one for each height - actually I think they
%sum to the cloud fraction (isn't a <-55 bin)



if ino_dbz_cfads==0



N_CFAD_bins = 15; %number of dbZ bins we want to load (starting from the highest dbZ)
Nmax=15; %total no. bins

cfad_alts_edges = load_am3('csatindx','(:)',nc_inst,gcm_idays); %height of CFAD bins (metres)

dz=mean(diff(cfad_alts_edges));
cfad_alts_edges = [cfad_alts_edges(1:end)-dz/2 cfad_alts_edges(end)+dz/2];
CFAD_dbz_edges_all = [-50:5:25];
CFAD_dbz_edges = CFAD_dbz_edges_all(Nmax-N_CFAD_bins+1:Nmax+1);


CFAD_temp = eval(['load_am3(''cloudsatcfad_' num2str(Nmax) ''',''(gcm_idays,:,:,:)'',nc_inst,gcm_idays,1);']); %
gcm_cloudsat_CFAD = NaN*ones([size(CFAD_temp,1) size(CFAD_temp,2) size(CFAD_temp,3) size(CFAD_temp,4) N_CFAD_bins]);

gcm_cloudsat_CFAD(:,:,:,:,N_CFAD_bins) = CFAD_temp;

for icfad=Nmax-1 : -1 : Nmax-N_CFAD_bins+1
    gcm_cloudsat_CFAD(:,:,:,:,icfad-Nmax+N_CFAD_bins) =  eval(['load_am3(''cloudsatcfad_' num2str(icfad) ''',''(gcm_idays,:,:,:)'',nc_inst,gcm_idays,1);']);
end

sCFAD=size(gcm_cloudsat_CFAD);

% b=squeeze(meanNoNan(gcm_cloudsat_CFAD,1));
% dbz = repmat(0.5*(CFAD_dbz_edges(1:end-1)+CFAD_dbz_edges(2:end)),[40 1 48 48]);
% dbz = permute(dbz,[1 3 4 2]);
% 
% ilat=13;
% lon_inds=[25:40];
% mean_dbz = sum( dbz.*b ,4) ./ sum(b,4);
% c=squeeze(mean_dbz(:,ilat,lon_inds));c(end+1,:)=NaN;
% figure
% lons = gcm_Plon2D_edges(ilat,lon_inds(1)-1:lon_inds(end));
% pcolor(lons,cfad_alts_edges,c);

% all %
dbz = repmat(0.5*(CFAD_dbz_edges(1:end-1)+CFAD_dbz_edges(2:end)),[sCFAD(1) 1 sCFAD(2) sCFAD(3) sCFAD(4)]);
dbz = permute(dbz,[1 3 4 5 2]);


%ilat=13;
%lon_inds=[25:40];
%the 'sum' flag does the sum instead of the mean
dbz_thresh=-30;
ithresh = find(CFAD_dbz_edges>=dbz_thresh);

[cf_CFAD_dbz] = meanNoNan(gcm_cloudsat_CFAD(:,:,:,:,ithresh(1):end),5,'sum',0); %might need this for averaging
%cf_CFAD_dbz will be the cloud fraction for each height now.
[mean_CFAD_dbz] = meanNoNan( 10.^(dbz/10).*gcm_cloudsat_CFAD ,5,'sum',0); 
mean_CFAD_dbz = 10*log10( mean_CFAD_dbz ./ cf_CFAD_dbz ); %are normalising by CF, so that this
%is the in-cloud dbz. May want to weight by CF when averaging up to give a
%grid box mean



%LIDAR scattering ratio
%N_CFAD_bins = 11; %number of dbZ bins we want to load (starting from the highest dbZ)
Nmax=15; %total no. bins
N_CFAD_bins = Nmax;

CFAD_sr_edges_all = [-100 0.01 1.2 3 5 7 10 15 20 25 30 40 50 60 80 999];
CFAD_sr_edges = CFAD_sr_edges_all(Nmax-N_CFAD_bins+1:Nmax+1);


CFAD_temp = eval(['load_am3(''calipsosrcfad_' num2str(Nmax) ''',''(gcm_idays,:,:,:)'',nc_inst,gcm_idays,1);']); %
gcm_calipso_CFAD = NaN*ones([size(CFAD_temp,1) size(CFAD_temp,2) size(CFAD_temp,3) size(CFAD_temp,4) N_CFAD_bins]);

gcm_calipso_CFAD(:,:,:,:,N_CFAD_bins) = CFAD_temp;

for icfad=Nmax-1 : -1 : Nmax-N_CFAD_bins+1
    gcm_calipso_CFAD(:,:,:,:,icfad-Nmax+N_CFAD_bins) =  eval(['load_am3(''calipsosrcfad_' num2str(icfad) ''',''(gcm_idays,:,:,:)'',nc_inst,gcm_idays,1);']);
end

sCFAD=size(gcm_calipso_CFAD);
dbz = repmat(0.5*(CFAD_sr_edges(1:end-1)+CFAD_sr_edges(2:end)),[sCFAD(1) 1 sCFAD(2) sCFAD(3) sCFAD(4)]);
dbz = permute(dbz,[1 3 4 5 2]);

%ilat=13;
%lon_inds=[25:40];

%Frequency over a certain range of dbzs/srs - N.B. the CFADs contain NaN values, so need
%to deal with those
%I think that this should give the cloud fraction above a certain sr
%threshold - which sr value to choose? Perhaps for thick Sc it is not so
%important.
sr_thresh=7;
ithresh = find(CFAD_sr_edges>=sr_thresh);

[cf_CFAD_sr,N2] = meanNoNan(gcm_calipso_CFAD(:,:,:,:,ithresh(1):end),5,'sum',0); %might need this for further averaging
mean_CFAD_sr = meanNoNan(10.^(dbz/10).*gcm_calipso_CFAD ,5,'sum',0); 
mean_CFAD_sr = 10*log10(mean_CFAD_sr ./ cf_CFAD_sr);
%actually the cf_CFAD_sr quantity might be the most useful thing here for
%CALIPSO since essentially it gives us the cloud fraction broken down
%vertically - the actual LIDAR sr values are perhaps not useful since the
%LIDAR beam is attenuated so quickly within a cloud. Maybe radar dbZs are
%useful, though. So prob best just to average cf_CFAD_sr.


% b=squeeze(meanNoNan(gcm_calipso_CFAD,1));
% sr = repmat(0.5*(CFAD_sr_edges(1:end-1)+CFAD_sr_edges(2:end)),[40 1 48 48]);
% sr = permute(sr,[1 3 4 2]);
% 
% ilat=13;
% lon_inds=[25:40];
% mean_sr = sum( sr.*b ,4 ) ./ sum(b,4);
% c=squeeze(mean_sr(:,ilat,lon_inds)); c(end+1,:)=NaN; c(:,end+1)=NaN;
% figure
% pcolor(lons,cfad_alts_edges,c); colorbar
% 
% pcfad = meanNoNan(gcm_phalf,1);
% ice_cfad = meanNoNan(gcm_iwc_av,1);
% ice_mean_slice = squeeze(ice_cfad(:,ilat,lon_inds)); ice_mean_slice(end+1,:)=NaN; ice_mean_slice(:,end+1)=NaN;
% p_mean_slice = squeeze(pcfad(:,ilat,lon_inds(1)-1:lon_inds(end))); 
% lons2D= repmat(lons,[size(p_mean_slice,1) 1]);
% figure
% pcolor(lons2D,p_mean_slice/100,ice_mean_slice); colorbar
% set(gca,'ydir','reverse');
% 
% lwc_cfad = meanNoNan(gcm_liq_av,1);
% lwc_mean_slice = squeeze(lwc_cfad(:,ilat,lon_inds)); lwc_mean_slice(end+1,:)=NaN; lwc_mean_slice(:,end+1)=NaN;
% figure
% pcolor(lons2D,p_mean_slice/100,lwc_mean_slice); colorbar
% set(gca,'ydir','reverse');


end

end

if iread_normal_vars==1
    
% -----------------------------------
%  standard variables
% -----------------------------------
gcm_cf = load_am3('cld_amt','(gcm_idays,:,:,:)',nc_inst,gcm_idays); %cloud fraction between 0 and 1
gcm_liq=load_am3('liq_wat','(gcm_idays,:,:,:)',nc_inst,gcm_idays) ./ gcm_cf; 
gcm_liq_av=load_am3('liq_wat','(gcm_idays,:,:,:)',nc_inst,gcm_idays);

gcm_iwc_av=load_am3('ice_wat','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
gcm_lwp_all_clouds=load_am3('LWP_all_clouds','(gcm_idays,:,:)',nc_inst,gcm_idays); %is this cloud averaged?
gcm_lwp=load_am3('LWP','(gcm_idays,:,:)',nc_inst,gcm_idays); %how to deal with CF?

gcm_iwp_all_clouds=load_am3('IWP_all_clouds','(gcm_idays,:,:)',nc_inst,gcm_idays); %is this cloud averaged?
gcm_iwp=load_am3('IWP','(gcm_idays,:,:)',nc_inst,gcm_idays); %how to deal with CF?

%convert precip to m/s to be consistent with CAM5 use h (m) = precip
%in (kg/m2) / rho_w
gcm_precT=load_am3('precip','(gcm_idays,:,:)',nc_inst,gcm_idays); %total precip rate kg/m2/s
gcm_precL=load_am3('prec_ls','(gcm_idays,:,:)',nc_inst,gcm_idays); %precip rate from stratiform cloud kg/m2/s
%(or the large scale rather, not convective)
gcm_precT = gcm_precT /1e3; %convert to m/s
gcm_precL = gcm_precL /1e3;

gcm_LSrain3D = load_am3('lscale_rain3d','(gcm_idays,:,:,:)',nc_inst,gcm_idays); %Rain fall rate from lscale  -3D 
% units: kg_h20 / m2 /s
gcm_LSrain3D = gcm_LSrain3D/1e3;

%gcm_LSsnow3D_snow = load_am3('lscale_snow3d','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
%gcm_LSsnow3D_snow = gcm_prec3D_snow/1e3;

gcm_Crain3D = load_am3('lscale_rain3d','(gcm_idays,:,:,:)',nc_inst,gcm_idays); %Rain fall rate from convection  -3D 
% units: kg_h20 / m2 /s
gcm_Crain3D = gcm_Crain3D/1e3;

%Don't read in yet, as we only have the rain fluxes and need qr....??
clear gcm_rain3D

%gcm_rain3D = [];
rwp=[]; rwp_isccp_low=[]; rwp_isccp_mid=[]; rwp_isccp_high=[]; 


%gcm_rain3D = gcm_LSrain3D + gcm_Crain3D;

%gcm_Csnow3D_snow = load_am3('lscale_Csnow3d','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
%gcm_Csnow3D_snow = gcm_prec3D_Csnow/1e3;




% -----------------------------------
%  droplet conc, LWP and cloud fraction
% -----------------------------------

%divide by cloud fractions to get in-cloud averages - beware values when
%have low CF, as will be inflated by zero cloud-fractions
gcm_drop_read=load_am3('liq_drp','(gcm_idays,:,:,:)',nc_inst,gcm_idays) ./ gcm_cf .*gcm_rho; %multiply by density to convert to per m3 (from per kg)
gcm_drop_read2=load_am3('liq_drp','(gcm_idays,:,:,:)',nc_inst,gcm_idays).*gcm_rho; % #/kg
%gcm_strat_Nd = nc_inst{'strat_droplet_number'}(gcm_idays,:,:,:) 

gcm_strat_Nd = load_am3('strat_droplet_number','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
gcm_strat_Nd(gcm_strat_Nd<-998)=NaN;
gcm_strat_Nd = gcm_strat_Nd .*gcm_rho; %multiply by density to convert to per m3 (from per kg);

gcm_REFFL = load_am3('strat_size_drop','(gcm_idays,:,:,:)',nc_inst,gcm_idays);
gcm_REFFL(gcm_REFFL<-998)=NaN;
gcm_REFFL = gcm_REFFL /2e6; %/2 for radius to diameter. And convert to metres.
%Think is already an in-cloud average

%ilowcf = find(gcm_cf<0.01);
%ilowcf = find(gcm_cf<0.8);
%gcm_REFFL(ilowcf) = NaN;

gcm_lwc=gcm_liq .* gcm_rho; %convert to kg/m3
gcm_lwc_av = gcm_liq_av.* gcm_rho; %convert to kg/m3;


LW_surf_down = load_am3('lwdn_sfc','(gcm_idays,:,:)',nc_inst,gcm_idays); %
LW_surf_up = load_am3('lwup_sfc','(gcm_idays,:,:)',nc_inst,gcm_idays); % 
LW_surf_net = LW_surf_down -  LW_surf_up; %presuming that for CAM net is positive downwards
%And want to be consistent with that (CAM only provides down and net for
%surface)
SW_surf_down = load_am3('swdn_sfc','(gcm_idays,:,:)',nc_inst,gcm_idays); %
SW_surf_up = load_am3('swup_sfc','(gcm_idays,:,:)',nc_inst,gcm_idays); %
SW_surf_net = SW_surf_down - SW_surf_up;

LW_TOA_net = load_am3('olr','(gcm_idays,:,:)',nc_inst,gcm_idays); % 
SW_TOA_up = load_am3('swup_toa','(gcm_idays,:,:)',nc_inst,gcm_idays); % 
SW_TOA_down = load_am3('swdn_toa','(gcm_idays,:,:)',nc_inst,gcm_idays); % 
SW_TOA_net = SW_TOA_down - SW_TOA_up;

albedo = SW_TOA_up ./ (SW_TOA_net + SW_TOA_up); %assuming that albedo = SWup/SWdown and SWnet = SWdown - SWup




end

% -----------------------------------
%  lat lon and time indices
% -----------------------------------
Plat=gcm_lat;
Plon=gcm_lon;

i180=find(Plon>180);
Plon(i180)=Plon(i180)-360;

% a = 0.5*(gcm_lat(:,2:end) + gcm_lat(:,1:end-1));
% b = 0.5*(a(2:end,:) + a(1:end-1,:));
% dlat_x = diff(gcm_lat(:,1:2));
% gcm_Plat2D_edges 

gcm_Plon2D=Plon;
gcm_Plat2D=Plat;

gcm_Plon2D_AM3=Plon;
gcm_Plat2D_AM3=Plat;

%not taking the last index as Matlab pcolor etc. decides to ignore the last strip of
%data and so needs grid arrays of the same size as the data....
%In the doc for surfc (on which pcolor is based) it says that can have data
%that is one dimension less than the X,Y data, but I can't get it to work?
%Have made a wrapper script called dpcolor or m_dpcolor that adds a row and
%column of NaNs to the C data to make it the same size as X and Y edges.
%Could start to use this if necessary
Plat=gcm_lat2(1:end-1,1:end-1);
Plon=gcm_lon2(1:end-1,1:end-1);

%Plat=gcm_lat2(1:end,1:end);
%Plon=gcm_lon2(1:end,1:end);


i180=find(Plon>180);
Plon(i180)=Plon(i180)-360;

gcm_Plon2D_edges=Plon;
gcm_Plat2D_edges=Plat;

gcm_Plon2D_edges_AM3=Plon;
gcm_Plat2D_edges_AM3=Plat;

%when plotting using m_pcolor it uses a plot-field that is the same size as
%the Plon,Plat grid - these are the cell edges in the plot so it only uses
%the field P(1:end-1,1:end-1) in the plot (ignores last row and column).
%So I don't include the end rows of gcm_lat2 for the edge arrays.

%dlon = 1.87, dlat = 2.5. Both constant throughout the grid for AM3

disp('Done read AM3 COSP');

  clear ioverride_read_am3
catch am3_error
    clear ioverride_read_am3
    rethrow(am3_error)
end






