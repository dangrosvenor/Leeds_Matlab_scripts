try
    
savemem=1;  %flag so that certain fields are not read in to save memory (mainly 3D fields)



comp='UWchallenger';

if ~exist('ioverride_read_am3') | ioverride_read_am3==0

    %Need to keep these inside of the override flag
    nc_dir = '/home/disk/eos8/d.grosvenor/CPT/CAM5/COSP_output/';
    nc_dir2 = '/home/disk/eos8/d.grosvenor/CPT/CAM5/';
    savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';

    nc_inst_file = ['cam5_1degcosp.cam.h0.0001-01.nc']; gcm_str='CAM5_COSP';
    %nc_inst_file = ['cam5_1degcosp.cam.h0.0002-01.nc']; gcm_str='CAM5_COSP';

    %nc_grid=netcdf([nc_dir nc_grid_file],'nowrite');
    nc_inst=netcdf([nc_dir nc_inst_file],'nowrite');

    %have monthly output for COSP - each file is one month
    gcm_idays=[-1];

    inew_h2_style = 0; %if want to set it to one then will need to set up nc_cosp_file filename for below
    if inew_h2_style==0
        nc_inst = nc_grid; %am calling it nc_grid to be consist with AM3, but is the h1 files
        %h2 files are nc_inst
    else   
        nc_inst = netcdf([nc_cosp_file],'nowrite');
    end
    
    %Select the lat lon range to extract
       %VOCALS region:-
    lat_range = [-40 10];
    lon_range = [-140 -50]+360;

end

%gcm_phalf_ref=nc_grid{'phalf'}(:);
%most arrays have these dimensions
        gcm_lat_full = nc_inst{'lat'}(:); %runs from -90 to 90 (192 pts for 1 degree data)
        gcm_lon_full = nc_inst{'lon'}(:); %runs from 0 to 358.75 (288 pts for 1 deg)

        gcm_lat_full = gcm_lat_full(1:end-2);
        gcm_lon_full = gcm_lon_full(1:end-2);

        gcm_slat_full = nc_inst{'slat'}(:); %runs from -89.5 to 89.5 (191 pts for 1 degree data)
        gcm_slon_full = nc_inst{'slon'}(:); %runs from -0.6250 (or =359.3750) to 358.1250 (288 pts for 1 deg)

        %gcm_lat_full([1 end]) ---> ans = [-90 90]  &  size(gcm_lat_full) = [96 1], i.e. cell edges
        %gcm_slat_full([1 end]) ---> ans= [-89.0526 89.0526]  & size(gcm_lat_full)= [95 1], i.e. cell centers

        %gcm_lon_full([1 end]) --> ans = [0 357.500]  &  size(gcm_lon_full) = [144 1], likely cell edges
        %gcm_slon_full([1 end]) --> ans = [-1.2500 356.2500]  &  size(gcm_lon_full) = [144 1], likely cell centers

        %decided that the data is at the cell centres, and we just ignore the first
        %and last values in latitude - gcm_slat_full then become the cell edges and
        %gcm_lat_full(2:end-1) the cell faces
        
        %in case we are loading the regional files (from 45S to 45N)
        %provided after the diurnal issues
        if length(gcm_lat_full)==0
            gcm_lat_full = nc_inst{'LAT_45s_to_45n'}(:); %runs from -44.76 to +44.76 (96 pts for 1 degree data)
            gcm_lon_full = nc_inst{'LON_150w_to_55w'}(:); %runs from 0 to 358.75 (288 pts for 1 deg)            
            %if extend these using the constant dlat = 0.9424 then can see
            %that runs from -90 to +90 like the "lat" above (not slat).
            %Same for lon.
            %So, assume that these are the cell centres as above
            %Now make cell edges (slat and slon)
            dlat=mean(diff(gcm_lat_full));
            dlon=mean(diff(gcm_lon_full));            
            
            gcm_slat_full = [gcm_lat_full(1)-dlat/2 : dlat : gcm_lat_full(end)+dlat/2];
            gcm_slon_full = [gcm_lon_full(1)-dlon/2 : dlon : gcm_lon_full(end)+dlon/2];            
        end
        



        ilat=find(gcm_slat_full>=lat_range(1) & gcm_slat_full<=lat_range(2));
        %go one either side
        %ilat = [max([1 ilat(1)-1]); ilat; min([length(gcm_lat_full) ilat(end)+1])];
        %ilat = [max([1 ilat(1)-1]):min([length(gcm_lat_full) ilat(end)+1])];

        %ilat = [ilat(1):ilat(end)];

        ilon=find(gcm_slon_full>=lon_range(1) & gcm_slon_full<=lon_range(2));
        %go one either side
        %ilon = [max([1 ilon(1)-1]):min([length(gcm_lon_full) ilon(end)+1])];
        %ilon = [ilon(1):ilon(end)];
        gcm_slat = gcm_slat_full(ilat);
        gcm_slon = gcm_slon_full([ilon]);

        %ilat = ilat(1:end-1);
        %ilon = ilon(1)-1:ilon(end)-2;

        ilat = ilat(1)+1:ilat(end);
        ilon = ilon(1):ilon(end)-1;

        
        gcm_lat = gcm_lat_full(ilat);
        gcm_lon = gcm_lon_full(ilon);


if length(gcm_idays)==1
    inotkeep=0;
else
    inotkeep=1;
end

gcm_time_read = nc_inst{'time'}(gcm_idays); %For CAM5 this is days since 1st Jan 0001
%But varies depending on data - can get this info from the units attribute
time_units=nc_inst{'time'}.units(:); %e.g. = days since 2006-01-01 00:00:00
istr = findstr(time_units,'days since '); istr_loc = istr + 11;
datestart = datenum(time_units(istr + 11:end));
gcm_time_matlab = gcm_time_read + datestart; %convert to Matlab time
[Y,MO,D,H,MI,S] = datevec(gcm_time_matlab); %this outputs numbers for the date components

gcm_time_days = D;

%will be consistent with the actual calendar day rather than nearest day

%the hour of the day
gcm_time_UTC = H;
%the month
gcm_month = MO;

gcm_decimal_days = gcm_time_read - gcm_time_read(1);
%days of year since beginning of the year
daynum_timeseries3 = floor(gcm_time_matlab - datenum(Y(1),1,1) + 1);
eval(['daynum_timeseries3_' gcm_str ' = daynum_timeseries3;']);
%zero difference means day 1




% -----------------------------------
%  now for the pressure 
% -----------------------------------
%pfull(i,j,k) = (phalf(i,j,k+1)-phalf(i,j,k))) / ln(phalf(i,j,k+1)/phalf(i,j,k))  ;
%gcm_pfull_ref = (gcm_phalf_ref(2:end)-gcm_phalf_ref(1:end-1)) ./ log(gcm_phalf_ref(2:end)./gcm_phalf_ref(1:end-1));

%gcm_siglev = nc_inst{'lev'}(:); %sigma co-ordinates. Don't think we need these.

%All the 3D variables are stored on lev(=30) levels. These are the
%mid-point levels. Level 30 is just above the surface. Level 31 of the
%interface index is the surface.


%gcm_p0 = nc_inst{'P0'}(:); %reference pressure (constant) - doesn't seem
%to be there - use 1000 hPa for now
gcm_p0 = 1000e2;
%fprintf(1,'\nUsing reference pressure of 1000 hPa - find the proper one! ***\n');
%actually 1000hPa is the correct reference pressure! (according to Andy
%Gettelman).

gcm_ps = squeeze_keep1(nc_inst{['PS' suffix]}(gcm_idays,ilat,ilon),inotkeep); %surface pressure (time,lat,lon), Pa


gcm_tsurf = squeeze_keep1(nc_inst{['TS' suffix]}(gcm_idays,ilat,ilon),inotkeep); %radiative surface temperature (time,lat,lon), K
gcm_zsurf = nc_inst{['PHIS' suffix]}(gcm_idays,ilat,ilon)/9.80616; %surface geopotential divided by g should be surface height
gcm_landmask = squeeze(nc_inst{['LANDFRAC' suffix]}(1,ilat,ilon)); %radiative surface temperature (time,lat,lon), K

gcm_hyai = nc_inst{'hyai'}(:); %a-coefficient for the model (interface) levels
gcm_hybi = nc_inst{'hybi'}(:); %b-coeff. %These just vary vertically

gcm_pki = gcm_hyai*gcm_p0;   
gcm_bki = gcm_hybi;

gcm_hyam = nc_inst{'hyam'}(:); %a-coefficient for the model levels
gcm_hybm = nc_inst{'hybm'}(:); %b-coeff. %These just vary vertically

gcm_pkm = gcm_hyam*gcm_p0; 
gcm_bkm = gcm_hybm;

%Use the formula gcm_phalf = pk + ps.*bk  to calculate the pressure
gcm_phalf = am3_make_pressure(gcm_ps,gcm_pki,gcm_bki); %
gcm_pfull = am3_make_pressure(gcm_ps,gcm_pkm,gcm_bkm); %
%(Pa)

%AM3 does it like this - so pfull refers to lev (mid-point indices) and phalf to ilev
%(interface indices) gcm_pfull = (gcm_phalf(:,2:end,:,:)-gcm_phalf(:,1:end-1,:,:)) ./ log(gcm_phalf(:,2:end,:,:)./gcm_phalf(:,1:end-1,:,:));


% ---------------------------------------------
%  Calculate the density
% ---------------------------------------------
gcm_temp = squeeze_keep1(nc_inst{['T' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %Air temperature in K
gcm_rho = density(gcm_pfull,gcm_temp); %kg/m3

%LTS
gcm_theta = gcm_temp.*( (1000e2./gcm_pfull).^0.286 ); %calculate the potential temp
gcm_theta700 = squeeze(lininterp1f_multidim_RUN(gcm_pfull,gcm_theta,700e2,2)); %interpolate to find theta at 700mb
gcm_theta0 = squeeze(gcm_theta(:,end,:,:)); %theta at first model level
gcm_LTS = gcm_theta700 - gcm_theta0;

gcm_theta1000 = squeeze(lininterp1f_multidim_RUN(gcm_pfull,gcm_theta,1000e2,2)); %interpolate to find theta at 700mb
gcm_LTS1000 = gcm_theta700 - gcm_theta1000;

%qv
gcm_RH = squeeze_keep1(nc_inst{['RELHUM' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %RH (%)
f=1e6*28.97/18;
gcm_qsat = SatVapPress(gcm_temp,'goff','liq',gcm_pfull,1)/f; %kg/kg
gcm_qv = gcm_RH./100 .* gcm_qsat;
gcm_qv700 = squeeze(lininterp1f_multidim_RUN(gcm_pfull,gcm_qv,700e2,2)); %interpolate to find qv at 700mb

LW_surf_down = squeeze_keep1(nc_inst{['FLDS' suffix]}(gcm_idays,ilat,ilon),inotkeep);
LW_surf_net = squeeze_keep1(nc_inst{['FLNS' suffix]}(gcm_idays,ilat,ilon),inotkeep); %presumably positive downwards?
LW_TOA_net = squeeze_keep1(nc_inst{['FLNT' suffix]}(gcm_idays,ilat,ilon),inotkeep);
SW_surf_down = squeeze_keep1(nc_inst{['FSDS' suffix]}(gcm_idays,ilat,ilon),inotkeep);
SW_surf_net = squeeze_keep1(nc_inst{['FSNS' suffix]}(gcm_idays,ilat,ilon),inotkeep);
SW_TOA_net = squeeze_keep1(nc_inst{['FSNTOA' suffix]}(gcm_idays,ilat,ilon),inotkeep);
SW_TOA_up = squeeze_keep1(nc_inst{['FSUTOA' suffix]}(gcm_idays,ilat,ilon),inotkeep);

albedo = SW_TOA_up ./ (SW_TOA_net + SW_TOA_up); %assuming that albedo = SWup/SWdown and SWnet = SWdown - SWup

% ---------------------------------------------
%  Calculate the height, assuming hydrostatic
% ---------------------------------------------
%the calculation below takes a long time - have saved the resulting .mat
%variables in the file below. N.B - can calculate LWP, etc using dp and the
%hydrostatic equation
%[h_full,h_half]=am3_calc_height(gcm_pfull,gcm_temp,gcm_phalf);
%mat_file='height_grid_cam5_1_17_dat.cam.h1.0001-01-01-00000.mat';
%mat_file='height_grid_cam5_1_17_dat.cam.h1.0001-01-01-00000_new_grid_March30th2012.mat';
%save([nc_dir mat_file],'h_full','h_half','gcm_rho','-V7.3');
%load([nc_dir2 mat_file],'h_full','h_half','gcm_rho');
%h_half(:,end,:,:)=0; %first level of hhalf is the surface

% -----------------------------------
%  droplet conc, LWP and cloud fraction
% -----------------------------------
gcm_cf = squeeze_keep1(nc_inst{['CLOUD' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %cloud fraction between 0 and 1

%divide by cloud fractions to get in-cloud averages - beware values when
%have low CF, as will be inflated by zero cloud-fractions
%gcm_drop_read=squeeze_keep1(nc_inst{'NUMLIQ'}(gcm_idays,:,ilat,ilon),inotkeep); %in per m3 - "Prognostic in-cloud water number conc"
gcm_drop_read=squeeze_keep1(nc_inst{['AWNC' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %see below
%am using AWNC because NUMLIQ is only present as a monthly average (may be
%what it is...?) for the latest COSP data (h0,h1&h2 files from Pete)
%AWNC is thresholded for a reasonable LWC and CF (looks like ICWNC is not). So it
%is zero when there is no reasonable cloud
gcm_nice=squeeze_keep1(nc_inst{['NUMICE' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %in per m3 - "Prognostic in-cloud ice number conc"
gcm_ccn3=squeeze_keep1(nc_inst{['CCN3' suffix]}(gcm_idays,:,ilat,ilon),inotkeep)*1e6; %converted from per cm3 to per m3 for consistency
%CCN concentration at 0.1% - grid-box average I guess (since CCN in-cloud
%average seems unlikely)
gcm_awnc=squeeze_keep1(nc_inst{['AWNC' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %in per m3 - "Average cloud water number conc"
gcm_icwnc=squeeze_keep1(nc_inst{['ICWNC' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %in per m3 - "Average cloud water number conc"
gcm_REFFL=1e-6*nc_inst{['REL' suffix]}(gcm_idays,:,ilat,ilon); %in um (convert to metres). Average cloud top effective radius liquid.



%CLDLIQ = Grid-box averaged water mixing ratio (kg/kg)
gcm_liq_av=squeeze_keep1(nc_inst{['CLDLIQ' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg
gcm_liq=gcm_liq_av ./ gcm_cf;








% **********************************************
%   Skip here for now
% **********************************************
iskip=0;
if iskip==0
    gcm_lwc_av=gcm_liq_av .* gcm_rho;
    gcm_lwc=gcm_liq .* gcm_rho;
    gcm_iwc_av=squeeze_keep1(nc_inst{['CLDICE' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg

    gcm_lwp=squeeze_keep1(nc_inst{['TGCLDLWP' suffix]}(gcm_idays,ilat,ilon),inotkeep); %Total grid-box LWP (kg/m2)
    gcm_iwp=squeeze_keep1(nc_inst{['TGCLDIWP' suffix]}(gcm_idays,ilat,ilon),inotkeep); %IWP

    gcm_precL=squeeze_keep1(nc_inst{['PRECL' suffix]}(gcm_idays,ilat,ilon),inotkeep); %large scale (stable) precip rate (m/s)
    gcm_precT=squeeze_keep1(nc_inst{['PRECT' suffix]}(gcm_idays,ilat,ilon),inotkeep); %TOTAL (large scale+convective) precip rate (m/s)
    gcm_rain3D=squeeze_keep1(nc_inst{['AQRAIN' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %3D rain mixing ratio (kg/kg) ("Average rain mixing ratio")
    gcm_rain3D = gcm_rain3D .* gcm_rho; %convert to kg/m3
    %also have
    if exist('gcm_Qrain3D') & prod(size(gcm_Qrain3D))>0
        gcm_Qrain3D=squeeze_keep1(nc_inst{['QRAIN' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %"Diagnostic grid-mean rain mixing ratio"
        gcm_Qrain3D = gcm_Qrain3D .* gcm_rho; %convert to kg/m3
    else
        gcm_Qrain3D =[];
    end

    if savemem==0
        gcm_omega3D=squeeze_keep1(nc_inst{['OMEGA' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %Vertical pressure velocity (Pa/s)
    end
    gcm_omega500=squeeze_keep1(nc_inst{['OMEGA500' suffix]}(gcm_idays,ilat,ilon),inotkeep); %Vertical pressure velocity at 500 hPa (Pa/s)

    if savemem==0
        gcm_U3D=squeeze_keep1(nc_inst{['U' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %Zonal wind speed (m/s)
    end
    gcm_U10=squeeze_keep1(nc_inst{['U10' suffix]}(gcm_idays,ilat,ilon),inotkeep); %10m zonal wind speed (m/s)

    if savemem==0
        gcm_V3D=squeeze_keep1(nc_inst{'V'}(gcm_idays,:,ilat,ilon),inotkeep); %Meridional wind speed (m/s)
    end
    gcm_V10=squeeze_keep1(nc_inst{'V10'}(gcm_idays,ilat,ilon),inotkeep); %10m meridional wind speed (m/s)


end %if iskip==1

cllcalipso = squeeze_keep1(nc_grid{['CLDLOW_CAL' suffix]}(gcm_idays,ilat,ilon),inotkeep); cllcalipso(cllcalipso<-1e29)=NaN;
clmcalipso = squeeze_keep1(nc_grid{['CLDMED_CAL' suffix]}(gcm_idays,ilat,ilon),inotkeep); clmcalipso(clmcalipso<-1e29)=NaN;
clhcalipso = squeeze_keep1(nc_grid{['CLDHGH_CAL' suffix]}(gcm_idays,ilat,ilon),inotkeep); clhcalipso(clhcalipso<-1e29)=NaN;

liqCF_modis = squeeze_keep1(nc_grid{['CLWMODIS' suffix]}(gcm_idays,ilat,ilon),inotkeep); liqCF_modis(liqCF_modis<-1e29)=NaN;
%Re and Tau below come as Re*CLWMODIS. Since CLWMODIS (liquid MODIS CF) is
%in % we divide them by 100 below to be consistent with other fields (i.e.
%the average over a GCM grid box). But do this after the -1e30 screening
liqRe_modis = squeeze_keep1(nc_grid{['REFFCLWMODIS' suffix]}(gcm_idays,ilat,ilon),inotkeep); liqRe_modis(liqRe_modis<-1e29)=NaN; liqRe_modis=liqRe_modis*1e-2;
liqTau_modis = squeeze_keep1(nc_grid{['TAUWMODIS' suffix]}(gcm_idays,ilat,ilon),inotkeep); liqTau_modis(liqTau_modis<-1e29)=NaN; liqTau_modis=liqTau_modis*1e-2;



%CPCT must be the cloud percentage (up to 100). This is the only way that
%the values make sense - so they are grid box means, not in-cloud means

%CFADs - each value is the fractional area with that dbz for the given
%height. So they should sum to one for each height - actually I think they
%sum to the cloud fraction (isn't a <-55 bin)


% ****************************************
%             CFADS
% ****************************************

cfad_alts_edges = squeeze_keep1(nc_grid{'cosp_ht'}(:),inotkeep); cfad_alts_edges(cfad_alts_edges<-1e29)=NaN;
%cfad_alts_edges = cfad_alts_edges';
dz=mean(diff(cfad_alts_edges));
%cfad_alts_edges = [cfad_alts_edges(1:end)-dz/2 cfad_alts_edges(end)+dz/2];
acat = cfad_alts_edges(1:end)-dz/2;
catdim=1;
if size(acat,1)==1; catdim=2; end
cfad_alts_edges = cat(catdim,acat,cfad_alts_edges(end)+dz/2);

if ino_dbz_cfads==0

N_CFAD_bins = 15; %number of dbZ bins we want to load (starting from the highest dbZ)
Nmax=15; %total no. bins

CFAD_dbz_edges_all = squeeze_keep1(nc_grid{'cosp_dbze_bnds'}(:),1); CFAD_dbz_edges_all(CFAD_dbz_edges_all<-1e29)=NaN;
%CFAD_dbz_edges_all = [-50:5:25];
CFAD_dbz_edges_all = (cat(1,CFAD_dbz_edges_all(:,1),CFAD_dbz_edges_all(end,2) ) )';




CFAD_dbz_edges = CFAD_dbz_edges_all(Nmax-N_CFAD_bins+1:Nmax+1);


gcm_cloudsat_CFAD = squeeze_keep1(nc_grid{'CFAD_DBZE94_CS'}(gcm_idays,:,:,ilat,ilon),inotkeep); gcm_cloudsat_CFAD(gcm_cloudsat_CFAD<-1e29)=NaN;
gcm_cloudsat_CFAD = permute(gcm_cloudsat_CFAD,[1 2 4 5 3]);
sCFAD=size(gcm_cloudsat_CFAD);

dbz = repmat(0.5*(CFAD_dbz_edges(1:end-1)+CFAD_dbz_edges(2:end)),[sCFAD(1) 1 sCFAD(2) sCFAD(3) sCFAD(4)]);
dbz = permute(dbz,[1 3 4 5 2]);

%the 'sum' flag does the sum instead of the mean
dbz_thresh=-30;
ithresh = find(CFAD_dbz_edges>=dbz_thresh);

[cf_CFAD_dbz] = meanNoNan(gcm_cloudsat_CFAD(:,:,:,:,ithresh(1):end),5,'sum',0); %might need this for averaging
[mean_CFAD_dbz] = meanNoNan( 10.^(dbz/10).*gcm_cloudsat_CFAD ,5,'sum',0); 
mean_CFAD_dbz = 10*log10( mean_CFAD_dbz ./ cf_CFAD_dbz ); %are normalising by CF, so that this
%is the in-cloud dbz. May want to weight by CF when averaging up to give a
%grid box mean

end

%  ***   LIDAR scattering ratio ***
%N_CFAD_bins = 11; %number of dbZ bins we want to load (starting from the highest dbZ)
Nmax=15; %total no. bins
N_CFAD_bins = Nmax;

CFAD_sr_edges_all = squeeze_keep1(nc_grid{'cosp_sr_bnds'}(:),1); CFAD_sr_edges_all(CFAD_sr_edges_all<-1e29)=NaN;
%CFAD_sr_edges_all = [-100 0.01 1.2 3 5 7 10 15 20 25 30 40 50 60 80 999];
CFAD_sr_edges_all = (cat(1,CFAD_sr_edges_all(:,1),CFAD_sr_edges_all(end,2) ) )';


CFAD_sr_edges = CFAD_sr_edges_all(Nmax-N_CFAD_bins+1:Nmax+1);

gcm_calipso_CFAD = squeeze_keep1(nc_grid{['CFAD_SR532_CAL' suffix]}(gcm_idays,:,:,ilat,ilon),inotkeep); gcm_calipso_CFAD(gcm_calipso_CFAD<-1e29)=NaN;
gcm_calipso_CFAD = permute(gcm_calipso_CFAD,[1 2 4 5 3]);

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


% ****************************************************************
%             New variables from Pre/Post Mphys for LWC runs
% ****************************************************************
  %CLDLIQ = Grid-box averaged water mixing ratio (kg/kg)
%gcm_lwcBP_av=gcm_rho.*squeeze_keep1(nc_inst{['CLDLIQBP' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg - before mphys
%gcm_lwcAP_av=gcm_rho.*squeeze_keep1(nc_inst{['CLDLIQAP' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg - after mphys
if exist('nc_prepostLWP_file')
    gcm_lwcBP_av=gcm_rho.*squeeze_keep1(nc_prepostLWP{['CLDLIQBMG' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg - before mphys
    gcm_lwcAP_av=gcm_rho.*squeeze_keep1(nc_prepostLWP{['CLDLIQAMG' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg - after mphys
end
gcm_lwcSEDTEN_av=gcm_rho.*squeeze_keep1(nc_inst{['QCSEDTEN' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg/s sedimentation tendency
gcm_lwcEVAPTEN_av=gcm_rho.*squeeze_keep1(nc_inst{['QCSEVAP' suffix]}(gcm_idays,:,ilat,ilon),inotkeep); %kg/kg/s rate of evap of falling cloud water
if prod(size(gcm_lwcBP_av))>0
    ilwcAPBP=1;
else
    ilwcAPBP=0;
end
% -----------------------------------
%  lat lon and time indices
% -----------------------------------
%cell edges (contained in gcm_slat)
Plat=gcm_slat;
Plon=gcm_slon;

Plat=gcm_lat_full;
Plon=gcm_lon_full;

Plat=gcm_slat_full;
Plon=gcm_slon_full;

Plat=gcm_slat;
Plon=gcm_slon;


i180=find(Plon>180);
Plon2=Plon;
Plon2(i180)=Plon2(i180)-360;

%Plon=[Plon2(1:end); Plon2(1)];
Plon=Plon2;



[gcm_Plon2D_edges,gcm_Plat2D_edges]=meshgrid(Plon,Plat);


%cell centres - contained in gcm_lat - seems that the fields are stored as point values at gcm_lat, which runs to lat=+/-90
% - so gcm_lat are effectively the cell face positions. And we ignore the
% first and last lat value (since cells can't extend beyond 90 degree lat)
Plat=gcm_lat;
Plon=gcm_lon;

i180=find(Plon>180);
Plon(i180)=Plon(i180)-360;

[gcm_Plon2D,gcm_Plat2D]=meshgrid(Plon,Plat);

%dlon = 1.87, dlat = 2.5. Both constant throughout the grid for CAM5

% dlat = diff(gcm_lat);
% dlon = diff(gcm_lon);
% Plat=gcm_lat+[dlat; dlat(end)]/2;
% Plon=gcm_lon+[dlon; dlon(end)]/2;
% 
% i180=find(Plon>180);
% Plon(i180)=Plon(i180)-360;
% 
% [gcm_Plon2D_edges,gcm_Plat2D_edges]=meshgrid(Plon,Plat);


fprintf(1,'\n Done read CAM5\n');


clear ioverride_read_am3
catch read_cam_error
clear ioverride_read_am3    
    rethrow(read_cam_error);
end




