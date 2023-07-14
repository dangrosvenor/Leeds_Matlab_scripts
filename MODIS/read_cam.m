try
    
savemem=1;  %flag so that certain fields are not read in to save memory (mainly 3D fields)


savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';

comp='UWchallenger';

if ~exist('ioverride_read_am3') | ioverride_read_am3==0
    nc_dir = '/home/disk/eos8/d.grosvenor/CPT/CAM5/';
    %nc_grid_file = ['19900101.grid_spec.tile' tile_no '.nc']; gcm_str='AM3';
    nc_inst_file = ['cam5_1_17_dat.cam.h1.0001-01-01-00000.nc']; gcm_str='CAM5'; am3_dataset = 'CAM5_2deg';

%    nc_dir = '/home/disk/eos5/d.grosvenor/CPT/CAM5/new_output_2ndJuly_2012/';
%    nc_inst_file = ['camclubb51_AMIP_1deg.cam.h1.2007-01-29-64800.nc']; gcm_str='CAM5_CLUBB'; am3_dataset = 'CAM5_CLUBB';

    
    %have output every 27 hours (325 in total for cam5_1_17_dat.cam.h1.0001-01-01-00000.nc).
    gcm_idays=[1:325];
    gcm_idays=[1:10];
    
    %Select the lat lon range to extract
       %VOCALS region:-
       lat_range = [-40 10];
       lon_range = [-140 -50]+360;

end

low_cf_thresh=[0 1];
icf_low_only=0;
i_ice_screen=-1;

                
                

%nc_grid=netcdf([nc_dir nc_grid_file],'nowrite');
nc_inst=netcdf([nc_dir nc_inst_file],'nowrite');

latlon_method = 'experiment';
switch latlon_method
    case 'old'

        %gcm_phalf_ref=nc_grid{'phalf'}(:);
        gcm_lat_full = nc_inst{'lat'}(:);
        gcm_lon_full = nc_inst{'lon'}(:);

        gcm_slat_full = nc_inst{'slat'}(:);
        gcm_slon_full = nc_inst{'slon'}(:);

        %gcm_lat_full([1 end]) ---> ans = [-90 90]  &  size(gcm_lat_full) = [96 1], i.e. cell edges
        %gcm_slat_full([1 end]) ---> ans= [-89.0526 89.0526]  & size(gcm_lat_full)= [95 1], i.e. cell centers

        %gcm_lon_full([1 end]) --> ans = [0 357.500]  &  size(gcm_lon_full) = [144 1], likely cell edges
        %gcm_slon_full([1 end]) --> ans = [-1.2500 356.2500]  &  size(gcm_lon_full) = [144 1], likely cell centers


        ilat=find(gcm_slat_full>=lat_range(1) & gcm_slat_full<=lat_range(2));
        %go one either side
        %ilat = [max([1 ilat(1)-1]); ilat; min([length(gcm_lat_full) ilat(end)+1])];
        %ilat = [max([1 ilat(1)-1]):min([length(gcm_lat_full) ilat(end)+1])];

        %ilat = [ilat(1):ilat(end)];

        ilon=find(gcm_slon_full>=lon_range(1) & gcm_slon_full<=lon_range(2));
        %go one either side
        %ilon = [max([1 ilon(1)-1]):min([length(gcm_lon_full) ilon(end)+1])];
        %ilon = [ilon(1):ilon(end)];



        gcm_lat = gcm_lat_full(ilat);
        gcm_lon = gcm_lon_full(ilon);
        gcm_slat = gcm_slat_full(ilat);
        gcm_slon = gcm_slon_full([ilon; ilon(end)+1]);

        ilat=[ilat(2:end); ilat(end)+1];
        ilon=[ilon; ilon(end)+1];

    case 'experiment' %have matched this with the landmask file
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

        lat_range = [-40 10];
        lon_range = [-140 -50]+360;

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


end



gcm_time_read = nc_inst{'time'}(gcm_idays); %For CAM5 this is days since 1st Jan 0001
%gcm_time_matlab = gcm_time_read + datenum('01-Jan-0001'); %convert to Matlab time
%is different for different runs it seems...
gcm_time_matlab = gcm_time_read + datenum('01-Jan-2006'); %convert to Matlab time
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

gcm_ps = nc_inst{'PS'}(gcm_idays,ilat,ilon); %surface pressure (time,lat,lon), Pa
gcm_tsurf = nc_inst{'TS'}(gcm_idays,ilat,ilon); %radiative surface temperature (time,lat,lon), K
gcm_zsurf = nc_inst{'PHIS'}(gcm_idays,ilat,ilon)/9.80616; %surface geopotential divided by g should be surface height
gcm_landmask = squeeze(nc_inst{'LANDFRAC'}(1,ilat,ilon)); %radiative surface temperature (time,lat,lon), K

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
%  Calculate the density, assuming hydrostatic
% ---------------------------------------------
gcm_temp = nc_inst{'T'}(gcm_idays,:,ilat,ilon); %Air temperature in K
%will load in gcm_rho
gcm_rho = density(gcm_pfull,gcm_temp); %kg/m3

% ---------------------------------------------
%  Calculate the height, assuming hydrostatic
% ---------------------------------------------
%the calculation below takes a long time - have saved the resulting .mat
%variables in the file below.
%[h_full,h_half]=am3_calc_height(gcm_pfull,gcm_temp,gcm_phalf);
%mat_file='height_grid_cam5_1_17_dat.cam.h1.0001-01-01-00000.mat';
mat_file='height_grid_cam5_1_17_dat.cam.h1.0001-01-01-00000_new_grid_March30th2012.mat';
%save([nc_dir mat_file],'h_full','h_half','gcm_rho','-V7.3');
%load([nc_dir mat_file],'h_full','h_half','gcm_rho');
%h_half(:,end,:,:)=0; %first level of hhalf is the surface

% -----------------------------------
%  droplet conc, LWP and cloud fraction
% -----------------------------------
gcm_cf = nc_inst{'CLOUD'}(gcm_idays,:,ilat,ilon); %cloud fraction between 0 and 1

%divide by cloud fractions to get in-cloud averages - beware values when
%have low CF, as will be inflated by zero cloud-fractions
gcm_drop_read=nc_inst{'ICWNC'}(gcm_idays,:,ilat,ilon); %in per m3 - "Prognostic in-cloud water number conc"
gcm_nice=nc_inst{'ICINC'}(gcm_idays,:,ilat,ilon); %in per m3 - "Prognostic in-cloud ice number conc"
gcm_ccn3=nc_inst{'CCN3'}(gcm_idays,:,ilat,ilon)*1e6; %converted from per cm3 to per m3 for consistency
%CCN concentration at 0.1% - grid-box average I guess (since CCN in-cloud
%average seems unlikely)
gcm_awnc=nc_inst{'AWNC'}(gcm_idays,:,ilat,ilon); %in per m3 - "Average cloud water number conc"

gcm_REFFL=1e-6*nc_inst{'REL'}(gcm_idays,:,ilat,ilon); %in um (convert to metres). Effective radius liquid from MG.
%could compare this to that calculated from LWC and Nd assuming monomodal
%droplets

%CLDLIQ = Grid-box averaged water mixing ratio (kg/kg)
gcm_liq_av=nc_inst{'CLDLIQ'}(gcm_idays,:,ilat,ilon); %kg/kg
gcm_liq=gcm_liq_av ./ gcm_cf;
gcm_lwc_av=gcm_liq_av .* gcm_rho; 
gcm_lwc=gcm_liq .* gcm_rho; 
gcm_iwc_av=nc_inst{'CLDICE'}(gcm_idays,:,ilat,ilon); %kg/kg

gcm_lwp=nc_inst{'TGCLDLWP'}(gcm_idays,ilat,ilon); %Total grid-box LWP (kg/m2)
gcm_iwp=nc_inst{'TGCLDIWP'}(gcm_idays,ilat,ilon); %IWP

gcm_precL=nc_inst{'PRECL'}(gcm_idays,ilat,ilon); %large scale (stable) precip rate (m/s)
gcm_precT=nc_inst{'PRECT'}(gcm_idays,ilat,ilon); %TOTAL (large scale+convective) precip rate (m/s)

if savemem==0
    gcm_qr3D=nc_inst{'AQRAIN'}(gcm_idays,:,ilat,ilon); %3D rain mixing ratio (kg/kg)
end

if savemem==0
    gcm_omega3D=nc_inst{'OMEGA'}(gcm_idays,:,ilat,ilon); %Vertical pressure velocity (Pa/s)
end
gcm_omega500=nc_inst{'OMEGA500'}(gcm_idays,ilat,ilon); %Vertical pressure velocity at 500 hPa (Pa/s)

if savemem==0
    gcm_U3D=nc_inst{'U'}(gcm_idays,:,ilat,ilon); %Zonal wind speed (m/s)
end
gcm_U10=nc_inst{'U10'}(gcm_idays,ilat,ilon); %10m zonal wind speed (m/s)

if savemem==0
    gcm_V3D=nc_inst{'V'}(gcm_idays,:,ilat,ilon); %Meridional wind speed (m/s)
end
gcm_V10=nc_inst{'V10'}(gcm_idays,ilat,ilon); %10m meridional wind speed (m/s)

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

 clear ioverride_read_am3
catch am3_error
    clear ioverride_read_am3
    rethrow(am3_error)
end










