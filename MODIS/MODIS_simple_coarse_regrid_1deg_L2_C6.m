thresh_CTH=[0 3.2];
thresh_SZA=[0 65];

imodis_file_override = 1;
open_L2_C6_MODIS_file_01

%first filter for >3km and SZA>65
N37_filtered = N37;

cth_1km_from_5km = griddata(Plat_L2_5km,Plon_L2_5km,cth_5km,Plat2_L2,Plon2_L2);
sza_1km_from_5km = griddata(Plat_L2_5km,Plon_L2_5km,solar_zenith,Plat2_L2,Plon2_L2);

inan = find(cth_1km_from_5km/1e3<thresh_CTH(1) | cth_1km_from_5km/1e3>thresh_CTH(2));
N37_filtered(inan) = NaN;
inan = find(sza_1km_from_5km<thresh_SZA(1) | sza_1km_from_5km>thresh_SZA(2));
N37_filtered(inan) = NaN;

%coarse grain
N=100; M=100;
[mat_new]=reduce_matrix_subsample_mean(N37_filtered,N,M);
[lat_new]=reduce_matrix_subsample_mean(Plat2_L2,N,M);
[lon_new]=reduce_matrix_subsample_mean(Plon2_L2,N,M);

%Re-Grid to 1x1 deg

MLAT=[-89.5:89.5];
MLON=[-179.5:179.5];
[MLAT2d,MLON2d]=meshgrid(MLAT,MLON);
[MLON2d,MLAT2d]=meshgrid(MLON,MLAT);
Nd37_1deg = griddata(lat_new,lon_new,mat_new,MLAT2d,MLON2d);

