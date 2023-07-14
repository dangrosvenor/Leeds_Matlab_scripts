filepath='/home/disk/eosftp/pub/d.grosvenor/JessHo/BL_nuc_run_u-ca232/';

var = 'nuc_mass_H2SO4_ukca';
nc=netcdf([filepath var '/ca232a.pm2017jan_201701161200_' var '_slice_at_iz=0_z=0_saved.nc']);
nuc=nc{var}(:);
nc=close(nc);

var = 'aitken_mass_H2SO4_ukca';
nc=netcdf([filepath var '/ca232a.pm2017jan_201701161200_' var '_slice_at_iz=0_z=0_saved.nc']);
aitken=nc{var}(:);
nc=close(nc);

var = 'accum_mass_H2SO4_ukca';
nc=netcdf([filepath var '/ca232a.pm2017jan_201701161200_' var '_slice_at_iz=0_z=0_saved.nc']);
accum=nc{var}(:);
nc=close(nc);

var = 'coarse_mass_H2SO4_ukca';
nc=netcdf([filepath var '/ca232a.pm2017jan_201701161200_' var '_slice_at_iz=0_z=0_saved.nc']);
coarse=nc{var}(:);
nc=close(nc);

tot = nuc + aitken + accum + coarse;

qpcolor(tot);

