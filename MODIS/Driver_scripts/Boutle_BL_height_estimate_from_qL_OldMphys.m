%Estimate the boundary layer height using mean qL profile for Old Mphys run

%Old mphys
%nc_qL=netcdf('/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/xmmz-n/xmmzn_qL_.pp.nc');
%Control run
nc_qL=netcdf('/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/xmmz-u/xmmzu_qL_.pp.nc');

qL_prof = zeros([70 1]);

nt=97;
for it=1:nt
    qL = nc_qL{'qL'}(it,:,:,:);
    qL_prof = qL_prof + meanNoNan(meanNoNan(qL,2),2)/97;        
end

clear nc_qL

%Shows that the highest level qL reaches is model level 18, which
%corresponds to a height of 1133m. This similar, or even a bit less than
%for the other runs (based on the radar reflectivity profiles).
%Could check those with this technique.