file_dir = '/home/disk/eos15/d.grosvenor/ESGF/';

models = {};
%models{end+1} = {'BCC-ESM1'};
% models{end+1} = {'CNRM-CM6-1'};
% models{end+1} = {'CNRM-ESM2-1'};
 %models{end+1} = {'EC-Earth3-AerChem'};
% %models{end+1} = {'MPI-ESM-1-2-HAM'};
% models{end+1} = {'IPSL-CM6A-LR'};
% models{end+1} = {'MIROC6'};
% models{end+1} = {'MIROC-ES2L'};
% %models{end+1} = {'MRI-ESM2-0'};
% models{end+1} = {'CESM2'};
% models{end+1} = {'CESM2-FV2'};
% models{end+1} = {'CESM2-WACCM'};
% models{end+1} = {'CESM2-WACCM-FV2'};
% models{end+1} = {'NorESM2-MM'};
% models{end+1} = {'NIMS-KMAUKESM1-0-LL'};
 models{end+1} = {'NOAA-GFDL/GFDL-CM4'};
% models{end+1} = {'GFDL-ESM4'};
 models{end+1} = {'MOHC/UKESM1-0-LL'};

var = 'cldnvi';

%UKESM
nc_file = [file_dir '/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Emon/cldnvi/gn/latest/cldnvi_Emon_UKESM1-0-LL_historical_r1i1p1f2_gn_185001-194912.nc3'];
nc01 = netcdf(nc_file);    
cldnvi_ukesm = nc01{var}(:);
nc01=close(nc01);
inan=find(cldnvi_ukesm>1e19);
cldnvi_ukesm(inan)=NaN;
dat_ukesm = meanNoNan(meanNoNan(cldnvi_ukesm(1,:,:),2),2);


%GFDL-CM4
nc_file = [file_dir 'NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/cldnvi/cldnvi_Eday_GFDL-CM4_historical_r1i1p1f1_gr2_18500101-18691231.nc3'];
nc01 = netcdf(nc_file);    
cldnvi_gfdl_cm4 = nc01{var}(:);
nc01=close(nc01);
%inan=find(cldnvi_ukesm>1e19);
%cldnvi_ukesm(inan)=NaN;
dat_gfdl_cm4 = meanNoNan(meanNoNan(cldnvi_gfdl_cm4(1,:,:),2),2);

%CESM2
nc_file = [file_dir 'NCAR/CESM2/historical/r1i1p1f1/Emon/cldnvi/gn/latest/merged.nc'];
nc01 = netcdf(nc_file);    
cldnvi_cesm2 = nc01{var}(:);
nc01=close(nc01);
%inan=find(cldnvi_ukesm>1e19);
%cldnvi_ukesm(inan)=NaN;
dat_cesm2 = meanNoNan(meanNoNan(cldnvi_cesm2(1,:,:),2),2);

var2='clw';
%CESM2
nc_file = [file_dir 'NCAR/CESM2/historical/r1i1p1f1/clw/clw_t0.nc3'];
nc01 = netcdf(nc_file);    
clw_cesm2 = nc01{var2}(:);
nc01=close(nc01);
%inan=find(cldnvi_ukesm>1e19);
%cldnvi_ukesm(inan)=NaN;
dat_cesm2 = meanNoNan(meanNoNan(cldnvi_cesm2(1,:,:),2),2);



