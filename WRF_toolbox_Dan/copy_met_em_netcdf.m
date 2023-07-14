
%dire = 'Y:\WRF\';

%rundir='ant_jan06_ecmwf_Nov08'; fileWRF(1).file=['met_em_d01-2006-01-05_00-00-00'];

dir_read = 'Y:\WRF\ant_jan06_sfRUC_v3\';
dir_netcdf = 'Y:\WRF\ant_jan06_ecmwf_Nov08\';
dir_ecmwf = 'Y:\WRF\jan06_met_orig\';

%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08_met_em/';

%file_ncep = 'met_em_d03-06_18-00-00';
%file_ecmwf = 'met_em_d03-01-06_18_ecmwf';

files = dir([dir_ecmwf '/met_em*']);

for i=1:length(files)
    cd(dir_read);
	file_read = [files(i).name]
	file_ecmwf = [dir_ecmwf files(i).name]   
    file_netcdf = [dir_netcdf files(i).name]
	nc_read = netcdf(file_read); %open ncep file for reading
%    nc_ecmwf = netcdf(file_ecmwf); %open ncep file for reading

varname='SKINTEMP';
vardata=nc_read{varname}(:);
nc_varput(file_netcdf,varname,vardata);       



% 
% varstruct.Dimension={'Time', 'z-dimension0024', 'south_north', 'west_east'};
% 		varstruct.Name='LANDUSEF';
%         varstruct.FieldType = int32(104) ;
% 		varstruct.MemoryOrder = 'XYZ';
% 		varstruct.units = 'category';
% 		varstruct.description = '24-category USGS landuse';
% 		varstruct.stagger = 'M';
%         varstruct.Nctype = nc_float;
%         
%         nc_varrename(file_netcdf,varstruct.Name,'DEL_LANDUSEF') 
%         
%         netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
        
        
        
    
    %    close(nc_ncep);
        close(nc_read);
end

    
