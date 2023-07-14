%rename on connolly to delete SM and ST for ecmwf files
%copy other fields (except ST and SM)
%add the new SM and ST on pc using the new dimension
%load SM and ST from ncep data on pc
%use nc_varput to put into ecmwf file


%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_sfRUC_v3/';

%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_sfRUC_v3/';
%dir_ecmwf = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08/';
dir_ncep = 'Y:\WRF\ant_jan06_sfRUC_v3\';
dir_ecmwf = 'Y:\WRF\ant_jan06_ecmwf_Nov08\';


%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08_met_em/';

%file_ncep = 'met_em_d03-06_18-00-00';
%file_ecmwf = 'met_em_d03-01-06_18_ecmwf';

files = dir([dir_ecmwf '/met_em*d01*05_18*']);

for i=1:length(files)
	file_ncep = [dir_ncep files(i).name]
	file_ecmwf = [dir_ecmwf files(i).name]
    file_netcdf= file_ecmwf;
	nc_ncep = netcdf(file_ncep); %open ncep file for reading
%	nc_ecmwf = netcdf(file_ecmwf,'write'); %open ecmwf file for writing


    nc_attput(file_netcdf,nc_global,'FLAG_SM000010', int32(1) );
    nc_attput(file_netcdf,nc_global,'FLAG_SM010200', int32(1) );
    nc_attput(file_netcdf,nc_global,'FLAG_ST000010', int32(1) );
    nc_attput(file_netcdf,nc_global,'FLAG_ST010200', int32(1) );

    nc_attput(file_netcdf,nc_global,'FLAG_SM100255', int32(0) );
    nc_attput(file_netcdf,nc_global,'FLAG_SM028100', int32(0) );
    nc_attput(file_netcdf,nc_global,'FLAG_SM007028', int32(0) );
    nc_attput(file_netcdf,nc_global,'FLAG_SM000007', int32(0) );
    
    nc_attput(file_netcdf,nc_global,'FLAG_ST100255', int32(0) );
    nc_attput(file_netcdf,nc_global,'FLAG_ST028100', int32(0) );
    nc_attput(file_netcdf,nc_global,'FLAG_ST007028', int32(0) );
    nc_attput(file_netcdf,nc_global,'FLAG_ST000007', int32(0) );
    
    nc_attput(file_netcdf,nc_global,'FLAG_SST', int32(0) );

    i_all=1;
if i_all==1
    
        
    ST=nc_ncep{'ST'}(:);
    SM=nc_ncep{'SM'}(:);
                
    close(nc_ncep);
   
    nc_add_dimension(file_ecmwf,'2num_sm_levels' , 2);
    
   varstruct.Dimension={'Time', '2num_sm_levels', 'south_north', 'west_east'};
		varstruct.Name='SM';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = '';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;        
    netcdf_add_variable(file_ecmwf,varstruct);
    
       varstruct.Dimension={'Time', '2num_sm_levels', 'south_north', 'west_east'};
		varstruct.Name='ST';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = '';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;        
    netcdf_add_variable(file_ecmwf,varstruct);

    
    nc_varput(file_ecmwf,'ST',ST);
    nc_varput(file_ecmwf,'SM',SM);
    
end
    
end
    
    
