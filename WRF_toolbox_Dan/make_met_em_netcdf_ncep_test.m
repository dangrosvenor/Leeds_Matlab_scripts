
%dire = 'Y:\WRF\';

%rundir='ant_jan06_ecmwf_Nov08'; fileWRF(1).file=['met_em_d01-2006-01-05_00-00-00'];

dir_ncep = 'Y:\WRF\ant_jan06_sfRUC_v3\';
%dir_netcdf = 'Y:\WRF\ant_jan06_ecmwf_Nov08\';
dir_netcdf = 'C:\Documents and Settings\dan\My Documents\WRF\ant_jan06_met_ncep_test\';
dir_ecmwf = 'Y:\WRF\jan06_met_orig\';

%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08_met_em/';

%file_ncep = 'met_em_d03-06_18-00-00';
%file_ecmwf = 'met_em_d03-01-06_18_ecmwf';

files = dir([dir_ecmwf '/met_em*']);

for i=1:length(files)
	file_ncep = [dir_ncep files(i).name]
	file_ecmwf = [dir_ecmwf files(i).name]   
    file_netcdf = [dir_netcdf files(i).name]
	nc_ncep = netcdf(file_ncep); %open ncep file for reading    
    nc_ecmwf = netcdf(file_ecmwf); %open ncep file for reading
%    nc_ecmwf = nc_ncep;
    
    icreate=1;
    if icreate==1

        %mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
        %nc_create_empty ( file_netcdf, mode);
        nc_create_empty ( file_netcdf);

    end
    
    data = nc_ecmwf{'PRES'}(:);
    west_east = (size(data,2));
    south_north = (size(data,3));
    num_metgrid_levels = (size(data,1));
    
    
    % global attributes - add these first as otherwise takes ages to add them to the end
        var='TITLE';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='SIMULATION_START_DATE';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='WEST-EAST_GRID_DIMENSION';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='SOUTH-NORTH_GRID_DIMENSION';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); south_north+1 ;
		var='BOTTOM-TOP_GRID_DIMENSION';val=nc_attget(file_ecmwf,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); num_metgrid_levels ;
		var='WEST-EAST_PATCH_START_UNSTAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='WEST-EAST_PATCH_END_UNSTAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); west_east ;
		var='WEST-EAST_PATCH_START_STAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='WEST-EAST_PATCH_END_STAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); west_east+1 ;
		var='SOUTH-NORTH_PATCH_START_UNSTAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='SOUTH-NORTH_PATCH_END_UNSTAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); south_north ;
		var='SOUTH-NORTH_PATCH_START_STAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='SOUTH-NORTH_PATCH_END_STAG';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); south_north+1 ;
		var='GRIDTYPE';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='DX';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='DY';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1875.;
		var='DYN_OPT';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 2 ;
		var='CEN_LAT';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='CEN_LON';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='TRUELAT1';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='TRUELAT2';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='MOAD_CEN_LAT';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='STAND_LON';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val);
		var='POLE_LAT';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='POLE_LON';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='corner_lats';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='corner_lons';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='MAP_PROJ';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 2 ;
		var='MMINLU';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 
		var='ISWATER';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 16 ;
		var='ISICE';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 24 ;
		var='ISURBAN';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='ISOILWATER';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 14 ;
		var='grid_id';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 3 ;
		var='parent_id';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 2 ;
		var='i_parent_start';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 145 ;
		var='j_parent_start';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 145 ;
		var='i_parent_end';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 256 ;
		var='j_parent_end';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 256 ;
		var='parent_grid_ratio';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 4 ;
		var='FLAG_METGRID';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_SNOW';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_PSFC';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_SM000010';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_SM010200';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_ST000010';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_ST010200';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_SLP';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_SOILHGT';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
		var='FLAG_MF_XY';val=nc_attget(file_ncep,nc_global,var); nc_attput(file_netcdf,nc_global,var,val); 1 ;
        
              

	nc_add_dimension(file_netcdf,'Time',1); %zero means to make unilimited
    nc_add_dimension(file_netcdf,'DateStrLen' , 19);
	nc_add_dimension(file_netcdf,'west_east' , west_east);
	nc_add_dimension(file_netcdf,'south_north' , south_north);
	nc_add_dimension(file_netcdf,'num_metgrid_levels' , num_metgrid_levels);
	nc_add_dimension(file_netcdf,'num_sm_levels' , 2);
	nc_add_dimension(file_netcdf,'num_st_levels' , 2);
	nc_add_dimension(file_netcdf,'south_north_stag' , south_north+1);
	nc_add_dimension(file_netcdf,'west_east_stag' , west_east+1);
	nc_add_dimension(file_netcdf,'z-dimension0012' , 12);
	nc_add_dimension(file_netcdf,'z-dimension0016' , 16);
	nc_add_dimension(file_netcdf,'z-dimension0024', 24);
  
    clear varstruct
    
   varstruct.Dimension={'Time','DateStrLen'};
		varstruct.Name='Times';
        varstruct.Nctype = nc_char;

         nc_addvar(file_netcdf,varstruct);
         Times=nc_ncep{'Times'}(:);
         nc_varput(file_netcdf,'Times',Times);

    
	varstruct.Dimension={'Time', 'num_metgrid_levels', 'south_north', 'west_east'};
		varstruct.Name='PRES';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = '';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
        
         
	varstruct.Dimension={'Time', 'num_sm_levels', 'south_north', 'west_east'};
		varstruct.Name='SM';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = '';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
        
	varstruct.Dimension={'Time', 'num_st_levels', 'south_north', 'west_east'};
		varstruct.Name='ST';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = '';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'num_metgrid_levels', 'south_north', 'west_east'};
		varstruct.Name='GHT';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'm';
		varstruct.description = 'Height';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SOILHGT';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'm';
		varstruct.description = 'Terrain field of source analysis';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SNOW';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'kg m-2';
		varstruct.description = 'Water equivalent snow depth';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SKINTEMP';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'K';
		varstruct.description = 'Skin temperature (can use for SST also)';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='LANDSEA';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '';
		varstruct.description = 'Land/Sea flag';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_int;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SEAICE';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '';
		varstruct.description = 'Ice flag';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='ST010200';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'K';
		varstruct.description = 'T of 10-200 cm ground layer';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='ST000010';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'K';
		varstruct.description = 'T of 0-10 cm ground layer';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SM010200';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '';
		varstruct.description = 'Soil Moisture of 10-200 cm ground layer';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SM000010';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '';
		varstruct.description = 'Soil Moisture of 0-10 cm ground layer';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='PSFC';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'Pa';
		varstruct.description = 'Surface Pressure';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SPECHUMD';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'kg kg-1';
		varstruct.description = 'Specific Humidity';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'num_metgrid_levels', 'south_north', 'west_east'};
		varstruct.Name='RH';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = 'Relative Humidity';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'num_metgrid_levels', 'south_north_stag', 'west_east'};
		varstruct.Name='VV';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'm s-1';
		varstruct.description = 'V';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'num_metgrid_levels', 'south_north', 'west_east_stag'};
		varstruct.Name='UU';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'm s-1';
		varstruct.description = 'U';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'num_metgrid_levels', 'south_north', 'west_east'};
		varstruct.Name='TT';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'K';
		varstruct.description = 'Temperature';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='PMSL';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'Pa';
		varstruct.description = 'Sea-level Pressure';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ecmwf);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SLOPECAT';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'category';
		varstruct.description = 'Dominant category';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SNOALB';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'percent';
		varstruct.description = 'Maximum snow albedo';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'z-dimension0012', 'south_north', 'west_east'};
		varstruct.Name='GREENFRAC';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'fraction';
		varstruct.description = 'Monthly green fraction';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'z-dimension0012', 'south_north', 'west_east'};
		varstruct.Name='ALBEDO12M';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'percent';
		varstruct.description = 'Monthly surface albedo';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'z-dimension0016', 'south_north', 'west_east'};
		varstruct.Name='SOILCBOT';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'category';
		varstruct.description = '16-category bottom-layer soil type';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'z-dimension0016', 'south_north', 'west_east'};
		varstruct.Name='SOILCTOP';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'category';
		varstruct.description = '16-category top-layer soil type';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SOILTEMP';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'Kelvin';
		varstruct.description = 'Annual mean deep soil temperature';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north_stag', 'west_east'};
		varstruct.Name='HGT_V';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'meters MSL';
		varstruct.description = 'Topography height';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east_stag'};
		varstruct.Name='HGT_U';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'meters MSL';
		varstruct.description = 'Topography height';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SLPY';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '-';
		varstruct.description = 'df/dy';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SLPX';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '-';
		varstruct.description = 'df/dx';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='HGT_M';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'meters MSL';
		varstruct.description = 'Topography height';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='LU_INDEX';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'category';
		varstruct.description = 'Dominant category';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'z-dimension0024', 'south_north', 'west_east'};
		varstruct.Name='LANDUSEF';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = 'category';
		varstruct.description = '24-category USGS landuse';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='LANDMASK';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Landmask varstruct. 1=land, 0=water';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_int; 
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='COSALPHA';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Cosine of rotation angle';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='SINALPHA';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Sine of rotation angle';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='F';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '-';
		varstruct.description = 'Coriolis F parameter';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='E';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = '-';
		varstruct.description = 'Coriolis E parameter';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east_stag'};
		varstruct.Name='MAPFAC_UY';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor (y-dir) on U grid';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
        
	varstruct.Dimension={'Time', 'south_north_stag', 'west_east'};
		varstruct.Name='MAPFAC_VY';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor (y-dir) on V grid';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='MAPFAC_MY';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor (y-dir) on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east_stag'};
		varstruct.Name='MAPFAC_UX';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor (x-dir) on U grid';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north_stag', 'west_east'};
		varstruct.Name='MAPFAC_VX';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor (x-dir) on V grid';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='MAPFAC_MX';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor (x-dir) on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east_stag'};
		varstruct.Name='MAPFAC_U';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor on U grid';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north_stag', 'west_east'};
		varstruct.Name='MAPFAC_V';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor on V grid';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='MAPFAC_M';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'none';
		varstruct.description = 'Mapfactor on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='CLONG';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees longitude';
		varstruct.description = 'Computational longitude on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='CLAT';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees latitude';
		varstruct.description = 'Computational latitude on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east_stag'};
		varstruct.Name='XLONG_U';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees longitude';
		varstruct.description = 'Longitude on U grid';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east_stag'};
		varstruct.Name='XLAT_U';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees latitude';
		varstruct.description = 'Latitude on U grid';
		varstruct.stagger = 'U';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north_stag', 'west_east'};
		varstruct.Name='XLONG_V';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees longitude';
		varstruct.description = 'Longitude on V grid';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north_stag', 'west_east'};
		varstruct.Name='XLAT_V';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees latitude';
		varstruct.description = 'Latitude on V grid';
		varstruct.stagger = 'V';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='XLONG_M';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees longitude';
		varstruct.description = 'Longitude on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
         
	varstruct.Dimension={'Time', 'south_north', 'west_east'};
		varstruct.Name='XLAT_M';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XY ';
		varstruct.units = 'degrees latitude';
		varstruct.description = 'Latitude on mass grid';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);



        close(nc_ncep);
        close(nc_ecmwf);
end

    
