
%dire = 'Y:\WRF\';

%rundir='ant_jan06_ecmwf_Nov08'; fileWRF(1).file=['met_em_d01-2006-01-05_00-00-00'];

dir_ncep = 'Y:\WRF\ant_jan06_sfRUC_v3\';
dir_netcdf = 'Y:\WRF\ant_jan06_ecmwf_Nov08\';
dir_ecmwf = 'Y:\WRF\jan06_met_orig\';

%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08_met_em/';

%file_ncep = 'met_em_d03-06_18-00-00';
%file_ecmwf = 'met_em_d03-01-06_18_ecmwf';

files = dir([dir_ecmwf '/met_em*d01*05_00*']);

for i=1:length(files)
	file_ncep = [dir_ncep files(i).name]
	file_ecmwf = [dir_ecmwf files(i).name]   
    file_netcdf = [dir_netcdf files(i).name]
	nc_ncep = netcdf(file_ncep); %open ncep file for reading
    nc_ecmwf = netcdf(file_ecmwf); %open ncep file for reading
    
    
    tic
    
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
    
    
    % global attributes
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
%        netcdf_add_variable(file_netcdf,varstruct);
         
	varstruct.Dimension={'Time', 'num_sm_levels', 'south_north', 'west_east'};
		varstruct.Name='SM';
        varstruct.FieldType = int32(104) ;
		varstruct.MemoryOrder = 'XYZ';
		varstruct.units = '';
		varstruct.description = '';
		varstruct.stagger = 'M';
        varstruct.Nctype = nc_float;
        
        netcdf_add_variable_data(file_netcdf,varstruct,nc_ncep);
%         netcdf_add_variable(file_netcdf,varstruct);
        

 %        Times=nc_ncep{'Times'}(:);
 %        nc_varput(file_netcdf,'Times',Times);
                          
%         data=nc_ecmwf{'PRES'}(:);
%         nc_varput(file_netcdf,'PRES',data);
                  
 %        data=nc_ncep{'SM'}(:);
 %        nc_varput(file_netcdf,'SM',data);



        close(nc_ncep);
        close(nc_ecmwf);
end

toc

    
