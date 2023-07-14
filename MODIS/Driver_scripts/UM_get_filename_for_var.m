function [dat,file_UM]=UM_get_filename_for_var(dirUM_i,fileUM,VAR_NAME_STR_file,VAR_NAME_STR,UM_datatype)

dat.dat=NaN;

    file_UM = remove_character([dirUM_i fileUM],'VAR_NAME',VAR_NAME_STR_file);
    switch UM_datatype
        case '2D time mean'
            file_UM = remove_character(file_UM,'.pp.nc','.pp.nc_timeseries.mat');
            dat = load(file_UM,'time_mean_2D');
        case 'All data'
            file_UM = remove_character(file_UM,'.pp.nc','.pp.nc.mat');
             load(file_UM,VAR_NAME_STR);
             eval(['dat.dat=' VAR_NAME_STR ';']);
             
    end
    
  
    