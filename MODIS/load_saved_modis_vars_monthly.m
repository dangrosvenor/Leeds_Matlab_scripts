function [output,dat_size_out,MLAT,MLON,daynum_timeseries3,modisyear_timeseries3,aqua_terra_timeseries3] = load_saved_modis_vars_func(data_type,modis_varname,savedir_var,modis_data_case,iknow_size,dat_size,field_string,high_lat_correction,now_str)
 
%once have run this once we know the size of the time dimension and then
%can allocate the array to save memory
if iknow_size==1
    output = NaN*ones(dat_size);
end
        
daynum_timeseries3=[];
modisyear_timeseries3=[];
aqua_terra_timeseries3=[];
itime=1;
for imod=1:length(modis_data_case)  %loop through the files to load from

    % runs this script to determine the filename
    filename_choose_saved_MODIS_files
    
    
    switch data_type
        case {'L3 processed data','Nd vs SZA profiles at constant LAT'}
            %load in the variable required
            switch modis_varname  %these variables are contained in a separate file
                case {'Cloud_Top_Pressure_Day_Maximum','Cloud_Top_Pressure_Day_Mean','Cloud_Top_Pressure_Day_Minimum','Cloud_Top_Pressure_Day_Standard_Deviation','Cloud_Top_Temperature_Day_Maximum','Cloud_Top_Temperature_Day_Mean','Cloud_Top_Temperature_Day_Minimum','Cloud_Top_Temperature_Day_Standard_Deviation','Retrieved_Temperature_Profile_Mean'}
                    load(savevarname{2},modis_varname);
                otherwise %other variables
                    load(savevarname{1},modis_varname);
            end


        case 'Monthly averages etc from timeseries3'
            load(savevarname{3},modis_varname);
            modisyear_timeseries3 = modis_year;
    end
  
    %load lat and lon each time
    MLAT = load(savevarname{1},'MLAT');
    MLON = load(savevarname{1},'MLON');
    
    %coarsen resolution as a function of latitude
    if high_lat_correction==1
        eval_str = [modis_var{ivar} field_string{ivar} '_coarse=reduce_high_lats(' modis_var{ivar} field_string{ivar} ',LAT,''reduce'')'];
        eval(eval_str);
        save([savevarname{1} now_str],[modis_var{ivar} field_string{ivar} '_coarse'],'-APPEND');
    end

    

    
    switch data_type
        case 'L3 processed data'
                size_dat_new = eval( ['size( ' modis_varname field_string ');'] );
                output(:,:,itime:itime+size_dat_new(3)-1) = eval( [modis_varname field_string ';'] );


        case 'Nd vs SZA profiles at constant LAT'
            size_dat_new = eval( ['size( ' modis_varname ');'] );
            output(:,:,itime:itime+size_dat_new(3)-1) = eval( [modis_varname ';'] );
            
            load(savevarname{1},'LAT_sza');
            load(savevarname{1},'mid_sza');
            
        case 'Monthly averages etc from timeseries3'
            size_dat_new = eval( ['size( ' modis_varname ');'] );
            output = eval( [modis_varname ';'] );

    end
   
    
    itime=itime+size_dat_new(3);    
    
    clear savevarname %to make sure that it is set each time
    

end

dat_size_out=size(output);

switch data_type
    case 'L3 processed data'
%        LAT = eval( [modis_varname '.timeseries3_LAT;'] );
%        LON = eval( [modis_varname '.timeseries3_LON;'] );
        
    case 'Nd vs SZA profiles at constant LAT'       
        MLAT=LAT_sza;        
        MLON=mid_sza;
        
    otherwise
%        MLAT=0;
%        MLON=0;
        
end




