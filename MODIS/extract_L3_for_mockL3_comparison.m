%extract a few days of data from a year of L3 for comparison to the mock L3

operation='load';
%operation='extract and save';

indices_str2 = '(:,:,164:181)'; %extract 13-30th June, 2007

savedir_var='/home/disk/eos8/d.grosvenor/saved_data_L2/';

switch operation
    case 'load'
        filepath = [savedir_var 'actual_L3_for_13-30Jun_2007_20120509T164441.mat'];
    case 'extract and save';
        filepath = [savedir_var 'actual_L3_for_13-30Jun_2007_' datestr(now,30) '.mat'];              
end

istring=1;
clear modis_var field_string indices_str

 modis_var{istring}='Cloud_Fraction_Liquid'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
modis_var{istring}='Cloud_Fraction_Liquid_Pixel_Counts'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=21)  %N.B. this is exactly the same as totN
modis_var{istring}='Cloud_Fraction_Combined'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=440)

modis_var{istring}='Cloud_Effective_Radius_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)

modis_var{istring}='Cloud_Optical_Thickness_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=337)

%modis_var{istring}='Sensor_Zenith_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 


%clear modis_var; istring=1;
modis_var{istring}='Cloud_Top_Temperature_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 

%modis_var{istring}='Cloud_Top_Pressure_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Pressure_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Pressure_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
        
modis_var{istring}='Solar_Zenith_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 

for ivar=1:length(modis_var)
    indices_str{ivar}=indices_str2;
end

modis_var{istring}='MLAT'; field_string{istring}=''; indices_str{istring}=''; istring=istring+1; 
modis_var{istring}='MLON'; field_string{istring}=''; indices_str{istring}=''; istring=istring+1;



    

        for ivar=1:length(modis_var)
           

            switch operation
                case 'extract and save'                    
                     fprintf(1,'\nCopying & saving variable %d of %d',ivar,length(modis_var));

                    eval_str = [modis_var{ivar} '_L32007_Arctic' field_string{ivar} '=' modis_var{ivar} field_string{ivar} indices_str{ivar} ';'];
                    eval(eval_str);

                    if ivar>1
                        app_str='''-APPEND'',';
                    else
                        app_str='';
                    end
                    
                    eval_str = ['save(''' filepath ''',''' modis_var{ivar} '_L32007_Arctic'',' app_str '''-V7.3'');'];
                    eval(eval_str);

                case 'load'
                    fprintf(1,'\nLoading variable %d of %d',ivar,length(modis_var));
                     
                    app_str='';
                    eval_str = ['load(''' filepath ''',''' modis_var{ivar} '_L32007_Arctic''' app_str ');'];
                    eval(eval_str);
            end



        end
        
        
%          switch operation
%                 case 'extract and save'  
%                     eval_str = ['save(''' filepath ''',''MLAT_L32007_Arctic'',' app_str '''-V7.3'');'];
%                     eval(eval_str);                    
%          end
%                     
                    
                    

