fprintf(1,'\n loading L2 data');

clear modis_saved_file_list
istring=1;

savedir_var=['/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/'];
%aqua
modis_saved_file_list{istring}=[savedir_var 'L2_data_20110607T212113_AQUA_sf_correction']; istring=istring+1;
%terra
modis_saved_file_list{istring}=[savedir_var 'L2_data_20110607T174144_TERRA_sf_correction']; istring=istring+1;

%aqua summer NH
%modis_saved_file_list{istring}=[savedir_var 'L2_data_20110609T174900']; istring=istring+1;%
%terra summer NH
%modis_saved_file_list{istring}=[savedir_var 'L2_data_20110613T111336']; istring=istring+1;


%put the required variable names into modis_var{i}
set_MODIS_L2_variable_names

app_str=')';

for ifile=1:length(modis_saved_file_list)

    savevarname = modis_saved_file_list{ifile};

    for iread_modis=1:length(modis_var)

        fprintf(1,'\n loading %d out of %d',iread_modis,length(modis_var));

        MODIS_varname = [modis_var{iread_modis}];

        eval(['load(''' savevarname ''',''' MODIS_varname '_tim''' app_str]);
        
        if ifile==1 %no need to append for the first file
            eval([MODIS_varname '_tim_ALL = ' MODIS_varname '_tim;']);
        else %append the newly loaded data
            cat_dim = eval(['length(size(' MODIS_varname '_tim))']);
            eval([MODIS_varname '_tim_ALL = cat(cat_dim,' MODIS_varname '_tim_ALL,' MODIS_varname '_tim);']);
        end
        
        %on last file put data back into array name without the ALL for consistency
        if ifile==length(modis_saved_file_list) 
            eval([MODIS_varname '_tim = ' MODIS_varname '_tim_ALL;']);
        end
            

    end


    load(savevarname,'timLAT');
    load(savevarname,'timLON');


end


