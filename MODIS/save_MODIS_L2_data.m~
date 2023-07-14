%this saves the timeseries of the variables in modis_var
%(i.e. the name with _tim added)
fprintf(1,'\n saving L2 data in modis_var array...');


app_str=')';

for iread_modis=1:length(modis_var)
    
    fprintf(1,'\n saving %d out of %d',iread_modis,length(modis_var));

    MODIS_varname = modis_var{iread_modis};

    %double quotes tells matlab to insert a ' in the string (without
    %closing the quote for the sring prcessing)
    
    switch action
        case {'make mock L3 data','make mock L3 data daily'}
            eval(['save(''' savevarname '.mat'',''' MODIS_varname '''' app_str]);
        otherwise
            eval(['save(''' savevarname '.mat'',''' MODIS_varname '_tim''' app_str]);
    end
    
    app_str=',''-APPEND'')';

end


save([savevarname  '.mat'],'timLAT','-APPEND');
save([savevarname  '.mat'],'timLON','-APPEND');
%MLAT = timLAT;
%MLON = timLON;
save([savevarname  '.mat'],'MLAT','-APPEND');
save([savevarname '.mat'] ,'MLON','-APPEND');
%save(savevarname,'scantime_matlab_tim','-APPEND');

fprintf(1,'\nDone\n');


