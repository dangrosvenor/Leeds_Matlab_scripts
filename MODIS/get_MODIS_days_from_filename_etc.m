clear days
for imod_read=1:nMOD_av
    imr=files_to_read(imod_read);
    
    file_name_h5 = files_mod(imr).name;
    iday=findstr(file_name_h5,'.A');
    modis_year_str=file_name_h5(iday+2:iday+5);
    modis_day_str=file_name_h5(iday+6:iday+8);
    
    days(imod_read) = str2num(modis_day_str);
    
end