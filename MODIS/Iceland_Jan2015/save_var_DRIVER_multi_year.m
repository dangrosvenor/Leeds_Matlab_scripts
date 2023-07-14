function save_var_DRIVER_multi_year(save_file_DRIVER,inew_file,dat,str)

eval([str '=dat;']);

if inew_file==1
    save(save_file_DRIVER,str);
else
    save(save_file_DRIVER,str,'-APPEND');    
end