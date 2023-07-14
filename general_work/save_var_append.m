function save_var_append(file,var_in,var_out)

eval([var_out '=var_in;']); %in case want to change names

save(file,var_out,'-APPEND','-v7.3');