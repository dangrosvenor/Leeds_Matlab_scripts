function save_vars_mat_func(filename,var_name,dat,iappend)

if iappend==1
    app_str=[',''-APPEND'''];
else
    app_str='';
end

eval_str=[var_name '=dat;'];
eval(eval_str);

i=findstr('.',var_name);
if i>0
   var_name=var_name(1:i-1); 
end

eval_str = ['save(filename,''' var_name ''',''-V7.3''' app_str ')'];
eval(eval_str);
