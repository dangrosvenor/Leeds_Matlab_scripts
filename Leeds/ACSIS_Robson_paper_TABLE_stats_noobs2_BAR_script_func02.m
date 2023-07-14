%Need to make sure we include the exp part - think this applies to the rest
%of the script too?
model_str=model; %default
switch model
    case 'ukesm'
        model_str='';
end

eval(['d' var '_' model '_' period_str ' = table_vals.' var model_str 'trend' period_str...
    ' * 10.^table_vals.' var model_str 'trend' period_str 'exp;']);

switch var
    case 'TS'
        eval(['d' var 'global_' model '_' period_str ' = table_vals_region0.' var model_str 'trend' period_str...
    ' * 10.^table_vals_region0.' var model_str 'trend' period_str 'exp;']);
        
end
    
