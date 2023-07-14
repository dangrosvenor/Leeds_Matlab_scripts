
eval( ['dat = table_vals.' var_str_tab  run_str 'trend' period_str ';'] );
eval( ['dat_un = table_vals.' var_str_tab  run_str 'trend' period_str 'un;']);

eval( [run_str ' = dat;';] );
eval( [run_str '_un = dat_un;']);

if (abs(dat) - dat_un > 0) dat_sig=1; else dat_sig=0; end
eval( [run_str '_sig = dat_sig;'] );

