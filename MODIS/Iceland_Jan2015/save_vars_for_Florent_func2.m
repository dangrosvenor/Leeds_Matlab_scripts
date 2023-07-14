function save_vars_for_Florent(savefile,Psave,Nsave,Nd_std_dev_save,ind,varname)

%P=Psave{ind};
%N=Nsave{ind};

eval([varname ' = Psave{ind};']);
eval([varname '_Ndays = Nsave{ind};']);
eval([varname '_std_dev = Nd_std_dev_save{ind};']);

save(savefile,varname,[varname '_Ndays'],[varname '_std_dev'],'-APPEND');
