function save_vars_for_Florent(savefile,Psave,Nsave,ind,varname)

%P=Psave{ind};
%N=Nsave{ind};

eval([varname ' = Psave{ind};']);
eval([varname '_Ndays = Nsave{ind};']);


save(savefile,varname,[varname '_Ndays'],'-APPEND');
