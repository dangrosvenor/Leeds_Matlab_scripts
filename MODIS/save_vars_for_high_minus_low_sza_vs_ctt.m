tag = 'highSZA';
%tag = 'lowSZA';
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt']; %re
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Nd']; %Nd
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau']; %tau

filedir_savevars = '/home/disk/eos1/d.grosvenor/';
filename_savevars = [filedir_savevars filename '.mat'];

save_or_load = 'save';


ivar=1; clear vars_to_save


vars_to_save{ivar}='xdat'; ivar=ivar+1;

vars_to_save{ivar}='ydat'; ivar=ivar+1;
vars_to_save{ivar}='errordatU'; ivar=ivar+1;
vars_to_save{ivar}='errordatL'; ivar=ivar+1;
vars_to_save{ivar}='thresh_str'; ivar=ivar+1;

%run the save script.
save_vars_mat