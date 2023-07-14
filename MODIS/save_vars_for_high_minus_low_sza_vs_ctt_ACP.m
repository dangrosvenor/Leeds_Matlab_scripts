
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt']; %re
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Nd']; %Nd
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau']; %tau
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau_for_ACP']; %tau
filename = ['saved_vars_for_high_minus_low_sza_vs_gammaTau_Tau_for_ACP']; %tau

%filedir_savevars = '/home/disk/eos1/d.grosvenor/';
filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
filename_savevars = [filedir_savevars filename '.mat'];

save_or_load = 'save';


ivar=1; clear vars_to_save
vars_to_save{ivar}='xdat'; ivar=ivar+1;
vars_to_save{ivar}='ydat'; ivar=ivar+1;
vars_to_save{ivar}='errordatU'; ivar=ivar+1;
vars_to_save{ivar}='errordatL'; ivar=ivar+1;
vars_to_save{ivar}='thresh_str'; ivar=ivar+1;


xdat_bk=xdat;
ydat_bk=ydat;
errordatU_bk = errordatU;
errordatL_bk = errordatL;
%thresh_str_bk = thresh_str;

isza=1;
xdat = xdat_bk(isza);
ydat_lowSZA = ydat_bk(isza);
errordatU = errordatU_bk(isza);
errordatL = errordatL_bk(isza);
thresh_str = thresh_str;

tag = 'lowSZA';
%run the save script.
save_vars_mat

isza=2;
xdat = xdat_bk(isza);
ydat = ydat_bk(isza);
errordatU = errordatU_bk(isza);
errordatL = errordatL_bk(isza);
%thresh_str_highSZA = thresh_str;

tag = 'highSZA';
%run the save script.
save_vars_mat
