
save_file='/home/disk/eos10/d.grosvenor/UM/Rosenfeld_analysis/ACSIS_global_mean_dSW_dNd.mat';
nudged=load(save_file,'Nd_PD_ALL','Nd_PI_ALL','SW_TOA_PD_mean','SW_TOA_PI_mean','dSW_TOA','dNd');
%region24_SW_trends{1}

dat_ukesm.dat_annual = dNd;
