Y_mean_file_1 = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved_cw_1e-6.mat'; icalc_21um=1; %Smaller test value of 1e-6
%test using 1.81, but v3 of ODran's data
%Y_mean_file_1 = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved_cw_1.81e-6.mat'; icalc_21um=0;
Y_mean_file_1pt8 = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved.mat'; %original cw=1.8e-6 value





dat_1 = load(Y_mean_file_1);
dat_1pt8 = load(Y_mean_file_1pt8);


Y_mean_file_1 = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved_reff_param_cw_1e-6.mat'; %Smaller test value of 1e-6
%test using 1.81, but v3 of ODran's data
%Y_mean_file_1 = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved_reff_param_cw_1.81e-6.mat'; icalc_21um=0;
Y_mean_file_1pt8 = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Y_mean_Odran_saved_reff_param.mat'; %original cw=1.8e-6 value

dat_reff_1 = load(Y_mean_file_1);
dat_reff_1pt8 = load(Y_mean_file_1pt8);

if icalc_21um==1

mean_tau_diff_21 = 100 * mean( (dat_1.Y_mean_all_Nd_21_mum - dat_1pt8.Y_mean_all_Nd_21_mum) ./ dat_1pt8.Y_mean_all_Nd_21_mum)
max_tau_diff_21 = 100 * max(abs( (dat_1.Y_mean_all_Nd_21_mum - dat_1pt8.Y_mean_all_Nd_21_mum) ./ dat_1pt8.Y_mean_all_Nd_21_mum ))

mean_reff_diff_21 = 100 * mean( (dat_reff_1.Y_mean_all_Nd_21_mum - dat_reff_1pt8.Y_mean_all_Nd_21_mum) ./ dat_reff_1pt8.Y_mean_all_Nd_21_mum)
max_reff_diff_21 = 100 * max(abs( (dat_reff_1.Y_mean_all_Nd_21_mum - dat_reff_1pt8.Y_mean_all_Nd_21_mum) ./ dat_reff_1pt8.Y_mean_all_Nd_21_mum ))

end

%3.7um
mean_tau_diff_37 = 100 * mean( (dat_1.Y_mean_all_Nd_37_mum - dat_1pt8.Y_mean_all_Nd_37_mum) ./ dat_1pt8.Y_mean_all_Nd_37_mum)
max_tau_diff_37 = 100 * max(abs( (dat_1.Y_mean_all_Nd_37_mum - dat_1pt8.Y_mean_all_Nd_37_mum) ./ dat_1pt8.Y_mean_all_Nd_37_mum ))

mean_reff_diff_37 = 100 * mean( (dat_reff_1.Y_mean_all_Nd_37_mum - dat_reff_1pt8.Y_mean_all_Nd_37_mum) ./ dat_reff_1pt8.Y_mean_all_Nd_37_mum)
max_reff_diff_37 = 100 * max(abs( (dat_reff_1.Y_mean_all_Nd_37_mum - dat_reff_1pt8.Y_mean_all_Nd_37_mum) ./ dat_reff_1pt8.Y_mean_all_Nd_37_mum ))