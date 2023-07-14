savefile = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs_highsens.mat';
savefile = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs.mat';
savefile = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs_low_sigCTT.mat';
savefile = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs_high_sigCTT.mat';

isave=0;
%% Tau
ioverride_122=1;
multi_plot_case_122 = 'Tau vs Tau';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

%% Reff
ioverride_122=1;
multi_plot_case_122 = 'Re vs Re';
ioverride_122_save=1;
xvar_save_str = 'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

ioverride_122=1;
multi_plot_case_122 = 'Re vs Re';
ioverride_122_save=1;
xvar_save_str = 'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

ioverride_122=1;
multi_plot_case_122 = 'Re vs Re';
ioverride_122_save=1;
xvar_save_str = 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

%%Nd
ioverride_122=1;
multi_plot_case_122 = 'Nd vs Nd';
ioverride_122_save=1;
xvar_save_str = 'Nd_{1.6} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

ioverride_122=1;
multi_plot_case_122 = 'Nd vs Nd';
ioverride_122_save=1;
xvar_save_str = 'Nd from grid vals timeseries3';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

ioverride_122=1;
multi_plot_case_122 = 'Nd vs Nd';
ioverride_122_save=1;
xvar_save_str = 'Nd_{3.7} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7';
man_choose_water_graph=1;
graph=131;
waterVapourMay2005
save_TauRe_attribution_PDF_data

%save the data to file
isave=1;
save_pdfs_filename = savefile;
save_TauRe_attribution_PDF_data


