%% 

%band_str='16'; %string to define which Re band to use
%band_str='21';
%band_str='37';

fid = fopen('/home/disk/eos1/d.grosvenor/matlab/work/MODIS/Nd_attr_table.txt','wt');
fid2 = fopen('/home/disk/eos1/d.grosvenor/matlab/work/MODIS/sza_means_stds.txt','wt');

%load the low sens, all sigCTT data for low and high sza
save_pdfs_filename='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs.mat';
load(save_pdfs_filename)
run_str = 'low \\sens';
run_str2 = 'All $\\sigma_{CTT}$, low \\sens';
make_sza_rel_increases_table_run

save_pdfs_filename='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs_highsens.mat';
load(save_pdfs_filename)
run_str = 'high \\sens';
run_str2 = 'All $\\sigma_{CTT}$, high \\sens';
make_sza_rel_increases_table_run

save_pdfs_filename='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs_low_sigCTT.mat';
load(save_pdfs_filename)
run_str = 'Low $\\sigma_{CTT}$';
run_str2 = 'Low $\\sigma_{CTT}$, low \\sens';
make_sza_rel_increases_table_run

save_pdfs_filename='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/saved_Arctic_tau_re_Nd_PDFs_high_sigCTT.mat';
load(save_pdfs_filename)
run_str = 'High $\\sigma_{CTT}$';
run_str2 = 'High $\\sigma_{CTT}$, low \\sens';
make_sza_rel_increases_table_run

fclose(fid);
fclose(fid2);

                            
                            
                            
                            
                            