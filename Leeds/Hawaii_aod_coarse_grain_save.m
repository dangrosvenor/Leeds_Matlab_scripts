function Hawaii_aod_coarse_grain_save(save_file,dat,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var)

%% Save AOD data in .mat file
%aod_save_file = [UM_base_dir um_case_PD '/aod550_total.mat'];


save_var_list = {'time_out_var','dat','gcm_Plat2D_UM','gcm_Plon2D_UM','gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM'};

M_coarse_grain=1;
N_coarse_grain=1;
save_str = '_1x1';
%eval_str=['clear dat' save_str]; eval(eval_str);
for it=1:size(dat,3)
    eval_str=['dat' save_str '(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);']; eval(eval_str);
end
eval_str=['gcm_Plat2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['gcm_Plon2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['[gcm_Plat2D_edges' save_str ',gcm_Plon2D_edges' save_str '] = get_edges_lat_lon(gcm_Plat2D' save_str ',gcm_Plon2D' save_str ');']; eval(eval_str);
save_var_list{end+1} = ['dat' save_str];
save_var_list{end+1} = ['gcm_Plat2D' save_str];
save_var_list{end+1} = ['gcm_Plon2D' save_str];
save_var_list{end+1} = ['gcm_Plat2D_edges' save_str];
save_var_list{end+1} = ['gcm_Plon2D_edges' save_str];

%Coarse grain the data to approximately match MODIS AOD (10km).
%Model resolution is 4km so do both 2x2 and 3x3 (8 and 12km)
%Could also do 2x3 and 3x2 since this would match the km^2 value of a
%10x10km pixel since 2*4 * 3*4 = 96 km^2
M_coarse_grain=2;
N_coarse_grain=2;
save_str = '_2x2';
%eval_str=['clear dat' save_str]; eval(eval_str);
for it=1:size(dat,3)
    eval_str=['dat' save_str '(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);']; eval(eval_str);
end
eval_str=['gcm_Plat2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['gcm_Plon2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['[gcm_Plat2D_edges' save_str ',gcm_Plon2D_edges' save_str '] = get_edges_lat_lon(gcm_Plat2D' save_str ',gcm_Plon2D' save_str ');']; eval(eval_str);
save_var_list{end+1} = ['dat' save_str];
save_var_list{end+1} = ['gcm_Plat2D' save_str];
save_var_list{end+1} = ['gcm_Plon2D' save_str];
save_var_list{end+1} = ['gcm_Plat2D_edges' save_str];
save_var_list{end+1} = ['gcm_Plon2D_edges' save_str];


M_coarse_grain=3;
N_coarse_grain=3;
save_str = '_3x3';
%eval_str=['clear dat' save_str]; eval(eval_str);
for it=1:size(dat,3)
    eval_str=['dat' save_str '(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);']; eval(eval_str);
end
eval_str=['gcm_Plat2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['gcm_Plon2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['[gcm_Plat2D_edges' save_str ',gcm_Plon2D_edges' save_str '] = get_edges_lat_lon(gcm_Plat2D' save_str ',gcm_Plon2D' save_str ');']; eval(eval_str);
save_var_list{end+1} = ['dat' save_str];
save_var_list{end+1} = ['gcm_Plat2D' save_str];
save_var_list{end+1} = ['gcm_Plon2D' save_str];
save_var_list{end+1} = ['gcm_Plat2D_edges' save_str];
save_var_list{end+1} = ['gcm_Plon2D_edges' save_str];

M_coarse_grain=2;
N_coarse_grain=3;
save_str = '_2x3';
%eval_str=['clear dat' save_str]; eval(eval_str);
for it=1:size(dat,3)
    eval_str=['dat' save_str '(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);']; eval(eval_str);
end
eval_str=['gcm_Plat2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['gcm_Plon2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['[gcm_Plat2D_edges' save_str ',gcm_Plon2D_edges' save_str '] = get_edges_lat_lon(gcm_Plat2D' save_str ',gcm_Plon2D' save_str ');']; eval(eval_str);
save_var_list{end+1} = ['dat' save_str];
save_var_list{end+1} = ['gcm_Plat2D' save_str];
save_var_list{end+1} = ['gcm_Plon2D' save_str];
save_var_list{end+1} = ['gcm_Plat2D_edges' save_str];
save_var_list{end+1} = ['gcm_Plon2D_edges' save_str];

M_coarse_grain=3;
N_coarse_grain=2;
save_str = '_3x2';
eval_str=['clear dat' save_str]; eval(eval_str);
for it=1:size(dat,3)
    eval_str=['dat' save_str '(:,:,it) = reduce_matrix_subsample_mean(dat(:,:,it),M_coarse_grain,N_coarse_grain);']; eval(eval_str);
end
eval_str=['gcm_Plat2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['gcm_Plon2D' save_str ' = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);']; eval(eval_str);
eval_str=['[gcm_Plat2D_edges' save_str ',gcm_Plon2D_edges' save_str '] = get_edges_lat_lon(gcm_Plat2D' save_str ',gcm_Plon2D' save_str ');']; eval(eval_str);
save_var_list{end+1} = ['dat' save_str];
save_var_list{end+1} = ['gcm_Plat2D' save_str];
save_var_list{end+1} = ['gcm_Plon2D' save_str];
save_var_list{end+1} = ['gcm_Plat2D_edges' save_str];
save_var_list{end+1} = ['gcm_Plon2D_edges' save_str];

    
%     %save(save_file,'-V7.3','dat','gcm_Plat2D_UM','gcm_Plon2D_UM','gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM','time_out_var');
% eval_str=['save(save_file,''-V7.3'',''' save_var_list{1} ''');']; eval(eval_str);
% for ivar=2:length(save_var_list)
%     eval_str=['save(save_file,''-APPEND'',''' save_var_list{ivar} ''');']; eval(eval_str);    
% end

save(save_file,'-V7.3'); %just save everything
    
display('Finished Hawaii_aod_coarse_grain_save');