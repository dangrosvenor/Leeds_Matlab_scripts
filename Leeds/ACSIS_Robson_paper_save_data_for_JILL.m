
save_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/totCF_data_for_Jill.mat']

clear y time

year = trend_dat_box_ens{1,1,1}.x;

for iens=1:size(trend_dat_box_ens,3)    
    y(:,iens) = trend_dat_box_ens{1,1,iens}.y;        
end

save(save_file,'year','y');
mat2nc_Dan(save_file,[save_file '.nc']);