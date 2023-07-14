function [dat_out,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM)

nT = length(time_inds);

switch load_type
    case 'mat'
        dat_out = NaN * ones([size(gcm_Plat2D_UM,1) size(gcm_Plat2D_UM,2) nT]);
        it=0;
        for it_global_diff=time_inds
            it=it+1;
            dat_out(:,:,it) = eval(['dat_global.' var_UM '_ALL{it_global_diff};']);
        end

    case {'merged netCDF','cam6'}
        if length(size(dat_global.dat))>2
            dat_out = squeeze(dat_global.dat(time_inds,:,:));
        else
           dat_out =  squeeze(dat_global.dat(:,:));
        end
        if length(size(dat_out))==3
            dat_out = permute(dat_out,[2 3 1]);
        end
end

