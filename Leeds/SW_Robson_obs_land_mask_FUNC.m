function [dat_out,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_tmp,gcm_Plat2D,gcm_Plon2D,lmask_in_optional,screen_type_optional)

if ~exist('land_ocean')
    land_ocean = 'ocean only';    %default to screening out land.    
end

if ~exist('lmask_in_optional')
    switch screen_type_optional
        case 'land+ocean'
            lmask_out = ones(size(gcm_Plat2D));
        otherwise
            land_mask_MODIS=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
            lmask_MODIS = flipdim(land_mask_MODIS.amsre_land_mask,1);
            lmask_MODIS = lmask_MODIS + 1; %Make it ones where have ocean
            lmask_out = griddata(land_mask_MODIS.gcm_Plat2D_AMSRE,land_mask_MODIS.gcm_Plon2D_AMSRE,lmask_MODIS,gcm_Plat2D,gcm_Plon2D,'nearest');
    end
else
    lmask_out = lmask_in_optional;
end

switch land_ocean
    case 'land+ocean'
        lmask = ones(size(lmask_out));
    case 'ocean only'
        lmask = lmask_out;
    case 'land only'
        land_only_mask = ones(size(lmask_out));
        land_only_mask(lmask_out==1)=NaN;
        lmask = land_only_mask;
end

if length(size(dat_tmp))==2
    dat_out = squeeze(dat_tmp).*lmask;
else
    
    for it=1:size(dat_tmp,1)
        dat_out(it,:,:) = squeeze(dat_tmp(it,:,:)).*lmask;
    end
    
end
