function [dat_out,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_tmp,gcm_Plat2D,gcm_Plon2D,optional);

if exist('optional')==1
    %Convert all of the variable names in the input structure to actual names
    %for ease of use
    name_struc='optional'; %The name of the structure
    names = eval(['fieldnames(' name_struc ');']);
    for i=1:length(names)
        eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
        eval(eval_str);
    end
end

if ~exist('screen_type_optional')
    screen_type_optional = 'ocean only';    %default to screening out land.    
end

if ~exist('lmask_in_optional')
    switch screen_type_optional
        case 'land+ocean'
            lmask_out = ones(size(gcm_Plat2D));
        otherwise
            lmask_type = 'UM'; %Makes it consistent with the model this way
            switch lmask_type
                case 'modis'
                    land_mask_MODIS=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
                    lmask_MODIS = flipdim(land_mask_MODIS.amsre_land_mask,1);
                    lmask_MODIS = lmask_MODIS + 1; %Make it ones where have ocean
                    lmask_out = griddata(land_mask_MODIS.gcm_Plat2D_AMSRE,land_mask_MODIS.gcm_Plon2D_AMSRE,lmask_MODIS,gcm_Plat2D,gcm_Plon2D,'nearest');
                    
                case 'UM'
                    load_type = 'merged netCDF';
                    var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
                    %um_case='u-bf666'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
                    %dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
                    um_case='u-cj725'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
                    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case '/glm/'];                    
                    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
                    % Re-grid the UM land mask to the obs land mask.
                    lmask_out = griddata(dat_global.gcm_Plat2D_UM,dat_global.gcm_Plon2D_UM,dat_global.dat,gcm_Plat2D,gcm_Plon2D,'nearest');
                    lmask_out(lmask_out==1)=NaN;
                    lmask_out = lmask_out + 1;
                    
            end                        
            
    end
else
    lmask_out = lmask_in_optional;
end

switch screen_type_optional
    case 'land+ocean'
        lmask = ones(size(lmask_out));
    case 'ocean only'
        lmask = lmask_out;
    case 'ocean only, no sea-ice'
        lmask = lmask_out;
        max_surf_dir = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/'];
        max_surf_name = [max_surf_dir 'max_surf_albedo_DAMIP_annual'];
        surf_mask = load(max_surf_name);
        %Re-grid the sufrace alebdo to the obs grid
        surf_out = griddata(surf_mask.gcm_Plat2D_UM,surf_mask.gcm_Plon2D_UM,surf_mask.max_surf,gcm_Plat2D,gcm_Plon2D,'nearest');        
        surf_alb_thresh = 0.2; %Max value of surface albedo allowed (higher assumes sea-ice)
        lmask(surf_out > surf_alb_thresh) = NaN;
        
    case 'land only'
        land_only_mask = ones(size(lmask_out));
        land_only_mask(lmask_out==1)=NaN;
        lmask = land_only_mask;
end

if length(size(dat_tmp))==2
    dat_out = squeeze(dat_tmp).*lmask;
else    
    lmask = repmat(lmask,[1 1 size(dat_tmp,1)]);
    lmask = permute(lmask,[3 1 2]);    
    
    dat_out = dat_tmp.*lmask;        
end
