function [mon_me_SO] = monthly_means_restrict_lat_lon(years,dat,gcm_Plat2D_UM,gcm_Plon2D_UM,thresh_LAT,thresh_LON,iscreen_land,land_mask_in,iarea_weight,obs_area)

if ~exist('iscreen_land') | iscreen_land==0
    iscreen_land=0;
    land_mask = ones(size(dat{1}));
else
    land_mask = repmat(land_mask_in,[1 1 size(dat{1},3)]);
end

if ~exist('iarea_weight')
    iarea_weight = 0;    
end

if iarea_weight==0
    obs_area = ones([size(dat{1},1) size(dat{2},2) ]);
end


ilat = find(gcm_Plat2D_UM(:,1)>thresh_LAT(1) & gcm_Plat2D_UM(:,1)<thresh_LAT(2));
ilon = find(gcm_Plon2D_UM(1,:)>thresh_LON(1) & gcm_Plon2D_UM(1,:)<thresh_LON(2));


for iy=1:length(years)  
    dat2 = NaN*ones(size(dat{iy}));
    dat2(ilat,ilon,:) = dat{iy}(ilat,ilon,:);
    dat2 = dat2 .* land_mask;
%     dat_tmp = permute(dat2(ilat,ilon,:),[3 1 2]);   
%     area_rep = repmat(obs_area(ilat,ilon),[1 1 size(dat2,3)]);
%     area_rep = permute(area_rep,[3 1 2]);
    
        
    %dat{iy}(inan_lat,inan_lon,:)=NaN;
    for im=1:12        
        %mon_me_SO{iy}(im) = meanNoNan(meanNoNan(dat2(:,:,im),1),1);
        dat_tmp = dat2(ilat,ilon,im);
        area_tmp = obs_area(ilat,ilon);
        mon_me_SO{iy}(im) = meanNoNan(dat_tmp(:),1,'',0,1,area_tmp(:));
    end
end
