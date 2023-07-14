function [obs_monthly_box,obs_annual_box,obs_annual_map] =...
    ACSIS_Robson_paper_process_obs(obs_dat,lat2d,lon2d,LAT_val,LON_val,years_obs,area_obs,season)

ilat = find(lat2d(:,1)>LAT_val(1) & lat2d(:,1)<LAT_val(2));
ilon = find(lon2d(1,:)>LON_val(1) & lon2d(1,:)<LON_val(2));
area_tmp = area_obs(ilat,ilon);

for it=1:size(obs_dat,1)
    dat_tmp = squeeze(obs_dat(it,ilat,ilon));    
    obs_monthly_box(it) = meanNoNan(dat_tmp(:),1,'',0,1,area_tmp(:));
    %dat_tmp = meanNoNan(meanNoNan(dat.dat_annual_ens(:,it,ilat,ilon),4),2);
    %dat_annual_box_ukesm_ens_std(it) = std(dat_tmp); %std dev across the ensemble
end

switch season
    case 'Annual'
        
        istart=1; iend=12;
        for it=1:length(years_obs)
            obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
            obs_annual_map(it,:,:) = meanNoNan(obs_dat(istart:iend,:,:),1);
            istart=istart+12;
            iend=iend+12;
        end
        
    case 'DJF'
        
        %Seasonal means
        %DJF
        istart=0; iend=istart+2;
        for it=1:length(years_obs)
            if it==1
                obs_annual_box(it) = NaN;
                obs_annual_map(it,:,:) = NaN*ones(size(obs_dat(1,:,:))); %set to NaN until have supplied the previous year
                %Or just use Jan and Feb
            else
                obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
                obs_annual_map(it,:,:) = meanNoNan(obs_dat(istart:iend,:,:),1);
            end
            istart=istart+12;
            iend=iend+12;
        end
        
    case 'MAM'
        
        %MAM
        istart=3; iend=istart+2;
        for it=1:length(years_obs)
            obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
            obs_annual_map(it,:,:) = meanNoNan(obs_dat(istart:iend,:,:),1);
            istart=istart+12;
            iend=iend+12;
        end
        
    case 'JJA'
        
        %JJA
        istart=6; iend=istart+2;
        for it=1:length(years_obs)
            obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
            obs_annual_map(it,:,:) = meanNoNan(obs_dat(istart:iend,:,:),1);
            istart=istart+12;
            iend=iend+12;
        end
        
    case 'SON'
        
        %SON
        istart=9; iend=istart+2;
        for it=1:length(years_obs)
            obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
            obs_annual_map(it,:,:) = meanNoNan(obs_dat(istart:iend,:,:),1);
            istart=istart+12;
            iend=iend+12;
        end
        
    case 'JFM'
        
        %JFM
        istart=1; iend=istart+2;
        for it=1:length(years_obs)
            obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
            obs_annual_map(it,:,:) = meanNoNan(obs_dat(istart:iend,:,:),1);
            istart=istart+12;
            iend=iend+12;
        end        
        
        
end
