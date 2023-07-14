iload=0;
icalc_sza=0;

if iload==1
    load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat';
    load(load_file_goes);
end



if icalc_sza==1
    %calculate for the middle of the domain for each day and then adjust using
%     %longtiude
%     mean_lat = meanNoNan(gcm_Plat2D_GOES(:),1);
%     mean_lon = meanNoNan(gcm_Plon2D_GOES(:),1);
%     sza_centre = sun_pos(times_GOES_save,mean_lat,mean_lon);
%     
%     clear sza_ALL
%     for it=1:length(times_GOES_save)
%         it
%         for ipos=1:length(gcm_Plon2D_GOES(:))
%             sza_ALL{it}(ipos) = sun_pos(times_GOES_save(it),gcm_Plat2D_GOES(ipos),gcm_Plon2D_GOES(ipos));
%         end
%     end
    
    
    ilat=[1:floor(size(gcm_Plat2D_GOES,1)/10):size(gcm_Plat2D_GOES,1)];
    ilon=[1:floor(size(gcm_Plon2D_GOES,2)/10):size(gcm_Plon2D_GOES,2)];   
    
   
    for it=1:length(times_GOES_save)
        it
        clear sza_coarse
        for ilat2=1:length(ilat)
            for ilon2=1:length(ilon)
                sza_coarse(ilat2,ilon2) = sun_pos(times_GOES_save(it),gcm_Plat2D_GOES(ilat(ilat2),ilon(ilon2)),gcm_Plon2D_GOES(ilat(ilat2),ilon(ilon2)));
            end
        end        
        
        sza_coarse2 = inpaint_nans(sza_coarse);
        lat2 = inpaint_nans(gcm_Plat2D_GOES(ilat,ilon));
        lon2 = inpaint_nans(gcm_Plon2D_GOES(ilat,ilon));
%        sza_ALL{it} = griddata(gcm_Plat2D_GOES(ilat,ilon),gcm_Plon2D_GOES(ilat,ilon),sza_coarse,gcm_Plat2D_GOES,gcm_Plon2D_GOES);
        sza_ALL{it} = griddata(lat2,lon2,sza_coarse,gcm_Plat2D_GOES,gcm_Plon2D_GOES);   
%        inan = find( isnan(gcm_Plat2D_GOES) == 1 | isnan(gcm_Plon2D_GOES) == 1 );
%        sza_ALL{it}(inan) = NaN;
    end
    
    sza_notes = 'Calculated using find_max_Nd_GOES_12thNov_paper.m';
       
    
    save(load_file_goes,'sza_ALL','sza_notes','-V7.3','-APPEND');

end

prctiles = [0:5:100];
max_Nd = NaN * ones(length(goes_Nd_multi));
prc_Nd = NaN * ones([length(goes_Nd_multi) length(prctiles)]);
for it=1:length(times_GOES_save)
    max_dat = goes_Nd_multi{it};
    isza = find(sza_ALL{it} <= 65 );  
    if length(isza)>0
        max_Nd(it) = maxALL(max_dat(isza));
        prc_Nd(it,:) = prctile(max_dat(:),prctiles);
    end
end

