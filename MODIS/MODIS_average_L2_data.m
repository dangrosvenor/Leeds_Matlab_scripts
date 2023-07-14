%works on an individual MODIS L2 swath to extract data for specified regions
%usually called from MODIS_process_multiple_L2_files


%lat_case = 'choose here';
lat_case = 'VOCALS cloud segments (Chris T)';

switch lat_case
    case 'choose here'
        %timLAT=[71+16.1/60];
        %timLON=[-156+48.3/60];

        timLAT=[0:2:90];
        timLON=[-156.75 -153 -149];
    case 'VOCALS cloud segments (Chris T)'
        VOCALS_cloud_leg_info  %(external script)
        
        time_match_cloud_legs_inds = time_match_cloud_legs{imod_read}(:);
        
        timLAT = lat_voc(time_match_cloud_legs_inds);
        timLON = lon_voc(time_match_cloud_legs_inds);
        
%        timLAT = lat_voc;
%        timLON = lon_voc;
end
    
eval(['Timeseries_L2_LON = timLON;']);
eval(['Timeseries_L2_LAT = timLAT;']);

%iLON=findheight_nearest(lon,timLON);
%iLAT=findheight_nearest(lat,timLAT);

%dlat=0.5;

clear inds dists
for iLON=1:length(timLON)
    switch lat_case
        case 'VOCALS cloud segments (Chris T)'
            LONS = timLON(iLON);
            LATS = timLAT(iLON);
        otherwise
            %make array of the longitude for all lats
            LONS = ones(size(timLAT))*timLON(iLON);
            LATS = timLAT

    end
    
    [ilat,ilon,dist_min]=getind_latlon_quick(lat,lon, LATS, LONS);
    %ilat,ilon now reference the lat,lon, etc. arrays
    %make 2D array of these indices (for each timLAT,timLON location)
    inds(:,iLON) = sub2ind(size(lat),ilat,ilon);
    dists(:,iLON) = dist_min;
end

dist_thresh = 10; %km - max distance allowed between the swath point and the desired
dist_thresh = 200;
%lat lon point

%find all the locations out of the timLAT,timLON locations that were
%within the swath (was a swath pixel closer then dist_thresh)
iyes=find(dists<=dist_thresh);
%so iyes indices are for the (timLAT,timLON) sized arrays

switch lat_case
    case 'VOCALS cloud segments (Chris T)'
        if imod_read==1
             clear VOCALS_leg_number_matches 
             nvoc_match_file = 0;
        end
            
        if length(iyes)>0 
           nvoc_match_file=nvoc_match_file+1;
           
           inds_voc = time_match_cloud_legs_inds(iyes);
           
           VOCALS_leg_number_matches(nvoc_match_file).time_inds = inds_voc;
           VOCALS_leg_number_matches(nvoc_match_file).filename = filename_h5;
           VOCALS_leg_number_matches(nvoc_match_file).dist_thresh = dist_thresh;
           VOCALS_leg_number_matches(nvoc_match_file).time_thresh = dtime_restrict;                      
           VOCALS_leg_number_matches(nvoc_match_file).leg_hours = hours_voc(inds_voc);
           VOCALS_leg_number_matches(nvoc_match_file).leg_mins = mins_voc(inds_voc);
           VOCALS_leg_number_matches(nvoc_match_file).sat_date_str = date_str;
           VOCALS_leg_number_matches(nvoc_match_file).sat_hour_str = modis_time_str(1:2);
           VOCALS_leg_number_matches(nvoc_match_file).sat_min_str = modis_time_str(3:4);  
           VOCALS_leg_number_matches(nvoc_match_file).distances = dists(iyes);
                                 
           times_flight = datenum(year_voc(inds_voc),month_voc(inds_voc),day_voc(inds_voc),hours_voc(inds_voc),mins_voc(inds_voc),secs_voc(inds_voc));
           times_swath = datenum([date_str ' ' modis_time_str(1:2) ':' modis_time_str(3:4)]);
            
            %difference between all the vocals cloud leg times in hours
            dtime = 24*60*(times_swath - times_flight);                                     
            VOCALS_leg_number_matches(nvoc_match_file).time_diff_mins = dtime;
            VOCALS_leg_number_matches(nvoc_match_file).time_diff_mins_txt = 'Time differences are sat time minus flight time';
        end
           
           
    otherwise

        for imodvar=1:length(modis_var)

            %clear the XXX_tim variable if is the first file
            if imod_read==1
                eval(['clear ' modis_var{imodvar} '_tim']);
            end

            eval_str = [modis_var{imodvar} '_tim'];

            %set all values to NaN as default
            eval([eval_str '(:,imod_read)=NaN*ones(size(inds(:)));']);

            %add data for the locations contained within the loaded L2 swath
            %inds(i) is index within the lat,lon,modis_var arrays for location i
            eval([eval_str '(iyes,imod_read)=' modis_var{imodvar} '(inds(iyes));']);

            %reshape from a [N,nMOD_av] to a [length(timLAT),length(timLON),nMOD_av]
            %array
            if imod_read==nMOD_av
                eval([eval_str '=reshape(' eval_str ',size(inds,1),size(inds,2),imod_read);']);
            end

        end




end



