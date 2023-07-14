% The time that we want to track back to.
time_end = datenum('13-Nov-2008 12:15:00');

% -- Load the file with the saved winds in

% Using nblock_lon = 75; nblock_lat = 42;
%load_file_track = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_tracked_winds_20151127T075942.mat';
% New file using the differences to the next GOES image even if it is an
% hour (but ignorning no.77). Not dividing by two here either.
load_file_track = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_tracked_winds_20151130T042813.mat';

load(load_file_track);

% run regions_of_joint_histo_space_12Nov2008_case
% in order to choose a region in the LWP-Nd phase space to track points for
% Or load in a .mat file with the lats and lons in
%load_file_roi = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_roi_Joint_Histo_phase_space_20151130T034948.mat';

%High Nd, low LWP
load_file_roi = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_roi_Joint_Histo_phase_space_20151201T032516.mat';
load(load_file_roi);

%% Load the GOES multi file
file_goes_multi_whole_dom = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat';
load(file_goes_multi_whole_dom);

%Could use inpaint_nans instead of this
%ilat_no_nan = 2:size(lat_track,1)-1;
%ilon_no_nan = 2:size(lat_track,2)-1;

lat_no_nan = inpaint_nans(lat_track);
lon_no_nan = inpaint_nans(lon_track);

%same for the full array - takes a while, so have added it to the .mat file
%gcm_Plat2D_GOES_nonan = inpaint_nans(gcm_Plat2D_GOES);
%gcm_Plon2D_GOES_nonan = inpaint_nans(gcm_Plon2D_GOES);

% So, now have the offsets in terms of indices that the features should
% move by (forwards in time) from one GOES image to the next. 
% Times are listed in time_offxy_save{i}, indices within the _multi array
% as it_inds_offxy 


%% Find the index for the time from which the roi was chosen for the
%% offsets
istart = find(time_offxy_save==gcm_time_matlab_GOES_roi);
iend = find(time_offxy_save==time_end);
% So the offset values from the index before should take us back to the
% position before.


    
  clear lat_new lon_new lat_old lon_old
%    LON_LAT_track = [-77.1,-20.6]; %itest=itest+1;
%    it_inds = [67:75 78:82];
     if istart>iend
        time_direc = 'backwards';
        it_inds_loop = [istart:-1:iend+1];
        it_inds = [istart:-1:iend];        
     else
         time_direc = 'forwards';
         it_inds_loop = [istart:iend-1];   %so we finish at the correct time and not one further on         
         it_inds = [istart:iend];   %so we finish at the correct time and not one further on
     end
     
   lat_new = NaN*ones([length(it_inds_loop)+1 length(lat_in_roi)]);
   lon_new = NaN*ones([length(it_inds_loop)+1 length(lat_in_roi)]);
   lat_old = NaN*ones([length(it_inds_loop)+1 length(lat_in_roi)]);
   lon_old = NaN*ones([length(it_inds_loop)+1 length(lat_in_roi)]);
   LWP_track = NaN*ones([length(it_inds_loop)+1 length(lat_in_roi)]);   
   Nd_track = NaN*ones([length(it_inds_loop)+1 length(lat_in_roi)]);      
     
   iloop=0;
   for it0=it_inds_loop
        iloop=iloop+1
        
        
        switch time_direc
            case 'forwards'
                ind_save = it_inds(iloop);
                off_x = off_x_LWP_save{ind_save};
                off_y = off_y_LWP_save{ind_save};
            case 'backwards'
                %Use the index before (in terms of time) and make the
                %offsets negative since we are moving backwards
                ind_save = it_inds(iloop+1);
                off_x = -off_x_LWP_save{ind_save};
                off_y = -off_y_LWP_save{ind_save};

        end
        ind_multi = it_inds_offxy(ind_save);  %indices for the full *multi* array


        
%% Write in the starting values
        if it0==it_inds_loop(1)
            ind_01 = it_inds_offxy(it_inds(1)); %First index
            LWP_track(iloop,:) = griddata(gcm_Plat2D_GOES_nonan,gcm_Plon2D_GOES_nonan,goes_LWP_multi{ind_01},lat_in_roi,lon_in_roi);
            Nd_track(iloop,:) = griddata(gcm_Plat2D_GOES_nonan,gcm_Plon2D_GOES_nonan,goes_Nd_multi{ind_01},lat_in_roi,lon_in_roi);
            lat_new(iloop,:) = lat_in_roi;
            lon_new(iloop,:) = lon_in_roi;            
        end
            


        
        for iroi=1:length(lat_in_roi)
            
            if it0==it_inds_loop(1)
                LON_LAT_track_in = [ lon_in_roi(iroi) lat_in_roi(iroi) ];
            else
                LON_LAT_track_in = [lon_new(iloop-1,iroi) lat_new(iloop-1,iroi)];
            end
        
          
       
            [lat_new(iloop+1,iroi),lon_new(iloop+1,iroi),lat_old(iloop+1,iroi),lon_old(iloop+1,iroi)] = track_motion_GOES_func(LON_LAT_track_in,lat_no_nan,lon_no_nan,off_x,off_y,gcm_Plat2D_GOES,gcm_Plon2D_GOES);
            
        end
        
 % Is much quicker to do lots of points at once than one at a time using
   % griddata
   LWP_track(iloop+1,:) = griddata(gcm_Plat2D_GOES_nonan,gcm_Plon2D_GOES_nonan,goes_LWP_multi{ind_multi},lat_new(iloop+1,:),lon_new(iloop+1,:));
   Nd_track(iloop+1,:) = griddata(gcm_Plat2D_GOES_nonan,gcm_Plon2D_GOES_nonan,goes_Nd_multi{ind_multi},lat_new(iloop+1,:),lon_new(iloop+1,:));        
 
    
   end

  
   LWP_track_mean=meanNoNan(LWP_track,2);
   Nd_track_mean=meanNoNan(Nd_track,2);
   
   


run_feature_track_test=0;

if run_feature_track_test==1

    itest=0;
    clear lat_new

    %Features to mark to test tracking

    % Interpolate the x and y offset (no. pixels)
    %can't have NaNs for griddata, so chop off the first column - need a more
    %consistent way to do this
    %Function to find the new lat and lon
    % LON_LAT_track = [-74.13,-20.47]; itest=itest+1;
    % [lat_new(itest),lon_new(itest),lat_old(itest),lon_old(itest)] = track_motion_GOES_func(LON_LAT_track,lat_track(ilat_no_nan,ilon_no_nan),lon_track(ilat_no_nan,ilon_no_nan),off_x_LWP(ilat_no_nan,ilon_no_nan),off_y_LWP(ilat_no_nan,ilon_no_nan),gcm_Plat2D_GOES,gcm_Plon2D_GOES);
    %
    % LON_LAT_track = [-76.35,-18.6]; itest=itest+1;
    % [lat_new(itest),lon_new(itest),lat_old(itest),lon_old(itest)] = track_motion_GOES_func(LON_LAT_track,lat_track(ilat_no_nan,ilon_no_nan),lon_track(ilat_no_nan,ilon_no_nan),off_x_LWP(ilat_no_nan,ilon_no_nan),off_y_LWP(ilat_no_nan,ilon_no_nan),gcm_Plat2D_GOES,gcm_Plon2D_GOES);
    %
    % LON_LAT_track = [-78.35,-18.55]; itest=itest+1;
    % [lat_new(itest),lon_new(itest),lat_old(itest),lon_old(itest)] = track_motion_GOES_func(LON_LAT_track,lat_track(ilat_no_nan,ilon_no_nan),lon_track(ilat_no_nan,ilon_no_nan),off_x_LWP(ilat_no_nan,ilon_no_nan),off_y_LWP(ilat_no_nan,ilon_no_nan),gcm_Plat2D_GOES,gcm_Plon2D_GOES);
    %
    % LON_LAT_track = [-77.2,-21.3]; itest=itest+1;
    % [lat_new(itest),lon_new(itest),lat_old(itest),lon_old(itest)] = track_motion_GOES_func(LON_LAT_track,lat_track(ilat_no_nan,ilon_no_nan),lon_track(ilat_no_nan,ilon_no_nan),off_x_LWP(ilat_no_nan,ilon_no_nan),off_y_LWP(ilat_no_nan,ilon_no_nan),gcm_Plat2D_GOES,gcm_Plon2D_GOES);


    %LON_LAT_track = [-77.1,-20.6]; itest=itest+1;
    %[lat_new(itest),lon_new(itest),lat_old(itest),lon_old(itest)] = track_motion_GOES_func(LON_LAT_track,lat_track(ilat_no_nan,ilon_no_nan),lon_track(ilat_no_nan,ilon_no_nan),off_x_LWP(ilat_no_nan,ilon_no_nan),off_y_LWP(ilat_no_nan,ilon_no_nan),gcm_Plat2D_GOES,gcm_Plon2D_GOES);

    clear lat_new lon_new lat_old lon_old
    LON_LAT_track = [-77.1,-20.6]; %itest=itest+1;
    it_inds = [67:75 78:82];
    iloop=0;
    for it0=it_inds
        iloop=iloop+1;

        if it0==it_inds(1)
            LON_LAT_track_in = LON_LAT_track;
        else
            LON_LAT_track_in = [lon_new(it_inds(iloop-1)) lat_new(it_inds(iloop-1))];
        end
        [lat_new(it0),lon_new(it0),lat_old(it0),lon_old(it0)] = track_motion_GOES_func(LON_LAT_track_in,lat_track(ilat_no_nan,ilon_no_nan),lon_track(ilat_no_nan,ilon_no_nan),off_x_LWP_save{it0}(ilat_no_nan,ilon_no_nan),off_y_LWP_save{it0}(ilat_no_nan,ilon_no_nan),gcm_Plat2D_GOES,gcm_Plon2D_GOES);

        ioverride_goes_multi=1;
        data_source = 'goes_multi_array'; imulti_goes=it0;
        POC_26Oct2008_GOES_general_maps_20141125T032943

        if it0==it_inds(1)
            m_plot(LON_LAT_track(1),LON_LAT_track(2),'ko','markersize',9,'markerfacecolor','w');
            m_plot(LON_LAT_track(1),LON_LAT_track(2),'k+','markersize',9);


        else

            %    if it0==it_inds(end)
            m_plot(LON_LAT_track_in(1),LON_LAT_track_in(2),'ko','markersize',9,'markerfacecolor','w');
            m_plot(LON_LAT_track_in(1),LON_LAT_track_in(2),'kx','markersize',9);
            %    end
        end
    end





    %for itest=1:length(lat_new)
    m_plot(lon_old(itest),lat_old(itest),'ko','markersize',9,'markerfacecolor','w');
    m_plot(lon_old(itest),lat_old(itest),'k+','markersize',9);
    m_plot(lon_new(itest),lat_new(itest),'ko','markersize',9,'markerfacecolor','w');
    m_plot(lon_new(itest),lat_new(itest),'kx','markersize',9);
    %end

    itest=it1;
    ioverride_goes_multi=1;
    data_source = 'goes_multi_array'; imulti_goes=it1;
    POC_26Oct2008_GOES_general_maps_20141125T032943

    %for itest=1:length(lat_new)
    m_plot(lon_old(itest),lat_old(itest),'ko','markersize',9,'markerfacecolor','w');
    m_plot(lon_old(itest),lat_old(itest),'k+','markersize',9);
    m_plot(lon_new(itest),lat_new(itest),'ko','markersize',9,'markerfacecolor','w');
    m_plot(lon_new(itest),lat_new(itest),'kx','markersize',9);
    %end

end