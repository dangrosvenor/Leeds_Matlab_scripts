%Can pick out regions on a 2D plot - e.g. in joint histogram phase space.
%Then can plot them on, say, a spatial map using for e.g. :-
% m_plot(lon_in_roi(:),lat_in_roi(:),'wx')

data_type = 'GOES';
data_type = 'UM';

save_matfile_pick_region = ['/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_roi_Joint_Histo_phase_space_' datestr(now,30) '.mat']

h=impoly(gca,[]); %starts an interactive tool to pick out a region by hand

%roi = getPosition(h); %returns [N 2] vector of x and y positions of the
%polygon - only works on newer Matlab versions. For older ones use :-
api = iptgetapi(h);
roi = api.getPosition();

%Split the region vertically into nYroi equal slices 

% nYroi=1000;
% yroi = linspace(min(roi(:,2)),max(roi(:,2)),nYroi);
% 
% for iroi=1:nYroi
%     inds_roi = find(X>
% end

% roi =
% 
%    32.6507  256.4010
%    62.4930  162.5201
%    95.9827  125.3220
%   135.7724  103.1803
%   158.6515   88.1240
%   156.6620   77.4959
%    76.7510   79.2673
%    61.1667   83.6956
%    35.9665  127.0934
%    27.6770  180.2335

%Store and save the time and the region used to make the joint histogram
LAT_val_DRIVER_roi = LAT_val_DRIVER;
LON_val_DRIVER_roi = LON_val_DRIVER;


        lat_roi =  eval( ['gcm_Plat2D_' gcm_str ';'] );
        lon_roi =  eval( ['gcm_Plon2D_' gcm_str ';'] );
        
        gcm_time_matlab_GOES_roi = eval( ['gcm_time_matlab_' gcm_str ';'] );

switch data_type
    case 'GOES'
        %
        LWP_roi=goes_LWP_save; %from pdf2D*
        Nd_roi=goes_Nd_save;        

    case 'UM'
        %
        LWP_roi = pdf2D_X_save; %from pdf2D*
        Nd_roi = pdf2D_Y_save;        

end
% 
in_roi = inpolygon(LWP_roi,Nd_roi,roi(:,1),roi(:,2)); %returns logical (zero or one) for all points in X and Y that are inside of on edge of roi
% %points in the polygon are now available as e.g. 
LWP_in_roi = LWP_roi(in_roi);
Nd_in_roi = Nd_roi(in_roi);
lat_in_roi =lat_roi(in_roi);
lon_in_roi = lon_roi(in_roi);

%can plto these on a map using e.g. 
%m_plot(lon_in_roi,lat_in_roi,'wo');



save(save_matfile_pick_region,'LWP_roi','Nd_roi','in_roi','LWP_in_roi','Nd_in_roi','lat_in_roi','lon_in_roi',...
    'gcm_time_matlab_GOES_roi','LAT_val_DRIVER_roi','LON_val_DRIVER_roi');

