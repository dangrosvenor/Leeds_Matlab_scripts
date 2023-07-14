%Looks like can't get later versions of Matlab to work stably, so will just
%write my own script using xcorr2

%Have done this - is called track_motion_Dan
%Have also made a .mat file using read_GOES_vocals_netcdf_files_MULTI.m
% that contains all of the GOES domain. Can then use sub-regions as
% requried. File is :- /home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat
%

file_goes_multi_whole_dom = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat';
load(file_goes_multi_whole_dom);

save_file_track = ['/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_tracked_winds_' datestr(now,30) '.mat'];


% The resolution of GOES appears to be around 2.4 x 4.4 km in lon x lat (?)
% So might want to use more pixels in lon direction.

% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
LAT_val = [-20.5 -17.5]; LON_val = [-78.75 -73.25];

LAT_val = [-22.70 -17.28]; LON_val =[-78.93 -73.08]; %FULL UM domain for 12th Nov

LAT_val =[]; LON_val=[];

if length(LAT_val)==0
    %This selects the full domain
    ilat=[1:size(gcm_Plat2D_GOES,1)];
    ilon=[1:size(gcm_Plat2D_GOES,2)];

else

    %[iregion_lin,iregion_lin_edges] = get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_GOES,gcm_Plon2D_GOES,gcm_Plat2D_edges_GOES,gcm_Plon2D_edges_GOES);

    %Even though the GOES lat and lon grids are not quite orthogonal and
    %regular will pick out a square region in terms of pixels
    %average the lats and lons over the rows/columns to get 1d vectors
    goes_lat = meanNoNan(gcm_Plat2D_GOES,2);
    goes_lon = meanNoNan(gcm_Plon2D_GOES,1);

    ilat = find( goes_lat>=LAT_val(1) & goes_lat<LAT_val(end) );
    ilon = find( goes_lon>=LON_val(1) & goes_lon<LON_val(end) );

    dlon=0.04; dlat=dlon;
    [lon_reg,lat_reg] = meshgrid([LON_val(1):dlon:LON_val(2)],[LAT_val(1):dlat:LAT_val(2)]);
    lat_reg = flipdim(lat_reg,1);


    sza = sun_pos(times_GOES_save,mean(LAT_val),mean(LON_val));

end


% SZA gets to 62.5 by no. 83
% and 64.1 for no. 67
% Is a gap between file 79 and file 80 (normally are every 30 mins, the gap
% there is 1 hour)
% Also between 74 and 75
% So, should only use 67-83
% datestr(times_GOES_save(67:83))
% 
% ans =
% 
% 13-Nov-2008 12:15:00 % 67
% 13-Nov-2008 12:45:00
% 13-Nov-2008 13:15:00
% 13-Nov-2008 13:45:00
% 13-Nov-2008 14:15:00
% 13-Nov-2008 14:45:00
% 13-Nov-2008 15:15:00
% 13-Nov-2008 15:45:00  % 74
% 13-Nov-2008 16:45:00  %hour gap - could do an average of the before and after
% 13-Nov-2008 17:15:00  % 76            % velocities
% 13-Nov-2008 17:45:00  %dodgy file - hour gap (no. 77)
% 13-Nov-2008 18:15:00
% 13-Nov-2008 18:45:00  % 79
% 13-Nov-2008 19:45:00  %hour gap
% 13-Nov-2008 20:15:00
% 13-Nov-2008 20:45:00
% 13-Nov-2008 21:15:00  % 83




%it0=75; it1=76; 
%it0=74; it1=75; 

%The files that are ok to use (and the range we want)
it_inds_offxy = [67:76 78:83]; %will run to one before the last index (and do the difference betweeen that and the last index)

i_it = 0;
clear off_y_LWP_save off_x_LWP_save time_offxy_save
for it0=[it_inds_offxy(1):it_inds_offxy(end-1)]  %[67:75 78:82]  %avoiding no. 77
    i_it = i_it + 1;
    fprintf(1,'it0=%d\n',it0);
    it1 = it_inds_offxy(i_it+1);  %so e.g. for no.76 should use 78 as the next one and void 77

max_dN_lon = 20; %25; %for 1.5 hrs (1.5*3600*10/1e3 = 54km to 10 m/s)
max_dN_lat = 10; %15;

nblock_lon = 50; %to split domain into about 5 equal regions
nblock_lat = 28;
%nblock_lon = 100; %
%nblock_lat = 56;
nblock_lon = 75; %to split domain into about 5 equal regions
nblock_lat = 42;



LWP_t0 = goes_LWP_multi{it0}(ilat,ilon); 
%LWP_t02 = goes_LWP_multi{it0}; %(ilat,ilon); 
%LWP_t02(isnan(LWP_t0))=0;
%a=NaN*ones(size(LWP_t0));
%a(iregion_lin)=0;
%LWP_t0 = LWP_t0 + a;

%lat_2d = gcm_Plat2D_GOES; lat_2d(isnan(lat_2d))=0;
%lon_2d = gcm_Plon2D_GOES; lon_2d(isnan(lon_2d))=0;
%LWP_t02 = griddata(lat_2d,lon_2d,LWP_t02,lat_reg,lon_reg);

LWP_t1 = goes_LWP_multi{it1}(ilat,ilon); %LWP_t1 = LWP_t1 + a;   
%LWP_t12 = goes_LWP_multi{it1}; %(ilat,ilon); 
%LWP_t12(isnan(LWP_t12))=0;
%LWP_t12 = griddata(lat_2d,lon_2d,LWP_t12,lat_reg,lon_reg);

% Teff_t0 = goes_Teff_multi{it0}(ilat,ilon); %Teff_t0 = Teff_t0 + a;    
% Teff_t1 = goes_Teff_multi{it1}(ilat,ilon); %Teff_t1 = Teff_t1 + a;  
% 
% Vis_t0 = goes_Vis_multi{it0}(ilat,ilon); %Vis_t0 = Vis_t0 + a;    
% Vis_t1 = goes_Vis_multi{it1}(ilat,ilon); %Vis_t1 = Vis_t1 + a; 
% 
% Nd_t0 = goes_Nd_multi{it0}(ilat,ilon); %Nd_t0 = Nd_t0 + a;    
% Nd_t1 = goes_Nd_multi{it1}(ilat,ilon); %Nd_t1 = Nd_t1 + a; 


% Getting outputs as [off_y,off_x...] so that y now corresponds to lat and
% x to lon
[off_y_LWP,off_x_LWP,max_corr_LWP,lat_track,lon_track]=track_motion_Dan(LWP_t0,LWP_t1,max_dN_lat,max_dN_lon,nblock_lat,nblock_lon,gcm_Plat2D_GOES(ilat,ilon),gcm_Plon2D_GOES(ilat,ilon));
%[off_y_Teff,off_x_Teff,max_corr_Teff,lat_track,lon_track]=track_motion_Dan(Teff_t0,Teff_t1,max_dN_lat,max_dN_lon,nblock_lat,nblock_lon,gcm_Plat2D_GOES,gcm_Plon2D_GOES);
%[off_y_Vis,off_x_Vis,max_corr_Vis,lat_track,lon_track]=track_motion_Dan(Vis_t0,Vis_t1,max_dN_lat,max_dN_lon,nblock_lat,nblock_lon,gcm_Plat2D_GOES,gcm_Plon2D_GOES);
%[off_y_Nd,off_x_Nd,max_corr_Nd,lat_track,lon_track]=track_motion_Dan(Nd_t0,Nd_t1,max_dN_lat,max_dN_lon,nblock_lat,nblock_lon,gcm_Plat2D_GOES,gcm_Plon2D_GOES);


%can plot using 
%m_vec(10,lon_track,lat_track,off_x_Vis,off_y_Vis,'k')

%off_x_LWP = off_x_LWP * 2.4; %convert to km (approx)
%off_y_LWP = off_y_LWP * 4.4; %convert to km (approx)
%off_x_Teff = off_x_Teff * 2.4; %convert to km (approx)
%off_y_Teff = off_y_Teff * 4.4; %convert to km (approx)
%off_x_Vis = off_x_Vis * 2.4; %convert to km (approx)
%off_y_Vis = off_y_Vis * 4.4; %convert to km (approx)
%off_x_Nd = off_x_Nd * 2.4; %convert to km (approx)
%off_y_Nd = off_y_Nd * 4.4; %convert to km (approx)

%fprintf(1,'LWP: off_x=%f off_y=%f max_corr=%f\n',off_x_LWP,off_y_LWP,max_corr_LWP);
%fprintf(1,'Teff: off_x=%f off_y=%f max_corr=%f\n',off_x_Teff,off_y_Teff,max_corr_Teff);
%fprintf(1,'Vis: off_x=%f off_y=%f max_corr=%f\n',off_x_Vis,off_y_Vis,max_corr_Vis);
%fprintf(1,'Nd: off_x=%f off_y=%f max_corr=%f\n',off_x_Nd,off_y_Nd,max_corr_Nd);


off_x_LWP_save{i_it} = off_x_LWP;
off_y_LWP_save{i_it} = off_y_LWP;

time_offxy_save(i_it) = times_GOES_save(it0);



end

% % for 74 will use the one before as that is closest
% %off_x_LWP_save{74} = off_x_LWP_save{73};
% %off_y_LWP_save{74} = off_y_LWP_save{73};
% 
% off_x_LWP_save{74} = off_x_LWP_save{74}/2; %since is an hour gap
% off_y_LWP_save{74} = off_y_LWP_save{74}/2; %since is an hour gap
% 
% 
% % for 76 will use the one before as that is closest
% off_x_LWP_save{76} = off_x_LWP_save{75};
% off_y_LWP_save{76} = off_y_LWP_save{75};
% 
% % for 77 will use the one after as that is closest
% off_x_LWP_save{77} = off_x_LWP_save{78};
% off_y_LWP_save{77} = off_y_LWP_save{78};
% 
% off_x_LWP_save{79} = off_x_LWP_save{79}/2; %since is an hour gap
% off_y_LWP_save{79} = off_y_LWP_save{79}/2; %since is an hour gap

% -- save the winds/offsets in a .mat file

save(save_file_track,'off_x_LWP_save','off_y_LWP_save','lat_track','lon_track','time_offxy_save','it_inds_offxy',...
    ,'gcm_Plat2D_GOES_nonan','gcm_Plat2D_GOES_nonan','-V7.3');

%% --- pick some feature points to mark on plots to test tracking

for it0=[67:82]  %

    ioverride_goes_multi=1;
    data_source = 'goes_multi_array'; imulti_goes=it0;
    POC_26Oct2008_GOES_general_maps_20141125T032943
    
    m_vec(8,lon_track,lat_track,off_x_LWP_save{it0},-off_y_LWP_save{it0},'k');
    
    saveas_ps_fig_emf(gcf,[savename],'',0,1);
    
    
end



run_example='yes';
run_example='no';

switch run_example 
    case 'yes'
        
%Using Blockmatcher (need to add vision toolbox and the /usr/local/MATLAB/R2013a/toolbox/shared/ to the path)
% Looks like the values returned in motion are imaginary - the real parts
% are x velocities and the imag. part y velocities.        

img1 = im2double(rgb2gray(imread('onion.png')));
htran = vision.GeometricTranslator('Offset', [5 5], 'OutputSize', 'Same as input image');
hbm = vision.BlockMatcher('ReferenceFrameSource', 'Input port', 'BlockSize', [35 35]);
hbm.OutputValue = 'Horizontal and vertical components in complex form';
halphablend = vision.AlphaBlender;

% Offset the first image by [5 5] pixels to create second image
img2 = step(htran, img1);

% Compute motion for the two images
motion = step(hbm, img1, img2);

% Blend two images
img12 = step(halphablend, img2, img1);

% Use quiver plot to show the direction of motion on the images
[X Y] = meshgrid(1:35:size(img1, 2), 1:35:size(img1, 1));
imshow(img12); hold on;
quiver(X(:), Y(:), real(motion(:)), imag(motion(:)), 0); hold off;

end