ioverride_plotglobal_loc=0;

% ----------------------------------------------------------------------------------
irestrict_domain=0 %whether to restrict the domain or not
% ----------------------------------------------------------------------------------

gcm_time_of_day_select=0; %is set in the plot type specification below. Set to 2 for time
%screening based on local time (but will need to change the way that
%time_inds_average is used

inew_cticks=0; %

 %type of projection
 iset_min_clim=0;  %flags for setting the colorscale limits or not
 iset_max_clim=0;

 proj_type='polar';
 proj_type='global oval'; %select this for anything but polar
 %        stereographic
 
 iplot_mgrid_lines=0; %whether to plot the grid lines for maps using m_grid.

 %        data_select='latest_modis_data';
 data_select='specific_modis_data';
 ifilter_ndays=0 %flag for whether to cut out points for which there aren't many days
 icontour=0;
 cont_col_str='c'; %Colour for the contour lines
 %N.B. - 'w' colour doesn't show up when saving for some reason?
 %Use 'c' (cyan) instead
 
 MODIS_varname2_plot='plot_global_maps plot';
 title_info='';
 units_str_plot='';

 noplot=0;
 
 savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';
 
 inew_figure=1;
 supress_colorbar=0; %Setting to one stops the colorbar from appearing
 i_increase_font_size_map_figures_OFF = 0;  % setting to one stops this from running

 plot_num_datapoints=0;


 
 
 
 screen_type ='none';
 
 
 

 
  thresh_ndays=0; %threshold no. days
%        thresh_ndays=5; %threshold no. days
%        thresh_ndays=1; %threshold no. days  
%        thresh_ndays=100; %threshold no. days          

        thresh_SZA=[65 90];
        thresh_SZA=[0 65];         
%        thresh_SZA=[0 90];        
        thresh_CF=0.8;
        thresh_CF=[0.8 1.00001];
%        thresh_CF=[0.1 1.00001];        
%        thresh_CF=[0 1.00001];        
%        thresh_CF=[0.0 0.4];        
%        thresh_CF=[-0.01 1.00001];


%        thresh_CF=[-0.01 10e-2];   %trying to look at AMSRE clear-sky bias
        
        thresh_NP=50;
%        thresh_NP=0;
        
        thresh_sensZA=45;
        thresh_sensZA=50;
        thresh_sensZA=40;
        thresh_sensZA=[0 41.4]; 
        thresh_sensZA=[0 90];         
        thresh_CTT = 273;
        thresh_CTT = [273-100 273+100];        
%        thresh_CTT = [273 273+100];    
%        thresh_CTT = [268 273+100];   
%        thresh_CTT = [273-10 273+100];  
        
        thresh_reff = [0 30];
        thresh_reff = [0 12];        
        
        thresh_CTP = [730 1000]; %Cloud top pressure (hPa)
        
        thresh_relAZ = [0 180];

        %    thresh_sensZA=80;
        thresh_maxSZA=[-1 81.4];
        thresh_minSZA=45;        
        thresh_stdSZA=[-1 0.5e9];
        thresh_stdSZA=[-1 0.5];
        
        thresh_dSZA=[-1 1e9];        
       thresh_dSZA=[-1 1];
        
        thresh_Nd_per_error = 100;
        thresh_Reff_per_error = 50;
        thresh_Reff_abs_error = 4;
        thresh_Reff = 30;
        
        gcm_CTT_thresh = [273.15 400];
        
        thresh_CTH = [-0.01 1e9]; %km
%        thresh_CTH = [-0.01 3.2]; %km  - the approx height of 680 hPa (low cloud pressure threshold for ISCCP, CALIPSO, etc)
%        thresh_CTH = [3.2 20];
%        thresh_CTH = [6.5 20];        

        thresh_zeroCF = 0.05;
        
        thresh_stdW = [0 1e9];
%        thresh_stdW = [0 5];
        
        thresh_sigCTT = [0 1e9];
%        thresh_sigCTT = [0 1];  


      minfrac_CF = 0.9; %minimum fraction of the sampled points that had successful cloudy/clear/phase
      %determination (i.e. Npix/Nptot_mockL3 =
      %Cloud_Fraction_Liquid_Pixel_Counts./Cloud_Fraction_Liquid./Total_pixels
      % - restriction (2) as presented in the paper

      minfrac_NpNd = 0.9;        
       %Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3./Cloud_Fraction_Liquid_Pixel_Counts.timeseries3
       %Fraction of points that remain after all previous filtering for
       %which we have an Nd retrieval. Restriction (4) in the SZA paper.
            
      thresh_NP_Nd = 50; %min no. of pixels required for an Nd, re, tau, etc measurement to count
      %(Uses Cloud_Fraction_Liquid_Pixel_Counts2). For this screening
      %usually
      %thresh_NP is the number of pixels that the swath must have covered -
      %i.e. the total number of pixels available. Only a portion of those 
      
        
%% Latitudes
        thresh_LAT = [70 74];
        %                      thresh_LAT = [72 73];


        %       thresh_LON = [-160 -135];
        %       thresh_LON = [-160 -140];
        %                       thresh_LON = [-152.5 -150];
        thresh_LAT = [67 75]; thresh_LON = [-180 -135]; %MPACE domain
        %         thresh_LAT = [70.33 71.33]; thresh_LON = [-154 -150]; %MPACE domain - zoomed into flight track
        % %        thresh_LAT = [70.33 71.33]; thresh_LON = [-158 -148]; %MPACE domain - zoomed into flight track even closer
        % %        thresh_LAT = [69 72]; thresh_LON = [-160 -146]; %MPACE
        %         thresh_LAT = [70 71.5]; thresh_LON = [-160 -146]; %MPACE
        %         thresh_LAT = [70 71.5]; thresh_LON = [-151 -146]; %MPACE - profiles to the east
        %         thresh_LAT = [70 71.5]; thresh_LON = [-153 -148]; %MPACE - profiles
        %         thresh_LAT = [70 71.5]; thresh_LON = [-153.5 -151.5]; %MPACE - profiles
        thresh_LAT = [60 64]; thresh_LON = [25 29]; %Sami's station (Puijo)
        thresh_LAT = [62.5 63.5]; thresh_LON = [25 29]; %Sami's station   (@ 62.909 N, 27.656 E)
%        thresh_LAT = [50 70]; thresh_LON = [10 40]; %Sami's station   
        
%                thresh_LAT = [-40 -10];  thresh_LON = [-100 -60];
                %VOCALS CAPT:- One used for CPT stuff in general 10th Jan
                %2013
                thresh_LAT = [-40 10];  thresh_LON = [-140 -50]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
                
                %To match Fig. 13 of Toniazzo 2011 paper 
%                 thresh_LAT = [-35 -10];  thresh_LON = [-90 -65]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
                 
               %GOES domain to roughly match the area of the VOCALS GOES
               %images
%                thresh_LAT = [-33 -7];  thresh_LON = [-103 -67]; lon_ticks=[-10:5:-65]; lat_ticks=[-30:5:-5];%VOCALS CAPT
 
                %new CPT domain extended up to California
%                thresh_LAT = [-45 45];  thresh_LON = [-150 -55]; 
                
                
%                thresh_LAT = [-60 60];  thresh_LON = [-180 -50]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT                %
                
%                thresh_LAT = [-50 50];  thresh_LON = [-160 -30]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
                
 %               thresh_LAT = [0 65]; thresh_LON = [-160 -60]; %Eastern US,
        %        UK western Europe
%                thresh_LAT = [0 85]; thresh_LON = [-180 180]; %N. hemisphere - whole globe    
        %        thresh_LAT = [-85 0]; thresh_LON = [-180 180]; %S. hemisphere - whole globe    
                
%                thresh_LAT = [0 65]; thresh_LON = [-180 180]; %N. hemisphere - whole globe
%                 thresh_LAT = [-65 0]; thresh_LON = [-180 180]; %S.hemisphere - whole globe
%                   thresh_LAT = [36 52]; thresh_LON = [-30 0]; %UK, Atlantic off SW coast of UK
%            thresh_LAT = [44 51]; thresh_LON = [-30 -6]; %UK, Atlantic off SW coast of UK

%                    thresh_LAT = [20 60]; thresh_LON = [-60 30]; %UK, Atlantic off SW coast of UK
                   
        %    thresh_LAT = [44 52]; thresh_LON = [-15 0]; %UK, Atlantic off
        %    SW coast of UK
        %    thresh_LAT = [36 65]; thresh_LON = [-15 60]; %mainland Europe

%                 thresh_LAT = [-10 60];  thresh_LON = [-36 60]; %Europe
%                 thresh_LAT = [30 80];  thresh_LON = [-90 80]; %UK and Iceland

%                 thresh_LAT = [-80 -40]; thresh_LON = [-200 -30];
%                 %Antarctica
%                 thresh_LAT = [-75 -60]; thresh_LON = [-80 -50];
%                 %Antarctica   
%                 thresh_LAT = [-72 -63]; thresh_LON = [-72 -55];
%                 %Antarctica    
%                 thresh_LAT = [-70 -65]; thresh_LON = [-72 -55]; %Used for paper. Antarctica Peninsula and Larsen C
%also set  detailed_BAS_coast=1 (run convert_Antarctic_esri_file_lat_lon_mmap to convert a shapefile to m-map lat/lon co-ords)                 
%                 thresh_LAT = [-69 -65]; thresh_LON = [-69 -64]; %Antarctica Peninsula - close up on Rothera 

%                 thresh_LAT = [-75 -40]; thresh_LON = [-80 -50]; %Antarctica   
%                 thresh_LAT = [-75 -40]; thresh_LON = [-90 -20];
%                 %Antarctica                    

%                 thresh_LAT = [-90 -50]; thresh_LON = [-180 180];  %Southern pole


%requested region was 70-80N, -20 to 60E
%thresh_LAT = [70 80]; thresh_LON = [-20 60]; %Arctic summer requested
%region
%thresh_LAT = [65 85]; thresh_LON = [-20 60]; %Arctic summer requested region

%thresh_LAT = [-23 -16]; thresh_LON = [-90 -80]; %

%thresh_LAT = [-60 -10]; thresh_LON = [110 155]; %Australia
%thresh_LAT = [-69 10]; thresh_LON = [90 180]; %

%thresh_LAT = [0 60]; thresh_LON = [90 180]; %China

%thresh_LAT = [-90 90]; thresh_LON = [-180 180]; %global

%thresh_LAT = [-65 -35]; thresh_LON = [-180 180]; %Southern Ocean
%thresh_LAT = [-60 -45]; thresh_LON = [50 100]; %Southern Ocean box
%thresh_LAT = [-65 -30]; thresh_LON = [0 240]; %Southern Ocean box

%Cape Grimm is 41S, 144.5E

% Hawaii - Amy's domain
%thresh_LAT = [12 25]; thresh_LON = [-170 -150]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Downwind of Hawaii
%thresh_LAT = [0 25]; thresh_LON = [-180 -150]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Downwind of Hawaii
%wider picture
%thresh_LAT = [-25 25]; thresh_LON = [-180 180]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Downwind of Hawaii

if length(thresh_CF)==1
    thresh_CF(2)=1;
end
if length(thresh_CTP)==1
    thresh_CTP(2)=1100;
end


 