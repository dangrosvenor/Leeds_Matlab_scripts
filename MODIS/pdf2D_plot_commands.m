% modisL3_screening_timeseries3   -- script to do the screening
% time_inds_modisL3_timeseries3   -- script to do the time selection
% CTRL + up and down arrows moves between blocks

% N.B. for a plto of mean values insteqad of a histogram set ndims_hist=3
% and set teh variable for the mean in Z

try %Catching errors so that flags are reset even if the program
    %aborts due to an error. Commands after the catch statement are
    %executed
    %if an error now occurs. The flags are then reset and the error is
    %"rethrown"
    
% *** IMPORTANT - make sure that (ilat,ilon,itime) is used for the arrays - even if all of the array 
% is being used since the itime indices are not in order, which means the
% order will be different from that of the screening arrays if itime is not
% used. ***

clear X Y Z

% -----------------------------------------------------------------------
disp('***** DON''T forget to set datatype *******');

if ~exist('ioverride_pdf') | ioverride_pdf==0

        if max(strcmp(gcm_str,{'MODIS','AMSRE','POLDER','AMSRE_time3','switchable'}))==1
            datatype='timeseries3'; gcm_time_of_day_select=0; %e.g. for MODIS
        elseif max(strcmp(gcm_str,{'GOES','UM','AMSRE2'}))==1
            datatype='gcm_data';    gcm_time_of_day_select=2;
        else
            datatype='timeseries3';    gcm_time_of_day_select=0; %e.g. for MODIS
            %datatype='mock L3'; %lat-lon 1 degree cells
            datatype='gcm_data';    gcm_time_of_day_select=2;
            %datatype = 'other'; %works for POLDER
            datatype = 'makeshift';
        end
        
           daynight='Daytime';
%   daynight='Nighttime'; %for LWP normalised by model CF
%    daynight='ALL times of day';

            ireduce_res_lwp=0

        
end

if ~exist('extra_title_info')
    extra_title_info='';
end
extra_title_info=''; %this was getting too long due to self-concatenation


if ~exist('ireduce_res_lwp')
    ireduce_res_lwp=0;
end


   


 % -----------------------------------------------------------------------


if ~exist('time_series_type')
    time_series_type = 'daily';
end

%select the timeseries type required
timeseries_case=''; %for timeseries
timeseries_case='3'; %for timeseries3

%LAT=eval(['Cloud_Optical_Thickness_Liquid_Mean.timeseries' timeseries_case '_LAT']);
%LON=eval(['Cloud_Optical_Thickness_Liquid_Mean.timeseries' timeseries_case '_LON']);

 %set some default values (leaving these unchanged will result in the bins
 %going from the min to the max value in nXpdf bins
    minXbins=-9e99
    maxXbins=-9e99;

    minYbins=-9e99;
    maxYbins=-9e99;
    
    minZbins=-9e99;
    maxZbins=-9e99;
    
    ikeep_X_above_zero=0;
    
    post_plotTime_commands={''};
   
    
%    iplot_mean_XY = 'x';
%    iplot_mean_XY = 'y';    

if ~exist('ioverride_pdf') | ioverride_pdf==0

    iplot_mean_XY = ''; %set to 'x' or 'y' for x or y means
    
    ndims_hist=2;
%    ndims_hist=3;

    %select the number of PDF bins
    nXpdf=10000;
%    nXpdf = 12;
   nXpdf = 150;
%   nXpdf = 1e3;

    nYpdf=10;
%    nYpdf=50; 
%    nYpdf=50;
%    nYpdf=25;
%    nYpdf=12;
%    nYpdf=8;    
    
    nZpdf=100;



   

   
    
    ichoose_Xbins=0; %flags to allow the direct specification of Xbins, Ybins, Zbins
    ichoose_Ybins=0;
    ichoose_Zbins=0;
    
    ipost_plotTime_commands = 0;
    
    
else
    %    clear ioverride_pdf
end

%% Choose lat and lon
    
    if ~exist('ioverride_location_selection') | ioverride_location_selection==0

            %select the LAT and LON values to make up the
            %points included - these will be the EDGES of the box required
            LAT_val = [-0.5000];
            %LAT_val = [-60.5000];
            %LAT_val = [-22.5000];
            LAT_val = [-60.5000];
            LAT_val = [-22.5:-18.5];
            LAT_val = [35.5:37.5];            
            %LAT_val = [-59.5000];
            %LAT_val = [-89.5:1:89.5];
            %LAT_val = [-30.5:1:30.5];
            %LAT_val = [-20.5:1:20.5];
            %LAT_val = [-10.5:1:10.5];
            %LAT_val = [-30.5:1:-10.5];
            %LAT_val = [-12.5000];

            %LON_val = [-80:-70];
            %LON_val = 0;


            LAT_val = [0.5:89.5]; LON_val = [-179.5:179.5]; %global
            % LAT_val = [36:65]; LON_val = [-15:60]; %mainland Europe
            LAT_val = [44:51]; LON_val = [-30:-6]; %SW UK/Atlantic - ocean only
            
            %LON_val = [-76.5 -72]; LON_val=[-81 -76.5]; LON_val = [-85.5-81];
            %LON_val = [-90 -85.5];

%            LAT_val = [-22 -18]; LON_val = [-76.5 -72];

%            LAT_val = [-72:-64]; LON_val = [-120:-90]; %Antarctica, random region
%            LAT_val = [-71:-70]; LON_val = [-120:-90]; %Antarctica, random region
%            LAT_val = [-80:-40]; LON_val = [-120:-90]; %Antarctica, whole latitude range
            %             LAT_val = [-65:-64]; LON_val = [-120:-90]; %Antarctica, random region
            %             LAT_val = [-77:-76]; LON_val = [-120:-90]; %Antarctica, random region

%            LAT_val = [72:75]; LON_val = [0:50]; %Arctic
            %             LAT_val = [40:90]; LON_val = [0:50]; %Arctic whole region
            


%            LAT_val = [72 75]; LON_val = [-3:48]; %Arctic summer box region
%            LAT_val = [72 75]; LON_val = [-3:22]; %Arctic summer box western half
%            LAT_val = [72 75]; LON_val = [23:48]; %Arctic summer box eastern half            
%            LAT_val = [-30.5 -19.5]; LON_val = [-73.25 -72.25]; %VOCALS one lon
%            LAT_val = [-30.5 -19.5]; LON_val = [-77.25 -76.25];
%            LAT_val = [-29.5 -10.5]; LON_val = [-103.25 -72.25]; %one used for AGU talk 2012
%            LAT_val = [-25.5 -14.5]; LON_val = [-140.25 -70.25]; %Terry's VOCALS region
%            LAT_val = [-30.5 10.5]; LON_val = [-140.25 -70.25]; %Wider VOCALS region including low lats
            LAT_val = [-40.5 10.5]; LON_val = [-140 -50]; %VOCALS CAPT (whole map for CPT paper plots)
            LAT_val = [-40.5 -30.5]; LON_val = [-140 -100]; %Smaller SW region of VOCALS where have the issue with high clouds, etc.
%            LAT_val = [-25.5 -15.5]; LON_val = [-80 -70]; %Region near coast where have Sc
            
%            LAT_val = [-29.5 -10.5]; LON_val = [-103.25 -81]; %like above, but avoiding all land            
%            LAT_val = [-29.5 -20.5]; LON_val = [-103.25 -81]; %       
%            LAT_val = [-29.5 -18.5]; LON_val = [-103.25 -72.25];            
%            LAT_val = [-20.5 10.5]; LON_val = [-130 -100]; %Approx region to the west where MODIS shows negative LWP bias
%            LAT_val = [-25.5 -14.5]; LON_val = [-130 -77]; %Matching the POLDER slide         
%            LAT_val = [-25.5 -14.5]; LON_val = [-120 -77]; %Matching the POLDER slide                     
            LAT_val = [-25.5 0.5]; LON_val = [-120 -77]; %Extending the POLDER result to larger lat range for more samples            

%            LAT_val = [-25.5 0.5]; LON_val = [-95 -77]; %Easternmost region only
%            LAT_val = [-25.5 0.5]; LON_val = [-120 -105]; %Westernmost region only  
            
%            LAT_val = [-24.5 -15.44]; LON_val = [-86.93 -77.08]; %GOES regin for UM comparison xkqk 26thOct POC
            
%           LAT_val = [25 34]; LON_val = [122:129]; %China sea

%            LAT_val = [-60 -57]; LON_val = [50 100]; %restritcted latitude Southern Ocean
%            LAT_val = [-60 -45]; LON_val = [50 100]; %Southern Ocean, 50-100
%            LAT_val = [-60 -45]; LON_val = [-180 180]; %Southern Ocean, all lons   
%            LAT_val = [-60 -40]; LON_val = [-180 180]; %Southern Ocean, all lons   
%            LAT_val = [-47 -35]; LON_val = [-180 180]; %Southern Ocean, all lons              
% %           LAT_val = [-47 -31]; LON_val = [-180 180]; %Southern Ocean, DJF lat region for points 
%            %affected by the <65 SZA restriction
%            %            LAT_val = [-60 -57]; LON_val = [50 100]; %Southern Ocean, restricted latitude       
% %            LAT_val = [-45 -44]; LON_val = [50 100]; %Southern Ocean, restricted latitude  
% %            LAT_val = [-60 -59]; LON_val = [50 75]; %Southern Ocean, restricted latitude              
% %            LAT_val = [-60 -59]; LON_val = [75 100]; %Southern Ocean, restricted latitude                  
% %            LAT_val = [-60 -59]; LON_val = [50 100]; %Southern Ocean, restricted latitude                              
% %            LAT_val = [-47 -44]; LON_val = [50 100]; %Southern Ocean, restricted latitude             
% %            LAT_val = [-43 -40]; LON_val = [143 144.5]; %Tasmania, Boers, OJRMS, 1998 - Southern Ocean - too small            
% %            LAT_val = [-43 -40]; LON_val = [110 144]; %Tasmania, Boers, OJRMS, 1998 - Southern Ocean -extended
% %LAT_val = [-60 -40]; LON_val = [40 144];
% %           LAT_val = [62 63]; LON_val = [27 28]; %Puijo (Sami) : 62.91N, 27.66 E
% %            LAT_val = [62 64]; LON_val = [26 29]; %Puijo (Sami) : 62.91N, 27.66 E 
% 
% %           LAT_val = [30 60]; LON_val = [-180 180]; %
% 
%            LAT_val = [-90 90]; LON_val = [-180 180]; %


%            LAT_val = [-1e9 1e9]; LON_val = [-1e9 1e9];
           
            %*** these will be the EDGES of the box required  ***

    else
        %        clear ioverride_location_selection
    end


    if LAT_val(1)>LAT_val(end)
        fprintf(1,'\n***  LAT_val needs to run from lowest to highest!  ***\n');
        return
    end

    if LON_val(1)>LON_val(end)
        fprintf(1,'\n***  LON_val needs to run from lowest to highest!  ***\n');
        return
    end








    
    switch datatype
        case 'makeshift'
        case {'gcm_data'}
            if exist(['gcm_lwp_' gcm_str])

                %data for screening
                screen_lwp = eval(['1e3*gcm_lwp_' gcm_str]); %g m^{-2}
                screen_lwp = permute(screen_lwp,[2 3 1]);

                if exist(['gcm_TLWP_' gcm_str])
                    screen_tlwp = eval(['1e3*gcm_TLWP_' gcm_str]); %g m^{-2}
                    screen_tlwp = permute(screen_tlwp,[2 3 1]);
                else
                    screen_tlwp = eval(['1e8*ones(size(screen_lwp));']); %g m^{-2}
                end

                screen_precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s
                screen_precip = screen_precip *3600*24;  %mm/day
                screen_precip = permute(screen_precip,[2 3 1]);

                screen_cf = eval(['cf_isccp_low_' gcm_str]);  %Clouf fraction (range 0 to 1)
                screen_cf = permute(screen_cf,[2 3 1]);

            end

            %other data
            Plat2D = eval(['gcm_Plat2D_' gcm_str]);
            Plon2D = eval(['gcm_Plon2D_' gcm_str]);

            Plat2D_edges = eval(['gcm_Plat2D_edges_' gcm_str]);
            Plon2D_edges = eval(['gcm_Plon2D_edges_' gcm_str]);

            dlat=meanNoNan(meanNoNan(diff(Plat2D),1),1);
            dlon=meanNoNan(meanNoNan(diff(Plon2D,1,2),1),1);  %for writing the actual boundaries (assuming equal spacing)
            %


            stime = eval(['size(gcm_time_matlab_' gcm_str ')']);
            Plat3D = repmat(Plat2D,[1 1 stime]);
            Plon3D = repmat(Plon2D,[1 1 stime]);
            Plat3D_edges = repmat(Plat2D_edges,[1 1 stime]);
            Plon3D_edges = repmat(Plon2D_edges,[1 1 stime]);

            gcm_time_UTC = eval(['gcm_time_UTC_' gcm_str]);
            daynum_timeseries3 = eval(['daynum_timeseries3_' gcm_str]);
            modisyear_timeseries3 = eval(['modisyear_timeseries3_' gcm_str]);


            switch gcm_str
                case 'CALIPSO_monthly'
                    Plat3D = permute(Plat3D,[3 1 2]);
                    Plon3D = permute(Plon3D,[3 1 2]);
                    Plat3D_edges = permute(Plat3D_edges,[3 1 2]);
                    Plon3D_edges = permute(Plon3D_edges,[3 1 2]);
            end

            iregion_lin = find(Plat3D>=LAT_val(1) & Plat3D<LAT_val(end) & Plon3D>=LON_val(1) & Plon3D<LON_val(end));
            iregion_lin_edges = find(Plat3D_edges>=LAT_val(1) & Plat3D_edges<LAT_val(end) & Plon3D_edges>=LON_val(1) & Plon3D_edges<LON_val(end));
            iregion_lin2D = find(Plat2D>=LAT_val(1) & Plat2D<LAT_val(end) & Plon2D>=LON_val(1) & Plon2D<LON_val(end));
            iregion_lin2D_edges = find(Plat2D_edges>=LAT_val(1) & Plat2D_edges<LAT_val(end) & Plon2D_edges>=LON_val(1) & Plon2D_edges<LON_val(end));

            if length(iregion_lin)==0
                fprintf(1,'\n*** WARNING - region is outside of Plat3D and Plon3D in pdf2d_plot_commands.m *** \n\n')
            end

            switch gcm_str
                case {'CALIPSO_monthly'}
                    time_series_type = 'CALIPSO';

                case {'GOES','UM','AMSRE2'}
                    %do nothing

                otherwise
                    s3D = size(Plat3D);
                    iocean = eval(['find(gcm_landmask_' gcm_str '==0);']);
                    [i_ocean,j_ocean,k_ocean] = ind2sub(s3D,iocean);
                    k2_ocean = repmat([1:s3D(3)],[length(i_ocean) 1]);
                    i2_ocean = repmat(i_ocean,[1 s3D(3)]);
                    j2_ocean = repmat(j_ocean,[1 s3D(3)]);
                    iregion_ocean_lin = sub2ind(s3D,i2_ocean(:),j2_ocean(:),k2_ocean(:));
                    % LAT=Plat2D;
                    % LON=Plon3D;

            end
            %%
        otherwise

            switch gcm_str
                case 'CLOUDSAT_PRECIP'
                    LAT = Plat2D_matt_centres(:,1);
                    LON = Plon2D_matt_centres(1,:);
            end

            if ~exist('LAT')
                LAT = eval(['LAT_' gcm_str]);
                LON = eval(['LON_' gcm_str]);
            end

            %LAT and LON are the mid-points

            %            ilat=findheight_nearest(LAT,LAT_val);
            %            ilon=findheight_nearest(LON,LON_val);
            dlat=prctile(diff(LAT),40);
            dlon=prctile(diff(LON),50);  %for writing the actual boundaries (assuming equal spacing)

            ilat = find(LAT>=LAT_val(1) & LAT<LAT_val(end));
            ilon = find(LON>=LON_val(1) & LON<LON_val(end));

    end

switch datatype
    case 'timeseries3'
        switch gcm_str
            case 'POLDER'
                daynum_timeseries3 = daynum_timeseries3_POLDER;
                modisyear_timeseries3 = modisyear_timeseries3_POLDER;   
            case 'AMSRE'
                daynum_timeseries3 = daynum_timeseries3_AMSRE;
%                modisyear_timeseries3 = modisyear_timeseries3_AMSRE;  
            case 'AMSRE_time3'
                daynum_timeseries3 = daynum_timeseries3_time3;
            otherwise
               daynum_timeseries3 = eval(['daynum_timeseries3_' gcm_str ';']);
               modisyear_timeseries3 = eval(['modisyear_timeseries3_' gcm_str ';']);
        end
end


nstart = 1;
%nstart = 32; %1st Feb
%nstart = 60; %1st March = day 60 (except in leap years when =61)
%nstart = 91; %1st March = day 60 (except in leap years when =61)
%nstart = 152; %1st June
%nstart = 213; %1st Aug
nstart = 244; %1st Sep
nstart = 274; %1st Oct
%nstart = 305; %1st Nov
%nstart = 334; %=335 1st Dec, but use 334 otherwise not enough days for 32 days
nstart = 355; %21st Dec 2005
%nstart = 1+96*3;
%nstart = 365-96+1;





switch datatype
    case 'timeseries3'
%        size_tim3_ALL = size(Cloud_Fraction_Liquid.timeseries3);
%Don't think this is needed anymore?
end


naim = 3/12*365; %=91.5 - so the nearest divisible by 16 is 96
naim = 32;
naim = 8;
naim = 1;

itime_orig=[nstart:nstart-1+16*round(naim/16)];
itime_orig=[nstart:nstart-1+naim];

%itime_orig=9:16;

%itime=[1:32*3];
%itime=[366:size_tim3_ALL(3)];
%itime=[1:size_tim3_ALL(3)];
%itime=1:365;

year_ends=[365 365];
%year_ends=[365 365 365 365 365 365]; %2005-2007 aqua and terra
%N.B. - not going to use the last one of these, but can put it in anyway

year_ends_cum=cumsum(year_ends);



itime=itime_orig;
for it=1:length(year_ends)-1
    itime=[itime itime_orig+year_ends_cum(it)]; %to include the aqua data as well
end

%% time screening
switch datatype
    case 'makeshift'
    otherwise
% ************************ time - external script ****************
time_inds_modisL3_timeseries3  %this is re-run later for some variables
% **** ***********************************************************
itime = time_inds_average;


date_range_str = [datestr(itime_orig(1)+datenum('01-Jan-2005')-1,'dd-mmm') ' to ' datestr(itime_orig(end)+datenum('01-Jan-2005')-1,'dd-mmm') ];

end


switch datatype
    case 'timeseries3'
        switch gcm_str
            case 'POLDER'
                size_time_dim = length(daynum_timeseries3_POLDER);
                size_tim3 = [length(ilat) length(ilon) length(itime)];
                %POLDER X values are permuted later since they are [time
                %lat lon]
            case {'AMSRE'}
                size_time_dim = size(lwp_amsre,3);
                size_tim3 = size(lwp_amsre(ilat,ilon,itime,1));
            case 'AMSRE_time3'
                size_time_dim = size(lwp_amsre_time3,3);
                size_tim3 = size(lwp_amsre_time3(ilat,ilon,itime,1));
            case 'GOES_time3'
                size_time_dim = size(lwp_goes_time3,3);
                size_tim3 = size(lwp_goes_time3(ilat,ilon,itime,1));
            case 'GOES'
                size_time_dim = size(goes_LWP,3);
                size_tim3 = size(goes_LWP(ilat,ilon,itime,1));
            case 'UM'
                size_time_dim = length(daynum_timeseries3_UM);
                size_tim3 = [ilat ilon itime];  
            case 'AMSRE2'
                size_time_dim = length(daynum_timeseries3_AMSRE2);
                size_tim3 = [ilat ilon itime];
            case 'switchable'
                size_time_dim = eval(['length(daynum_timeseries3_' gcm_str ')';]);
                size_tim3 = [length(ilat) length(ilon) length(itime)];
            case 'CLOUDSAT_PRECIP'    
                size_time_dim = eval(['length(daynum_timeseries3_' gcm_str ')';]);
                size_tim3 = [length(ilat) length(ilon) length(itime)];
            otherwise
                size_time_dim = size(Cloud_Fraction_Liquid.timeseries3,3);
                size_tim3 = size(Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime));
        end
end

if ~exist('iswap_xy')
    iswap_xy=0; %flag to swap the defined x and y
             %  -- is also set below --
end

%% Choose x variables

if ~exist('ioverride_pdf_varchoose') | ioverride_pdf_varchoose==0

iswap_xy=0; %flag to swap the defined x and y

x_axis_vals = 'Dummy data';
x_axis_vals = 'General GCM-style x-axis';
%x_axis_vals = 'Nd dist skewness lat-lon cells';
% x_axis_vals = 'dNd/dSZA overall';
%x_axis_vals = 'Nd from grid vals timeseries3';
%x_axis_vals = 'Nd_{1.6} from grid vals timeseries3';
%x_axis_vals = 'Nd_{3.7} from grid vals timeseries3';
%x_axis_vals = 'Nd from individual pixels timeseries3';
%x_axis_vals = 'Nd from individual pixels 1.6\mum timeseries3';
%x_axis_vals = 'Nd from individual pixels 3.7\mum timeseries3';
%x_axis_vals = 'Nd uncertainty from grid vals timeseries3';
%x_axis_vals = 'Nd % uncertainty from grid vals timeseries3';
%x_axis_vals = 'Re % uncertainty from grid vals timeseries3';
%x_axis_vals = 'Optical Depth % uncertainty from grid vals timeseries3';
%x_axis_vals = 'H from grid vals timeseries3';
%x_axis_vals = 'LWP grid-cell mean from grid vals timeseries3'; %MODIS LWP - can select MOD06 ot MOD35 scaling
%x_axis_vals = 'MODIS LWP in-cloud timeseries3';
%x_axis_vals = 'MODIS minus AMSRE LWP bias';
%x_axis_vals = 'LWP 1.6um MODIS grid-cell mean from mockL3'; %MODIS LWP - can select MOD06 or MOD35 scaling
%x_axis_vals = 'LWP 2.1um MODIS grid-cell mean from mockL3'; %MODIS LWP - can select MOD06 or MOD35 scaling
%x_axis_vals = 'LWP 3.7um MODIS grid-cell mean from mockL3'; %MODIS LWP - can select MOD06 or MOD35 scaling
%x_axis_vals = 'MODIS 1.6um minus AMSRE from mockL3'; %MODIS LWP - can select MOD06 or MOD35 scaling
%x_axis_vals = 'MODIS 2.1um minus AMSRE from mockL3'; %MODIS LWP - can select MOD06 or MOD35 scaling
%x_axis_vals = 'MODIS 3.7um minus AMSRE from mockL3'; %MODIS LWP - can select MOD06 or MOD35 scaling
%x_axis_vals = 'MockL3 reff 1.6um timeseries3';
%x_axis_vals = 'MockL3 reff 2.1um timeseries3';
%x_axis_vals = 'MockL3 reff 3.7um timeseries3';
%x_axis_vals = 'POLDER minus AMSRE LWP';
%x_axis_vals = 'POLDER reff daymean';

%x_axis_vals = 'CTP std dev timeseries3, x-axis';
%x_axis_vals = 'Nd from histogram vals timeseries';
%x_axis_vals = 'Nd from histogram vals timeseries3';
%x_axis_vals = 'Tau'; %timeseries3
%x_axis_vals = 'Tau COSP GCM'; %timeseries3
%x_axis_vals = 'LWP COSP GCM'; %timeseries3
%x_axis_vals = 'Re COSP GCM'; %timeseries3
%x_axis_vals = 'Nd maxliq GCM'; %timeseries3
%x_axis_vals = 'R_{eff 2.1 \mum} (\mum)';  %pick this one if using L3 data (i.e. no 1.6 and 3.7)
%x_axis_vals = 'R_{eff 1.6 \mum} (\mum)';
%x_axis_vals = 'R_{eff 3.7 \mum} (\mum)';
%x_axis_vals = 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%x_axis_vals = 'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%x_axis_vals = 'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7';
%x_axis_vals = 'Nd factor';
%x_axis_vals = 'Reff^{-5/2}'; %timeseries3
%x_axis_vals = 'Tau^{1/2}'; %timeseries3
%                x_axis_vals = 'Std. dev of Nd from histogram timeseries';
%x_axis_vals = 'Std. dev of Nd from histogram timeseries3';
%x_axis_vals = 'Mean SZA timeseries3';
%x_axis_vals = 'Relative Azimuth Angle timeseries3';
%x_axis_vals = 'Mean Sensor ZA timeseries3';
%x_axis_vals = 'Cloud fraction timeseries3';
%x_axis_vals = 'Std. dev of Reff from histogram timeseries';
%x_axis_vals = 'Std. dev of Reff';
%x_axis_vals = 'Minimum Reff timeseries3';
%x_axis_vals = 'Maximum Tau timeseries3';
%x_axis_vals = 'Sensor ZA timeseries3';
%x_axis_vals = 'Latitude';  %for MODIS or GCMs 
%x_axis_vals = 'Longitude';  %for MODIS or GCMs
%x_axis_vals='Mean CTT timeseries3';
%x_axis_vals = 'Cloud Top Temp standard deviation, liquid pixels';
%x_axis_vals = 'Nd dist skewness lat-lon cells';
%x_axis_vals = 'Homogeneity Parameter timeseries3 using W from mean tau and Re';
%x_axis_vals = 'Homogeneity Parameter timeseries3 using mean W for each pixel';
%x_axis_vals = 'Homogeneity Parameter timeseries3 using Cahalan log mean W (pixel level)';
%x_axis_vals = 'Homogeneity Parameter Cahalan Optical Depth (Seethala)';
%x_axis_vals = 'CF all'; %ice + liquid _+ undetermined cloud fraction
%x_axis_vals = 'Percentage Nd difference allSZA vs lowSZA';
%x_axis_vals = 'CDR Polder'; %uses the 4 orbits as separate datapoints 
   %- so can have more than one per day
%x_axis_vals = 'CDR Polder2'; %uses the mean over the orbits - only one datapoint per day
%x_axis_vals = 'CDR Polder Colocated with MODIS'; %
%x_axis_vals = 'CDR Polder select region'; %
%x_axis_vals = 'Dispersion Polder';
%x_axis_vals = 'Std Deviation Polder';
%x_axis_vals = 'In-cloud averaged AMSRE TLWP'; % (by CFliq)
%x_axis_vals = 'In-cloud averaged (using daytime MOD35) AMSRE TLWP';
x_axis_vals = 'Grid-box mean AMSRE TLWP';




%x_axis_vals = 'In-cloud LWP GCM norm by COSP-MODIS-CF'; %timeseries3
%x_axis_vals = 'In-cloud LWP GCM norm by Model CF'; %timeseries3
%x_axis_vals = 'In-cloud TLWP GCM norm by COSP-MODIS-CF'; %timeseries3
%x_axis_vals = 'qv700 GCM'; %timeseries3
%x_axis_vals = 'LTS GCM'; %timeseries3
%x_axis_vals = 'CDP Nd'; %timeseries3
%x_axis_vals = 'DMPS Nacc'; %timeseries3
%x_axis_vals = 'DMPS Nd estimate';
%x_axis_vals = 'DMPS Nd estimate CDP NaN match';
%x_axis_vals = 'CDP Nd all data';
%x_axis_vals = 'CF COSP-MODIS GCM';
%x_axis_vals = 'CF GCM';
%x_axis_vals = 'CF COSP-CALIPSO GCM';
%x_axis_vals = 'MOD35 Cloud Fraction from grid vals timeseries3'; %MOD35 CF
%                   - not necessarily liquid, though, or just low cloud
%x_axis_vals = 'Cloud Fraction from grid vals timeseries3'; %Liquid cloud fraction
%x_axis_vals = 'CALIPSO cloud fraction';
%x_axis_vals = 'LWP GCM'; %LWP divided by MODIS COSP CF
%x_axis_vals = 'LWP GCM grid-box mean'; %
%x_axis_vals = 'Pre-mphys LWP GCM grid-box mean'; %
%x_axis_vals = 'Post-mphys LWP GCM grid-box mean'; %
%x_axis_vals = 'Post minus pre-mphys LWP GCM grid-box mean'; %
%x_axis_vals = 'Fraction of LWP lost to mphys GCM grid-box mean';
%x_axis_vals = 'LWP removal rate due to mphys';
%x_axis_vals = 'LWP+RWP GCM grid-box mean'; 
%x_axis_vals = 'LWP+precip removed LWP GCM grid-box mean';
%x_axis_vals = 'In-cloud LWP normalised by GCM CF'; %LWP divided by model CF
%x_axis_vals = 'In-cloud LWP normalised by COSP-CALIPSO CF'; %LWP divided by COSP-CALIPSO CF
%x_axis_vals = 'In-cloud LWP GCM min COSP MODIS CF'; %LWP divided by an estimated minimum possible MODIS COSP CF
%x_axis_vals = 'Precip rate';
%x_axis_vals = 'LWP removed by precip rate';
x_axis_vals = 'GOES LWP x-axis2';
%x_axis_vals = 'Dummy data';

%% Choose y variable

y_axis_vals ='Nd lat-lon cells';
y_axis_vals ='Time of Day';
% y_axis_vals='WMOD';
% y_axis_vals='CF';
% %y_axis_vals='Scattering angle';
% y_axis_vals='SZA std dev timeseries';
%%%y_axis_vals='Mean SZA timeseries';
%y_axis_vals='Mean SZA timeseries3';
%y_axis_vals='Mean CTT timeseries3, y-axis';
%y_axis_vals='Min CTT timeseries3, y-axis';
y_axis_vals='Mean CTH timeseries3, y-axis';
%y_axis_vals = 'SST - Mean CTT timeseries3, y-axis';
%y_axis_vals='Mean CTP timeseries3, y-axis';
%y_axis_vals='Mean Sensor Zenith Angle timeseries3';
% %y_axis_vals='Max Sensor Zenith Angle timeseries3';
% %y_axis_vals='SZA std dev timeseries2';
% %y_axis_vals='SZA std dev timeseries3';
% %y_axis_vals = 'Nd from histogram vals timeseries';
%y_axis_vals = 'Nd from histogram vals timeseries3';
%y_axis_vals = 'LWP std dev from histogram vals timeseries3';
% y_axis_vals='SZA max timeseries3';
%y_axis_vals='SZA min timeseries3';
%y_axis_vals='SZA min max difference';
%y_axis_vals='Mean Scattering Angle timeseries3';
%y_axis_vals='Min Scattering Angle timeseries3';
%y_axis_vals='Max Scattering Angle timeseries3';
%y_axis_vals = 'SZA max from composite';
y_axis_vals = 'Dummy data for 1D';
%y_axis_vals = 'Latitude';
%y_axis_vals = 'Longitude';
%y_axis_vals = 'Nd from grid vals timeseries3';
%y_axis_vals = 'Nd from grid vals re*0.8 timeseries3';
%y_axis_vals = 'Nd from individual pixels timeseries3';
%y_axis_vals = 'Nd from individual pixels 1.6\mum timeseries3';
%y_axis_vals = 'Nd from individual pixels 3.7\mum timeseries3';
%y_axis_vals = 'Cloud Fraction from grid vals timeseries3';
y_axis_vals = 'MOD35 Cloud Fraction from grid vals timeseries3';
%y_axis_vals = 'CF COSP-MODIS GCM';
%y_axis_vals = 'CF COSP-CALIPSO GCM';
%y_axis_vals = 'CF GCM';
%y_axis_vals = 'LWP GCM Grid-box mean pre-mphys';
%y_axis_vals = 'LWP removal timescale';
%y_axis_vals = 'Precip rate';
%y_axis_vals = 'Tau^{1/2}'; %timeseries3
%y_axis_vals = 'Tau'; %timeseries3
%y_axis_vals = 'Homogeneity Parameter timeseries3 using W from mean tau and Re';
%y_axis_vals = 'Homogeneity Parameter timeseries3 using mean W for each pixel';
%y_axis_vals = 'Homogeneity Parameter timeseries3 using Cahalan log mean W (pixel level)';
% y_axis_vals = 'Heterogeneity Parameter Cahalan Optical Depth (Seethala)';
%y_axis_vals = 'LWP standard deviation';
%y_axis_vals = 'Normalised LWP standard deviation';
%y_axis_vals = 'Homogeneity Parameter Cahalan Optical Depth (Seethala)';
%y_axis_vals = 'Cloud Top Temp standard deviation (normalised by mean CTT), liquid pixels';
%y_axis_vals = 'Cloud Top Temp standard deviation, liquid pixels';
%y_axis_vals = 'Cloud Top Temp standard deviation, all pixels';
%y_axis_vals = 'CTT homogeneity, liquid pixels';
%y_axis_vals = 're_{3.7} - re_{2.1} (\mum) reduced dataset Re_1.6 Re_3.7';
%y_axis_vals = 'Cloud Top Height all clouds';
%y_axis_vals = 'R_{eff 2.1 \mum} (\mum)';  %pick this one if using L3 data (i.e. no 1.6 and 3.7)
%y_axis_vals = 'R_{eff 2.1 \mum} (\mum) minus 20%';  %pick this one if using L3 data (i.e. no 1.6 and 3.7)
%y_axis_vals = 'CDR Polder2'; %uses the mean over the orbits - only one
%datapoint per day
%y_axis_vals = 'Re COSP GCM'; %timeseries3
%y_axis_vals = 'Nd GCM';
%y_axis_vals = 'TLWP DAYTIME GCM, grid-box average';
%y_axis_vals = 'TLWP NIGHTTIME GCM, grid-box average';
%y_axis_vals = 'LWP DAYTIME GCM, grid-box average';
%y_axis_vals = 'LWP NIGHTTIME GCM, grid-box average';
%y_axis_vals = 'AMSRE TLWP DAYTIME';
%y_axis_vals = 'AMSRE TLWP NIGHTTIME';
%y_axis_vals = 'MODIS LWP, grid-box average using MOD06';
%y_axis_vals = 'MODIS LWP, grid-box average using MOD35';
%y_axis_vals = 'MODIS LWP, grid-box average using MOD35 minus 15%';
%y_axis_vals = 'SZA Polder Colocated with MODIS'; %uses the mean over the orbits - only one datapoint per day
%y_axis_vals = 'SZA Polder select region'; %uses the mean over the orbits - only one datapoint per day
%y_axis_vals = 'CDP Nd'; %timeseries3
%y_axis_vals = 'TOA albedo';
%y_axis_vals = 'Precip rate Comstock using MODIS';
%y_axis_vals = 'Precip rate Comstock using AMSRE';
%y_axis_vals = 'POLDER reff daymean';
%y_axis_vals = 'MODIS minus POLDER mock L3';
%y_axis_vals = 'GOES LWP';
y_axis_vals = 'GOES Nd2';
%y_axis_vals = 'UM LWP';
%y_axis_vals = 'Precip rate DAYTIME GCM';
%y_axis_vals = 'Precip rate NIGHTTIME GCM';

%y_axis_vals = 'General y-axis no ilat';

% ---- Z variables ---------------------------------------------------------
z_axis_vals = 'Nd from grid vals timeseries3';
%z_axis_vals = 'Nd absolute uncertainty from grid vals timeseries3';
z_axis_vals = 'Tau';

else

%    clear ioverride_pdf_varchoose

end

thresh=150;
thresh=0.3;
thresh=0;

%cf = Cloud_Fraction_Liquid.data;
%WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
%tau = Cloud_Optical_Thickness_Liquid_Mean.data;
%                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres



%cf_time=Cloud_Fraction_Liquid.timeseries(ilat,:);
%NP_time=Cloud_Fraction_Liquid_Pixel_Counts.timeseries(ilat,:)./cf_time;

%tau_time = Cloud_Optical_Thickness_Liquid_Mean.timeseries(ilat,:);
%reff_time = Cloud_Effective_Radius_Liquid_Mean.timeseries(ilat,:)*1e-6; %convert to metres
%Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)

%[N_time,H,W,k,Q,cw]=MODIS_N_H_func(tau_time,reff_time,Wflag,0);

switch datatype
    case 'timeseries33'
        %timeseries3
        cf_time3=Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3(ilat,ilon,itime)./cf_time3;
        %sensZA_time3 = Sensor_Zenith_Mean.timeseries3(ilat,ilon,itime);
        max_sensZA_time3 = Sensor_Zenith_Maximum.timeseries3(ilat,ilon,itime);
        maxSZA_time3 = Solar_Zenith_Maximum.timeseries3(ilat,ilon,itime);
        minSZA_time3 = Solar_Zenith_Minimum.timeseries3(ilat,ilon,itime);

        tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
        reff_time3 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime)*1e-6; %convert to metres
        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)

        [N_time3,H,W,k,Q,cw]=MODIS_N_H_func(tau_time3,reff_time3,Wflag,0);



end



%%%%
switch datatype
    case 'timeseries3'
        %extra_title_info = [' for day ' modis_day_str ', Y' modis_year_str];
        %        extra_title_info = [' for Y' modis_year_str2];

        if length(itime)==size_time_dim
            tim_str='';
        else
            tim_str=[' for days ' num2str(itime_orig(1)) ' to ' num2str(itime_orig(end))];
        end

        tim_str = [' ' date_range_str];
        tim_str = [' ' time_mean_str];        
        
    otherwise
        tim_str='';
        modis_year_str2='';
end

switch datatype
    case 'makeshift'
        LAT_str='';
        LON_str='';
        
    case 'gcm_data'
        if length(iregion_lin)==1
            LAT_str=num2str(Plat3D(iregion_lin),'%2.2f');
        else
            LAT_str=[num2str(min(Plat3D(iregion_lin))-dlat/2,'%2.2f') ' to ' num2str(max(Plat3D(iregion_lin))+dlat/2,'%2.2f')];
        end

         if length(iregion_lin)==1
            LON_str=num2str(Plon3D(iregion_lin),'%2.2f');
        else
            LON_str=[num2str(min(Plon3D(iregion_lin))-dlon/2,'%2.2f') ' to ' num2str(max(Plon3D(iregion_lin))+dlon/2,'%2.2f')];
        end
        

    otherwise

        if length(ilat)==1
            LAT_str=num2str(LAT(ilat),'%2.2f');
        else
            LAT_str=[num2str(LAT(ilat(1))-dlat/2,'%2.2f') ' to ' num2str(LAT(ilat(end))+dlat/2,'%2.2f')];
        end

        if length(ilon)==1
            LON_str=num2str(LON(ilon));
        else
            LON_str=[num2str(min(LON(ilon))-dlon/2) ' to ' num2str(max(LON(ilon))+dlon/2)];
        end

end




%% switch for x-variable

switch x_axis_vals

     case 'General GCM-style x-axis simple2'
         
%         switch xlabelstr
%             case 'Dummy data'
%                 xlabelstr = 'x values from outside script';
%         end
    
    %Simple case where specify everything outside in the driver script with
    %no screening, region or time selection, etc.
      
        
%        X = X(ilat,ilon,itime);
        X = X_driver;    %choose only the requested lat lon points  
%         if length(size(Plat3D))==3
%             X=permute(X,[2 3 1]);
%         end
%        a=NaN*ones(size(X));
%        a(iregion_lin)=0;
%        X=X+a;
       
       pdf2D_X_save = X;
        
  
        
        Xbins = Xbins_driver;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
%        ikeep_X_above_zero=1;               
 
%want to do this afterwards
%        ipost_plotTime_commands = 0; post_plotTime_commands = {''};
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};


case 'General GCM-style x-axis simple'
    
    %Simple case where specify everything outside in the driver script with
    %no screening, region or time selection, etc.



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = X_driver;    %choose only the requested lat lon points  
        if length(size(Plat3D))==3
            X=permute(X,[2 3 1]);
        end
       a=NaN*ones(size(X));
       a(iregion_lin)=0;
       X=X+a;
       
       pdf2D_X_save = X;
        
  
        
        Xbins = Xbins_driver;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
%        ikeep_X_above_zero=1;               
 
%want to do this afterwards
%        ipost_plotTime_commands = 0; post_plotTime_commands = {''};
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};


case 'General GCM-style x-axis simple 4D'
    
    %Simple case where specify everything outside in the driver script with
    %no screening or time selection, etc.
    % But does now do region selection



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = X_driver;    %choose only the requested lat lon points  
%         if length(size(Plat3D))==3
%             X=permute(X,[2 3 1]);
%         end
%        a=NaN*ones(size(X));
%        a(iregion_lin)=0;
%        X=X+a;

      a = spatial_screening_4D(X,LAT_val,LON_val,Plat2D,Plon2D);
      X = X + a;
       
       pdf2D_X_save = X;
        
  
        
        Xbins = Xbins_driver;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
%        ikeep_X_above_zero=1;               
 
%want to do this afterwards
%        ipost_plotTime_commands = 0; post_plotTime_commands = {''};
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};

case 'General reg lat lon screening x-axis simple 4D'
    
    %Simple case where specify everything outside in the driver script with
    %no screening or time selection, etc.
    % But does now do region selection



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = X_driver;    %choose only the requested lat lon points  
%         if length(size(Plat3D))==3
%             X=permute(X,[2 3 1]);
%         end
%        a=NaN*ones(size(X));
%        a(iregion_lin)=0;
%        X=X+a;

%This method works for irregular lat long grids, but uses a lot of memory%
%      a = spatial_screening_4D(X,LAT_val,LON_val,Plat2D,Plon2D)%;
%      X = X + a;
      
% So, will use this to work with regulat lat/lon grids when need to save
% memory
 ilat = find(Plat2D(:,1)>=LAT_val(1) & Plat2D(:,1)<LAT_val(end));
 ilon = find(Plon2D(1,:)>=LON_val(1) & Plon2D(1,:)<LON_val(end));

 X = X(ilat,ilon,:,:);      
       
       pdf2D_X_save = X;
        
  
        
        Xbins = Xbins_driver;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
%        ikeep_X_above_zero=1;               
 
%want to do this afterwards
%        ipost_plotTime_commands = 0; post_plotTime_commands = {''};
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};



    case 'GOES LWP x-axis'   
        %Use read_GOES_vocals_netcdf_files.m to get the goes data
        
        xlabelstr = ['LWP (g m^{-2})'];
        

         X = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2

        
%        Y(Y<30)=NaN;

         X(isnan(X))=0;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        goes_LWP_save = X;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     



        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
        Xbins = [-0.01 30:10:2500]; ichoose_Ybins=1;      
        Xbins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Xbins=1;
%        Xbins = [-0.01 10.^[log10(30):0.1:log10(500)]]; ichoose_Xbins=1;   

        Xbins = Xbins_DRIVER; ichoose_Xbins=1;
        
        

        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
    case 'GOES LWP x-axis2'   
        %Use read_GOES_vocals_netcdf_files.m to get the goes data
        
        xlabelstr = ['LWP (g m^{-2})'];
        

%         X = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2
         
         X = goes_LWP_multi{79};

        
%        Y(Y<30)=NaN;

         X(isnan(X))=0;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        goes_LWP_save = X;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     



        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
        Xbins = [-0.01 30:10:2500]; ichoose_Xbins=1;      
        Xbins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Xbins=1;
        Xbins = [-0.01 10.^[log10(5):0.05:log10(300)]]; ichoose_Xbins=1;        
        

        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        
     case 'Cloud Top Temp standard deviation, liquid pixels'
        xlabelstr = 'Cloud Top Temp standard deviation, liquid pixels (K)';
        xlabelstr = '\sigma_{CTT} (K)';
        X = (Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3(ilat,ilon,itime));

        minXbins=0;
        maxXbins=20; 
%        maxYbins=6;

Xbins=[0:0.25:3]; ichoose_Xbins=1;
%Xbins=[0:0.25:2]; ichoose_Xbins=1;        

nXpdf = 125;
nXpdf = 25;

        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        
%        het = 1 - homog_tau_Cahalan(ilat,ilon,itime);
%        iremove = find(~(het < 0.05));
%        X(iremove)=NaN;

        
     case 'Precip rate'
        xlabelstr = 'Precipitation rate (mm day^{-1})';
        CF_gcm_thresh=[-0.01 1.01];
        %CF_gcm_thresh=0.8;
        
      %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
%               LWP = eval(['gcm_lwp_' gcm_str]);  %kg/m2
            precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s
            %only calcuate when precip greater than a threshold - otherwise
            %the zero precip times will swamp the answer
            %Perhaps better plotted vs precip rate? Or a PDF might be
            %better
%            precip_thresh = 1e-10; max_timescale = 4*1;
%            iprecip = find(precip<precip_thresh);
%            precip(iprecip)=precip_thresh;

max_bin_val = 0.2*24;

            X = precip*3600*24;  %mm/day
            X(X>max_bin_val)=max_bin_val;
            X=permute(X,[2 3 1]);
        
        a=NaN*ones(size(X));                
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
 
dbin=0.002*24;
Xbins = [-0.00001:dbin:max_bin_val-0.00001+dbin]; ichoose_Ybins=1; 


  case 'LWP removed by precip rate'
        xlabelstr = 'LWP removed in 30 mins (g m^{-2})';
        CF_gcm_thresh=[-0.01 1.01];
        %CF_gcm_thresh=0.8;
        
      %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
%               LWP = eval(['gcm_lwp_' gcm_str]);  %kg/m2
            precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s
            %only calcuate when precip greater than a threshold - otherwise
            %the zero precip times will swamp the answer
            %Perhaps better plotted vs precip rate? Or a PDF might be
            %better
%            precip_thresh = 1e-10; max_timescale = 4*1;
%            iprecip = find(precip<precip_thresh);
%            precip(iprecip)=precip_thresh;

max_bin_val = 0.2*0.5*1e3;
dbin=0.002*0.5*1e3;

            X = precip*3600*0.5*1e3;  %mm/hr * 30 mins * 1e3 =  g/m2 removed in 30 mins
            X(X>max_bin_val)=max_bin_val;
            X=permute(X,[2 3 1]);
        
        a=NaN*ones(size(X));                
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
 

Xbins = [-0.00001:dbin:max_bin_val-0.00001+dbin]; ichoose_Ybins=1; 



    case 'CALIPSO cloud fraction'
        daynight_str = 'Daytime';
%        daynight_str = 'Nighttime';
        daynight_str = 'DayNight average';
        
        xlabelstr = ['CALISPO' daynight_str 'Cloud Fraction'];
        
        switch daynight_str
            case 'Daytime'
                %X = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
                X = cllcalipso_monthly/100;
            case 'DayNight average'
                X =  cllcalipso_monthly_AVERAGE/100;
        end
        
        
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;

        
        
        %         minYbins = 0;
        %         maxYbins = 1;

        %Xbins = [0:0.1:1.0 1.00001]; ichoose_Xbins=1;
        Xbins = [0:0.1:1.00001]; ichoose_Xbins=1;
        Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1; %the problem here is that the bin widths do not reflect the actual
        %bin widths, since cannot have CF<0 or >1

        Xbins = [-0.0001:0.1:1.0001]; ichoose_Xbins=1;
        Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1;
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
      case 'MOD35 Cloud Fraction from grid vals timeseries3'         
         xlabelstr = 'MOD35 Daytime Cloud Fraction';        
         X = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
         
%         minYbins = 0;
%         maxYbins = 1;

%Xbins = [0:0.1:1.0 1.00001]; ichoose_Xbins=1;
Xbins = [0:0.1:1.00001]; ichoose_Xbins=1;
Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1; %the problem here is that the bin widths do not reflect the actual
        %bin widths, since cannot have CF<0 or >1
        
Xbins = [-0.0001:0.1:1.0001]; ichoose_Xbins=1;
Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1;

 ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data

    
    case 'CF COSP-CALIPSO GCM'
        xlabelstr = 'COSP-CALIPSO GCM Cloud Fraction';
        CF_gcm_thresh=0.01;
        %CF_gcm_thresh=0.8;
        
        
            
        X = eval(['cllcalipso_' gcm_str '/100']);
        X=permute(X,[2 3 1]);
        
        a=NaN*ones(size(X));                
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
%        Y(Y<CF_gcm_thresh)=NaN;

%COSP CFs are multiples of 0.1 (0, 0.1, 0.2 ... 1.0)
Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1;
%Ybins = [-0.0999:0.1:0.9999 1.0001]; ichoose_Ybins=1; %[-0.0999    0.0001    0.1001    0.2001    0.3001    0.4001    0.5001    0.6001    0.7001    0.8001    0.9001    1.0001]
%this will ensure that the each CF value falls into one bin only (0, 0.1,
%0.2, etc.)

Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1; %this makes sure that the bin widths actually represent
%the ranges possible. This makes sense for MODIS. But for COSP I'm not sure
%what the singlar values actually mean - does 0.1 mean 0.05 to 0.15? Or
%exactly 0.1??
% so here we have 0 <= CF < 0.05, 0.05 <= CF < 0.1, ....  0.95 <= CF <= 1.0






    case 'CF GCM'
        dat_modis = eval(['cf_isccp_low_' gcm_str]);    
        
        switch daynight
            case 'Daytime'            
                xlabelstr = 'Max low cloud fraction (no screening) DAYTIME';
                times_required = [12:15];  %Aqua daytime
            case 'Nighttime'
                xlabelstr = 'Max low cloud fraction (no screening) NIGHTTIME';
                times_required = [0:3]; %Aqua nighttime
        end
        
        
         var_gcm_thresh=[-0.01 1e9]; %g m-2                

        xlabelstr = [xlabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];                
       var = eval(['1e3*gcm_lwp_' gcm_str]);  
       inan = find(var<var_gcm_thresh(1) | var>=var_gcm_thresh(2));
                      
        X = eval(['cf_isccp_low_' gcm_str]);
        X(inan) = NaN;
        X=permute(X,[2 3 1]);
        
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
         ioverride_time_selection=1;         
        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        X=X+permute(time_inds_average2,[2 3 1]);
        
        

        Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1;
        Xbins = [-.0001:0.1:0.8999 1.0001]; ichoose_Xbins=1;
        Xbins = [-0.1 0.001:0.1:0.9001 0.999 1.1]; ichoose_Xbins=1;
        Xbins = [0:0.05:0.95 1.001]; ichoose_Xbins=1; 
        Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1; %same as for model COSP bins
        
        Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1; %this makes sure that the bin widths actually represent
%the ranges possible. This makes sense for MODIS. But for COSP I'm not sure
%what the singlar values actually mean - does 0.1 mean 0.05 to 0.15? Or
%exactly 0.1??
% so here we have 0 <= CF < 0.05, 0.05 <= CF < 0.1, ....  0.95 <= CF <= 1.0

        extra_title_info = [extra_title_info ' ' gcm_str];
        
        
    case 'CF COSP-MODIS GCM'
        xlabelstr = 'COSP-MODIS GCM Cloud Fraction';
        CF_gcm_thresh=[-0.01 1.01];
        %CF_gcm_thresh=0.8;
        
        %N.B. - don't think there is a COSP MOD35 and MOD06 - seems
        %unlikely that COSP can determine edge pixels since it assumes a
        %horizontally homogeneous cloud
            
        X = eval(['liqCF_modis_' gcm_str])/100;
        X=permute(X,[2 3 1]);
        
        a=NaN*ones(size(X));                
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
 
%COSP CFs are multiples of 0.1 (0, 0.1, 0.2 ... 1.0)
Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1;
%Ybins = [-0.0999:0.1:0.9999 1.0001]; ichoose_Ybins=1; %[-0.0999    0.0001    0.1001    0.2001    0.3001    0.4001    0.5001    0.6001    0.7001    0.8001    0.9001    1.0001]
%this will ensure that the each CF value falls into one bin only (0, 0.1,
%0.2, etc.)

Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1; %this makes sure that the bin widths actually represent
%the ranges possible. This makes sense for MODIS. But for COSP I'm not sure
%what the singlar values actually mean - does 0.1 mean 0.05 to 0.15? Or
%exactly 0.1??
% so here we have 0 <= CF < 0.05, 0.05 <= CF < 0.1, ....  0.95 <= CF <= 1.0

    case 'LWP removal timescale'
        xlabelstr = 'LWP removal timescale (hrs)';
        CF_gcm_thresh=[-0.01 1.01];
        %CF_gcm_thresh=0.8;
        
      %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
               LWP = eval(['gcm_lwp_' gcm_str]);  %kg/m2
            precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s
            %only calcuate when precip greater than a threshold - otherwise
            %the zero precip times will swamp the answer
            %Perhaps better plotted vs precip rate? Or a PDF might be
            %better
            precip_thresh = 1e-10; max_timescale = 24*7;
%            precip(precip<precip_thresh)=precip_thresh;
            
            X = LWP./precip/3600; 
            X(precip<precip_thresh)=max_timescale;
            X=permute(X,[2 3 1]);
        
        a=NaN*ones(size(X));                
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
 

Xbins = [-0.00001:0.2:max_timescale+0.00001]; ichoose_Xbins=1; 




    case 'qv700 GCM'

        
        xlabelstr = '700 mb water vapour mixing ratio (g kg^{-1})';
        
%                 CF_gcm_thresh=0.01;
% %        CF_gcm_thresh=0.8;
% 
%         cf = eval(['liqCF_modis_' gcm_str])/100;
%         cf(cf<CF_gcm_thresh)=NaN;


%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];

        X = eval(['1e3*gcm_qv700_' gcm_str]);    %convert to g/kg form kg/kg
        %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        %keep only the points that are ocean AND are in our required region
        inds_keep  = intersect(iregion_lin,iregion_ocean_lin);
        a(inds_keep)=0;
        X=X+a;
        
        Xbins = [0:0.1:20]; ichoose_Xbins=1;  
%        Xbins = [0:0.1:1]; ichoose_Xbins=1;          
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
    case 'LTS GCM'

        
        xlabelstr = 'LTS (K)';
        
%                 CF_gcm_thresh=0.01;
% %        CF_gcm_thresh=0.8;
% 
%         cf = eval(['liqCF_modis_' gcm_str])/100;
%         cf(cf<CF_gcm_thresh)=NaN;


%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];

        X = eval(['gcm_LTS_' gcm_str]);    %convert to g/kg form kg/kg
        %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        Xbins = [0:1:30]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        
     case 'CDP Nd'

        
        xlabelstr = 'CDP Nd (cm^{-3})';       
        X = squeeze(Nd_CDP.timeseries3(ilat,ilon,itime));    %use the indices to remove the data at the end, which
        %is blank
        
        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data    
        
        
 case 'CDP Nd all data'        
        xlabelstr = 'CDP Nd all data (cm^{-3})';      
        
        [Y,M,D] = datevec(MatlabTime_Puijo);                    
        X = squeeze(Nd_Puijo);
        CDP_all_month = 'Oct'; month_choose=10;
        X(M~=month_choose)=NaN;
        
        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data  
        
        
 case 'days for CDP Nd all data'        
        xlabelstr = 'CDP Nd all data (cm^{-3})';               
        date_str = datestr(MatlabTime_Puijo);
        X = day_of_year_from_date_func(date_str);
        

        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data          
        
        
    case 'CDP Nd NaN match'
        
        inanCDP = find(isnan(Nd_CDP.timeseries3(ilat,ilon,itime))==1);
        inanMODIS = find(isnan(Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime))==1);
        inanDMPS = find(isnan(Nacc(ilat,ilon,itime))==1);       


        
        xlabelstr = 'CDP Nd (cm^{-3}) - reduced dataset';       
        
        X = squeeze(Nd_CDP.timeseries3(ilat,ilon,itime));    %use the indices to remove the data at the end, which
        %is blank
        
        X(inanMODIS)=NaN;
        X(inanDMPS)=NaN;
        
        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data    
        
        
     case 'DMPS Nacc'
        
        xlabelstr = 'Accumulation mode aerosol conc. (cm^{-3})';       
        X = squeeze(Nacc(ilat,ilon,itime));    %use the indices to remove the data at the end, which
        %is blank
        
        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data   
        
        
      case 'log DMPS Nacc'
        
        xlabelstr = 'log(N_{acc}) (cm^{-3})';       
        X = log(squeeze(Nacc(ilat,ilon,itime)));    %use the indices to remove the data at the end, which
        %is blank
        
        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data   
        
        
     case 'DMPS Nd estimate'
        
        xlabelstr = 'N_d calculated from Nacc (cm^{-3})';   
        Nacc_temp = squeeze(Nacc(ilat,ilon,itime));
       switch DMPS_type
            case 'Irshad new'
                a = 30.13;  b = 0.36; %Irshad values
            case 'Irshad'
                a = 45.08;  b = 0.31; %Irshad values    
            case 'marine'
                a = 8.75;  b = 0.61; %Portin, 2009, marine
            case 'continental'
                a = 4.06;  b = 0.7; %Portin, 2009, continental
            case 'all'
                a = 10.2;  b = 0.54; %Portin, 2009, all
           case 'Dan'
                a = 15.1781; b = 0.4408; %from a pfit to the cat 2 data with nlayers>1 and CTH>800 removed
        end
        
        
        Nd_DMPS.timeseries3 = a.*Nacc_temp.^b;  
        
        X = Nd_DMPS.timeseries3;
        %is blank
        
        Xbins = [0:25:1000]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data      
        
     case 'DMPS Nd estimate CDP NaN match'
         
        inanCDP = find(isnan(Nd_CDP.timeseries3(ilat,ilon,itime))==1);
        inanMODIS = find(isnan(Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime))==1);
        
        
        xlabelstr = 'N_d calculated from Nacc (cm^{-3}) - reduced dataset';   
        Nacc_temp = squeeze(Nacc(ilat,ilon,itime));
        Nacc_temp(inanCDP) = NaN; %make it so that the same time periods are used for the CDP and DMPS 
        Nacc_temp(inanMODIS) = NaN; %make it so that the same time periods are used for the CDP and DMPS 
        
        switch DMPS_type
            case 'Irshad'
                a = 45.08;  b = 0.31; %Irshad values
            case 'marine'
                a = 8.75;  b = 0.61; %Portin, 2009, marine
            case 'continental'
                a = 4.06;  b = 0.7; %Portin, 2009, continental
            case 'all'
                a = 10.2;  b = 0.54; %Portin, 2009, all
        end
        
        Nd_DMPS.timeseries3 = a.*Nacc_temp.^b;  
        
        X =  Nd_DMPS.timeseries3;
        %is blank
        
        Xbins = [0:5:1000]; ichoose_Xbins=0;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data      
        

        
        
    case 'Percentage Nd difference 1.6um allSZA vs lowSZA'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Percentage Nd difference';
        
%        if ~exist('ioverride_pdf')
 %           band_str='21';
            band_str='16';
 %           band_str='37';
 %       end
            
            switch band_str
                case '21'
                    band_str2='2.1 \mum';
                case '16'
                    band_str2='1.6 \mum';
                case '37'
                    band_str2='3.7 \mum';
            end
            
            title_nice = [LAT_str ', ' time_mean_str ', ' band_str2];

        X=eval(['squeeze(100*(Nd_' band_str '_allSZA(ilat,ilon,itime) - Nd_' band_str '_lowSZA(ilat,ilon,itime))./Nd_' band_str '_allSZA(ilat,ilon,itime));']);
        
%        dsza = eval(['squeeze(Solar_Zenith_Maximum_allSZA.timeseries3(ilat,ilon,itime) - Solar_Zenith_Maximum_lowSZA.timeseries3(ilat,ilon,itime));']);
%        X(dsza<1e-3)=NaN;
        
        x_axis_vals = [', ' band_str2];
        
%        Xbins = [-50:2:50]; ichoose_Xbins=1;
        
case 'Percentage Nd difference 2.1um allSZA vs lowSZA'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Percentage Nd difference';
        
%        if ~exist('ioverride_pdf')
            band_str='21';
%            band_str='16';
 %           band_str='37';
 %       end
            
            switch band_str
                case '21'
                    band_str2='2.1 \mum';
                case '16'
                    band_str2='1.6 \mum';
                case '37'
                    band_str2='3.7 \mum';
            end
            
            title_nice = [LAT_str ', ' time_mean_str ', ' band_str2];

        X=eval(['squeeze(100*(Nd_' band_str '_allSZA(ilat,ilon,itime) - Nd_' band_str '_lowSZA(ilat,ilon,itime))./Nd_' band_str '_allSZA(ilat,ilon,itime));']);
        
%        dsza = eval(['squeeze(Solar_Zenith_Maximum_allSZA.timeseries3(ilat,ilon,itime) - Solar_Zenith_Maximum_lowSZA.timeseries3(ilat,ilon,itime));']);
%        X(dsza<1e-3)=NaN;
        
        x_axis_vals = [', ' band_str2];
        
%        Xbins = [-50:2:50]; ichoose_Xbins=1;
        
case 'Percentage Nd difference 3.7um allSZA vs lowSZA'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Percentage Nd difference';
        
%        if ~exist('ioverride_pdf')
 %           band_str='21';
%            band_str='16';
           band_str='37';
 %       end
            
            switch band_str
                case '21'
                    band_str2='2.1 \mum';
                case '16'
                    band_str2='1.6 \mum';
                case '37'
                    band_str2='3.7 \mum';
            end
            
            title_nice = [LAT_str ', ' time_mean_str ', ' band_str2];

        X=eval(['squeeze(100*(Nd_' band_str '_allSZA(ilat,ilon,itime) - Nd_' band_str '_lowSZA(ilat,ilon,itime))./Nd_' band_str '_allSZA(ilat,ilon,itime));']);
        
%        dsza = eval(['squeeze(Solar_Zenith_Maximum_allSZA.timeseries3(ilat,ilon,itime) - Solar_Zenith_Maximum_lowSZA.timeseries3(ilat,ilon,itime));']);
%        X(dsza<1e-3)=NaN;
        
        x_axis_vals = [', ' band_str2];
        
        %Xbins = [-50:2:50]; ichoose_Xbins=1;
        
        
        
    case 'Percentage Nd difference allSZA vs lowSZA'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Percentage Nd difference';
        
        if ~exist('ioverride_pdf')
            band_str='21';
            %            band_str='16';
            band_str='37';
        end
            
            switch band_str
                case '21'
                    band_str2='2.1 \mum';
                case '16'
                    band_str2='1.6 \mum';
                case '37'
                    band_str2='3.7 \mum';
            end
            
title_nice = [LAT_str ', ' time_mean_str ', ' band_str2];

        X=eval(['squeeze(100*(Nd_' band_str '_allSZA(ilat,ilon,itime) - Nd_' band_str '_lowSZA(ilat,ilon,itime))./Nd_' band_str '_allSZA(ilat,ilon,itime));']);
        
        dsza = eval(['squeeze(Solar_Zenith_Maximum_allSZA.timeseries3(ilat,ilon,itime) - Solar_Zenith_Maximum_lowSZA.timeseries3(ilat,ilon,itime));']);
%        X(dsza<1e-3)=NaN;
        
        x_axis_vals = [', ' band_str2];
        
        Xbins = [-50:2:50]; ichoose_Xbins=1;
        
    case 'CF ice'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Ice Cloud Fraction';
        X = Cloud_Fraction_Ice.timeseries3(ilat,ilon,itime);

        
     case 'CF undet'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Undetermined Cloud Fraction';
        X = Cloud_Fraction_Undetermined.timeseries3(ilat,ilon,itime);
 
        
    case 'CF all'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Cloud Fraction (all phases)';
        X = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime)+Cloud_Fraction_Ice.timeseries3(ilat,ilon,itime)+Cloud_Fraction_Undetermined.timeseries3(ilat,ilon,itime);
%        Cloud_Fraction_Combined.timeseries(ilat,ilon,itime) %think
%        Combined is the same
       
        
    case 'CTP std dev timeseries3, x-axis'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Cloud Top Pressure std dev (hPa)';
        X = Cloud_Top_Pressure_Day_Standard_Deviation.timeseries3(ilat,ilon,itime);

        %        minXbins=100;
        %        maxXbins=1080;

    case 'Nd dist skewness lat-lon cells'
        ihtot = [1:prod(size(Nd_all(:,:,1)))]; thresh_str='xxx'; %all the data
        xlabelstr = 'Skewness';
        X = skewness(Nd_all,1,3);
    case 'dNd/dSZA overall'
        ihtot = [1:prod(size(grad_sza(:,:,itime)))]; thresh_str='xxx'; %all the data
        xlabelstr = 'dNd/dSZA overall';
        X = grad_sza(:,:,itime);
    case 'Sensor ZA timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Mean Sensor ZA';
        X = Sensor_Zenith_Mean.timeseries3(ilat,ilon,itime);
    case 'Latitude'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Latitude';
        X = squeeze(repmat(LAT(ilat),[1 1 length(ilon) size(Solar_Zenith_Standard_Deviation.timeseries3,3)]));

    case 'Longitude'
        xlabelstr = 'Longitude';

        switch datatype
            case 'gcm_data'
                X = Plon3D;  %Note that Plon3D doesn't need to be permuted
                a=NaN*ones(size(Plon3D)); a(iregion_lin)=0;
                X = X + a;
                       
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
                %MODIS lat and lons are the bin centres, so choose bins to be the edges between them
                %Tryin to make these consistent between plots - make them
                %odd numbers
                Xbins = [1+2*floor(minALL(X)/2):2:-1+2*ceil(maxALL(X)/2)]; ichoose_Xbins=1;
                
                lon_bin = 5;
                Xbins = [lon_bin*floor(minALL(X)/lon_bin):lon_bin:lon_bin*ceil(maxALL(X)/lon_bin)]; ichoose_Xbins=1;                
            otherwise
                X = squeeze(repmat(LON(ilon),[length(ilat) 1 length(itime)]));
                %MODIS lat and lons are the bin centres, so choose bins to be the edges between them
%                Xbins = [LON(ilon(1))-0.5:2:LON(ilon(end))+0.5]; ichoose_Xbins=1;
                
                lon_bin = 5;
                Xbins = [lon_bin*floor(minALL(X)/lon_bin):lon_bin:lon_bin*ceil(maxALL(X)/lon_bin)]; ichoose_Xbins=1;   
        end
        
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data

    case 'Std. dev of Reff'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Std Dev of Reff (\mum)';
        X = Cloud_Effective_Radius_Liquid_Standard_Deviation.timeseries3(ilat,ilon,itime);

    case 'Minimum Reff timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Minimum Reff (\mum)';
        X = Cloud_Effective_Radius_Liquid_Minimum.timeseries3(ilat,ilon,itime);

    case 'Maximum Tau timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Maximum Optical Thickness';
        X = Cloud_Optical_Thickness_Liquid_Maximum.timeseries3(ilat,ilon,itime);

    case 'Cloud fraction timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Mean Cloud Fraction';
        X = cf_time3;

        %        minXbins=-1e-10;
        %        maxXbins=1+1e-10;

    case 'Mean SZA timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Solar Zenith Angle';
        X = Solar_Zenith_Mean.timeseries3(ilat,ilon,itime);
        
    case 'Relative Azimuth Angle timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Relative Azimuth Angle';
        X = Relative_Azimuth_Mean.timeseries3(ilat,ilon,itime);
        
    case 'Mean Sensor ZA timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Sensor Zenith Angle';
        X = Sensor_Zenith_Mean.timeseries3(ilat,ilon,itime);    

    case 'Mean CTT timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Mean Cloud Top Temperature (K)';
        X = Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime);

        minXbins=200;
        maxXbins=310;


    case 'Mean CTP timeseries3'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        xlabelstr = 'Mean Cloud Top Pressure (hPa)';
        X = Cloud_Top_Pressure_Day_Mean.timeseries3(ilat,ilon,itime);

        minXbins=100;
        maxXbins=1080;




    case 'Nd from grid vals'
        xlabelstr = 'N_d (cm^{-3})';
        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
        xdat(1).x = N;

    case 'Nd from grid vals timeseries3'
        xlabelstr = 'N_d (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_time3(ilat,ilon,itime);
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=600;

Xbins = [0:10:1000]; ichoose_Xbins=1;

    case 'Nd from individual pixels 1.6\mum timeseries3'
         xlabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels using re_{1.6}';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Droplet_Number_Concentration_16.timeseries3(ilat,ilon,itime);
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=600;

%Xbins = [0:10:1000]; ichoose_Xbins=0;

 case 'Nd from individual pixels timeseries3'
         xlabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=600;

%Xbins = [0:10:1000]; ichoose_Xbins=1;

 case 'Nd from individual pixels 3.7\mum timeseries3'
         xlabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels using re_{3.7}';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime);
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=600;

%Xbins = [0:10:1000]; ichoose_Xbins=0;

case 'days for Nd from individual pixels 3.7\mum timeseries3'
         xlabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels using re_{3.7}';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        inotnan = find(isnan(Date_Time_Swath.timeseries3(ilat,ilon,itime))==0);
        X = NaN * ones(size(Date_Time_Swath.timeseries3(ilat,ilon,itime)));
        date_str = datestr(Date_Time_Swath.timeseries3(inotnan));        
        X(inotnan) = day_of_year_from_date_func(date_str);
                
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=600;

%Xbins = [0:10:1000]; ichoose_Xbins=0;
        




     case 'Nd from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'
        xlabelstr = 'N_d (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_time3(ilat,ilon,itime);
        %        X = N_time3;

        minXbins=0;
        maxXbins=600;       
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;
        
        
        
    case 'Nd_{1.6} from grid vals timeseries3'
        xlabelstr = 'N_{d Re 1.6} (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_time3_16(ilat,ilon,itime);
        %        X = N_time3;

        minXbins=0;
        maxXbins=600;
        
    case 'Nd_{1.6} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'
        xlabelstr = 'N_{d Re 1.6} (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_time3_16(ilat,ilon,itime);
        %        X = N_time3;

        minXbins=0;
        maxXbins=600;    
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;
        
     case 'Nd_{3.7} from grid vals timeseries3'
        xlabelstr = 'N_{d Re 3.7} (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_time3_37(ilat,ilon,itime);
        %        X = N_time3;
        
        std_dev = Standard_Deviation_Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime); 

        minXbins=0;
        maxXbins=600;  
        
    case 'Nd_{3.7} from grid vals timeseries3 reduced dataset Re_1.6 Re_3.7'
        xlabelstr = 'N_{d Re 3.7} (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_time3_37(ilat,ilon,itime);
        %        X = N_time3;

        minXbins=0;
        maxXbins=600;  
        
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;
        
    case 'Nd factor'
        xlabelstr = 'Pre-factor for Nd';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
      
        %        X = N_time3;

        minXbins=0;
        maxXbins=600;  
        
        
        Tau = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime); 
        
        N = N_time3_16(ilat,ilon,itime);
        
        %not dependent on whether use 1.6, 2.1 or 3.7 um since only depends
        %on the CTT
        X = N./(Tau.^0.5 .* Re_16.^(-2.5));
        
        %make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;
        
        
        
 case 'Nd uncertainty from grid vals timeseries3'
        xlabelstr = 'N_d uncertainty (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_un_time3(ilat,ilon,itime);
        %        X = N_time3;
        
        %N.B. - N_un_time3 = N_time3/100 .* sqrt(
        %(0.5*Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3).^2 + (-5/2*Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3).^2 );
                
        minXbins=0;
        maxXbins=600;
        
 case 'Nd 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)'
        xlabelstr = 'N_d uncertainty (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_un_time3(ilat,ilon,itime) .* N_time3_16(ilat,ilon,itime)./N_time3(ilat,ilon,itime);
%just need to multiply by the ratio of Nd for the different wavelengths
%since are assuming the same relative error in Re, since:- 
%N_un_time3 = N_time3/100 .* sqrt(
%(0.5*Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3).^2 + (-5/2*Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3).^2 );
%(i.e. relative Nd error is the same)

        minXbins=0;
        maxXbins=600;
        
case 'Nd 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)'
        xlabelstr = 'N_d uncertainty (cm^{-3})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = N_un_time3(ilat,ilon,itime) .* N_time3_37(ilat,ilon,itime)./N_time3(ilat,ilon,itime);
%just need to multiply by the ratio of Nd for the different wavelengths
%since are assuming the same relative error in Re, since:- 
%N_un_time3 = N_time3/100 .* sqrt(
%(0.5*Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3).^2 + (-5/2*Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3).^2 );
%(i.e. relative Nd error is the same)

        minXbins=0;
        maxXbins=600;        
        
        
 case 'Nd % uncertainty from grid vals timeseries3'
        xlabelstr = 'N_d uncertainty (%)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = 100*N_un_time3(ilat,ilon,itime)./N_time3(ilat,ilon,itime);
        %        X = N_time3;

        minXbins=0;
        maxXbins=200;      
        
case 'Re % uncertainty from grid vals timeseries3'
        xlabelstr = 'R_e uncertainty (%)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        %Re and Tau uncertainties are already in percent
        X = Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3(ilat,ilon,itime);

        minXbins=0;
        maxXbins=100;   
        
case 'Re uncertainty from grid vals timeseries3'
        xlabelstr = 'R_e uncertainty (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        %Re and Tau uncertainties are already in percent
        X = Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3(ilat,ilon,itime).*Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime)/100;

        minXbins=0;
        maxXbins=10;    
        
case 'Re 1.6 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)'
        xlabelstr = 'R_e uncertainty (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        %Re and Tau uncertainties are already in percent
        X = Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3(ilat,ilon,itime).*Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime)/100;

        minXbins=0;
        maxXbins=10;   
        
case 'Re 3.7 \mum uncertainty from grid vals timeseries3 (assumed 2.1 \mum percentage error)'
        xlabelstr = 'R_e uncertainty (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        %Re and Tau uncertainties are already in percent
        X = Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3(ilat,ilon,itime).*Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime)/100;

        minXbins=0;
        maxXbins=10;         
        
case 'Optical Depth % uncertainty from grid vals timeseries3'
        xlabelstr = 'Optical Depth Uncertainty (%)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        %Re and Tau uncertainties are already in percent
        X = Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3(ilat,ilon,itime);

        minXbins=0;
        maxXbins=200;  
        
    case 'Optical Depth uncertainty from grid vals timeseries3'
        xlabelstr = 'Optical Depth Uncertainty';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        %Re and Tau uncertainties are already in percent
        X = Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3(ilat,ilon,itime).*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime)/100;

        minXbins=0;
        maxXbins=200;          

        
     case 'H from grid vals timeseries3'
        xlabelstr = 'Cloud Depth (m)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = H_time3(ilat,ilon,itime);
        X = sqrt(1/1.15/0.7) * H_time3(ilat,ilon,itime); xlabelstr=[xlabelstr ' Painemal']; 
        %        X = N_time3;

        minXbins=0;
        maxXbins=1500; 
        
     case 'LWP grid-cell mean from grid vals timeseries3'
        cf_type = 'MOD06';
        cf_type = 'MOD35';
         
         xlabelstr = 'LWP (g m^{-2})';
        
        switch cf_type
            case 'MOD06'
                cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
            case 'MOD35'
                cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
        end
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        X = (1/1.15)*1e3*W_time3(ilat,ilon,itime).*cf; xlabelstr=[xlabelstr ' minus 15%, CF=' cf_type];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

Xbins = [0:10:500]; ichoose_Xbins=1;   
Xbins = [0:1:500]; ichoose_Xbins=1;   


        case 'LWP 1.6um MODIS grid-cell mean from mockL3'
            cf_type = 'MOD06';
            cf_type = 'MOD35';

             xlabelstr = 'LWP (g m^{-2})';

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_MODIS  = 5/9*1e3*Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime)*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);

            X = (1/1.15)*1e3*LWP_MODIS.*cf; xlabelstr=[xlabelstr ' minus 15%, CF=' cf_type];        
            %        X = N_time3;

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [0:10:500]; ichoose_Xbins=1;   
            Xbins = [0:1:500]; ichoose_Xbins=1;
            
        case 'LWP 2.1um MODIS grid-cell mean from mockL3'
            cf_type = 'MOD06';
            cf_type = 'MOD35';

             xlabelstr = 'LWP (g m^{-2})';

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_MODIS  = 5/9*1e3*Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime)*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);

            X = (1/1.15)*1e3*LWP_MODIS.*cf; xlabelstr=[xlabelstr ' minus 15%, CF=' cf_type];        
            %        X = N_time3;

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [0:10:500]; ichoose_Xbins=1;   
            Xbins = [0:1:500]; ichoose_Xbins=1;
            
         case 'LWP 3.7um MODIS grid-cell mean from mockL3'
            cf_type = 'MOD06';
            cf_type = 'MOD35';

             xlabelstr = 'LWP (g m^{-2})';

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_MODIS  = 5/9*1e3*Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime)*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);

            X = (1/1.15)*1e3*LWP_MODIS.*cf; xlabelstr=[xlabelstr ' minus 15%, CF=' cf_type];        
            %        X = N_time3;

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [0:10:500]; ichoose_Xbins=1;   
            Xbins = [0:1:500]; ichoose_Xbins=1;
            
         case 'MODIS 1.6um minus AMSRE from mockL3'
            cf_type = 'MOD06';
            cf_type = 'MOD35';
            
            colocate_POLDER=1; %flag to determine whether to match the data samping to that of POLDER or not.
            MODIS_minus_15pc=0;

             xlabelstr = 'LWP (g m^{-2})';
             xlabelstr=[xlabelstr ', CF=' cf_type];

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_MODIS  = 1e3.*cf.*5/9.*1e3.*Cloud_Effective_Radius_16_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime).*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
            if MODIS_minus_15pc==1
                LWP_MODIS = LWP_MODIS*(1/1.15);
                xlabelstr=[xlabelstr ' minus 15%'];
            end
            
            LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1);
            X = LWP_MODIS - LWP_AMSRE_PDF; 
            
            
            %        X = N_time3;
            
           
            
            if colocate_POLDER==1
                X(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
            end
       
                       

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [-100:10:100]; ichoose_Xbins=1;

             
  
            
             

            
         case 'MODIS 2.1um minus AMSRE from mockL3'
            cf_type = 'MOD06';
            cf_type = 'MOD35';
            
            colocate_POLDER=1; %flag to determine whether to match the data samping to that of POLDER or not.
            MODIS_minus_15pc=0;

             xlabelstr = 'LWP (g m^{-2})';
             xlabelstr=[xlabelstr ', CF=' cf_type];

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_MODIS  = 1e3.*cf.*5/9.*1e3.*Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime).*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
            if MODIS_minus_15pc==1
                LWP_MODIS = LWP_MODIS*(1/1.15);
                xlabelstr=[xlabelstr ' minus 15%'];
            end
            
            LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1);
            X = LWP_MODIS - LWP_AMSRE_PDF; 
            
            
            %        X = N_time3;
            
           
            
            if colocate_POLDER==1
                X(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
            end
       
                       

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [-100:10:100]; ichoose_Xbins=1;
             
             
          
            
         case 'MODIS 3.7um minus AMSRE from mockL3'
            cf_type = 'MOD06';
            cf_type = 'MOD35';
            
            colocate_POLDER=1; %flag to determine whether to match the data samping to that of POLDER or not.
            MODIS_minus_15pc=0;

             xlabelstr = 'LWP (g m^{-2})';
             xlabelstr=[xlabelstr ', CF=' cf_type];

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_MODIS  = 1e3.*cf.*5/9.*1e3.*Cloud_Effective_Radius_37_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime).*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
            if MODIS_minus_15pc==1
                LWP_MODIS = LWP_MODIS*(1/1.15);
                xlabelstr=[xlabelstr ' minus 15%'];
            end
            
            LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1);
            X = LWP_MODIS - LWP_AMSRE_PDF; 
            
            
            %        X = N_time3;
            
           
            
            if colocate_POLDER==1
                X(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
            end
       
                       

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [-100:10:100]; ichoose_Xbins=1;
            
          case 'POLDER minus AMSRE LWP'
            cf_type = 'MOD06';
            cf_type = 'MOD35';

             xlabelstr = 'LWP (g m^{-2})';

            switch cf_type
                case 'MOD06'
                    cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
                case 'MOD35'
                    cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
            end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

            LWP_POLDER  = 1e3.*cf.*5/9.*1e3.*daymean_Par2_CDR(ilat,ilon,itime).*1e-6.*Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
            LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1);
            X = LWP_POLDER - LWP_AMSRE_PDF; xlabelstr=[xlabelstr ', CF=' cf_type];        
            %        X = N_time3;

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [-100:10:100]; ichoose_Xbins=1;    
            
            
        case 'POLDER reff daymean'
%            cf_type = 'MOD06';
%            cf_type = 'MOD35';

             xlabelstr = 'POLDER r_e (\mum)';

%             switch cf_type
%                 case 'MOD06'
%                     cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
%                 case 'MOD35'
%                     cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
%             end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);


            X = daymean_Par2_CDR(ilat,ilon,itime);
            %        X = N_time3;

    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [0:0.1:31]; ichoose_Xbins=1;    
            Xbins = [0:1:31]; ichoose_Xbins=1;        
            
        case 'POLDER from structure reff daymean'
%            cf_type = 'MOD06';
%            cf_type = 'MOD35';

             xlabelstr = 'POLDER r_e (\mum)';

%             switch cf_type
%                 case 'MOD06'
%                     cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
%                 case 'MOD35'
%                     cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
%             end
            ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
    %        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);


%            X = daymean_Par2_CDR(ilat,ilon,itime);
            eval(['X = ' switchable_str{idat_multi} 'daymean_Par2_CDR(ilat,ilon,itime);']);


    %        minXbins=0;
    %        maxXbins=1500; 

            Xbins = [0:0.1:31]; ichoose_Xbins=1;    
            Xbins = [0:1:31]; ichoose_Xbins=1;                
            

    case 'MockL3 reff 1.6um timeseries3'
        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end
        
         xlabelstr = 'r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        X = Cloud_Effective_Radius_16_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        
        if colocate_POLDER==1
            X(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
        end
        
    case 'MockL3 reff 2.1um timeseries3'

        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end
         
        xlabelstr = 'r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        X = Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        Xbins = [0:1:31]; ichoose_Xbins=1;
        
        if colocate_POLDER==1
            X(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
        end
        
    case {'MockL3_no_conf reff 2.1um timeseries3','MockL3_no_conf reff 1.6um timeseries3','MockL3_no_conf reff 3.7um timeseries3'}

        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end
         
        xlabelstr = 'r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%       X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

%        dat = flipdim(mockL3_no_conf.Cloud_Effective_Radius_Liquid_Mean.timeseries3,1);

        switch x_axis_vals
            case 'MockL3_no_conf reff 1.6um timeseries3'
                X = mockL3_no_conf.Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];
            case 'MockL3_no_conf reff 2.1um timeseries3'
                X = mockL3_no_conf.Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];
            case 'MockL3_no_conf reff 3.7um timeseries3'
                X = mockL3_no_conf.Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];
        end

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        Xbins = [0:1:31]; ichoose_Xbins=1;
        
        if colocate_POLDER==1
             ilat_POL = find(LAT_MODIS>=LAT_val(1) & LAT_MODIS<LAT_val(end));
             ilon_POL = find(LON_MODIS>=LON_val(1) & LON_MODIS<LON_val(end));
             days = daynum_timeseries3(itime);
             years = modisyear_timeseries3(itime);
             itime_POL = NaN*ones([1 length(days)]);
             for ifind=1:length(days)
                 itime_POL(ifind) = find(daynum_timeseries3_MODIS==days(ifind) & modisyear_timeseries3_MODIS==years(ifind));
             end
%            dat = flipdim(daymean_Par2_CDR,1);

            %dat = daymean_Par2_CDR(ilat_POL,ilon_POL,itime_POL);            
            eval(['dat = ' switchable_str{idat_multi} 'daymean_Par2_CDR(ilat_POL,ilon_POL,itime_POL);']);            
            
             %Now X and dat should be the same except that dat is flipped
             %in the lat dimension (1st dim)
            dat = flipdim(dat,1);
            X(isnan(dat))=NaN;

        end
        
        
    case 'MockL3 reff 1.6um timeseries3 re.LT.20um'

        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end
         
        xlabelstr = 'r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%       X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

%        dat = flipdim(mockL3_no_conf.Cloud_Effective_Radius_Liquid_Mean.timeseries3,1);
        X = mockL3_conf_re16_re37_LT_20um.Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        Xbins = [0:1:31]; ichoose_Xbins=1;
        
        if colocate_POLDER==1
            %This script finds the required indices
            find_isnan_polder_mockL3_flipped
            X(inan_POL)=NaN;
        end
        
     case 'MockL3 reff 3.7um timeseries3 re.LT.20um'

        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end
         
        xlabelstr = 'r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%       X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

%        dat = flipdim(mockL3_no_conf.Cloud_Effective_Radius_Liquid_Mean.timeseries3,1);
        X = mockL3_conf_re16_re37_LT_20um.Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        Xbins = [0:1:31]; ichoose_Xbins=1;
        
        if colocate_POLDER==1
            %This script finds the required indices
            find_isnan_polder_mockL3_flipped
            X(inan_POL)=NaN;
        end    
        
     case 'MockL3 reff 2.1um timeseries3 re.LT.20um'

        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end
         
        xlabelstr = 'r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%       X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

%        dat = flipdim(mockL3_no_conf.Cloud_Effective_Radius_Liquid_Mean.timeseries3,1);
        X = mockL3_conf_re16_re37_LT_20um.Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        Xbins = [0:1:31]; ichoose_Xbins=1;
        
        if colocate_POLDER==1
            %This script finds the required indices
            find_isnan_polder_mockL3_flipped
            X(inan_POL)=NaN;
        end    
        
        
    case 'MockL3 reff 3.7um timeseries3'
        
        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end

        
        xlabelstr = 'MODIS 3.7\mum r_e (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        X = Cloud_Effective_Radius_37_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];        
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [0:0.1:31]; ichoose_Xbins=1;
        Xbins = [0:1:31]; ichoose_Xbins=1;

        if colocate_POLDER==1
            X(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
        end
        
        
    case 'MODIS minus AMSRE LWP bias'
        xlabelstr = 'LWP bias (g m^{-2})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
%        X = 1e3*W_time3(ilat,ilon,itime).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
%MOD35       
cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);        

%MOD06
%cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

LWP_MODIS_PDF = (1/1.15)*1e3*W_time3(ilat,ilon,itime).*cf; xlabelstr=[xlabelstr ' MODIS minus 15% '];        
        LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1);
        X = LWP_MODIS_PDF - LWP_AMSRE_PDF;

%        minXbins=0;
%        maxXbins=1500; 

        Xbins = [-100:10:100]; ichoose_Xbins=1;
        

     case 'Cloud Fraction from grid vals timeseries3'
        xlabelstr = 'Cloud Fraction';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        %        X = N_time3;

%        minXbins=0;
%        maxXbins=1;   

        Xbins = [-0.05:0.1:1.05]; ichoose_Xbins=1; %the problem here is that the bin widths do not reflect the actual
        %bin widths, since cannot have CF<0 or >1
        Xbins = [-0.0001:0.1:0.89999 1.0001]; ichoose_Xbins=1;
        Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1;


    case 'Nd from histogram vals timeseries'
        xlabelstr = 'N_d (cm^{-3}) from histogram timeseries';
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data

        [histo_output]=...
            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,ihtot) );
        xdat(1).x=histo_output.N_histo_mean;
        %                        xdat(1).x=histo_output.N_histo_std;

    case 'Nd from histogram vals timeseries3'
        xlabelstr = 'N_d (cm^{-3}) from histogram timeseries3';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        %        [histo_output]=...
        %            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );
        %                        xdat(1).x=histo_output.N_histo_mean;
        X = Nd_timeseries.mean(ilat,ilon,itime);

        %        minXbins=0;
        %        maxXbins=500;
   
        
    case 'R_{eff 2.1 \mum} (\mum)'
        if ~exist('ioverride_pdf') | ioverride_pdf==0
            colocate_POLDER=1;
        end

        xlabelstr = 'R_{eff 2.1 \mum} (\mum)';
        
        X = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);   
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
         ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
         
          Xbins = [0:0.5:50]; ichoose_Xbins=1; 
          
          
          
          
          if colocate_POLDER==1

              dat = daymean_Par2_CDR(ilat,ilon,itime);
              X(isnan(dat))=NaN;

          end
        
    case 'R_{eff 2.1 \mum} (\mum) reduced dataset Re_1.6 Re_3.7'

        xlabelstr = 'R_{eff 2.1 \mum} (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;

        

    case 'R_{eff 1.6 \mum} (\mum)'

        xlabelstr = 'R_{eff 1.6 \mum} (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
    case 'R_{eff 1.6 \mum} (\mum) reduced dataset Re_1.6 Re_3.7'

        xlabelstr = 'R_{eff 1.6 \mum} (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;

    case 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7'

        xlabelstr = 'R_{eff 3.7 \mum} (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;
        
 case 'R_{eff 3.7 \mum} (\mum)'

        xlabelstr = 'R_{eff 3.7 \mum} (\mum)';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
        
        
    case 'Reff^{-5/2}'

        xlabelstr = 'Reff ^{-5/2} (\mum ^{-2.5})';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = (Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime)).^(-5/2);        

    case 'Tau'

        xlabelstr = 'Mean Optical Depth';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
        
        Xbins = [0:0.2:50]; ichoose_Xbins=1;
        
    case 'Tau COSP GCM'                
        a = eval(['liqTau_modis_' gcm_str]);
        b = eval(['liqCF_modis_' gcm_str])/100;
        
        xlabelstr = 'Mean Optical Depth';


        CF_gcm_thresh=0.01;
%        CF_gcm_thresh=0.8;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];

        a(b<CF_gcm_thresh)=NaN;
        a=a./b; %divide by the cloud fraction to get the in-cloud Tau
        a=permute(a,[2 3 1]);

        
        
%        X = X(ilat,ilon,itime);
        X = a;    %choose only the requested lat lon points   
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
         ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data    
        
        Xbins = [0:0.2:50]; ichoose_Xbins=1;
        
    case 'In-cloud averaged AMSRE TLWP'
        cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.99 1.01];
        cf(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2))=NaN;
        
        X = 1e3*lwp_amsre_time3(ilat,ilon,itime,1)./cf; %last index is ascending(=1)/descending(=2)
               
        xlabelstr = 'AMSRE In-Cloud TLWP (weighted by MODIS CF)';
             
%        Xbins = [0:20:300]; ichoose_Xbins=1;      
        Xbins = [0:20:300]; ichoose_Xbins=1;      
        Xbins = [-60:10:300]; ichoose_Xbins=1;  
         
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        
case 'In-cloud averaged (using daytime MOD35) AMSRE TLWP'
        cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
%        CF_thresh=0.01;
        CF_gcm_thresh=[0.99 1.01];
        CF_gcm_thresh=[0.05 0.35];
        CF_gcm_thresh=[0.35 0.45];        
        
%        CF_gcm_thresh=[0.65 0.95];
%        CF_gcm_thresh=[0.95 1.01];        
        
%        CF_gcm_thresh=[0.399 0.901];   
        CF_gcm_thresh=[-0.01 1.01]; 
        
        cf(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2))=NaN;
        
        X = 1e3*lwp_amsre_time3(ilat,ilon,itime,1)./cf;  %last index is ascending(=1)/descending(=2)
               
        xlabelstr = 'AMSRE In-Cloud TLWP (weighted by MOD35 Daytime CF)';
             
%        Xbins = [0:20:300]; ichoose_Xbins=1;      
        Xbins = [0:20:300]; ichoose_Xbins=1;      
        Xbins = [-60:10:500]; ichoose_Xbins=1;  
        Xbins = [-200:10:1500]; ichoose_Xbins=1;  
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data        
        
        
case 'Grid-box mean AMSRE TLWP'    
        cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
%        CF_thresh=0.01;
        CF_gcm_thresh=[0.99 1.01];
        CF_gcm_thresh=[0.05 0.35];
        CF_gcm_thresh=[0.35 0.45];        
        
%        CF_gcm_thresh=[0.65 0.95];
%        CF_gcm_thresh=[0.95 1.01];        
        
%        CF_gcm_thresh=[0.399 0.901];        
        CF_gcm_thresh=[-0.01 1.01];   

%        cf(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2))=NaN;

amsre_daynight = 'daytime';
%amsre_daynight = 'nighttime';

amsre_datatype = 'time3';
amsre_datatype = 'normal';

switch amsre_daynight
    case 'daytime'
        X = 1e3*lwp_amsre_time3(ilat,ilon,itime,1);  %last index is ascending(=1)/descending(=2)
        xlabelstr = 'AMSRE Grid-box mean Daytime TLWP (g m^{-2})';
    case 'nighttime'
        X = 1e3*lwp_amsre_time3(ilat,ilon,itime,2);  %last index is ascending(=1)/descending(=2)
        xlabelstr = 'AMSRE Grid-box mean Nighttime TLWP (g m^{-2})';
end
   
             
%        Xbins = [0:20:300]; ichoose_Xbins=1;      
        Xbins = [0:20:300]; ichoose_Xbins=1;      
        Xbins = [-60:10:500]; ichoose_Xbins=1;  
        Xbins = [-0.0001 10.^([0:0.04:3])]; ichoose_Xbins=1;
        
        ikeep_X_above_zero=1;
         
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data        
        
%want to do this afterwards
        ipost_plotTime_commands = 1;
        post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};
        
        
        
 case 'Re COSP GCM'
        
        cf = eval(['liqCF_modis_' gcm_str])/100;        
%        tau = eval(['liqTau_modis_' gcm_str './cf']);
        re = eval(['liqRe_modis_' gcm_str './cf']);

        
        xlabelstr = 'COSP Reff';

        CF_gcm_thresh=0.01;
        CF_gcm_thresh=0.8;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];

        extra_title_info = [thresh_str ' '];
%        tau(cf<CF_gcm_thresh)=NaN;
%        tau=permute(tau,[2 3 1]);
        re(cf<CF_gcm_thresh)=NaN;
        re=permute(re,[2 3 1]);
        
        
%        X = X(ilat,ilon,itime);
        X = 1e6*re;    %choose only the requested lat lon points   
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        Xbins = [0:1:50]; ichoose_Xbins=1;        
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
    
    case 'LWP COSP GCM'
        
        
        tau = eval(['liqTau_modis_' gcm_str]);
        re = 1.6*eval(['liqRe_modis_' gcm_str]);
        cf = eval(['liqCF_modis_' gcm_str])/100;
        
        xlabelstr = 'COSP LWP';
        ihtot = [1:prod(size(b))]; thresh_str='xxx'; %all the data

        CF_gcm_thresh=0.01;
%        CF_gcm_thresh=0.8;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];

        tau(cf<CF_gcm_thresh)=NaN;
        tau=permute(tau,[2 3 1]);
        re(b<CF_gcm_thresh)=NaN;
        re=permute(re,[2 3 1]);
        
        
%        X = X(ilat,ilon,itime);
        X = 1000*5/9*1000.*re.*tau;    %choose only the requested lat lon points   
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        Xbins = [0:10:400]; ichoose_Xbins=1;        
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
     
     case 'Nd maxliq GCM'
        
        xlabelstr = 'Droplet Concentration at maxliq (cm^{-3})';

       
        cf = eval(['liqCF_modis_' gcm_str])/100;
        cf(cf<CF_gcm_thresh)=NaN;

        CF_gcm_thresh=0.01;
%        CF_gcm_thresh=0.8;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = Nd;    %calculate using gcm_process_lite  
%        X = eval(['gcm_Nd_maxliq_' gcm_str])/100; %note already restricted to CF>0.8
        
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        Xbins = [0:10:1000]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data    
        
        
        
    case 'In-cloud LWP normalised by GCM CF'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
 %       CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];    
%        CF_gcm_thresh=[0.399 0.901];
        CF_gcm_thresh=[0.65 0.95];
%        CF_gcm_thresh=[0.05 0.35];
        
        xlabelstr = ['In-cloud mean Liquid Water Path (g m^{-2}) for model CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
        cf = eval(['cf_isccp_low_' gcm_str]);
        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str './cf']);    %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
         %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
 case 'In-cloud LWP normalised by COSP-CALIPSO CF'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.99 1.01];
 %       CF_gcm_thresh=[0.01 1.01];  
        CF_gcm_thresh=[0.099 0.301];    
%        CF_gcm_thresh=[0.399 0.901];

        xlabelstr = ['In-cloud mean Liquid Water Path (g m^{-2}) for model CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
        cf = eval(['cllcalipso_' gcm_str '/100']);
        
        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str './cf']);    %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
         %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        


case 'LWP GCM'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];  
        CF_gcm_thresh=[0.65 0.95];   
        CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
%        CF_gcm_thresh=[0.05 0.35];
        
        xlabelstr = ['In-cloud mean Liquid Water Path (g m^{-2}) for COSP CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
        cf = eval(['liqCF_modis_' gcm_str])/100;
        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str './cf']);    %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        

case 'LWP GCM grid-box mean'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];  
        CF_gcm_thresh=[0.65 0.95];   
        CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
%        CF_gcm_thresh=[0.05 0.35];
        CF_gcm_thresh=[-0.01 1.01];       
%        CF_gcm_thresh=[0.05 1.01];               

        xlabelstr = ['Grid-box mean Liquid Water Path (g m^{-2}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
%       cf = eval(['liqCF_modis_' gcm_str])/100;
       cf = eval(['cf_isccp_low_' gcm_str]);  
       inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
%        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
% 

%LWP_postMP = eval(['lwpAP_isccp_low_' gcm_str]) + eval(['lwpAP_isccp_mid_' gcm_str]) + eval(['lwpAP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points
%LWP_preMP = eval(['lwpBP_isccp_low_' gcm_str]) + eval(['lwpBP_isccp_mid_' gcm_str]) + eval(['lwpBP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points

%LWP_postMP = LWP_postMP *1e3;  %convert to g/m2
%LWP_preMP = LWP_preMP *1e3;
%dLWP = LWP_postMP - LWP_preMP; %Change due to mphys (g/m2)
        
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str]);    %choose only the requested lat lon points  
%        X = LWP_postMP;
        X(inan) = NaN;
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        Xbins = [0:5:500]; ichoose_Xbins=1;    
        
        Xbins = [-0.0001 10.^([0:0.04:3])]; ichoose_Xbins=1;
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        ikeep_X_above_zero=1;
 
%want to do this afterwards
        ipost_plotTime_commands = 1;
        post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};

%        
%        
        

    case {'Pre-mphys LWP GCM grid-box mean','Post-mphys LWP GCM grid-box mean','Post minus pre-mphys LWP GCM grid-box mean','Fraction of LWP lost to mphys GCM grid-box mean'}
%         CF_gcm_thresh=0.01;
%         CF_gcm_thresh=[0.95 1.01];
% %        CF_gcm_thresh=[0.01 1.01];  
% %        CF_gcm_thresh=[0.099 0.301];  
%         CF_gcm_thresh=[0.65 0.95];   
%         CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
%        CF_gcm_thresh=[0.05 0.35];
%        CF_gcm_thresh=[-0.01 1.01];       
%        CF_gcm_thresh=[0.05 1.01];               

       
        
%screening is done at the end now 

%       cf = eval(['liqCF_modis_' gcm_str])/100;
%       cf = eval(['cf_isccp_low_' gcm_str]);  
%       inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
%        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];


        LWP_postMP = eval(['lwpAP_isccp_low_' gcm_str]) + eval(['lwpAP_isccp_mid_' gcm_str]) + eval(['lwpAP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points 
        LWP_preMP = eval(['lwpBP_isccp_low_' gcm_str]) + eval(['lwpBP_isccp_mid_' gcm_str]) + eval(['lwpBP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points 
  
        LWP_postMP = LWP_postMP *1e3;  %convert to g/m2      
        LWP_preMP = LWP_preMP *1e3;
        dLWP = LWP_postMP - LWP_preMP; %Change due to mphys (g/m2)
        % 

%         
        
%        X = X(ilat,ilon,itime);
%        X = eval(['1e3*gcm_lwp_' gcm_str]);    %choose only the requested
%        lat lon points  

        Xbins = [0:10:500]; ichoose_Xbins=1;   
        Xbins = [0:5:500]; ichoose_Xbins=1;    
        
        Xbins = [-0.0001 10.^([0:0.04:3])]; ichoose_Xbins=1;
        
        
        

switch x_axis_vals
    case 'Pre-mphys LWP GCM grid-box mean'
        X = LWP_preMP;
        xlabelstr = ['Pre-mphys grid-box mean Liquid Water Path (g m^{-2}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
    case 'Post-mphys LWP GCM grid-box mean'
        X = LWP_postMP;
        xlabelstr = ['Post-mphys grid-box mean Liquid Water Path (g m^{-2}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
    case 'Post minus pre-mphys LWP GCM grid-box mean'
        X = (LWP_postMP-LWP_preMP);
        xlabelstr = ['Post minus pre-mphys grid-box mean Liquid Water Path (g m^{-2}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];        
        dbin = 0.1;
        Xbins = [-100-dbin/2:0.1:100+dbin/2]; ichoose_Xbins=1; %This allows a bit centred on exactly zero
    case 'Fraction of LWP lost to mphys GCM grid-box mean'
        X = 100 * dLWP./LWP_preMP;
        xlabelstr = ['Percentage change of GCM grid-box mean LWP during mphys'];
        dbin = 10;
        Xbins = [-100-dbin/2:dbin:200+dbin/2]; ichoose_Xbins=1; %This allows a bit centred on exactly zero        
end

%        X(inan) = NaN;   %should try to do screening at the end
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        

        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        ikeep_X_above_zero=0;
 
%want to do this afterwards
        ipost_plotTime_commands = 0; post_plotTime_commands = {''};
%       ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};

%        




    
    case 'LWP removal rate due to mphys'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];  
        CF_gcm_thresh=[0.65 0.95];   
        CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
%        CF_gcm_thresh=[0.05 0.35];
        CF_gcm_thresh=[-0.01 1.01];       
%        CF_gcm_thresh=[0.05 1.01];               

        xlabelstr = ['LWP removal rate due to mphys (mm day^{-1}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
%       cf = eval(['liqCF_modis_' gcm_str])/100;
       cf = eval(['cf_isccp_low_' gcm_str]);  
       inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
%        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        LWP_postMP = eval(['lwpAP_isccp_low_' gcm_str]) + eval(['lwpAP_isccp_mid_' gcm_str]) + eval(['lwpAP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points 
        LWP_preMP = eval(['lwpBP_isccp_low_' gcm_str]) + eval(['lwpBP_isccp_mid_' gcm_str]) + eval(['lwpBP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points 
%        precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2
        %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
              
        time_step = 30/60/24; %model timestep = 30 mins. Convert to days.
        X = (LWP_preMP - LWP_postMP) / time_step;
%        X = 1e3*(LWP_postMP);        
        X(inan) = NaN;
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        Xbins = [0:5:500]; ichoose_Xbins=1;    
        
        Xbins = [-0.0001 10.^([0:0.04:3])]; ichoose_Xbins=1;
        
%        max_bin_val = 0.4*24;
%        dbin=0.002*24;
%        Xbins = [-max_bin_val-0.00001:dbin:max_bin_val-0.00001+dbin]; ichoose_Xbins=0;
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        ikeep_X_above_zero=1;
 
%want to do this afterwards
        ipost_plotTime_commands = 1;
        %Need to comment out the post_plot commands below if not using them
        post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};

%        
%        
        
        


case 'LWP+RWP GCM grid-box mean'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];  
        CF_gcm_thresh=[0.65 0.95];   
        CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
        CF_gcm_thresh=[0.05 1.01];
        CF_gcm_thresh=[-0.01 1.01];

        xlabelstr = ['Grid-box mean LWP+RWP (g m^{-2}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
%        cf = eval(['liqCF_modis_' gcm_str])/100;
%        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
%        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_TLWP_' gcm_str]);    %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        Xbins = [0:5:500]; ichoose_Xbins=1;    
        
        Xbins = [-0.0001 10.^([0:0.04:3])]; ichoose_Xbins=1;
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        ikeep_X_above_zero=1;               
 
%want to do this afterwards
        ipost_plotTime_commands = 0; post_plotTime_commands = {''};
        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};
        
case 'LWP+precip removed LWP GCM grid-box mean'
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];  
        CF_gcm_thresh=[0.65 0.95];   
        CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
%        CF_gcm_thresh=[0.05 0.35];
        CF_gcm_thresh=[-0.01 1.01];

        xlabelstr = ['Grid-box mean LWP+precip (g m^{-2}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        
%        cf = eval(['liqCF_modis_' gcm_str])/100;
%        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
%        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;        



        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        dt_gcm = 0.5*3600; %model timestep - assume 30 mins. But is it for CLUBB?
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str]);    %LWP in g/m2     
        precip = eval(['1e3*gcm_precT_' gcm_str]); %is in m/s originally
        %multiply by 1e3 to give mm/s = kg/m2/s
        X = X + dt_gcm.*1e3*precip; %1e3 here to covert to g/m2/s
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        Xbins = [0:5:500]; ichoose_Xbins=1;    
        
        Xbins = [-0.0001 10.^([0:0.04:3])]; ichoose_Xbins=1;
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        ikeep_X_above_zero=1;
        
        %want to do this afterwards
        ipost_plotTime_commands = 1;
        post_plotTime_commands = {'set(gca,''xscale'',''log'');','set(gca,''xlim'',[0.9 300])'};
        
        
        
    case 'In-cloud LWP GCM min COSP MODIS CF'                
%Choose a MODIS COSP CF range - then assume that the rest of the sub-grid box contains very low
%LWC cloud (call LWC=0). And that the COSP-MODIS fraction has tau=0.3. Then
%calc whether the mean tau would then be < 0.03 and therefore be undetected
%by CALISPO - since we know that in all cases the COSP-CALIPSO CF was the
%same as the COSP-MODIS CF. If so, then use the COSP-MODIS CF for
%normalisation
% Need to only select points that have a model CF, but no COSP CF - i.e.
% COSP removes their CF due to them being too thin - for both MODIS and
% CALIPSO
        
xlabelstr = ['In-cloud mean Liquid Water Path (g m^{-2}) for COSP CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];  
        CF_gcm_thresh=[0.65 0.95];   
        CF_gcm_thresh=[0.35 0.45];           
%        CF_gcm_thresh=[0.399 0.901];
%        CF_gcm_thresh=[0.05 0.35];

CF_mean = mean(CF_gcm_thresh);

cf_MODIS = eval(['liqCF_modis_' gcm_str])/100;            
cf_CALIPSO = eval(['cllcalipso_' gcm_str '/100']);
%model CF
cf = eval(['cf_isccp_low_' gcm_str]);
%points that were not removed by COSP when there was some model CF
inot_COSP_remove = find(~(cf_MODIS<0.05 & cf_CALIPSO<0.05 & cf > 0.05)) ;
        

inan = find(cf<=CF_mean); %need the model CF to be at least as big as our desired CF
cf(inan) = NaN;               
        
%        cf = eval(['liqCF_modis_' gcm_str])/100;
%        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
%        cf(inan) = NaN;
%        cf(cf<CF_gcm_thresh(1))=NaN;    

tau_thin = 0;
tau_modis = 0.3;
tau_calipso = 0.03;
tau_calipso = 0.2;

tau_mean = (tau_thin.*(cf - CF_mean) + tau_modis.*CF_mean) ./ cf;
ireject = find(tau_mean>tau_calipso);




        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];



% 
%         
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str './CF_mean']);    %choose only the requested lat lon points   
        X(ireject) = NaN; %reject those that were visible by CALIPSO
        X(inot_COSP_remove) = NaN; %reject those that would not have been removed by COSP
        %  since we only have the MODIS/CALIPSO constaints on these
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        Xbins = [0:10:500]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        
  case 'In-cloud TLWP GCM norm by COSP-MODIS-CF'
        
        xlabelstr = 'In-cloud mean LWP+RWP (g m^{-2})';
        

        
         CF_gcm_thresh=0.01;
        CF_gcm_thresh=[0.95 1.01];
%        CF_gcm_thresh=[0.01 1.01];  
%        CF_gcm_thresh=[0.099 0.301];    
%        CF_gcm_thresh=[0.399 0.901];  
        CF_gcm_thresh=[0.65 0.95];
%        CF_gcm_thresh=[0.05 0.35];        
        
        cf = eval(['liqCF_modis_' gcm_str])/100;
        inan = find(cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));
        cf(inan) = NaN;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_TLWP_' gcm_str './cf']);    %choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
         %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        X = X + ts; %makes the times we don't want NaN
        
        
%        Xbins = [-60:10:300]; ichoose_Xbins=1;                   
        Xbins = [0:10:500]; ichoose_Xbins=1;  
               
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
    case 'In-cloud LWP GCM norm by Model CF'
                
        switch daynight
            case 'Daytime'            
                xlabelstr = 'DAYTIME In-cloud mean LWP (norm by model CF) (g m^{-2})';
                times_required = [12:15];  %Aqua daytime
            case 'Nighttime'
                xlabelstr = 'NIGHTTIME In-cloud mean LWP (norm by model CF) (g m^{-2})';
                times_required = [0:3]; %Aqua nighttime
        end
        
        
        
         CF_gcm_thresh=0.01;
%        CF_gcm_thresh=0.8;

        cf = eval(['cf_isccp_low_' gcm_str]);
        cf(cf<CF_gcm_thresh)=NaN;

       


        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str './cf']);    
        %bit below choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        ioverride_time_selection=1;
        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        X=X+permute(time_inds_average2,[2 3 1]);
        
        
        
        Xbins = [-60:10:300]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        
    case 'In-cloud LWP GCM norm by COSP-MODIS-CF'
        
        xlabelstr = 'In-cloud mean LWP (norm by COSP-MODIS-CF) (g m^{-2})';
        
         CF_gcm_thresh=0.01;
%        CF_gcm_thresh=0.8;

        cf = eval(['liqCF_modis_' gcm_str])/100;
        cf(cf<CF_gcm_thresh)=NaN;

       


        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  
        
%        X = X(ilat,ilon,itime);
        X = eval(['1e3*gcm_lwp_' gcm_str './cf']);    
        %bit below choose only the requested lat lon points   
        X=permute(X,[2 3 1]);
        a=NaN*ones(size(X));
        a(iregion_lin)=0;
        X=X+a;
        
        Xbins = [-60:10:300]; ichoose_Xbins=1;   
        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data

        
        
     case 'Tau reduced dataset Re_1.6 Re_3.7'

        xlabelstr = 'Mean Optical Depth';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  
        
        %make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        X(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        X(re_nan) = NaN;
        
        
 case 'Tau^{1/2}'

        xlabelstr = 'Mean Optical Depth ^{1/2}';
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data
        X = (Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime)).^0.5;        


    case 'Std. dev of Nd from histogram timeseries'

        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data

        [histo_output]=...
            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,ihtot) );

        xlabelstr = 'Std dev N_d (cm^{-3}) from histogram timeseries';
        xdat(1).x=histo_output.N_histo_std;

        xlabelstr = 'Normalised std dev N_d from histogram timeseries';
        xdat(1).x=histo_output.N_std_norm;

    case 'Std. dev of Nd from histogram timeseries3'

        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        %        [histo_output]=...
        %            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );

        xlabelstr = 'Std dev N_d (cm^{-3}) from histogram timeseries';
        X = Nd_timeseries.std_dev(ilat,ilon,itime);


        xlabelstr = 'Normalised std dev N_d from histogram timeseries';
        X = Nd_timeseries.std_dev(ilat,ilon,itime)./Nd_timeseries.mean(ilat,ilon,itime);

    case 'Std. dev of Reff from histogram timeseries'

        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        %        [histo_output]=...
        %            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );

        xlabelstr = 'Std dev Reff (\mum) from timeseries';
        X = Cloud_Effective_Radius_Liquid_Standard_Deviation.timeseries3(ilat,ilon,itime);
        

    case 'Homogeneity Parameter timeseries3 using W from mean tau and Re'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        xlabelstr = 'Homogeneity Parameter using W from mean tau and Re';
        X = homog_time3_meanW(ilat,ilon,itime);

        minXbins=0;
        maxXbins=20;

    case 'Homogeneity Parameter timeseries3 using mean W for each pixel'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        xlabelstr = 'Homogeneity Parameter using mean W for each pixel';
        X = homog_time3_W(ilat,ilon,itime);        
%        minXbins=0;
%        maxXbins=20;

    case 'Homogeneity Parameter timeseries3 using Cahalan log mean W (pixel level)'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        xlabelstr = 'Homogeneity Parameter using Cahalan log mean W (pixel level)';
        X = homog_time3_logW_W(ilat,ilon,itime);

%        minXbins=0;
%        maxXbins=20;

    case 'Homogeneity Parameter Cahalan Optical Depth (Seethala)'
        ihtot = [1:prod(size_tim3)]; thresh_str='xxx'; %all the data

        xlabelstr = 'Homogeneity Parameter using Cahalan log mean tau';       
        X = homog_tau_Cahalan(ilat,ilon,itime);

%        minYbins=0;
%        maxYbins=20; 


        
case 'CDR Polder'
        
        xlabelstr = 'CDR (\mum)';

        thresh_str=[''];

        X = Par2_CDR_ALL(:,ilat,ilon,:);    %choose only the requested lat lon points 
        %N.B. - here we are using the data with teh orbit dimension kept so
        %that every data point is a separate datapoint. I.e. could have up
        %to 4 datapoints for one location on a particular day.
        %CDR Polder2 below uses the mean of the 4 orbits so that only have
        %one measurement per day
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        Xbins = [0.5:1:20.5]; ichoose_Xbins=1;
        
case 'CDR Polder2'
        
        xlabelstr = 'CDR (\mum)';

        thresh_str=[''];

        %ilat and ilon are as for a martix of size 
        X = daymeanALL_Par2_CDR(itime,ilat,ilon);    %choose only the requested lat lon points        
        X = permute(X,[2 3 1]); %change to be the same way around as MODIS (for time indices etc).
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        Xbins = [0.5:1:20.5]; ichoose_Xbins=1;    
        
    case 'CDR Polder select region'        
        xlabelstr = 'POLDER R_eff (\mum)';

        thresh_str=[''];

        %ilat and ilon are as for a martix of size 
        X = Par2_CDR_cut; %(ilat,ilon,itime);    %choose only the requested lat lon points        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        Xbins = [0.5:1:20.5]; ichoose_Xbins=1; 
        Xbins = [0.5:1:18.5]; ichoose_Xbins=1;     
        
    case 'CDR Polder Colocated with MODIS'
        %Arctic box colocation for POLDER
        %  '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_CDR_colocated_Arctic.mat'
        
        xlabelstr = 'POLDER R_{eff} (\mum)';

        thresh_str=[''];

        %ilat and ilon are as for a martix of size 
        X = Par2_CDR_coloc(ilat,ilon,itime);    %choose only the requested lat lon points        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        Xbins = [0.5:1:20.5]; ichoose_Xbins=1; 
        Xbins = [0.5:1:17.5]; ichoose_Xbins=1;         
        
        
       
        
        
    case 'Dispersion Polder'
        
        xlabelstr = '\sigma_{re}/re_{mean}';

        thresh_str=[''];

        X = Par2_CDRstd_ALL(:,ilat,ilon)./Par2_CDR_ALL(:,ilat,ilon);    %choose only the requested lat lon points           
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        Xbins = [0.0:0.05:3]; ichoose_Xbins=1;  
        
        
 case 'Std Deviation Polder'
        
        xlabelstr = '\sigma_{re}';

        thresh_str=[''];

        X = Par2_CDRstd_ALL(:,ilat,ilon);    %choose only the requested lat lon points           
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
        Xbins = [0.0:0.1:5]; ichoose_Xbins=1;          




end


%extra_title_info = [', LAT=' LAT_str ', LON=' LON_str tim_str ' ' modis_year_str2];
extra_title_info = [extra_title_info ', LAT=' LAT_str ', LON=' LON_str tim_str];


%% Y values
switch y_axis_vals
    case 'Dummy data for 1D'
        %just make Y the same as X
        Y = X;
        ylabelstr = xlabelstr;
        Ybins = Xbins;
        
    case 'General y-axis'      
%        ylabelstr = 'R_{eff 2.1 \mum} (\mum)';
        Y = Y_driver(ilat,ilon,itime);
%        Ybins = [0:0.5:30]; ichoose_Ybins=1;

    case 'General y-axis no ilat'                       
        Y = Y_driver;
        
        bin_size = 10;
        
%        Ybins = [0:0.1:maxALL(Y_driver)+0.1]; ichoose_Ybins=1;
        Ybins = [0:bin_size:maxALL(Y_driver)+0.1]; ichoose_Ybins=1;        
        
        ylabelstr = 'Precip rate (mm/hr)';
        
    case 'General y-axis no ilat simple'
        Y = Y_driver;
        

        
%        Ybins = [0:0.1:maxALL(Y_driver)+0.1]; ichoose_Ybins=1;
%        Ybins = [0:bin_size:maxALL(Y_driver)+0.1]; ichoose_Ybins=1;        
        
%        ylabelstr = '';
        
    case 'POLDER reff daymean'
        ylabelstr = 'POLDER r_e (\mum)';
        Y = daymean_Par2_CDR(ilat,ilon,itime);
        Ybins = [0:1:31]; ichoose_Ybins=1;
            

    case 'MODIS minus POLDER mock L3'
        MOD_band = '1.6\mum';        
        MOD_band = '2.1\mum';
%        MOD_band = '3.7\mum';        
        
        ylabelstr = ['MODIS ' MOD_band ' minus POLDER r_e (\mum)'];
        POL = daymean_Par2_CDR(ilat,ilon,itime);
        switch MOD_band
            case '1.6\mum'
                MOD = Cloud_Effective_Radius_16_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];
            case '2.1\mum'
                MOD = Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];
            case '3.7\mum'
                MOD = Cloud_Effective_Radius_37_Liquid_Mean_mockL3.timeseries3(ilat,ilon,itime); %xlabelstr=[xlabelstr ' minus 15% '];

        end
        
        Y = MOD - POL;
        
%        Ybins = [0:1:31]; ichoose_Ybins=1;            
    
    case 'GOES LWP'   
        %Use read_GOES_vocals_netcdf_files.m to get the goes data
        
        ylabelstr = ['LWP (g m^{-2})'];
        

         Y = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2

        
%        Y(Y<30)=NaN;

         Y(isnan(Y))=0;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        goes_LWP_save = Y;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     



        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;      
%        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1;
        
    case 'GOES Nd'   
        %Use read_GOES_vocals_netcdf_files.m to get the goes data
        
        ylabelstr = ['Nd (cm^{-3})'];
        
        Y=MODIS_justN_func(goes_Tau,goes_Reff*1e-6,'calc',0,goes_Teff,'N');
        
         goes_LWP = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2

        
        Y(goes_LWP<thresh_LWP_DRIVER)=NaN;

%         Y(isnan(Y))=0;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        goes_Nd_save = Y;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     

        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;      
%        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1;

 Ybins = [0:10:290 300:50:950 1000:500:6000]; ichoose_Ybins=1;
 Ybins = Ybins_DRIVER; ichoose_Ybins=1;     
        
   case 'GOES Nd2'   
        %Use read_GOES_vocals_netcdf_files.m to get the goes data
        
        ylabelstr = ['Nd (cm^{-3})'];
        
%        Y=MODIS_justN_func(goes_Tau,goes_Reff*1e-6,'calc',0,goes_Teff,'N');
        Y=goes_Nd_multi{79};
%         goes_LWP = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2
        goes_LWP = goes_LWP_multi{79};

        
        Y(goes_LWP<thresh_LWP_DRIVER)=NaN;

%         Y(isnan(Y))=0;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        goes_Nd_save = Y;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     

        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;      
%        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1;

 Ybins = [0:10:290 300:50:950 1000:500:6000]; ichoose_Ybins=1;
 Ybins = [0:10:290 300:50:550]; ichoose_Ybins=1;        
        
    case 'General GCM-style 4D'   
        %
        
%        ylabelstr = ['LWP (g m^{-2})')];
        
%        Y = goes_LWP;
%         Y = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2        
%        Y(Y<30)=NaN;

% Y_driver is set outside the script
        Y = Y_driver;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
%        a=NaN*ones(size(Y));
%        a(iregion_lin)=0;

% function to do the stuff below, but might use more memory?
%        a = spatial_screening_4D(Y,LAT_val,LON_val,Plat2D,Plon2D);

%         siz = size(Y);
%         Plat4D = repmat(Plat2D,[1 1 siz(3) siz(4)]);
%         Plon4D = repmat(Plon2D,[1 1 siz(3) siz(4)]);
%         iregion_lin = find(Plat4D>=LAT_val(1) & Plat4D<LAT_val(end) & Plon4D>=LON_val(1) & Plon4D<LON_val(end));
%         %iregion_lin_edges = find(Plat3D_edges>=LAT_val(1) & Plat3D_edges<LAT_val(end) & Plon3D_edges>=LON_val(1) & Plon3D_edges<LON_val(end));        
% 
%         a = NaN*ones(size(Y));
%         a(iregion_lin) = 0;
%         Y=Y+a;

% Hi-jacking this to work with regulat lat/lon grids only since was using
% to mcuch memory otherwise...
 ilat = find(Plat2D(:,1)>=LAT_val(1) & Plat2D(:,1)<LAT_val(end));
 ilon = find(Plon2D(1,:)>=LON_val(1) & Plon2D(1,:)<LON_val(end));

 Y = Y(ilat,ilon,:,:);

        
        pdf2D_Y_save = Y;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     



        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
%        Ybins = [0:10:500]; ichoose_Ybins=1;              


    case 'General GCM-style'   
        %
        
%        ylabelstr = ['LWP (g m^{-2})'];
        
%        Y = goes_LWP;
%         Y = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2        
%        Y(Y<30)=NaN;

% Y_driver is set outside the script
        Y = Y_driver;
% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        
        if length(size(Plat3D))==3
            Y=permute(Y,[2 3 1]);
        end              
    
        times_required = [0:24]; %



        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;  %iregion_lin based on Plat3D
        Y=Y+a;
        
        pdf2D_Y_save = Y;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     



        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
%        Ybins = [0:10:500]; ichoose_Ybins=1;              
        
    case 'UM LWP'
        %Use read_GOES_vocals_netcdf_files.m to get the goes data
        
        ylabelstr = ['LWP (g m^{-2})'];
        
%        Y = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2
 
        Y = lwp_UM_n5;
        
        
%        Y(Y<30)=NaN;
        
              
    
        times_required = [0:24]; %

% Regional screening for data with 2D lat lon grids (non-regular)   
        %Change the order to [lat lon time] to fit with Plat3D
        % May not need to do this!!
%        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
%         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
% Time screening taking into account local time variation with longitude
%         ioverride_time_selection=1;         
%         %do the time screening again with the override
%         time_inds_modisL3_timeseries3
%         %calculates time_inds_average2, which is NaN at times we don't want
%         Y=Y+permute(time_inds_average2,[2 3 1]);
%     



        extra_title_info = [extra_title_info ' ' gcm_str];        
%        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;      
%        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1;
%        Ybins = [-0.01 10.^[log10(10):0.1:log10(2500)]]; ichoose_Ybins=1;        
        

    case 'Precip rate Comstock using MODIS'                
        ylabelstr = 'Precipitation Rate using MODIS LWP (mm hr^{-1})';
        
        cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);        

%MOD06
%cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

        LWP_MODIS_PDF = (1/1.15)*1e3*W_time3(ilat,ilon,itime); ylabelstr=[ylabelstr ' minus 15% '];        
%        LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1)./cf;
        
        N  = N_time3(ilat,ilon,itime);
        
        Y = precip_rate_Wood_2008(LWP_MODIS_PDF,N);
            %function P = precip_rate_Wood_2008(W,N)
            % Returns precipitation rate in mm/hour given the LWP (W) in g/m2
            % and N in cm^-3 using the formula
            
            Ybins = [0:0.001:1]; ichoose_Ybins=1;
            Ybins = [0:0.01:1]; ichoose_Ybins=1;
            
     case 'Precip rate Comstock using AMSRE'                
        ylabelstr = 'Precip Rate using AMSRE LWP (scaled w/ MOD35) (mm hr^{-1})';
        
        cf = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);        

%MOD06
%cf = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);

%        LWP_MODIS_PDF = (1/1.15)*1e3*W_time3(ilat,ilon,itime); ylabelstr=[ylabelstr ' MODIS LWP minus 15% '];        
        LWP_AMSRE_PDF = 1e3*lwp_amsre_time3(ilat,ilon,itime,1)./cf;
        
        N  = N_time3(ilat,ilon,itime);
        
        Y = precip_rate_Wood_2008(LWP_AMSRE_PDF,N);
            %function P = precip_rate_Wood_2008(W,N)
            % Returns precipitation rate in mm/hour given the LWP (W) in g/m2
            % and N in cm^-3 using the formula
            
            Ybins = [0:0.001:1]; ichoose_Ybins=1;
            Ybins = [0:0.01:1]; ichoose_Ybins=1;           
        

 case 'TOA albedo'
        CF_gcm_thresh=[0.35 0.45];
%        CF_gcm_thresh=[-0.01 1.01];
        
        LWP_gcm_thresh=[-0.01 2]; % g/m2
%        LWP_gcm_thresh=[10 20]; % g/m2        
%        LWP_gcm_thresh=[-0.01 5e3]; % g/m2        
        
        ylabelstr = ['TOA Albedo for LWP.GTE.' num2str(LWP_gcm_thresh(1)) '.AND.LT.' num2str(LWP_gcm_thresh(2)) 'and for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2))];
        
        LWP = eval(['1e3*gcm_lwp_' gcm_str]); %g/m2
        LWP=permute(LWP,[2 3 1]);
        CF = eval(['cf_isccp_low_' gcm_str]);
        CF=permute(CF,[2 3 1]);
        irem = find(LWP<=LWP_gcm_thresh(1) | LWP>LWP_gcm_thresh(2) | CF<=CF_gcm_thresh(1) | CF>CF_gcm_thresh(2));
            
        Y = eval(['albedo_' gcm_str]);
        Y=permute(Y,[2 3 1]);
        Y(irem) = NaN; %remove the points we don't want to include for LWP thresholding
        
%        SW_down =  eval(['SW_TOA_net_' gcm_str '+SW_TOA_up_' gcm_str]);
%        SW_down =  permute(SW_down,[2 3 1]);
%        Y(SW_down<10) = NaN;
        
        a=NaN*ones(size(Y));                
        a(iregion_lin)=0;
        Y=Y+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        Y = Y + ts; %makes the times we don't want NaN
        
        
%        Y(Y<CF_gcm_thresh)=NaN;



%Xbins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Xbins=1; %
Ybins = [0:0.02:1]; ichoose_Ybins=1; %

% ipost_plotTime_commands = 1;
% post_plotTime_commands = {'set(gca,''ylim'',[0 0.6])'};

extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];
        
    case 'SZA Polder Colocated with MODIS'
        %Arctic box colocation for POLDER
        %  '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_CDR_colocated_Arctic.mat'
        
        ylabelstr = 'POLDER SZA (degrees)';

        %ilat and ilon are as for a martix of size 
        Y = Par2_sza_coloc(ilat,ilon,itime);    %choose only the requested lat lon points  
        
    case 'SZA Polder select region'
        ylabelstr = 'POLDER SZA (degrees)';

        %ilat and ilon are as for a martix of size 
        Y = Par2_sza_cut; %(ilat,ilon,itime);    %choose only the requested lat lon points 
        
       
        
    case 'AMSRE TLWP DAYTIME'
        ylabelstr = 'Daytime AMSRE TLWP (g m^{-2})';        
        Y = 1e3*lwp_amsre_time3(ilat,ilon,itime,1); %the one at the end is for the ascending node (daytime)   
        
        if ireduce_res_lwp==1
            %Reduce the resolution by averaging over boxes of size M_lwp x
            %N_lwp
            Y1 = reduce_matrix_subsample_mean(Y(:,:,1),M_lwp,N_lwp);
            siz_Y1 = size(Y1);
            siz_tim = size(Y,3);
            Y2 = NaN*ones([siz_Y1 siz_tim]);
            for ired=1:siz_tim
                Y2(:,:,ired) = reduce_matrix_subsample_mean(Y(:,:,ired),M_lwp,N_lwp);
            end
            Y=Y2;
        end        

        if ~exist('ioverride_Ybins') | ioverride_Ybins==0
            Ybins = [-60:10:2000]; ichoose_Ybins=1;
        end

    case 'AMSRE TLWP NIGHTTIME'
        ylabelstr = 'Nighttime AMSRE TLWP (g m^{-2})';        
        Y = 1e3*lwp_amsre_time3(ilat,ilon,itime,2); %the one at the end is for the ascending node (daytime)                          
        Ybins = [-60:10:2000]; ichoose_Ybins=1;  
        
        if ireduce_res_lwp==1
            %Reduce the resolution by averaging over boxes of size M_lwp x
            %N_lwp
            Y1 = reduce_matrix_subsample_mean(Y(:,:,1),M_lwp,N_lwp);
            siz_Y1 = size(Y1);
            siz_tim = size(Y,3);
            Y2 = NaN*ones([siz_Y1 siz_tim]);
            for ired=1:siz_tim
                Y2(:,:,ired) = reduce_matrix_subsample_mean(Y(:,:,ired),M_lwp,N_lwp);
            end
            Y=Y2;
        end 
        
    case 'Cloudsat precip DAYTIME'
        ylabelstr = 'Daytime Cloudsat precip rate (mm hr^{-1})';      
        
%Haven't set up the correct time selection options yet...        
         dat_modis = rain_warm_ocean(:,ilat,ilon); %mm/hr
        
                                    
                    Y = dat_modis(1:2:end,:,:);
%                    time_inds_average = time_inds_average(1:2:end);

            
                    
%        Y = 1e3*lwp_amsre_time3(ilat,ilon,itime,1); %the one at the end is for the ascending node (daytime)     
        
        
        Ybins = [-60:10:2000]; ichoose_Ybins=0;  

    case 'Cloudsat precip NIGHTTIME'
        ylabelstr = 'Nighttime AMSRE TLWP (g m^{-2})';   
        dat_modis = rain_warm_ocean(:,ilat,ilon);
        Y = dat_modis(2:2:end,:,:);
%                    time_inds_average = time_inds_average(2:2:end);        
 
%        Y = 1e3*lwp_amsre_time3(ilat,ilon,itime,2); %the one at the end is for the ascending node (daytime)                          
        Ybins = [-60:10:2000]; ichoose_Ybins=0;         

  case 'TLWP DAYTIME GCM, grid-box average'        
        ylabelstr = 'Daytime grid-box average LWP+RWP (g m^{-2})';
        
%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  thresh_str='';
        
        Y = eval(['1e3*gcm_TLWP_' gcm_str ';']);    %choose only the requested lat lon points   
        Y=permute(Y,[2 3 1]);
        a=NaN*ones(size(Y));
        %iregion_lin is based on Plat3D which is of size [lat lon time] -
        %this is why we permute Y above
        a(iregion_lin)=0;
        Y=Y+a; 
        

        
        Ybins = [0:10:2000]; ichoose_Ybins=1;   
        
        ioverride_time_selection=1;
         times_required = [12:15];  %Aqua daytime
%                        times_required = [0:3]; %Aqua nighttime

        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);

        
       



 case 'TLWP NIGHTTIME GCM, grid-box average'        
     ylabelstr = 'Nighttime grid-box average LWP+RWP (g m^{-2})';
        
%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  thresh_str='';
        
        Y = eval(['1e3*gcm_TLWP_' gcm_str ';']);    %choose only the requested lat lon points   
        Y=permute(Y,[2 3 1]);
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        Ybins = [0:10:2000]; ichoose_Ybins=1;   
        
        ioverride_time_selection=1;
%         times_required = [12:15];  %Aqua daytime
         times_required = [0:3]; %Aqua nighttime
        
        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        

        
       
        
        
  case 'Precip rate DAYTIME GCM'        
        ylabelstr = 'Daytime precip rate (mm day^{-1})';
        
%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
        thresh_str='';
  
        precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s

        max_bin_val = 0.4*24;

        Y = precip*3600;  %mm/hr
        

        Y=permute(Y,[2 3 1]);
        a=NaN*ones(size(Y));
        %iregion_lin is based on Plat3D which is of size [lat lon time] -
        %this is why we permute Y above
        a(iregion_lin)=0;
        Y=Y+a;        
        
%        Ybins = [0:10:2000]; ichoose_Ybins=1;   
        
        ioverride_time_selection=1;
         times_required = [12:15];  %Aqua daytime
%                        times_required = [0:3]; %Aqua nighttime

        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        
%        dbin=0.002*24;
%        Ybins = [-0.00001:dbin:max_bin_val-0.00001+dbin]; ichoose_Ybins=1;
        
  case 'Precip rate NIGHTTIME GCM'        
        ylabelstr = 'Nighttime precip rate (mm day^{-1})';
        
%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
        thresh_str='';
  
        precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s

        max_bin_val = 0.4*24;

        Y = precip*3600;  %mm/hr
        

        Y=permute(Y,[2 3 1]);
        a=NaN*ones(size(Y));
        %iregion_lin is based on Plat3D which is of size [lat lon time] -
        %this is why we permute Y above
        a(iregion_lin)=0;
        Y=Y+a;        
        
%        Ybins = [0:10:2000]; ichoose_Ybins=1;   
        
        ioverride_time_selection=1;
         %times_required = [12:15];  %Aqua daytime
         times_required = [0:3]; %Aqua nighttime

        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        
%        dbin=0.002*24;
%        Ybins = [-0.00001:dbin:max_bin_val-0.00001+dbin]; ichoose_Ybins=1;


  case 'LWP DAYTIME GCM, grid-box average'        
        ylabelstr = 'Daytime grid-box average LWP (g m^{-2})';
        
%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  thresh_str='';
        
        Y = eval(['1e3*gcm_lwp_' gcm_str ';']);    %choose only the requested lat lon points   
        Y=permute(Y,[2 3 1]);
        a=NaN*ones(size(Y));
        %iregion_lin is based on Plat3D which is of size [lat lon time] -
        %this is why we permute Y above
        a(iregion_lin)=0;
        Y=Y+a;        
        
        Ybins = [0:10:2000]; ichoose_Ybins=1;   
        
        ioverride_time_selection=1;
         times_required = [12:15];  %Aqua daytime
%                        times_required = [0:3]; %Aqua nighttime

        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        
        
        


 case 'LWP NIGHTTIME GCM, grid-box average'        
     ylabelstr = 'Nighttime grid-box average LWP (g m^{-2})';
        
%        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
  thresh_str='';
        
        Y = eval(['1e3*gcm_lwp_' gcm_str ';']);    %choose only the requested lat lon points   
        Y=permute(Y,[2 3 1]);
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        Ybins = [0:10:2000]; ichoose_Ybins=1;   
        
        ioverride_time_selection=1;
%         times_required = [12:15];  %Aqua daytime
         times_required = [0:3]; %Aqua nighttime
        
        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        
       

    case 'MODIS LWP, grid-box average using MOD06';
        ylabelstr = 'MODIS LWP gridbox average using MOD06 (g m^{-2})';        
        Y = 1e3*W_time3(ilat,ilon,itime,1).*Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime,1); %the one at the end is for the ascending node (daytime)                          
        Ybins = [-60:10:2000]; ichoose_Ybins=1; 
        
    case 'MODIS LWP, grid-box average using MOD35';
        ylabelstr = 'MODIS LWP gridbox average using MOD35 (g m^{-2})';        
        Y = 1e3*W_time3(ilat,ilon,itime,1).*Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime,1); %the one at the end is for the ascending node (daytime)                          
        Ybins = [-60:10:2000]; ichoose_Ybins=1;     
        
    case 'MODIS LWP, grid-box average using MOD35 minus 15%';
        ylabelstr = 'MODIS LWP gridbox average using MOD35 minus 15% (g m^{-2})';        
        Y = (1/1.15).*1e3*W_time3(ilat,ilon,itime,1).*Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime,1); %the one at the end is for the ascending node (daytime)                          
        Ybins = [-60:10:2000]; ichoose_Ybins=1;         
        
        
    case 'Re COSP GCM'        
        cf = eval(['liqCF_modis_' gcm_str])/100;        
%        tau = eval(['liqTau_modis_' gcm_str './cf']);
        re = eval(['liqRe_modis_' gcm_str './cf']);
        
        ylabelstr = 'COSP Reff';

        CF_gcm_thresh=0.01;
        CF_gcm_thresh=0.8;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
        extra_title_info = [extra_title_info ' ' thresh_str ' '];

        re(cf<CF_gcm_thresh)=NaN;
        re=permute(re,[2 3 1]);
               
        Y = 1e6*re;    
        %choose only the requested lat lon points   
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;  %so if a is NaN then Y will become NaN
        
        Ybins = [0:1:30]; ichoose_Ybins=1;                

        
    case 'R_{eff 2.1 \mum} (\mum)'        
        ylabelstr = 'R_{eff 2.1 \mum} (\mum)';
        Y = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Ybins = [0:0.5:30]; ichoose_Ybins=1;
        
case 'R_{eff 2.1 \mum} (\mum) minus 20%'        
        ylabelstr = 'R_{eff 2.1 \mum} (\mum) minus 20%';
        Y = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Y = Y*0.8;
        Ybins = [0:0.5:30]; ichoose_Ybins=1;        
        
case 'CDR Polder2'
        
        ylabelstr = 'CDR (\mum)';
        %ilat and ilon are as for a matrix of size 
        Y = daymeanALL_Par2_CDR(itime,ilat,ilon);    %choose only the requested lat lon points   
        Y = permute(Y,[2 3 1]); %change to be the same way around as MODIS (for time indices etc).
        % Note - here am not choosing the time dimension
        Ybins = [0.5:1:30.5]; ichoose_Ybins=1;        

    case 'Nd GCM'
%        cf = eval(['liqCF_modis_' gcm_str])/100;        
%        tau = eval(['liqTau_modis_' gcm_str './cf']);
        Nd = eval(['gcm_Nd_max_screen_' gcm_str]);
        
        ylabelstr = 'Nd (cm^{-3})';

        CF_gcm_thresh=0.01;
        CF_gcm_thresh=0.8;

        thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];

 %       Nd(cf<CF_gcm_thresh)=NaN;
        Y=permute(Nd,[2 3 1]);
               
        %choose only the requested lat lon points   
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;  %so if a is NaN then Y will become NaN
        
        Ybins = [0:10:1000]; ichoose_Ybins=1;                


        
        
      case 'Time of Day'
        ylabelstr = 'Local Time';
        %the remainder from the previous whole day is the fraction of days
        %from 0 UTC
        modis_UTC = ( Date_Time_Swath.timeseries3 - floor(Date_Time_Swath.timeseries3) )*24;
        
        lon = squeeze(repmat(LON(ilon),[length(ilat) 1 size_tim3(3)]));
        Y = modis_UTC(ilat,ilon,itime) + lon*24/360;   
        Y = mod(Y,24);
        
         Ybins = [0:24]; ichoose_Ybins=1;
         Ybins = [0.5:1:23.5]; ichoose_Ybins=1;         
        
    case 'Nd lat-lon cells'
        ylabelstr = 'Nd (cm^{-3})';
        Y = meanNoNan(Nd_all,3);
    case 'Latitude'
        ylabelstr = 'Latitude';
        
        switch datatype
            case 'gcm_data'
                Y = Plat3D;  %Note that Plon3D doesn't need to be permuted
                a=NaN*ones(size(Plat3D)); a(iregion_lin)=0;
                Y = Y + a;
                       
                %also add the time screening
                ts = permute(time_inds_average2,[2 3 1]);
                Y = Y + ts; %makes the times we don't want NaN
        
        
                %MODIS lat and lons are the bin centres, so choose bins to be the edges between them
                %Trying to make these consistent between plots - make them
                %odd numbers
                Ybins = [1+2*floor(minALL(Y)/2):5:-1+2*ceil(maxALL(Y)/2)]; ichoose_Ybins=1;
                
                lat_bin = 5;
                Ybins = [lat_bin*floor(minALL(Y)/lat_bin):5:lat_bin*ceil(maxALL(Y)/lat_bin)]; ichoose_Ybins=1;                
                
            otherwise
                
                Y = squeeze(repmat(LAT(ilat),[1 1 length(ilon) size_tim3(3)]));
        
        end
        
  case 'Longitude'
        ylabelstr = 'Longitude';
        
        switch datatype
            case 'gcm_data'
                Y = Plat3D;  %Note that Plon3D doesn't need to be permuted
                a=NaN*ones(size(Plon3D)); a(iregion_lin)=0;
                Y = Y + a;
                       
                %also add the time screening
                ts = permute(time_inds_average2,[2 3 1]);
                Y = Y + ts; %makes the times we don't want NaN
        
        
                %MODIS lat and lons are the bin centres, so choose bins to be the edges between them
                %Trying to make these consistent between plots - make them
                %odd numbers
                Ybins = [1+2*floor(minALL(Y)/2):5:-1+2*ceil(maxALL(Y)/2)]; ichoose_Ybins=1;
                
              
                            
                
            otherwise                
                Y = squeeze(repmat(LON(ilon),[length(ilat) 1 length(itime)]));
        end
        
         lat_bin = 5;
         Ybins = [lat_bin*floor(minALL(Y)/lat_bin):5:lat_bin*ceil(maxALL(Y)/lat_bin)]; ichoose_Ybins=1;   
        
        
    case 'Min Scattering Angle timeseries3'
        ylabelstr = 'Min Scattering Angle';
        Y = Scattering_Angle_Maximum.timeseries3(ilat,ilon,itime);
    case 'Max Scattering Angle timeseries3'
        ylabelstr = 'Max Scattering Angle';
        Y = Scattering_Angle_Minimum.timeseries3(ilat,ilon,itime);

    case 'Mean Scattering Angle timeseries3'
        ylabelstr = 'Mean Scattering Angle';
        Y = Scattering_Angle_Mean.timeseries3(ilat,ilon,itime);

    case 'SZA min timeseries3'
        ylabelstr = 'Minimum SZA';
        Y = Solar_Zenith_Minimum.timeseries3(ilat,ilon,itime);

    case 'SZA max timeseries3'
        ylabelstr = 'Maximum SZA';
        Y = Solar_Zenith_Maximum.timeseries3(ilat,ilon,itime);

    case 'SZA min max difference'
        ylabelstr = 'Min Max SZA difference';
        Y = Solar_Zenith_Maximum.timeseries3(ilat,ilon,itime) - Solar_Zenith_Minimum.timeseries3(ilat,ilon,itime);


    case 'WMOD'
        ylabelstr = 'W times five sixths (kg m^{-2})';
        Y = WMOD(:); %MOD06

    case 'CF'
        ylabelstr = 'Cloud Fraction';
        Y = cf(:); %MOD06 CF
        
        
    case 'LWP GCM Grid-box mean pre-mphys'
          var_gcm_thresh=[1 1e9]; %screening set up for pre mphys LWP
          var_gcm_thresh=[-0.01 1e9]; %screening set up for pre mphys LWP          
%          var_gcm_thresh=[-0.01 1]; %g m-2       

          precip_gcm_thresh = [0.1 1e9]; %mm day^{-1}
          precip_gcm_thresh = [-0.01 1e9]; %mm day^{-1}

                     
        
        switch daynight
            case 'Daytime'            
                ylabelstr = 'Pre Mphys LWP (g m^{-2}), DAYTIME';
                times_required = [12:15];  %Aqua daytime
            case 'Nighttime'
                ylabelstr = 'Pre Mphys LWP (g m^{-2}), NIGHTTIME';
                times_required = [0:3]; %Aqua nighttime
            case 'ALL times of day'
                ylabelstr = 'Pre Mphys LWP (g m^{-2}), ALL times of day';
                times_required = [0:24]; %
        end
        
        LWP_postMP = eval(['lwpAP_isccp_low_' gcm_str]) + eval(['lwpAP_isccp_mid_' gcm_str]) + eval(['lwpAP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points 
        LWP_preMP = eval(['lwpBP_isccp_low_' gcm_str]) + eval(['lwpBP_isccp_mid_' gcm_str]) + eval(['lwpBP_isccp_high_' gcm_str]);  %kg/m2  %choose only the requested lat lon points 
  
        LWP_postMP = LWP_postMP *1e3;  %convert to g/m2      
        LWP_preMP = LWP_preMP *1e3;
        dLWP = LWP_postMP - LWP_preMP; %Change due to mphys (g/m2)
        
%        Y = eval(['1e3*gcm_lwp_' gcm_str]);
         Y = LWP_preMP;       %g/m2
        
        var = LWP_preMP;
        inan = find(var<var_gcm_thresh(1) | var>=var_gcm_thresh(2));
        Y(inan) = NaN;
        
        
%        var = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s
%        var = var *3600*24;  %mm/day
%        inan = find(var<precip_gcm_thresh(1) | var>=precip_gcm_thresh(2));
%        Y(inan) = NaN;
        
        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        
%        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
         ylabelstr = [ylabelstr ' for pre_mphysLWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2))];
         
         ioverride_time_selection=1;         
        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        
        
        
%        Ybins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Ybins=1; %
        Ybins = [0:1:1000]; ichoose_Ybins=1; %        



        extra_title_info = [extra_title_info ' ' gcm_str];
        
        ipost_plotTime_commands = 1; post_plotTime_commands = {'set(gca,''yscale'',''log'');','set(gca,''ylim'',[0.9 300])'};

        
        
        
    case 'CF GCM'
          var_gcm_thresh=[-0.01 1e9]; %g m-2        
%          var_gcm_thresh=[-0.01 1]; %g m-2       

          precip_gcm_thresh = [0.1 1e9]; %mm day^{-1}
          precip_gcm_thresh = [-0.01 1e9]; %mm day^{-1}

                     
        
        switch daynight
            case 'Daytime'            
                ylabelstr = 'Max low cloud fraction (no screening) DAYTIME';
                times_required = [12:15];  %Aqua daytime
            case 'Nighttime'
                ylabelstr = 'Max low cloud fraction (no screening) NIGHTTIME';
                times_required = [0:3]; %Aqua nighttime
            case 'ALL times of day'
                ylabelstr = 'Max low cloud fraction (no screening) ALL times of day';
                times_required = [0:24]; %
        end
        
        Y = eval(['cf_isccp_low_' gcm_str]);
       
        
        var = eval(['1e3*gcm_lwp_' gcm_str]);
        inan = find(var<var_gcm_thresh(1) | var>=var_gcm_thresh(2));
        Y(inan) = NaN;
        
        
        var = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s
        var = var *3600*24;  %mm/day
        inan = find(var<precip_gcm_thresh(1) | var>=precip_gcm_thresh(2));
        Y(inan) = NaN;
        
        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));
        a(iregion_lin)=0;
        Y=Y+a;
        
        
        ylabelstr = [ylabelstr ' for LWP.GTE.' num2str(var_gcm_thresh(1)) '.AND.LT.' num2str(var_gcm_thresh(2)) ' AND.precip.GTE.' num2str(precip_gcm_thresh(1)) '.AND.LT.' num2str(precip_gcm_thresh(2)) ' mm day^{-1}'];
        
         ioverride_time_selection=1;         
        %do the time screening again with the override
        time_inds_modisL3_timeseries3
        %calculates time_inds_average2, which is NaN at times we don't want
        Y=Y+permute(time_inds_average2,[2 3 1]);
        
        

        Ybins = [-0.05:0.1:1.05]; ichoose_Ybins=1;
        Ybins = [-.0001:0.1:0.8999 1.0001]; ichoose_Ybins=1;
        Ybins = [-0.1 0.001:0.1:0.9001 0.999 1.1]; ichoose_Ybins=1;
        Ybins = [0:0.05:0.95 1.001]; ichoose_Ybins=1; 
        Ybins = [-0.05:0.1:1.05]; ichoose_Ybins=1; %same as for model COSP bins
        
        Ybins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Ybins=1; %this makes sure that the bin widths actually represent
%the ranges possible. This makes sense for MODIS. But for COSP I'm not sure
%what the singlar values actually mean - does 0.1 mean 0.05 to 0.15? Or
%exactly 0.1??
% so here we have 0 <= CF < 0.05, 0.05 <= CF < 0.1, ....  0.95 <= CF <= 1.0



        extra_title_info = [extra_title_info ' ' gcm_str];


        
 case 'CF COSP-MODIS GCM'
        ylabelstr = 'COSP-MODIS GCM Cloud Fraction';
        CF_gcm_thresh=[-0.01 1.01];
        %CF_gcm_thresh=0.8;
        
        
            
        Y = eval(['liqCF_modis_' gcm_str])/100;
        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));                
        a(iregion_lin)=0;
        Y=Y+a;
        
        
%        Y(Y<CF_gcm_thresh)=NaN;

%COSP CFs are multiples of 0.1 (0, 0.1, 0.2 ... 1.0)
Ybins = [-0.05:0.1:1.05]; ichoose_Ybins=1;
Ybins = [-0.0999:0.1:0.9999 1.0001]; ichoose_Ybins=1; %[-0.0999    0.0001    0.1001    0.2001    0.3001    0.4001    0.5001    0.6001    0.7001    0.8001    0.9001    1.0001]
%this will ensure that the each CF value falls into one bin only (0, 0.1,
%0.2, etc.)
% so here we have 0 <= CF < 0.05, 0.05 <= CF < 0.1, ....  0.95 <= CF <= 1.0
Ybins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Ybins=1;
  
extra_title_info = [extra_title_info ' ' gcm_str];
        

case 'CF COSP-CALIPSO GCM'
        ylabelstr = 'COSP-CALIPSO GCM Cloud Fraction';
        CF_gcm_thresh=0.01;
        %CF_gcm_thresh=0.8;
        
        
            
        Y = eval(['cllcalipso_' gcm_str '/100']);
        Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));                
        a(iregion_lin)=0;
        Y=Y+a;
        
        
%        Y(Y<CF_gcm_thresh)=NaN;

%COSP CFs are multiples of 0.1 (0, 0.1, 0.2 ... 1.0)
Ybins = [-0.05:0.1:1.05]; ichoose_Ybins=1;
%Ybins = [-0.0999:0.1:0.9999 1.0001]; ichoose_Ybins=1; %[-0.0999    0.0001    0.1001    0.2001    0.3001    0.4001    0.5001    0.6001    0.7001    0.8001    0.9001    1.0001]
%this will ensure that the each CF value falls into one bin only (0, 0.1,
%0.2, etc.)

Ybins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Ybins=1; %this makes sure that the bin widths actually represent
%the ranges possible. This makes sense for MODIS. But for COSP I'm not sure
%what the singlar values actually mean - does 0.1 mean 0.05 to 0.15? Or
%exactly 0.1??
% so here we have 0 <= CF < 0.05, 0.05 <= CF < 0.1, ....  0.95 <= CF <= 1.0

extra_title_info = [extra_title_info ' ' gcm_str];


    case 'LWP removal timescale'
        ylabelstr = 'LWP removal timescale (hrs)';
        CF_gcm_thresh=[-0.01 1.01];
        %CF_gcm_thresh=0.8;
        
      %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
               LWP = eval(['gcm_lwp_' gcm_str]);  %kg/m2
            precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s
            %only calcuate when precip greater than a threshold - otherwise
            %the zero precip times will swamp the answer
            %Perhaps better plotted vs precip rate? Or a PDF might be
            %better
            precip_thresh = 1e-10; max_timescale = 4*1;
            iprecip = find(precip<precip_thresh);
            precip(iprecip)=precip_thresh;
            
            Y = LWP./precip/3600; 
            Y(Y>max_timescale)=max_timescale;
            Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));                
        a(iregion_lin)=0;
        Y=Y+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        Y = Y + ts; %makes the times we don't want NaN
        
        
 
dhrs=0.05;
Ybins = [-0.00001:dhrs:max_timescale-0.00001+dhrs]; ichoose_Ybins=1; 

extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];


    case 'Precip rate'

        CF_gcm_thresh=[-0.01 1.01];
%        CF_gcm_thresh=[0.05 1.01];        
    
   
%       cf = eval(['liqCF_modis_' gcm_str])/100;
       cf = eval(['cf_isccp_low_' gcm_str]);  
      
       
        lwp_gcm_thresh=[-0.01 1e9]; %g m-2        
%        lwp_gcm_thresh=[-0.01 1]; %g m-2   
%        lwp_gcm_thresh=[1 1e9]; %g m-2           

%what to use for the LWP screening
LWP_or_TLWP = 'LWP'; 
LWP_or_TLWP = 'TLWP'; %Or LWP+RWP

switch LWP_or_TLWP
    case 'LWP'
        lwp = eval(['1e3*gcm_lwp_' gcm_str]);
    case 'TLWP'
        lwp = eval(['1e3*gcm_TLWP_' gcm_str]);
end

inan = find(lwp<lwp_gcm_thresh(1) | lwp>=lwp_gcm_thresh(2) | cf<CF_gcm_thresh(1) | cf>=CF_gcm_thresh(2));


ylabelstr = ['Precipitation rate (mm day^{-1}) for CF.GTE.' num2str(CF_gcm_thresh(1)) '.AND.LT.' num2str(CF_gcm_thresh(2)) ', ' LWP_or_TLWP '.GTE.' num2str(lwp_gcm_thresh(1)) '.AND.LT.' num2str(lwp_gcm_thresh(2))];
        
      %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
%               LWP = eval(['gcm_lwp_' gcm_str]);  %kg/m2
            precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s = mm/s
            %only calcuate when precip greater than a threshold - otherwise
            %the zero precip times will swamp the answer
            %Perhaps better plotted vs precip rate? Or a PDF might be
            %better
%            precip_thresh = 1e-10; max_timescale = 4*1;
%            iprecip = find(precip<precip_thresh);
%            precip(iprecip)=precip_thresh;

precip(inan) = NaN;

max_bin_val = 0.4*24;

            Y = precip*3600*24;  %mm/day
            Y(Y>max_bin_val)=max_bin_val;
            Y=permute(Y,[2 3 1]);
        
        a=NaN*ones(size(Y));                
        a(iregion_lin)=0;
        Y=Y+a;
        
        %also add the time screening
        ts = permute(time_inds_average2,[2 3 1]);
        Y = Y + ts; %makes the times we don't want NaN
        
        
 
dbin=0.002*24;
Ybins = [-0.00001:dbin:max_bin_val-0.00001+dbin]; ichoose_Ybins=1; 

extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];



        
    case 'CF all'
        ylabelstr = 'Cloud Fraction (all phases)';
        Y = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime)+Cloud_Fraction_Ice.timeseries3(ilat,ilon,itime)+Cloud_Fraction_Undetermined.timeseries3(ilat,ilon,itime);
%        Cloud_Fraction_Combined.timeseries(ilat,ilon,itime) %think
%        Combined is the same
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];    
        
    case 'CF ice'
        ylabelstr = 'Ice Cloud Fraction';
        Y = Cloud_Fraction_Ice.timeseries3(ilat,ilon,itime);
%        Cloud_Fraction_Combined.timeseries(ilat,ilon,itime) %think
%        Combined is the same
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];    
        
     case 'CF undet'
        ylabelstr = 'Undetermined Cloud Fraction';
        Y = Cloud_Fraction_Undetermined.timeseries3(ilat,ilon,itime);
%        Cloud_Fraction_Combined.timeseries(ilat,ilon,itime) %think
%        Combined is the same
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];        

    case 'Scattering angle'
        ylabelstr = 'Scattering angle (degrees)';
        Y = sangle(:);

    case 'SZA std dev timeseries'
        ylabelstr = 'Std dev of SZA';
        Y = Solar_Zenith_Standard_Deviation.timeseries(ilat,:);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Mean SZA timeseries'
        ylabelstr = 'Mean SZA';
        Y = Solar_Zenith_Mean.timeseries(ilat,:);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Mean SZA timeseries3'
        ylabelstr = 'Solar Zenith Angle';
        Y = Solar_Zenith_Mean.timeseries3(ilat,ilon,itime);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        
    case 'Mean CTT timeseries3, y-axis'
        ylabelstr = 'Mean Cloud Top Temperature (K)';
        Y = Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=200;
        maxYbins=310;
        
    case 'Min CTT timeseries3, y-axis'
        ylabelstr = 'Minimum Cloud Top Temperature (K)';
        Y = Cloud_Top_Temperature_Day_Minimum.timeseries3(ilat,ilon,itime);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

%        minYbins=200;
%        maxYbins=310;

 case 'Mean CTH_max_CTP_min timeseries3 flag zeroCF, y-axis'
          %When the cloud fraction is very low (but not NaN) the CTH is set to a
     %value to flag it as such (-1e9) - since if the CF was zero there could
     %(presumably? Does it get set to NaN? Yes, seems to.) be no CTT value.
        ylabelstr = 'Mean Cloud Top Height (km)';
        
        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
%        subtract_val = -2.5; %Works well using mean CTT from L3 for the region     LAT_val = [-40.5 -30.5]; LON_val = [-140 -100]; 
        
        Y = CTH_hybrid_std_atmos_CTP.timeseries3(ilat,ilon,itime);
        
        cf_min_thresh = 0.01; %0.01;
        ilow_cf = find( Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime) < cf_min_thresh);
        Y(ilow_cf) = Ybins(1)+1e-30; %set this as a marker that is not NaN, but make it within range of Ybins - will then discard
        %the first bin
        % It seems likely that the inclusion of points where the CF is very
        % small (and so no CTH reading) is not what is important here - it
        % is more the inclusion of the zeros at all heights above teh
        % observed CTH that makes a difference, I think.
        
        %Need to somehow add an option to have a CTH if there is no cloud
        %(CF=0) - did something for this in the screening routine.
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20; 
        
        %CF profiles for each location are created in case 961 of watervap.
        %(run using case 131 with  multi_plot_case_122 = 'CF_MOD35 vs CTH
        %MOD35').
        % This sets the CF for each profile to one at the CTH and zero
        % above this (with NaNs below). FOr cases where we have no CTH we
        % assume the CF is zero for all heights.
      
        


     

 case 'Mean CTH_max timeseries3 flag zeroCF, y-axis'
     %When the cloud fraction is very low (but not NaN) the CTH is set to a
     %value to flag it as such (-1e9) - since if the CF was zero there could
     %(presumably? Does it get set to NaN? Yes, seems to.) be no CTT value.
        ylabelstr = 'Mean Cloud Top Height (km)';
        
        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
%        subtract_val = -2.5; %Works well using mean CTT from L3 for the region     LAT_val = [-40.5 -30.5]; LON_val = [-140 -100]; 
        
        Y = CTH_max.timeseries3(ilat,ilon,itime);
        
        cf_min_thresh = 0.01; %0.01;
        ilow_cf = find( Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime) < cf_min_thresh);
        Y(ilow_cf) = Ybins(1)+1e-30; %set this as a marker that is not NaN, but make it within range of Ybins - will then discard
        %the first bin
        % It seems likely that the inclusion of points where the CF is very
        % small (and so no CTH reading) is not what is important here - it
        % is more the inclusion of the zeros at all heights above teh
        % observed CTH that makes a difference, I think.
        
        %Need to somehow add an option to have a CTH if there is no cloud
        %(CF=0) - did something for this in the screening routine.
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20; 
        
        %CF profiles for each location are created in case 961 of watervap.
        %(run using case 131 with  multi_plot_case_122 = 'CF_MOD35 vs CTH
        %MOD35').
        % This sets the CF for each profile to one at the CTH and zero
        % above this (with NaNs below). FOr cases where we have no CTH we
        % assume the CF is zero for all heights.
      
        

        
        
 case 'Mean CTH_Ryan timeseries3 flag zeroCF, y-axis'
     %When the cloud fraction is very low (but not NaN) the CTH is set to a
     %value to flag it as such (-1e9) - since if the CF was zero there could
     %(presumably? Does it get set to NaN? Yes, seems to.) be no CTT value.
        ylabelstr = 'Mean Cloud Top Height (km)';
        
        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
%        subtract_val = -2.5; %Works well using mean CTT from L3 for the region     LAT_val = [-40.5 -30.5]; LON_val = [-140 -100]; 
        
        Y = CTH_Ryan(ilat,ilon,itime);
        
        cf_min_thresh = 0.01; %0.01;
        ilow_cf = find( Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime) < cf_min_thresh);
        Y(ilow_cf) = Ybins(1)+1e-30; %set this as a marker that is not NaN, but make it within range of Ybins - will then discard
        %the first bin
        % It seems likely that the inclusion of points where the CF is very
        % small (and so no CTH reading) is not what is important here - it
        % is more the inclusion of the zeros at all heights above teh
        % observed CTH that makes a difference, I think.
        
        %Need to somehow add an option to have a CTH if there is no cloud
        %(CF=0) - did something for this in the screening routine.
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20; 
        
        %CF profiles for each location are created in case 961 of watervap.
        %(run using case 131 with  multi_plot_case_122 = 'CF_MOD35 vs CTH
        %MOD35').
        % This sets the CF for each profile to one at the CTH and zero
        % above this (with NaNs below). FOr cases where we have no CTH we
        % assume the CF is zero for all heights.
      
        
 case 'Mean CTH timeseries3 flag zeroCF, y-axis'
     %When the cloud fraction is very low (but not NaN) the CTH is set to a
     %value to flag it as such (-1e9) - since if the CF was zero there could
     %(presumably? Does it get set to NaN? Yes, seems to.) be no CTT value.
        ylabelstr = 'Mean Cloud Top Height (km)';
        
        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
%        subtract_val = -2.5; %Works well using mean CTT from L3 for the region     LAT_val = [-40.5 -30.5]; LON_val = [-140 -100]; 
        
        Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) - subtract_val) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
        
        cf_min_thresh = 0.01; %0.01;
        ilow_cf = find( Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime) < cf_min_thresh);
        Y(ilow_cf) = Ybins(1)+1e-30; %set this as a marker that is not NaN, but make it within range of Ybins - will then discard
        %the first bin
        % It seems likely that the inclusion of points where the CF is very
        % small (and so no CTH reading) is not what is important here - it
        % is more the inclusion of the zeros at all heights above teh
        % observed CTH that makes a difference, I think.
        
        %Need to somehow add an option to have a CTH if there is no cloud
        %(CF=0) - did something for this in the screening routine.
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20; 
        
        %CF profiles for each location are created in case 961 of watervap.
        %(run using case 131 with  multi_plot_case_122 = 'CF_MOD35 vs CTH
        %MOD35').
        % This sets the CF for each profile to one at the CTH and zero
        % above this (with NaNs below). FOr cases where we have no CTH we
        % assume the CF is zero for all heights.
        
        
 case 'Mean CTH timeseries3, y-axis'
        ylabelstr = 'Mean Cloud Top Height (km)';
        
        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
        
        Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) - subtract_val) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
        
        %Need to somehow add an option to have a CTH if there is no cloud
        %(CF=0) - did something for this in the screening routine.
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20;        
        
 case 'Mean CTH_no_subtract timeseries3, y-axis'       
        subtract_val = 2.35; %original Zuidema factor - inversion strength?
        subtract_val = -4;
%        subtract_val = 0;
        
        ylabelstr = ['Mean Cloud Top Height (km), subtract val=' num2str(subtract_val)];
        
        Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) - subtract_val) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20;        
        
        
     case 'Cloud Top Height all clouds'
        ylabelstr = 'Mean Cloud Top Height (km)';
        
%        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
        
%        Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) - subtract_val) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
 
Y = CTH_all.timeseries3(ilat,ilon,itime);

        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=-2;
        maxYbins=12; 
        nYpdf = 20;  
        
        
     case 'CTT all clouds'
        ylabelstr = 'Mean Cloud Top Temperature (K)';
        
%        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
        
%        Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) - subtract_val) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
 
Y = Cloud_Top_Temperature_Day_ice_liq_Mean.timeseries3(ilat,ilon,itime);

        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

%         minYbins=-2;
%         maxYbins=12; 
%         nYpdf = 20;      
        
    case 'SST - Mean CTT timeseries3, y-axis'
        ylabelstr = 'SST - Mean CTT (K)';
%        ylabelstr = 'SST - Min CTT (K)';        
        
%        subtract_val = 2.35; %original Zuidema factor - inversion strength?
%        subtract_val = 0;
        
%        Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) - subtract_val) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
         Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime) ) ;      
%         Y = (273.15 + sst_amsre_time3(ilat,ilon,itime) - Cloud_Top_Temperature_Day_Minimum.timeseries3(ilat,ilon,itime) ) ;
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

%        minYbins=-2;
%        maxYbins=12; 
%        nYpdf = 20;    

    case 'SST timeseries3, y-axis'
        ylabelstr = 'SST (K)';
         Y = sst_amsre_time3(ilat,ilon,itime);
        
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

%        minYbins=-2;
%        maxYbins=12; 
%        nYpdf = 20;  
        


    case 'Mean CTP timeseries3, y-axis'
        ylabelstr = 'Mean Cloud Top Pressure (hPa)';
        Y = Cloud_Top_Pressure_Day_Mean.timeseries3(ilat,ilon,itime);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str ' for ' time_mean_str];

        minYbins=100;
        maxYbins=1080;

    case 'Mean Sensor Zenith Angle timeseries3'
        ylabelstr = 'Mean Sensor ZA';
        Y = Sensor_Zenith_Mean.timeseries3(ilat,ilon,itime);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Max Sensor Zenith Angle timeseries3'
        ylabelstr = 'Max Sensor ZA';
        Y = Sensor_Zenith_Maximum.timeseries3(ilat,ilon,itime);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'SZA std dev timeseries2'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'Std dev of SZA';
        Y = Solar_Zenith_Standard_Deviation.timeseries(ilat,:);
        %                        xdat(1).x = N_time(ihtot);

    case 'SZA std dev timeseries3'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'Std dev of SZA';
        Y= Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,itime);
        %                        xdat(1).x = N_time(ihtot);

    case 'Nd from histogram vals timeseries'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'N_d from histo timeseries (cm^{-3})';

        [histo_output]=...
            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,:) );
        Y=histo_output.N_histo_mean;

    case 'Nd from histogram vals timeseries3'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'N_d from histo timeseries3 (cm^{-3})';

        %        [histo_output]=...
        %            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );
        %        Y=histo_output.N_histo_mean;

        Y = Nd_timeseries.mean(ilat,ilon,itime);

    case 'LWP std dev from histogram vals timeseries3'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'Std dev LWP from histo timeseries3 (normalised)';

        Y = W_timeseries.std_dev_norm(ilat,ilon,itime);        

    case 'SZA max from composite'
        extra_title_info = [''];
        ylabelstr = 'Max SZA';

        Y = mid_sza(:,:,itime);
        
     case 'Nd from grid vals timeseries3'
         extra_title_info = [''];
        ylabelstr = 'N_d (cm^{-3})';
       
        Y = N_time3(ilat,ilon,itime);
        
        std_dev = Standard_Deviation_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime); 
        
%            nYpdf=50;

%Ybins = [0:10:90 100:20:280 300:50:1000]; ichoose_Ybins=1;
Ybins = [0:10:1000]; ichoose_Ybins=1;

 case 'CDP Nd'
     extra_title_info = [''];
        
        ylabelstr = 'CDP Nd (cm^{-3})';       
        Y = squeeze(Nd_CDP.timeseries3(ilat,ilon,itime));    %use the indices to remove the data at the end, which
        %is blank
        
        Ybins = [0:5:1000]; ichoose_Ybins=1;   
        
        
        
        

    case 'Nd_{3.7} from grid vals timeseries3'
        ylabelstr = 'N_{d Re 3.7} (cm^{-3})';
       
        Y = N_time3_37(ilat,ilon,itime);
        %        X = N_time3;
        
        std_dev = Standard_Deviation_Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime); 

        Ybins = [0:10:1000]; ichoose_Ybins=1;
        
 case 'Nd_{1.6} from grid vals timeseries3'
        ylabelstr = 'N_{d Re 1.6} (cm^{-3})';
       
        Y = N_time3_16(ilat,ilon,itime);
        %        X = N_time3;
        
        std_dev = Standard_Deviation_Droplet_Number_Concentration_16.timeseries3(ilat,ilon,itime); 

        Ybins = [0:10:1000]; ichoose_Ybins=1;        
        

     case 'Nd from grid vals re*0.8 timeseries3'
        extra_title_info = [''];
        ylabelstr = 'N_d (cm^{-3})';
       
        [N_time3_re80]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,0.8*Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');
        
        Y = N_time3_re80(ilat,ilon,itime);
        
%            nYpdf=50;

%Ybins = [0:10:90 100:20:280 300:50:1000]; ichoose_Ybins=1;
Ybins = [0:10:1000]; ichoose_Ybins=1;


         case 'Nd from individual pixels timeseries3'
             extra_title_info = [''];
             ylabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels using re_{2.1}';
             Y = Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
             
             std_dev = Standard_Deviation_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
             instrument_err = Absolute_Combined_Error_Mean_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
             
             Ybins = [0:25:1000]; ichoose_Ybins=1;
             
          case 'log Nd from individual pixels timeseries3'
             extra_title_info = [''];
             ylabelstr = 'log MODIS N_d (cm^{-3}) - from individual pixels using re_{2.1}';
             Y = log(Droplet_Number_Concentration.timeseries3(ilat,ilon,itime));
             
             std_dev = Standard_Deviation_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
             instrument_err = Absolute_Combined_Error_Mean_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
             
             Ybins = [0:25:1000]; ichoose_Ybins=1;
             
             
         case 'Nd from individual pixels 1.6\mum timeseries3'
             extra_title_info = [''];
             ylabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels using re_{1.6}';
             Y = Droplet_Number_Concentration_16.timeseries3(ilat,ilon,itime);    
             
             std_dev = Standard_Deviation_Droplet_Number_Concentration_16.timeseries3(ilat,ilon,itime);
             instrument_err = Absolute_Combined_Error_Mean_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);
             
        case 'Nd from individual pixels 3.7\mum timeseries3'
             extra_title_info = [''];
             ylabelstr = 'MODIS N_d (cm^{-3}) - from individual pixels using re_{3.7}';
             Y = Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime);       
             
             std_dev = Standard_Deviation_Droplet_Number_Concentration_37.timeseries3(ilat,ilon,itime);  
             instrument_err = Absolute_Combined_Error_Mean_Droplet_Number_Concentration.timeseries3(ilat,ilon,itime);

     case 'Cloud Fraction from grid vals timeseries3'
         extra_title_info = [''];
         ylabelstr = 'Cloud Fraction';        
         Y = Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
         
%         minYbins = 0;
%         maxYbins = 1;

%Ybins = [0:0.1:1.0 1.00001]; ichoose_Ybins=1;
Ybins = [0:0.1:1.00001]; ichoose_Ybins=1;
Ybins = [-0.05:0.1:1.05]; ichoose_Ybins=1;
Ybins = [-0.0001:0.1:0.89999 1.0001]; ichoose_Ybins=1;

extra_title_info = [', LAT=' LAT_str ', LON=' LON_str tim_str];


  case 'MOD35 Cloud Fraction from grid vals timeseries3'
         extra_title_info = [''];
         ylabelstr = 'MOD35 Daytime Cloud Fraction';        
         Y = Cloud_Fraction_Day_Mean.timeseries3(ilat,ilon,itime);
         
%         minYbins = 0;
%         maxYbins = 1;

%Ybins = [0:0.1:1.0 1.00001]; ichoose_Ybins=1;
Ybins = [0:0.1:1.00001]; ichoose_Ybins=1;
Ybins = [-0.05:0.1:1.05]; ichoose_Ybins=1;
%Ybins = [-0.0001:0.1:1.0001]; ichoose_Ybins=1;
Ybins = [-0.00001 0.04999:0.1:0.94999 1.00001]; ichoose_Ybins=1;

extra_title_info = [', LAT=' LAT_str ', LON=' LON_str tim_str];



    case 'Tau^{1/2}'

        ylabelstr = 'Mean Optical Depth ^{1/2}';
        Y = (Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime)).^0.5;
        
 case 'Tau'

        ylabelstr = 'Mean Optical Depth';
        Y = (Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime));    
        
        Ybins = [0:1:50]; ichoose_Ybins=1;
        
    case 'Homogeneity Parameter timeseries3 using W from mean tau and Re'
        ylabelstr = 'Homogeneity Parameter using W from mean tau and Re';
        Y = homog_time3_meanW(ilat,ilon,itime);

        minYbins=0;
        maxYbins=20;

    case 'Homogeneity Parameter timeseries3 using mean W for each pixel'
        ylabelstr = 'Homogeneity Parameter using mean W for each pixel';
        Y = homog_time3_W(ilat,ilon,itime);

%        minYbins=0;
%        maxYbins=20;

    case 'Homogeneity Parameter timeseries3 using Cahalan log mean W (pixel level)'
        ylabelstr = 'Homogeneity Parameter using Cahalan log mean W (pixel level)';
        Y = homog_time3_logW_W(ilat,ilon,itime);

%        minYbins=0;
%        maxYbins=20;     
        
   

    case 'Homogeneity Parameter Cahalan Optical Depth (Seethala)'
        ylabelstr = 'Homogeneity Parameter using Cahalan log mean tau';
        Y = homog_tau_Cahalan(ilat,ilon,itime);

%        minYbins=0;
%        maxYbins=20;   

nYpdf = 125;

 case 'Heterogeneity Parameter Cahalan Optical Depth (Seethala)'
        ylabelstr = 'Heterogeneity Parameter using Cahalan log mean tau';
        ylabelstr = '\gamma_\tau';
        Y = 1 - homog_tau_Cahalan(ilat,ilon,itime);
        
        Ybins=1-fliplr([0.720001:0.02:1.00001]); ichoose_Ybins=1;
%        Ybins=1-fliplr([0.720001:0.0002:1.00001]); ichoose_Ybins=1; %for means better to have lots of bins
        
%        minYbins=0;
%        maxYbins=20;   

        Ybins=1-fliplr([0.720001:0.002:1.00001]); ichoose_Ybins=1;
        
nYpdf = 125;


%tau_ctt = (Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3(ilat,ilon,itime));
%iremove = find(~(tau_ctt < 1.5));
%Y(iremove)=NaN;



        if ~exist('ioverride_pdf') | ioverride_pdf==0
         colocate_POLDER=1;
        end       
        
        if colocate_POLDER==1
            Y(isnan(daymean_Par2_CDR(ilat,ilon,itime)))=NaN;
        end
        
        
        

    case 'Cloud Top Temp standard deviation (normalised by mean CTT), liquid pixels'
        ylabelstr = 'Cloud Top Temp standard deviation (normalised by mean CTT), liquid pixels';
        Y = (Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3(ilat,ilon,itime)./Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime));

%        minYbins=0;
%        maxYbins=20; 

nYpdf = 125;

    case 'Cloud Top Temp standard deviation, liquid pixels'
        ylabelstr = 'Cloud Top Temp standard deviation, liquid pixels (K)';
        ylabelstr = '\sigma_{CTT} (K)';
        Y = (Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3(ilat,ilon,itime));

        minYbins=0;
        maxYbins=20; 
%        maxYbins=6;
        

nYpdf = 125;
nYpdf = 25;



 case 'Cloud Top Temp standard deviation, all pixels'
        ylabelstr = 'Cloud Top Temp standard deviation, ALL pixels (K)';
        Y = (Cloud_Top_Temperature_Day_ice_liq_Standard_Deviation.timeseries3(ilat,ilon,itime));

        minYbins=0;
        maxYbins=20; 

nYpdf = 125;

 case 'CTT homogeneity, liquid pixels'
        ylabelstr = 'Cloud Top Temp homogeneity, liquid pixels';
        Y = (Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3(ilat,ilon,itime)./Cloud_Top_Temperature_Day_Mean.timeseries3(ilat,ilon,itime));

        Y = (1./Y).^2;
        maxval = 1e7;
        Y(Y>maxval)=maxval;
        
        Ybins = 1e7*[0:0.003:0.03 0.032:0.02:0.065 0.075 0.1 0.75 1]; ichoose_Ybins=1;
        Ybins = 1e7*[0:0.0002:0.01 0.013:0.015:0.03 0.032:0.02:0.065 0.075 0.1 0.75 1]; ichoose_Ybins=1;        
        
        Ybins = 10.^([3:0.2:7]); ichoose_Ybins=1; 
%        minYbins=0;
%        maxYbins=20; 

    case 'LWP standard deviation'
        ylabelstr = 'LWP standard deviation (g m^{-2})';
        Y = 1e3*Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3(ilat,ilon,itime);

%        minYbins=0;
%        maxYbins=1000;        
       

 case 'Normalised LWP standard deviation'
        ylabelstr = 'Normalised LWP standard deviation';
        Y = Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3(ilat,ilon,itime) ./ W_time3(ilat,ilon,itime);

%        minYbins=0;
%        maxYbins=1000;  

case 're_{3.7} - re_{2.1} (\mum) reduced dataset Re_1.6 Re_3.7'

        ylabelstr = 're_{3.7} - re_{2.1} (\mum)';       

        
        Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime);
        Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilat,ilon,itime);  

        Y = Re_37 - Re_21;
                
        %make = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
        %slightly different NaN values)
        re_nan = find(isnan(Re_16)==1);
        Y(re_nan) = NaN;
        re_nan = find(isnan(Re_21)==1);
        Y(re_nan) = NaN;
        re_nan = find(isnan(Re_37)==1);
        Y(re_nan) = NaN;
       
end

%in case we want x to be dummy data
switch x_axis_vals
    case 'Dummy data'
               
        X = Y;
        xlabelstr = 'Dummy data';  


        
        ikeep_X_above_zero=1;        
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data        
                        
end

%% End of y-variable switch



if  ndims_hist>2
    switch z_axis_vals
        case 'Z from outside script'
            Z = Z_driver;
%            Zbins = Zbins_driver;
        case 'Nd from grid vals timeseries3'
            zlabelstr = 'N_d (cm^{-3})';
            Z = N_time3(ilat,ilon,itime);

            minZbins=0;
            maxZbins=600;
            
        case 'Tau' %z-axis
            zlabelstr = 'O';
            Z = (Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime));    
        
%        Zbins = [0:1:50]; ichoose_Zbins=1;

            minZbins=0;
            maxZbins=600;    

        case 'Nd absolute uncertainty from grid vals timeseries3'
            zlabelstr = 'N_d absolute uncertainty (cm^{-3})';
            Z = N_un_time3(ilat,ilon,itime);

            minZbins=0;
            maxZbins=600;
    end

end


switch x_axis_vals
    case 'Dummy data for 1D'
       X=Y;     
end


if iswap_xy==1
    y_bk = ylabelstr;
    ylabelstr=xlabelstr;
    xlabelstr=y_bk;

    y_bk = Y;
    Y=X;
    X=y_bk;
    
    Ybins_bk=Ybins;
    Ybins=Xbins;
    Xbins=Ybins_bk;
    
    ichoose_Ybins_bk = ichoose_Ybins;
    ichoose_Xbins = ichoose_Ybins;
    ichoose_Ybins = ichoose_Ybins_bk;

end



%%%%%%%%   now X and Y data are in place. Restrict the data if required


% set_MODIS_NH_flags=1; %gets reset every time
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
% MODIS_N_H_calc %runs the routine to calculate Nd
% N_H_calc_histo
%
% thresh=0.8;






%[MLAT2d,MLON2d]=meshgrid(MLON,MLAT); %these are the mid points e.g. MLAT=89.5 means 89-90 I think
%ihtot = find(abs(MLAT2d)<=30 & cf>=thresh);  thresh_str='LAT LTE 30 & CF GTE ';

%ihtot = find(abs(Solar_Zenith_Mean.data)>=70 & cf>=thresh);  thresh_str='SZA GTE 70 & CF GTE ';
%                        ihtot = find(cf>=thresh);  thresh_str='CF GTE';




%%%%%  choosing only certain points (lat-lon grid)
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
%                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
%                 ihtot = find(WMOD>=thresh);  thresh_str='W times
%                 five-sixths';


%%%%%  choosing only certain points (timeseries)
%                ihtot = find(NP_time>50 & cf_time>0.8);
%                thresh_str='NP.GT.50 AND CF.GT.0.8';

%% Choose thresholds for screening

if ~exist('ioverride_pdf') | ioverride_pdf==0 | ( exist('inotoverride_screening') & inotoverride_screening==1 )

    
    %cf screening is >thresh_CF(1) & <=thresh_CF(2)
    
thresh_SZA=[75 90]; %upper range for SZA paper
thresh_SZA=[50 55];
%thresh_SZA=[45 55];
thresh_SZA=[-1 90];
%thresh_SZA=[0 70];
%thresh_SZA=[70 90];
%thresh_SZA=[0 65];
%thresh_SZA=[50 55];  %lower range for SZA paper
%thresh_SZA=[48.78 52.82];
%thresh_SZA=[75 82];
%thresh_SZA=[46 56];
%thresh_SZA=[55 60];
%thresh_SZA=[50 55];
%thresh_SZA=[45 60];
%thresh_SZA=[45 85];

thresh_CF=[0.799 1.0000001];
%thresh_CF=[0.47 1.0000001];
%thresh_CF=[0.99 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.1 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.99 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.01 0.8];
thresh_CF=[-0.1 1.000001];

thresh_NP=0;
%thresh_NP=10;
thresh_NP=25;
thresh_NP=50;

thresh_ndays=15;

thresh_zeroCF = 0.05;



thresh_sensZA=45;
%thresh_sensZA=50;
thresh_sensZA=[0 20];
thresh_sensZA=[0 41.4];
%thresh_sensZA=[41.4 90];
thresh_sensZA=[-1 90];
%thresh_sensZA=[0 65];
%thresh_sensZA=[55 65];
%thresh_sensZA=[30 70];
%thresh_sensZA=58;

thresh_CTH = [-20 1e9];
%thresh_CTH = [-20 2];
thresh_CTH = [-0.01 3.2];

thresh_AZ = [50 130];

thresh_relAZ = [-1 181];
%thresh_relAZ = [0 80];
%thresh_relAZ = [100 180];

thresh_CTT = [273 273+100];
thresh_CTT = [268 273+100];
thresh_CTT = [173 273+100];
%thresh_CTT = [273 273+100];
%thresh_CTT = [273 273+100];
%thresh_CTT = [273-100 273+100];
%thresh_CTT = [265 273+100];
%thresh_CTT = [250 273+100];
thresh_CTP = 800; %Cloud top pressure (hPa)

thresh_reff = [0 30];
%thresh_reff = [0 14];

thresh_sigCTT = [0 1e9];

%    thresh_sensZA=80;
thresh_maxSZA=200;

thresh_Nd_per_error = 100;
thresh_Reff_per_error = 50;
thresh_Reff_abs_error = 4;
thresh_Reff = 30;

thresh_tau = [-1 300];

%this is set up for homog = (LWP/std_LWP).^2  - so is HOMOGENEITY. Lower
%numbers are more heterogeneous
thresh_stdW = [25 1e9];
thresh_stdW = [0 1e9];


minfrac_CF = 0.9; 
%minfrac_CF = 0.99;  
%minfrac_CF = 0.79;  
%minfrac_CF = 0.47;  
%minfrac_CF = 0.1;  
%minfrac_CF = 0.0;  
%minimum fraction of the sampled points that had successful cloudy/clear/phase
      %determination (i.e. Npix/Nptot_mockL3 =
      %Cloud_Fraction_Liquid_Pixel_Counts./Cloud_Fraction_Liquid./Total_pixels
      % - restriction (2) as presented in the paper
  

minfrac_NpNd = 0.9;
%minfrac_NpNd = 0.99;
%minfrac_NpNd = 0.47;
%minfrac_NpNd = 0.1;
%minfrac_NpNd = 0.0;
%Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3./Cloud_Fraction_Liquid_Pixel_Counts.timeseries3
%Fraction of points that remain after all previous filtering for
%which we have an Nd retrieval. Restriction (4) in the SZA paper.

thresh_dreff = 0.6; %maximum allowed difference between 1.6, 2.1 and 3.7 um retrievals
thresh_dreff = 1e9; %maximum allowed difference between 1.6, 2.1 and 3.7 um retrievals - is set as a percentage at the moment
%thresh_dreff = 10; %maximum allowed difference between 1.6, 2.1 and 3.7 um retrievals - is set as a percentage at the moment

%max allowed height of the upper cloud layer
thresh_maxlayerH = [0 1e5];
%thresh_maxlayerH = [1800 1e5];
thresh_maxlayerH = [0 800];
thresh_nlayers = [1 1]; %max number of allowed layers
%thresh_nlayers = [1 4]; %max number of allowed layers

thresh_calipso_highCF = [-0.01 0.3]; %desired range of Calipso mid+high cloud

%% GCM screening

thresh_gcm_lwp = [-0.01 1]; %g m^{-2}
thresh_gcm_lwp = [-0.01 1e9]; %g m^{-2}
%thresh_gcm_lwp = [1 1e9]; %g m^{-2}

thresh_gcm_tlwp = [-0.01 1]; %g m^{-2}
thresh_gcm_tlwp = [-0.01 1e9]; %g m^{-2}

thresh_gcm_precip = [0.1 1e9]; %mm day^{-1}
thresh_gcm_precip = [-0.01 1e9]; %mm day^{-1}

thresh_gcm_cf = [0.05 1.01]; %
thresh_gcm_cf = [-0.05 1.01]; %

ocean_only_flag = 'Land and ocean';
%ocean_only_flag = 'Ocean only';
ocean_only_flag = 'None';  %only select 'None' for non-GCM plots
%calculate gcm_landmask_rep outside of the override switch - need size of Y


%% Choose screening type

%screen_type='NP + CF + MAX sensZA';
%                                    screen_type='NP + MAX sensZA';
screen_type='NP + CF + MEAN sensZA';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT + noice';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT';
%screen_type ='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';


%                                   screen_type='NP + CF + MAX solarZA';
%screen_type='NP + CF + MEAN solarZA';
%screen_type = 'NP + CF + MEAN solarZA + MEAN sensZA';
screen_type='NP + CF';
%the one used in case 122 for stdCTT plots, etc. for the SZA paper
screen_type='NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';
%screen_type='NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff + seaice';
%screen_type='NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff';

screen_type =  'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff';

%with min CTT over the considered area
%screen_type='NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH';

%with min CTT over the considered area & threshold for dreff (max diff
%between 1.6, 2.1 and 3.7 um)
%screen_type='NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH + dreff';
%when have Irshad ceilometer data
%screen_type='NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH + dreff + ceil_nlayers + ceil_layer_max_height';


%best for level-3 data? Actually the CTT is a bit dubious since it is the
%CTT for all cloud in the 1x1 box - not just the liquid CTT used for the Nd
%So will screen out a lot of cloud
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH, no homog screening';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT';

%screen_type='NP + CF_L3'; %choose this for L3 data (is the same, but put into the 'new' case of the switch)
  % (i.e. for L3 as opposed to mockL3)
%screen_type='NP + CF mockL3';
%                                    screen_type='NP + CF + warm';
%                                    screen_type='NP + CF + min pressure';
%                                    screen_type='NP + CF + min pressure + min temp';
%                                    screen_type='NP + CF + mean pressure';
%                                    screen_type='NP + CF + mean pressure + mean temp';
%                                    screen_type='NP + warm';
%                                  screen_type='NP';
%
%                                    screen_type='CF';

%screen_type = 'am3 lat-lon';

screen_type = 'calipso_high_cloud_ocean_only';
%screen_type = 'gcm_screening';

% --------  NO screening case  -----------------
screen_type ='none';





%%%%%  choosing only certain points (timeseries3)- comment all out for all
%%%%%  points
%timeseries3
%ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF)];
%ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF & tau_time3<thresh_tau); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF) '.AND.tau.LT.' num2str(thresh_tau)];
%ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF & maxSZA_time3<=thresh_maxSZA); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF) ' AND maxSZA.LTE.' num2str(thresh_maxSZA)];

%ihtot = find(NP_time3>thresh_NP);thresh_str=['NP.GT.' num2str(thresh_NP)];
%ihtot = find(cf_time3>thresh_CF); thresh_str=['CF.GT.' num2str(thresh_CF)];
%ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF & sensZA_time3<thresh_SensZA  & tau_time3<thresh_tau); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF) ' AND SensZA.LT.' num2str(thresh_SensZA)  '.AND.tau.LT.' num2str(thresh_tau)];
%ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF & max_sensZA_time3<thresh_SensZA  & tau_time3<thresh_tau); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF) ' AND Max SensZA.LT.' num2str(thresh_SensZA)  '.AND.tau.LT.' num2str(thresh_tau)];
%ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF &
%max_sensZA_time3<thresh_SensZA  & tau_time3<thresh_tau & abs(maxSZA_time3-minSZA_time3)<=thresh_minmax); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF) ' AND SensZA.LT.' num2str(thresh_SensZA)  '.AND.tau.LT.' num2str(thresh_tau) '.AND.min max SZA diff .LT.' num2str(thresh_minmax)];

end



switch screen_type
    case 'calipso_high_cloud_ocean_only'
        switch daynight
            case 'Daytime'
                cll = meanNoNan(cllcalipso_monthly,1)/100;
                clh2 = meanNoNan(clhcalipso_monthly,1)/100;
                clm2 = meanNoNan(clmcalipso_monthly,1)/100;
                clh = clm2 + clh2 - clm2.*clh2;  %random overlap assumption
            case 'Nighttime'
                cll = meanNoNan(cllcalipso_monthly_NIGHTTIME,1)/100;
                clh2 = meanNoNan(clhcalipso_monthly_NIGHTTIME,1)/100;
                clm2 = meanNoNan(clmcalipso_monthly_NIGHTTIME,1)/100;
                clh = clm2 + clh2 - clm2.*clh2;  %random overlap assumption
        end
        
        switch gcm_str
            case 'MODIS'
                Plon2D = repmat(LON(ilon),[length(LAT(ilat)) 1]);
                Plat2D = (repmat(LAT(ilat),[length(LON(ilon)) 1]))';                
                %Don't really need the landmask if using AMSRE
                load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
                
                %Interpolate to teh correct grid
                gcm_landmask_MODIS = griddata(MLON_standard,MLAT_standard,amsre_land_mask,Plon2D,Plat2D);

                
                
        end
        
        %interpolate the CALIPSO data onto the same grid as the GCM
        clh_int = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,clh,Plon2D,Plat2D); %interpolate GCM data onto same grid
        %replicate over the time dimension to match the GCM data
        clh_int_rep = repmat(clh_int,[1 1 size(Y,3)]);
        
        gcm_landmask = eval(['gcm_landmask_' gcm_str ';']);
        gcm_landmask_rep = repmat(gcm_landmask,[1 1 size(Y,3)]);
        
end



switch ocean_only_flag
    case 'Ocean only'
         gcm_landmask = eval(['gcm_landmask_' gcm_str ';']);
         gcm_landmask_rep = repmat(gcm_landmask,[1 1 size(Y,3)]);
    case 'Land and ocean'
         gcm_landmask_temp = eval(['gcm_landmask_' gcm_str ';']);
         gcm_landmask_temp(:)=0; %remove all land so that all points are selected
         gcm_landmask_rep = repmat(gcm_landmask_temp,[1 1 size(Y,3)]);         
end




% *** external script *** %
iplot_global=0;
modisL3_screening_timeseries3;
% ************************************  %

%% End of screening

%only select the points required
X2=X;
X=NaN*ones(size(X2));
X(ihtot)=X2(ihtot);
Y2=Y;
Y=NaN*ones(size(Y2));
Y(ihtot)=Y2(ihtot);
%in case the data is [1 1 N]
X=squeeze(X);
Y=squeeze(Y);

if ndims_hist==3
    Z=Z(ihtot);
    Z=squeeze(Z);
end


if ireduce_res_lwp==1
    %Reduce the resolution by averaging over boxes of size M_lwp x
    %N_lwp
    Y1 = reduce_matrix_subsample_mean(Y(:,:,1),M_lwp,N_lwp);
    siz_Y1 = size(Y1);
    siz_tim = size(Y,3);
    Y2 = NaN*ones([siz_Y1 siz_tim]);
    for ired=1:siz_tim
        Y2(:,:,ired) = reduce_matrix_subsample_mean(Y(:,:,ired),M_lwp,N_lwp);
    end
    Y=Y2;
    X=Y2; %Will only use this when x is dummy data
end





%plot_type
short_plot_name=[x_axis_vals ' vs ' y_axis_vals ];
tit(1).tit=[short_plot_name extra_title_info];
tit(1).tit=[tit(1).tit '_' gcm_str];
savename=[savedir tit(1).tit];







if ichoose_Xbins==0
    if minXbins<-8.9e99
        Xbins=make_PDF_bins(X,nXpdf); %if not set limits then use the default of minALL(X):dX:maxALL(X)
    else
        Xbins=make_PDF_bins(X,nXpdf,minXbins,maxXbins);
    end
end

if ichoose_Ybins==0
    if minYbins<-8.9e99
        Ybins=make_PDF_bins(Y,nYpdf);
    else
        Ybins=make_PDF_bins(Y,nYpdf,minYbins,maxYbins);
        %    Ybins=[-0.01:0.1:0.89 0.99 1.09];
    end
end

if ndims_hist==3
    if ichoose_Zbins==0
        if minZbins<-8.9e99
            Zbins=make_PDF_bins(Z,nZpdf);
        else
            Zbins=make_PDF_bins(Z,nZpdf,minZbins,maxZbins);
        end
    end
end

%Adding this becauase Ybins may not have existed when this was done above
switch x_axis_vals
    case 'Dummy data'
        Xbins = Ybins; ichoose_Xbins=1;
end


%clear flags in case of no errors
clear ioverride_pdf ioverride_pdf_varchoose ioverride_location_selection inotoverride_screening
iswap_xy=0; %flag to swap the defined x and y

%and in case of errors too
catch pdf2D_ERROR
    %clear the flags in case of a runtime error
    clear ioverride_pdf ioverride_pdf_varchoose ioverride_location_selection inotoverride_screening
    iswap_xy=0; %flag to swap the defined x and y
    rethrow(pdf2D_ERROR); %"re-issue" the error (also creates the links to the error etc)
end

