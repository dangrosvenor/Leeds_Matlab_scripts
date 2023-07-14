% 29th May, 2015 -- Script to produce the Nd dataset for Stefan Kinne using
% the mockL3 data. 
% Have this for one year only (Dec 2006-Nov 2007)

%So have saved an all-CF and a CF>80% dataset as produced by this script, here:-
% all CF :-
%  /home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_minCTT_173_ice_allowed_SZA_65//Stefan_Kinne_dataset_monthly_SZA65_allCF_Dec2006-Nov2007.mat
% CF>80% :-
% /home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65//Stefan_Kinne_dataset_monthly_SZA65_CF80_Dec2006-Nov2007.mat.nc

iplot=1;
isave=1;

Date_text = '1st Dec 2006 to 30th Nov 2007';
Notes =['N.B. Months 1-12 refer to Jan-Dec, although the December is actually from 2006 and the other months from 2007.'...
    'Global droplet concentrations from MODIS using technique described in Grosvenor & Wood, ACP,2014, doi:10.5194/acp-14-7291-2014. '...
    'Data is given for the three MODIS near_infrared bands :- _16 refers to 1.6um, _21 = 2.1um and _37 = 3.7um. 3.7um is perhaps more trusted. '...
    'Daily data was aggregated from MODIS Terra Collection 5.1 L2 swath data at 1km resolution into 1x1 degree regions. Only data with SZAs<65 degrees is included. ' ...
    'Regions containing sea-ice have been excluded using Nimbus-7 daily seaice data from NSDIC from http://nsidc.org/data/docs/daac/nsidc0051_gsfc_seaice.gd.html '...
    'For the *_ocean variables, a landmask has also been applied to give ocean only regions.'...
    'The daily 1x1 degree datapoints were averaged into monthly values. The number of such days for each monthly mean is in ndays_XX. '...
    'Only locations for which ndays_XX >4 have been included. Fill values are set to NaN. '...
    'Standard deviations of the ndays_XX daily values comprising each mean are given in std_Nd_XX_ocean. '...
    ];

halo=10; %size of halo to add to landmask and seaice for smoothing
mask_size = 3; %NxN filter for averaging

choose_CF = 'all CF';
%choose_CF = 'CF_gt_80';
switch choose_CF
    case 'CF_gt_80'

        %File given to Dan McCoy for Part2 of his first paper - has CF>0.8
        %screening, but not other screening.
        % CF > 80%
        file01 = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65/mockL3_saved_data_20121017T115935.mat';
        savedir01 = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65/';
        savefile01 = [savedir01 '/Stefan_Kinne_dataset_monthly_SZA65_CF80_Dec2006-Nov2007.mat'];
        
        Notes = [Notes ' Only 1x1 degree daily datapoints with a mean liquid cloud fraction greater than 80% have been included in this data. '...
            'This is due to concerns about the retrieval of effective radius and optical depth in regions of heterogeneous clouds (e.g. 3D radiative effects). '...
            'Also, this tends to restrict the retrievals to stratocumulus regions where many of the assumptions have been well validated. '...
            'E.g. the adibaticity of the LWC profiles are likely to be close to those assumed by the Nd retrieval technique. '...
            ]
        

    case 'all CF'
        % all CF
        file01 = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_minCTT_173_ice_allowed_SZA_65/mockL3_saved_data_20150603T051150.mat';
        savedir01 = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_minCTT_173_ice_allowed_SZA_65/';
        savefile01 = [savedir01 '/Stefan_Kinne_dataset_monthly_SZA65_allCF_Dec2006-Nov2007.mat'];
        Notes = [Notes ' No filtering based on liquid cloud fraction has been performed on this data. '...
            'There are concerns about the retrieval of effective radius and optical depth in regions of heterogeneous clouds (e.g. 3D radiative effects). '...
            'Consider using the CF>80% data. Also, this tends to restrict the retrievals to stratocumulus regions where many of the assumptions have been well validated. '...
            'E.g. the adibaticity of the LWC profiles are likely to be close to those assumed by the Nd retrieval technique. '...
            'However, the global coverage is much less as a result than for the no CF filtering dataset. '...
            ]        
end

Nd_CF80 = load(file01,'N_time3_16','N_time3','N_time3_37','MLAT','MLON','modisyear_timeseries3','daynum_timeseries3','seaice_max_time3','seaice_time3');
% *** Have combined the NH seaice here and saved the results in the file01
% above - the steps follwed are below:-
% %% Combine NH and SH seaice (is only SH in the file)
% seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2007_all_NHemisphere_20150529T103046.mat';
% load(seaice_fileload_1deg_NH);
% %Load these to do the time matching
% load(file01,'modisyear_timeseries3','daynum_timeseries3');
% seaice_match_times
% %So will now have NaN values in both SH and NH files. Will set the NaNs to zero for
% %now and then make the places where both were NaN.
% % -- Actually the seaice_time3 files don't seem to contain any NaNs...
% seaice_time3_NH = seaice_time3;
% seaice_max_time3_NH = seaice_max_time3;
% %inanNH = isnan(seaice_time3_NH); seaice_time3_NH(inanNH)=0;
% %inanSH = isnan(Nd_CF80.seaice_time3); Nd_CF80.seaice_time3(inanSH) = 0;
% %inanNH_max = isnan(seaice_max_time3_NH); seaice_max_time3_NH(inanNH_max)=0;
% %inanSH_max = isnan(Nd_CF80.seaice_max_time3); Nd_CF80.seaice_max_time3(inanSH_max)=0;
% 
% %iboth = inanNH & inanSH; %When both are NaN - need to use the logical since they are not just integers
% %iboth_max = inanNH_max & inanSH_max; %When both are NaN - need to use the logical since they are not just integers
% 
% seaice_time3 = seaice_time3_NH + Nd_CF80.seaice_time3;
% %seaice_time3(iboth) = NaN;
% 
% seaice_max_time3 = seaice_max_time3_NH + Nd_CF80.seaice_max_time3;
% %seaice_max_time3(iboth_max) = NaN;
% 
% save(file01,'-APPEND','-V7.3','seaice_time3','seaice_max_time3');

%days of month funciton will be useful to give days required for each month
%- e.g. 32:60 for Feb.
% [days_required_for_mean,time_mean_str] = days_of_month(month_no)


landmask_load= load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
landmask = repmat(flipdim(landmask_load.amsre_land_mask,1),[1 1 12]);
%Smooth out the land mask to make it more aggressive - seemed to be
%allowing some small bits of land otherwise.
%Also need to do this for the sea ice

fs=fspecial('average',mask_size); %Make the averaging filter
landmask_tmp = add_halo(landmask_load.amsre_land_mask,halo); 
landmask_smooth = filter2(fs,landmask_tmp,'same');
landmask_smooth = remove_halo(landmask_smooth,halo);
landmask2 = repmat(flipdim(landmask_smooth,1),[1 1 12]);


%Land mask is zero where there is land and NaN where not - so can just add to the
%array to make the non-land regions NaN

%Doing for months 1:12, but bear in mind that the data actually runs from
%Dec2006 to Nov2007

clear idays Nd_16 Nd_21 Nd_37 sea_ice_max_tmp sea_ice_tmp



for imonth=1:12
    [days_required_for_mean,time_mean_str] = days_of_month(imonth);

    
    icount=0;
    for i=1:length(days_required_for_mean)
        ii = find(Nd_CF80.daynum_timeseries3 == days_required_for_mean(i));
        if length(ii)>0
            icount=icount+1;
            idays(icount) = ii;
        end
    end
    
    [ Nd_16(:,:,imonth), ndays_16(:,:,imonth), std_Nd_16(:,:,imonth) ] = meanNoNan(Nd_CF80.N_time3_16(:,:,idays),3);
    [ Nd_21(:,:,imonth), ndays_21(:,:,imonth), std_Nd_21(:,:,imonth) ] = meanNoNan(Nd_CF80.N_time3(:,:,idays),3);
    [ Nd_37(:,:,imonth), ndays_37(:,:,imonth), std_Nd_37(:,:,imonth) ] = meanNoNan(Nd_CF80.N_time3_37(:,:,idays),3);  


%    count_ndays = ones(size(Nd_CF80.N_time3_37(:,:,idays)));
%    count_ndays( isnan(Nd_CF80.N_time3_37(:,:,idays)) )= 0;
%    ndays_data(:,:,imonth) = sum(count_ndays,3);

    
    %Sea-ice is only present for the SH
    %The max value is that from the time of the daily data +/- one week
    sea_ice_tmp(:,:,imonth) = meanNoNan(Nd_CF80.seaice_time3(:,:,idays),3);
    sea_ice_max_tmp(:,:,imonth) = meanNoNan(Nd_CF80.seaice_max_time3(:,:,idays),3);
    
    %Make a longitude halo to allow smoothing over the edge of the map
    sea_ice_tmp2 = add_halo(sea_ice_tmp(:,:,imonth),halo);
    sea_ice_max_tmp2 = add_halo(sea_ice_max_tmp(:,:,imonth),halo);
    
    %Smooth the sea-ice spatially to be conservative
    sea_ice_tmp3 = filter2(fs,sea_ice_tmp2,'same');
    sea_ice_max_tmp3 = filter2(fs,sea_ice_max_tmp2,'same');
    
    %Remove halo
    sea_ice(:,:,imonth) = remove_halo(sea_ice_tmp3,halo);
    sea_ice_max(:,:,imonth) = remove_halo(sea_ice_max_tmp3,halo);
    
    time_mean_str_save{imonth}=time_mean_str;
    
   
end

 %find points where there is any sea-ice using the threshold below
    sea_ice_thresh = 1e-5;
%    sea_ice_thresh = 1e9;
    isea_high = find(sea_ice>sea_ice_thresh);
%    isea_high = find(sea_ice_max>sea_ice_thresh);    
    Nd_16(isea_high)=NaN;
    Nd_21(isea_high)=NaN;
    Nd_37(isea_high)=NaN;
    
 %Exclude points below a threshold number of days included
    thresh_ndays=5;
    
    iNdays = find(ndays_16<thresh_ndays); Nd_16(iNdays)=NaN;
    iNdays = find(ndays_21<thresh_ndays); Nd_21(iNdays)=NaN;
    iNdays = find(ndays_37<thresh_ndays); Nd_37(iNdays)=NaN;    
    
    
 %Make the land regions NaN
    Nd_16_ocean = Nd_16 + landmask2;
    Nd_21_ocean= Nd_21 + landmask2;    
    Nd_37_ocean = Nd_37 + landmask2;
    std_Nd_16_ocean = std_Nd_16 + landmask2;
    std_Nd_21_ocean= std_Nd_21 + landmask2;    
    std_Nd_37_ocean = std_Nd_37 + landmask2;    
    

    
    for imonth=[1:12]

        if iplot==1
            figure
            pcolor(Nd_CF80.MLON,Nd_CF80.MLAT,Nd_37(:,:,imonth)); shading flat; caxis([0 200]); colorbar
%            pcolor(Nd_CF80.MLON,Nd_CF80.MLAT,std_Nd_37_ocean(:,:,imonth)); shading flat; caxis([0 200]); colorbar            
            title(time_mean_str_save{imonth})
        end

    end
    
    
        clear vars_to_save
        i=1;
        vars_to_save{i}='Date_text'; i=i+1;   
        vars_to_save{i}='Notes'; i=i+1;           
        vars_to_save{i}='MLAT'; i=i+1;
        vars_to_save{i}='MLON'; i=i+1;
        vars_to_save{i}='Nd_16_ocean'; i=i+1;
        vars_to_save{i}='Nd_21_ocean'; i=i+1;
        vars_to_save{i}='Nd_37_ocean'; i=i+1;
        vars_to_save{i}='std_Nd_16_ocean'; i=i+1;
        vars_to_save{i}='std_Nd_21_ocean'; i=i+1;
        vars_to_save{i}='std_Nd_37_ocean'; i=i+1;
        vars_to_save{i}='Nd_16'; i=i+1;
        vars_to_save{i}='Nd_21'; i=i+1;
        vars_to_save{i}='Nd_37'; i=i+1;
        vars_to_save{i}='std_Nd_16'; i=i+1;
        vars_to_save{i}='std_Nd_21'; i=i+1;
        vars_to_save{i}='std_Nd_37'; i=i+1;        
        vars_to_save{i}='ndays_16'; i=i+1;        
        vars_to_save{i}='ndays_21'; i=i+1;        
        vars_to_save{i}='ndays_37'; i=i+1;                
        
 %Save here
    if isave==1
        filename_savevars = savefile01;
        tag='';
        save_vars_mat
        mat2nc_Dan(filename_savevars,[filename_savevars '.nc']);
    end
