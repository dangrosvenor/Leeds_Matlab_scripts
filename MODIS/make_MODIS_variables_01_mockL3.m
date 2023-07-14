%variables to make for the mock L3 data
%Called from MODIS_process_multiple_L2_files.m
%See also mock_L3_from_L2.m

inds_str = '(inds_swath)'; %the string expression for below
inds_str_cell = '{inds_swath}'; %the string expression for below

%use sensor_zenith as an indicator of the region that this swath covered in
%terms of L3 1x1 deg boxes
%will increment the swath counter only at these locations
%and will only put data in array for these points too (will already
%be NaN so don't need to worry about not over-writing them)
inotnan = find(isnan(meanSensorZA_mockL3)==0);

if length(inotnan)==0
   return     
end
%size of the 2D data
sdat=size(meanSensorZA_mockL3);
%size of the 3D data (2D * No. swaths)
sdat2=size(Sensor_Zenith_Mean.timeseries3);
[inotnanX,inotnanY]=ind2sub(sdat,inotnan);

switch action2
    case {'make mock L3 data','make mock L3 all days'}
%        inds_str = '(:,:,iswath)';
%make the 3rd index equal to iswath, the swath number if are storing all of
%them
        ISWATH = iswath*ones([length(inotnanX) 1]);
        inds_swath=sub2ind(sdat2,inotnanX,inotnanY,ISWATH);
        
        
    case 'make mock L3 data daily'

%Want to create arrays of size [ndays nlat nlon Nswath_dim]
%(will want to pre-allocate). ndays can be worked out from the filenames
%Nswath_dim is the max number of orbits possible in one day
%Then will also need to put each swath into the right day, which should be
%easy. The tricky part will be putting in the right orbit position without
%looping through all the lat lons - would need to know the next free
%location for each lat lon for each time. Then can prob use sub2ind, e.g.
%a=magic([5 5]); b=[50 50 50 50 50]; i=[1 1 2 2 5]; ii=[1:5]; 
%inds=sub2ind(size(a),ii,i); a(inds)=b; puts b into a with different
%indices (i) for each row
%will also want to include the sea ice data



nswath(inotnan)=nswath(inotnan)+1;
%Hack to prevent script failures when nswath>20 - just tell it to write in
%no. 20 again. Could increase arrays to >20, but might cause other issues.
%Think only goes >20 swaths on day 167 of 2010.
i20=find(nswath>sdat2(3));
nswath(i20)=sdat2(3);



%nswath is the position to place data at each lat lon
%only select the positions that have data for this swath, otherwise will be
%overwriting with NaNs or nswath will be zero
ISWATH = nswath(inotnan);
%ISWATH = ISWATH(:);

%ix = [1:sdat(1)]; ix=ix(inotnanX);
%iy = [1:sdat(2)]; iy=iy(inotnanY);

%make a set of indices for lat and lon for all locations in the 2D grid
%[IX,IY]=meshgrid(inotnanX,inotnanY);
%transpose because the (:) operation moves through the rows fastest
%but the meshgrid will create an IX that does the opposite
%IX=IX'; IY=IY';
%find the linear index of the combinations of the 2D index and the position
%in the 3rd dimension as given by nswath
inds_swath=sub2ind(sdat2,inotnanX,inotnanY,ISWATH);
%now we can just use dat(inds) = dat2D to insert the 2D map in the desired
%location

end



eval(['Date_Time_Swath.timeseries3' inds_str '= modis_date_time;']); %the Matlab based numerical datenum
%value for the swath date and time

%MOD06 Cloud Fractions
%Npmask_mockL3 is the number of pixels for which there was a successful determinations of whether 
%the pixel was cloudy or not. So CF = n_cloudy./Npmask_mockL3
%Cloud_Fraction_Liquid_old.timeseries3(:,:,iswath) = Nliq_mockL3 ./ Npmask_mockL3; %liquid only
%Cloud_Fraction_Liquid_Pixel_Counts_old.timeseries3(:,:,iswath) = Npmask_mockL3; %no, pixels used in CF determination

%basing on how it is actually done by the MODIS people (using just the
%phase determination and the primary cloud retrieval outcome flags)
eval(['Npix = Nliq_mockL3_new(inotnan) + Nice_mockL3_new(inotnan) + Nundet_mockL3_new(inotnan) + N_clear_new(inotnan);']); 


%eval(['Cloud_Fraction_Liquid_Pixel_Counts.timeseries3' inds_str ' = Nliq_mockL3_new(inotnan) + Nice_mockL3_new(inotnan) + Nundet_mockL3_new(inotnan) + N_clear_new(inotnan);']); 
%no, pixels for which there was a successful phase determination, or which
%were successfully determined to be clear

%no pixels which were liquid phase:-
eval(['Cloud_Fraction_Liquid_Pixel_Counts.timeseries3' inds_str ' = Nliq_mockL3_new(inotnan);']);
%below is actual number of points that had successful Nd retrievals - can
%use this to calculate a (more stringent) CF:-
eval(['Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3' inds_str ' = Np_mockL3(inotnan);']); 

eval(['Total_pixels.timeseries3' inds_str ' = Nptot_mockL3(inotnan);']); 
eval(['Pixel_fraction_CF.timeseries3' inds_str ' = Nliq_mockL3_new(inotnan)./Nptot_mockL3(inotnan);']); 
eval(['Pixel_fraction_Nd.timeseries3' inds_str ' = Np_mockL3(inotnan)./Nptot_mockL3(inotnan);']); 

eval(['Np_single_liquid_layer.timeseries3' inds_str ' = N_single_liquid_layer_mockL3(inotnan);']); 



%Npix = eval(['Cloud_Fraction_Liquid_Pixel_Counts.timeseries3' inds_str ';']);

%if the above is low due to the swath cutting across a grid box then may
%want to reject
%Cloud_Fraction_Ice.timeseries3' inds_str ' = Nice_mockL3 ./ Npmask_mockL3;
%Cloud_Fraction_Undetermined.timeseries3' inds_str ' = Nundet_mockL3 ./ Npmask_mockL3;
%Cloud_Fraction_Combined.timeseries3' inds_str ' = CF_mockL3; %note, this doesn't equal the sum of the other CFs

%Cloud fractions used to try to match Rosenfeld (2019) - less conservative
%masks
% I think this is the MOD35 cloud fraction :-
eval(['Cloud_Fraction_MOD35_L2.timeseries3' inds_str ' = CF_mockL3(inotnan);']); %CF directly from L2 file ('Cloud_Fraction' variable from HDF file)
% These are the fractions from the probably cloudy and confident cloudy
% pixels from the cloud mask. Would need to add the two to get all cloud considered as probable or confident.
eval(['Cloud_Fraction_mask_probable.timeseries3' inds_str ' = N_cloudy_probable_mockL3(inotnan) ./ Npmask_mockL3(inotnan);']); %
eval(['Cloud_Fraction_mask_confident.timeseries3' inds_str ' = N_cloudy_confident_mockL3(inotnan) ./ Npmask_mockL3(inotnan);']); %
%Looking at the results of this, it looks like the MOD35 cloud
%mask is very close to the probable+confident masks.

eval(['Cloud_Fraction_NoOpt.timeseries3' inds_str ' = N_cloudy_no_opt_mockL3(inotnan) ./ Npmask_mockL3(inotnan);']);
%above is the fraction of pixels where the cloud mask identified cloud, but no optical
%retrievals could be made (no tau confidence, etc)

eval(['Cloud_Fraction_Liquid.timeseries3' inds_str ' = Nliq_mockL3_new(inotnan) ./ Npix;']); %liquid only
eval(['Cloud_Fraction_Liquid2.timeseries3' inds_str ' = Np_mockL3(inotnan) ./ Npix;']); %
eval(['Cloud_Fraction_Ice.timeseries3' inds_str ' = Nice_mockL3_new(inotnan) ./ Npix;']);
eval(['Cloud_Fraction_Undetermined.timeseries3' inds_str ' = Nundet_mockL3_new(inotnan) ./ Npix;']);
eval(['Cloud_Fraction_Combined.timeseries3' inds_str ' = (Nliq_mockL3_new(inotnan) + Nice_mockL3_new(inotnan) + Nundet_mockL3_new(inotnan)) ./ Npix;']); %





eval(['Cloud_Effective_Radius_Liquid_Mean.timeseries3' inds_str ' = meanRe_mockL3(inotnan);']);
eval(['Cloud_Effective_Radius_16_Liquid_Mean.timeseries3' inds_str ' = meanRe_mockL3_16(inotnan);']);
eval(['Cloud_Effective_Radius_37_Liquid_Mean.timeseries3' inds_str ' = meanRe_mockL3_37(inotnan);']);
eval(['Cloud_Effective_Radius_Liquid_Standard_Deviation.timeseries3' inds_str ' = stdRe_mockL3(inotnan);']);
eval(['Cloud_Effective_Radius_16_Liquid_Standard_Deviation.timeseries3' inds_str ' = stdRe_mockL3_16(inotnan);']);
eval(['Cloud_Effective_Radius_37_Liquid_Standard_Deviation.timeseries3' inds_str ' = stdRe_mockL3_37(inotnan);']);
eval(['Cloud_Effective_Radius_Liquid_Minimum.timeseries3' inds_str ' = minRe_mockL3(inotnan);']);
eval(['Cloud_Effective_Radius_16_Liquid_Minimum.timeseries3' inds_str ' = minRe_mockL3_16(inotnan);']);
eval(['Cloud_Effective_Radius_37_Liquid_Minimum.timeseries3' inds_str ' = minRe_mockL3_37(inotnan);']);
eval(['Cloud_Effective_Radius_Liquid_Maximum.timeseries3' inds_str ' = maxRe_mockL3(inotnan);']);
eval(['Cloud_Effective_Radius_16_Liquid_Maximum.timeseries3' inds_str ' = maxRe_mockL3_16(inotnan);']);
eval(['Cloud_Effective_Radius_37_Liquid_Maximum.timeseries3' inds_str ' = maxRe_mockL3_37(inotnan);']);
eval(['Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3' inds_str ' = meanRe_un_mockL3(inotnan);']);
eval(['Cloud_Effective_Radius_Liquid_Log_Mean_Uncertainty.timeseries3' inds_str ' = meanRe_logun_mockL3(inotnan);']);
%optical thickness
eval(['Cloud_Optical_Thickness_Liquid_Mean.timeseries3' inds_str ' = meanTau_mockL3(inotnan);']);
eval(['Cloud_Optical_Thickness_Liquid_Standard_Deviation.timeseries3' inds_str ' = stdTau_mockL3(inotnan);']);
eval(['Cloud_Optical_Thickness_Liquid_Minimum.timeseries3' inds_str ' = minTau_mockL3(inotnan);']);
eval(['Cloud_Optical_Thickness_Liquid_Maximum.timeseries3' inds_str ' = maxTau_mockL3(inotnan);']);
eval(['Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3' inds_str ' = meanTau_un_mockL3(inotnan);']);
eval(['Cloud_Optical_Thickness_Liquid_Log_Mean_Uncertainty.timeseries3' inds_str ' = meanTau_logun_mockL3(inotnan);']); 

eval(['Cloud_Water_Path_Liquid.timeseries3' inds_str ' = meanW_mockL3(inotnan);']);
eval(['Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3' inds_str ' = stdW_mockL3(inotnan);']);
eval(['Cloud_Water_Path_Liquid_Log_Mean.timeseries3' inds_str ' = meanlogW_mockL3(inotnan);']);

eval(['Cloud_Optical_Thickness_Liquid_Log_Mean.timeseries3' inds_str ' = meanTau_log_mockL3(inotnan);']); 
eval(['Cloud_Effective_Radius_Liquid_Log_Mean.timeseries3' inds_str ' = meanRe_log_mockL3(inotnan);']);
eval(['Cloud_Effective_Radius_16_Liquid_Log_Mean.timeseries3' inds_str ' = meanRe_log_mockL3_16(inotnan);']);
eval(['Cloud_Effective_Radius_37_Liquid_Log_Mean.timeseries3' inds_str ' = meanRe_log_mockL3_37(inotnan);']);

eval(['Droplet_Number_Concentration.timeseries3' inds_str ' = meanNd_mockL3(inotnan);']);
eval(['Droplet_Number_Concentration_16.timeseries3' inds_str ' = meanNd_mockL3_16(inotnan);']);
eval(['Droplet_Number_Concentration_37.timeseries3' inds_str ' = meanNd_mockL3_37(inotnan);']);
eval(['Standard_Deviation_Droplet_Number_Concentration.timeseries3' inds_str ' = stdNd_mockL3(inotnan);']);
eval(['Standard_Deviation_Droplet_Number_Concentration_16.timeseries3' inds_str ' = stdNd_mockL3_16(inotnan);']);
eval(['Standard_Deviation_Droplet_Number_Concentration_37.timeseries3' inds_str ' = stdNd_mockL3_37(inotnan);']);
eval(['Percent_Error_Mean_Droplet_Number_Concentration.timeseries3' inds_str ' = meanNd_un_mockL3(inotnan);']);
eval(['Absolute_Combined_Error_Mean_Droplet_Number_Concentration.timeseries3' inds_str ' = meanNd_un_combined_mockL3(inotnan);']);

eval(['Cloud_Top_Temperature_Day_Mean.timeseries3' inds_str ' = meanCTT_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_Minimum.timeseries3' inds_str ' = minCTT_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_Maximum.timeseries3' inds_str ' = maxCTT_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3' inds_str ' = stdCTT_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_ice_liq_Mean.timeseries3' inds_str ' = meanCTT_ALL_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_ice_liq_Standard_Deviation.timeseries3' inds_str ' = stdCTT_ALL_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_CF_Mean.timeseries3' inds_str ' = meanCTT_CF_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Day_CF_Standard_Deviation.timeseries3' inds_str ' = stdCTT_CF_mockL3(inotnan);']);


%CTT properties for the ice cloud fraction only
eval(['Cloud_Top_Temperature_Ice_Day_Mean.timeseries3' inds_str ' = meanCTT_ice_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Ice_Day_Minimum.timeseries3' inds_str ' = minCTT_ice_mockL3(inotnan);']);
eval(['Cloud_Top_Temperature_Ice_Day_Maximum.timeseries3' inds_str ' = maxCTT_ice_mockL3(inotnan);']);



%not bothering with min,max and std dev for SZA, etc. since we won't be including
%multiple overpasses
eval(['Solar_Zenith_Mean.timeseries3' inds_str ' = meanSolarZA_mockL3(inotnan);']);
eval(['Sensor_Zenith_Mean.timeseries3' inds_str ' = meanSensorZA_mockL3(inotnan);']);
eval(['Solar_Azimuth_Mean.timeseries3' inds_str ' = meanSolarAz_mockL3(inotnan);']);
eval(['Sensor_Azimuth_Mean.timeseries3' inds_str ' = meanSensorAz_mockL3(inotnan);']);

if length(inds_swath)>0
    %or could use the deal function here - e.g. a=cell([2 3]); [a{[2 5]}] = deal('xx','yx')
    for icell = 1:length(inds_swath)
        
%    file_str_inds = eval(['repmat({''' file_name_h5 '''},[1 size(inds_swath)]);']);
    MODIS_swath_filename.timeseries3{inds_swath(icell)} = file_name_h5;
    end
end
        
%CDP and Puijo related variables
if exist('Nd_Puijo')
        
        eval(['weather_index.timeseries3' inds_str ' = iweather;']);
        eval(['Nd_CDP.timeseries3' inds_str ' = Nd_Puijo(time_match_cloud_legs{imod_read});']);
        eval(['stdNd_CDP.timeseries3' inds_str ' = stdNd_Puijo(time_match_cloud_legs{imod_read});']);       
        eval(['LWC_CDP.timeseries3' inds_str ' = LWC_Puijo(time_match_cloud_legs{imod_read});']);               
end

%these will  now be the box centres
MLAT = 0.5*(LATS(1:end-1) + LATS(2:end) );
MLON = 0.5*(LONS(1:end-1) + LONS(2:end) );

%MLAT=LATS;
%MLON=LONS;




