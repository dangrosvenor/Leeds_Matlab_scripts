function [success]=L2_test_file(filename_h5)
success=1;

INFO = hdfinfo(filename_h5);
SD_id = hdfsd('start',filename_h5,'read'); %open the file
[solar_zenith,dimsizes_solar_zenith]=get_hdf_data_dan('Solar_Zenith',SD_id,INFO); %int16, min=-32767

i=find(solar_zenith>81.5 | isnan(solar_zenith));
if length(i)==prod(size(solar_zenith))
    success = 0;
end


%evalc('[lat_10km,dimsizes_lat,sampling_along_10km,sampling_across_10km,success]=get_hdf_data_dan(''Latitude_10km'',SD_id,INFO);');
%evalc('[tau,dimsizes_tau,tmp,tmp,success_tau]=get_hdf_data_dan(''Cloud_Optical_Thickness'',SD_id,INFO);'); %int16 (>0)

%success = min([success success_tau]);

%close the hdf file
status = hdfsd('end',SD_id);
%success is returned by the above

