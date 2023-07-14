sat_case = 'AQUA';
%sat_case = 'TERRA';

switch sat_case
    case 'AQUA'
        VOCALS_write_file = '/home/disk/eos1/d.grosvenor/modis_work/VOCALS_AQUA_L2_flight_leg_matches.txt';
        VOCALS_read_file = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/VOCALS_10-40S_100-60W_12Oct-17Nov_2008_AQUA_L2_VOCALS_matches_20111116T093256.mat';
    case 'TERRA'
        VOCALS_write_file = '/home/disk/eos1/d.grosvenor/modis_work/VOCALS_TERRA_L2_flight_leg_matches.txt';
        VOCALS_read_file = '/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/VOCALS_10-40S_100-60W_12Oct-17Nov_2008_TERRA_L2_VOCALS_matches_20111115T162311.mat';
end

%load in the VOCALS flight leg info
VOCALS_cloud_leg_info  %(external script)

load(VOCALS_read_file);        
        
fid=fopen(VOCALS_write_file,'wt');

fprintf(fid,'%s\t%s\t%s\t%s\t%s\n','Leg_No.','Sat_time_minus_flight_time_(mins)','Distance_diff_(km)','Leg_LAT','Leg_LON');
for i=1:length(VOCALS_leg_number_matches)
    fprintf(fid,'%s %s %s:%s :-\n',VOCALS_leg_number_matches(i).filename(73:end),VOCALS_leg_number_matches(i).sat_date_str,VOCALS_leg_number_matches(i).sat_hour_str,VOCALS_leg_number_matches(i).sat_min_str);
    for j=1:length(VOCALS_leg_number_matches(i).time_inds)
        voc_ind = VOCALS_leg_number_matches(i).time_inds(j);
        fprintf(fid,'%d\t%f\t%f\t%f\t%f\n',VOCALS_leg_number_matches(i).time_inds(j),VOCALS_leg_number_matches(i).time_diff_mins(j),VOCALS_leg_number_matches(i).distances(j),lat_voc(voc_ind),lon_voc(voc_ind));
%        fprintf(1,'\n');
    end

end

fclose(fid);

