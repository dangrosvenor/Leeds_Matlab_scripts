filedir = '/home/disk/eos8/d.grosvenor/MAS_data/temp/';

file_name_h5 = 'ace_mas_980530.017';

filename_h5=[filedir file_name_h5];


INFO = hdfinfo(filename_h5);






SD_id = hdfsd('start',filename_h5,'read'); %open the file


[lat,dimsizes_lat,sampling_along,sampling_across]=get_hdf_data_dan('Latitude',SD_id,INFO);