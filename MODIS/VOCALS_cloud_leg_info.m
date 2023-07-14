filedir_vocals = '/home/disk/eos1/d.grosvenor/modis_work/VOCALS/';
filename_vocals = 'VOCALS_C130_IncloudLegs.txt';

% fid_voc=fopen(filename_vocals,'rt');
% fgetl(fid_voc)
% [dat_VOCALS]=fscanf(fid_voc,'%f',[5 Inf]);


[month_voc,day_voc,hours_voc,mins_voc,secs_voc,lat_voc,lon_voc]=read_vocals_info_files_func(filedir_vocals,filename_vocals);



year_voc=2008*ones(size(month_voc));