

seaice_dir = '/home/disk/eos8/d.grosvenor/sea_ice_data/';

seaice_fileload_1deg_NH = [seaice_dir 'north/saved_seaice_1deg_2000-2015_all_NHemisphere_20161223T184753.mat'];
seaiceNH = load(seaice_fileload_1deg_NH);
seaice_fileload_1deg_SH = [seaice_dir 'south/saved_seaice_1deg_2000-2015_all_SHemisphere_20170110T094057.mat'];
seaiceSH = load(seaice_fileload_1deg_SH);

[Y,M,D]=datevec(seaiceNH.seaice_array_1deg_datenum);
[Y2,M2,D2]=datevec(seaiceSH.seaice_array_1deg_datenum);

years = unique([Y;Y2]);



for iy=1:length(years)
    year = years(iy);
    year_str = num2str(year);

    tinds = find(Y==year);
    seaice_array_1deg = seaiceNH.seaice_array_1deg(:,:,tinds);
    seaice_array_1deg_max = seaiceNH.seaice_array_1deg_max(:,:,tinds);
    seaice_array_1deg_datenum = seaiceNH.seaice_array_1deg_datenum(tinds);      
    Notes = ['Separated into individual years  from ' seaice_fileload_1deg_NH ' using sea_ice_separate_years.m'];
    savename = [seaice_dir 'north/saved_seaice_1deg_' year_str '_all_NHemisphere_' datestr(now,29)];
    save(savename,'seaice_array_1deg','seaice_array_1deg_max','seaice_array_1deg_datenum','Notes','-V7.3');
    
    tinds = find(Y2==year);
    seaice_array_1deg = seaiceSH.seaice_array_1deg(:,:,tinds);
    seaice_array_1deg_max = seaiceSH.seaice_array_1deg_max(:,:,tinds);
    seaice_array_1deg_datenum = seaiceSH.seaice_array_1deg_datenum(tinds);     
    Notes = ['Separated into individual years  from ' seaice_fileload_1deg_SH ' using sea_ice_separate_years.m'];
    savename = [seaice_dir 'south/saved_seaice_1deg_' year_str '_all_SHemisphere_' datestr(now,29)];
    save(savename,'seaice_array_1deg','seaice_array_1deg_max','seaice_array_1deg_datenum','Notes','-V7.3');

end