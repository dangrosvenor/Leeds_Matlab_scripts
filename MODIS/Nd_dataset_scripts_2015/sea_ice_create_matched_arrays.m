%Loads in the sea-ice data and produces arrays the same size at MODIS data
%and at same times.
%Uses Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3 for
%matching

iload_seaice=1;

if iload_seaice==1

%% Combine NH and SH seaice (is only SH in the file)
%seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2007_all_NHemisphere_20150529T103046.mat';
%seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2008_JJA_all_NHemisphere_20161222T222241.mat';
seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2000-2015_all_NHemisphere_20161223T184753.mat';
seaiceNH = load(seaice_fileload_1deg_NH);
%seaice_fileload_1deg_SH = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2008_JJA_all_SHemisphere_20161221T095015.mat';
seaice_fileload_1deg_SH = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2000-2015_all_SHemisphere_20170110T094057.mat';
seaiceSH = load(seaice_fileload_1deg_SH);

end


[seaice_time3,seaice_max_time3] = sea_ice_create_matched_arrays(Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3,seaiceNH,seaiceSH)
