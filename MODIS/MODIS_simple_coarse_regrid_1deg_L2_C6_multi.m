%Had a problem with this script after doing a certain number of files (31)
%- had to reset Matlab - not sure why...?? Think was because was not
%closing the HDF files - this is done now, so should be ok.

thresh_CTH=[0 3.2];
thresh_SZA=[0 65];

%include the / at the end of filedir
%filedir = '/home/disk/eos15/d.grosvenor/eos8/MOD_L2/NAtlantic/'; %Aqua
filedir = '/home/disk/eos10/d.grosvenor/MOD_L2/NAtlantic/01Aug2016/'; %Terra
fnames = dir([filedir '*.hdf']);

start_file = 1; %in case want to skip some of the files, etc.
end_file = length(fnames);
files = [start_file:end_file];
%files=[1];

for ifname=files

    file_name_h5 = fnames(ifname).name;
    
    
    %% run script - loads the file, filters for CTH and SZA, coarse grains
    %% to around 100km and regrids to 1x1deg regular grid
    MODIS_simple_coarse_regrid_1deg_L2_C6
    
    savename = [filedir '1deg/' file_name_h5 '_1deg.mat'];
    save(savename,'Nd37_1deg','MLAT','MLON','-V7.3');
    
    
end


%	combine_overlapping_data_with_Nans.m - can use this to combine the data
%	from several days into one map (averaging when have more than one
%	data point).
