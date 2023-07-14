%save some POLDER fields for Dan McCoy

savefile_POLDER = '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_Reff_2007-2008.mat';

save(savefile_POLDER,'daymeanALL_Par2_CDR','daynum_timeseries3_POLDER','modisyear_timeseries3_POLDER','-V7.3');
