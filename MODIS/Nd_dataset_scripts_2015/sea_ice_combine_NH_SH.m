
%Do the time matching Load these to do the time matching (already loaded in
%MODIS* script)
[seaice_time3_NH,seaice_max_time3_NH] = seaice_match_times_FUNC(seaiceNH,Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3);
[seaice_time3_SH,seaice_max_time3_SH] = seaice_match_times_FUNC(seaiceSH,Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3);
%So will now have NaN values in both SH and NH files. Will set the NaNs to zero for
%now and then make NaN the places where both were NaN.
inanNH = isnan(seaice_time3_NH); seaice_time3_NH(inanNH)=0;
inanSH = isnan(seaice_time3_SH); seaice_time3_SH(inanSH)=0;
inanNH_max = isnan(seaice_max_time3_NH); seaice_max_time3_NH(inanNH_max)=0;
inanSH_max = isnan(seaice_max_time3_SH); seaice_max_time3_SH(inanSH_max)=0;

iboth = inanNH & inanSH; %When both are NaN - need to use the logical since they are not just integers
iboth_max = inanNH_max & inanSH_max; %When both are NaN - need to use the logical since they are not just integers

seaice_time3 = seaice_time3_NH + seaice_time3_SH;
seaice_time3(iboth) = NaN;

seaice_max_time3 = seaice_max_time3_NH + seaice_max_time3_SH;
seaice_max_time3(iboth_max) = NaN;
