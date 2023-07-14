function [dat_out] = GLACIES_temp_range_FUNC(dat,temp,T0,T1)

dat_temp = NaN*ones([size(temp)]);
inds = find(temp>T0 & temp<=T1);
dat_temp(inds) = dat(inds);
dat_out = meanNoNan(dat_temp,2); %average over 2nd dimension (height).
dat_out = meanNoNan(dat_out,1); %average over time

