function sst_amsre_time32 = amsre_block_av_time_function2(sst_amsre_time3,N,years,months,days)
%avearge over N day windows in order to help remove NaNs due to orbit gaps

[s1,s2,s3]=size(sst_amsre_time3);

time_mat = datenum(years,months,days);

%find when the days we have are continuous
diffs = diff(time_mat);
i=find(diffs>1);
istarts = [1 i+1]; %the starts of all the continuous periods
iends = [i length(time_mat)]; %the ends

sst_amsre_time32 = NaN*ones([s1 s2 s3]);

for it=1:length(istarts)
    inds = istarts(it):iends(it);
    sst_amsre_time32(:,:,inds) = amsre_block_av_time_function(sst_amsre_time3(:,:,inds),N);
end

