function sst_amsre_time32 = amsre_block_av_time_function(sst_amsre_time3,N)
%avearge over N day windows in order to help remove NaNs due to orbit gaps

[s1,s2,s3]=size(sst_amsre_time3);


%sst_amsre_time3 = dat;
clear dat

it_inds = [1:N:s3];

if diff([it_inds(end) s3]) < 4
    it_inds(end) = s3;
end
   
sst_amsre_time32 = NaN*ones([s1 s2 s3]);

for it=1:length(it_inds)-1
    
    it_inds2 = it_inds(it):it_inds(it+1);
    
    val = meanNoNan( sst_amsre_time3(:,:,it_inds2),3 );
    
    for it2=1:length(it_inds2)
        sst_amsre_time32(:,:,it_inds2(it2)) = val;
    end
    
end