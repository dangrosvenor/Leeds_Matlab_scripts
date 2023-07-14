file_goes_multi_whole_dom = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151117T234455.mat';

load(file_goes_multi_whole_dom,'gcm_Plat2D_GOES','gcm_Plon2D_GOES');

d=diff(gcm_Plat2D_GOES,[],1);
dlat_GOES = meanNoNan(meanNoNan(d,1),1);



d=diff(gcm_Plon2D_GOES,[],2);
dlon_GOES = meanNoNan(meanNoNan(d,1),1);
