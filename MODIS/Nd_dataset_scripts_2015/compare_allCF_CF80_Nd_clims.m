
tit_str='ANNUAL';
r=(Nd_SZA_lt_70_CF_gt_80_AND_CTT_gt_273_clim_ANNUAL./Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_ANNUAL);
dat=r; qpcolor(dat); caxis([1 1.5]); title([tit_str ' relative difference CF80/allCF max = ' num2str(maxALL(dat))]);

tit_str='DJF';
rDJF=(Nd_SZA_lt_70_CF_gt_80_AND_CTT_gt_273_clim_DJF./Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_DJF);
dat=rJJA; qpcolor(dat); caxis([1 1.5]); title([tit_str ' relative difference CF80/allCF max = ' num2str(maxALL(dat))]);

tit_str='MAM';
rJJA=(Nd_SZA_lt_70_CF_gt_80_AND_CTT_gt_273_clim_MAM./Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_MAM);
dat=rMAM; qpcolor(dat); caxis([1 1.5]); title([tit_str ' relative difference CF80/allCF max = ' num2str(maxALL(dat))]);

tit_str='JJA';
rJJA=(Nd_SZA_lt_70_CF_gt_80_AND_CTT_gt_273_clim_JJA./Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_JJA);
dat=rJJA; qpcolor(dat); caxis([1 1.5]); title([tit_str ' relative difference CF80/allCF max = ' num2str(maxALL(dat))]);

tit_str='SON';
rSON=(Nd_SZA_lt_70_CF_gt_80_AND_CTT_gt_273_clim_SON./Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_SON);
dat=rSON; qpcolor(dat); caxis([1 1.5]); title([tit_str ' relative difference CF80/allCF max = ' num2str(maxALL(dat))]);