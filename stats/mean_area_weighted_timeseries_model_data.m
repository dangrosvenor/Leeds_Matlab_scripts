function [me]=mean_area_weighted_timeseries_model_data(dat,var,lats,lons)


i = find(dat.gcm_Plat2D_UM >= lats(1) & dat.gcm_Plat2D_UM < lats(2) & ...
    dat.gcm_Plon2D_UM >= lons(1) & dat.gcm_Plon2D_UM < lons(2));


eval_str = ['siz = size(dat.' var ');']; eval(eval_str);

for it=1:siz(1)
    eval_str = ['dat2 = dat.' var '(it,:,:);']; eval(eval_str);
    me(it) = meanNoNan(dat2(i),1,'',1,1,dat.gcm_area_UM(i));
    
end

