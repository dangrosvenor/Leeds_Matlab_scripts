function [out]=calc_stats(prc_thresh,true_data,test_data)



nT=size(true_data,3);

out.aux.prc_thresh = prc_thresh;


abs_bias = test_data - true_data;    
prc_bias = 100*(test_data - true_data) ./true_data;
prc_bias(true_data<prc_thresh)=NaN;



%Make data with consistent NaNs - might want to consider a switch for this?
inan = find(isnan(true_data)==1 | isnan(test_data)==1);

true_data2 = true_data;
test_data2 = test_data;
true_data2(inan)=NaN; test_data2(inan)=NaN;

true_data3 = true_data;
test_data3 = test_data;
true_data3(inan)=[]; test_data3(inan)=[];

%correlation coefficient
if length(true_data3(:))>0 & length(test_data3(:))>0
    out.rsq = corr(true_data3(:),test_data3(:));
else
    out.rsq=NaN;
end

[me_true,N,std_dev] = meanNoNan(true_data2(:),1);
[me_test,N,std_dev] = meanNoNan(test_data2(:),1);
out.me_true = me_true;
out.me_test = me_test;

%Also output the non-NaN-matched means - these match the ones calculated in
%plot_global_maps that appear on the titles of the figures.
[me_true_non_matched,N,std_dev_non_matched] = meanNoNan(true_data(:),1);
[me_test_non_matched,N,std_dev_non_matched] = meanNoNan(test_data(:),1);
out.me_true_non_matched = me_true_non_matched;
out.me_test_non_matched = me_test_non_matched;

% [i_lin,i_lin_edges,abs_bias_regional]=get_lat_lon_irregular_with_time(1,lat_bnds,lon_bnds,lat2d,lon2d,lat2d_edges,lon2d_edges,abs_bias);
% out.regional_mean_bias = meanNoNan(abs_bias_regional(:),1);
% 
% out.regional_mean_prc_bias = 100*out.regional_mean_bias./meanNoNan(true_data(:),1);
% %regional mean of % values
% [i_lin,i_lin_edges,prc_regional]=get_lat_lon_irregular_with_time(1,lat_bnds,lon_bnds,lat2d,lon2d,lat2d_edges,lon2d_edges,prc_bias);
% out.regional_mean_prc_bias_area_av = meanNoNan(prc_regional(:),1);


%RMSE
sq = (true_data2(:)-test_data2(:)).^2;
out.rmse = sqrt ( meanNoNan( sq , 1) );
out.rmse_div_mean = out.rmse ./ me_true;




%% Normalised mean bias factor (BNMBF) - see Gustafson, ASL, 2012
if me_test >= me_true
    out.BNMBF = 100 * (me_test./me_true - 1); %covert to percentage
else
    out.BNMBF = 100* (1 - me_true./me_test);
end

%% Normalised mean absolute error factor (ENMAEF) - see Gustafson, ASL, 2012
% Bit like RMSE
ab_vals = abs(true_data2(:) - test_data2(:));
[sum_ab,N,std_dev] = meanNoNan(ab_vals(:),1,'sum');
[sum_true,N,std_dev] = meanNoNan(true_data2(:),1,'sum');
[sum_test,N,std_dev] = meanNoNan(test_data2(:),1,'sum');

if me_test >= me_true
    out.ENMAEF = 100* (sum_ab ./ sum_true); %covert to percentage
else
    out.ENMAEF = 100* (sum_ab ./ sum_test); %covert to percentage
end


%% Revised normalised mean bias factor (BNMBF) - see Gustafson, ASL, 2012
% revised in the above paper to deal with when have negative data
if abs(me_test) >= abs(me_true) & sign(me_test)==sign(me_true)
    out.BNMBFrev = 100* (abs(me_test./me_true) - 1); %covert to percentage
elseif sign(me_test)==sign(me_true)
    out.BNMBFrev = 100* (1 - abs(me_true./me_test)); %covert to percentage
else
    out.BNMBFrev = NaN; %Not defined if signs are different.
end

%% Revised Normalised mean absolute error factor (EMAEF) - see Gustafson, ASL, 2012
% revised in the above paper to deal with when have negative data
% Bit like RMSE
ab_vals = abs(true_data2(:) - test_data2(:));
EMAGE = meanNoNan(ab_vals(:),1);

if abs(me_test) >= abs(me_true) & sign(me_test)==sign(me_true)
    out.ENMAEFrev = 100 * (EMAGE ./ abs(me_true)); %covert to percentage
elseif sign(me_test)==sign(me_true)
    out.ENMAEFrev = 100 * (EMAGE ./ abs(me_test)); %covert to percentage
else
    out.ENMAEFrev = NaN; %Not defined if signs are different.
end




