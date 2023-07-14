function [LS_rain_bins] = UM_Iceland_rain_bins_FUNC(LS_rain_PI_ALL,bin_vals,min_val,nbins,deltalog)

itidy_bins=0;
switch bin_vals    %Note - the values in the LS_rain_PI_ALL array are in kg/m2/s (mm/s) not mm/hr
    case 'choose'

    case 'Log10 Percentiles' %Works better in log space
        %min_val=0.05;        
        %min_val=0.005; 
        %nbins=7;
        dprc=100/nbins;
        prcs_vals=[0 dprc:dprc:100];
        lwp_for_prc = LS_rain_PI_ALL(:); %convert to mm/hr     
        lwp_for_prc(lwp_for_prc<min_val)=NaN;
        lwp_for_prc = log10(lwp_for_prc); %do prctiles in log space - to ensure enough points in each bin 
        prcs = prctile(lwp_for_prc,prcs_vals);        
        dlog_bin = prcs(2) - prcs(1); %The first bin will contain all those with values less than min_val - will make the
        %(arbitary) lower boundary of this bin equal to the upper bound minus the
        %width of the second bin to make the 2d histo plot look nicer.
        min_bin = prcs(1) - dlog_bin; %min_val/5;
        LS_rain_bins = 10.^[min_bin prcs]; 
        
        
    case 'log10'
        %deltalog=0.3;
        rmax = maxALL(LS_rain_PI_ALL);
        maxlog = ceil(log10(rmax)/deltalog)*deltalog;
        logmin=log10(min_val);
        LS_rain_bins = 10.^[logmin:deltalog:maxlog]; %Start with say around 0.01 mm/hr up to max
        itidy_bins=1;
end


X_driver = LS_rain_PI_ALL; %convert to mm/hr
X_driver(X_driver<LS_rain_bins(1)) = LS_rain_bins(1)*1.01; %Force values less than a min bin to be = to the min bin so can have a log x scale.
%Y_driver = LWP_PD_ALL - LWP_PI_ALL;

clear Y_me Y_N Y_std
for i=1:length(LS_rain_bins)-1
   ime = find(X_driver>=LS_rain_bins(i) & X_driver<LS_rain_bins(i+1));
   %length(ime)
   %[Y_me(i),Y_N(i),Y_std(i)] = meanNoNan(Y_driver(ime),1); 
   Y_N(i) = length(ime);
end
if itidy_bins==1
    inan=find(Y_N == 0);
    Y_N(inan)=[]; LS_rain_bins(inan)=[];
    Ycum_rev = flipdim(cumsum(flipdim(Y_N,2)) / sum(Y_N) , 2);
    min_prc=0.03;
    imin = find(Ycum_rev<min_prc);
    LS_rain_bins(imin+1)=LS_rain_bins(end);
    LS_rain_bins(imin+2:end)=[];
end

