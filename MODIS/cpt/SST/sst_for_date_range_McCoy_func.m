function [sst,dates,diff_days] = sst_for_date_range_McCoy_func(sst_in,sst_time,date_range,flag)
%Matches the dates supplied to those from Dan McCoy's SST file

if ~exist('flag')
    flag='';
end

if length(date_range)==1
    date_range(2) = date_range(1);
end

%Load data - note, there is also an SST2.mat file - not sure what this is?
load('/home/disk/eos7/dtmccoy/research/SST/SST.mat')

%The Matlab file SST.mat contains the read in version with an SST array that is 360x180x1253. The 1253 is 1253/52 = 24.0962 years, which makes sense since is labelled as 1990 to present (i.e. 24 years if present is 2014).
%This is the same size as for the .nc file (sst.wkmean.1990-present.nc). Time convention is 
%  days since 1800-1-1 00:00:00
% lat and lon are 1x1 deg.

switch flag
    case 'nearest'
        [diff_days(1),it_01] = min( abs(sst_time-date_range(1)) );
        [diff_days(2),it_02] = min( abs(sst_time-date_range(2)) );      
        it=[it_01:it_02];
    otherwise
        it = find(sst_time >= date_range(1) & sst_time <= date_range(2) );
        diff_days(1:2)=0;
end

if length(it)>0

    sst = sst_in(it,:,:);
    dates = sst_time(it);
    
else
    
    sst = NaN;
    dates = NaN;
    diff_days = NaN;

end

