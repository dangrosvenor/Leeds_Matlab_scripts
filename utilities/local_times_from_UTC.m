function [LT] = local_times_from_UTC(time_utc,lon_2d)

if maxALL(lon_2d>181)
    lon_2d(lon_2d>180) = lon2d(lon_2d>180) - 360;
end

time_ALL = repmat(time_utc,[1 size(lon_2d,1) size(lon_2d,2)]);
lon_ALL = repmat(lon_2d,[1 1 length(time_utc)]);
%lon_ALL = permute(lon_ALL,[3 1 2]);
time_ALL = permute(time_ALL,[2 3 1]);


LT = mod(time_ALL + lon_ALL/15 , 24); %limit to less than 24 hours.
