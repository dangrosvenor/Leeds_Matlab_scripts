function dat_out = squeeze_keep1(dat,inotkeep)
if nargin<2
    inotkeep=0; %default is to keep/add a dimension so end up with [1 M N P etc]
end

%squeeze to start with
dat_out = squeeze(dat);


if inotkeep==0
%add one singleton dimension back
   dat_out = shiftdim(dat_out,-1);
end
