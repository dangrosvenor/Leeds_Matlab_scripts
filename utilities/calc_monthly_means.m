function [me,years] = calc_monthly_means(dat,time)
%Data needs to be in form [lat lon time] - or rather just the time dimension
%needs to be last and it needs to be a 3D array.
siz=size(dat);

%Make sure that the time and resulting M and Y vectors are the correct
%orientation.
if size(time,1)==1
    time = time'
end

[Y,M,D] = datevec(time);
M2 = repmat(M,[1 siz(1) siz(2)]);
M2 = permute(M2,[2 3 1]);
Y2 = repmat(Y,[1 siz(1) siz(2)]);
Y2 = permute(Y2,[2 3 1]);

years=unique(Y);

for iy=1:length(years)
for im=1:12
    i=find(M2 ~= im | Y2 ~= years(iy) );
    dat2 = dat;
    dat2(i)=NaN;
    me{iy}(:,:,im) = meanNoNan(dat2,3);        
end
end
   