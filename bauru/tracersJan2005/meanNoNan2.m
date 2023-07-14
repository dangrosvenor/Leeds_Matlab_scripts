function me=meanNoNan2(dat,ind)
ndims=length(size(dat));

if ndims==1
    ind=1;
    mean_ind=1;
else
    mean_ind=2;
end

x=permute(dat,[ndims-ind+1 ind]); %reoder matrix so dim dimension is the first

for i=1:size(x,1)
    inan=~isnan(x(i,:));
    me(i)=mean(x(i,inan),mean_ind); %multiplication by inan removes all NaN numbers from sum
end

