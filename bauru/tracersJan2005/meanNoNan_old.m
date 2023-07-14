function me=meanNoNan_old(dat)

inan=~isnan(dat);
nnums=sum(inan,2); %number of non-NaN numbers in each row
ib=find(inan==0);
dat(ib)=0;

me=sum(dat,2)./nnums; %multiplication by inan removes all NaN numbers from sum

