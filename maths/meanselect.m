function me=meanselect(dat,rule)
%me=meanselect(dat,rule) finds mean of only certain points governed by rule
%rule is a string of the form e.g. 'dat>0' or 'dat==30'
%for 2-d matrices - finds mean for column 1

%ilnb=eval(['find(' rule ');']); %evaluate rule

ilnb=eval(['find(~(' rule '));']); %evaluate opposite of rule

dat2=dat;
dat2(ilnb)=NaN; %label all the points not required for mean as NaN so that the others can be easily counted
inan=~isnan(dat2); %find all points required
nps=sum(inan,2); %count number of points in each row

dat(ilnb)=0; %set all points that don't require to zero

me=sum(dat,2)./nps; %calculate sum and divide by number of points in data required