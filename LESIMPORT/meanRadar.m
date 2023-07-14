%calc means and rms of radar stats
 rrates=[.15 .3 .6 1.3 2.7 5.6 12 24 49 100 205]; %rain rates corresponding to radar reflectivities
 reflect=[10 15 20 25 30 35 40 45 50 55 60]; %reflectivites for above
 values=[13.5:1.5:60];

ratevals=interp1(reflect,rrates,values,'spline');
summ=sum(sum(npp,1).*values)/sum(sum(npp));    %mean of all dbZ values at all times
summ2=sum(sum(npp,1).*ratevals)/sum(sum(npp));    %mean of all rate values at all times


for i=1:size(npp,1)
    meany(i)=sum(npp(i,:).*values)/sum(npp(i,:)); %mean of each time (dbZ)
    meany2(i)=sum(npp(i,:).*ratevals)/sum(npp(i,:)); %mean of each time (rates)
end

funct
    