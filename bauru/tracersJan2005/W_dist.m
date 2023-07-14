clear f f2

step=2; %bin size m/s

ww=w_prctiles(1).w(2:21,:,:);

wb=[minALL(ww):step:ceil(maxALL(ww))];
%dw=wb(2:end)-wb(1:end-1);


ntot=size(ww,2)*size(ww,3); %total number of points in population, each is 5% of that point in terms of number

for i=2:length(wb)
    a=find(ww>=wb(i-1) & ww<wb(i));
    f(i)=length(a)/20/ntot/step; %each point in time-height space is 5%/ntot of the total number of w values (20 5th percentiles)
end

figure
plot(wb,log10(f));
hold on;

ww=w_prctiles(2).w(2:21,:,:);

wb=[minALL(ww):step:ceil(maxALL(ww))];
%dw=wb(2:end)-wb(1:end-1);

ntot=size(ww,2)*size(ww,3); %total number of points in population, each is 5% of that point in terms of number

for i=2:length(wb)
    a=find(ww>=wb(i-1) & ww<wb(i));
    f2(i)=length(a)/20/ntot/step; %each point in time-height space is 5%/ntot of the total number of w values (20 5th percentiles)
end


plot(wb,log10(f2),'r');
'done'