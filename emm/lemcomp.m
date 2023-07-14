
figure
h=hemm(:,1);
wmaxprof=max(wemm,[],2);
plot(wmaxprof,h,'b');
hold on;

maxprof=max(MaxW(1).w(:,1:18),[],2);
plot(maxprof,GridDan(1).Z/1000,'r');

figure;
[maxwt maxwti]=max(MaxW(1).w(:,1:18),[],1);
maxh=GridDan(1).Z(maxwti)/1000; %height of the max w (thermal centre)
plot(GridDan(1).t(1:18)+3,maxh);
hold on;

tmaxemm=max(wemm,[],1);
for i=1:size(wemm,2)
    ii=find(wemm(:,i)==tmaxemm(i));
    maxhemm(i)=(h(ii(1))+h(ii(end)) )/2;
end

plot(temm(1,:),maxhemm,'r');



