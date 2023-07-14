readsound;

degs=1;

t1=pr(2).p(1,3);
th1=theta(t1,pr(2).p(1,2))+degs;

t=th1./(1000./pr(2).p(:,2)).^0.286;

diff=t-pr(2).p(:,3)-273.15;

a=find(diff>=0);
pr(3).p=pr(2).p;
pr(3).p(1:a(end),3)=t(1:a(end))-273.15;