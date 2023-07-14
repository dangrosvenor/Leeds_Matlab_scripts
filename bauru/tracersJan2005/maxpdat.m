function [max,h1,t1,min,h2,t2]=maxpdat(pdat,zz,time,pend)

[max b]=maxALL(pdat(1).p(:,1:pend));
h1=zz(1).z(b(1))+620;
t1=time(1).t(b(2));

[min b]=minALL(pdat(1).p(:,1:pend));
h2=zz(1).z(b(1))+620;
t2=time(1).t(b(1));