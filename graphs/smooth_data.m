function [xdat,ydat]=smooth_data(xdat,ydat,nfilter)

bfilter=ones([1 nfilter])*1/nfilter;
ydat=filter(bfilter,1,ydat);
ydat(end-nfilter+1:end)=[]; 
xdat(end-nfilter+1:end)=[];
ydat(1:nfilter)=[];
xdat(1:nfilter)=[];