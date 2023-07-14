function [temm,Temm,arr]=readEMMf(fdir,fname,col);
% [temm,Temm,arr]=readEMMf(fdir,fname,col);
% col is the column in which the time is stored
iwczt=dlmread([fdir fname],' ');
dt=iwczt(2:end,col)-iwczt(1:end-1,col);
a=find(dt>0);
da=a(2:end)-a(1:end-1);
nt=size(iwczt,1);
aim=(length(a)+1)*da(1);
if (nt<aim)
    iwczt(nt+1:aim,:)=NaN;
end

iwczt(aim+1:nt,:)=[];

arr=reshape(iwczt,[da(1) length(a)+1 size(iwczt,2)]);
Temm=arr(:,1,2);
temm=arr(1,:,3);