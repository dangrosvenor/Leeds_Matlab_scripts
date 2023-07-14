 ztot=Radar(GridDan(idir),TwoD,izmin,izmax);
ztot=10*log10(ztot(:,xinds)); 
[a b]=find(ztot>=10); %all points with >=10 dBz
ihs=unique(a); %all the different height indices with points > 10 dBz

n10dbz(j).n(1:length(Grid.Z),jj)=0;
for iradar=1:length(ihs)
    n10dbz(j).n(ihs,jj)=length(find(a==ihs(iradar))); %number of points > 10 dBz at each height index with > 10 dBz
end