ihm=3;
thresh=0.1e-3;

for k=1:length(GridDan(idir).Z)
    row=sum(TwoDDan(idir).Q(k,:,ihm),3);
	ith=find(row>=thresh);
    me(k)=mean(row(ith));
end
