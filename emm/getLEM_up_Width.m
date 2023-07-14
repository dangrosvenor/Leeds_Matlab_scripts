wdeb=2;

for k=1:length(GridDan(1).Z)
    w=TwoD(idir).W(k,:);
    [maxw iimax]=max(w);
    ii2=find(w<wdeb);
    idiff=ii2-iimax;
    idiff_gt0=idiff;
    idiff_lt0=idiff;
    idiff_gt0(idiff<0)=1e30;
    idiff_lt0(idiff>0)=-1e30;
    [min_i_diffpos iright]=min(idiff_gt0); %index of closest 2 m/s point to right of max 
    [min_i_diffneg ileft]=max(idiff_lt0); %to left
    width(j).w(k,jj)=GridDan(idir).Y1(ii2(iright))-GridDan(idir).Y1(ii2(ileft));
end
