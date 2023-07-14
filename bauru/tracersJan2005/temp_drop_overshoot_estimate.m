ih=findheight(GridDan(1).Z+620,15.75e3);
ih2=findheight(GridDan(1).Z+620,16.05e3);
th=mean ( TwoD.TH1(ih2,:) + GridDan(1).THREF(ih2) , 2);
p=GridDan(1).PREFN(ih2);

T=th * (p/1e5)^(0.286)

dp= - diff(GridDan(1).PREFN([ih2 ih]));

dT= - dp * th * 0.286 * (p/1e5)^(-0.286-1) /1e5 

drho=p.*28.97e-3/8.3144./T/T *dT

