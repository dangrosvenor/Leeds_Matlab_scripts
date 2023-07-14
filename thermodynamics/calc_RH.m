function rh=calc_RH(qv,p,T)
%function rh=calc_RH(qv,p,T)

f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone. Or f=1e6/eps where eps=0.622 - MR = ppmv/f
qs = SatVapPress(T,'goff','liq',p,1) / f;

rh = qv./qs;