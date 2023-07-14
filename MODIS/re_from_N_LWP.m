function re=re_from_N_LWP(N,W,k,CTT)
%function re=re_from_N_LWP(N,W,k,CTT)
%Returns re from Nd and LWP
%N in m^-3, W in kg/m2, CTT in K

cw=MODIS_justN_func(999,999,'calc',0,CTT,'cw');
re= (9/8 * cw .* W .* (1/(pi*k.*N*1e3)).^2 ).^(1/6);

