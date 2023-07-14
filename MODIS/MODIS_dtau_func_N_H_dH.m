function [dtau,tauc,cw]=MODIS_dtau_func_N_H_dH(H,Nm3,dH,CTT,k,P,fad)
%W and Nm3 in kg/m2 and /m3

if ~exist('k') | length(k)==0 | isnan(k)==1
    k=0.8;
end
if ~exist('CTT')
    CTT = 278;
end
if ~exist('P')
    P=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
%the pressure dependence is not very strong
end
if ~exist('fad')
    fad = 0.8;
    fad = 1.0;    
end

[tauc,W,cw]=MODIS_tau_func_N_H(H,Nm3,CTT,k,P,fad);
[tau_h,W]=MODIS_tau_func_N_H(H-dH,Nm3,CTT,k,P,fad);

dtau = tauc - tau_h;



