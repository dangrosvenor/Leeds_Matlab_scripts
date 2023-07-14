function [tau,W,cw]=MODIS_tau_func_N_H(H,Nm3,CTT,k,P,fad)
%W and Nm3 in kg/m2 and /m3

if ~exist('k')
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

Q=2;
rhow=1000;


Parr=ones(size(CTT))*P; %
[cw]=adlwcgm2_just_robs(CTT,Parr) * 1e-3; 
%H = (2*W./cw).^0.5; %cloud depth
W = 0.5.*fad.*cw.*H.^2;

B = 2*sqrt(10)*rhow.^2*fad*cw.^(0.5).*(5/9).^(5/2) ./ (k.*pi.*Q.^3);
tau = (Nm3.*W.^(5/2)./B).^(1/3);