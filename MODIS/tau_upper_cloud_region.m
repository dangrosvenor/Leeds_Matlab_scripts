function tau = tau_upper_cloud_region(H,N,z_lwc_max,f_ad,CTT,k,Q)

if ~exist('k')
    k=0.88;
end
if ~exist('Q')
    Q=2;
end



%Cloud top pressure
Parr=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
%the pressure dependence is not very strong

cw = adlwcgm2_just_robs(CTT,Parr) * 1e-3;

Lmax = f_ad.*cw.*H*z_lwc_max; %calculate the max LWC

%% Now ramp down LWC and Nd linearlly from Lmax to zero from z_lwc_max to H
dz = 0.05; %metres
% z goes from top downwards; z=0 is the top
z = [0:dz:300*(1-z_lwc_max)];
nz=length(z);
N2 = [0:N/(nz-1):N];
L2 = [0:Lmax/(nz-1):Lmax];

beta = extinction_func_N_L(N2,L2,k,Q); 

tau = dz.* trapz(beta);
    





