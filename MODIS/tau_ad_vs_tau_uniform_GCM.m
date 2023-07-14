%If assume W, H and N about a GCM constant level cloud:-
W=100e-3; %kg/m2
H=100;
N=100e6; %m^-3
k=0.8;
CTT = 273+10;

%Can calc LWC = W/H
L = W./H;

% LWC = N *4/3 *pi*R^3 * rho_w
% So R = re*k^(1/3) = (3.*L ./ (4*N*pi*1000) ).^(1/3)
rec = ( 3.*L ./ (4*k*N*pi*1e3) ).^(1/3);

%Formula for tau for a constant LWC cloud :-
tau_c = 3/2 .*W ./(1e3.*rec);

% P=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
% %the pressure dependence is not very strong
% Parr=ones(size(CTT))*P;
% [cw]=adlwcgm2_just_robs(CTT,Parr) * 1e-3;

re_ad = re_from_N_LWP(N,W,k,CTT);

%Had = sqrt( 2*W./cw );

tau_ad = 9/5 .*W ./(1e3*re_ad); %Works out as giving a higher tau for the adiabatic model - so for a given GCM model LWP
%the real cloud with the equivalent LWP would be brighter. So would mean
%that would have to tune the model LWP upwards to match the albedo of a
%real cloud.

tau_ad./tau_c