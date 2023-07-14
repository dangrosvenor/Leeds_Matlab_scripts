D = 200; % length of dehydrated air in 2-d case 
D = 1150;

L=1000; %length of dehyd zone (km)
dndt=24 / (1000e3)^2;  %rate of events per km2 / month
A = pi * D^2 / 4;   % area of dehydrated air produced by model
V = 5 * 1e3 / (3600*24*30);  % horiz processing speed
qdef=1;
qdef=0.5;

dQ = dndt * A * qdef * L / V  %calculation of total dq change for air parcel traversing dehydration zone
                              %(assuming no uplift) of length L at speed V km/month
                              
% = 0.4 ppmv for qdef=1, L=1000, v=5 m/s, D=200 km, dndt = 24 /month over 1000x1000 km                              
%derived from:
% mean reduction of mixing ratio over an area A2 of MR qenv due to the presence of dehydrated air of 
%q ppmv covering area A is:
% dq = (A*q + (A2 - A)*qenv ) / A2  -  qenv
% dq = A*(q - qenv)/A2   define q - qenv = qdef

% total number of events per month over A2 area = dndt * A2 where dndt is rate of events occuring per unit area
% then rate over area A2 dqdt = A*qdef*dndt
% so total change in q of air experiencing this for time of T = L/V is
% dQ = A*qdef*dndt * L/V

%should really apply this interactively so that the proper qdef is calculated from evolving qenv
% *** no need to consider uplift !! *** but obviously only going to affect certain height range 
% so uplift will need to be considered after air has passed through the region
% e.g. if work out that 1 km layer reduced by 2 ppmv for reasonable passage time then
% uplift rate will determine how often air will need to experience such a processing through one of these
% regions before moist air starts to follow it into stratosphere.
