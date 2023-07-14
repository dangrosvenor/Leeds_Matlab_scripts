% solving equation for stratosphere source sink model
% dq/dt = dq/dt(source) - dq/dt(sink) = pdq/dt0 - aq
% where p is some factor increase of the new source relative to the original one
% (e.g. assuming sudden rise in source of by, say, 8.6 % : then p=1.086)
% are also assuming that the sink is proportional to the amount of q in 
% the stratosphere (perhaps justified if assume that the "exit flow rate" remains
% the same but just q changes. Are also assuming here that the mixing time is very quick throughout
% the stratosphere so that input of q is instantly mixed to give a new mean q.
% Then the ODE becomes dq/dt + aq = c  (where c = pdq/dt0)
% Can solve this by mulitplying each side by exp(at) and noting that
% d/dt(exp(at)q) = exp(at)*(dq/dt + aq) for LHS of equation
% This leads to q = q0(1-p)exp(-at) + pq0, where q0 is the original q value at t=0 and a = 1/q0 * dq/dt0(original input rate)


%need estimate of stratospheric mass to calculate dq/dt
%know from global_est_strat_percentage.m that:
%area of globe from -L to +L latitude
L=15; %latitude for tropics (degrees)
R=6400e3; %radius of the Earth (m)
A=4*pi*R^2*sin(pi*L/180);
%so can now approximate the volume using the height of the strat (=15e3?) and the density (=0.2 kg/m3?)
%really need just the mass involved in the Brewer Dobson circulation
H=15e3;
rho=0.2;
M=A*H*rho;

%now work out flux of water vapour
q0=3.8; %mean mixing ratio of air at trop (ppmv) Dessler 1998 (assumed value at start of 50 year period)
f=1e6*28.97/18; %conversion factor for kg/kg to ppmv
%NOTE - Holton (1995) quotes an annual average mass flux across 100 hPa of 85e8 kg/s (for L=-15 to +15)
%This is based on Rosenlof and Holton (1993), based on UKMO measurements (radiosondes I think)
dmdt = 85e8;  %total mass of air crossing 100 hPa in a year based on Rosenlof and Holton (1995) 
dqdt0 = dmdt * qvap/f; %flux of vapour mass (kg/s)

%now calc values for the ODE
a = 1/q0 * dqdt0;
per = 8.6; % percentage increase required
p = 1 + per/100;

%now solve for a give t
tyears=50;
y2secs=365*24*3600; %conversion factor from years to seconds
tsecs=tyears*y2secs; %convert from years to seconds
npoints=200;
t=[0:tsecs/npoints:tsecs];

q = q0.*(1-p).*exp(-at) + p*q0;

figure;
plot(t/y2secs,q);






