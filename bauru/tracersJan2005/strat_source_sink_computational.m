% computational model to solve for gradually (linearly) increasing vapour flu
% into strat

%need estimate of stratospheric mass to calculate dq/dt
%know from global_est_strat_percentage.m that:
%area of globe from -L to +L latitude
clear diff

L=15; %latitude for tropics (degrees)
R=6400e3; %radius of the Earth (m)
A=4*pi*R^2*sin(pi*L/180);
A=4*pi*R^2;
%so can now approximate the volume using the height of the strat (=15e3?) and the density (=0.2 kg/m3?)
%really need just the mass involved in the Brewer Dobson circulation
H=15e3;
rho=0.2;
M=A*H*rho ;

%set up time grid
tyears=45;
y2secs=365*24*3600; %conversion factor from years to seconds
tsecs=tyears*y2secs; %convert from years to seconds
npoints=2000;
time=[0:tsecs/npoints:tsecs];
dt=diff(time(1:2));

%now work out initial flux of water vapour
f=1e6*28.97/18; %conversion factor for kg/kg to ppmv
q0=3.8/f; %mean mixing ratio of air at trop (kg/kg) Dessler 1998 (=3.8 ppmv) (assumed value at start of 50 year period)
q0=4.5/f; %value assumed from Rosenlof 2001 GRL paper for annual increase of 0.045 ppmv and 0.5% per year over 45 year
%       - i.e. 1.45*x=x+(0.045*45). x=4.5 rising to 6.5 pppmv after 45 years.
%NOTE - Holton (1995) quotes an annual average mass flux across 100 hPa of 85e8 kg/s (for L=-15 to +15)
%This is based on Rosenlof and Holton (1993), based on UKMO measurements (radiosondes I think)
dmdt = 85e8;  %total mass of air crossing 100 hPa in a year based on Rosenlof and Holton (1995) 
dmvdt = dmdt * q0; %flux of vapour mass (kg/s)
dqdt0 = dmvdt / M; %rate of change of stratospheric vapour mixing ratio

dqdt0 = dqdt0 * 2; %if assuming that methane oxidation contributes the same source as Brewer


% equation is of form 
% dq/dt = dq/dt(source) - dq/dt(sink) = ( dq/dt0 + (T/50)*dq/dt0*(p-1) ) - aq
%                           dq/dt     = dq/dt0 * ( 1 + T/50 * (p-1) ) - aq
% the bit before the -aq is just a linear increase in dq/dt(source) from dq/dt0 to p*dq/dt0 at T=50 years (e.g. p=1.086 for 8.6% rise)
% as before
a = 1/q0 * dqdt0; %for times two are assuming inital methane source equals vaopur source
per = 8.6; % percentage increase required
per = 22.5;
per = 45;
p = 1 + per/100;


add=0.045/2/f/y2secs; %additional term for dq/dt to represent the increase due to methane increases of 0.045/2 ppmv per year
%add=0;
dqdt_noadd=zeros(size(time));
dqdt_sink=zeros(size(time));
delay=8; %delay in years


% now solve with Euler method
q(1)=q0;
for t=2:length(time)
    tdelay = round ( t - delay*y2secs/dt ) - 1;%time index for delay years ago
    tdelay = max([1 tdelay]); %make sure not less than one
    td(t)=tdelay;
%    dqdt_noadd(t) = dqdt0 * ( 1 + time(t)/(45*y2secs) * (p-1) ) ;
     dqdt_noadd(t) = dqdt0*p;
     dqdt_sink(t) = a*q(t-1) ;   %original sink before put in tdelay stuff
%    dqdt_sink = a*q(tdelay) ;
%    dqdt_sink = a*q(1)*( 1 + 0.45*tdelay/length(time) ) ;
%   dqdt_sink = dqdt0;
%    add = dqdt0 * p;

    q(t) = q(t-1) + (   add +  dqdt_noadd(t) - dqdt_sink(t) )*dt; %solve analytical solution for each t and convert to ppmv
%    reldiff(t)= ( (dqdt_noadd(t) - dqdt_sink(t) ) ) /add;
    
    if tdelay==1
%        q(t) = q(1)*(1 + 0.45*t*dt/45/y2secs);
    end
    
end
q=q*f; %convert to ppmv

%figure;
plot(time/y2secs,q,'kx-');






