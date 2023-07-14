function [P,h]=hyd_solve(Tprof,Pprof,PSPAN)
%[P,h]=hyd_solve(Tprof,Pprof,PSPAN)
%this solves the hydrostatic equation for height as a function
%of pressure over the pressure interval defined by PSPAN
%when have profiles as functions of pressure (Tprof vs Pprof)

%[P,h] = ODE45(@hydrostatic2,PSPAN,0,[],Pprof,Tprof); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		

%[P,h] = ODE45_dan(@hydrostatic2,PSPAN(1),PSPAN(2),0,[],Pprof,Tprof); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
   

%[P,h] = ODE45_dan(@ode_test,PSPAN(1),PSPAN(2),0); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		


% C-mex version of ODE solver
%[a,b]=ode45_dan_hyd('hydrostatic2',100e2,1e5,17000,1e-3,PP,TT)
[P,h] = ode45_dan_hyd('hydrostatic3',-PSPAN(1),-PSPAN(end),0,1e-3,Pprof,Tprof);
P=-P; h=-h;  %doing this because the C_mex function needs the array in increasing order
%so used the limits -PSPAN(1) and -PSPAN(end) (hydrostatic3 makes negative
%p and z). So need to do the negative at the end
