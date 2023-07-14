%ho2.m
%Calculation of position, velocity, and acceleration for a body in
%free fall with air resistance versus time. 
%The equations of motion are used for small time intervals
clear;
%NPTS=200;TMAX=20.0;%example Maximum number of points and maximum time
TTL=input(' Enter the title name TTL:','s');%string input
NPTS=input(' Enter the number calculation steps desired NPTS: ');
TMAX=input(' Enter the run time TMAX: ');
NT=NPTS/10;%to print only every NT steps
%G=9.8;M=1.0;C=0.05;y0=0;v0=110;% example Parameters
G=input(' Enter value of gravity G: ');
M=input(' Enter the object mass M: ');
C=input(' Enter the drag coefficient C: ');
y0=input(' Enter the initial height y0: ');  % Initial Conditions
v0=input(' Enter the initial velocity v0: ');% Initial Conditions
FLAG=input(' Enter 0 (v drag) or 1 (v^2 drag) FLAG: ');
t0=0.0;% start at time t=0
dt=TMAX/NPTS;%time step size
if FLAG ==0
   F=-M*G-C*v0;           % initial force - case 1
   vt=abs(M*G/C);         % terminal velocity
   elseif FLAG==1
     F=-M*G-C*v0*abs(v0); % initial force - case 2
     vt=sqrt(M*G/C);      % terminal velocity
end;
%dt,FLAG, and vt used
fprintf(' FLAG=%1i, Time step dt=TMAX/NPTS=%5.2f, vt=%5.2f\n',FLAG,dt,vt);
a0=F/M;% initial acceleration
fprintf('    t       y       v       a\n');%output column labels
v(1)=v0;
y(1)=y0;
a(1)=a0;
t(1)=t0;
fprintf('%7.4f %7.4f %7.4f %7.4f\n',t(1),y(1),v(1),a(1));%print initial values
for i=1:NPTS
    v(i+1)=v(i)+a(i)*dt;                   %new velocity
    y(i+1)=y(i)+v(i+1)*dt;                 %new position
    t(i+1)=t(i)+dt;                        %new time
    if FLAG ==0
       F=-M*G-C*v(i+1);                    %new force - case 1
       elseif FLAG==1
         F=-M*G-C*v(i+1)*abs(v(i+1));      %new force - case 2
    end;
    a(i+1)=F/M;                            %new acceleration
% print only every NT steps
    if(mod(i,NT)==0)
        fprintf('%7.4f %7.4f %7.4f %7.4f\n',t(i+1),y(i+1),v(i+1),a(i+1));
    end;
end;
    plot(t,y,'k-',t,v,'b:',t,a,'r-.');
    ylabel('y (m), v (m/s), a (m/s^2)','FontSize',14);
    xlabel('time','FontSize',14);
    title(TTL,'FontSize',14);
    h=legend('position','velocity','acceleration',0); set(h,'FontSize',14)
