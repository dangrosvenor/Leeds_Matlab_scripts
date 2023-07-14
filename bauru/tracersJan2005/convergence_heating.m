clear diff

dT=2;
w=5;
H=1.5e3; %height of influx
X=5e3; %distance over which are calculating influx (updraught width?)

%dhdt=w*14e3*2.5e3*1.012e3*dT/A/1e3 %heating rate equivalent of mass flux at 10 m/s of dT temperature difference

dhdt=w*H/X  * 1.012e3 * dT / 1e3

x=[0:0.1:7e3];

x=[0:0.1:7e3];

dx=diff(x);

A2=pi*7e3^2;
sumth=sum( ( cos(pi.*x(2:end)/14e3) ).^2 *2*pi.*x(2:end).*dx)/A2; %3d radially symmetric average
            %=0.298
sumth=sum( ( cos(pi.*x(2:end)/14e3) ).^2 .*dx)/sum(dx); %in 2-d average is integral(f(x) dx) / integral(dx)
            %=0.5
thav=2.1e-2*2.5e3*1.012e3*sumth /1e3 %in kW/m2

dedt_up=w*1*1.012e3*dT /1e3  %w*rho*CdT  using updraught flux to calculate amount of heat being input
                             %into cloud base if assuming that low level convergence supplies
                             %the updraught with air dT K above environment
                             
                             
                             
%if consider updraught flux then dqdt=wdQ/dz = rate at which q is replaced in layer dz deep, dQ=qnew-qold for input of 
%moist air of qnew relaceing air of qold. Continual replacement justified as otherwise "dry" air in rest of domain would flow to 
%replace air moved upwards. so dQ=dqdt * dZ /w
dqdt=5.025e-6; %moistening rate applied in model
w=3;
dz=2.5e3;

dQ=dqdt*dz/w

dQ/17e-3

