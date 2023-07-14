z=[-0.5:0.01:0.5]*1000;
dZ=z-z(1);
pot0=295;
%assume constant lapse rate temperature starting at pot0
T0=-11;
Tend=-18;
T=273.15+[T0 : (Tend-T0)/(length(z)-1) : Tend];
press0=((273+T0)/pot0).^(1/0.286) * 1000e2;

ZSPAN=[z(1) z(end)];
[h,P] = ODE45(@hydrostatic,ZSPAN,press0,[],z,T); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H		
T_h = interp1(z,T,h);

pot=T_h.*(1000e2./P).^0.286;
%N=0.01; %assume constant stratifcation
N=sqrt( 9.81/pot(1) * (pot(end)-pot(1))/(h(end)-h(1)) )


x=[0:0.01:20];
y=1000*sin(x/3)/4;

potx0=interp1(h,pot,y(1)); %starting potemp
%now work out temp at each displacement in y
py = interp1(h,P,y); %pressure for all displacements
Ty = potx0 ./ (1000e2./py).^0.286 - 273.15;

%can create different sin waves and then interpolate onto a regular grid

y_upper=150+1*y;
Ty_upper=Twave(y_upper,h,pot,P);

y_lower=-150+1*y;
potx0=interp1(h,pot,y_lower(1)); %starting potemp
%now work out temp at each displacement in y
Ty_lower=Twave(y_lower,h,pot,P);

y_upper2=200+1*y;
Ty_upper2=Twave(y_upper2,h,pot,P);

y_upper3=220+1*y;
Ty_upper3=Twave(y_upper3,h,pot,P);


YY=[y_upper3; y_upper2; y_upper; y; y_lower];
TT=[Ty_upper3; Ty_upper2; Ty_upper; Ty; Ty_lower];

for i=1:length(x)
    T_reg(:,i)=interp1(YY(:,i),TT(:,i),z);
end
% 
% ZZ=100 + 0.35*y_upper; %create a shallow sin path through for the aircraft
% T_path=interp2(x,z,T_reg,x,ZZ); %interpolate from the 2D field along the path

y_upper3_old=220+0.3*y;

ZZ=-150 + y_upper3_old; %create a shallow sin path through for the aircraft
T_path=interp2(x,z,T_reg,x,ZZ);

ZZ2=-200 + y_upper3_old; %create a shallow sin path through for the aircraft
T_path2=interp2(x,z,T_reg,x,ZZ2);

%ZZ2=100 + 0.25*y_upper; %create a shallow sin path through for the aircraft
%T_path2=interp2(x,z,T_reg,x,ZZ2); %interpolate from the 2D field along the path

%have chosen ZZ3 to be of the same amplitude as T_upper3 as this
%had an amplitude of 1.5 degrees, similar to what the temperature
%oscillation would be for air that moved at the same amplitude as that
%which the aircraft moved (based on the pressure of the aircraft)
ZZ3=-100 + y_upper3_old; %create a shallow sin path through for the aircraft
T_path3=interp2(x,z,T_reg,x,ZZ3); %interpolate from the 2D field along the path

%same, but higher up
ZZ4=-50 + y_upper3_old; %create a shallow sin path through for the aircraft
T_path4=interp2(x,z,T_reg,x,ZZ4); %interpolate from the 2D field along the path


figure
pcolor(x,z,T_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
plot(x,y_lower);
plot(x,y);
plot(x,y_upper);
plot(x,y_upper2);
plot(x,y_upper3);
title('Temperature (^{o}C) field of gravity waves with aircraft paths marked');
xlabel('X (km)');
ylabel('Z (m)');

figure
plot(x,T_path,'b'); set(gca,'ydir','reverse'); hold on
plot(x,T_path2,'r');
plot(x,T_path3,'r--');
plot(x,T_path4,'g--');
ylabel('T ^{o}C');
xlabel('X (km)');
title('Temperature change along aircraft paths')


