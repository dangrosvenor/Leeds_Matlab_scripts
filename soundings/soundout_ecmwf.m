%prepare pr from ecmwf output combined with LES for writing sounding
%load in first LEM dump first

clear p
iz1=49;
iztop=60;
%izlem=findheight(GridDan(1).Z,z(iz1)); %z=ecmwf height

p(:,1)=GridDan(1).Z+620;
p(:,2)=GridDan(1).PREFN;
T=TempLES(GridDan);
p(:,3)=T-273.15;
p(:,10)=GridDan(1).OLQBAR(:,1);

p(:,8)=GridDan(1).VBAR;
p(:,9)=135;

iend=length(T);


hyd_pressure_height;

p(iend+1:iend+iztop-iz1+1,1)=z(iz1:iztop)';
p(iend+1:iend+iztop-iz1+1,2)=parr(iz1:iztop)*100;
p(iend+1:iend+iztop-iz1+1,3)=Tarr(iz1:iztop)-273.15;
p(iend+1:iend+iztop-iz1+1,10)=ecmwf(1).q(iz1:iztop,1,2);

u=ecmwf(1).u(iz1:iztop,1,2);
v=ecmwf(1).v(iz1:iztop,1,2);
bear=bearing(u,v);
mag=sqrt(u.^2+v.^2);
p(iend+1:iend+iztop-iz1+1,8)=mag;
p(iend+1:iend+iztop-iz1+1,9)=bear;


pr(3).p=p;

%now run writesound