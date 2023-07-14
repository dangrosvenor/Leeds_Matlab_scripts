% cycloid3d_der.m 
function derivs = cycloid3d_der( t, r, flag, qm,B,E )
% Cycloid3d_der: returns the derivatives needed by cycloid3d 
% for a charged particle in an electromagnetic field
% arrays are used as follows:
% B=[Bx;By;Bz], E=[Ex;Ey;Ez] - magnetic, electric fields
% r(1)=x, r(2)=dx/dt=vx
% r(3)=y, r(4)=dy/dt=vy
% r(5)=z, r(6)=dz/dt=vz
% qm - charge to mass ratio
%derivs calculates the derivatives of r. For example
%dx/dt=vx->r(2), dvx/dt=q*(vy*Bz-vz*By+Ex)/m->dr(2)/dt
derivs = [ r(2); qm*(r(4)*B(3)-r(6)*B(2)+E(1)); ...
           r(4); qm*(r(6)*B(1)-r(2)*B(3)+E(2)); ...
           r(6); qm*(r(2)*B(2)-r(4)*B(1)+E(3))];
           
