function Z=ptz(P,T,p_start,z_start)
% work out height from pressure
R_gas=287;
grav=9.8;

dp=-100;

imax=fix((p_start-P(end))./-dp);
z_grid=zeros(1,imax);
p=p_start;
z=z_start;
    z_grid(1)=z;
    p_grid(1)=p;
for i=1:imax
    dz=-dp.*R_gas.*interp1(P,T,p)./p./grav;
    z=z+dz;
    p=p+dp;
    i=i+1;
    z_grid(i)=z;
    p_grid(i)=p;
end

% now interpolate to grid to get Z
Z=interp1(p_grid,z_grid,P);