function z = hyd_pressure_height(parr,Tarr,z0)
%finds the height for a given pressure and temp profile starting at z0 (m) using hydrostatic equation
%usage: z = hyd_pressure_height(parr,Tarr,z0)

z(1)=z0;
% parr=GridDan(1).PREFN;
% Tarr=tempLES(GridDan);


M = 28*1.67E-27;
k = 1.38E-23;
G = 9.81;

%T(1)=Tarr(1);
%p(1)=parr(1);




for i=2:length(parr)
    rho=parr(i-1)*M/k/Tarr(i-1);
   % dz=zarr(i)-zarr(i-1);
    dp=parr(i)-parr(i-1);
    F=-1/(rho*G);
    z(i)=z(i-1)+F*dp;
  %  T(i)=interp1(parr,Tarr,zarr(i));
end

% clear p2
% p2(1)=parr(1);
% zarr=GridDan(1).Z;

% for i=2:length(parr)
%     rho=p2(i-1)*M/k/Tarr(i-1);
%     dz=zarr(i)-zarr(i-1);
%    % dp=parr(i)-parr(i-1);
%     F=-(rho*G);
%     p2(i)=p2(i-1)+F*dz;
%   %  T(i)=interp1(parr,Tarr,zarr(i));
% end