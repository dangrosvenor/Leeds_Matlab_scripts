function [x_p,y_p] = um_xy_nc(x,y,polelon,polelat)

% um_xy_nc.m
% calculates real latitude and lontitude variables for the UM
% lat = y_p(i,j); lon = x_p(i,j)
% from Jeff Chagnon, 2006
if polelat >= 0
    sin_phi_pole = sin(pi/180*polelat);
    cos_phi_pole = cos(pi/180*polelat);
else
    sin_phi_pole = -sin(pi/180*polelat);
    cos_phi_pole = -cos(pi/180*polelat);
end
Nx = length(x);
Ny = length(y);
Exi=x(:);Eyj=y(:);
for i = 1:Nx
    for j = 1:Ny
        %convert to radians

        E_x=pi/180*Exi(i);
        E_y=pi/180*Eyj(j) ;
        %scale eq long to range from -180 to 180
        if E_x > pi; E_x = E_x - 2*pi; elseif E_x < -pi; E_x = E_x + 2*pi;end

        % Compute latitude using equation (4.7)

        arg=cos_phi_pole*cos(E_x)*cos(E_y) +sin(E_y)*sin_phi_pole ;
        arg=min(arg, 1.0)  ;
        arg=max(arg,-1.0) ;
        a_phi=asin(arg) ;
        y_p(i,j)=180/pi*a_phi  ;

        % Compute longitude using equation (4.8)

        term1 =(cos(E_x)*cos(E_y)*sin_phi_pole -sin(E_y)*cos_phi_pole);
        term2=cos(a_phi);
        if abs(term2) < 1e-5
            a_lambda=0.0 ;
        else
            arg=term1/term2;
            arg=min(arg, 1.0);
            arg=max(arg,-1.0);
            a_lambda=180/pi*acos(arg);
            a_lambda=a_lambda*sign(Exi(i)*pi/180-2*pi);%*sign(a_lambda);
        end
        x_p(i,j)=a_lambda+polelon-180;

    end
end
