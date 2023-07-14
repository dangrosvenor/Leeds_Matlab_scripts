ivert_slope=1; %flag to use the vertical slope concept
vert_slope_angle = 68.2; %angle (degrees) of vertical slope 
%vert_slope_angle = 45; %angle (degrees) of vertical slope 

clear W_af U_af

heights_map = [2400:10:4200]; %heights to make a 2D (u,v) map for

for ip=1:length(W2)
    z = Z2(ip);
    zz = heights_map-z;
    
    dx = - zz*tan(vert_slope_angle*pi/180); %if z is positive then the dx will be negative
    %this assumes that the flight level was the same everywhere as at X2(ip)
    %so really should calculate where the line with the given angle intercepts the flight track
    dz = (X2-X2(ip))./tan(vert_slope_angle*pi/180);
    
    for iz=1:length(zz)
        zline = heights_map(iz) + dz;
        [a,b]=min(abs(zline-Z2)); %find the closest point where the line crosses the flight path
        if b==1 & zline(1)>Z2(1) %if have negative z and line is not crossing the flight path
            W_af(ip,iz)=NaN;
            U_af(ip,iz)=NaN;
        elseif b==length(X2) & zline(end)<Z2(end) %if have negative z and line is not crossing the flight path
            W_af(ip,iz)=NaN;
            U_af(ip,iz)=NaN;
        else            
            %W_af(ip,iz) = interp1(X2,W2,dx+X2(ip));            
            %U_af(ip,iz) = interp1(X2,U2,dx+X2(ip)); 
            W_af(ip,iz) = W2(b);
            U_af(ip,iz) = U2(b);            
        end
    end
    
    
end

