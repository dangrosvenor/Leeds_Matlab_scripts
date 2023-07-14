function [dir]=wind_dir(u,v,lat2d,lon2d,ilat,ilon,iloc)

jnorth = ilat(iloc) + 10;
lons_north = lon2d.var(jnorth,:);
[temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

%angle of the local north line relative to the grid
thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );




for iuv=1:length(u)

    theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

    if u(iuv)==0 & v(iuv)==0
        dir(iuv) = 0;
    elseif u(iuv)>=0 & v(iuv)>=0
        dir(iuv) = theta2;
    elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
        dir(iuv) = 180 + theta2;
    elseif u(iuv)<=0 & v(iuv)<=0
        dir(iuv) = 180 + theta2;
    elseif u(iuv)<0 & v(iuv)>0
        dir(iuv) = 360 + theta2; %theta2 is negative
    end

end

dir = dir' + 180 - thetaN*180/pi; %add 180 to make it the direction wind is coming from
% take away thetaN to give direction relative to north
i360 = find(dir>=360);
dir(i360) = dir(i360) - 360;