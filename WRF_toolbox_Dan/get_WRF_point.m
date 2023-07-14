function [ydat]=get_WRF_point(nc,ih_wrf,ilat,ilon,var,iloc)


%NOTE that choosing the single values from nc{xxx} didn't work
%so this is a profile is first taken out
switch var
    case 'Pressure'
        p=nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
        ydat = p(ih_wrf);
    case 'Temp'
        T=nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
        ydat = T(ih_wrf) - 273.15;
    case 'Vapour'
        f=1e6*28.97/18;
        rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
        T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
        P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
        qsat = satvappress(T(ih_wrf),'goff','liq',P(ih_wrf),1)/f;
        ydat = rh(ih_wrf)/100 .* qsat;
    case 'Wind'
        u=nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
        v=nc{'VV'}(1,:,ilat(iloc),ilon(iloc));

        ydat= sqrt( u(ih_wrf).^2 + v(ih_wrf).^2 );
    case 'Wind dir'
        u=nc{'UU'}(ih_wrf,ih_wrf,ilat(iloc),ilon(iloc));
        v=nc{'VV'}(ih_wrf,ih_wrf,ilat(iloc),ilon(iloc));


        jnorth = ilat(iloc) + 10;
        lons_north = lon2d.var(jnorth,:);
        [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

        %angle of the local north line relative to the grid
        thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );




        for iuv=1:length(u)

            theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

            if u(iuv)==0 & v(iuv)==0
                ydat(iuv) = 0;
            elseif u(iuv)>=0 & v(iuv)>=0
                ydat(iuv) = theta2;
            elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                ydat(iuv) = 180 + theta2;
            elseif u(iuv)<=0 & v(iuv)<=0
                ydat(iuv) = 180 + theta2;
            elseif u(iuv)<0 & v(iuv)>0
                ydat(iuv) = 360 + theta2; %theta2 is negative
            end

        end

        ydat = ydat' + 180 - thetaN*180/pi; %add 180 to make it the direction wind is coming from
        % take away thetaN to give direction relative to north
        i360 = find(ydat>=360);
        ydat(i360) = ydat(i360) - 360;



    case 'cloud'
        cloud=nc{'QCLOUD'}(time,ih_wrf,ilat(iloc),ilon(iloc));
        cloud=cloud+nc{'QICE'}(time,ih_wrf,ilat(iloc),ilon(iloc));
        cloud=cloud+nc{'QSNOW'}(time,ih_wrf,ilat(iloc),ilon(iloc));
        cloud=cloud+nc{'QGRAUP'}(time,ih_wrf,ilat(iloc),ilon(iloc));
        ydat = cloud;
end