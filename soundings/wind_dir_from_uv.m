function bearing=wind_dir_from_uv(u,v,ilat,ilon)
%%% calculate wind direction from u and v


    bearing = 180/pi * atan ( u ./ v );

    bearing( u==0 & v==0 ) = 0;
    
    i=find(u>0 & v<0);
        bearing(i) = 180 + bearing(i);
    
    i=find(u<=0 & v<=0);
        bearing(i) = 180 + bearing(i);
        
    i=find(u<0 & v>0);
        bearing(i) = 360 + bearing(i); 
        
    bearing=bearing+180;  %add 180 to make the direction the wind comes FROM
    i = find(bearing>=360);
    bearing(i) = bearing(i) - 360;  %keep within the 0-360 range
    
    

    
    
%     if u(iuv)==0 & v(iuv)==0
%         ydat(idat).y(iuv) = 0;
%     elseif u(iuv)>=0 & v(iuv)>=0
%         ydat(idat).y(iuv) = theta2;
%     elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
%         ydat(idat).y(iuv) = 180 + theta2;
%     elseif u(iuv)<=0 & v(iuv)<=0
%         ydat(idat).y(iuv) = 180 + theta2;
%     elseif u(iuv)<0 & v(iuv)>0
%         ydat(idat).y(iuv) = 360 + theta2; %theta2 is negative
%     end

end