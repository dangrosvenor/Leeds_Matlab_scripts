%what we really want is to know the path that air would take given an
%initial bearing and windspeed. However, we can assume the constant bearing
%or windspeed for long anyway due to likely pressure and Coriolis forces.
%So is proabably best just to use orthogonal lat-lon lines locally to
%calculate a change in lat and lon proportional to the wind speed in each
%directoin

R=6378.140e3;


dx=distlatlon(30,30,30,31);



[X,Y,Z] = SPH2CART(30*pi/180,30*pi./180,R);
[X2,Y2,Z2] = SPH2CART(30*pi/180,31*pi./180,R);

D=sqrt((X2-X).^2 + (Y2-Y).^2 + (Z2-Z).^2)/1e3

[TH,PHi,R]=CART2SPH(X+100e3,Y,Z);
TH*180/pi
PHi*180/pi


%do on a lambertian projection?
h=m_plot([-163 -162],[70 70]);
xd=get(h,'Xdata')
yd=get(h,'Ydata')
sc=m_scale;

d=sc*sqrt( (xd(1)-xd(2)).^2 + (yd(1)-yd(2)).^2 )/1e3
distlatlon(70,-163,70,-162)

%albers is equal-area
%gnomonic means that straight lines are great circles
%m_ll2xy converts from lon,lat to x,y
m_proj('gnomonic','lon',[-160],'lat',70,'rad',10);
m_grid
[xd,yd]=m_ll2xy([-163 -162],[70 70]);
%sc=m_scale;
%xd=xd*sc; yd=yd*sc;
d=sqrt( (xd(1)-xd(2)).^2 + (yd(1)-yd(2)).^2 )/1e3
d2=distlatlon(70,-163,70,-162)
sc2=d2/d;


latA1=73;
latA2=74;
lonB1 = -170;
lonB2 = -171;
[xd,yd]=m_ll2xy([lonB1 lonB2],[latA1 latA2]);
xd=xd*sc2; yd=yd*sc2;
d=sqrt( (xd(1)-xd(2)).^2 + (yd(1)-yd(2)).^2 )/1e3
d2=distlatlon(latA1,lonB1,latA2,lonB2)


%N.B. distlatlon and lldist use Haversine formula - this is the great
%circle distance!

LATS=[66:74];
LONS=[-170:-150];
[lat2d,lon2d]=meshgrid(LATS,LONS);
[lon2d2,lat2d2]=meshgrid(LONS,LATS);
[X,Y]=m_ll2xy(lon2d2,lat2d2);





