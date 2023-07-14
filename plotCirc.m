th=1:360;
r=240;

yd=r.*sin(th.*pi./180);
xd=r.*cos(th.*pi./180);

ix=63+ceil(xd./10);
iy=38+ceil(yd./10);

hold on;

plot(ix,iy,'kx');