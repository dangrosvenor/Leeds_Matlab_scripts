 hold on;plot(-49.03,-22.36,'ko','MarkerSize',8,'Linewidth',2,'markerfacecolor','k'); %plot location of Bauru

lon=ecmwf(1).lon;
lon=ecmwf(1).lon;

ih=1;

fw=2;
quiver(lon-360,lat,squeeze(ecmwf(1).v(it,ih,:,:)),squeeze(ecmwf(1).u(it,ih,:,:)),fw,'k');

