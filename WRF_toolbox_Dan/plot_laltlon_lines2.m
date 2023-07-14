hold on
[lines_con_a, lines_con_b]=contour(timesTH(1).t(1:pend),zz(1).z,lon2d.var,[-55 -60 -65 -70 -75],'color','w');
clabel(lines_con_a, lines_con_b,'labelspacing',200,'color','w'); %default 144

[lines_con_a, lines_con_b]=contour(timesTH(1).t(1:pend),zz(1).z,lat2d.var,[-55 -60 -65 -70 -75],'color','w');
clabel(lines_con_a, lines_con_b,'labelspacing',200,'color','w'); %default 144