%plot3d - script to plot 3D isosurfaces for of clouds
figure

idir=1;
X=GridDan(idir).X1(2:end-1)/1000;
Y=GridDan(idir).Y1(2:end-1)/1000;
ih=findheight(GridDan(idir).Z+620,20e3);
Z=GridDan(idir).Z(1:ih)/1000;





% FV = ISOSURFACE(X,Y,Z,V,ISOVALUE) computes isosurface geometry for
%     data V at isosurface value ISOVALUE. Arrays (X,Y,Z) specify the points
%     at which the data V is given. The struct FV contains the faces and
%     vertices of the isosurface and can be passed directly to the PATCH
%     command.

%rain
minMR=0.1;
V=wrap3d(1e3*ThreeDDan(idir).Q(2:end-1,2:end-1,1:ih,2));  % g/kg
p=patch(isosurface(X,Y,Z,V,minMR));
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');



% %liq
% %minMR=0.5;
% V=wrap3d(1e3*ThreeDDan(idir).Q(2:end-1,2:end-1,1:ih,1));  % g/kg
% p=patch(isosurface(X,Y,Z,V,minMR));
% set(p, 'FaceColor', [0.3 0.3 1], 'EdgeColor', 'none');



%ice
V=wrap3d(1e3*sum(ThreeDDan(idir).Q(2:end-1,2:end-1,1:ih,5),4));  % g/kg
p=patch(isosurface(X,Y,Z,V,minMR));
set(p, 'FaceColor', [0.5 1 1], 'EdgeColor', 'none');

%graupel
V=wrap3d(1e3*sum(ThreeDDan(idir).Q(2:end-1,2:end-1,1:ih,3),4));  % g/kg
p=patch(isosurface(X,Y,Z,V,minMR));
set(p, 'FaceColor','yellow', 'EdgeColor', 'none');

%snow
V=wrap3d(1e3*sum(ThreeDDan(idir).Q(2:end-1,2:end-1,1:ih,4),4));  % g/kg
p=patch(isosurface(X,Y,Z,V,minMR));
%set(p, 'FaceColor', [0.6 1 1], 'EdgeColor', 'none');
set(p, 'FaceColor','red', 'EdgeColor', 'none');

%daspect([1 1 1])

axis tight
camlight right; camlight left; lighting gouraud
view(3)
axis([[X(1) X(end) Y(1) Y(end)]/1.8 0 20]);

xlabel('X (km)','fontsize',18);
ylabel('Y (km)','fontsize',18);
zlabel('Height (km)','fontsize',18);

set(gca,'fontsize',16);