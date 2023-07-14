nx=4; %multiplication factor for each dimension
ny=4; %16 is too high! 4 is ok.

x_fine=x_grid(1):DX/nx:x_grid(end);
y_fine=y_grid(1):DY/ny:y_grid(end);

[X,Y]=meshgrid(x_fine,y_fine); %makes X and Y matrices so that X(:),Y(:) would
%be all the x,y pairings


lat2d_fine=interp2(x_grid,y_grid,squeeze(lat2d.var),...
    X,Y);

lon2d_fine=interp2(x_grid,y_grid,squeeze(lon2d.var),...
    X,Y);





    
