%Make a box with lots of small lines so that it curves along a line of
%constant latitude

%lat_box = [72 75];
%lon_box = [0 50.0];

lat_box = [72 75];
lon_box = [-3 48];

%lat_box = LAT_val([1 end]);
%lon_box = LON_val([1 end]);

%bottom part of box
N = 200;
lats_bot = lat_box(1)*ones([1 N]);
lons_bot = lon_box(1):(lon_box(2)-lon_box(1))/(N-1):lon_box(2);

lats_top = lat_box(2)*ones([1 N]);
lons_top = lon_box(1):(lon_box(2)-lon_box(1))/(N-1):lon_box(2);


xvec = [lon_box(1) lon_box(2); lon_box(1) lon_box(2)];
yvec = [lat_box(1) lat_box(1); lat_box(2) lat_box(2)];

color_vec  = [0.95 0.95 0.95];

m_line(xvec,yvec,'linewidth',2,'color',color_vec);

m_line(lons_bot,lats_bot,'linewidth',2,'color',color_vec);
m_line(lons_top,lats_top,'linewidth',2,'color',color_vec);

