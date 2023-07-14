function qpcolor(dat)
%quick pcolor of the data recieved (opens a new figure and squeezes the data)

figure
[h,hback]=dpcolor(squeeze(dat));
shading flat;
colorbar;
set(hback,'facecolor',[0 0 0]); %the background (NaN) colour
