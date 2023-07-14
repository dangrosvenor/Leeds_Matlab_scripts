%NOTE - doesn't work very well for plots with the pole in as it plots a line at -180 degrees for all the 
%latitude lines. One workaround might be to make the line continue to the longitude 180 degress opposite
%by changing e.g. the -60+180 = 120 degree longitudes temporarily to -60 so that the line continues
%the problem would be with labelling as it would label the line on both sides of the poles
%could switch to manual labelling or use text command

nLAT = 3;
nLAT = 9;
nLON = 5;

%probably should scale it according to the domain size rather than just
%d02 and d03
if strfind(fileWRF(1).file,'d03')

    nLAT=36;
    nLON=30;
    
    nLAT=18;
    nLON=15;

elseif strfind(fileWRF(1).file,'d02')

    nLAT=18;
    nLON=15;

else

    nLAT=72;
    nLON=60;

end

sp_LAT = [0.1 0.2 0.5 1 2 5 10 18];
sp_LON = [0.1 0.2 0.5 1 2 5 10 20 40 50 60];

minLAT=minALL(lat2d(1).var);
maxLAT=maxALL(lat2d(1).var);
dLAT = (maxLAT-minLAT)/nLAT;
[a,b] = min(abs(sp_LAT - dLAT));
dLAT2 = sp_LAT(b);

minLON=minALL(lon2d(1).var);
maxLON=maxALL(lon2d(1).var);
dLON = (maxLON-minLON)/nLON;
[a,b] = min(abs(sp_LON - dLON));
dLON2 = sp_LON(b);

conts_LAT = (dLAT2*(ceil(minLAT/dLAT2)) : dLAT2 : dLAT2*(floor(maxLAT/dLAT2)) );
conts_LON = (dLON2*(ceil(minLON/dLON2)) : dLON2 : dLON2*(floor(maxLON/dLON2)) );


%looks like you can't use three number vector specification of colour, just strings:
%r = Red
%g = Green
%b = Blue
%c = Cyan
%m = Magenta
%y = Yellow
%k = Black
%w = White
colour_string = 'w';
%colour_string = 'b';
colour_string = 'k';
colour_string = 'g';
%colour_string = 'r';

if ~exist('xinds')
    xinds=1:length(x_grid);
    yinds=1:length(y_grid);
end

hold on
% xinds2=[xinds(1):10:xinds(end)];
% yinds2=[yinds(1):10:yinds(end)];
% 


xinds2=xinds;
yinds2=yinds;

 dat_lon=lon2d.var(yinds2,xinds2);
 dat_lat=lat2d.var(yinds2,xinds2);

% nfilter=50;
% bfilter=ones([nfilter nfilter])*1/nfilter/nfilter;
% lat_fil=filter2(bfilter,lat2d.var);
% lon_fil=filter2(bfilter,lon2d.var);
% 
% dat_lon=lon_fil(yinds2,xinds2);
% dat_lat=lat_fil(yinds2,xinds2);



%[lines_con_a, lines_con_b]=contour(timesTH(1).t(1:pend),zz(1).z,lon2d.var(yinds,xinds),conts_LON,colour_string);
[lines_con_a, lines_con_b]=contour(x_grid(xinds2),y_grid(yinds2),dat_lon,conts_LON,colour_string);
clabel(lines_con_a, lines_con_b,'labelspacing',200,'color',colour_string); %default 144

%[lines_con_a, lines_con_b]=contour(timesTH(1).t(1:pend),zz(1).z,lat2d.var(yinds,xinds),conts_LAT,colour_string);
[lines_con_a, lines_con_b]=contour(x_grid(xinds2),y_grid(yinds2),dat_lat,conts_LAT,colour_string);
clabel(lines_con_a, lines_con_b,'labelspacing',200,'color',colour_string); %default 144