%read Nimbus-7 sea ice data from NSDIC
%data is stored as bytes in a flat map of a sterographic plot
%see http://nsidc.org/data/docs/daac/nsidc0051_gsfc_seaice.gd.html#format


shapefile_dir = '/home/disk/eos1/d.grosvenor/matlab/work/m_map1.4/topography files/';
esri_shapefile = 'cst00_linestring';

%read the shapefile
S=shaperead([shapefile_dir esri_shapefile]);

%set the map details of the input file (need to know beforehand)
lat_shape = [-60 -60 -61.18470 -90];
lon_shape = [0 -90 -54.05071 0];
%x_extent = [-2800000 -1000000];
%y_extend = [200000 2000000];

D=3.5e6;
D=3.3357e+06;
%x_extent = [0 -3e6];
%y_extent = [3e6 0];

x_extent = [0 D];
y_extent = [D 0];

%y_extent = [0 1.874e6];
%x_extent = [-2.589e6 0];


LX_shp = abs(diff(x_extent));
LY_shp = abs(diff(y_extent));

lon_stereo = 0; %rotation of the stereo projection - this corresponds to the 6 o'clock position
%Getting the longitude rotation right is critical
lat_stereo = 60; %latitude to which the plot extends
lat_pole=-90; %-90 for south pole

figure
m_proj('stereographic','lat',lat_pole,'lon',lon_stereo,'rad',90-abs(lat_stereo)); 
%N.B. the latitude extension doesn't matter here since we allow the map
%tool to go beyond the edge of the map
%is 180+lon_stereo because m_map puts the specified longitude at the bottom
%of the map



if lat_pole<0
    m_grid('xaxislocation','top','yaxislocation','middle','tickdir','out','linest','-','backcolor',[0 0 0],'color','c');   %white doesn't come out when saved for some reason so can use cyan ('c')
else
    m_grid('xaxislocation','bottom','yaxislocation','middle','tickdir','out','linest','-','fontsize',14,'backcolor',[0 0 0],'color','k');
end


%stereo_ylims = get(gca,'ylim');
%stereo_xlims = get(gca,'xlim');





        [stereo_x,stereo_y] = m_ll2xy([lon_shape],[lat_shape],'clip','off');

        
%pick the two reference points we want to use
  stereo_xlims = stereo_x([1 2]);
  stereo_ylims = stereo_y([1 2]);
 
 LX = abs(diff(stereo_xlims));
 LY = abs(diff(stereo_ylims));
 
 %now run through all of the line vectors and convert from x,y to lat lon
 %using projection that the coastline was drawn on
% for k=1:length(M.ncst)
% x=M.ncst{k}(:,1);
% y=M.ncst{k}(:,2);
% 
% 
% end

%get the size of the vector
NN=0;
for i=1:length(S)
    NN=NN+length(S(i).X(:));
end
% 
% 
 X=NaN*ones([1 NN]);
 Y=NaN*ones([1 NN]);
 
 LATshp=NaN*ones([1 NN]);
 LONshp=NaN*ones([1 NN]);
%ipos=1;
for i=1:length(S)
x=S(i).X(:);
y=S(i).Y(:);

SHP(i).X = x * LX/LX_shp;
SHP(i).Y = y * LY/LY_shp;
% 
% X(ipos:ipos+length(x)-1) = SHP(i).X;
% Y(ipos:ipos+length(x)-1) = SHP(i).Y;
% 
% ipos = ipos + length(x);

%convert from x,y co-ords to lat-lon
[SHP(i).lon,SHP(i).lat]=m_xy2ll(SHP(i).X,SHP(i).Y);

%LATshp(ipos:ipos+length(x)-1) = SHP(i).lat;
%LONshp(ipos:ipos+length(x)-1) = SHP(i).lon;

end

% X2 = X * LX/LX_shp;
% Y2 = Y * LY/LY_shp;


% seaice_xgrid = linspace(stereo_xlims(1),stereo_xlims(2),nrows);
% seaice_ygrid = linspace(stereo_ylims(1),stereo_ylims(2),ncols);
% [seaice_Xgrid,seaice_Ygrid] = meshgrid(seaice_xgrid,seaice_ygrid);

%do some plotting


for i=1:length(SHP)
  m_line(SHP(i).lon,SHP(i).lat,'color','w');
end




% shp_lat = reshape(shp_lat,[ncols nrows]);
% shp_lon = reshape(shp_lon,[ncols nrows]);
% 
% m_pcolor(shp_lon,shp_lat,seaice_data); shading interp
% 
if lat_pole<0
    m_grid('xaxislocation','top','yaxislocation','middle','tickdir','out','linest','-','backcolor',[0 0 0],'color','c');   %white doesn't come out when saved for some reason so can use cyan ('c')
else
    m_grid('xaxislocation','bottom','yaxislocation','middle','tickdir','out','linest','-','fontsize',14,'backcolor',[0 0 0],'color','k');
end
% 
% hcoast=m_coast('line','linewidth',2,'color',[0.2 0.2 0.2]);
% 
% seaice_year = str2num(char(HEADER(103:108)'));
% seaice_day = str2num(char(HEADER(109:114)'));


%HEADER = fread(fid,hsize,'integer');
%fread(fid,hsize,'float=>double');





% Bytes 	Description
% 1-6 	Missing data integer value
% 7-12 	Number of columns in polar stereographic grid
% 13-18 	Number of rows in polar stereographic grid
% 19-24 	Unused/internal
% 25-30 	Latitude enclosed by polar stereographic grid
% 31-36 	Greenwich orientation of polar stereographic grid
% 37-42 	Unused/internal
% 43-48 	J-coordinate of the grid intersection at the pole
% 49-54 	I-coordinate of the grid intersection at the pole
% 55-60 	Five-character instrument descriptor (SMMR, SSM/I)
% 61-66 	Two descriptors of two characters each that describe the data;
% (for example, 07 cn = Nimbus-7 ice concentration)
% 67-72 	Starting Julian day of grid data
% 73-78 	Starting hour of grid data (if available)
% 79-84 	Starting minute of grid data (if available)
% 85-90 	Ending Julian day of grid data
% 91-96 	Ending hour of grid data (if available)
% 97-102 	Ending minute of grid data (if available)
% 103-108 	Year of grid data
% 109-114 	Julian day of grid data
% 115-120 	Three-digit channel descriptor (000 for ice concentrations)
% 121-126 	Integer scaling factor
% 127-150 	24-character file name
% 151-230 	80-character image title
% 231-300 	70-character data information (creation date, data source, etc.)

%fseek(fid,hsize,'cof');