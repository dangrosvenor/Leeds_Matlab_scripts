function [iregion_lin,iregion_lin_edges,Y_out]=get_lat_lon_irregular_with_time(stime,LAT_val,LON_val,Plat2D,Plon2D,Plat2D_edges,Plon2D_edges,Y)
%function [iregion_lin,iregion_lin_edges,Y_out]=get_lat_lon_irregular_with_time(stime,LAT_val,LON_val,Plat2D,Plon2D,Plat2D_edges,Plon2D_edges,Y)
%If Y is supplied then it makes the data in Y NaN outside of the lat lon
%values that are input.
%Otherwise returns the indices for a 3D matrix (with 3rd index as time) for regions
%within LAT_val and LON_val. Then can do for e.g.
%        a=NaN*ones(size(Y));
%        a(iregion_lin)=0;
%        Y=Y+a;
% To make the regions we don't want NaN in Y.

Plat3D = repmat(Plat2D,[1 1 stime]);
Plon3D = repmat(Plon2D,[1 1 stime]);
Plat3D_edges = repmat(Plat2D_edges,[1 1 stime]);
Plon3D_edges = repmat(Plon2D_edges,[1 1 stime]);

iregion_lin = find(Plat3D>=LAT_val(1) & Plat3D<LAT_val(end) & Plon3D>=LON_val(1) & Plon3D<LON_val(end));
iregion_lin_edges = find(Plat3D_edges>=LAT_val(1) & Plat3D_edges<LAT_val(end) & Plon3D_edges>=LON_val(1) & Plon3D_edges<LON_val(end));   

if nargin>7
    a=NaN*ones(size(Y));
    a(iregion_lin)=0;
    Y_out = Y;
    Y_out=Y_out+a;
end

