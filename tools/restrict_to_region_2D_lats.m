function [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges,iLAT,iLON,iT,iLAT2,iLON2,iT2] = restrict_to_region_2D_lats(LAT_val,LON_val,Plat2D,Plon2D,Plat2D_edges,Plon2D_edges,stime)


dlat=mean(mean(diff(Plat2D)));
dlon=mean(mean(diff(Plon2D,1,2)));  %for writing the actual boundaries (assuming equal spacing)
%


Plat3D = repmat(Plat2D,[1 1 stime]);
Plon3D = repmat(Plon2D,[1 1 stime]);
Plat3D_edges = repmat(Plat2D_edges,[1 1 stime]);
Plon3D_edges = repmat(Plon2D_edges,[1 1 stime]);






iregion_lin = find(Plat3D>=LAT_val(1) & Plat3D<LAT_val(end) & Plon3D>=LON_val(1) & Plon3D<LON_val(end));
iregion_lin_edges = find(Plat3D_edges>=LAT_val(1) & Plat3D_edges<LAT_val(end) & Plon3D_edges>=LON_val(1) & Plon3D_edges<LON_val(end));
iregion_lin2D = find(Plat2D>=LAT_val(1) & Plat2D<LAT_val(end) & Plon2D>=LON_val(1) & Plon2D<LON_val(end));
iregion_lin2D_edges = find(Plat2D_edges>=LAT_val(1) & Plat2D_edges<LAT_val(end) & Plon2D_edges>=LON_val(1) & Plon2D_edges<LON_val(end));

%Spatial arrays are ordered lat, lon
%Get the i, j and time indices for the 2D array
[iLAT,iLON,iT] = ind2sub(size(Plat3D),iregion_lin);
%Find min and max of these for cases where can assume regularity for all
%times, etc.
% I.e. these are just the ilon and ilat indices that can be used for
% looping
iLAT2 = [minALL(iLAT) maxALL(iLAT)];
iLON2 = [minALL(iLON) maxALL(iLON)];
iT2 = [minALL(iT) maxALL(iT)];


