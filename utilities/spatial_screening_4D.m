function [mask_arr] = spatial_screening_4D(arr,LAT_val,LON_val,Plat2D,Plon2D)
%function [Plat4D,Plon4D,Plat4D_edges,Plon4D_edges] = spatial_screening_4D(Plat2D,Plon2D,N3,N4)

 siz = size(arr);
 Plat4D = repmat(Plat2D,[1 1 siz(3) siz(4)]);
 Plon4D = repmat(Plon2D,[1 1 siz(3) siz(4)]);
        
 %Plat4D_edges = repmat(Plat2D_edges,[1 1 N3 N4]);
 %Plon4D_edges = repmat(Plon2D_edges,[1 1 N3 N4]);
 
 
iregion_lin = find(Plat4D>=LAT_val(1) & Plat4D<LAT_val(end) & Plon4D>=LON_val(1) & Plon4D<LON_val(end));
%iregion_lin_edges = find(Plat3D_edges>=LAT_val(1) & Plat3D_edges<LAT_val(end) & Plon3D_edges>=LON_val(1) & Plon3D_edges<LON_val(end));        

 mask_arr = NaN*ones(size(arr));
 mask_arr(iregion_lin) = 0;
% arr = arr + mask_arr;