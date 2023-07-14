function [gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,gcm_area_UM] = UM_lat_lon_intialise_FUNC(lat,lon)
% lat and lon inputs can be eithe 1D or 2D arrays.
% lon can be either 0-360 or -180 to 180 :- this converts them to -180 to
% +180 form.

%Deal with single lat and lon vectors
if size(lat,2)==1
    [gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(lon,lat);
else
    gcm_Plat2D_UM = lat;
    gcm_Plon2D_UM = lon;
end

% N.B. - Values from Python scripts are already rotated to proper lat and lons!
%         switch run_type
%             case 'nested'
%                 [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(gcm_Plat2D_UM,gcm_Plo
%                 n2D_UM,pole_lat,pole_lon);
%        end


i180 = find(gcm_Plon2D_UM>180);
gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;

[gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

[gcm_area_UM] = calc_area_lat_lon2d(gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);




