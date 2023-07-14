function [y] = map_flight_onto_sat(X,Y,sat_dat,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace_5km)  
%function [y] = map_flight_onto_sat(X,Y,sat_dat,mpace_lon_mapped,mpace_lat_mapped,box_type,ilinear_mpace_5km)  
switch box_type
    case 'NxN pixel square';
%        [ilat,ilon,dist_min,ilin]=getind_latlon_quick(Plat,Plon,mpace_lat_mapped,mpace_lon_mapped);
%will have loaded the relevant ilinear_mpace_5km etc data (in
%read_MPACE_met_files)
        y = sat_dat(ilinear_mpace_5km);
    otherwise
       y = interp2(X,Y,sat_dat,mpace_lon_mapped,mpace_lat_mapped);
end