function CF_out = UM_calc_coarsen_calc_CF(v)

% %Convert all of the variable names in the input structure to actual names
% %for ease of use
% names = fieldnames(UM_time_in);
% for i=1:length(names)
%     eval_str = [names{i} ' = UM_time_in.' names{i} ';'];
%     eval(eval_str);
% end

if v.icoarsen==1
    % Coarse grain the LWP data to GOES resolution (4km) for comparison
    % to that
    [lwp_coarse,lat2,lon2,lat_edges2,lon_edges2] = coarse_grain(v.dat,v.lat,v.lon,v.dlat_target_coarsen,v.dlon_target_coarsen);
end

% --- Now calc the cloud fraction

%Work out the number of points to average over to get 0.25x0.25 degree
%resolution

d=diff(lat2,[],1);
dlat_UM = meanNoNan(meanNoNan(d,1),1);
N = ceil(abs(v.target_dlat/dlat_UM));

d=diff(lon2,[],2);
dlon_UM = meanNoNan(meanNoNan(d,1),1);
M = ceil(abs(v.target_dlon/dlon_UM));

%Cloud fraction will be defined as fraction of points within each N*M
%box with an LWP greater than a threshold. Total no. points =N*M
%Make an array of ones and make the ones that we don't want to
%count zero
Nlwp = zeros(size(lwp_coarse));
Nlwp(lwp_coarse>=v.thresh_LWP_driver)=1;
%Now coarse grain (avearge) the Nlwp array - this will now be our
%cloud fraction
CF_out.cf= reduce_matrix_subsample_mean(Nlwp,N,M);

%Also coarse grain the lat and lon to the same grid as cf
CF_out.lat = reduce_matrix_subsample_mean(lat2,N,M);
CF_out.lon = reduce_matrix_subsample_mean(lon2,N,M);
CF_out.lat_edges = reduce_matrix_subsample_mean(lat_edges2,N,M);
CF_out.lon_edges = reduce_matrix_subsample_mean(lon_edges2,N,M);
