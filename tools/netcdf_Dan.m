function [dat] = netcdf_Dan(nc_out,var_str,ind_str)
% For old version use as nc_out = nc file ID, var_str=e.g. 'time'
%ind_str needs to be something like '(:)' or '(1,:,:)', etc.
% For new version nc_out is the filename, ind_str is as for above (but
% optional)
        
version_mat=version;
switch version_mat
    case 'Old version?' 
        if ~exist('ind_str')
            ind_str='(:)'
        end
        dat = eval(['nc_out{var_str}' ind_str ';']);
    otherwise
%       if ~exist('ind_str')
%           ind_str = ones([1 100]);
           
        dat = double(ncread(nc_out,var_str));
        if exist('ind_str')
            dat = eval(['dat' ind_str ';']);
        end
end     
        
        
    