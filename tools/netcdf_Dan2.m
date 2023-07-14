function [dat] = netcdf_Dan(varargin)
%function [dat] = netcdf_Dan(nc_out,var_str,start_vals,end_vals)
% For old version use as nc_out = nc file ID, var_str=e.g. 'time'
%ind_str needs to be something like '(:)' or '(1,:,:)', etc.
% For new version nc_out is the filename, ind_str is as for above (but
% optional)
        
version_mat=version;
vfind = strfind(version_mat,'R2007b');    
if length(vfind)>0; version_mat='R2007b'; end

switch version_mat
    case 'R2007b' 
        if ~exist('ind_str')
            ind_str='(:)'
        end
        if length(varargin)==4   
            ind_str = '(';
            for i_inds=1:length(varargin{3})
               if  varargin{4}(i_inds)==Inf
                   inds{i_inds} = [num2str(varargin{3}(i_inds)) ':end'];
               else
                   inds{i_inds} = [num2str(varargin{3}(i_inds)) ':' num2str(varargin{4}(i_inds))];
               end
               
               ind_str=[ind_str inds{i_inds} ','];               
            end
            ind_str=[ind_str(1:end-1) ')'];  
%            dat = double(ncread(varargin{1},varargin{2},varargin{3},varargin{4}));  %ncread(SOURCE,VARNAME,START, COUNT, STRIDE)

            eval_str = [ 'varargin{1}{' varargin{2} '}(' ];
%            inds{1} ',' inds{2} ',' inds{3} ',' inds{iv} ',' 
        else
            dat = double(ncread(varargin{1},varargin{2}));
        end
        dat = eval(['nc_out{var_str}' ind_str ';']);
    otherwise
%       if ~exist('ind_str')
%           ind_str = ones([1 100]);
        if length(varargin)==4   
            dat = double(ncread(varargin{1},varargin{2},varargin{3},varargin{4}));  %ncread(SOURCE,VARNAME,START, COUNT, STRIDE)
        else
            dat = double(ncread(varargin{1},varargin{2}));
        end
        if size(dat,2)~=1 % if is a vector then is the same way aroudn as old NetCDF method
            ndims = length(size(dat));
            per = [ndims:-1:1];
            dat = permute(dat,per);
        end
        
end     
        
        
    