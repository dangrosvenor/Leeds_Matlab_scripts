function [arr_out] = apply_inds_from_one_dim_in_3d_array(inds, arr2)
%function [arr_out] = apply_inds_from_one_dim_in_3d_array(inds, arr2)
% returns the values of the array (arr2) at the indices for the FIRST
% dimension as specified by inds. So if arr2 is of size [M N P] then inds
% should be of size [N P] and refer to the location in the M dimension for
% each N,P gridbox.
% inds need to refer to FIRST dimension of arr2 and the size of the array
% used to find inds should be the same size as arr2

arr_out=NaN;
siz=size(arr2);

%[min_vals,inds] = min(arr,[],1); 

%find the indices for 2nd/original array at the position of the min value within the first array
inds2 = sub2ind([siz(1) siz(2)*siz(3)],inds(:),[1:siz(2)*siz(3)]');

if nargin==2
    arr_out = reshape(arr2(inds2),siz(2),siz(3)); %
end
