function [inds,arr_out,min_vals] = min_column_inds(arr, arr2)
%function [inds,arr_out] = min_column_inds(arr, arr2)
% returns the values of the second array (arr2) at the min of the first
% (arr) over the 1st dimension. Also returns the indices (inds). 
%Can leave out arr2 if just want the indices
% D0es the min/max over the FIRST index
arr_out=NaN;
siz=size(arr);

[min_vals,inds] = min(arr,[],1); 

%find the indices for 2nd/original array at the position of the min value within the first array
inds2 = sub2ind([siz(1) siz(2)*siz(3)],inds(:),[1:siz(2)*siz(3)]');

if nargin==2
    arr_out = reshape(arr2(inds2),siz(2),siz(3)); %
end
