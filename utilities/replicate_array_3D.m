function [arr3D_out] = replicate_array_3D(arr1D,arr3D)

arr3D_out=repmat(arr1D,[1 size(arr3D,1) size(arr3D,2)]); 
arr3D_out=permute(arr3D_out,[2 3 1]);