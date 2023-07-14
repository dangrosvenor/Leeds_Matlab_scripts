function [arr_out] = indices_apply_2D_to_3D_all_column(inds,arr3D)
% given a set of indcies from a 2D array (e.g., specifying x,y locations)
% pick out all of the columns from a 3D array for those x,y positions
% If arr3D is of size [M N P] then inds were calclated from an array of
% size [N P]

arr3D = arr3D(:,:); %If was of size [M N P] will now be [M N*P]
arr_out = arr3D(:,inds);

