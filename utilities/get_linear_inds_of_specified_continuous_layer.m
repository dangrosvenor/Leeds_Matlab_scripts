function [out]=get_linear_inds_of_specified_continuous_layer(siz,iCTH,iCTB,K,T,I,J)
% *** Decided to abandon this - using the max CF within low, middle and
% high categories should be fine (according to Rob).
%From arrays of data in the form [K,T,I,J] (height, time, x and y)
%gets a list of linear indices for continuous data in the K dimension.
%Start and end points of the continuous data are already specified by an array
%of starting indices and ending indices  (a long list of linear indices, with corresponding i,j,k,t)
%For the layers that we want to identify we have indices for a point within
%those layers. These will be in the form of one K index for each (t,i,j) -
%Will convert to linear to make a long list of the required indices (that
%reference the main array).
%Need to match these with a start and end index for the
%layer.


%Since the array being referenced is ordered [K,T,I,J] the linear indices
%go through the K dimension first - so adjacent K indices are continuous in
%height unless they belong to the next T index


all_inds3D = [1:prod(siz(2:end))];
[T,I,J] = ind2sub(siz(2:end),all_inds3D); %these should now be the size of K
%TK = repmat(T(:),[length(K) 1]);
%IK = repmat(I(:),[length(K) 1]);
%JK = repmat(J(:),[length(K) 1]);


%linear indices of the required layers
layer_indices = sub2ind(siz,K,T,I,J);

pos=1;
%Loop through each indenfied layer location
for ic = 1:length(layer_indices)
     %find the next lowest (or equal) index in the iCTH index array - this will be the
     %top of the layer - N.B. the K index here is ordered with top of atmos
     %as K=1 and surface as K=end
    imatch = find(iCTH<=layer_indices(ic));
    iCTH_layer = iCTH(imatch(end));
    
    %find the next highest (or equal) index in the iCTB index array - this will be the
     %bottom of the layer
    imatch = find(iCTB>=layer_indices(ic));
    iCTB_layer = iCTB(imatch(1));
    
    %Now write out a list of the indices of the identified layer (linear
    %indices)
    list = iCTH_layer:iCTB_layer;
    N = length(list);
    out(pos:pos+N-1) = list;
    pos = pos + N;
end