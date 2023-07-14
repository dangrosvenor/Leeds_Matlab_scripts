 function YI = lininterp1f_multidim_RUN(X,Y,XI,dim)
%  F = lininterp1f_multidim_RUN(X,Y,XI,Ydefault,dim)
%    Performs 1D interpolation on arrays of vectors using a fast C
%    1D interpolation routine (running the 1D routine for arrays using for loops from
%    Matlab is painfully slow!)
%    >>> X or Y are arrays of vectors, e.g. of size [N M P Qetc]
%    At the moment it operates assuming that either X is an array (size [N M P Qetc]) and Y is a vector (size [1 N]; 
%    i.e. the same Y-values for all of the X-vectors), OR 
%    vice versa (size(X)=[1 N]; size(Y)=[N M P Qetc]; the same X-values for a series of Y-vectors)
%     OR that both are arrays
%    >>> dim specifies the X dimension in the array (does not have to be the
%    first dimension, e.g. can have [M P N Q], if specify dim=3)
%    The X vectors must be increasing monotonically, but do not have to be linearly spaced
%    >>> XI = the x values at which the interpolated values are required
%    (1D vector at the moment)
%
%    Returns NaN values where the requested values are out of range
%    -----  Daniel Grosvenor (1D C routine written by Umberto Picchini) ----



% lininterp1f_multidim requires that  the vectors for the arrays are ordered so that the interpolation data is 
%  stored in the first dimension. So, need to permute the arrays to make
%  this so.

sX = size(X);
sY = size(Y);

if length(sX)==2 & min(sX)==1 & length(sY)==2 & min(sY)==1 %both vectors
    %but need them to both be the same size
    if sX~=sY
%        error('X and Y are both vectors, but they are not the same size - unsure whether dim applies to X or Y!!');
         if sY(1)==1  %if Y is of size [1 N] then flip
             Y=Y'; 
         else %otherwise flip X
             X=X';
         end         
         dim=1;
         %N.B. - will enter the first if below too
    end
end

if length(sX)==2 & min(sX)==1 % X is a vector, assume Y is the array
    Y = shiftdim(Y,-1); %make sure the array has at least 3 dimensions
    sY2 = size(Y);

    ndims=[1:length(sY2)];
    %create a vector with all of the dimensions except the required one
    ndims2=ndims;
    ndims2(dim+1)=[]; % +1 because we added an extra dimension
    %perumte the dimensions so that the required interpolation dimension
    %(dim) is first
    Y = permute(Y,[dim+1 ndims2]); % +1 because we added an extra dimension

elseif length(sY)==2 & min(sY)==1 % Y is a vector, assume X is the array
    X = shiftdim(X,-1); %make sure the array has at least 3 dimensions
    sX2 = size(X);

    ndims=[1:length(sX2)];
    %create a vector with all of the dimensions except the required one
    ndims2=ndims;
    ndims2(dim+1)=[]; % +1 because we added an extra dimension
    %perumte the dimensions so that the required interpolation dimension
    %(dim) is first
    X = permute(X,[dim+1 ndims2]); % +1 because we added an extra dimension
        
elseif length(sY)==length(sX)
    %both are arrays
    X = shiftdim(X,-1); %make sure the array has at least 3 dimensions
    Y = shiftdim(Y,-1);
    sX2 = size(X);

    ndims=[1:length(sX2)];
    %create a vector with all of the dimensions except the required one
    ndims2=ndims;
    ndims2(dim+1)=[]; % +1 because we added an extra dimension
    %perumte the dimensions so that the required interpolation dimension
    %(dim) is first
    X = permute(X,[dim+1 ndims2]); % +1 because we added an extra dimension
    Y = permute(Y,[dim+1 ndims2]); % +1 because we added an extra dimension
else
    error('Either X must be a vector and Y an array or vice versa, or both arrays');
end

%run the 1D interpolation C routine
YI = lininterp1f_multidim(X,Y,XI,9e99);

%remove the added singleton dimension
dims = [1:length(size(YI))];
dims(2) =[]; %will shift the 2nd dimension to the end as this is the singleton that was created
YI = permute(YI,[dims 2]); % - will have the effect of removing it, but keeping other singleton dimensions, so
%that the shape of the original array is preserved

%Replace the out of range values with NaN
YI(YI>8.9e99)=NaN;

