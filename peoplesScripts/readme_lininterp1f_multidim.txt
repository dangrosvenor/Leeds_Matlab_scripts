Steps (from Matlab terminal) :

1) Compile the c-file from the Matlab terminal using "mex
   lininterp1f_multidim.c" (cd to the containing directory first)
2) Type "help lininterp1f_multidim_RUN" to get the info on how to run
   the script

%  F = lininterp1f_multidim_RUN(X,Y,XI,Ydefault,dim)
%    Performs 1D interpolation on arrays of vectors using a fast C
%    1D interpolation routine (running the 1D routine for arrays using for loops from
%    Matlab is painfully slow!)
%    >>> X or Y are arrays of vectors, e.g. of size [N M P Qetc]
%    At the moment it operates assuming that either X is an array (size [N M P Qetc]) and Y is a vector (size [1 N]; 
%    i.e. the same Y-values for all of the X-vectors), OR 
%    vice versa (size(X)=[1 N]; size(Y)=[N M P Qetc]; the same X-values for a series of Y-vectors)
%    >>> dim specifies the X dimension in the array (does not have to be the
%    first dimension, e.g. can have [M P N Q], if specify dim=3)
%    The X vectors must be increasing monotonically, but do not have to be linearly spaced
%    >>> XI = the x values at which the interpolated values are required
%    (1D vector at the moment)
%
%    Returns NaN values where the requested values are out of range
%    -----  Daniel Grosvenor (1D C routine written by Umberto Picchini) ----