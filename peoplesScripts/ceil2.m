function b = ceil2(a, n)
%FIX2 Round to a specified number of decimals towards zero.
%
%   FIX2(A) returns the same as FIX(A).
%   Y = FIX2(A, N) fixes the elements of A to decimals specified in N.
%
%   See also: FIX, FLOOR, CEIL, ROUND, ROUND2.


% simpler case(round towards nearest integers)
if nargin < 2 | ~any(n) | isempty(a)
    b = ceil(a); return
end

% general case(round to n decimals )
%if ~isscalar(a) & ~isscalar(n) & ~isequal(size(a),size(n))
%    error('Non-scalar inputs must be the same size')
%end
t2n = 10.^floor(n);
b = ceil(a.*t2n)./t2n;