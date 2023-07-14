% calculate the autocorrelation function of A, A must be a column vector
% Author: Sheng Liu
% Email: ustc.liu@gmail.com
% Date: 7/16/2015
% *** N.B. - NEEDS A to be the values with the mean subtracted (y-mean(y))
% !! ***
function x = autocorrel(A)
% get the size of A
[row,col] = size(A);
if (row ~= 1 && col ~= 1)
    error('The input should be a vector, not a matrix!');
end
if row == 1
    A = A';
end
N = length(A);
x = zeros(N,1);
x(1) = sum(A.*A);
for ii = 2:N
    nshift = (ii-1);
    B = circshift(A,-nshift); %Shift B backwards by number (ns) varying from 1 upwards
    inds = [1:N-nshift]; %Only use the first N-ns numbers
    B = B(inds); 
    x(ii) = sum(B.*A(inds));
end
x = x/x(1);


