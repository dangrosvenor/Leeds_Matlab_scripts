function x = autocorrel_Dan(A)
% function x = autocorrel_Dan(A)
% Calculates the autocorrelation coefficients for different lags (0 to
% length(A)).
% Based on :- https://www.itl.nist.gov/div898/handbook/eda/section3/eda35c.htm
% Box and Jenkins 1976
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
me = mean(A);
for ii = 1:N %ii-1 is the lag of autocorr
    nshift = (ii-1);
    B = circshift(A,-nshift); %Shift B backwards by number (ns) varying from 1 upwards
    inds = [1:N-nshift]; %Only use the first N-ns numbers
    B = B(inds); 
    
    ck(ii) = sum( (A(inds)-me).*(B-me) );
end
c0 = sum( (A-me).^2 ); 

x = ck/c0;


