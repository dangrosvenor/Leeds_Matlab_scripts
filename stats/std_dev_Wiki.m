function std = std_dev_Wiki(data)
%Based on https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
% - computing shifted data section
% This uses a shift in the x values to avoid potential rounding error
% issues


data=data(:); %Linear roll out of possible array of data

%Use the first datapoint as the K value through which the data is shifted
K = data(1);

n = 0;
Sum = 0;
Sum_sqr = 0;

for i=1:length(data)
    x = data(i);
    n = n + 1;
    Sum = Sum + x - K;
    Sum_sqr = Sum_sqr + (x - K) * (x - K);
end

std = sqrt( (Sum_sqr - (Sum * Sum)/n)/(n - 1) );
!use n instead of (n-1) if want to compute the exact variance of the given data
!use (n-1) if data are samples of a larger population