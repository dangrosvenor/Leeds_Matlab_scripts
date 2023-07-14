function [diff]=diff_dan(z,q)

diff=z(2:end)-z(1:end-1);
