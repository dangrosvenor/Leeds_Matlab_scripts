function [diff]=diff(z,q)

diff=z(2:end)-z(1:end-1);
