function [diff]=diffrep(z,q)

diff=z(2:end)-z(1:end-1);
diff=repmat(diff,[1 size(q,2)]);