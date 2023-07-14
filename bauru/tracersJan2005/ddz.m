function [ddz,dz]=ddz(z,q)

dz=z(2:end)-z(1:end-1);
%dz=repmat(dz,[1 size(q,2)]);
%ddz=( q(2:end,:)-q(1:end-1,:) ) ./ dz;