function rwp = calc_LWP(qL,rho,Z)
% Gets or calculates LWP or RWP field from the UM
% At the moment it assumes that Z is a 1D vector and that the height
% dimension of qL is the first dimension

if length(Z) == size(qL,1)
    Z = cat(1,0,Z);  % CAT(DIM,A,B) concatenates the arrays A and B along
                     % the dimension DIM.
end
    
dz = diff(Z);

dz2 = repmat(dz,[1 size(qL,2) size(qL,3)]);

rwp = squeeze( sum(qL.*dz2.*rho,1) );

