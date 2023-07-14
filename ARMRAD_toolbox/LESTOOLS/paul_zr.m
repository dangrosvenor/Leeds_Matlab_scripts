function z=paul_zr(cx,NX0,LAMBDA,alphax,dx)

%z=0.19*(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./...
%    (RHO.*LAMBDA.^(1+alphax+2.*dx)))./((1E-3).^6));

z=(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./...
    (LAMBDA.^(1+alphax+2.*dx)))./((1E-3).^6));

% M(D) = cx*D^dx

% dBz is then just 10*log10(z)
