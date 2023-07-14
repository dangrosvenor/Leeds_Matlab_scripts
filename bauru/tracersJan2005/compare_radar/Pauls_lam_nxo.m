function [LAMBDA,NX0]=Pauls_lam_nxo(Q,NQ,alphax,cx,dx,RHO)

LAMBDA=((NQ.*cx.*gamma(1+alphax+dx)./(Q.*gamma(1+alphax))).^(1./dx));
NX0=(RHO.*NQ.*(LAMBDA.^(1+alphax))./gamma(1+alphax));