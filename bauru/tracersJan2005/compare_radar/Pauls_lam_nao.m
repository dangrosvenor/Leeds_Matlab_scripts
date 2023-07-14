[LAMBDA NXO]=Pauls_lam_nxo(Q,NQ,alphax,dx)

LAMDA=real((NQ.*cx.*gamma(1+alphax+dx)./(Q.*gamma(1+alphax))).^(1./dx));
NX0=real(RHO.*NQ.*(LAMDA.^(1+alphax))./gamma(1+alphax));