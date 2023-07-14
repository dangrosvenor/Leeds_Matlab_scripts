nit=5;
W=10;


clear serr nerr

for ia=1:nit
    smax(ia)=fzeroDan(@nenes,[min(snew) max(snew)],[],10e4,283,W,NnewT,snew);
    serr(ia)=(smax(ia)-Sacc(ia))*100/Sacc(ia);
    [a,b,c,Nn(ia),is,NnewT]=nenes(smax(ia),10e4,283,W,NnewT,snew);
    nerr(ia)=100*(Nn(ia)-Nacc(ia))/Nacc(ia);
end