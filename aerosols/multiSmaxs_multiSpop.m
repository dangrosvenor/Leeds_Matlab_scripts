

W=[2:2:50];

multipop %creates an aerosol distribution composed of an accumulation and nuclei mode

for i=1:length(W)
    sm=fzerodan(@nenes,[min(snew2) max(snew2)],[],10e4,T,W(i),Ntot,snew2);
    [I(i),sp(i),smax(i)]=nenes(sm,10e4,283,W(i),Ntot,snew2);
%    is=findheight(snew2,smax(i));
    is=findheight_nearest(snew2,smax(i));    
    Ndrops(i)=sum(Ntot(1:is));
end



figure;plot(W,Ndrops);