id=1;
for jj=1:88
 prcs=[0:5:100];
 vap_prctiles.v(1:250,jj,21)=(max(vap.v(:,:,jj),[],2));
end