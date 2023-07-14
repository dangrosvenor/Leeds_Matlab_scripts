terr_fine = -1*ones(size(ufine(:,:,1)));
for iz=1:izmax
    inan=find(ufine(:,:,iz)>-1e9 & terr_fine==-1);
    terr_fine(inan)=zfine(iz);
end
    