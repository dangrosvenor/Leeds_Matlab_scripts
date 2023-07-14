clear depM depN

q39norm=NaeroDetQ(:,:,39)/NnewT(39);
depN(:,:,39)=NaeroDetQ(:,:,39);
depM(:,:,39)=depN(:,:,39)*mmav2(39);

for iq=1:38
    depN(:,:,iq)=q39norm*NnewT(iq);
    depM(:,:,iq)=depN(:,:,iq)*mmav2(iq);
end

