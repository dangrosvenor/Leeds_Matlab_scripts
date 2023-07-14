clear MDetQ flux

'LOAD a run into memory first for dgstrDan'

for j=1:length(diag)

for iq=1:14 %54
    if iq<10
        dgs=strcat('ALL_WQ0',num2str(iq));
    else
        dgs=strcat('ALL_WQ',num2str(iq));

    end
    [MDetQ(j).m(:,iq,:),flux(j).m(:,iq,:)]=detrain(diag(j),dgs,times,dgstrDan(j));
end

end
