summ=0;
sumn=0;
for i=2:40 %droplets now stored in QNloss01
    if i<10
        dgs=strcat('QNloss0',num2str(i));
    else
        dgs=strcat('QNloss',num2str(i));
    end
    dgfind=findhead(dgs,dgstrDan(1).dg);
    summ=summ+TimeAv.DGAV(:,dgfind(1))*mmav2(i-1);
    sumn=sumn+TimeAv.DGAV(:,dgfind(1));
end

figure;plot(summ,Grid.Z);