function plotemmtimHf(temm,Temm,dat,col,savedir,nam)

figure;
pcolor(temm,Temm,dat(:,:,col));
shading flat;
set(gca,'ydir','reverse');
colorbar;
set(gcf,'name',nam);
picname=[savedir nam];
print(gcf,picname,'-dmeta');
close(gcf);
