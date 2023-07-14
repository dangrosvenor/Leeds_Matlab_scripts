%plots graphs of averages (avall), maximums (maxall) and total integrated (totall) against ccn conc

clear ccns xx yy;
for i=1:nplots
    size(direcDan(i).dir);
    ccns(i)=str2num(direcDan(i).dir(43:ans(2)));
end;
xx=ccns;

%figure('name',figlab,'Position',posit);
if maxflag==1
    yy=maxall;
elseif totflag==1
    yy=totall;
else
    yy=avall;
end

[a I]=sort(xx);
yy=yy(I);
xx=sort(xx);
subplot(rowtot,1,rowno);
hs(rowno).h=newplot;


yy=yy./divfac; %factor to change from e.g. m to km

hhh=plot(xx./1e6,yy,'-b',xx/1e6,yy,'dk');
set(hhh,'markersize',5);
axis([0 max(xx)/1e6 min(yy)/1.01 max(yy)*1.01]);
cons=int2str(con);
if rowno==rowtot;
    xlabel('Droplet concentration (/cm^3)','Fontsize',14);
end;



ylabel(ytit,'Fontsize',11);
%title(tit22);

%set(gca,'Fontsize',10);



%text(0,-1.5,direcDan(1).dir,'units','centimeters');

