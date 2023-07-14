%plots graphs of averages (avall), maximums (maxall) and total integrated (totall) against ccn conc

clear ccns xx yy;

if ~exist('dirs')
    dirs=[1:nplots];
end

for i=1:length(dirs)
    id=dirs(i);
    siccn=findstr('ccn',direcDan(id).dir);
    
    if str2num(direcDan(id).dir(siccn+3:end))
        ccns(i)=str2num(direcDan(id).dir(siccn+3:end));
    else
        ccns(i)=i;
    end
end;
xx=ccns;



yy=avall;

[a I]=sort(xx);
yy=yy(I);
xx=sort(xx);
subplot(3,1,1);
h1=newplot;
plot(xx,yy,'-b',xx,yy,'xk');
if ( max(xx)>0 & max(yy)+abs(max(yy))*0.01>min(yy) )
axis([0 max(xx) min(yy) max(yy)+abs(max(yy))*0.01]);
end
%cons=int2str(con);
title(strcat('Average Value of Column ',cons));

yy=maxall;
yy=yy(I);
subplot(3,1,2);
h2=newplot;
plot(xx,yy,'-b',xx,yy,'xk');
if ( max(xx)>0 & max(yy)+abs(max(yy))*0.01>min(yy) )
axis([0 max(xx) min(yy) max(yy)+abs(max(yy))*0.01]);
end
%cons=int2str(con);
title(strcat('Max Value of Column ',cons));

yy=totall;
yy=yy(I);
subplot(3,1,3);
h3=newplot;
plot(xx,yy,'-b',xx,yy,'xk');
if ( max(xx)>0 & max(yy)+abs(max(yy))*0.01>min(yy) )
axis([0 max(xx) min(yy) max(yy)+abs(max(yy))*0.01]);
end
%cons=int2str(con);
title(strcat('Total Integrated Value of Column ',cons));

text(0,-1.5,direcDan(1).dir,'units','centimeters');

