%conic1.m - script for possible conic section curves with rmim and
% eccentricity e.
clear;
L=14;                   %for viewing
th=[0:0.001:2*pi];      %angle range
rmin=1; N=7;            %any rmin, and curves to draw
emin=-1.8;emax=1.8;     %range of e
de=(emax-emin)/(N-1);   %step size
e=[emin:de:emax];
c=get(gca,'ColorOrder');%get available colors to use
stp=cat(2,'o*vsdph');   %plot symbols (more are .ox+*sdv^><ph ) 
cz=size(c);
for i=1:1:N
r=rmin*(1+e(i))./(1+e(i)*cos(th));
x=r.*cos(th);y=r.*sin(th);m=length(x);ms=5;%count points to plot
j=mod(i,cz(1))+1;       %only cz(1) colors available so 1 < j < 7
plot(x(1:ms:m),y(1:ms:m),stp(j),'MarkerSize',2.5,'Color',c(j,:))%in color order
axis ([-L L -L L])
hold on
str=cat(2,'e(',num2str(i),')=',num2str(e(i)));%curve label
sz=size(str); sc(i,1:sz(2))=str;%store labels in sc array according to length
end
str2=cat(2,'plot: r = rmin*(1+e)/(1+e*cos(\theta)) vs \theta, varying e');
title(str2,'FontSize',14)
h=legend(sc,-1);        % place all labels as legend
set(h,'FontSize',12)
xlabel('x=r*cos(\theta)','FontSize',14)
ylabel('y=r*sin(\theta)','FontSize',14)
