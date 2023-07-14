%euler_ang.m  - shows Euler angles: ph, th, ps; planes, line of nodes
clear,rc=pi/180; v=1+.1; ph=15; th=25; ps=35; %rc=deg-rad, v=multiplier, angles
phi=ph*rc; the=th*rc; psi=ps*rc;              %angles ph,th,ps in radians
axis ([-.8,.8,-.8,.8,-.8,.8]), view([1,1,1]), hold on; %window view 
x=1;y=1;z=1; r=[[x;0;0;],[0;y;0;],[0;0;z]];   %starting cartesian axes
line([0,r(1,1)],[0,r(1,2)],[0,r(1,3)],'color', 'b','LineStyle','-')
line([0,r(2,1)],[0,r(2,2)],[0,r(2,3)],'color', 'b','LineStyle','-')
line([0,r(3,1)],[0,r(3,2)],[0,r(3,3)],'color', 'b','LineStyle','-')
text(r(1,1)*v,r(1,2)*v,r(1,3)*v,'x','FontSize',11)
text(r(2,1)*v,r(2,2)*v,r(2,3)*v,'y','FontSize',11)
text(r(3,1)*v,r(3,2)*v,r(3,3)*v,'z','FontSize',11)
xlabel('x','FontSize',14,'Position',[0 1 -.8]);
ylabel('y','FontSize',14,'Position',[1 0 -.8]);
zlabel('z','FontSize',14)
%--- phi-rotation about z by ph degrees
rphi=[[cos(phi);-sin(phi);0;],[sin(phi);cos(phi);0;],[0;0;1;]];%z matrix rotation
rp=rphi*r;
line([0,rp(1,1)],[0,rp(1,2)],[0,rp(1,3)],'color', 'k','LineStyle','--')
line([0,rp(2,1)],[0,rp(2,2)],[0,rp(2,3)],'color', 'k','LineStyle','--')
line([0,rp(3,1)],[0,rp(3,2)],[0,rp(3,3)],'color', 'k','LineStyle','--')
text(rp(1,1)*v,rp(1,2)*v,rp(1,3)*v,'x\prime','FontSize',11)
text(rp(2,1)*v,rp(2,2)*v,rp(2,3)*v,'y\prime','FontSize',11)
text(rp(3,1)*v,rp(3,2)*v,rp(3,3)*v,'z\prime=z','FontSize',11)
x1=[0;.5;rp(1,1)/2]; y1=[0;0;rp(1,2)/2]; z1=[0;0;rp(1,3)/2];
h(1)=fill3(x1,y1,z1,'b'); %color ph angle from x-axis
x1=[0;0;rp(2,1)/2]; y1=[0;.5;rp(2,2)/2]; z1=[0;0;rp(2,3)/2];
fill3(x1,y1,z1,'b');      %color ph angle from y-axis
%--- theta-rotation about x' by th degrees from y-axis
rthe=[[1;0;0;],[0;cos(the);-sin(the);],[0;sin(the);cos(the);]];%x' matrix rotation
rpp=rthe*rp;
line([0,rpp(1,1)],[0,rpp(1,2)],[0,rpp(1,3)],'color', 'r','LineStyle','-.')
line([0,rpp(2,1)],[0,rpp(2,2)],[0,rpp(2,3)],'color', 'r','LineStyle','-.')
line([0,rpp(3,1)],[0,rpp(3,2)],[0,rpp(3,3)],'color', 'r','LineStyle','-.')
text(rpp(1,1)*v,rpp(1,2)*v,rpp(1,3)*v,'x\prime\prime','FontSize',11)
text(rpp(2,1)*v,rpp(2,2)*v,rpp(2,3)*v,'y\prime\prime','FontSize',11)
text(rpp(3,1)*v,rpp(3,2)*v,rpp(3,3)*v,'z\prime\prime','FontSize',11)
x1=[0;rp(3,1)/2;rpp(3,1)/2]; y1=[0;rp(3,2)/2;rpp(3,2)/2]; 
z1=[0;rp(3,3)/2;rpp(3,3)/2];
h(2)=fill3(x1,y1,z1,'k'); %color th angle from - z axis
x1=[0;rp(2,1)/2;rpp(2,1)/2]; y1=[0;rp(2,2)/2;rpp(2,2)/2]; 
z1=[0;rp(2,3)/2;rpp(2,3)/2];
fill3(x1,y1,z1,'k');      %color th angle from - y' axis
%--- psi-rotation about z'' by ps degrees
rpsi=[[cos(psi);-sin(psi);0;],[sin(psi);cos(psi);0;],[0;0;1;]];%z'' matrix-rotation
rppp=rpsi*rpp;
line([0,rppp(1,1)],[0,rppp(1,2)],[0,rppp(1,3)],'color', 'g','LineStyle','--')
line([0,rppp(2,1)],[0,rppp(2,2)],[0,rppp(2,3)],'color', 'g','LineStyle','--')
line([0,rppp(3,1)],[0,rppp(3,2)],[0,rppp(3,3)],'color', 'g','LineStyle','--')
text(rppp(1,1)*v,rppp(1,2)*v,rppp(1,3)*v,'x\prime\prime\prime=1','FontSize',11)
text(rppp(2,1)*v,rppp(2,2)*v,rppp(2,3)*v,'y\prime\prime\prime=2','FontSize',11)
text(rppp(3,1)*v,rppp(3,2)*v,rppp(3,3)*v,...
    'z\prime\prime\prime=z\prime\prime=3','FontSize',11)
x1=[0;rpp(1,1)/2;rppp(1,1)/2];y1=[0;rpp(1,2)/2;rppp(1,2)/2];
z1=[0;rpp(1,3)/2;rppp(1,3)/2];
h(3)=fill3(x1,y1,z1,'r'); %color ps angle from x''-axis
x1=[0;rpp(2,1)/2;rppp(2,1)/2];y1=[0;rpp(2,2)/2;rppp(2,2)/2];
z1=[0;rpp(2,3)/2;rppp(2,3)/2];
fill3(x1,y1,z1,'r');      %color ps angle from y'' axis
str1=cat(2,'\phi=',num2str(ph,2),'^o'); str2=cat(2,'\theta=',num2str(th,2),'^o');
str3=cat(2,'\psi=',num2str(ps,2),'^o'); 
h=legend(h,str1,str2,str3,1); set(h,'FontSize',11)
%Draw and rotate planes (rotation uses a minus sign due to internal construction)
[xx,yy,zz]=ellipsoid(0,0,0,x,y,0,50); %flat unit ellipsoids are the desired planes
mesh(xx,yy,zz,'EdgeColor',[0 0 1],'FaceAlpha',[0],'EdgeAlpha',[0.1]); % draw plane
%draw, rotate planes: first about z axis by (-ph), use rphi rot. Mat. & reverse order
xp=rthe(1,1)*xx-rthe(1,2)*yy-rthe(1,3)*zz; %minus sign on off diagonals due to -ph
yp=-rthe(2,1)*xx+rthe(2,2)*yy-rthe(2,3)*zz;
zp=-rthe(3,1)*xx-rthe(3,2)*yy+rthe(3,3)*zz;
%next rotation: about x' axis by (-th) - use rthe rot. matrix, then draw plane
xpp=rphi(1,1)*xp-rphi(1,2)*yp-rphi(1,3)*zp;%minus sign on off diagonals due to -th
ypp=-rphi(2,1)*xp+rphi(2,2)*yp-rphi(2,3)*zp;
zpp=-rphi(3,1)*xp-rphi(3,2)*yp+rphi(3,3)*zp;
mesh(xpp,ypp,zpp,'EdgeColor',[1 0 0],'FaceAlpha',[0],'EdgeAlpha',[.1]); %draw plane