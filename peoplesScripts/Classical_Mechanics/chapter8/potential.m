%potential.m - script to plot an attractive potential, energy, etc
%for a body under a central force
m=1; ri=1; vi=5; K=-45; thi=75; %inputs
th=thi*pi/180; vt=vi*sin(th); vr=vi*cos(th);%velocity components
L=m*ri*vt;                      %angular momentum, vt=tangential velocity
E=m*vr^2/2+L^2/m/ri^2/2+K/ri;   %E is constant
e=sqrt(1+2*E*L^2/m/K^2);        %eccentricity
r=[0.2:0.01:5];                 %variable
V=K./r; CT=L^2./r.^2/2/m;       %V(r), and Ang. Momentum term
Veff=V+CT;                      %Effective potential
rmin=1/(-m*K/L^2+sqrt((m*K/L^2)^2+2*m*E/L^2));
rmax=1/(-m*K/L^2-sqrt((m*K/L^2)^2+2*m*E/L^2));
Ec=-m*K^2/L^2/2; rc=-L^2/m/K;
plot(r,V,'k--',r,CT,'-.',r,Veff,'r')
line([rmin;rmax],[E;E],'Color','b','LineStyle',':')
axis([0 4 -60 50])
h=legend('V(r)','L^2/2mr^2','Veff(r)','E',4);
set(h,'FontSize',14)
xlabel('r','FontSize',14); ylabel('Potentials','FontSize',14); hold on;
plot(rmin,E,'k.',rmax,E,'k.',rc,Ec,'k.');%dot rmin,rmax,rc positions
str1=cat(2,'r_{min}=',num2str(rmin,2));
str2=cat(2,'r_{max}=',num2str(rmax,2));
str3=cat(2,'r_c=',num2str(rc,2),', E_c=',num2str(Ec,3));
text(rmin*(1+0.1),E*(1-0.2),str1,'FontSize',10); %post rmin
text(rmax*(1-0.1),E*(1-0.2),str2,'FontSize',10); %post rmax
text(rc*(1-0.8),Ec*(1+0.1),str3,'FontSize',9);   %post rc,Ec
str4=cat(2,' m=',num2str(m),', r_i=',num2str(ri,2),', v_i=',num2str(vi,2),...
  ', \theta_i=',num2str(thi,3),'^0',', K=',num2str(K,3),', L=',num2str(L,3),...
  ', E=',num2str(E,3),', e=',num2str(e,3));
text(0.5,40,str4,'FontSize',11);                 %post the parameters used
title('Potential Plot','FontSize',14)