%kepler3rd.m - This script is designed to look at the solar system planets
%and to look at the tau versus a relationship.
clear;
planet=['Mercury ';'Venus   ';'Earth   ';'Mars    ';'Jupiter ';...
        'Saturn  ';'Uranus  ';'Neptune ';'Pluto   '];
p=[0.241;0.615;1.000;1.881;11.86;29.46;84.01;164.8;248.5];%period
a=[0.387;0.723;1.000;1.524;5.203;9.539;19.191;30.061;39.529];%semimajor axis
x=a.^3; y=p.^2;% define the variables
[c,e]=polyfit(x,y,1);%Linear fit of y(x) can give coefficient b, error e(3)
[cp,ep]=polyfit(y,x,1);%Linear fit of x(y) needed to calculate r-squared
rr=c(1)*cp(1); %definition of r_squared (Bevington's Data Reduction...)
ybest1 = c(1).*x + c(2);%used to compare predicted and actual
x2=[min(x):(max(x)-min(x))/9:max(x)*(1+0.15)];%more points for the prediction
ybest2 = c(1).*x2 + c(2);                  %general predicted
subplot(1,2,1);plot(x,y,'k.',x,ybest1,'o',x2,ybest2,'r:')
xlabel('a^3 (AU^3)','FontSize',14); ylabel('\tau^2 (yr^2)','FontSize',14)
title('\tau^2 Versus a^3','FontSize',14)
h=legend('Planets','Comparison',' Fit',4); set(h,'FontSize',12)
str1=cat(2,' Line Fit=',num2str(c(1),4),'*x + ',num2str(c(2),4));
str2=cat(2,' r^2=',num2str(rr,10));
text(x(6),y(9)*(1+0.05),str1,'FontSize',9);%post rmin
text(x(6),y(9),str2,'FontSize',9);         %post rmin
%---- tau = a^{3/2} figure
xth=[0:0.1:max(a)];                        %variable for pth
pth=xth.^(3/2);                            %theoretical curve
subplot(1,2,2);semilogy(a,p,'k.',xth,pth,'r:')%use semilog for better viewing
for i=1:9
text(a(i),p(i),planet(i,:),'Color','b','FontSize',12);%The planets
end
str3=cat(2,'\tau = a^{3/2}');
h=legend('Planets',str3,4); set(h,'FontSize',14)
xlabel('a (AU)','FontSize',14); ylabel('Period (yr)','FontSize',14);
title('\tau versus a','FontSize',14)