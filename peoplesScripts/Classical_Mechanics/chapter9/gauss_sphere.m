%gauss_sphere.m - plots the gravitational field for a sphere of mass M and 
%inner radius a and outer radius b.
clear;
a=3; b=6; N=250; %a, b in meters, and points to plot
rmin=0; rmax=2*b; dr=(rmax-rmin)/N; %ranges
for i=1:N
    r(i)=rmin+(i-1)*dr; g(i)=0;
    if r(i)>=a & r(i)<b g(i)=(r(i)-a)/(b-a)/r(i)^2;
    else if r(i)>=b g(i)=1/r(i)^2;
        end
    end
end
[c,j]=max(g);%get maximum of g and its index
plot(r,g,'k',a,0,'r.',b,0,'r.')
text(a*(1+0.03),0.0008,'a','FontSize',14)
text(b*(1+0.03),0.0008,'b','FontSize',14)
line([a,a],[0,c/2],'Color','r','LineStyle','--')
line([b,b],[0,c/2],'Color','r','LineStyle','--')
text(a/2,c/3,'I','FontSize',14)
text((a+b)/2,c/3,'II','FontSize',14)
text(b+(rmax-b)/2,c/3,'III','FontSize',14)
xlabel('r(m)','FontSize',14);ylabel('-g/GM','FontSize',14);
title('Gravitational Field from Sphere with radii a,b','FontSize',14)
