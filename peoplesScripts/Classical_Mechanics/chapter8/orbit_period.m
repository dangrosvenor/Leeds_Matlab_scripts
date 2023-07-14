%orbit_period.m
%Script made to obtain the time it takes to go from rmin to rmax
%in an orbit due to a force of the form -a*r^p
clear;
a=108.0;             %force strength in N/m^p
p=1;                 %power of r
r0=1; v0=6.0; th0=0; %initial position(m),tangential velocity(m/s), angle(rad)
%ra is an approximate value of rmin needed to estimate vr in the next line
ra=0.4813; vr=v0*sqrt((r0/ra)^2-1); %vr=radial velocity component needed in E
m=1;                 %object's mass in kg
L=m*v0*r0;           %angular momentum, v0=tangential velocity
% Below E is the constant energy, I is the integrand definition
if p ==-1 
    E=0.5*m*vr^2+L^2/m/r0^2/2+a*log(r0);
    I = inline('1./sqrt(2*(E-a*log(r)-L^2/m./r.^2/2)/m)',...
        'r','a','p','E','m','L');
else
    E=0.5*m*vr^2+L^2/m/r0^2/2+a*r0^(p+1)/(p+1);
    I = inline('1./sqrt(2*(E-a*r.^(p+1)/(p+1)-L^2/m./r.^2/2)/m)',...
        'r','a','p','E','m','L'); 
end
% search the min, max turning points
st=0.0002; rlim=10*r0; % search step size and upper limit of search
for r=[r0:-st:0]
    if real(I(r,a,p,E,m,L)) <= 0, rmin=r+st; break;
    end
end
for r=[r0:st:rlim]
    if real(I(r,a,p,E,m,L)) <=0,  rmax=r-st; break;
    end
end
T = 4*quad(I,rmin,rmax,[],[],a,p,E,m,L);%period is 4 times area
%r=[rmin:0.01:rmax]; plot(r,I(r,a,p,E,m,L))%integrand plot
r=[rmin:0.01:rmax]; 
area(r,I(r,a,p,E,m,L),'FaceColor',[0.7 0.8 0.9])%area plot
str1=cat(2,'Central Force (-a*r^p',', p=',num2str(p,3),')',...
    ' Integrand & Period');
str2=cat(2,'(rmin,rmax)=(',num2str(rmin,3),...
    ',',num2str(rmax,3),'),',' \tau = ',num2str(T,3),' sec');
axis([rmin rmax 0 0.6]);xlabel('r(m)','FontSize',14);
ylabel('Integrand & Area','FontSize',14);title(str1,'FontSize',14)
text(0.6,0.1,str2,'FontSize',14)