%drive_power.m
%plots the power supplied by the driving force versus frequency
clear;
m=0.5;                                %mass
k=0.5;                                %spring constant
F0=0.5;                               %driving force amplitude
wo=sqrt(k/m);                         %SHO natural frequency
wmin=0.01;                            %minimum frequency
wmax=3;                               %maximum frequency
NPTS=200;                             %number of points
dw=(wmax-wmin)/(NPTS-1);              %step for w
hold on                               %get ready to superimpose plots
cmin=0.2;                             %minimum value of c
cmax=1;                               %maximum c
cstep=(cmax-cmin)/3;                  %c step size
for c=cmin:cstep:cmax,                %loop over the drag coefficient
  gam=c/2/m;                          %find gamma
  for i=1:NPTS
    w(i)=wmin+(i-1)*dw;
    desc=(2*gam*w(i))^2+(wo^2-w(i)^2)^2; 
    A=F0/m/sqrt(desc);         %The driven ho amplitude
    den=wo^2-w(i)^2;
    if den==0, den=1.e-3; end
    if(w(i) <= wo)
    ph=atan(2*gam*w(i)/den);   %the phase difference
    else
    ph=pi+atan(2*gam*w(i)/den);%shift by pi needed if w > wo
    end
    power(i)=0.5*F0*A*w(i)*sin(ph);
  end
  plot(w,power)
  [p,j]=max(power);            %point where the power is maximum
  str=cat(2,'\gamma=',num2str(gam,3));
  text(w(j),power(j)+0.02,str,'FontSize',10,'Color','red');    
end
title('Average Power Supply vs Drive Frequency','FontSize',14)
ylabel('Power','FontSize',14);
xlabel('\omega_D','FontSize',14);
