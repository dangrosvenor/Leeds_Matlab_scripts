%drive_amp.m
%plots the amplitude of the solution for a driven HO
clear;
m=0.5;                                %mass
k=0.5;                                %spring constant
F0=0.5;                               %driving force amplitude
wo=sqrt(k/m);                         %SHO natural frequency
wmin=0.1;                             %minimum frequency
wmax=2;                               %maximum frequency
NPTS=200;                             %number of points
w=[wmin:wmax/NPTS:wmax];              %w array
hold on                               %get ready to superimpose plots
cmin=0.2;                             %minimum value of c
cmax=2*m*wo/sqrt(2);                  %maximum c so that om_res is real
cstep=(cmax-cmin)/5;                  %c step size
for c=cmin:cstep:cmax,                %loop over the drag coefficient
    gam=c/2/m;                        %find gamma
    desc=(2*gam*w).^2+(wo^2-w.^2).^2; 
    A=F0/m./sqrt(desc);               %The driven ho amplitude
    plot(w,A);                        %plot amplitude
    om_res=sqrt(wo^2-2*gam^2);        %resonant frequency
    Amax=F0/2/m/gam/sqrt(wo^2-gam^2); %Maximum amplitude at om_res
                                      %next, draw point at position of Amax
    line([om_res;om_res],[Amax;Amax],'Color','red','Marker','.');
                                      %num2str(c,p) c to string with p digits
                                      %cat(2,'a','b') concatenates a and b
    str=cat(2,'\gamma=',num2str(gam,2));
    text(om_res+0.1,Amax-0.06,str,'FontSize',10,'Color','red');
end
title('Amplitude versus Motor Frequency','FontSize',14)
ylabel('A(\omega_D)','FontSize',14);
xlabel('\omega_D','FontSize',14);
