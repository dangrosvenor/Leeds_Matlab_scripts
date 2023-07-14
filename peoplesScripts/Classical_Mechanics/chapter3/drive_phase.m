%drive_phase.m
%plots the phase difference between the driving force and
%the solution for a driven HO
clear;
m=0.5;                                %mass
k=0.5;                                %spring constant
wo=sqrt(k/m);                         %SHO natural frequency
wmin=0.0;                             %minimum frequency
n=2.5;                                %used to increase wmax and NPTS
wmax=n*wo;                            %maximum frequency in terms of wo
NPTS=n*33+1;                          %w points
wstep=(wmax-wmin)/NPTS;               %w step size
hold on                               %get ready to superimpose plots
cmin=0.01;                            %minimum value of c
cmax=1;                               %maximum value of c
cstep=(cmax-cmin)/5;                  %c step size
for c=cmin:cstep:cmax,                %loop over the drag coefficient
 gam=c/2/m;                           %gamma
 for i=1:NPTS
   w(i)=wmin+(i-1)*wstep;
   den=wo^2-w(i)^2;
     if(w(i)<=wo)
       ph(i)=atan(2*gam*w(i)/den);    %The phase difference    
     else
       ph(i)=pi+atan(2*gam*w(i)/den); %shift by pi needed if w>wo
     end
 end
 plot(w,ph);        %plot phi
                    %num2str(c,p) converts c to a string with p digits
                    %cat(2,'a','b') concatenates a and b
 str=cat(2,'\gamma=',num2str(gam,2));
 text(w(20),ph(20)*(1+0.05),str,'FontSize',8,'Color','red');
 text(w(78),ph(78),str,'FontSize',8,'Color','red');
end
title('Phase Difference Between Driving Force and Solution','FontSize',14)
ylabel('\phi(\omega_D)','FontSize',14);
xlabel('\omega_D','FontSize',14);
