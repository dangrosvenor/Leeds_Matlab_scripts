%least_action.m - simulates Hamilton's Least Action principle for a
%particle under the action of gravity. The trajectory is compared
%with what is expected analytically
clear; tol=1.e-7;   %clear, and tolerance reasonably small
v0=5; g=9.8; m=1;   %initial speed, gravity and mass, given
y0=0; yf=1;         %initial, final height, given
t0=0; N=15; dy=.01; %initial time, # of points, max change in y allowed
tf=v0/g+sqrt((v0/g)^2-2*(yf-y0)/g); %final time estimate from analytic soln
if (v0/g)^2-2*(yf-y0)/g < 0,        %work with proper v0, y0, yf in this case
    disp 'v0 is too small, or yf is too large, stopped'
    return
end
dt=(tf-t0)/(N-1);              %time step
t=t0:dt:tf;                    %time array
yan=y0+v0*t-g*t.^2/2;          %exact analytic trajectory
yg=y0+(yf-y0)*(t-t0)/(tf-t0);  %trajectory guess: interpolate end points
%Kinetic and Potential energies follow
T(1:N-1)=m*((yg(2:N)-yg(1:N-1))/dt).^2/2; V(1:N-1)=m*g*yg(1:N-1);
L=T-V;             %initial Lagrangian (Kinetic minus Potential energies)
S1=dt*sum(L);      %initial action, rectangular rule - use N-1 terms
yn=yg;             %the initial trajectory is the guess
rand('state',0)    %reset random numbers for reproducible runs , else comment
dS=abs(S1); n=0;   %change in the action, initially it is equal to S1
je=0;              %exchange counter for each trial
subplot(2,1,1)
while dS > tol,    %simulation stops when the change in the action is tiny
 cla;
 n=n+1;            %trial number
 yp=yn;            %trajectory memory, use the saved trajectory
 for k=2:N-1       %y(1), y(N) don't change so start at k=2, end at k=N-1
 yp(k)=yp(k)+2*dy*(rand(size(yp(k)))-0.5);%modify kth step in trajectory 
                                          %memory, use Monte-Carlo step
 T(1:N-1)=m*((yp(2:N)-yp(1:N-1))/dt).^2/2; V(1:N-1)=m*g*yp(1:N-1);%energies
 L=T-V;       %Lagrangian (Kinetic minus Potential energies)
 S2=dt*sum(L);%new action integral by rectangular rule - use N-1 terms
  if S2 < S1       %accept step if the action S decreases  
     dS=abs(S2-S1);%variation of the action (must be small before stopping)
     S1=S2;        %save the lower value of the action
     yn=yp;        %modify trajectory based on accepted step, and keep it
     je=je+1;      %count exchanges made
  else
     yp=yn;        %if the step is not accepted, discard memory changes made,
                   %and refresh the trajectory memory
  end  
 end               
S(n)=S1;           %keep track of the action for each trial
pen=100*sqrt(sum((yn-yan).^2./(yan.^2+.001^2)))/N;%percent error estimate
hold on
plot(t,yan,'k--',t,yg,'b-.',t,yn(:),'r.')%plot analytic,guess, and simulated
str1=cat(2,'Trials: ',num2str(n,4),', %E_n=',num2str(pen,3));
str2=cat(2,'Exchanges:',num2str(je,4));
str3=cat(2,'dS:',num2str(dS,4));
text(t0,max(yan)*(1-0.01),str1);
text(t0,max(yan)*(1-0.13),str2);
text(t0,max(yan)*(1-0.25),str3);
pause (.01)
end
h=legend('Analytic','Init Guess','Least Action',-1); set(h,'FontSize',11)
xlabel('t','FontSize',14),ylabel('y(t)','FontSize',14)
title(['MC: Hamilton''s Least Action Principle',' (tolerance=',... %MC=Monte-Carlo
        num2str(tol,2),')'],'FontSize',12)
subplot(2,1,2), plot(1:n,S(1:n))         %plot the action for each trial
xlabel('trial','FontSize',14),ylabel('Action (S)','FontSize',14)
str4=cat(2,'(y_0, y_f, dy; v_0; g; m; N) = (',num2str(y0,2),', ',...
    num2str(yf,2),', ',num2str(dy,2),'m; ',...
    num2str(v0,2),'m/s; ',num2str(g,2),'m/s^2; ',num2str(m,2),'kg; ',...
    num2str(N,2),')');
text(round(n/40),S(round(n/40)),str4);
