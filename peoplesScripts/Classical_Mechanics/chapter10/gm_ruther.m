%gm_ruther.m - program to do a comparison between the Geiger and Marsden 1913
%experimental data and the Rutherford Scattering formula for the number of particles
clear; warning off;
m=1;                    %projectile mass in units of alpha particle mass
za=2;                   %za=projectile
zt=[47 79];             %zt=target charges (47-> silver, 79->gold)
K=za*zt;                %dimensionless force constant
vb=0.01965;             %velocity units (a_b/tau_b) (in inits of c=light speed)
ma=3730e6;              %alpha particle mass energy in eV
Ene=1e6;                %assume this energy
v0=sqrt(2*Ene/m/ma)/vb; %initial speed in units of vb
factor=K.^2/4/m^2/v0^4; %regular cross-section factor
cr=[0.019 0.01];        %fit coef for comparison with experiment
%================== Experimental Data ========================
%Actual Data from H. Geiger & E. Marsden, Phil. Mag. Vol.25, 605 (1913)
%for angle in degrees and number of scintillations obtained for Ag and Au
dth=[15 22.5 30 37.5 45 60 75 105 120 135 150]; % Exp. angle in Degrees from paper
rth=dth*2*pi/360;       %convert angle to radians
rsig(1,:)=[105400 20300 5260 1760  989 320 136 47.3 33.0 27.4 22.2];%Silver 
rsig(2,:)=[132000 27300 7800 3300 1435 477 211 69.5 51.9 43.0 33.1];%Gold
%==============================================================
N=50;thmin=min(rth-0.05);thmax=max(rth+0.05);ths=(thmax-thmin)/(N-1);
th=thmin:ths:thmax;     %uses experimental angle range in radians, N points
thd=th*360/2/pi;        %angle in degrees, use for plotting
for i=1:2               %do formula on silver (i=1) and gold (i=2)
   for j=1:N
       Nth(i,j)=cr(i)*factor(i)/sin(th(j)/2)^4;
   end
end
semilogy(thd(:),Nth(1,:),'r:',dth(:),rsig(1,:),'bo') %Ag
hold on
semilogy(thd(:),Nth(2,:),'k:',dth(:),rsig(2,:),'md') %Au
xlabel('\Theta (degrees)','FontSize',14)
ylabel('N(\Theta)','FontSize',14)
title('Comparison between Rutherford Formula and Actual 1913 Experiment')
h=legend('Ag-Theory','Ag: Exp. Geiger & Marsden (1913)','Au-Theory',...
    'Au: Exp. Geiger & Marsden (1913)',1);
set(h,'FontSize',13)
