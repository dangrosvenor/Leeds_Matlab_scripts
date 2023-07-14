% Platnick (2000) vertical penetration paper.
%Calculations for the representative optical depth for re retreivals using
%the re* values in Table 3 and assuming the adiabatic profile type (Profile
%B in Table 2).Equations for this are given in Section 3.1

%They test 4 clouds with tau_tot = 5, 8, 10 and 15.
clear tauc_Plat tau_star21 tau_star37
%The 4 clouds in Table 3a
%Using the wm estimates for profile B:-
i=1;
tauc_Plat(i)=15; re_bot(i)=4; re_top(i)=10; re_star21(i)=9.3; re_star37(i)=9.9; i=i+1;
tauc_Plat(i)=10; re_bot(i)=6; re_top(i)=15; re_star21(i)=13.5; re_star37(i)=14.4; i=i+1;
tauc_Plat(i)=8; re_bot(i)=5; re_top(i)=12; re_star21(i)=10.6; re_star37(i)=11.4; i=i+1;
tauc_Plat(i)=5; re_bot(i)=8; re_top(i)=12; re_star21(i)=10.7; re_star37(i)=11.1; i=i+1;

%Using the retrievals for profile B
i=1;
tauc_Plat(i)=15; re_bot(i)=4; re_top(i)=10; re_star21(i)=9.3; re_star37(i)=9.9; i=i+1;
tauc_Plat(i)=10; re_bot(i)=6; re_top(i)=15; re_star21(i)=13.6; re_star37(i)=14.5; i=i+1;
tauc_Plat(i)=8; re_bot(i)=5; re_top(i)=12; re_star21(i)=10.7; re_star37(i)=11.4; i=i+1;
tauc_Plat(i)=5; re_bot(i)=8; re_top(i)=12; re_star21(i)=10.7; re_star37(i)=11.2; i=i+1;

x=-3; %For re(tau) propto tauc_Plat - tau
x=1; %for adiabatic cloud, profile B


clear re
%Loop over clouds
for i=1:length(tauc_Plat)
   n = (2*x+3)/x;
   a0=re_top(i).^n;
   a1=a0-re_bot(i).^n;
   
   %range of tau to consider
   tau=[0:0.01:tauc_Plat(i)];
    
   re{i} = (a0 - a1.*tau./tauc_Plat(i)).^(1/n); %see below Eqn.(2) at the end of the page
   tau_star21(i) = interp1(re{i}(:),tau,re_star21(i));    
   tau_star37(i) = interp1(re{i}(:),tau,re_star37(i));      
end

%% Plots
lwidth=3;
fsize=16;
ylim_set=[0 5];


figure
plot(tauc_Plat,tau_star21,'bo-','linewidth',lwidth,'markerfacecolor','b');
hold on
plot(tauc_Plat,tau_star37,'rs--','linewidth',lwidth,'markerfacecolor','r');
xlabel('\tau_c');
ylabel('d\tau');
legend({'2.1 \mum','3.7 \mum'},'Location','NorthWest');
grid on
fontsize_figure(gcf,gca,fsize);

set(gca,'ylim',ylim_set);

savedir = '/home/disk/eos1/d.grosvenor/modis_work/vert_plots/';
savename = [savedir 'vert_pen_dtau_error_vs_tau'];

% relative error

% figure
% err_constant = ( tauc_Plat ./ (tauc_Plat-3.3) ).^0.5;
% err_var = ( tauc_Plat ./ (tauc_Plat-tau_star21) ).^0.5
% plot(tauc_Plat,err_constant,'bo-');
% hold on
% plot(tauc_Plat,err_var,'rx--');
% xlabel('\tau_c');
% ylabel('N/N_{corr}');
% legend({'d\tau=3.3','d\tau=\tau^*'});
% title('2.1 \mum');
% 
% figure
% err_constant = ( tauc_Plat ./ (tauc_Plat-2.0) ).^0.5;
% err_var = ( tauc_Plat ./ (tauc_Plat-tau_star37) ).^0.5
% plot(tauc_Plat,err_constant,'bo-');
% hold on
% plot(tauc_Plat,err_var,'rx--');
% xlabel('\tau_c');
% ylabel('N/N_{corr}');
% legend({'d\tau=2.0','d\tau=\tau^*'});
% title('3.7 \mum');

%% relative error plot

lwidth=3;
fsize=16;
ylim_set=[1 1.5];

err_var21 = ( tauc_Plat ./ (tauc_Plat-tau_star21) ).^0.5;
err_var37 = ( tauc_Plat ./ (tauc_Plat-tau_star37) ).^0.5;


figure
plot(tauc_Plat,err_var21,'bo-','linewidth',lwidth,'markerfacecolor','b');
hold on
plot(tauc_Plat,err_var37,'rs--','linewidth',lwidth,'markerfacecolor','r');
xlabel('\tau_c');
ylabel('N_{standard} / N_{re(h*)}');
legend({'2.1 \mum','3.7 \mum'},'Location','NorthWest');
grid on
fontsize_figure(gcf,gca,fsize);

set(gca,'ylim',ylim_set);

savedir = '/home/disk/eos1/d.grosvenor/modis_work/vert_plots/';
savename = [savedir 'vert_pen_Nd_error_var_dtau_vs_tau'];

