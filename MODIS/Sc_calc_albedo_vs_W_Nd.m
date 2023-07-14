%Also see albedo_sensitivity.m
% In the above  made a function that follows the steps below to make the
% albedo :-
% [Ac,tau] = albedo_cloudy_func_Seinfeld(W0,N0);

k=0.8;
Q=2;
cw=500;
rhow=1000;

CTT = 278;
P=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
%the pressure dependence is not very strong
Parr=ones(size(CTT))*P; %
[cw]=adlwcgm2_just_robs(CTT,Parr) * 1e-3; 

B = 2*sqrt(10)*rhow.^2*cw.^(0.5).*(5/9).^(5/2) ./ (k.*pi.*Q.^3);


W=[10:300]*1e-3;
H = (2*W./cw).^0.5;

%% Plot albedo vs LWP

figure
clear labs
hold on

N=10;Nm3=N*1e6;

tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
plot(W*1e3,A,'r');
labs{1} = ['Nd=' num2str(N) ' cm^{-3}'];

N=100; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
plot(W*1e3,A,'b');
labs{2} = ['Nd=' num2str(N) ' cm^{-3}'];

N=1000; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
plot(W*1e3,A,'g');
labs{3} = ['Nd=' num2str(N) ' cm^{-3}'];


xlabel('LWP (g m^{-2})');
ylabel('Albedo');
legend(labs,'location','southeast');
grid on


%% Plot albedo vs cloud thickness

figure
clear labs
hold on

N=10; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
plot(H,A,'r');
labs{1} = ['Nd=' num2str(N) ' cm^{-3}'];

N=100; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
plot(H,A,'b');
labs{2} = ['Nd=' num2str(N) ' cm^{-3}'];

N=1000; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
plot(H,A,'g');
labs{3} = ['Nd=' num2str(N) ' cm^{-3}'];



xlabel('Cloud thickness (m)');
ylabel('Albedo');
legend(labs,'location','southeast');
grid on


figure
plot(H,W*1e3);
xlabel('Cloud thickness (m)');
ylabel('LWP (g m^{-2})');
%legend(labs,'location','southeast');
grid on


%% dA/dW
figure
clear labs
hold on

N=10;Nm3=N*1e6;

tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
dAdW = (A(3:end)-A(1:end-2)) ./ (W(3:end)-W(1:end-2));
plot(W(2:end-1)*1e3,dAdW,'r');
labs{1} = ['Nd=' num2str(N) ' cm^{-3}'];

N=100; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
dAdW = (A(3:end)-A(1:end-2)) ./ (W(3:end)-W(1:end-2));
plot(W(2:end-1)*1e3,dAdW,'b');
labs{2} = ['Nd=' num2str(N) ' cm^{-3}'];

N=1000; Nm3=N*1e6;
tau = (Nm3.*W.^(5/2)./B).^(1/3);
A = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
dAdW = (A(3:end)-A(1:end-2)) ./ (W(3:end)-W(1:end-2));
plot(W(2:end-1)*1e3,dAdW,'g');
labs{3} = ['Nd=' num2str(N) ' cm^{-3}'];



xlabel('LWP (g m^{-2})');
ylabel('dA/dLWP');
legend(labs,'location','northeast');
grid on
set(gca,'yscale','log');

%% Now look at albedo vs Nd for fixed LWP
figure
clear labs
hold on

N_range = [10:0.1:1500]*1e6; %per m3

W0=75e-3; %kg/m2
[Ac_W75,tau_W75] = albedo_cloudy_func_Seinfeld(W0,N_range);
plot(N_range/1e6,Ac_W75,'r','linewidth',3);
labs{1} = ['LWP=' num2str(W0*1e3) ' g m^{-2}'];


W0=125e-3; %kg/m2
[Ac_W125,tau_W125] = albedo_cloudy_func_Seinfeld(W0,N_range);
plot(N_range/1e6,Ac_W125,'b','linewidth',3);
labs{2} = ['LWP=' num2str(W0*1e3) ' g m^{-2}'];


W0=250e-3; %kg/m2
[Ac_W250,tau_W250] = albedo_cloudy_func_Seinfeld(W0,N_range);
plot(N_range/1e6,Ac_W250,'g','linewidth',3);
labs{3} = ['LWP=' num2str(W0*1e3) ' g m^{-2}'];


xlabel('N_d (cm^{-3})');
ylabel('A_c');
legend(labs,'location','northeast');
set(gca,'xlim',[0 200]);
grid on

%% dA/Nd for fixed LWP
figure
clear labs
hold on

N_range = [10:0.1:1500]*1e6; %per m3

W0=75e-3; %kg/m2
[Ac_W75,tau_W75] = albedo_cloudy_func_Seinfeld(W0,N_range);
dAdN_W75 = (Ac_W75(3:end)-Ac_W75(1:end-2)) ./ (N_range(3:end)-N_range(1:end-2));
plot(N_range(2:end-1)/1e6,1e6*dAdN_W75,'r','linewidth',3);
labs{1} = ['LWP=' num2str(W0*1e3) ' g m^{-2}'];


W0=125e-3; %kg/m2
[Ac_W125,tau_W125] = albedo_cloudy_func_Seinfeld(W0,N_range);
dAdN_W125 = (Ac_W125(3:end)-Ac_W125(1:end-2)) ./ (N_range(3:end)-N_range(1:end-2));
plot(N_range(2:end-1)/1e6,1e6*dAdN_W125,'b','linewidth',3);
labs{2} = ['LWP=' num2str(W0*1e3) ' g m^{-2}'];


W0=250e-3; %kg/m2
[Ac_W250,tau_W250] = albedo_cloudy_func_Seinfeld(W0,N_range);
dAdN_W250 = (Ac_W250(3:end)-Ac_W250(1:end-2)) ./ (N_range(3:end)-N_range(1:end-2));
plot(N_range(2:end-1)/1e6,1e6*dAdN_W250,'g','linewidth',3);
labs{3} = ['LWP=' num2str(W0*1e3) ' g m^{-2}'];


xlabel('N_d (cm^{-3})');
ylabel('dA/dN_d (cm^3)');
legend(labs,'location','northeast');
set(gca,'xlim',[0 200]);
grid on
%set(gca,'yscale','log');

%Calculate some values for a 100 per cc cloud
xx = N_range(2:end-1)/1e6;

Nval=100;
dAdNd_N100_W75 = interp1(xx,1e6*dAdN_W75,Nval)  %convert from m3 to cm3
dAdNd_N100_W125 = interp1(xx,1e6*dAdN_W125,Nval)
dAdNd_N100_W250 = interp1(xx,1e6*dAdN_W250,Nval)

Nval=50;
dAdNd_N50_W75 = interp1(xx,1e6*dAdN_W75,Nval)
dAdNd_N50_W125 = interp1(xx,1e6*dAdN_W125,Nval)
dAdNd_N50_W250 = interp1(xx,1e6*dAdN_W250,Nval)

Nval=20;
dAdNd_N20_W75 = interp1(xx,1e6*dAdN_W75,Nval)
dAdNd_N20_W125 = interp1(xx,1e6*dAdN_W125,Nval)
dAdNd_N20_W250 = interp1(xx,1e6*dAdN_W250,Nval)


