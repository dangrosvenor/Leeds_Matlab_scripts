N0=[10:10:2000]*1e6;  %per m3
W0 = 100/1e3; %kg/m2

%calculates cloud albedo based on tau/(tau+7.7) and some assumptions about
%cloud temperature etc. for adiabaticity
[Ac,tau] = albedo_cloudy_func_Seinfeld(W0,N0);


dAdW = 5/6 * Ac.*(1-Ac) ./ W0;
figure
plot(N0/1e6,dAdW);
xlabel('N_d (cm^{-3})');
ylabel('dA/dLWP (m^2 kg^{-1})');


%% Check that the equation is right by calculating the gradient of dA/dW at
%% N=50 and N=1000 per cc at W=100 g/m2 - seems to be accurate

figure
N=50e6;
dW=10;
W=[10:dW:400];
[Ac50,tau50] = albedo_cloudy_func_Seinfeld(W/1e3,N);
plot(W,Ac50);
hold on

N=1000e6;
W=[10:dW:400];
[Ac1000,tau1000] = albedo_cloudy_func_Seinfeld(W/1e3,N);
plot(W,Ac1000,'r');

g50=gradient(Ac50,dW/1e3);
g1000=gradient(Ac1000,dW/1e3);

[minval,i]=min(abs(W-100));
g50(i)
g1000(i)

[minval,i]=min(abs(N0-50e6));
dAdW(i)
[minval,i]=min(abs(N0-1000e6));
dAdW(i)


