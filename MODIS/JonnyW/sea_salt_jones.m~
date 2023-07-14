%Jones et al, JGR, 2001 sea salt param
savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';

u=[0:0.1:15]; %m/s
Nmin=5e6;
ASO4 = 0;

Af = NaN*ones(size(u));
Aj = NaN*ones(size(u));

i2 = find(u<=2);
    Af(i2) = 3.856e6*(1-exp(-0.736.*u(i2)));
    Aj(i2) = 0.671e6*(1-exp(-1.351.*u(i2)));
i17 = find(u>2 & u<=17.5);
    Af(i17) = 10.^(0.095.*u(i17)+6.283);
    Aj(i17) = 10.^(0.0422.*u(i17)+5.7122);
imax = find(u>17.5);
    Af(imax) = 1.5e8*(1 - 97.87*exp(-0.313.*u(imax)));
    Aj(imax) = 3.6e6*(1 - 103.926*exp(-0.353.*u(imax)));
   
A = ASO4+Af+Aj;

N = 3.75e8 * (1 - exp(-2.5e-9.*A));
N2 = 1e-6*max(N,Nmin); %per cc

z=[0:0.2:200];
zf=exp(-z/900);
Af2D = repmat(Af,[length(z) 1]);
Aj2D = repmat(Aj,[length(z) 1]);
z2D = ( repmat(zf,[length(Aj) 1]) )';
Afz=Af2D.*z2D;
Ajz=Aj2D.*z2D;

Az = ASO4+Afz+Ajz;
N2D = 1e-6 * 3.75e8 * (1 - exp(-2.5e-9.*Az));
N2D(N2D<5) = 5;


%% Clarke - based on sea-salt (SS) flux
D = 3e-6; %Divergence (per sec)
% get N = F_SS / (D*zi), where S is the flux of sea-salt particles
N_clarke = 214.*u.^3.41 / (D*1500) / 1e6;  %zi = inversion height (=1000m?)


%% plotting
figure
plot(u,N2);
hold on
plot(u,N_clarke,'r');
xlabel('U (m s^{-1})');
ylabel('N_d (cm^{-3})');
savename=[savedir 'NvsU_Jones2001'];
grid
legend('Jones','Clarke');


figure
Us = [2 5 10 15];
iu = findheight_nearest(u,Us);
    
plot(N2D(:,iu),z);
hold on
legend('u=2','u=5','u=10','u=15')

xlabel('z (m)');
xlabel('N_d (cm^{-3})');
grid
set(gca,'xlim',[0 60]);
savename=[savedir 'NvsZ_Jones2001'];


figure
plot(u,Af)
hold on
plot(u,Aj,'k--')
xlabel('U (m/s)');
ylabel('A');
legend('Af','Aj');
grid
savename=[savedir 'Aj_Af_vsU_Jones2001'];












