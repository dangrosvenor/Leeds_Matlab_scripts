clear nccn_all 

if exist('r_single')==0
    r_single = 117e-9; %Default value - used for volcano as of 17th Nov 2016
end

if exist('N_type')==0
        N_type = 'prescribed number';
        N_type = 'prescribed single aerosol mass';
end

%N_type = 'prescribed number';
%N_type = 'prescribed single aerosol mass';
        
        
% air mass properties
Tc = 10;         % deg C
T = Tc +273.15;  % K
p = 90000;       % Pa
dw = 0.05;
w = [0.01:dw:1]; % m/s

% constants
Mw = 0.18015e-1;  % Molecular weight of water  [kg mol-1]
zetasa = 0.8e-1;  % Surface tension at solution-air interface
Ru = 8.314472;    % Universal gas constant
Rd = 287.05;      % gas constant for dry air
Rv =  461.5;      % gas constant for water vapour
eps = 1.6077;     % (Rv/Rd)
rhow = 997.0;     % water density
g = 9.8;          % gravitational acceleration ms-2
Lv = 0.2501e7;    % Latent heat of vapourization
cp = 1005.0;
Dv = 0.226e-4;    % diffusivity of water vapour in air
ka = 0.243e-1;    % thermal conductivity of air

% aerosol properties
nmodes = 1;   % number of modes
vantHoff = [3.0,3.0];
density = [1777,1777];

m_aero_single = density .* 4/3*pi.*r_single.^3;

massMole = [132.0e-3,132.0e-3];
M_aero0 = [5.e-8 , 5.843200e-10] * 2;
%M_aero0 = [1e-9]; %For comparison with Fig. 3 of Ghan (2011)
M_aero0 = [3e-8]; %For comparison with Fig. 4 of Ghan (2011). That was for r_mode=50nm
  %Need to work out how to calc the mode size. Actually, think they mean
  %median diameter, which is what is calculated here as rd

switch N_type
    case 'prescribed number'
        N_aero0 = [5e+09 , 8.6e+08] * 2;
    case 'prescribed single aerosol mass'
        N_aero0 = M_aero0 ./ m_aero_single;
end

sigma = [1.5 , 1.5];
sigma = [2 , 2]; %sigma=2.0 for the Ghan (2011) paper
%sigma = [2.5 , 2.5];
df=0.01;
factor = [df:df:5];
%factor=1;


for k=1:length(factor)  %in range(0,len(factor)):
    nccn=zeros([length(w) nmodes]);

    M_aero = factor(k)*M_aero0;
    N_aero = factor(k)*N_aero0;
    rd = (3.0.*M_aero.*exp(-4.5.*log(sigma).^2.)./(4.0.*N_aero.*pi.*density)).^(1.0/3.0);
    % print factor[k],rd[0]

    for n=1:length(w)
        Ak=2.0*Mw*zetasa/(Ru*T*rhow);
        alpha=g*(Lv/(eps*cp*T)-1)/(T*Rd);
        es=(100.0*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc));
        gamma=eps*p/es+Lv.^2/(Rv*T.^2*cp);
        bigG=1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)));
        zeta=(2.0/3.0)*Ak*(w(n)*alpha/bigG).^0.5;
        rsmax2=0.0;

        s_cr = zeros([1 nmodes]);
        for i=1:length(nmodes)
            Bk=vantHoff(i)*Mw*density(i)/(massMole(i)*rhow);
            s_cr(i)=(2.0/sqrt(Bk))*(Ak/(3.0*rd(i))).^1.5;
            eta=(w(n)*alpha/bigG).^1.5/(2.0*pi*rhow*gamma*N_aero(i));
            f1=0.5*exp(2.5*(log(sigma(i))).^2);
            f2=1.0+0.25*log(sigma(i));
            rsmax2=rsmax2+(f1*(zeta/eta).^1.5+f2*(s_cr(i)*s_cr(i)/(eta+3.0*zeta)).^.75)/(s_cr(i)*s_cr(i));

        end

        smax=sqrt(1.0/rsmax2);
        for i=1:length(nmodes)
            error_func=1.0-erf(2.0*log(s_cr(i)/smax)/(3.0*sqrt(2.0)*log(sigma(i))));
            nccn(n,i)=0.5*N_aero(i)*error_func;

            % Make sure we don't activate too many...
            nccn(n,i)=min([nccn(n,i) 0.999*N_aero(i)]);
        end

    end
    
    nccn_all(k,:) = nccn; 

    % print rd(i)*(s_cr(i)/smax).^(2./3.)
%    plot(w,nccn(:,1)*10e-6,'color','b','linewidth',2);
%    leg{k}=num2str(factor(k));
%    hold on
    
end

%    legend(leg);    
    figure
    
    w2=[w(1)-dw w w(end)-dw];
    w2=0.5*[w2(2:end) + w2(1:end-1)];
        
    f2=[factor(1)-df factor factor(end)-df];
%    num = factor.*Maero0;
    f2=[factor(1)-df factor factor(end)-df];    
    f2=0.5*[f2(2:end) + f2(1:end-1)];  
    
    dpcolor(w2,f2,nccn_all/1e6); %Converted to per mg (same as per cc for air density=1)
    shading interp
    colorbar
    xlabel('W (m s^{-1})','fontsize',15);
%    ylabel('Nc (cm^{-3})','fontsize',15);
    ylabel('factor','fontsize',15);
    title(['r=' num2str(r_single*1e9) ' nm, rd=' num2str(rd(1)*1e9,'%.0f') ' nm, N_{aero MAX}=' num2str(N_aero0(1)*factor(end)/1e6,'%.0f') ' mg^{-1}'],'fontsize',15);
    
    
    figure
    [temp,iw] = min(abs(w-0.5));
    tot_aero = factor.*N_aero0(1);
    plot(tot_aero/1e6,nccn_all(:,iw) ./ tot_aero');
    set(gca,'xscale','log');
    set(gca,'xlim',[100 10000]);
    set(gca,'ylim',[0 1]);
    xlabel('Number conc. (cm^{-3})');
    ylabel('Activated fraction');
