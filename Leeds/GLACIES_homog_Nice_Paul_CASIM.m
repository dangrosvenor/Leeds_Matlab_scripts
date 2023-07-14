function [Nice,Nice_lim] = GLACIES_homog_Nice_Paul_CASIM(w,P,Tk,Ndrops)
%w in m/s, P in Pa, T in K, Ndrops in per kg

Ra=287.1;             %#spec. gas const of air [J kg-1 K-1]
Rv=461.5;
Lv = 2.501e6;   % Latent heat of vapourization
Lf = 3.34e5; %Latent heat of fusion (freezing)
cp = 1004.67; 
Dv = 0.226e-4;    % diffusivity of water vapour in air
g = 9.81; %check that g is this

d0_homog=50e-6; %m - assumed diameter of frozen ice crystals
%d0_homog=10e-6; %m - assumed diameter of frozen ice crystals

rhoa = density(P,Tk);
%[qsw,esw]=qsatw_rob(Tk,P)
[Ew,qv] = SatVapPress(Tk,'goff','liq',P); %Assume that the vapour MR is at water saturation (kg/kg)
[Ei,qi] = SatVapPress(Tk,'goff','ice',P); %Calculate ice sat vap press (qi not needed)

% alternate limiter based on w - this is explicit w. Will overwrite dnumber
        %qv = qfields(k, i_qv)
        %Tk = th*exner(k)
        ka=((5.69+0.017*(Tk-273.15))*1e-5) * 418.6; 

        cap=1.0; %assume spheres and radius used
        rhoi=200.0;
            %#                   vapour
        min_homog_ni=1e2; %kg-1
        

        %Ei=qisaturation(Tk,pressure(k)/100.)*pressure(k)/(qisaturation(Tk,pressure(k)/100.)+0.62198) %vap press over ice [Pa]
        %Ew=qsaturation(Tk,pressure(k)/100.)*pressure(k)/(qsaturation(Tk,pressure(k)/100.)+0.62198) %vap press over liq [Pa]

        bm=1.0/(qv)+Lv*Lf/(cp*Rv*Tk.^2);
        Ai=1.0/(rhoi*Lf.^2/(ka*Rv*Tk.^2)+rhoi*Rv*Tk/(Ei*Dv));
        B0=4.0*pi*cap*rhoi*Ai/rhoa;
        Bis=bm*B0*(Ew/Ei-1.0);
        aw=g/(Ra*Tk)*(Lv*Ra/(cp*Rv*Tk)-1.0);
         
        %all ql mass onto ni new particles
        %m1=ni.^(4/3)*(dmass*3.0/4.0/pi/rhoi)
        %ustar=(Bis*m1/aw)
        %dnumber=ni
%        dnumberi=max(w(k),0.0)*(aw/Bis).^0.75*((pi*rhoi*4.0/dmass/3.0).^0.25).^0.75
%    try just assuming a starting size

        dt=1; %dummpy value for overall number
            %dnumberi=max(w(k),0.0)*(aw/Bis)/d0_homog  / dt  %convert to a rate
        dnumberi=max(w,0.0)*(aw/Bis)/d0_homog  / dt;  %convert to a rate
        
        %make max just below limit
        dnumber = Ndrops; %limit to 90% of number of droplets.
        dnumberi_lim = min(dnumber*0.90, max(min_homog_ni, dnumberi));
     

Nice = dnumberi*rhoa/1e6; %convert to per cc
Nice_lim = dnumberi_lim*rhoa/1e6;




%Karcher (2016) Fig. 5 values (starting with 100 per cc droplets)
wKar = [0.1 0.17 0.3 0.5 0.8 1 1.7 3 5 8 10];
NiKar = [0.1 0.2 0.6 1.2 3 4 9 27 40 80 90];
figure
set(gcf,'color','w');
plot(wKar,NiKar./100); %Fraction of 100 per cc droplets frozen
hold on
plot(wKar,wKar*0.1,'b--');
xlabel('Updraught speed (m s^{-1})');
ylabel('Fraction of droplets frozen');
legend('Karcher parcel model Fig. 5','w*0.1');

