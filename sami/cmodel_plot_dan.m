function f=cmodel_plot_dan(fig0,n_time,nbin,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry,n_bins)


if (nargin==0)
    fig0=1;
end
fig=fig0;
%
plot_RH_T=1;    % piirret‰‰nkˆ RH ja l‰mpˆtila ajan funktiona
plot_cgas=1;    % piirret‰‰nkˆ kaasujen pitoisuudet ajan funktiona
plot_distrib=1; % piirret‰‰nkˆ hiukkasten kokojakaumat
plot_kohler=1;  % piirret‰‰nkˆ Kˆhler k‰yr‰
plot_final_c=1; % piirret‰‰nkˆ kokojakauma lopussa
%

[ind_r_act,ind_t_act]=find_act(time,r_bin,nbin);
% ind_r_act     Ensimm‰isen aktivoituneen binin indeksi
% ind_t_act     Aktivoitumisen ajanhetken indeksi
%
if ~isempty(ind_r_act)
    % valitaan halutut aktivoituneiden hiukkasten indeksit. Hiukkanen on
    % aktivoitunut kun s‰de on riitt‰v‰n suuri, eli r_act
    r_act=r_bin(n_time,ind_r_act);
    ind_act=find(r_bin(n_time,:)>=r_act);
    ind_noact=find(r_bin(n_time,:)<r_act);
    %
    calc_act_frac(ind_r_act,ind_act,c_bin,nbin,max(RH))
end 
%
if plot_RH_T
    % Piirret‰‰n RH ja l‰mpˆtila samaan kuvaan
    fig=fig+1;
    figure(fig), clf
    [AX,H1,H2]=plotyy(time,RH.*100,time,T);
    xlabel('t (s)')
    ylabel('RH (%)')
    set(get(AX(2),'Ylabel'),'String','T (K)')
    title('RH and temperature as a function of time')
end
%
if plot_cgas
    % piirret‰‰n kaasujen pitoisuudet samaan kuvaan
    fig=fig+1;
    figure(fig), clf
    for i=1:n_gases
        subplot(n_gases,1,i), plot(time,c_gases(:,i),'k-')
        ylabel('c (mol/cm^3)')
        xlabel('t (s)')
        if i==1
            title('Gas concentrations as a function of time')
        end
    end
end
%
if and(~isempty(c_bin),plot_distrib) 
    % piirret‰‰n kokojakauma alussa ja lopussa
    fig=fig+1;
    figure(fig), clf, hold on
    for i=1:length(nbin)
        ii=sum(nbin(1:i-1))+1:1:sum(nbin(1:i));
        loglog(r_dry(ii),c_bin(ii),'b-',r_bin(1,ii),c_bin(ii),'k-',...
            r_bin(length(time),ii),c_bin(ii),'r-')
    end
    set(gca,'XScale','log')
    xlabel('radius (m)')
    ylabel('N (1/cm^3)')
    legend('Dry','Initial','Final')
    title('Aerosol size distribution')
end
%
if and(~isempty(ind_r_act),plot_kohler)
    % Piirret‰‰n pienimm‰n aktivoituneen hiukkasen Kˆhler-k‰yr‰
    fig=fig+1;
    figure(fig), clf
    ind=find(RH>0.998);
    plot(r_bin(ind,ind_r_act).*1.0e6,(RH(ind)-1).*100,'k+-')
    xlabel('radius (\mum)')
    ylabel('S (%)')
    title(['Kˆhler curve for initially dry ',num2str(r_dry(ind_r_act)*1e9),...
            ' nm partile'])
end
%
if and(length(nbin)>1,plot_final_c) 
    % piirret‰‰n kokojakauma lopussa
    distr=1;        % alku 
    distr=n_time;   % lopputulos
    %
    fig=fig+1;
    figure(fig), clf, hold on
    len=round(sum(nbin)/3);
    % jaetaan s‰de sopiviin v‰leihin (ensimm‰inen sis‰lt‰‰ kaikki pienet)
    rr=exp(linspace(log(min(r_bin(distr,:)))-1,log(max(r_bin(distr,:))),len));
    yy=[];
    for i=1:len-1
        ind=find(and(rr(i)<r_bin(distr,:),r_bin(distr,:)<=rr(i+1)));
        yy=[yy,sum(c_bin(ind))];
    end
    xx=(rr(1:len-1)+rr(2:len))./2;
    plot(xx,yy,'k-')
    set(gca,'XScale','log')
    xlabel('radius (m)')
    ylabel('N (1/cm^3)')
    title('Final aerosol size distribution')
end
%
% piirret‰‰n kokojakauman kehitys ajan funktiona. Jos binej‰ voi olla paljon,
% piirret‰‰n niist‰ vain osa.
nn=40;  % piirrett‰vien binien lukum‰‰r‰
figure(fig0), clf
if ~isempty(ind_r_act)
    % valitaan halutut aktivoituneiden hiukkasten indeksit. Hiukkanen on
    % aktivoitunut kun s‰de on riitt‰v‰n suuri
    if nn>=n_bins
        ii_act=1:length(ind_act);
        ii_noact=1:length(ind_noact);
    else
        ii_act=round(linspace(1,length(ind_act),ceil(nn*length(ind_act)/n_bins)));
        ii_noact=round(linspace(1,length(ind_noact),ceil(nn*length(ind_noact)/n_bins)));
    end
    ind_act=ind_act(ii_act);
    ind_noact=ind_noact(ii_noact);
    %
    % ensimm‰inen ja viimeinen aktivoitunut, sek‰ ensimm‰inen ja viimeinen
    % aktivoitumaton bini (s‰de)
    limits=[min(r_bin(1,:)),r_bin(1,ind_r_act-1),...
            r_bin(1,ind_r_act),max(r_bin(1,:))].*1e9;
    %
    a=semilogy(time,r_bin(:,ind_noact),'r-'); hold on
    b=semilogy(time,r_bin(:,ind_act),'k-');
    legend([a(1),b(1)],...
        ['R_p=',num2str(limits(1),'%9.2f'),' - ',num2str(limits(2),'%9.2f'),' nm'],...
        ['R_p=',num2str(limits(3),'%9.2f'),' - ',num2str(limits(4),'%9.2f'),' nm'],2)
else
    % ei aktivoitumista
    ind=round(linspace(1,n_bins,nn));
    semilogy(time,r_bin(:,ind),'k-')
end
title('Aerosol population as a function of time')
ylabel('r (m)')
xlabel('t (s)')
%--------------------------------------------------------------------------


function f=calc_act_frac(ind_r_act,ind_act,c_bin,nbin,Scrit)
% calc_act_frac(ind_r_act,c_bin), funktio laskee aktivoituneiden hiukkasten
% osuuden kaikista hiukkasista
% ind_r_act     ensimm‰inen aktivoitunut bini
% c_bin         hiukkasten pitoisuudet (1/cm^3)
%
if length(nbin)==1  % yksi jakauma
    % Koska tarkka aktivoitumis‰de lˆytyy jostain ensimm‰isen aktivoituneen
    % binin sis‰lt‰, ensimm‰inen aktivoitunut bini on hieman ep‰varma.
    % Oletetaan, ett‰ puolet ensimm‰isest‰ binist‰ todella aktivoituu,
    % jolloin ep‰varmuus on puolet ensimm‰isen binin m‰‰r‰st‰.
    n_act=sum(c_bin(ind_act))-c_bin(ind_r_act)/2; % oletusarvo
    dn_act=c_bin(ind_r_act)/2;          % +/- virhe
    %
    % fix:
    % Tarkka aiktivoitumiss‰de lˆytyy ensimm‰isen aktivoituneen ja
    % aktivoitumattoman binin v‰list‰, joten odotusarvona on, ett‰ koko
    % ensimm‰inen bini aktivoituu. Ep‰tarkkuus on puolet n‰iden binien
    % hiukkaslukum‰‰r‰st‰ (-N(act)/2 ... +N(act-1)/2)
    n_act=sum(c_bin(ind_act));
    dn_act=c_bin(ind_r_act)/2;      % +/- virheen likiarvo
elseif length(nbin)>1  % useita jakaumia
    % Oletetaan, ett‰ ensimm‰inen bini aktivoituu kokonaan
    dn_act=0;
    n_act=sum(c_bin(ind_act));
end
n_tot=sum(c_bin);
%
disp(' ')
disp(['S_max (%)             ',num2str(Scrit*100-100)])
disp(['Hiukkaspitoisuus      ',num2str(n_tot)])
disp(['Aktivoituneita        ',num2str(n_act)])
disp(['Ei aktivoituneita     ',num2str(n_tot-n_act)])
if dn_act==0
    disp(['Aktivoituneiden osuus ',num2str(n_act/n_tot)])
else
    disp(['Aktivoituneiden osuus ',num2str(n_act/n_tot),' +/- ',...
            num2str(dn_act/n_tot)])
end
disp(' ')
%--------------------------------------------------------------------------


function [ind_r_act,ind_t_act]=find_act(time,r_bin,nbin)
% [ind_r_act,ind_t_act]=find_act(time,r_bin), funktio etsii ensimm‰isen 
% aktivoituneen binin indeksin (ind_r_act) ja aktivoitumisen ajanhetken
% indeksen (ind_t_act). Jos aktivoitumista ei tapahdu, palautetaan tyhji‰ 
% matriiseja.
%
% time              aika-vektori
% r_bin(time,bin)   s‰de-matriisi
%
n_act=0;
ind_r_act=[];
%
[n_time,n_bin]=size(r_bin);
for i=1:n_bin-1
    if (r_bin(1,i)/r_bin(1,i+1) > r_bin(n_time,i)/r_bin(n_time,i+1)*10)
        n_act=n_act+1;
        ind_r_act=[ind_r_act,i+1];
    end
end
%
if n_act==1
    % Aktivoituminen tapahtuu silloin kuin ensimm‰inen aktivoitunut bini
    % kasvaa nopeiten
    [aa,ind_t_act]=max(r_bin(2:n_time,ind_r_act)-r_bin(1:n_time-1,ind_r_act));
elseif n_act<=length(nbin)
    % erilliset moodit, jotka voivat aktivoitua eri aikaan ja eri kokoiset
    % hiukkaset aktivoituvat eri aikana. Valitaan aktivoituneiksi
    % hiukkasiksi ne, joiden koko ylitt‰‰ pienimm‰n aktivoituneen hiukkasen
    % loppukoon
    [val,ii]=min(r_bin(n_time,ind_r_act));
    ind_r_act=ind_r_act(ii); ind_t_act=[];
else
    ind_r_act=[]; ind_t_act=[];
end
%--------------------------------------------------------------------------






function [n_time,n_bins,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry]=read_data2
% [time,r_bin,RH,T,c_gases,n_bins,n_gases,n_values]=read_data, funktio
% lukee aikaisemmin lasketun datan tekstitiedostosta.
% 
% Luetaan tekstitiedosto, miss‰ data on muodossa
% a=[t (s), [r_bin (m)], RH, T(K), [c_gas (mol/m^3)]]
%a=textread('C:\Tomi\fortran\cmodel\output\output.out','%s');
a=textread('C:\Tomi\fortran\cmodel_f90\output\output01.out','%s');
% Nyt vektori a sis‰lt‰‰ kaikki luvut, merkit ja merkkinot nerkkijonoksi
% muunnettuna. Selvitet‰‰n binien ja kaasujen lukum‰‰r‰ sek‰ yhden ajanhetken 
% tietojen kokonaispituus
if a{1}=='%'
    % alussa kommenttirivej‰, joissa annettu mm. binien lukum‰‰r‰ ja
    % lukum‰‰r‰pitoisuudet kussakin biniss‰. Esim:
    % 
    % output.dat
    %
    % Number of size bins =   INT
    % Number of gases =       INT
    % Record length =         INT
    % Simulation time =       DBLE
    % Vertical speed =        DBLE
    % Number concentrations (#/cm-3) for size bins
    %  DBLE ...
    % Dry radius (m) for size bins
    %  DBLE ...
    %
    % T‰ss‰ siis kommenttimerkkin‰ on '%'
    nconst=6;   % vakiolukujen lukum‰‰r‰ tiedoston alussa
    k=0;        % luettuja vakioarvoja
    c_bin=[]; nc=0;   % luettuja lukum‰‰r‰pitoisuuksia
    r_dry=[]; dr=0;   % luettuja kuivas‰teit‰
    for i=2:1:1000
        if k<nconst    % kaikki vakioita ei ole viel‰ luettu
            if and(a{i-1}=='=',k==0) % erillisten moodien lukum‰‰r‰
                n_modes=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==1) % binien lukum‰‰r‰ erillisiss‰ moodeissa
                n_bins=str2num(a{i});
                for j=2:n_modes
                    n_bins(j)=str2num(a{i+j-1});
                end
                k=k+1;
            elseif and(a{i-1}=='=',k==2) % kaasujen lukum‰‰r‰
                n_gases=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==3) % tietueen pituus
                rlen=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==4) % tietueen pituus
                stime=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==5) % tietueen pituus
                vspeed=str2num(a{i});
                k=k+1;
            end
        else
            % Kaikki vakiot jo luettu, joten luetaan lukum‰‰r‰pitoisuudet alussa
            % tai kuivas‰teet
            apu=str2num(a{i});
            if isempty(apu)
                % ei ole luku, ei tehd‰ mit‰‰n                
            elseif nc<sum(n_bins)
                % lˆydettiin lukum‰‰r‰pitoisuus
                c_bin=[c_bin,apu];
                nc=nc+1;
            elseif dr<sum(n_bins)
                % lˆydettiin kuivas‰de
                r_dry=[r_dry,apu];
                dr=dr+1;
            else
                % Kaikki pitoisuudet on luettu, joten t‰m‰ luku on
                % ensimm‰inen ajasta riippuva arvo
                a=a(i:length(a));       % ei huomioida kommenttirivej‰
                break
            end
        end
    end
else
    % ei kommenttimerkkej‰, joten tietoja joudutaan etsim‰‰n datan
    % perusteella
    c_bin=[]; r_dry=[];
    [n_bins,n_gases,rlen]=find_bins(a);
end
%
% Muutetaan teksti numeroiksi
arvot=[];
for i=1:rlen:length(a)
    apu=[];
    for j=0:1:rlen-1
        apu=[apu,str2num(a{i+j})];
    end
    arvot=[arvot;apu];
end
%
% T‰st‰ saadaan eri tiedot
time=arvot(:,1);
n_time=length(time);
r_bin=arvot(:,2:sum(n_bins)+1);
RH=arvot(:,sum(n_bins)+2);
T=arvot(:,sum(n_bins)+3);
c_gases=arvot(:,sum(n_bins)+4:rlen);
%--------------------------------------------------------------------------


function [n_bins,n_gases,rlen]=find_bins(vec)
% Funktio etsii binien ja kaasujen lukum‰‰r‰t, kun data on muodossa
%   vec=[t (s), [r_bin(1:n_bins) (m)], RH, T(K), [c_gas(1:n_gases) (mol/cm^3)]]
% rlen on yhden ajanhetken tietojen pituus
%
% vec(t)=[time,r_bin(n_bin),RH,T,c_gas(n_gas)]
for i=2:length(vec)
    aa=str2num(vec{i});
    if (aa>0.1) % RH
        n_bins=i-2
        break
    end
end
%
for i=0:length(vec)
    aa=str2num(vec{n_bins+4+i});
    if (aa>0.01)
        n_gases=i
        break
    end
end
rlen=n_gases+n_bins+3;
%--------------------------------------------------------------------------
