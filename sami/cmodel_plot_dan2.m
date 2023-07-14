


function f=calc_act_frac(ind_r_act,ind_act,c_bin,nbin,Scrit)
% calc_act_frac(ind_r_act,c_bin), funktio laskee aktivoituneiden hiukkasten
% osuuden kaikista hiukkasista
% ind_r_act     ensimm�inen aktivoitunut bini
% c_bin         hiukkasten pitoisuudet (1/cm^3)
%
if length(nbin)==1  % yksi jakauma
    % Koska tarkka aktivoitumis�de l�ytyy jostain ensimm�isen aktivoituneen
    % binin sis�lt�, ensimm�inen aktivoitunut bini on hieman ep�varma.
    % Oletetaan, ett� puolet ensimm�isest� binist� todella aktivoituu,
    % jolloin ep�varmuus on puolet ensimm�isen binin m��r�st�.
    n_act=sum(c_bin(ind_act))-c_bin(ind_r_act)/2; % oletusarvo
    dn_act=c_bin(ind_r_act)/2;          % +/- virhe
    %
    % fix:
    % Tarkka aiktivoitumiss�de l�ytyy ensimm�isen aktivoituneen ja
    % aktivoitumattoman binin v�list�, joten odotusarvona on, ett� koko
    % ensimm�inen bini aktivoituu. Ep�tarkkuus on puolet n�iden binien
    % hiukkaslukum��r�st� (-N(act)/2 ... +N(act-1)/2)
    n_act=sum(c_bin(ind_act));
    dn_act=c_bin(ind_r_act)/2;      % +/- virheen likiarvo
elseif length(nbin)>1  % useita jakaumia
    % Oletetaan, ett� ensimm�inen bini aktivoituu kokonaan
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
% [ind_r_act,ind_t_act]=find_act(time,r_bin), funktio etsii ensimm�isen 
% aktivoituneen binin indeksin (ind_r_act) ja aktivoitumisen ajanhetken
% indeksen (ind_t_act). Jos aktivoitumista ei tapahdu, palautetaan tyhji� 
% matriiseja.
%
% time              aika-vektori
% r_bin(time,bin)   s�de-matriisi
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
    % Aktivoituminen tapahtuu silloin kuin ensimm�inen aktivoitunut bini
    % kasvaa nopeiten
    [aa,ind_t_act]=max(r_bin(2:n_time,ind_r_act)-r_bin(1:n_time-1,ind_r_act));
elseif n_act<=length(nbin)
    % erilliset moodit, jotka voivat aktivoitua eri aikaan ja eri kokoiset
    % hiukkaset aktivoituvat eri aikana. Valitaan aktivoituneiksi
    % hiukkasiksi ne, joiden koko ylitt�� pienimm�n aktivoituneen hiukkasen
    % loppukoon
    [val,ii]=min(r_bin(n_time,ind_r_act));
    ind_r_act=ind_r_act(ii); ind_t_act=[];
else
    ind_r_act=[]; ind_t_act=[];
end
%--------------------------------------------------------------------------


function [n_time,n_bins,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry]=read_data(dirpath,out_folder)
% [time,r_bin,RH,T,c_gases,n_bins,n_gases,n_values]=read_data, funktio
% lukee aikaisemmin lasketun datan tekstitiedostosta.
% 
% Luetaan tekstitiedosto, miss� data on muodossa
% a=[t (s), [r_bin (m)], RH, T(K), [c_gas (mol/m^3)]]
%a=textread('C:\Tomi\fortran\cmodel\output\output.out','%s');
%a=textread('/home/romakkan//joninkoodi/joni1.3/glegg_vanha/output/output01.out','%s');

a=textread([dirpath out_folder 'output001.out'],'%s');


% Nyt vektori a sis�lt�� kaikki luvut, merkit ja merkkinot nerkkijonoksi
% muunnettuna. Selvitet��n binien ja kaasujen lukum��r� sek� yhden ajanhetken 
% tietojen kokonaispituus
if a{1}=='%'
    % alussa kommenttirivej�, joissa annettu mm. binien lukum��r� ja
    % lukum��r�pitoisuudet kussakin biniss�. Esim:
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
    % T�ss� siis kommenttimerkkin� on '%'
    nconst=6;   % vakiolukujen lukum��r� tiedoston alussa
    k=-1;        % luettuja vakioarvoja
    c_bin=[]; nc=0;   % luettuja lukum��r�pitoisuuksia
    r_dry=[]; dr=0;   % luettuja kuivas�teit�
    for i=2:1:1000
        if k<nconst    % kaikki vakioita ei ole viel� luettu
            if and(a{i-1}=='=',k==-1)     % termodynamiikkamallin nimi
                % ei lueta
                k=k+1;
            elseif and(a{i-1}=='=',k==0) % erillisten moodien lukum��r�
                n_modes=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==1) % binien lukum��r� erillisiss� moodeissa
                n_bins=str2num(a{i});
                for j=2:n_modes
                    n_bins(j)=str2num(a{i+j-1});
                end
                k=k+1;
            elseif and(a{i-1}=='=',k==2) % kaasujen lukum��r�
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
            % Kaikki vakiot jo luettu, joten luetaan lukum��r�pitoisuudet alussa
            % tai kuivas�teet
            apu=str2num(a{i});
            if isempty(apu)
                % ei ole luku, ei tehd� mit��n                
            elseif nc<sum(n_bins)
                % l�ydettiin lukum��r�pitoisuus
                c_bin=[c_bin,apu];
                nc=nc+1;
            elseif dr<sum(n_bins)
                % l�ydettiin kuivas�de
                r_dry=[r_dry,apu];
                dr=dr+1;
            else
                % Kaikki pitoisuudet on luettu, joten t�m� luku on
                % ensimm�inen ajasta riippuva arvo
                a=a(i:length(a));       % ei huomioida kommenttirivej�
                break
            end
        end
    end
else
    % ei kommenttimerkkej�, joten tietoja joudutaan etsim��n datan
    % perusteella
    c_bin=[]; r_dry=[];
    [n_bins,n_gases,rlen]=find_bins(a);
end
%
% Muutetaan teksti numeroiksi
arvot=[];
for i=1:rlen:length(a)-1
    arvot=[arvot;str2num(char(a(i:i+rlen-1)))'];
end
%
% T�st� saadaan eri tiedot
time=arvot(:,1);                        n=1;
r_bin=arvot(:,n+1:sum(n_bins)+1);       n=n+sum(n_bins);
RH=arvot(:,n+1);                        n=n+1;         
T=arvot(:,n+1);                         n=n+1;
c_gases=arvot(:,n+1:n+n_gases);         n=n+n_gases;
c_w_aq=arvot(:,n+1:n+sum(n_bins));      n=n+sum(n_bins);
C_hno3=arvot(:,n+1:n+sum(n_bins));      n=n+sum(n_bins);
% ei muuta
n_time=length(time);
%--------------------------------------------------------------------------



function [n_time,n_bins,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry]=read_data2
% [time,r_bin,RH,T,c_gases,n_bins,n_gases,n_values]=read_data, funktio
% lukee aikaisemmin lasketun datan tekstitiedostosta.
% 
% Luetaan tekstitiedosto, miss� data on muodossa
% a=[t (s), [r_bin (m)], RH, T(K), [c_gas (mol/m^3)]]
%a=textread('C:\Tomi\fortran\cmodel\output\output.out','%s');
a=textread('C:\Tomi\fortran\cmodel_f90\output\output01.out','%s');
% Nyt vektori a sis�lt�� kaikki luvut, merkit ja merkkinot nerkkijonoksi
% muunnettuna. Selvitet��n binien ja kaasujen lukum��r� sek� yhden ajanhetken 
% tietojen kokonaispituus
if a{1}=='%'
    % alussa kommenttirivej�, joissa annettu mm. binien lukum��r� ja
    % lukum��r�pitoisuudet kussakin biniss�. Esim:
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
    % T�ss� siis kommenttimerkkin� on '%'
    nconst=6;   % vakiolukujen lukum��r� tiedoston alussa
    k=0;        % luettuja vakioarvoja
    c_bin=[]; nc=0;   % luettuja lukum��r�pitoisuuksia
    r_dry=[]; dr=0;   % luettuja kuivas�teit�
    for i=2:1:1000
        if k<nconst    % kaikki vakioita ei ole viel� luettu
            if and(a{i-1}=='=',k==0) % erillisten moodien lukum��r�
                n_modes=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==1) % binien lukum��r� erillisiss� moodeissa
                n_bins=str2num(a{i});
                for j=2:n_modes
                    n_bins(j)=str2num(a{i+j-1});
                end
                k=k+1;
            elseif and(a{i-1}=='=',k==2) % kaasujen lukum��r�
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
            % Kaikki vakiot jo luettu, joten luetaan lukum��r�pitoisuudet alussa
            % tai kuivas�teet
            apu=str2num(a{i});
            if isempty(apu)
                % ei ole luku, ei tehd� mit��n                
            elseif nc<sum(n_bins)
                % l�ydettiin lukum��r�pitoisuus
                c_bin=[c_bin,apu];
                nc=nc+1;
            elseif dr<sum(n_bins)
                % l�ydettiin kuivas�de
                r_dry=[r_dry,apu];
                dr=dr+1;
            else
                % Kaikki pitoisuudet on luettu, joten t�m� luku on
                % ensimm�inen ajasta riippuva arvo
                a=a(i:length(a));       % ei huomioida kommenttirivej�
                break
            end
        end
    end
else
    % ei kommenttimerkkej�, joten tietoja joudutaan etsim��n datan
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
% T�st� saadaan eri tiedot
time=arvot(:,1);
n_time=length(time);
r_bin=arvot(:,2:sum(n_bins)+1);
RH=arvot(:,sum(n_bins)+2);
T=arvot(:,sum(n_bins)+3);
c_gases=arvot(:,sum(n_bins)+4:rlen);
%--------------------------------------------------------------------------


function [n_bins,n_gases,rlen]=find_bins(vec)
% Funktio etsii binien ja kaasujen lukum��r�t, kun data on muodossa
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
