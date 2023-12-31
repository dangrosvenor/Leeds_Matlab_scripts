function [n_time,n_bins,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry]=read_data(dirpath,out_folder)
% [time,r_bin,RH,T,c_gases,n_bins,n_gases,n_values]=read_data, funktio
% lukee aikaisemmin lasketun datan tekstitiedostosta.
% 
% Luetaan tekstitiedosto, missä data on muodossa
% a=[t (s), [r_bin (m)], RH, T(K), [c_gas (mol/m^3)]]
%a=textread('C:\Tomi\fortran\cmodel\output\output.out','%s');
%a=textread('/home/romakkan//joninkoodi/joni1.3/glegg_vanha/output/output01.out','%s');

a=textread([dirpath out_folder 'output001.out'],'%s');


% Nyt vektori a sisältää kaikki luvut, merkit ja merkkinot nerkkijonoksi
% muunnettuna. Selvitetään binien ja kaasujen lukumäärä sekä yhden ajanhetken 
% tietojen kokonaispituus
if a{1}=='%'
    % alussa kommenttirivejä, joissa annettu mm. binien lukumäärä ja
    % lukumääräpitoisuudet kussakin binissä. Esim:
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
    % Tässä siis kommenttimerkkinä on '%'
    nconst=6;   % vakiolukujen lukumäärä tiedoston alussa
    k=-1;        % luettuja vakioarvoja
    c_bin=[]; nc=0;   % luettuja lukumääräpitoisuuksia
    r_dry=[]; dr=0;   % luettuja kuivasäteitä
    for i=2:1:1000
        if k<nconst    % kaikki vakioita ei ole vielä luettu
            if and(a{i-1}=='=',k==-1)     % termodynamiikkamallin nimi
                % ei lueta
                k=k+1;
            elseif and(a{i-1}=='=',k==0) % erillisten moodien lukumäärä
                n_modes=str2num(a{i});
                k=k+1;
            elseif and(a{i-1}=='=',k==1) % binien lukumäärä erillisissä moodeissa
                n_bins=str2num(a{i});
                for j=2:n_modes
                    n_bins(j)=str2num(a{i+j-1});
                end
                k=k+1;
            elseif and(a{i-1}=='=',k==2) % kaasujen lukumäärä
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
            % Kaikki vakiot jo luettu, joten luetaan lukumääräpitoisuudet alussa
            % tai kuivasäteet
            apu=str2num(a{i});
            if isempty(apu)
                % ei ole luku, ei tehdä mitään                
            elseif nc<sum(n_bins)
                % löydettiin lukumääräpitoisuus
                c_bin=[c_bin,apu];
                nc=nc+1;
            elseif dr<sum(n_bins)
                % löydettiin kuivasäde
                r_dry=[r_dry,apu];
                dr=dr+1;
            else
                % Kaikki pitoisuudet on luettu, joten tämä luku on
                % ensimmäinen ajasta riippuva arvo
                a=a(i:length(a));       % ei huomioida kommenttirivejä
                break
            end
        end
    end
else
    % ei kommenttimerkkejä, joten tietoja joudutaan etsimään datan
    % perusteella
    c_bin=[]; r_dry=[];
    [n_bins,n_gases,rlen]=find_bins(a);
end
%
% Muutetaan teksti numeroiksi
arvot=[];
for i=1:rlen:length(a)-rlen+1
    arvot=[arvot;str2num(char(a(i:i+rlen-1)))'];
end
%
% Tästä saadaan eri tiedot
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