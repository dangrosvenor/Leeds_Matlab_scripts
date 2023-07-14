function [CAPE,CIN,HLCL,TLCL,PLCL]=plot_tephi_data2(t_dat,p_dat,q_dat,qsat_dat,alt_dat,dotted_flag)
% Plots data on a tephigram
% compute the psuedo adiabat
% compute LCL, CAPE and CIN
% estimate overshoot level by equating areas?
tephi_handle=gca;
global dp;dp=-10.;
p_end=10;
%p_start=110000.;

% PHYSICAL CONSTANTS
global A;A=2.53e11;
global B;B=5.42e3;
global epsilon;epsilon=0.622;
global k_gas;k_gas=0.286;
global latent;latent=2.5e6;
global cp_heat;cp_heat=1005.;

if nargin<6; dotted_flag=0; end %set this flag if not input

rh_dat=100.*q_dat./qsat_dat;
th_dat=t_dat.*(100000./p_dat).^k_gas;
% workout the dew point from humidity
for i=1:length(t_dat)
    if(~isnan(q_dat(i)))
        tdew_dat(i)=fzero(inline(sprintf('q_sat(x,%f)-%.20f;',p_dat(i),q_dat(i))),t_dat(i));
    else
        tdew_dat(i)=nan;
    end
end
size(tdew_dat);
[rr,cc]=size(p_dat);
tdew_dat=reshape(tdew_dat,rr,cc);
%index=find(~isnan(tdew_dat));
th_dew_dat=tdew_dat.*(100000./p_dat).^k_gas;
%th_dew_dat=th_dat

%[tdew_dat_r,th_dat_r]=rotate_coords(tdew_dat,th_dat,pi./4);
%[t_dat_r,th_dat_r]=rotate_coords(t_dat,th_dat,pi./4);
%plot(t_dat_r,th_dat_r,'b','linewidth',2);
%plot(tdew_dat_r,th_dat_r,'b','linewidth',2);

% compute the adiabat
[t_grid,th_grid,p_grid]=moist_adiabat(th_dat(1),p_dat(1),p_end,rh_dat(1));
%[t_grid,th_grid,p_grid]=moist_adiabat(t_dat(1),p_dat(1),p_end,100.*q_dat(2)./qsat_dat(1));
% rotate data
[t_psu_r,th_psu_r]=rotate_coords(t_grid,th_grid,pi./4);
[t_dat_r,th_dat_r]=rotate_coords(t_dat,th_dat,pi./4);
%[tdew_dat_r,th_dat_r2]=rotate_coords(tdew_dat,th_dew_dat,pi./4); % not sure?
[tdew_dat_r,th_dat_r2]=rotate_coords(tdew_dat,th_dew_dat,pi./4);

% plot data

if dotted_flag==0
		plot(t_psu_r,th_psu_r,'r','linewidth',2);
		plot(t_dat_r,th_dat_r,'b','linewidth',2);
		plot(tdew_dat_r,th_dat_r2,'k','linewidth',2);
	


% label psuedo-adiabat and Tdew, T curves
i=find(p_grid<25000);
ii=find(max(p_grid(i))==p_grid(i));
i=i(ii(1));
axes(tephi_handle);
text(t_psu_r(i),th_psu_r(i),'PSEUDO-ADIABAT','color','red','fontsize',20);
% tdew etc
i=find(p_dat<30000);
ii=find(max(p_dat(i))==p_dat(i));
i=i(ii(1));
axes(tephi_handle);
text(t_dat_r(i),th_dat_r(i),'T','color','blue','fontsize',20);
text(tdew_dat_r(i),th_dat_r2(i),'T_d','color','black','fontsize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% table of values
axes('position',[0.13 0.725 0.2 0.2]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);box on;
% workout cape etc
[CAPE,CIN,HLCL,TLCL,PLCL,CAPE_surf,CIN_surf,CAPEP,NAP]=calc_cape(p_dat,t_dat,q_dat,qsat_dat,alt_dat);
text(0.1,0.9,sprintf('Surface CAPE: %.0f J/Kg',CAPE_surf));
text(0.1,0.8,sprintf('CIN: %.0f J/Kg',CIN));
text(0.1,0.7,sprintf('P_{lcl}: %.0f mb',PLCL./100));
text(0.1,0.6,sprintf('T_{lcl}: %.1f C',TLCL));
text(0.1,0.5,sprintf('Alt_{lcl}: %.1f m',HLCL));


% find the integrated area under psuedo-adiabat
[I.pressure,I.integrated,I.dist]=overshoot(p_dat,t_dat,th_dat,p_grid,...
    t_grid,th_grid,100000,p_end);
% find the LNB 
index=find(I.pressure<40000);
inds=find(abs(I.dist(index))==min(abs(I.dist(index))));
LNB=I.pressure(index(inds(1)));
% estimate HNB
HNB=interp1(p_dat,alt_dat,LNB);
text(0.1,0.4,sprintf('P_{lnb}: %.1f mb',LNB./100));
text(0.1,0.3,sprintf('Alt_{lnb}: %.1f m',HNB));
% estmimate tropopause
t_interp=interp1(p_dat,t_dat,p_grid);
alt_interp=interp1(p_dat,alt_dat,p_grid);
%i=find(min(t_interp)==t_interp);

[TTROP,HTROP]=tropopause(t_interp,alt_interp); %finds tropopause height and temp

%text(0.1,0.2,sprintf('T_{trop}: %.1f C',TTROP));
%text(0.1,0.1,sprintf('Alt_{trop}: %.1f m',HTROP));


elseif dotted_flag==2    %if plotting extra lines then don't bother with labelling and new CAPE, etc values
		plot(t_psu_r,th_psu_r,'r--','linewidth',2);
		plot(t_dat_r,th_dat_r,'b--','linewidth',2);
		plot(tdew_dat_r,th_dat_r2,'k--','linewidth',2);
else %if don't want it dotted but want to add to a previous plot (or no CAPE, etc values)
        plot(t_psu_r,th_psu_r,'r-','linewidth',2);
		plot(t_dat_r,th_dat_r,'b-','linewidth',2);
		plot(tdew_dat_r,th_dat_r2,'k-','linewidth',2);
        
end %if dotted_flag==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [p_grid,summing,dist_xy]=overshoot(p_dat,t_dat,th_dat,p_psu,t_psu,th_psu,p_start,p_end)
% calculates the over shoot in the sounding

global dp;

p_grid=p_start:dp:p_end;

% varables on pressure grid
th_dat_grid=interp1(p_dat,th_dat,p_grid);
th_psu_grid=interp1(p_psu,th_psu,p_grid);
t_dat_grid=interp1(p_dat,t_dat,p_grid);
t_psu_grid=interp1(p_psu,t_psu,p_grid);
%figure
%plot(p_grid,th_dat_grid);hold;plot(p_dat,th_dat,'r');
dist_xy=-((th_dat_grid-th_psu_grid)./abs(th_dat_grid-th_psu_grid)).*...
    sqrt((th_dat_grid-th_psu_grid).^2+(t_dat_grid-t_psu_grid).^2);
summing=dist_xy;
for i=2:length(p_grid)
    summing(i)=sum(dist_xy(1:i));
end


function [t_grid,th_grid,p_grid]=moist_adiabat(th_start,p_start,p_end,rhstart)
% compute the moist adiabat
% run vars
global dp;


% PHYSICAL CONSTANTS
global A;
global B;
global epsilon;
global k_gas;
global latent;
global cp_heat;

p=p_start;
T=th_start.*(p./100000).^k_gas;
imax=fix((p_start-p_end)./-dp);


esg=A*exp(-B./T(1)); %es at ground
wg=rhstart./100.*esg.*epsilon./p_start;

% setup grid
t_grid=zeros(1,imax);t_grid(1)=T;
th_grid=zeros(1,imax);th_grid(1)=th_start;
p_grid=zeros(1,imax);p_grid(1)=p_start;
% initialise
t_grid(1)=T;
th_grid(1)=th_start;
p_grid(1)=p;

for i=2:imax
    dwsdp=-epsilon.*A.*exp(-B./T)./p./p;
    dwsdt=B.*epsilon.*A.*exp(-B./T)./T./T./p;
    es=A.*exp(-B./T);
    e=p.*wg./epsilon;
    %e=es; %% saturated adiabat
    if e>es
        dTdp=( (k_gas./p) - (latent.*dwsdp./T./cp_heat) ) ...
        ./( (1./T) + (latent.*dwsdt./T./cp_heat) );
        T=T+dTdp.*dp;
        p=p+dp;
        th=T.*(1e5./p).^k_gas;
        % set grid
    else
        p=p+dp;
        T=th_start.*(p./1e5).^k_gas;  %dry adiabat
        th=th_start;
    end
        
    t_grid(i)=T;
    th_grid(i)=th;
    p_grid(i)=p;
end


function [t,h]=tropopause(t_interp,alt_interp)

hposdtdz=1200;

dtdz=diff(t_interp)./diff(alt_interp);
posdt=find(dtdz>0);
goflag=1;
ipos=1;
while goflag==1
    ihpos=posdt(ipos);
    hpos=alt_interp(ihpos);
	ihpos2=findheight(alt_interp,hpos+hposdtdz); %find point hposdtdz km above the positive lapse rate
    ipos2=ipos+ihpos2-ihpos; %index in posdt if all points hposdtdz km above also have positive dtdz
    if ipos2>length(posdt) %if have reached end of positive lapse rate points
        goflag=0;
    else
        diffipos=diff(posdt(ipos:ipos2)); %all diffs should be one if all points positive
        idiffgt1=find(diffipos>1);
        %match=posdt(ipos:ipos2)==[hpos:hpos2];
        if length(idiffgt1)==0 %if differences in posdt indices for a hposdtdz km jump are all one - i.e. dtdz positive for 1.2km
            goflag=0; %stop search - altitude we need is hpos
        else        
            ipos=ipos+idiffgt1(1); %start searching again from height where lapse rate not positive
        end
    end
end
%TTROP=t_interp(min(i));
%HTROP=alt_interp(min(i));
%TTROP=t_interp(itrop);
%HTROP=alt_interp(itrop);
t=t_interp(posdt(ipos))-273.15;
h=hpos;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

