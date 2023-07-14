function I=plot_tephi_data(t_dat,p_dat,q_dat,qsat_dat)
% Plots data on a tephigram
% compute the psuedo adiabat
% comptute LCL, CAPE and CIN
% estimate overshoot level by equating areas?

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
tdew_dat=tdew_dat';
index=find(~isnan(tdew_dat));
th_dew_dat=tdew_dat(index).*(100000./p_dat(index)).^k_gas;

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
[tdew_dat_r,th_dat_r2]=rotate_coords(tdew_dat(index),th_dat(index),pi./4);

% plot data
plot(t_psu_r,th_psu_r,'r','linewidth',2);
plot(t_dat_r,th_dat_r,'b','linewidth',2);
plot(tdew_dat_r,th_dat_r2,'k','linewidth',2);

% find the integrated area under psuedo-adiabat
[I.pressure,I.integrated]=overshoot(p_dat,t_dat,th_dat,p_grid,...
    t_grid,th_grid,100000,p_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_grid,summing]=overshoot(p_dat,t_dat,th_dat,p_psu,t_psu,th_psu,p_start,p_end)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

