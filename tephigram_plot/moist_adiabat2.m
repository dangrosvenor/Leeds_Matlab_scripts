function [t_grid,th_grid,p_grid]=moist_adiabat2(th_start,p_start,p_end)
% compute the moist adiabat 
% function [t_grid,th_grid,p_grid]=moist_adiabat2(th_start,p_start,p_end)
% th_start = starting potemp
% p_start=start pressure (Pa)
% p_end = end pressure (Pa)

% PHYSICAL CONSTANTS
A=2.53e11;
B=5.42e3;
epsilon=0.622;
k_gas=0.286;
latent=2.5e6;
cp_heat=1005.;

dp=sign(p_end-p_start)*10.;

p=p_start;
T=th_start.*(p./100000).^k_gas;
imax=fix((p_start-p_end)./-dp);

% setup grid
t_grid=zeros(1,imax);t_grid(1)=T;
th_grid=zeros(1,imax);th_grid(1)=th_start;
p_grid=zeros(1,imax);p_grid(1)=p_start;

for i=2:imax
    dwsdp=-epsilon.*A.*exp(-B./T)./p./p;
    dwsdt=B.*epsilon.*A.*exp(-B./T)./T./T./p;
    es=A.*exp(-B./T);
    %e=p.*wg./epsilon;
    e=es; %% saturated adiabat
    dTdp=( (k_gas./p) - (latent.*dwsdp./T./cp_heat) ) ...
        ./( (1./T) + (latent.*dwsdt./T./cp_heat) );
    T=T+dTdp.*dp;
    p=p+dp;
    th=T.*(1e5./p).^k_gas;
    % set grid
    t_grid(i)=T;
    th_grid(i)=th;
    p_grid(i)=p;
end
