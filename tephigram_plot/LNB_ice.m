function [t_grid,th_grid,p_grid,lnb]=lnb_ice(th_start,p_start,penv,tenv,henv,qi,qv)
% compute the moist adiabat

% PHYSICAL CONSTANTS
A=2.53e11;
B=5.42e3;
epsilon=0.622;
k_gas=0.286;
latent=2.5e6;
cp_heat=1005.;
dpmag=10; %magnitude of the pressure step

p=p_start;
t_start=th_start.*(p./100000).^k_gas;
T=t_start;

qtot=qi+qv;
etot=p.*qtot./(qtot+epsilon); %vapour pressure of total water

Tenv_start=interp1(penv,tenv,p_start); %find environmental temperature at thermal pressure
dp=sign(Tenv_start-T)*dpmag;

% imax=fix((p_start-p_end)./-dp);

imax=1;
% setup grid
t_grid=T;
th_grid=th_start;
p_grid=p_start;

sig=sign(T-Tenv_start);
signew=sig;
finish=0;

ifin=[];
fin=zeros(size(sig));
answer=zeros(size(sig));
answer=answer*NaN;
sfinish=[];
nits=0;

while (sfinish~=prod(size(sig)) & nits<2000)
    nits=nits+1;
    dwsdp=-epsilon.*A.*exp(-B./T)./p./p;
    
    dwsdt=B.*epsilon.*A.*exp(-B./T)./T./T./p;
    es=A.*exp(-B./T);
    
    isat=find(etot>=es);
    Tsat=T(isat);
    psat=p(isat);
    %e=p.*wg./epsilon;
    e=es; %% saturated adiabat
    th=T.*(1e5./p).^k_gas; %update th for saturated adiabat but not for dry (constant potemp)
    
    dTdp(isat)=( (k_gas./psat) - (latent.*dwsdp(isat)./Tsat./cp_heat) ) ...
        ./( (1./Tsat) + (latent.*dwsdt(isat)./Tsat./cp_heat) );
    
    idry=find(etot<es);
    dTdp(idry)=th(idry).*k_gas.*1e-5.^k_gas.*p(idry).^(k_gas-1); %constant potemp
    T=T+dTdp.*dp;
    p=p+dp;
    
    % set grid
    t_grid=T;
    th_grid=th;
    p_grid=p;
    Tenv=interp1(penv,tenv,p_grid); %find new environmental temp for new thermal temp
    
    T(ifin)=t_start(ifin); %reset ones that have already been found to same as at start so that are not stored again
    Tenv(ifin)=Tenv_start(ifin); %reset so that signew will be same as original one
    
    signew=sign(T-Tenv);
    ifin=find(signew~=sig);
    fin(ifin)=1;
    answer(ifin)=p_grid(ifin); %store final pressure for those already finished
    
    inonzero=find(fin~=0);
    sfinish=prod(size(inonzero));

    
end

if nits==2001
    fprintf(1,'\n nits == 2001');
end


lnb=interp1(penv,henv,answer); %find height of LNB pressure level
