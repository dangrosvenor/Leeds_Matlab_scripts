function [t_grid,th_grid,p_grid,lnb]=lnb_ice3(th_start,z_start,penv,tenv,henv,qi,qv)
% compute the moist adiabat
%this one uses the lapse rate definition of the pseudoadiabat - prob best using lnb_ice2
% PHYSICAL CONSTANTS
A=2.53e11;
B=5.42e3;
epsilon=0.622;
Rv=461; %individual gas constant for water vapour. epsilon=Rd/Rv
Rd=epsilon*Rv;
k_gas=0.286;
latent=2.5e6;
latent=2.834e6;
cp_heat=1005.;

dzmag=10; %magnitude of the height step

p_start=interp1(henv,penv,z_start);
t_start=th_start.*(p_start./100000).^k_gas;

p=p_start;
T=t_start;




qtot=qi+qv;
etot=p.*qtot./(qtot+epsilon); %vapour pressure of total water

Tenv_start=interp1(penv,tenv,p_start); %find environmental temperature at thermal pressure

dz=-sign(Tenv_start-T)*dzmag;

% imax=fix((p_start-p_end)./-dp);

imax=1;
% setup grid
t_grid=T;
th_grid=th_start;
p_grid=p_start;
th=th_start;
z=z_start;

sig=sign(T-Tenv_start);
signew=sig;
finish=0;

ifin=[];
fin=zeros(size(sig));
answer=zeros(size(sig));
answer=answer*NaN;
sfinish=[];
nits=0;

dtdz_dry=-9.81/cp_heat;

nitsmax=2000;
nits=0;

while (sfinish~=prod(size(sig)) & nits<nitsmax)
    nits=nits+1;
    %dwsdp=-epsilon.*A.*exp(-B./T)./p./p;
    dwsdt=numdiff(@ws,T,p,1e-4,'x');
    %dwsdt=B.*epsilon.*A.*exp(-B./T)./T./T./p;
    dwsdp=numdiff(@ws,T,p,1e-2,'y');
    
    %es=A.*exp(-B./T);
    es=SatVapPress(T,'lem','ice',p);
    
    isat=find(etot>=es);
    Tsat=T(isat);
    psat=p(isat);

    th(isat)=T(isat).*(1e5./p(isat)).^k_gas; %update th for saturated adiabat but not for dry (constant potemp)
    

	wsat=ws(T,p);
	dTdz(isat)=dtdz_dry.*(1 + latent*wsat(isat)/Rd/T(isat))./(1 + latent*latent*epsilon*wsat(isat)./Rd./cp_heat./T(isat)./T(isat));
    
    idry=find(etot<es);
    dTdz(idry)=dtdz_dry;
    
    T=T+dTdz.*dz;
    p=interp1(henv,penv,z);
    z=z+dz;
    
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

if nits==nitsmax
    fprintf(1,'\n nits == nitsmax');
end


lnb=interp1(penv,henv,answer); %find height of LNB pressure level
