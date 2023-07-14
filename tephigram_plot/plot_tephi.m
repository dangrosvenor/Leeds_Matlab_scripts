function [TP,h]=plot_tephi(low_t,high_t,t_min,p_min)
% plots lines and labels for a tephigram
global k_gas;
figure('name','tephigram');
% compute tephigram
TP=tephigram;

% get rotated limits%%%%%%%%%%%%%%%%%%%%%
low_th=low_t.*(100000./105000).^k_gas;
high_thy=high_t.*(100000./5000).^k_gas;
[low_t_r,low_th_r]=rotate_coords(low_t,low_th,pi./4);% x limit1
[high_t_rx,high_th_rx]=rotate_coords(high_t,high_t,pi./4);%x_limit2
%[high_t_rx,high_th_rx]=rotate_coords(high_t,high_thy,pi./4);%x_limit2
% max y depends on p_min and t_min
high_th=t_min.*(100000./p_min).^k_gas;
[high_t_ry,high_th_ry]=rotate_coords(t_min,high_th,pi./4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(TP.t_piso_r',TP.th_piso_r','--g');hold on
for i=1:length(TP.p_iso) % label for p
    % interpolate to find the position to place the label in page coords cf max t_r
    x=high_t_rx-7;

    y=interp1(TP.t_piso_r(i,:),TP.th_piso_r(i,:),x);
    %[x,y]=rotate_coords(TP.t_dry(1,i),273+50,pi./4);
    %if(x>=low_t_r & x<high_t_rx & y>= low_th_r & y< high_th_ry)
    h(1).h(i)=text(x,y,sprintf('%.0f mb',TP.p_iso(i)./100));
        %end
end
set(h(1).h,'color','green');
  
plot(TP.t_dry_r,TP.th_dry_r,'k');
for i=1:length(TP.t_dry(1,:)) % label for t
    [x,y]=rotate_coords(TP.t_dry(1,i),273+50,pi./4);
    if(x>=low_t_r & x<high_t_rx & y>= low_th_r & y< high_th_ry)
        h(2).h(i)=text(x,y,sprintf('T = %.0f',TP.t_dry(1,i)-273),'rotation',45);
    else
        h(2).h(i)=nan;
    end
end
plot(TP.t_dry_r',TP.th_dry_r','k');
for i=1:length(TP.th_dry(:,1)) % label for theta
    [x,y]=rotate_coords(273-20,TP.th_dry(i,1),pi./4);
    if(x>=low_t_r & x<high_t_rx & y>= low_th_r & y< high_th_ry)
        h(3).h(i)=text(x,y,sprintf('\\theta = %.0f',TP.th_dry(i,1)-273),'rotation',-45);
    else
        h(3).h(i)=nan;
    end
end
h(4).h=plot(TP.t_psu_r',TP.th_psu_r','--b');
%for i=1:length(TP.th_psu(:,1)) % label for theta
%    [x,y]=rotate_coords(TP.t_psu(i,1),TP.th_psu(i,1),pi./4);
%    if(x>=low_t_r & x<high_t_rx & y>= low_th_r & y< high_th_ry)
%        text(x,y,sprintf('{\\theta}_w = %.0f',TP.th_psu(i,1)-273),'rotation',-90);
%    end
%end

plot(TP.t_smr_r',TP.th_smr_r','-.m');
for i=1:length(TP.smr_iso) % label for p
    % interpolate to find the position to place the label in page coords cf max t_r
    % match the position with axes
    % interpolate to find x at ymin, if greater than xmin use ymin as y cord and interpolate to find x
    % else use xmin as x cord and interpolate to find y
    x_test=interp1(TP.th_smr_r(i,:),TP.t_smr_r(i,:),low_th_r);
    if(x_test>=low_t_r)
        y=low_th_r;
        x=x_test;
    else
        x=low_t_r;
        y=interp1(TP.t_smr_r(i,:),TP.th_smr_r(i,:),low_t_r);
    end
    h(5).h(i)=text(x,y,sprintf('%.2f',TP.smr_iso(i).*1000));
end
set(h(5).h,'color','magenta');

% set rotated limits
set(gca,'xlim',[low_t_r,high_t_rx]);
set(gca,'ylim',[low_th_r,high_th_ry]);
% clear x,y ticks
set(gca,'xtick',[]);
set(gca,'ytick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TP=tephigram
% FUNCTION COMPUTES DRY ADIABATS, PSUEDO-ADIABATS,
% SMR LINES, ISOBARS AND ISOTHERMS OF A TEPHIGRAM

% RUN CONSTANTS
global dp;dp=-10.;
p_end=100;
%p_start=110000.;

% PHYSICAL CONSTANTS
global A;A=2.53e11;
global B;B=5.42e3;
global epsilon;epsilon=0.622;
global k_gas;k_gas=0.286;
global latent;latent=2.5e6;
global cp_heat;cp_heat=1005.;

% limits on th and t at 1000mb
low_th=-40+273;
high_th=40+273;
low_t=-40+273;
high_t=40+273;

% DRY ADIABATS
%th_min=-50+273;
%th_max=200+273;
th_dry=(-200:10:200)+273;
%t_min=-100+273;
%t_max=50+273;
t_dry=(-200:10:50)+273;
[TP.t_dry,TP.th_dry]=meshgrid(t_dry,th_dry);

% PSUEDO - ADIABATS
disp('Computing psuedo-adiabats...');
p_start=100000.;
th_moist=[-2:2:45]+273;
for i=1:length(th_moist) % all adiabats
    [t_grid(i,:),th_grid(i,:),p_grid(i,:)]=moist_adiabat(th_moist(i),p_start,p_end);
end 

TP.t_psu=t_grid;
TP.th_psu=th_grid;
TP.p_psu=p_grid;
disp('Done');

% isobars
t_min=-190;
t_max=330;
p_iso=[1050 1000 900 800 700 600 500 400 350 300 250 200 150 100 75 50].*100;
%p_iso=105000;
for i=1:length(p_iso)
    TP.t_piso(i,:)=(t_min:5:t_max);
    TP.th_piso(i,:)=TP.t_piso(i,:).*(100000./p_iso(i)).^k_gas;
end
TP.p_iso=p_iso;

% smr
smr_iso=[48 36 28 20 16 12 9 6 4 2.5 1.5 0.8 0.4 0.2 0.1 0.05 0.03 0.02 0.01]./1000;
TP.smr_iso=smr_iso;
disp('Calculating smr isolines');
p_start=110000;
jmax=fix((p_start-p_end)./-dp);
% allocate
TP.t_smr=zeros(length(smr_iso),jmax);
TP.th_smr=zeros(length(smr_iso),jmax);
for i=1:length(smr_iso)
    p=p_start;
    ws=smr_iso(i);
    for j=1:jmax
        %ws=epsilon.*A.*exp(-B./t)./(100000-A.*exp(-B./t));
        %log(((epsilon+ws).*A)./(100000.*ws))=(B./t)
        TP.t_smr(i,j)=B./log(((epsilon+ws).*A)./(p.*ws));
        TP.th_smr(i,j)=TP.t_smr(i,j).*(100000./p).^k_gas;
        p=p+dp;
    end 
end
disp('Done');
% also compute rotated coords
[TP.t_psu_r,TP.th_psu_r]=rotate_coords(TP.t_psu,TP.th_psu,pi./4);
[TP.t_dry_r,TP.th_dry_r]=rotate_coords(TP.t_dry,TP.th_dry,pi./4);
[TP.t_piso_r,TP.th_piso_r]=rotate_coords(TP.t_piso,TP.th_piso,pi./4);
[TP.t_smr_r,TP.th_smr_r]=rotate_coords(TP.t_smr,TP.th_smr,pi./4);


function [t_grid,th_grid,p_grid]=moist_adiabat(th_start,p_start,p_end)
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
