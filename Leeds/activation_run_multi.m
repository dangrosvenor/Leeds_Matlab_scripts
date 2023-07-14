clear nccn_all var_act 

dw = 0.1;
w = [dw:dw:1]; % m/s

%dw = 0.1;
w = [dw:0.5:10]; % m/s

%dw = 0.5;
%w = [dw:dw:20]; % m/s

df=0.1;
factor = [df:df:1];
factor = 10.^[-2:0.25:0];
factor = 10.^[-2:0.25:2];

%df=1;
%factor = [df:df:10];

var_act.N_type = N_type;
var_act.M_aero0 = M_aero0; % Total aerosol mass mixing ratio kg/kg
var_act.density = density; % aerosol density. E.g. ammonium sulphate = 1777 kg/m3 
var_act.sig = sig; %width of lognormal size dist
var_act.w = w; %updraft in m/s
var_act.T = T; %temperature in K
var_act.P = P; %Air pressure in Pa

for k=1:length(factor)  %in range(0,len(factor)):
    %    nccn=zeros([length(w) nmodes]);

    switch N_type
        case 'prescribed number'
            %Vary N_aero, but keep M_aero0 the same
            var_act.M_aero0 = M_aero0;
            var_act.N_aero0 = factor(k)*N_aero0
            % N_aero0 set outside
            tit_str = 'for varying number (constant mass)';

        case 'prescribed single aerosol mass'
            %Vary M_aero and keep r_single (translates to mass for monomodal dist as used in CASIM code to get number) the same
            %This single mass is then used to get he number, so number will
            %scale with the mass (i.e. by factor(k)).
            var_act.M_aero0 = factor(k)*M_aero0;
            %N_aero = 1; %dummy value - will be calculated in the function based on M_aero and r_single (size)
            tit_str = 'for varying mass and number (constant size)';

        case 'prescribed median radius (rd)'
            var_act.rd = r_single;
            var_act.M_aero0 = factor(k)*M_aero0;
            %N_aero0 will be calculated based on the mass and size
            tit_str = 'for varying mass and number (constant size)';
            
       case 'prescribed median radius (rd) vary accum'
            var_act.rd = r_single;
            var_act.M_aero0(1) = M_aero0(1);
            var_act.M_aero0(2) = factor(k)*M_aero0(2);            
            %N_aero0 will be calculated based on the mass and size
            tit_str = 'for varying mass and number (constant size)';     
    end
    
    for iroutine=1:length(act_routines)
        var_act.act_routine = act_routines{iroutine};

        %N.B. - think that the Nenes function needs w as cm/s, so multiplying
        %by 100 here
%        [nccn_all_Nenes(k,:),Dg,N_aero_out,r_single] = nenes_run_func(N_aero,M_aero,w*100,r_single,N_type);

        [nccn_all{iroutine}(k,:),Dg,N_aero_out,r_single_out] = activation_run_func(var_act);

    end

end

N_aero_out_Aitken=0;

switch N_type
    case 'prescribed median radius (rd) vary accum'
        N_aero_out_save = N_aero_out;
        N_aero_out = N_aero_out(2);
        N_aero_out_Aitken = N_aero_out_save(1);
        M_aero0_save = M_aero0;
        M_aero0 = M_aero0(2);
end

%% Make plots - pcolor plot first

for iroutine=1:length(act_routines)
    
%    legend(leg);    
    figure
    
    w2=[w(1)-dw w w(end)+dw];
    w2=0.5*[w2(2:end) + w2(1:end-1)];
        

%    num = factor.*Maero0;
    f2=[factor(1)-df factor factor(end)+df];    
    f2=0.5*[f2(2:end) + f2(1:end-1)];  
    
    dpcolor(w2,f2,nccn_all{iroutine}/1e6); %Converted to per mg (same as per cc for air density=1)
    shading flat %Note shading interp messes up the axes and gives a blank row and column...
    colorbar
    xlabel('W (m s^{-1})','fontsize',15);
%    ylabel('Nc (cm^{-3})','fontsize',15);
    ylabel(['factor ' tit_str],'fontsize',15);
    title([act_routines{iroutine} ':- rdmin=' num2str(Dg/2*1e9,'%.0f') ' nm, N_{aero MAX}=' num2str(N_aero_out(1)/1e6,'%.0f') ' mg^{-1}'],'fontsize',15);
        
    
end

w_find = 5;
w_find = 0.5;
w_find = 0.1;

ifig=0;
clear xsave_act ysave_act
%% Activated fraction plot
figure
ifig=ifig+1;
cols={'b','r'};
marks={'^','s'};
for iroutine=1:length(act_routines)    
    [temp,iw] = min(abs(w-w_find));
    wval=w(iw);
    tot_aero = factor.*N_aero_out;
    plot(tot_aero/1e6,nccn_all{iroutine}(:,iw) ./ tot_aero','marker',marks{iroutine},'color',cols{iroutine});
    xsave_act{ifig}(iroutine).x = tot_aero/1e6;
    ysave_act{ifig}(iroutine).y =  nccn_all{iroutine}(:,iw) ./ tot_aero';  
    hold on
    set(gca,'xscale','log');
    set(gca,'xlim',[100 10000]);
    set(gca,'ylim',[0 1]);
    xlabel('Aerosol number conc. (cm^{-3})');
    ylabel('Activated fraction');   
    title_str = ['w=' num2str(wval) ' m s^{-1}, median radius=' num2str(1e9*Dg/2) ' nm, sig=' num2str(sig)];
    title(title_str);
end
legend(act_routines);


%% Nd vs aerosol number plot
figure
ifig=ifig+1;
cols={'b','r'};
for iroutine=1:length(act_routines)
    [temp,iw] = min(abs(w-w_find));
    wval=w(iw);
    tot_aero = factor.*N_aero_out;
%    plot(tot_aero/1e6,nccn_all{iroutine}(:,iw) ./ tot_aero','color',cols{iroutine});
    plot(tot_aero/1e6,nccn_all{iroutine}(:,iw)/1e6,'marker',marks{iroutine},'color',cols{iroutine});    
    xsave_act{ifig}(iroutine).x = tot_aero/1e6;
    ysave_act{ifig}(iroutine).y = nccn_all{iroutine}(:,iw)/1e6;     
    hold on
    set(gca,'xscale','log');
    set(gca,'xlim',[100 11000]);
    set(gca,'ylim',[0 1500]);
    xlabel('Aerosol number conc. (cm^{-3})');
    ylabel('Activated CCN conc. (cm^{-3})');    
    title(['w=' num2str(wval) ' m s^{-1}, median radius=' num2str(1e9*Dg/2) ' nm, sig=' num2str(sig)]);
end
legend(act_routines);


%% Nd vs aerosol mass plot
figure
ifig=ifig+1;
cols={'b','r'};
for iroutine=1:length(act_routines)
    [temp,iw] = min(abs(w-w_find));
    wval=w(iw);
    tot_aero = factor.*M_aero0;
%    plot(tot_aero/1e6,nccn_all{iroutine}(:,iw) ./ tot_aero','color',cols{iroutine});
    plot(tot_aero,nccn_all{iroutine}(:,iw)/1e6,'marker',marks{iroutine},'color',cols{iroutine}); 
    xsave_act{ifig}(iroutine).x = tot_aero;
    ysave_act{ifig}(iroutine).y = nccn_all{iroutine}(:,iw)/1e6; 
    hold on
    set(gca,'xscale','log');
    set(gca,'xlim',[1e-10 1e-7]);
    set(gca,'ylim',[0 1500]);
    xlabel('Aerosol mass MR. (kg kg^{-1})');
    ylabel('Activated CCN conc. (cm^{-3})');    
    title(['w=' num2str(wval) ' m s^{-1}, median radius=' num2str(1e9*Dg/2) ' nm, sig=' num2str(sig)]);
end
legend(act_routines);




%% Nd vs updraft plot extra plot with Twomey curve too for specific aerosol
%% concentration
aero_find = 200e6; % per m3
aero_find = 600e6; % per m3
aero_find = 1000e6; % per m3

figure
cols={'b','r'};
for iroutine=1:length(act_routines)
    tot_aero = factor.*N_aero_out;
    [temp,iw] = min(abs(tot_aero-aero_find));   
    Naero = tot_aero(iw); %For use in Twomey approx below
    title_str = ['N_{aero}=' num2str(Naero/1e6) ' cm^{-3}, median radius=' num2str(1e9*Dg/2) ' nm, sig=' num2str(sig)];
%    plot(tot_aero/1e6,nccn_all{iroutine}(:,iw) ./ tot_aero','color',cols{iroutine});
    plot(w,nccn_all{iroutine}(iw,:)/1e6,'marker',marks{iroutine},'color',cols{iroutine});    
    hold on
    set(gca,'xscale','log');
%    set(gca,'xlim',[100 10000]);
%    set(gca,'ylim',[0 1500]);
    xlabel('Updraft (m s^{-1})');
    ylabel('Activated CCN conc. (cm^{-3})');    
end

leg = act_routines;
ileg=length(leg)+1;

%Add on the Aitken mode for the Twomey calc
Naero = Naero + N_aero_out_Aitken;

%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
k=0.4; Wmax=50;
N = Twomey_Nd_approx_using_Wmax(Naero,k,Wmax,w);
frac50 = N ./Naero;   
plot(w,N/1e6,'b--','linewidth',2);
hold on
leg{ileg} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero/1e6)]; ileg=ileg+1;


%Wmax=16m/s, Naero=1000;
%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
k=0.4; Wmax=16;
N = Twomey_Nd_approx_using_Wmax(Naero,k,Wmax,w);
frac16 = N ./Naero;   
plot(w,N/1e6,'r--','linewidth',2);
hold on
leg{ileg} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero/1e6)]; ileg=ileg+1;



%Wmax=16m/s, Naero=1000;
%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
k=0.4; Wmax=3;
N = Twomey_Nd_approx_using_Wmax(Naero,k,Wmax,w);
frac3 = N ./Naero;
plot(w,N/1e6,'k--','linewidth',2);
hold on
leg{ileg} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero/1e6)]; ileg=ileg+1;


%
title(title_str);
legend(leg);
set(gca,'xlim',[w(1) w(end)]);


%% As above, but for activated fraction
figure
cols={'b','r'};
for iroutine=1:length(act_routines)
    tot_aero = factor.*N_aero_out;
    [temp,iw] = min(abs(tot_aero-aero_find));   
    Naero = tot_aero(iw); %For use in Twomey approx below
%    plot(tot_aero/1e6,nccn_all{iroutine}(:,iw) ./ tot_aero','color',cols{iroutine});
    plot(w,nccn_all{iroutine}(iw,:)./Naero,'marker',marks{iroutine},'color',cols{iroutine});    
    hold on
    set(gca,'xscale','log');
%    set(gca,'xlim',[100 10000]);
%    set(gca,'ylim',[0 1500]);
    xlabel('Updraft (m s^{-1})');
    ylabel('Activated fraction');    
end

leg = act_routines;
ileg=length(leg)+1;

k=0.4; Wmax=50;
hold on
plot(w,frac50,'b--','linewidth',2);
leg{ileg} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero/1e6)]; ileg=ileg+1;


%Wmax=16m/s, Naero=1000;
%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
k=0.4; Wmax=16;
plot(w,frac16,'r--','linewidth',2);
hold on
leg{ileg} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero/1e6)]; ileg=ileg+1;



%Wmax=16m/s, Naero=1000;
%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
k=0.4; Wmax=3;
plot(w,frac3,'k--','linewidth',2);
hold on
leg{ileg} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero/1e6)]; ileg=ileg+1;


%
title(title_str);
legend(leg);
set(gca,'xlim',[w(1) w(end)]);