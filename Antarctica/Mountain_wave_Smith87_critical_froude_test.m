gd=1;
U=20;
N=0.02;

hhat=[0:0.01:1];

%continuous stratification case (Fig. 15b)
del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
H0_hat_strat=hhat - del_hat + acos(hhat./del_hat);
H0_crit=U/N*H0_hat_strat; %
F0_crit=1 ./ H0_hat_strat; %=U/NH=F0=1/H0_hat for Fig 15b


%single layer case (Fig. 15a)

H0_crit2 = [0:10:30000];
F0_crit2=U./sqrt(gd.*H0_crit2);
F0_crit3 = U/N./H0_crit2;
h_crit2 = H0_crit2 .* ( 1 + 0.5*F0_crit2.^2 - 1.5*F0_crit2.^(2/3) );
h_crit3 = H0_crit2 .* ( 1 + 0.5*F0_crit3.^2 - 1.5*F0_crit3.^(2/3) );
hhat_Houghton=h_crit3*N/U;


plot_case='constant_strat';
plot_case='inversion';
%plot_case='H vs h';

switch plot_case
    case 'constant_strat'

        figure
        % plot(hhat,F0_crit);
        % hold on

        %plotting here as for Fig. 15b
        plot(h_crit3*N/U,F0_crit3,'b'); %using Smith formula for F0 but with F0=U/(N*H0) - acheived by approximating
        %g' with N^2*H0 where H0=dz in N^2=g/theta * dtheta/dz
        hold on
        plot(hhat,F0_crit,'r'); %plot the original formulation for continuous stratification
        % are different but.... if multiply F0 by a factor of 0.62...
        plot(hhat_Houghton,F0_crit3*0.62,'g'); %they overlay each other almost exactly. Why 0.62???
        %is independant of N and U

        ylabel('F0=U/(N*H0)');
        xlabel('hN/U');


        set(gca,'xlim',[0 1.5]);
        set(gca,'ylim',[0 1.0]);

    case 'inversion'

        %now plot as for Fig. 15a (i.e. F0 vs. h/H0)
        figure

        %first give range of H0 values and calculate h using (16) of Smith and Sun
        plot(h_crit2./H0_crit2,F0_crit2,'b'); %using Smith formula for F0 for an inversion F0=U/sqrt(g'H0)
        %now give range of h' values and calculate H0 using (20). Note F0=1/H0' since H0' = H0*N/U
        hold on
        plot(2*hhat./H0_hat_strat,sqrt(2)*F0_crit,'r'); %plot the original formulation for continuous stratification with transformed variables
        % i.e. using Heff instead of H0 and F0=U/(N*Heff)=sqrt(2)U/(N*H0)
        
%        plot(1/(0.62^2)*hhat./H0_hat_strat,1/0.62*F0_crit,'k--');
        
        plot(hhat./H0_hat_strat/0.62,F0_crit/0.62,'g');
        
%        plot(hhat.*F0_crit,U./sqrt(gd*(1./F0_crit)),'r');
        %here are plotting hhat/H0_hat = h/H0 since H0_hat=1/F0 = N*H0/U
        % are different but.... if divide F0_crit by a factor of 0.62...
%        plot(hhat.*F0_crit/0.62,F0_crit/0.62,'g'); %they overlay each other almost exactly. Why 0.62???
        %is independant of N and U
        %

 %       plot(0.62*h_crit3./H0_crit2,F0_crit3*0.62,'k--');
        %so in order to get the critical line for the inversion case (like Fig. 15a) to match the critical line predicted using the
        %continuous stratification theory (as for 15b) the axes of 15a need to be multiplied by 0.62. So this should also be the case
        %for the critical line in Houghton (as the formulae are the same) - but then this should mean that the lines used for
        %predicting whether the jump is stationary or not and the location of the jumps could also be used for the
        %contiuous stratification case.
        
        
%         Hb_hat=H0_crit2*N/U;
%         heff_hat=Hb_hat/2 .*( 1 + 1./Hb_hat.^2 - 3/2 .*(sqrt(2)./Hb_hat).^(2/3) );
%         plot(2*heff_hat./Hb_hat,sqrt(2)*1./Hb_hat,'k--'); %is 2*heff_hat./Hb_hat because need h/Heff where Heff=Hb/2

        set(gca,'xlim',[0 1.5]);
        set(gca,'ylim',[0 2.0]);
        
        xlabel('h/H0');
        ylabel('F0');
        
        
        %%% so to transform constant stratification cases to plot in F0,h/H0 space (as for Houghton) need to halve the critical
        %%% height (height of the top of the stratification, Hb, i.e. height of dividing streamline) for the h/H0 part
        %%% and then use sqrt(2)*U/(N*Hb) as the Froude number on the y-axis. 

    case 'H vs h'

        %now plot as for Fig. 5 in Smith and Sun (i.e. Hhat vs. hhat)
        figure

        %first give range of H0 values and calculate h using (16) of Smith and Sun
        plot(h_crit2*N/U,H0_crit2*N/U,'b'); %using Smith formula for F0 for an inversion F0=U/sqrt(g'H)
        %now give range of h' values and calculate H0 using (20). Note F0=1/H0' since H0' = H0*N/U
        hold on
        %%plot(hhat.*F0_crit,F0_crit,'r'); %plot the original formulation for continuous stratification but assuming F0=U/NH
        plot(hhat,1./F0_crit,'r');
        %here are plotting hhat/H0_hat = h/H0 since H0_hat=1/F0 = N*H0/U
        % are different but.... if divide F0_crit by a factor of 0.62...
        
        
%        plot(hhat,1./(F0_crit/0.62),'g'); %not as good an overlay - but better than F0*sqrt(2)
        %BUT, H0hat-1.35 is better...
        
%%        plot(h_crit2*N/U,1/0.62*H0_crit2*N/U,'g'); 
        
        Hb=H0_crit2*N/U;
        heff_hat=Hb/2 .*( 1 + 1./Hb.^2 - 3/2 .*(sqrt(2)./Hb).^(2/3) );
%%        plot(heff_hat,Hb,'k--');
        
        Hc=H0_crit2-sqrt(2);
        Fc=U/N./Hc;
        hc = Hc .* ( 1 + 0.5*Fc.^2 - 1.5*Fc.^(2/3) );
%%        plot(hc*N/U,Hc*N/U,'b--');
        
        

%        plot(0.62*h_crit3./H0_crit2,F0_crit3*0.62,'k--');
        %so in order to get the critical line for the inversion case (like Fig. 15a) to match the critical line predicted using the
        %continuous stratification theory (as for 15b) the axes of 15a need to be multiplied by 0.62. So this should also be the case
        %for the critical line in Houghton (as the formulae are the same) - but then this should mean that the lines used for
        %predicting whether the jump is stationary or not and the location of the jumps could also be used for the
        %contiuous stratification case.

        set(gca,'xlim',[0 2.1]);
        set(gca,'ylim',[0 15.0]);


end




%by interpolating h for a given H0 are making the points fit the critical curve of Smith 87 - so is a bit pointless!
%However, the fact that the h values do not match the critical H0 values in Fig. 15b) of Smith 87 suggests that
%either that approximating those cases as continuos stratification ones is not correct (some are also plotted as a single layer
%with inversion case as in 15a ), or that the flow can choose a mountain height that it sees based on the forced H0
%for example due to a critical layer or wave breaking. One way to test whether the flow can choose a H0
%might be to examine the displacement values predicted
%by the Smith theory to see which H0 and h combination they match for.


idat=0;

%1) M6 case
idat=idat+1;
label_data(idat).dat='M6';
H_data(idat)=3500;
hhat_data(idat)=0.65;
h_data(idat)=800;
Ldat(idat)=hhat_data(idat)/h_data(idat);

%2) M7 case
idat=idat+1;
label_data(idat).dat='M7';
H_data(idat)=2200;
hhat_data(idat)=1.2;
h_data(idat)=800;
Ldat(idat)=hhat_data(idat)/h_data(idat);

%3) Boulder case
idat=idat+1;
label_data(idat).dat='Boulder';
H_data(idat)=9500;
hhat_data(idat)=0.45;
h_data(idat)=1600;
Ldat(idat)=hhat_data(idat)/h_data(idat);

%4) Windy Gap case
idat=idat+1;
label_data(idat).dat='Windy Gap';
H_data(idat)=700;
hhat_data(idat)=0.17;
h_data(idat)=100;
Ldat(idat)=hhat_data(idat)/h_data(idat);


for i=1:idat
    L=hhat_data(i)/h_data(i);
    hmax=interp1(HH,hh,H_data(i)*L)/L; %get hmax value in metres (max mountain height for given H0)
    H0=H_data(i)-hmax;
    H0 = H_data(i);
    H0=H_data(i)-1/L;  hmax=interp1(HH,hh,H0*L)/L; %get hmax value in metres (max mountain height for given H0)
    H0_max(i)=interp1(hh,HH,min([hhat_data(i) hh(end)]))/L;
    del_max(i)=interp1(hh,delL,min([hhat_data(i) hh(end)]))/L;

    F0(i)=1/L/H0;
    hx(i)=hmax*L;
    %     plot(hx(i),F0(i),'rx');
    %     text(hx(i),F0(i),label_data(i).dat);
    %     plot(hhat_data(i),1/L/H_data(i),'ro'); %original Smith formulation
    %     text(hhat_data(i),1/L/H_data(i),label_data(i).dat);
end





