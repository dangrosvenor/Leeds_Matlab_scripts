savedir_DRIVER = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/';

method = 'constant dtau';
method = 'Mean curves from Odran'; %Run the "Nd and error calcs" block of MODIS_vert_pen_VOCALS_tau_multi_region_PDFs.m first
%with dtau_method = 'function of tauc';
method = 'Mean curves from Odran reff ratio';

switch method
    case 'constant dtau'
        dtau21=3.3; %the penetration depth from cloud top that corresponds to the retrieved re, allowing for penetration depth.
        %Fig. 4 of Platnick (2000) sugggests that this is around 335 for 2.1um.
        %=2.0 for 3.7um
        dtau37=2.0;

        tau=[3:0.15:30];  %Fig. 6 of Painemal indicates that tau during VOCALS ranged from 2 to 20 with plenty of values aroud the 5
        % region (or less)
        
        leg01 = ['2.1 \mum; \Delta\tau = ' num2str(dtau21)];
        leg02 = ['3.7 \mum; \Delta\tau = ' num2str(dtau37)];
        
        
    case 'Mean curves from Odran'
        tau=[3:0.15:35];
        dtau21 = interp1(tauc,tau_star21,tau);
        dtau37 = interp1(tauc,tau_star37,tau);  
        
        leg01 = ['2.1 \mum'];
        leg02 = ['3.7 \mum'];
        
    case 'Mean curves from Odran reff ratio'
        %I guess we could just use tau andc and reff_ratio21 directly...
        tau=[3:0.15:35];
        reff_ratio21_new= interp1(tauc,reff_ratio21,tau);
        reff_ratio37_new = interp1(tauc,reff_ratio37,tau);  
        
        leg01 = ['2.1 \mum'];
        leg02 = ['3.7 \mum'];

end

switch method
    case {'constant dtau', 'Mean curves from Odran'}
        %Estimate of the relative difference between the standard assumption that
        %re=re(H) and a revised version with re at some penetration dtau
        Nerr = (tau ./ (tau-dtau21)).^0.5;
        Nerr2 = (tau ./ (tau-dtau37)).^0.5;
        xvals = tau;
        
        Nerr = (tauc ./ (tauc - tau_star21)).^0.5;
        Nerr2 = (tauc ./ (tauc- tau_star37)).^0.5;       
        xvals = tauc;
        
        ylab = 'N_{standard} / N(\tau^{*})';      
        
        ylim_set=[1 3];
        
    case {'Mean curves from Odran reff ratio'}
        %Estimate of the relative difference between the standard assumption that
        %re=re(H) and a revised version with re at some penetration dtau
        Nerr = (reff_ratio21_new).^2.5;
        Nerr2 = (reff_ratio37_new).^2.5;
        Nerr = (reff_ratio21).^2.5;
        Nerr2 = (reff_ratio37).^2.5;
        ylab = 'N_{standard} / N(r_e(H))';
        xvals = tauc;
        
        ylim_set=[1 1.8];

end

lwidth=4;
fsize=16;


figure
plot(xvals,Nerr,'b-','linewidth',lwidth);
xlabel('\tau_c');
ylabel(ylab);

%title(['For \Delta\tau = ' num2str(dtau)]);

grid on
hold on


plot(xvals,Nerr2,'r--','linewidth',lwidth);


legend({leg01,leg02});
fontsize_figure(gcf,gca,fsize);

set(gca,'ylim',ylim_set);
set(gca,'xlim',[0 35])

savename = [savedir_DRIVER 'vert_pen_Nd_error_vs_tau'];

value_tau5_21 = interp1(xvals,Nerr,5)
value_tau5_37 = interp1(xvals,Nerr2,5)

value_tau10_21 = interp1(xvals,Nerr,10)
value_tau10_37 = interp1(xvals,Nerr2,10)


