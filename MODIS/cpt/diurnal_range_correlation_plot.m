%diurnal range correlation plot
%use data stored in
%/home/disk/eos1/d.grosvenor/savefile_dirunal_ranges_all_models.mat
% Name                     Size            Bytes  Class     Attributes
% 
%   lwp_diurnal_save         1x1              4096  struct              
%   precip_diurnal_save      1x1              4118  struct              

plot_case = 'diurnal range';
%plot_case = 'mean_vals';
            

load('/home/disk/eos1/d.grosvenor/savefile_dirunal_ranges_all_models.mat');

figure

syms = {'*','o','^','d','p','o','o'};
cols={'r','y','b','m','k','w','c'};

nmodels=length(lwp_diurnal_save.labs);

for im=1:nmodels
    ilon=find(lwp_diurnal_save.xdat(im).x<80); %west of 80W
    switch plot_case
        case 'diurnal range'
            y = meanNoNan(2*lwp_diurnal_save.halfrange(im).dat(ilon),2);
            x = meanNoNan(24*2*precip_diurnal_save.halfrange(im).dat(ilon),2);
            xlab='Precip diurnal range (mm day^{-1})';
            ylab='LWP diurnal range (g m^{-2})';
            set(gca,'xlim',[0 1]);
            
        case 'mean_vals'
            ilon=3;
            y = meanNoNan(lwp_diurnal_save.ydat(im).y(ilon),2);
            x = meanNoNan(24*precip_diurnal_save.ydat(im).y(ilon),2);
            xlab='Mean Precip (mm day^{-1})';
            ylab='Mean LWP (g m^{-2})';
            set(gca,'xlim',[0 2]);
            title(['LON=' num2str(lwp_diurnal_save.xdat(im).x(ilon))]);
    end
    
    leg{im}=lwp_diurnal_save.labs(im).l;
    leg{1} = 'OBS';
    h=plot(x,y,[cols{im} syms{im}]);
    set(h,'markerfacecolor',cols{im},'markeredgecolor','k','markersize',12)
    text(x+24*0.0005,y-2,leg{im});
    hold on
end

xlabel(xlab);
ylabel(ylab);
fontsize_figure(gcf,gca,16);


savename=[savedir 'Diurnal_correlation_plot_CPT'];






