iloc=findheight(LAT,-67.67);

isave=1;

imethod='mean_at_lat';
switch imethod
        case 'mean_at_lat'
            figname = ['Mean contribution to melt energy at LAT=' num2str(LAT(iloc)) ' for ' filestr];
        case 'single_points'
            figname = ['Contribution to melt energy at LAT=' num2str(LAT(iloc)) 'LON=' num2str(LON(iloc)) ' for ' filestr];
    end
    titname = figname;
    savename = figname;
    hf=figure('name',figname,'Position',posit);

    ydat=[melt_dat(iloc) sw_dat(iloc) lw_dat(iloc) sh_dat(iloc)...
        lh_dat(iloc) sh_dat(iloc)+lh_dat(iloc) grd_dat(iloc)];

    ydat=[ydat; zeros([1 length(ydat)])]; %so that it draws the bars in different colours

    h=bar(ydat);

    legend('Melt','SW','LW','SH','LH','SH+LH','GRD');

    set(gca,'fontsize',20);

    set(gca,'xlim',[0.6 1.6]);
%    set(gca,'ylim',[minALL(ydat) maxALL(ydat)]);
    set(gca,'ylim',[-40 85]);
    %set(gca,'xscale','log');
    %set(gca,'yscale','log');


    set(gca,'xticklabel',''); %remove x-axis labels

    %xlabel('Critical Supersaturation (%)');
    ylabel('Melt rate (mm day^{-1})');
    ylabel('Total melting (mm w.e.)');
    
    title(titname,'fontsize',16);

    if isave==1
        multisaveplot
    end
    

