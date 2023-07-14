%plot lon-height slice from make_z_lon_slice

scrsz=get(0,'ScreenSize');
posit=[scrsz];

LAT_val2={[-32.74 -28],[-22.74 -18;],[-12.74 -8;]};
LON_val2={[-104 -71.25],[-104 -71.25],[-104 -71.25]};


times_required2 = {[0:3],[12:15]};
time_mean_str ='ALL';

for itime=1:length(times_required2)
    for ilat_slice = 1:length(LAT_val2)
        
        
        times_required = times_required2{itime};
        LAT_val = LAT_val2{ilat_slice};
        LON_val = LON_val2{ilat_slice};  
        
        %get the data for the slice
        make_z_lon_slice

%% now plot        
        
        figname = titlenam;

        hf=figure('name',figname,'Position',posit);
        dpcolor(x_zslice,z_zslice,dat_zslice);
        set(gca,'ydir','reverse');
        colorbar


        titlenam = remove_character(titlenam,'_','-');
        titlenamw=textwrap({titlenam},70);
        htit=title(titlenamw,'fontsize',10);


        set(gca,'ylim',[500 1020]);
        set(gca,'clim',[0 0.15]);

        xlabel('Longitude');
        ylabel('Pressure (hPa)');

        fontsize_figure(gcf,gca,15);

        figname = remove_character(figname,'.','pt');
        figname = remove_character(figname,'for comparison','');
        figname = remove_character(figname,'  UTC','');  
        figname = remove_character(figname,'   ','');
        figname = remove_character(figname,'  ','');        
        savename = [savedir figname];
        saveas_ps_fig_emf(gcf,savename,'',0,0);


    end
end

