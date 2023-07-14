isave_pdf = 0;

iplot_amip = 1;
iplot_obs = 1;

xval = 0.01;
xval = 0.03;
xval = 0.05;
xval = 0.07;
xval = 0.09;
xval = 0.11;
xval = 0.13;
xval = 0.15;
xval = 0.17;
xval = 0.19;

[minval,ibin] = min(abs(mid_Xbins_ens_mean{1} - xval));

%% Plot final PDF
figure('color','w');
set(gcf,'position',[3         197        856         422]);
ileg=1; clear leg_str
clear hline

ind_ens=0;
for iens=1:Nens
    if ind_ens==1
        figure('color','w');
        set(gcf,'position',[3         197        856         422]);
        ileg=1; clear leg_str
        clear hline
    end
    
    hline{iens}=plot(mid_Xbins_ens{iens},pdf2D_ens{iens}(ibin,1:end-1));
    
    if iens>1 & ind_ens==0
        %This line stops the line being listed in the legend :-
        set(get(get(hline{iens},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    elseif ind_ens==1
        leg_str{ileg}=['iens=' num2str(iens)]; ileg=ileg+1;
    else
        leg_str{ileg}='Ensemble members'; ileg=ileg+1;
    end
    
    set(hline{iens},'color',[0.5 0.5 0.5]);
    hold on
    
    if ind_ens==1 | iens==Nens
        
        hline2=plot(mid_Xbins_ens_mean{1},pdf2D_ens_mean{1}(ibin,1:end-1),'linewidth',3);        
        set(hline2,'color',[0 0 1]);
        leg_str{ileg}='Ensemble mean'; ileg=ileg+1;
        
        if iplot_amip==1
            hline3=plot(mid_Xbins_AMIP{1},pdf2D_AMIP{1}(ibin,1:end-1),'linewidth',3);
            set(hline3,'color',[1 0 0]);
            leg_str{ileg}='AMIP model'; ileg=ileg+1;
        end
        
        if iplot_obs==1
            hline4=plot(mid_Xbins_obs{1},pdf2D_obs{1}(ibin,1:end-1),'linewidth',3);
            set(hline4,'color',[0 0 0]);
            leg_str{ileg}='Observations'; ileg=ileg+1;
        end
        
        set(gca,'fontsize',18);
        
        grid on
        xlabel(ylabelstr);
        ylabel('No. datapoints');
        tit_str = [box_region_str ' ' land_ocean ',Xbin=' num2str(mid_Xbins_ens_mean{1}(ibin))];
        title(tit_str);
        legend(leg_str);
        
    end
    
end



if isave_pdf==1
    savename=[savedir_date 'PDF ' var_ukesm ' ' tit_str];
    %savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    
    savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
end




