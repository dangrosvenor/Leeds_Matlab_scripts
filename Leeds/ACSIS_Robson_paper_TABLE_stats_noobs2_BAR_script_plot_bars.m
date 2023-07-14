% Run from ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script.m

switch stacked_bar_or_error_bars
    case 'error_bar'
        %SW_paper_dSW_plot2
        %SW_paper_dSW_plot_horiz
        SW_paper_dSW_plot_horiz_just_vals
        
        
        
        
        
        
    case 'stacked_bar'
        
        ibar=ibar+1;
        
        ibar_start = ibar; %save this
        
        xlabs_str{ibar}='';
        
        iblock=1;
        bar_dat(ibar,iblock) = local_indirect; iblock=iblock+1;
        bar_dat(ibar,iblock) = local_direct; iblock=iblock+1;
        bar_dat(ibar,iblock) = non_local; iblock=iblock+1;
        bar_dat(ibar,iblock) = ghg; iblock=iblock+1;%GHGs
        bar_dat(ibar,iblock) = nat; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1; %to represent the actual model trend
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for net of first bar
        ibar=ibar+1;
        xlabs_str{ibar} = ''; %period_lab;
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = net; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for calcaulted trend
        ibar=ibar+1;
        xlabs_str{ibar} = period_lab;
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = had; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        if ilocal(ivar)==1
            %UKESM bar
            ibar=ibar+1;
            xlabs_str{ibar}=''; %'Period 1';
            iblock=1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = ukesm_local_indirect; iblock=iblock+1;
            bar_dat(ibar,iblock) = ukesm_local_direct; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        end
        
        %UKESM bar
        ibar=ibar+1;
        xlabs_str{ibar}=''; %'Period 1';
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = ukesm; iblock=iblock+1;
        
        
        
        
        
        ibar=ibar+1;
        xlabs_str{ibar}=''; %'Period 1';
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        %% Now repeat for uncertrainties
        ibar=ibar_start;
        iblock=1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = aer_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = ghg_un;  iblock=iblock+1;%GHGs
        bar_dat_UN(ibar,iblock) = nat_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;%to represent the actual model trend
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for net of first bar
        ibar=ibar+1;
        iblock=1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = net_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for calcaulted trend
        ibar=ibar+1;
        iblock=1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = had_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        if ilocal(ivar)==1
            iblock=1;
            ibar=ibar+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            
        end
        
        iblock=1;
        ibar=ibar+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = ukesm_un; iblock=iblock+1;
        
        
        iblock=1;
        ibar=ibar+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        
        
        % Define the colours oursleves as otherwise can get a weird coloru
        % order.
        %bar_cols = {'b','r','g','c','k','m',};
        bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]};
        offset_errbars=[0 0 1 1 1 0 0 0 0 0];
        
        
        if ilocal(ivar)==1
            if strcmp(var_str_tab,'SW')
                leg_str = {'Aerosol ACI forcing','Aerosol ARI forcing','Aerosol feedback','GHG','Natural','Net','HADGEM total','UKESM aerosol ACI forcing','UKESM aerosol ARI forcing','UKESM total'};
                %leg_str = {'Local aerosol','N/A','GHG','Natural','Net','HADGEM'};
                cut_inds=[];
            else
                leg_str = {'Aerosol forcing','Aerosol feedback','GHG','Natural','Net','HADGEM total','UKESM aerosol forcing','UKESM total'};
                cut_inds = [1 8]; %cut the ACI bar (rebadge as local aerossol)
            end
        else
            %leg_str = {'N/A','All aerosol','GHG','Natural','Net','HADGEM','N/A','UKESM'};
            %leg_str = {'All aerosol','GHG','Natural','Net','HADGEM'};
            %leg_str = {'All aerosol','N/A','GHG','Natural','Net','HADGEM'};
            
            %Cut out the local aerosol stacked bars and colours
            cut_inds = [1 2 8 9];
            leg_str = {'All aerosol','GHG','Natural','Net','HADGEM','UKESM'};
            
        end
        
        bar_dat(:,cut_inds)=[];
        bar_dat_UN(:,cut_inds)=[];
        bar_cols(cut_inds)=[];
        offset_errbars(cut_inds)=[];
        
        
        pad_list = [1:size(bar_dat,2)];
        
        bar_dat_UN_neg = bar_dat_UN;
        bar_dat_UN(bar_dat<0)=0;
        bar_dat_UN_neg(bar_dat>0)=0;
        
        %Matlab plots the stacked bars over the last one so that it covers them up
        %if the sign of the bars change. Here willl plot +ve and -ve separately and
        %then have a net bar.
        bar_dat_neg = bar_dat;
        bar_dat(bar_dat<0)=0;
        bar_dat_neg(bar_dat_neg>0)=0;
        
        
        
        
        
        
        figure
        set(gcf,'color','w');
        hbar = bar(bar_dat,'stacked'); %will give N bars each with M components for bar_dat[N,M]
        %manually set the colours to have control - otherwise not sure how to
        %get them as FaceColor is set to 'Flat'!
        for i=1:length(hbar)
            set(hbar(i),'FaceColor',bar_cols{i});
        end
        hold on
        hbar = bar(bar_dat_neg,'stacked');
        for i=1:length(hbar)
            set(hbar(i),'FaceColor',bar_cols{i});
        end
        %alpha(0.5); %Make transparent?
        increase_font_size_map_figures
        %title(['SW_{in}=' num2str(SW_in) ', A_{clear}=' num2str(A_clear)]);
        %if ioutput_abs==0
        %    ylabel('SW_{up TOA} Trend (W m^{-2} yr^{-1}) x 10^{-2}');
        %else
        %    ylabel('\Delta fc (W m^{-2})');
        %end
        if exponent~=0
            exp_str = [' x10^{' num2str(exponent) '}'];
        else
            exp_str='';
        end
        ylabel([ylabs{ivar} exp_str]);
        
        set(gca,'xticklabel',xlabs_str);
        
        % pad_list=[1 2];
        % for i=1:length(pad_list)
        %     iL=pad_list(i);
        %     leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
        % end
        
        for i=1:length(pad_list)
            iL=pad_list(i);
            leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
        end
        %L=legend(leg_str,'location','SouthWest','fontsize',14);
        L=legend(leg_str,'location','BestOutside','fontsize',14);
        %set(L,'string',leg_str)
        
        
        set(gca,'xlim',[0 ibar+1]);
        
        
        %cumsum_bar = cumsum(bar_dat,2);
        %cumsum_bar_neg = cumsum(bar_dat_neg,2);
        
        % Plot the error bars
        
        bar_dat_plot = bar_dat;
        bar_dat_UN_plot = bar_dat_UN;
        ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot_PLOT
        
        bar_dat_plot = bar_dat_neg;
        bar_dat_UN_plot = bar_dat_UN_neg;
        ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot_PLOT                                
        
        
end
