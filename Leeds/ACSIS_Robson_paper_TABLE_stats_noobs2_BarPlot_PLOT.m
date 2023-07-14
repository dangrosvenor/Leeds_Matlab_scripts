% Does the error bars for the bar plots - is run from ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot.m


cumsum_bar_plot = cumsum(bar_dat_plot,2);

% Do the actual bar plot
    
    for i=1:size(cumsum_bar_plot,1)
        %dat = [0 cumsum_bar_plot(i,:)]; %pad with a zero for mid points.
        %dat_mid = 0.5 * (dat(1:end-1) + dat(2:end));
        for j=1:size(cumsum_bar_plot,2) 
            if abs(bar_dat_plot(i,j)) > 1e-10
                %if j>=3 & j <= 5 %the stacked bars with error bars - offset the error bar x pos
                if offset_errbars(j)==1 %offset_errbars=[0 0 1 1 1 0 0 0 0 0]; or similar - indices where to offset
                    inds = find(offset_errbars==1);
                    k = min(inds)+1;                    
                    ipos = i + (j-k)*0.2; %so have i-dx, i and i+dx for j=2,3,4
                %elseif j>=7 & j <= 8 %the stacked bars - offset the error bar x pos
                %    ipos = i + (j-7)*0.2; %so have i-dx, i and i+dx for j=2,3,4
                else
                    ipos=i;
                end
                                
                if bar_cols{j}==[0 0 0]
                    bcol='g';
                else
                    bcol='k';
                end
                %if j>2 & j ~= 8 & j ~= 9 %don't want to plot for the local aerosol bars
                if bar_dat_UN_plot(i,j)>1e-10
                    herr=errorbarYY('vert',ipos,cumsum_bar_plot(i,j),bar_dat_UN_plot(i,j),gca,bcol,'o',4,0.01);
                    herr=errorbarYY('vert',ipos,cumsum_bar_plot(i,j),bar_dat_UN_plot(i,j),gca,bar_cols{j},'o',2,0.005);
                    plot(ipos,cumsum_bar_plot(i,j),[bcol 'o'],'markerfacecolor',bar_cols{j});
                end
            end
%            if bars_sig_neg(i,j)==1
               %plot(i,dat_mid(j),'ko','markerfacecolor','w');
            %end                            
        end
    end
    
    





