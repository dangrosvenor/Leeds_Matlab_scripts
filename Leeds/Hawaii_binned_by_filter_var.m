function [ymean,ystd,yn,ymean_overall,yn_overall,ystd_overall]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts)
    
    if isfield(opts,'ylims')==1
        i_ylim=1;
    else
        i_ylim=0;
    end        

    if min(size(var_filter) == size(var_test))==0
       error(['*** DPG - size of var_test (' num2str(size(var_test)) ') needs to be the same as the size of var_filter (' num2str(size(var_filter)) ')***']);
    end        
        
    for i=1:length(var_filter_bin_edges)-1
        %yvals(:,i)=Hawaii_plume_timeseries(dLWP,dNd,[var_filter_bin_edges(i) 1e99],time,title_str,ylab_str,ylims,opts);        
        imean=find(var_filter>=var_filter_bin_edges(i) & var_filter<var_filter_bin_edges(i+1));
        [ymean(i),yn(i),ystd(i)] = meanNoNan(var_test(imean),1);        
    end
    
    [ymean_overall,yn_overall,ystd_overall] = meanNoNan(var_test(:),1);
    
    %yvals_timemean = meanNoNan(yvals,1);
    figure
    %plot(var_filter_bin_edges,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(mid_points,ymean,'linewidth',3);
    xlabel(xlab_str);
    ylabel(['Time mean ' ylab_str]);   
    fontsize_figure(gcf,gca,18);
    grid on
    if i_ylim==1
       set(gca,'ylim',opts.ylims); 
    end
    
    figure
    %plot(var_filter_bin_edges,yvals_timemean,'linewidth',3);
    plot(mid_points,yn,'linewidth',3);
    %plot(var_filter_bin_edges,ystd./sqrt(yn),'linewidth',3);
    xlabel(xlab_str);
    ylabel('N_{samples}');
    fontsize_figure(gcf,gca,18);
    grid on
    set(gca,'yscale','log')  
    
    
    figure    
    plot(mid_points,ystd./sqrt(yn),'linewidth',3);
    %plot(var_filter_bin_edges,ystd./sqrt(yn),'linewidth',3);
    xlabel(xlab_str);
    ylabel(['\sigma/\surdN ' y_units_str]);
    fontsize_figure(gcf,gca,18);
    grid on
    %set(gca,'yscale','log')     
    
%


    