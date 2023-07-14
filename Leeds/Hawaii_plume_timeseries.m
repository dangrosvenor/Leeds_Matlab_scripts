function [dom_mean_var_plume,std_dev,N,inan]=Hawaii_plume_timeseries(var,var_filter,filter_vals,time,title_str,ylab_str,ylims,opts)
    
    if min(size(var_filter) == size(var))==0
       error(['*** DPG - size of var (' num2str(size(var)) ') needs to be the same as the size of var_filter (' num2str(size(var_filter)) ')***']);
    end

    inan = find(var_filter<filter_vals(1) | var_filter>filter_vals(2));
           
    var(inan)=NaN;
    %if var_filter contains NaN values then they won't satisfy the find for
    %inan above - so need to exclude those too.
    var(isnan(var_filter))=NaN;
    %N.B. - need to do this rather than mean over one dimension and then
    %the other since the latter will then weight all the means from one axis
    %equally when doing the the mean along the 2nd dir.
    var=permute(var,[3 1 2]);
    var=var(:,:);
    [dom_mean_var_plume,N,std_dev] = meanNoNan(var,2);           
    
    if isfield(opts,'no_plot')==0 | opts.no_plot~=1
        figure
        plot(time,dom_mean_var_plume,'linewidth',3);
        datetick('x','dd');
        
        %legend(leg_strs,'location','northwest');
        xlabel('Time (day in December)');
        %ylabel('SW surface forcing (W m^{-2})');
        ylabel(ylab_str);
        set(gca,'ylim',ylims);
        fontsize_figure(gcf,gca,18);
        grid on
        title(title_str);
        
    end
    