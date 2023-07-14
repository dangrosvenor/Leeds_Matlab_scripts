function UM_save_plot(gcf,isave_plot,savedir,save_str)

if isave_plot==1
    savename=[savedir save_str];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        opts.isavefig=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end
