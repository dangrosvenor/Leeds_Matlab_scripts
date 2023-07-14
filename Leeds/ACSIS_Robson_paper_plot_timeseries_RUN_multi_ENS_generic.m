%Plot lines for all ensemble members individually

multi_script = 'ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic';
%multi_script = 'ACSIS_Robson_paper_plot_timeseries_RUN_multi';
%multi_script = 'ACSIS_Robson_paper_plot_timeseries_RUN_multi_calipso_CCI';

for iens_plot=1:9
    eval(multi_script);
    %ACSIS_Robson_paper_plot_timeseries_RUN_multi_DeepC
end
