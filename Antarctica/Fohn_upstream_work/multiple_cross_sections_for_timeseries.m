%runs PlotTimeHeightVap3 several times in order to create multiple cross
%secitons and store the data for the timeseries plot
% This saves each cross section in a sperate .mat file.
% Need to next run make_cross_section_timeseries to extract a time series
% for a particular point. 
% Then case 144 of watervap. to plot



logflag=0;
dlogflag=0;
noplot=0; %flag to say to just do the calculations and not to plot
    
    
i577 = 'wrf_vert_cross_section';

%var_plot_multi = {'Potential temperature (K)','Wind speed (m s^{-1})','Component horizontal (UV) wind speed (m s^{-1})'...
%    ,'Wind direction (degrees)','Vertical wind speed (m s^{-1})'};
%var_plot_multi = {'Vertical wind speed (m s^{-1})'};
var_plot_multi={'Wind speed perpendicular component (m s^{-1})'};
var_plot_multi={'Relative humidity (%)'};

time_multi = [1:25];  %6th Jan:- 11=6UTC, 12=9UTC, 13=12UTC, 15=18UTC, 17=00UTC

for time2=1:length(time_multi)
    time=time_multi(time2);
    for ig_multi=1:length(var_plot_multi)
        var_plot = var_plot_multi{ig_multi};
        man_choose_plotTimeHeight_graph=1;
        plotTimeHeightVap3
        save_cross_section_data2
    end

end