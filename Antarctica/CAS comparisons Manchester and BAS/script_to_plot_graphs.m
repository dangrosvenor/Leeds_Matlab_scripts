%plotting the results for the BAS and Manchester CAS comparisons

isave_script_to_plot_graphs=1;

%case 1
time_plot=[10.4712 10.4842; 1+10.4400 1+10.4635]; %BAS CAS was running on summertime
%%time_plot=[10.4219 10.4426; 10.4400+1 10.4635+1]; %BAS CAS was running on summertime
%case 2
time_plot=[11.1316 11.1443; 1+11.1394 1+11.1800];
%case 3
time_plot=[12.8076 12.8484; 1+12.7673 1+12.7974];
%15 micron calibration
time_plot=[13+1/60+43/3600 13+01/60+48.5/3600; 1+12+54/60+13.5/3600 1+12+54/60+19.5/3600];
%20 micron calibration
time_plot=[13+3/60+3/3600 13+03/60+9.5/3600; 1+12+55/60+35.5/3600 1+12+55/60+42.5/3600];
%below is the second burst from the MAN CAS - v. similar distributions
%time_plot=[13+3/60+22/3600 13+03/60+25.5/3600; 1+12+55/60+35.5/3600 1+12+55/60+42.5/3600];
%30 micron calibration
time_plot=[13+5/60+7.5/3600 13+05/60+14.5/3600; 1+12+56/60+25.5/3600 1+12+56/60+33.5/3600];
%40 micron calibration
time_plot=[13+6/60+23/3600 13+06/60+33/3600; 1+12+57/60+15.5/3600 1+12+57/60+25.5/3600];
%1st burst of BAS CAS calibration (above was the second longer burst)
%time_plot=[13+6/60+23/3600 13+06/60+33/3600; 1+12+57/60+4/3600 1+12+57/60+9/3600];

%%% 29th July runs (in the chamber)
%order here is MAN, BAS, FSSP
%run1
time_plot=[12+53/60+0/3600 12+58/60+0/3600; 12+53/60+0/3600 12+58/60+0/3600];
time_plot=[12+53/60+9/3600 12+53/60+45/3600; 12+53/60+9/3600 12+53/60+45/3600];
%run2
time_plot=[13+24/60+0/3600 13+42/60+0/3600; 13+24/60+0/3600 13+42/60+0/3600; 3+24/60+0/3600 13+42/60+0/3600];

%50 um bead calibration
time_plot=[15.9145 15.9153; 15.8775 15.8785; 15.8490 15.8496];
time_plot=[15.9145 15.9153; 15.8751 15.8756; 15.8490 15.8496];
time_plot=[15.9145 15.9153; 15.8751 15.8785; 15.8490 15.8496];
time_plot=[15.9206 15.9214; 15.8751 15.8785; 15.8490 15.8496];
%time_plot=[15.9145 15.9150; 15.8751 15.8785; 15.8490 15.8496];
time_plot=[15.9145 15.9214; 15.8751 15.8785; 15.8490 15.8496];

  % for 0729346.min FSSP file
time_plot=[15.9145 15.9214; 15.8751 15.8785; 15.7815 15.7821];


%60 um bead calibration
%time_plot=[16.026 16.0276; 15.989 15.9904; 15.9664 15.9678];

%expt 1:-
time_plot=[12.885 12.962; 12.885 12.962; 12.885 12.962; 12.885 12.962];
%expt 2:-
time_plot=[13+25/60 13+32/60; 13+25/60 13+32/60; 13+25/60 13+32/60; 13+25/60 13+32/60];
time_plot=[13+32/60 13+42/60; 13+32/60 13+42/60; 13+32/60 13+42/60; 13+32/60 13+42/60];
%expt 3:-
%time_plot=[13+54/60 14+08/60; 13+54/60 14+08/60; 13+54/60 14+08/60; 13+54/60 14+08/60];
%time_plot=[14+1.5/60 14+08/60; 14+1.5/60 14+08/60; 13+54/60 14+08/60; 13+54/60 14+08/60]; %towards the end where LWC agreement is better
%expt 4:-
%time_plot=[14+22/60+7/3600 14+36/60+40/3600; 14+22/60+7/3600 14+36/60+40/3600; 14+22/60+7/3600 14+36/60+40/3600; 14+22/60+7/3600 14+36/60+40/3600];
%time_plot=[14+29/60+0/3600 14+36/60+0/3600; 14+29/60+0/3600 14+36/60+0/3600; 14+29/60+0/3600 14+36/60+0/3600; 14+29/60+0/3600 14+36/60+0/3600];
%expt 5:-
%time_plot=[15+21/60+18/3600 15+32/60+58/3600; 15+21/60+18/3600 15+32/60+58/3600; 15+21/60+18/3600 15+32/60+58/3600; 15+21/60+18/3600 15+32/60+58/3600];
%the above deliberately misses the initial burst of LWC seen in the LWC time-height trace as just wanted to sample the period where the cloud was steady
%time_plot=[15+18/60+0/3600 15+32/60+58/3600; 15+18/60+0/3600 15+32/60+58/3600; 15+18/60+0/3600 15+32/60+58/3600; 15+18/60+0/3600 15+32/60+58/3600];
%this uses the whole time period of experiment 5 (above)

% time_plot=[21.3125 21.3297]; %flight 102 section where CAS saw LWC, but not hotwire
% time_plot=[20.2795 20.2934]; %flight 102 section where CAS and hotwire saw LWC, 50-54 km
% 
% time_plot=[20.7894 20.7902; 20.7894 20.7902]; %flight 102 - where hotwire saw more than CAS - around 3.32 km.
% time_plot=[20.7790 20.7865; 20.7790 20.7865];

%time_plot=[0.8618 0.8623; 0.8684 0.8686]*24; %flight 102 - %comparing the size distributions 'good' CAS LWC points to the 'bad' ones as determined by Estimate_LWC_from_equiv_potemp
%with 20% tolerance - these were 'good' then 'bad'
%time_plot=[20+44/60+32/3600 20+44/60+35.7/3600; 20+42/60+11.5/3600 20+42/60+27.7/3600]; %flight 102 - good then bad
%time_plot=[20+38/60+24/3600 20+38/60+40.8/3600]; %flight 102    
%time_plot=[20+49/60+25/3600 20+49/60+36/3600; 20+50/60+33/3600 20+50/60+45/3600]; %flight 102 - good then bad


disp('*** Remember to set the correct time_plot values ***');

clear plot_graph_case
istring=1;
%plot_graph_case{istring}='Ratio of BAS to MAN counts as function of size'; istring=istring+1;
plot_graph_case{istring}='Size distributions MAN and BAS CAS on same plot'; istring=istring+1;
  %note, this can be used for single instrument size dists too
%plot_graph_case{istring}='Total number timseries'; istring=istring+1;

fsize=12;


for iplot_comparison=1:length(plot_graph_case)

    plot_graph_case2 = plot_graph_case{iplot_comparison};
    
switch plot_graph_case2
    case 'Need to set these'

        %N per bin plots 
        itype='timh';
        var_plot='N (cm^{-3})';
        i577 = 'Particle size dist vs time';
        imaxovr=1;
        iminovr=1;
        logflag=0;
        mincovOvr = 0;
        maxcovOvr = 100;

        %mincovOvr = 0;
        %maxcovOvr = 20;

        %logflag=1;
        %mincovOvr=1;



        %N per bin plots (BAS CAS)
        instrument='CAS no CIP';
        iset_colour_limits=1;
        iset_plot_and_instrument=1;
        man_choose_plotTimeHeight_graph=1;
        ichoose_itype_multisaveplot=1;
        if isave_script_to_plot_graphs==1
            multisaveplot;
        else
            plotTimeHeightVap3;
        end


        %LWC per bin plots 
        itype='timh';
        var_plot='LWC size dist (g cm^{-3})';
        i577 = 'Particle size dist vs time';
        imaxovr=1;
        iminovr=1;
        logflag=0;
        mincovOvr = 0;
        maxcovOvr = 0.05;

        %logflag=1;
        %mincovOvr=1;

        %(Manchester CAS)
        instrument='CAS MAN';
        iset_colour_limits=1;
        iset_plot_and_instrument=1;
        man_choose_plotTimeHeight_graph=1;
        ichoose_itype_multisaveplot=1;
        if isave_script_to_plot_graphs==1
            multisaveplot;
        else
            plotTimeHeightVap3;
        end

        %N per bin plots (BAS CAS)
        instrument='CAS no CIP';
        iset_colour_limits=1;
        iset_plot_and_instrument=1;
        man_choose_plotTimeHeight_graph=1;
        ichoose_itype_multisaveplot=1;
        if isave_script_to_plot_graphs==1
            multisaveplot;
        else
            plotTimeHeightVap3;
        end


            %'MAN CAS time-height size distribution'
            %(Manchester CAS)
            instrument='CAS MAN';
            iset_colour_limits=1;
            iset_plot_and_instrument=1;
            man_choose_plotTimeHeight_graph=1;
            ichoose_itype_multisaveplot=1;
            if isave_script_to_plot_graphs==1
                multisaveplot;
            else
                plotTimeHeightVap3;
            end



        case 'Size distributions MAN and BAS CAS on same plot'
        %size distributions
        itype='prof';
        graph=87;
        air_speed_type='constant';
        airspeed_constant=10;
        instrument_sd_all={'MAN CAS','BAS CAS','FSSP','Welas'};
        instrument_sd_all={'MAN CAS','BAS CAS','Welas'};
        instrument_sd_all={'BAS CAS','MAN CAS'};
%        instrument_sd_all={'BAS CAS','BAS CAS'};
%        instrument_sd_all={'BAS CAS'};        
        ylab='N (cm^{-3})';
        ylab='dN/dlogD (cm^{-3} \mum^{-1})';
%        ylab='dLWC/dlogD (g m^{-3} \mum^{-1}))'
%        ylab='LWC (g m^{-3})'; %prob best to just use LWC
        instrument_ratio=0;
        iytick_relabel=1;       
        y_axis_type='log10_matlab';
 %       y_axis_type='';        
        x_axis_type='log10_matlab';
%        x_axis_type='';


        man_choose_water_graph=1;
        ichoose_itype_multisaveplot=1;
        ichoose_times_size_dists=1;
        if isave_script_to_plot_graphs==1
            multisaveplot;
        else
            waterVapourMay2005;
        end


        case 'Ratio of BAS to MAN counts as function of size'
            %ratio of size distributions
            itype='prof';
            graph=87;
            air_speed_type='constant';
            airspeed_constant=10;
            instrument_sd_all={'MAN CAS','BAS CAS'};
            ylab='N (cm^{-3})';    
            instrument_ratio=0;
                        
            iytick_relabel=1;       
            y_axis_type='log10_matlab';
            y_axis_type='';        
            x_axis_type='log10_matlab';
    %        x_axis_type='';


            man_choose_water_graph=1;
            ichoose_itype_multisaveplot=1;
            ichoose_times_size_dists=1;
            instrument_ratio=1;
            if isave_script_to_plot_graphs==1
                multisaveplot;
            else
                waterVapourMay2005;
            end

        case 'Total number timseries'

            itype='prof';
            graph=46; %timeseries
            man_choose_itimser=1;
            itimser = 'CAS plots';
            air_speed_type='constant';
            airspeed_constant=17;
            instrument={'MAN CAS','BAS CAS','FSSP','Welas'};
            instrument={'MAN CAS','BAS CAS','Welas'};   
            instrument={'BAS CAS'};   
            time_graph = {'Total number CAS','Total number CAS','Total number CAS','Total number CAS'};
            time_graph = {'Total number CAS','Total number CAS','Total number CAS'};  
            time_graph = {'Total number CAS'};              
        %    time_graph = {'MVD cutoff','MVD cutoff','MVD cutoff'};    
            xlims=0;
            xlimits=[time_plot(1,1) time_plot(1,2)]/24;

            ifull_screen=0;
            ioverride_xlims=1;
            man_choose_flt_graph=1; %timeseries
            man_choose_water_graph=1;
            ichoose_itype_multisaveplot=1;
            if isave_script_to_plot_graphs==1
                multisaveplot;
            else
                waterVapourMay2005;
            end




    end


end





