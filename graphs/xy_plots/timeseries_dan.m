if ~exist('man_choose_itimser')
    %itimser='cumatHtotiN_resdiff';
    %itimser='cumatHtotiM_resdiff';
    %itimser='cumatHtotwAD_resdiff';
    %itimser='cumatHtotw_resdiff';
    %itimser='cumatHvap_resdiff'; %advective and microphysical source of vapour
    %itimser='cumatHchangeiceNC_resdiff';
    % itimser='gwave_k-spectra';
    % %itimser='vap_dist';
    % %itimser='av_w';
    % %itimser='cumatHdepsub_resdiff';
    % %itimser='fall_comp';
    % %itimser='low_tracer';
    % %itimser='ice_dep_rate';
    % %itimser='ice_mass';
    % %itimser='LNB_max_dqtot_path';
    % %itimser='av_rad';
    % %itimser='cumatHtoti'; %cumulative plots of advection and fall speed for total ice
    %itimser='cumatHtotw'; %cumulative plots of advection and fall speed for total water
    % itimser='cumatHvapflux'; %cumulative plots of advection of vapour
    % %itimser='min_icesatMR'; %cumulative plots of advection and fall speed for total water
    %itimser='ice_num'; %cumulative plots of advection and fall speed for total water
    %itimser='ice_dep_rate'; %
    % %itimser='ice_mic_rate'; %
    % %itimser='ice_mass'; %
    % itimser='emm_dcw';
    % itimser='emm_ncw';
    % itimser='emm_rwc'; %
    % %itimser='emm_nrwc'; %
    % %itimser='emm_lwc'; %
    % %itimser='emm_iwc'; %
    % %itimser='emm_isg'; %
    % %itimser= 'ice_proc_rates';
    % itimser= 'dqtotsum'; %sum of total water points below 5 ppmv
    %itimser= 'dqvapsum'; %sum of total water points below 5 ppmv
    % itimser= 'rhopert_vap'; %sum of density perts for points with vapour below 5 ppmv
    % %itimser='cumatH_tracer';
    %itimser= 'eddy_flux'; %timeseries of eddy flux th'w'
    % %itimser = 'height_dqvapmax';
    %itimser = 'mode_diam'; %mode diamter of ice mass
    %itimser = 'echo_top'; %echotop timeseries
    %itimser = 'max_supersat';
    %itimser = 'max_MR_hslice';
    %itimser = 'max_MR_totprc';
    %itimser = 'mean_vapour';
    %itimser = 'latent_heating';
    %itimser = 'av_up'; %check out cloud top too from profiles
    %itimser = 'max_hm_gm3';
    %itimser = 'max_mac3';
    itimser = 'antjan06_flt';
%    itimser = 'AWS_antjan06';
    %itimser = 'wrf_timser';
    itimser = 'CAS plots';
%    itimser = 'ACPIM timeseries';
else
    clear man_choose_itimser
end

isum_smooth=0; %flag to tell the smoothing function to sum rather than average in each window.

multi={'ice_dep_rate','ice_mass','ice_num'};

if subplotting==1
    itimser=multi{iplot};
end


idatetick=0; %flag to say the want the xaxis in proper time format rather than decimal time
%specify the type with datetick_type (see help datetick)


try
    npes=npess2(idir); %no. of processors

    izmin=3;
    izmax=length(GridDan(idir).Z);

    itmin=999;
    for itim=1:length(icediagsALL)
        itmin=min([itmin size(icediagsALL(itim).i,2)]);
    end

    dumprange=[1:itmin];
    dumprange=[1:44];
    %	dumprange=[1:15];
    dumprange=[1:33];


    ixtime=1; %**** make sure this is set to zero if not doing timeseries **** flag to say x axis is time so that hours above >=24 are renumbered

    logflag=0;

    izlim=0;
    %z=GridDan(idir).Z;





    time1=18;
    time2=19;

    t1=findheight(GridDan(idir).t+3,time1);
    t2=findheight(GridDan(idir).t+3,time2);

    %  t1=dumprange(1);
    %     t1=27;
    %  t2=dumprange(end);
    %t2=dumprange(62);



    xlims=1;
    xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];

    dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;

catch %in case GridDan(idir).t doesn't exist, if say doing a SER timeseries
end

xlab=['Time (UTC)'];
%xlab=['Time (hrs)'];


ititle=1;

switch itimser
     case 'template'

        %timeseries of data from the ACPIM model
        clear plot_cases
        iplot_cases=0;

        DT_acpim=0.5; %timestep in seconds
        xlab=['Time (s)'];

        iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Temperature ^{o}C';
        %iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Saturation Ratio Liquid';
        %iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Saturation Ratio Ice';


                iytick_relabel=0;    

                nmark=0;
                lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                ismooth=[0 0 0];

                izlim=0;
                zmin=0;
                zmax=0.5;

                idat=0;

        for idat=1:length(plot_cases)
                xdat(idat).x=[1:length(acpim.T)]*DT_acpim;        

                switch plot_cases{idat}
                    case 'Saturation Ratio Liquid'
                        ydat(idat).y = acpim.R;
                        ylab=['Saturation Ratio'];      
                        labs(idat).l='Saturation Ratio Liquid';

                    case 'Saturation Ratio Ice'
                        ydat(idat).y = acpim.Rice;
                        ylab=['Saturation Ratio'];      
                        labs(idat).l='Saturation Ratio Ice';    

                    case 'Temperature ^{o}C'
                        ydat(idat).y = acpim.Rice;
                        ylab=['Saturation Ratio'];      
                        labs(idat).l='Saturation Ratio Ice'; 
                end


        end
        
        titlenam = ['ACPIM plots'];
        figname=[titlenam '-' ylab];
        savename=figname;


    case 'ACPIM timeseries'

%timeseries of data from the ACPIM model
clear plot_cases
iplot_cases=0;

DT_acpim=0.5; %timestep in seconds
xlab=['Time (s)'];

%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Temperature ^{o}C';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Saturation Ratio Liquid';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Saturation Ratio Ice';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Ice number per litre';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Ice number per kg';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Ice mass';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Ice Diameter';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='Number of droplets';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='LWC';
iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='ACPIM Ice number vs CIP number';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='ACPIM Ice mean ice size vs CIP';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='ACPIM Simple Ice mean ice size vs CIP';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='ACPIM LWC vs aircraft';
%iplot_cases=iplot_cases+1; plot_cases{iplot_cases}='ACPIM droplet no. vs CAS';

        iytick_relabel=0;    
        
        nmark=0;
        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
        ismooth=[0 0 0];
                
        izlim=0;
        zmin=0;
        zmax=0.5;
        
        idat=0;
        
        rho=density(acpim.P,acpim.T);
        
        x_start = 146;
        x_end=156.9;
        %calculate the distances flown by the plane
        eval(['X_flt = X_flt' flight_no ';']);
        eval(['Y_flt = Y_flt' flight_no ';']);
        dist = [0; cumsum(sqrt( (diff(X_flt)).^2 + (diff(Y_flt)).^2 ))];

        [inds_1 inds_2] = findheight_nearest(dist,x_start,x_end);
        inds2=inds_1:inds_2;
        
        X_segment=dist(inds2);
        Z_segment = dat_flt(inds2,col_alt); %altitude (m)
        
        time_flt=time_flt102;
        
        time_CIP = CIP_time_Jonny2/3600;
        dist_CIP = interp1(time_flt,dist,time_CIP);
        [i1,i2]=findheight_nearest(dist_CIP,min(X_segment),max(X_segment));
        inds_CIP=[i1:i2];
        
        time_acpim=[1:size(acpim.Z,1)]*DT_acpim;   
        dist_acpim = 20*time_acpim/1000 + x_start;
        %height of aircraft on the acpim horiz axis
        Zseg_acpim = interp1(X_segment,Z_segment,dist_acpim);        
        

        
for idat=1:length(plot_cases)
     
        xdat(idat).x=time_acpim;
        
        switch plot_cases{idat}
            case 'Saturation Ratio Liquid'
                ydat(idat).y = acpim.R;
                ylab=['Saturation Ratio'];      
                labs(idat).l='Saturation Ratio Liquid';
                
            case 'Saturation Ratio Ice'
                ydat(idat).y = acpim.Rice;
                ylab=['Saturation Ratio'];      
                labs(idat).l='Saturation Ratio Ice';    
                
            case 'Temperature ^{o}C'
                ydat(idat).y = acpim.T-273.15;
                ylab=['Temperature (^{o}C)'];      
                labs(idat).l='Temperature'; 
                
            case 'Ice number per litre'
                ydat(idat).y = rho.*acpim.conci*1e-3; %convert from #/kg to #/litre
                ylab=['Ice Number Concentration (L^{-1})'];      
                labs(idat).l='Ice No.'; 
             
            case 'Ice number per kg'
                ydat(idat).y = acpim.conci; %convert from #/kg to #/litre
                ylab=['Ice Number Concentration (kg^{-1})'];      
                labs(idat).l='Ice No.';     
                
            case 'Ice mass'
                ydat(idat).y = rho.*acpim.massi*1e6; %convert from kg/kg to mg/m3
                ylab=['Ice Mass (mg m^{-3})'];      
                labs(idat).l='Ice Mass';     
                
            case 'Ice Diameter'
                ydat(idat).y = acpim.diami; %
                ylab=['Ice Diameter (\mum)'];      
                labs(idat).l='Ice Diameter';  
                
            case 'Number of droplets'
                ydat(idat).y = 1e-6*rho.*acpim.conc2; %
                ylab=['Number of droplets (cm^{-3})'];      
                labs(idat).l='No. droplets';     
                
             case 'LWC'
                ydat(idat).y = 1e3*rho.*acpim.qc; %
                ylab=['LWC (g m^{-3})'];      
                labs(idat).l='LWC';  
                
            case 'ACPIM Ice number vs CIP number'
                ylab=['Ice number (L^{-1})'];
                xlab=['Distance along flight track (km)'];                                

                ydat(1).y = 1000*ice_no_Jonny(inds_CIP);
                xdat(1).x = dist_CIP(inds_CIP);                   
                labs(1).l='CIP';  

                xdat(2).x = dist_acpim;                 
                %interpolate each profile for the value at the aircraft height
                for ix=1:size(rho,1)
                    ydat(2).y(ix)=interp1(acpim.Z(ix,:),rho(ix,:).*acpim.conci(ix,:)*1e-3,Zseg_acpim(ix));
                end
                labs(2).l='ACPIM';  

                
                lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'ACPIM Ice mean ice size vs CIP'
                ylab=['Ice size (\mum)'];
                xlab=['Distance along flight track (km)'];                                

                ydat(1).y = mean_ice_size(inds_CIP);
                xdat(1).x = dist_CIP(inds_CIP);                   
                labs(1).l='CIP';  

                xdat(2).x = dist_acpim;                 
                %interpolate each profile for the value at the aircraft height
                for ix=1:size(rho,1)
                    ydat(2).y(ix)=interp1(acpim.Z(ix,:),acpim.diami(ix,:),Zseg_acpim(ix));
                end
                labs(2).l='ACPIM';  

                
                lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

            case 'ACPIM Simple Ice mean ice size vs CIP'
                ylab=['Ice size (\mum)'];
                xlab=['Distance along flight track (km)'];                                

                ydat(1).y = mean_ice_size(inds_CIP);
                xdat(1).x = dist_CIP(inds_CIP);                   
                labs(1).l='CIP';  

                xdat(2).x = dist_acpim;  
                rho_ice=0.3e3; %kg/m3
                Vi_mean = acpim.massi./acpim.conci/rho_ice; %mean vol m3
                Di_mean = 1e6*(6/pi*Vi_mean).^(1/3); %um
                
                %interpolate each profile for the value at the aircraft height
                for ix=1:size(rho,1)
                    ydat(2).y(ix)=interp1(acpim.Z(ix,:),Di_mean(ix,:),Zseg_acpim(ix));
                end
                labs(2).l='ACPIM';  

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'ACPIM LWC vs aircraft'
                
                ylab=['LWC (g m^{-3})'];
                xlab=['Distance along flight track (km)'];                                
                
                LWC = interp1(CIP_time_all/3600,LWC_CAS_all,time_flt); %LWC (g/m3) - interpolate for dat_flt time base
                ydat(1).y = LWC(inds2);
                xdat(1).x = X_segment;
                labs(1).l='Hotwire';                 

                xdat(2).x = dist_acpim;                 
                
                %interpolate each profile for the value at the aircraft height
                for ix=1:size(rho,1)
                    ydat(2).y(ix)=interp1(acpim.Z(ix,:),1e3*rho(ix,:).*acpim.qc(ix,:),Zseg_acpim(ix));
                end
                labs(2).l='ACPIM';  
                
                
                LWC_cas = interp1(CAS_time_all/3600,LWC_dist_cas,time_flt); %CAS LWC (g/m3) - interpolate for dat_flt time base
                ydat(3).y = LWC_cas(inds2);
                xdat(3).x = X_segment;
                labs(3).l='CAS'; 

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'ACPIM droplet no. vs CAS'
                
                ylab=['Droplet no. (cm^{-3})'];
                xlab=['Distance along flight track (km)'];                                
                
                CAS_no = interp1(CAS_time_all/3600,CAS_total_number{icas_count},time_flt); %for sizes abouve cut_off_size
                ydat(1).y = CAS_no(inds2);
                xdat(1).x = X_segment;
                labs(1).l='CAS';  

                xdat(2).x = dist_acpim;                 
                
                %interpolate each profile for the value at the aircraft height
                for ix=1:size(rho,1)
                    ydat(2).y(ix)=interp1(acpim.Z(ix,:),1e-6*rho(ix,:).*acpim.conc2(ix,:),Zseg_acpim(ix));
                end
                labs(2).l='ACPIM';  

                lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
        end

        
end

 titlenam = ['ACPIM plots'];
 figname=[titlenam '-' ylab];
 savename=figname;
        
        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'CAS plots'

%timeseries of data from the CAS instrument 

%eval(['dat_flt = dat_flt' flight_no ';']);  %put the data for the required flight here
%eval(['time_flt = time_flt' flight_no ';']);
        iytick_relabel=0;    
        

        nmark=0;
        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
        
        if exist('override_binned_options') &  override_binned_options==1
        else
            ismooth=[0 0 0];
        end

        
        
        izlim=0;
        zmin=0;
        zmax=0.5;
        
        iaxis_square=0;
        iplot_error_bar=0;
        
        if ~exist('man_choose_flt_graph')
            
            time_graph = {'Total number CAS','Total number CAS'};
%            time_graph = {'LWC CDP Karl','LWC CAS Karl'};
            time_graph = {'Total number CAS'};   
%            time_graph = {'Total number CAS Karl'};    
%            time_graph = {'Total number CDP Karl','Total number CAS Karl'};                
%            time_graph = {'Total number normalised by bin width','Total number normalised by bin width','Total number normalised by bin width'};            
            
%            time_graph = {'Total number CIP'};
%            time_graph = {'Total large ice number CIP'};
%            time_graph = {'LWC_hotwire'}; %LWC as measured by the hotwire probe
%            time_graph = {'LWC_dist_CAS'}; %LWC as calculated from the size distribution
%            time_graph = {'LWC_dist_CIP'};
%            time_graph = {'LWC_dist_CAS','LWC_hotwire','LWC_dist_CAS_cutoff'};
%            time_graph = {'LWC_dist_CAS','LWC_hotwire'};            
%            time_graph = {'LWC_dist_CAS_cutoff','LWC_hotwire'};            
%            time_graph = {'LWC_dist_CAS','LWC_hotwire'};            
%            time_graph = {'Total number CPC'};
%            time_graph = {'Last bin number CAS'};
%             time_graph = {'Second bin number CIP'};
%             time_graph = {'Number in latter bins for CAS'};
%             time_graph = {'Second bin dN/dlogD CIP'};
%             time_graph = {'Reject diags'};
%             time_graph = {'PADS spreadsheet','PADS spreadsheet Dan''s values'};
%              time_graph = {'Mean ice diameter'};
%              time_graph = {'Ice number'};     
%              time_graph = {'Ice number size range'};     
%              time_graph = {'Ice mass'};  
%              time_graph = {'LWC mode'}; %the bin that contributes most to the LWC
%            time_graph = {'LWC_dist_CAS_cutoff','LWC_dist_CAS_cutoff','LWC_dist_CAS_cutoff'};
%            time_graph = {'LWC_dist_CAS_cutoff','LWC_dist_CAS_cutoff'};
            

%            instrument={'MAN CAS','BAS CAS','FSSP'};
%            instrument={'MAN CAS','BAS CAS'};            
            instrument={'BAS CAS','BAS CAS'};   
            instrument={'BAS CAS'};               
%            instrument={'MAN CAS'};  
%            instrument={'MAN CAS Karl'};  
%            instrument={'MAN CAS Karl','MAN CDP Karl'};              
%            instrument={'MAN CDP Karl','MAN CAS Karl'};
            
            %exist(x) returns 5 if it is a function and 1 if a variable
             if exist('icas_count')~=1
                icas_count=1;
            end

            %%% these are set in plotTimeHeightVap3.m - might want to override here
            % or to make sure they are the same if requried
            if exist('cut_off_size')~=1
                cut_off_size=1; %size (microns) below which to ignore counts for particle concentration
            end
            if exist('air_speed_type')~=1
                %air_speed_type = 'aircraft';
                %air_speed_type = 'constant 60m/s';
                air_speed_type = 'CIP probe';
            end
            if exist('airspeed_constant')~=1
                airspeed_constant=0;
            end
                
            
%            air_speed_type = 'constant 60m/s';
            if prod(size(CIP_counts_all))>0
                air_speed_type = 'CIP probe';
            else %no CIP probe
                air_speed_type='constant';
                airspeed_constant=100; %assumed 10 m/s for Manchester 19th June tests
                airspeed_constant=17; %measured 17 m/s for Manchester 19th July tests                
            end
        
        fsize=12;    
        y_axis_type='';
        x_axis_type='';
        clear times
        
        ix_distance=0; %flag to plot by distance along flight track rather than time    indsCAS = 1:length(CAS_time_all);

        else 
               clear man_choose_flt_graph  %reset the flag
        end
            
                    
        
        try
            eval(['time_flt=time_flt' flight_no ';']);
        end    
        
        switch air_speed_type
            case 'constant'
                disp('**** WARNING - USING CONSTANT AIRSPEED ****');
        end
            
            idatetick=1; %flag to say the want the xaxis in proper time format rather than decimal time
            datetick_type=15; %specify the type with datetick_type (see help datetick) 15= HH:MM 13=HH:MM:SS
            
           
            
            
            
            iMAN_CAS=0; %flag to choose to use the manchester CAS data in PACS format

            comparison_test=3;            
            switch comparison_test
                case 1
                    %10:28 flights
                    switch instrument{i_instr}
                        case 'BAS CAS'
                            times_mean = [11+26/60+24/3600 11+27/60+48/3600]; %BAS
                        case 'MAN CAS'
                            times_mean = [10+28/60+16/3600 10+29/60+3/3600];  %MAN
                        case 'FSSP'
                            %%%
                    end

                case 2
                    %11:08 flights
                    if iMAN_CAS==0
                        times_mean = [12+8/60+21/3600 12+10/60+48/3600]; %BAS
                    else
                        times_mean = [11+7/60+53/3600 11+8/60+39/3600];  %MAN
%                        times_mean = [11+8/60+0/3600 11+8/60+39/3600];  %MAN
                    end
                    
                case 3
                    %11:08 flights
                    if iMAN_CAS==0 %BAS file 5
                        times_mean = [13+46/60+2/3600 13+47/60+50/3600]; %BAS
                    else %MAN file 6 (last file)
                        times_mean = [12+48/60+27/3600 12+50/60+54/3600];  %MAN
%                        times_mean = [11+8/60+0/3600 11+8/60+39/3600];  %MAN
                    end    

            end
            

            
            
if exist('dat_flt')            
if prod(size(dat_flt))>0
    if exist('times')==1
        if length(times>0)       
            [inds0 inds1]=findheight_nearest(dat_flt(:,1)/1000/3600/24,times(1),times(2));
            inds=[inds0:inds1];
        else
            clear times
            inds=1:length(dat_flt(:,1));
        end
    else
        inds=1:length(dat_flt(:,1));
    end
end
end

            
for idat=1:length(time_graph)
    
    if ~exist('ioverride_lwc_cutoffs') | ioverride_lwc_cutoffs==0
    
    %%%    ------------------------------------   %%%
    %size outside of which to ignore LWC from CAS  
            CAS_LWC_cut_off_sizes = [0 15]; %    
            CAS_LWC_cut_off_sizes = [5 15]; %
%            CAS_LWC_cut_off_sizes = [15 25]; %           
 %           CAS_LWC_cut_off_sizes = [25 35]; %these 3 all coincide between teh FSSP and the CAS so are good choices
            %to avoid one or the other from missing numbers from part of the size range
            %(assuming they're all sizing ok....)
%            CAS_LWC_cut_off_sizes = [35 50];   %50 um doesn't coincide (47 micron max for FSSP)  
            CAS_LWC_cut_off_sizes = [0 20]; %
            
            CAS_LWC_cut_off_sizes = [0 10.2];
%            CAS_LWC_cut_off_sizes = [10.2 50];            
            CAS_LWC_cut_off_sizes = [6.5 7.2];
            CAS_LWC_cut_off_sizes = [2 15];            
            CAS_LWC_cut_off_sizes = [30 50]; %
%            CAS_LWC_cut_off_sizes = [10 30]; %
%            CAS_LWC_cut_off_sizes = [0.61 10]; %
            
%should actually normalise by the bin widths used...

    end


    
    if ix_distance==1
        eval(['X_flt = X_flt' flight_no ';']);
        eval(['Y_flt = Y_flt' flight_no ';']);
        dist_flt = [0; cumsum(sqrt( (diff(X_flt)).^2 + (diff(Y_flt)).^2 ))]; %distance from position (X_pos,Y_pos);
    end
    


    if exist('override_binned_options') &  override_binned_options==1
    else

        %            ismooth(idat)=1;
        %            nsmooth_steps=5; %number of data points over which to smooth data (running average) if smooth(idat)=1
                    nsmooth_steps=300;
                    nsmooth_steps=30;
%                    nsmooth_steps=30;     %using for LWC
        %nsmooth_steps=1;
    end

            
            
            cut_off_size = 1; %min diameter cut off for total number from CAS (microns)
%            cut_off_size = 5; %min diameter cut off for total number from CAS (microns)   
            
            

            
        
    
            
            
                
        switch instrument{idat}
                case 'BAS CAS'
                                                       
         %   if exist('dat_flt')                  
            
           %get the sample volume and total concentrations, plus air speed if required.           
%             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
%                 ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip...
%                 ,CAS_mode_diameter,CAS_mean_diameter,LWC_dist_cas_cutoff]...
%                 =cas_sample_volume_and_stats(dat_flt,CAS_time_all,...
%                 CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all...
%                 ,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes);   
            
            time_timeseries = CAS_time_all;
            
            [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                ,CAS_total_number_cutoff ...                                    
                ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                        ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off,LWC_dist_cas_cutoff_LOW,LWC_dist_cas_cutoff_HIGH]...
                =cas_sample_volume_and_stats2...
                (dat_flt,time_timeseries,...
                CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,airspeed_constant); 
            
            

            
        case 'MAN CAS'

            
            if exist('data_CAS_PACS')
                data_particle = data_CAS_PACS(41:41+29,:)';
                time_timeseries = data_CAS_PACS(1,:);

                [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_PACS',data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
            else
                disp('*** Load in the MAN CAS data! ***')                
            end
            
                                   
        case 'FSSP'
            disp('*** WARNING - applying FSSP flow speed correction factor ***');
            data_particle = data_FSSP(23:42,:)'/3.2906;
%            sample_volume_FSSP = 100*airspeed_constant * 0.24e-2; %0.24e-2 is the CAS laser area in cm^2
            
            time_timeseries = time_FSSP*3600;
            
            %just do at first to calculate the sample volume applied
            [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_FSSP',data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                
                %now compensate for the applied sample volume
                %since we already have concentrations for FSSP
                data_particle = data_particle.*sample_volume_CAS; 
                
             %and recalculate the other stuff   
             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_FSSP',data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                
                
            
          case 'Welas'

            data_particle = 1e-6*welas_dat.conc1';
%            sample_volume_FSSP = 100*airspeed_constant * 0.24e-2; %0.24e-2 is the CAS laser area in cm^2
            
            time_timeseries = welas_dat.time_of_day;
            bins_instrument = welas_dat.size(:,1);
            
            %just do at first to calculate the sample volume applied
            [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_instrument,data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                
                %now compensate for the applied sample volume
                %since we already have concentrations for FSSP
                data_particle = data_particle.*sample_volume_CAS; 
                
             %and recalculate the other stuff   
             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                    ,CAS_total_number_cutoff ...                                    
                    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (0,time_timeseries,...
                    bins_instrument,data_particle,[],[],[],air_speed_type,cut_off_size,[]...
                    ,CAS_LWC_cut_off_sizes,airspeed_constant);
                
            case 'MAN CAS Karl'
                
                time_timeseries = CAS_time_Karl; %seconds from 0 UTC
                data_particle = CAS_per_cc;              
                bins_instrument = CAS_bins_Karl;
                                                               
                sample_volume_CAS = ones(size(data_particle));
                
             case 'MAN CDP Karl'
                
                time_timeseries = CDP_time_Karl; %seconds from 0 UTC
                data_particle = CDP_per_cc;
                sample_volume_CDP = ones(size(data_particle));
                bins_instrument = CDP_bins_Karl;
                    
                
      
        end %switch instrument{idat}
        
        
        if exist('times')==1 %otherwise it might exist as a function
            indsCAS = find(time_timeseries/3600/24>=times(1) & time_timeseries/3600/24<=times(2));
            indsCIP = find(CIP_time_all/3600/24>=times(1) & CIP_time_all/3600/24<=times(2));
        else
            indsCAS = 1:length(time_timeseries);
            indsCIP = 1:length(CIP_time_all);
        end
        
                      
        
                            
                switch time_graph{idat}
                    
                    case 'LWC CAS Karl'
                        
                         iaxis_square=0;
                        
                        %                icas_count=1;
                        ylab=['LWC (g m^{-3})'];
                        

                         ydat(idat).y = CAS_LWC_3to50;

                         
%                        disp('**** WARNING - APPLYING 0.4 factor ****');
                        
                        xdat(idat).x = time_timeseries/3600;
      
                        labs(idat).l=[instrument{idat} ' for 3-50 \mum'];
                        if ismooth(idat)==1
                            labs(idat).l=[labs(idat).l ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
                        end
                        
                        nmark(idat)=0;
                                                                                                                   
                        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                    
                 case 'LWC CDP Karl'
                        
                         iaxis_square=0;
                        
                        %                icas_count=1;
                        ylab=['LWC (g m^{-3})'];
                        

                         ydat(idat).y = CDP_LWC_3to50;

                         
%                        disp('**** WARNING - APPLYING 0.4 factor ****');
                        
                        xdat(idat).x = time_timeseries/3600;
      
                        labs(idat).l=[instrument{idat} ' for 3-50 \mum'];
                        if ismooth(idat)==1
                            labs(idat).l=[labs(idat).l ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
                        end
                        
                        nmark(idat)=0;
                                                                                                                   
                        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                    
                        
                case 'Total number CDP Karl'
                     
                        iaxis_square=0;
                        
                        %                icas_count=1;
                        ylab=['CDP Total number Karl (cm^{-3})'];
                        

%                         ydat(idat).y = CDP_tot_0_06to50;
                         ydat(idat).y = CDP_tot_3to50;
%                         ydat(idat).y = CDP_tot_5to50;
                         
%                        disp('**** WARNING - APPLYTING 0.4 factor ****');
                        
                        xdat(idat).x = time_timeseries/3600;
      
                        labs(idat).l=[instrument{idat} ' for 3-50 \mum'];
                        if ismooth(idat)==1
                            labs(idat).l=[labs(idat).l ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
                        end
                        
                        nmark(idat)=0;
                                                                                                                   
                        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                        
                case 'Total number CAS Karl'
                     
                        iaxis_square=0;
                        
                       
                        
                        %                icas_count=1;
                        ylab=['CAS Total number Karl (cm^{-3})'];
                        
                        %14th CAS bin =3 um - the 14th per cc value is for <3 um
                        CAS_tot_3to50 = sum(CAS_per_cc(:,15:end),2);

%                         ydat(idat).y = CAS_tot_0_06to50;
                         ydat(idat).y = CAS_tot_3to50;
                         
%                        disp('**** WARNING - APPLYTING 0.4 factor ****');
                        
                        xdat(idat).x = time_timeseries/3600;
      
                        labs(idat).l=[instrument{idat} ' for 3-50 \mum'];
                        if ismooth(idat)==1
                            labs(idat).l=[labs(idat).l ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
                        end
                        
                        nmark(idat)=0;
                                                                                                                   
                        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                                                
                    
                    case 'Total number CAS'
                     
                        i_adjust_cas_LWC=1;   
                         ismooth(idat)=1;
                        
                        LWC_cas_adjust_factor = 1/2.5; %2.5 works very well
                        LWC_cas_adjust_factor = 1/1.4;                        
                                                

                        
                     
                        
                        %                icas_count=1;
                        ylab=['CAS Total number (cm^{-3})'];
%                        ydat(idat).y = CAS_total_number{icas_count};
                         ydat(idat).y = CAS_total_number_cutoff(indsCAS);                         
                         
%                        disp('**** WARNING - APPLYTING 0.4 factor ****');
                        
                            xdat(idat).x = time_timeseries(indsCAS)/3600;                                                      

                        labs(idat).l=[instrument{idat} ' for ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) '\mum'];
                        switch i_adjust_cas_LWC
                            case 1
                                ydat(idat).y = LWC_cas_adjust_factor * ydat(idat).y;
                                labs(idat).l=[labs(idat).l  ' * 1/' num2str(1/LWC_cas_adjust_factor) ' '];
                        end                        
                        if ismooth(idat)==1
                            labs(idat).l=[labs(idat).l ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
                        end
                        
                        nmark(idat)=0;
                        
                      
                        
                        
                   
                        
                        

                    lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                    
                    case 'MVD cutoff'
                        ismooth(idat)=0;
                        nsmooth_steps=2;
                        
                        
                        %                icas_count=1;
                        ylab=['MVD cutoff (\mum)'];
%                        ydat(idat).y = CAS_total_number{icas_count};
                         ydat(idat).y = MVD_cut_off;
                         
%                        disp('**** WARNING - APPLYTING 0.4 factor ****');
                        
                            xdat(idat).x = time_timeseries/3600;
      
                        labs(idat).l=[instrument{idat} ' for ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) '\mum'];
                        if ismooth(idat)==1
                            labs(idat).l=[labs(idat).l ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
                        end
                        
                        nmark(idat)=0;
                        
                   
                        
                        

                    lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                        
                    
                    
                case 'Total number normalised by bin width'
                        %                icas_count=1;
                        ylab=['dN/dlogD (cm^{-3} \mum^{-1})'];
%                        ydat(idat).y = CAS_total_number{icas_count};
                         ydat(idat).y = CAS_total_number_cutoff./( log10(bin_range(2)) - log10(bin_range(1)) );
                         
%                        disp('**** WARNING - APPLYTING 0.4 factor ****');
                        
                            xdat(idat).x = time_timeseries/3600;
      
                        labs(idat).l=['dN/dlogD for ' instrument{idat} ' for >' num2str(CAS_LWC_cut_off_sizes(1)) ' \mum and < ' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];
                        
                        ismooth(idat)=1;
                        nmark(idat)=0;
                        
                   
                        
                        

                    lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                                            
                        
                    case 'Last bin number CAS'
                        %                icas_count=1;
                        ylab=['CAS last bin number conc. (cm^{-3})'];                        
                        ydat(idat).y = CAS_counts_all(:,30)./sample_volume_CAS(:,30);
                        idatetick=0;
                        

                            xdat(idat).x = CAS_time_all/3600;
     
                        labs(idat).l='Total number';    
                        
                        ismooth(idat)=0;
                                           
                        

                    lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                    case 'Second bin number CIP'
                        %                icas_count=1;
                        ylab=['CIP second bin number conc. (cm^{-3})'];
                        ydat(idat).y = CIP_counts_all(:,2)./sample_volume_CIP(:,2);
                        idatetick=0;
                        

                            xdat(idat).x = CIP_time_all/3600;

                        labs(idat).l='Total number';    
                        
                        ismooth(idat)=0;   
                        
                 case 'Number in latter bins for CAS'
                        %                icas_count=1;
                        
                        dlog_bins = repmat(log10(diff(CAS_bins)),[1 size(CAS_counts_all,1)]);
                        ylab=['CAS dN/dlogD (cm^{-3} \mum^{-1})'];                        
                        ydat(idat).y = sum(CAS_counts_all(:,29:30)./sample_volume_CAS(:,29:30) ./ dlog_bins(28:29,:)' ,2);                        
%                       ydat(idat).y = sum(CAS_counts_all(:,29:30)./sample_volume_CAS(:,29:30) ,2);
                         
                        idatetick=0;
                        

                            xdat(idat).x = CAS_time_all/3600;

                        labs(idat).l='CAS';    
                        
                        ismooth(idat)=0;
                        
                 case 'Second bin dN/dlogD CIP'
                        %                icas_count=1;
                        ylab=['CIP dN/dlogD (cm^{-3} \mum^{-1})'];     
                        ydat(idat).y = CIP_counts_all(:,2)./sample_volume_CIP(:,2)/log10(25);
                        %CIP bins are 25 microns wide
                        idatetick=0;
                        

                            xdat(idat).x = CIP_time_all/3600;

                        labs(idat).l='CIP';    
                        
                        ismooth(idat)=0;   
                        
                        
                         
                        
                    case 'Total number CIP'
                        %                icas_count=1;
                        ylab=['Total CIP number (cm^{-3})'];
                        ydat(idat).y = CIP_total_number{icas_count};

                            xdat(idat).x = CIP_time_all/3600;

                        labs(idat).l='Total number';
                        
                        ismooth(idat)=0;
                        
                    case 'Total large ice number CIP'
                        %                icas_count=1;
                        thresh_ind=3;
                        D_thresh = 25*thresh_ind-12.5;                        
                        ylab=['Total CIP N_{D>' num2str(D_thresh) '\mum} (L^{-1})'];

                        ydat(idat).y =1000*sum(CIP_counts_all(indsCIP,thresh_ind:62)./sample_volume_CIP(indsCIP,thresh_ind:62),2);
                        %number in bin thresh_ind upwards - converted to per L from per cc
                        xdat(idat).x = CIP_time_all(indsCIP)/3600;

                        labs(idat).l='Total number of large ice';
                        
                        ismooth(idat)=0; 
                        
                        %N.B. - doing numbers like this will be overestimate as one large
                        %ice crystal will get counted as several if it extends sideways

                    case 'LWC_hotwire'
                        ylab=['Liquid water content (g m^{-3})'];
                        ydat(idat).y = LWC_CAS_all(indsCIP)';

                        xdat(idat).x = CIP_time_all(indsCIP)/3600;

                        labs(idat).l='LWC hotwire';
                        
                        ydat(idat).y(ydat(idat).y<-0.5)=NaN;
                        ydat(idat).y(ydat(idat).y>3)=NaN;
                        
                         ismooth(idat)=1;
%                         ismooth(idat)=0; 
                         
                        izlim=1;
                        zmin=0;
                        zmax=0.5;
                        
                        if adjust_STP==1
                            disp('**** WARNING - concentrations are adjusted to STP values *****');
                            %rho_factor was based on CIP_time_Jonny
                            rho_factor_cip = rho_CIP(indsCIP)'/rho_stp;
                            ydat(idat).y = ydat(idat).y./rho_factor_cip;

                            ylab = [ylab ' at STP'];
                        end
                        
                    case 'LWC_dist_CAS_LOW'
                        
                        i_use_cutoff=1;
                        i_adjust_cas_LWC=1;    
                        
                        LWC_cas_adjust_factor = 1/2.5; %2.5 works very well
                        LWC_cas_adjust_factor = 1/1.4;                        
                                                
                        switch i_use_cutoff
                            case 0
                                ydat(idat).y = LWC_dist_cas(indsCAS);
                                labs(idat).l=['LWC CAS'];
                            case 1
                                ydat(idat).y = (1/1)*LWC_dist_cas_cutoff_LOW(indsCAS);
                                labs(idat).l=[instrument{idat}  ' LOW ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];
                        end  
                        
                        switch i_adjust_cas_LWC
                            case 1
                                ydat(idat).y = LWC_cas_adjust_factor * ydat(idat).y;
                                labs(idat).l=[labs(idat).l  ' * 1/' num2str(1/LWC_cas_adjust_factor) ' '];
                        end

                        
                        ylab=['Liquid water content (g m^{-3})'];
                        
                        if adjust_STP==1
                            disp('**** WARNING - CAS LWC adjusted to STP values *****');
                            rho_factor_cas = rho_CAS(indsCAS)' / rho_stp;  %=Vol_stp / Vol_altitude
                            % for constant mass
                            ydat(idat).y = ydat(idat).y./rho_factor_cas;

                            ylab = [ylab ' at STP'];
                        end

                        
                        
                        
                        if iMAN_CAS==1
                            xdat(idat).x = data_CAS_PACS(1,:)'/3600;
                        else
                            xdat(idat).x = CAS_time_all(indsCAS)/3600;
                        end
                        
                         ismooth(idat)=1;
%                          ismooth(idat)=0;
                          
                          icalc_mean=1;
                          switch icalc_mean
                              case 1                                                     

                                  [a,b]=findheight_nearest(xdat(1).x,times_mean(1),times_mean(2));
                                  %N.B. times may be divided by 24 later if using datetick
                                  mean_lwc=mean(ydat(1).y(a:b));
                                  std_lwc = std(ydat(1).y(a:b));
                                  
                                  disp(['mean = ' num2str(mean_lwc) ' std = ' num2str(std_lwc)]);
                                  

                                  
                                  ylab = [ylab ' (Mean of ' num2str(mean_lwc) ' for ' datestr(times_mean(1)/24,13) ' - ' datestr(times_mean(2)/24,13)  ')'];
                                  
                          end
                          
                          if iMAN_CAS==1
                              ylab = [ylab ' MAN'];
                              labs(idat).l = [labs(idat).l ' (MAN)'];
                          end
                          
                    case 'LWC_dist_CAS_HIGH'
                        
                        i_use_cutoff=1;
                        i_adjust_cas_LWC=1;    
                        
                        LWC_cas_adjust_factor = 1/2.5; %2.5 works very well
                        LWC_cas_adjust_factor = 1/1.4;                        
                                                
                        switch i_use_cutoff
                            case 0
                                ydat(idat).y = LWC_dist_cas(indsCAS);
                                labs(idat).l=['LWC CAS'];
                            case 1
                                ydat(idat).y = (1/1)*LWC_dist_cas_cutoff_HIGH(indsCAS);
                                labs(idat).l=[instrument{idat}  ' HIGH ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];
                        end  
                        
                        switch i_adjust_cas_LWC
                            case 1
                                ydat(idat).y = LWC_cas_adjust_factor * ydat(idat).y;
                                labs(idat).l=[labs(idat).l  ' * 1/' num2str(1/LWC_cas_adjust_factor) ' '];
                        end

                        
                        ylab=['Liquid water content (g m^{-3})'];
                        
                        if adjust_STP==1
                            disp('**** WARNING - CAS LWC adjusted to STP values *****');
                            rho_factor_cas = rho_CAS(indsCAS)' / rho_stp;  %=Vol_stp / Vol_altitude
                            % for constant mass
                            ydat(idat).y = ydat(idat).y./rho_factor_cas;

                            ylab = [ylab ' at STP'];
                        end

                        
                        
                        
                        if iMAN_CAS==1
                            xdat(idat).x = data_CAS_PACS(1,:)'/3600;
                        else
                            xdat(idat).x = CAS_time_all(indsCAS)/3600;
                        end
                        
                         ismooth(idat)=1;
%                          ismooth(idat)=0;
                          
                          icalc_mean=1;
                          switch icalc_mean
                              case 1                                                     

                                  [a,b]=findheight_nearest(xdat(1).x,times_mean(1),times_mean(2));
                                  %N.B. times may be divided by 24 later if using datetick
                                  mean_lwc=mean(ydat(1).y(a:b));
                                  std_lwc = std(ydat(1).y(a:b));
                                  
                                  disp(['mean = ' num2str(mean_lwc) ' std = ' num2str(std_lwc)]);
                                  

                                  
                                  ylab = [ylab ' (Mean of ' num2str(mean_lwc) ' for ' datestr(times_mean(1)/24,13) ' - ' datestr(times_mean(2)/24,13)  ')'];
                                  
                          end
                          
                          if iMAN_CAS==1
                              ylab = [ylab ' MAN'];
                              labs(idat).l = [labs(idat).l ' (MAN)'];
                          end   
                          
                          %override some values for line plotting
                            lwidth=[3 3 1.5 1.5];
                            override_patterns=1;
                            patterns(1).p='-';
                            patterns(2).p=':';
                            patterns(3).p='-';
                            patterns(4).p='-';
                        
                

                    case 'LWC_dist_CAS'
                        
                        i_use_cutoff=1;
                        i_adjust_cas_LWC=1;    
                        
                        LWC_cas_adjust_factor = 1/2.5; %2.5 works very well
                        LWC_cas_adjust_factor = 1/1.4;                        
                                                
                        switch i_use_cutoff
                            case 0
                                ydat(idat).y = LWC_dist_cas(indsCAS);
                                labs(idat).l=['LWC CAS'];
                            case 1
                                ydat(idat).y = (1/1)*LWC_dist_cas_cutoff(indsCAS);
                                labs(idat).l=[instrument{idat}  ' ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];
                        end  
                        
                        switch i_adjust_cas_LWC
                            case 1
                                ydat(idat).y = LWC_cas_adjust_factor * ydat(idat).y;
                                labs(idat).l=[labs(idat).l  ' * 1/' num2str(1/LWC_cas_adjust_factor) ' '];
                        end

                        
                        ylab=['Liquid water content (g m^{-3})'];
                        
                        if adjust_STP==1
                            disp('**** WARNING - CAS LWC adjusted to STP values *****');
                            rho_factor_cas = rho_CAS(indsCAS)' / rho_stp;  %=Vol_stp / Vol_altitude
                            % for constant mass
                            ydat(idat).y = ydat(idat).y./rho_factor_cas;

                            ylab = [ylab ' at STP'];
                        end

                        
                        
                        
                        if iMAN_CAS==1
                            xdat(idat).x = data_CAS_PACS(1,:)'/3600;
                        else
                            xdat(idat).x = CAS_time_all(indsCAS)/3600;
                        end
                        
                         ismooth(idat)=1;
%                          ismooth(idat)=0;
                          
                          icalc_mean=1;
                          switch icalc_mean
                              case 1                                                     

                                  [a,b]=findheight_nearest(xdat(1).x,times_mean(1),times_mean(2));
                                  %N.B. times may be divided by 24 later if using datetick
                                  mean_lwc=mean(ydat(1).y(a:b));
                                  std_lwc = std(ydat(1).y(a:b));
                                  
                                  disp(['mean = ' num2str(mean_lwc) ' std = ' num2str(std_lwc)]);
                                  

                                  
                                  ylab = [ylab ' (Mean of ' num2str(mean_lwc) ' for ' datestr(times_mean(1)/24,13) ' - ' datestr(times_mean(2)/24,13)  ')'];
                                  
                          end
                          
                          if iMAN_CAS==1
                              ylab = [ylab ' MAN'];
                              labs(idat).l = [labs(idat).l ' (MAN)'];
                          end
                          
                          
                          
                        
                    case 'LWC_dist_CAS_cutoff'
                        
                        i_use_cutoff=1;
                                                
                        switch i_use_cutoff
                            case 0
                                ydat(idat).y = LWC_dist_cas;
                                labs(idat).l=['LWC CAS'];
                            case 1
                                ydat(idat).y = (1/1)*LWC_dist_cas_cutoff;
%                                labs(idat).l=['LWC CAS ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];
                                labs(idat).l=instrument{idat};
                        end
                        
                        ylab=['Liquid water content (g m^{-3}) ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];                        
                        xdat(idat).x = time_timeseries/3600;  
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

                        
                        
                   case 'LWC_dist_CIP'
                        ylab=['Liquid water content (g m^{-3})'];
                        ydat(idat).y = LWC_dist_cip;

                            xdat(idat).x = CIP_time_all/3600;

                        labs(idat).l='LWC CIP';  
                        
                  case 'Total number CPC'

                        ylab=['Total CPC number (cm^{-3})'];
                        ydat(idat).y = CPC_conc';

                            xdat(idat).x = CPC_time;

                        labs(idat).l='Total number';
                        
                        
                  case 'Reject diags'

                        ylab=['Counts'];
                        
                        idat2=1;
                        ydat(idat2).y = stats_CAS_all(:,10);                        
                        labs(idat2).l='Sum of particles';
                        
                        idat2=idat2+1;
                        ydat(idat2).y = stats_CAS_all(:,11);                       
                        labs(idat2).l='Sum of transit';
                        
                        idat2=idat2+1;
                        ydat(idat2).y = sum(CAS_counts_all,2);
                        labs(idat2).l='Sum of CAS bins';
                        
                        for idat2=1:length(ydat)

                                xdat(idat2).x = CAS_time_all/3600;

                        end
                        
                        y_axis_type='log10_matlab';
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        
                 case 'PADS spreadsheet'
                        ylab=['LWC (g m^{-3})'];
                        %run read_PACS_LWC.m and read_PACS_CAS.m (not actually necessary)
                        %but needed for lwc_PADS_comparison.m, which also needs to be run
                        
                        
                        ydat(idat).y = data_LWC_PACS(7,:); 
                        xdat(idat).x = PACS_TAS_time/3600;   
                        labs(idat).l=['PADS spreadsheet values'];
                        
                        ismooth(idat)=1;

                 case 'PADS spreadsheet Dan''s values'
                        ylab=['LWC (g m^{-3})'];
                        %run read_PACS_LWC.m and read_PACS_CAS.m (not actually necessary)
                        %but needed for lwc_PADS_comparison.m, which also needs to be run
                        
                        ydat(idat).y = LWC_calculated; 
                        xdat(idat).x = PACS_TAS_time/3600;   
                        labs(idat).l=['Dan''s values'];
                        
                        izlim=1;
                        zmin=0;
                        zmax=2;
                        
                        ismooth(idat)=1;
                        
                        

                    lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                                                
                    
                 case 'Mean ice diameter'
                        %                icas_count=1;
                        ylab=['Mean ice size (microns)'];
                        [i1,i2]=findheight_nearest(CIP_time_Jonny/3600,min(time_flt),max(time_flt));
                        inds=[i1:i2];
                        
                        ydat(idat).y = mean_ice_size(inds);

%                        disp('**** WARNING - APPLYTING 0.4 factor ****');

                        xdat(idat).x = CIP_time_Jonny(inds)/3600;

                        labs(idat).l=['Mean ice size'];

                        ismooth(idat)=0;
                        
                case 'Ice number'
                        %                icas_count=1;
                        ylab=['Ice number (L^{-1})'];
                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                        
                        ydat(idat).y = 1000*ice_no_Jonny(inds2);

%                        disp('**** WARNING - APPLYING 0.4 factor ****');

                        xdat(idat).x = CIP_time_Jonny2(inds2)/3600;

                        labs(idat).l=['Total ice number'];

                        ismooth(idat)=0;   
                        
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;

                        ice_accepts = sum(ice_no_particles_PSD(:,inds2),1);

                        num_dat(idat).dat = ice_accepts;
                        iadd_nums_above=1;
                    
                        
                    i0=find(num_dat(idat).dat==0);
%                    Aerr = ydat(idat).y.*sqrt(1./num_dat(1).dat +0.04); %combined absolute error due to Poisson and 20% in TAS
                    Aerr = ydat(idat).y.*sqrt(1./num_dat(idat).dat'); %combined absolute error due to Poisson and 20% in TAS                    
                    SV_mean = TAS_Jonny*160/1000; %mean sample volume in litres assuming 160 mm2 for area
                    NT=1;
                    Aerr(i0)= 1.8./SV_mean(i0)/NT; %set absolute error to 1.8 crystals when have no crystals
                    %then work out an error in the concentration from the mean sample volume                    
                                    
                    error_bar(1).dat=Aerr;
                    
                    iplot_error_bar=1;
                    
 

                case 'Round Number (Jonny)'
                        %                icas_count=1;
                        ylab=['Round ice number (L^{-1})'];
                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                        
                        

%                        disp('**** WARNING - APPLYING 0.4 factor ****');

                        

                        labs(idat).l=['Total round number'];

                        ismooth(idat)=0;   
                        
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;

                        xdat(idat).x = CIP_time_Jonny2(inds2)/3600;
                        ydat(idat).y = 1000*round_no_Jonny(inds2);
                        ice_accepts = sum(round_no_particles_PSD(:,inds2),1);
                        meanvals_tas=TAS_Jonny(inds2);  
                       
                        NT=30; %no. secs
                        iadd_nums_above=1;
                        iplot_error_bar=1;
                                                
                        if NT==0                            
                            num_dat(1).dat = ice_accepts;
                        else
                            time_bins = [xdat(1).x(1):NT/3600:xdat(1).x(end)];
                            [meanvals_Nice,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(xdat(1).x,ydat(1).y,time_bins);                           
                            [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(xdat(1).x,ice_accepts,time_bins);
                            num_dat(1).dat = sum_vals;
                            %for TAS
                            [meanvals_tas,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(xdat(1).x,tas_dat,time_bins);
                            
                            ydat(1).y = meanvals_Nice;
                            xdat(1).x=mid_points;
                        end
                                                                                              
                        i0=find(num_dat(1).dat==0);
%                        Aerr = meanvals_Nice.*sqrt(1./num_dat(1).dat +0.04); %combined absolute error due to Poisson and 20% in TAS
                        Aerr = ydat(1).y.*sqrt(1./num_dat(1).dat); %combined absolute error due to Poisson and 20% in TAS
                        SV_mean = meanvals_tas*160/1000; %mean sample volume in litres assuming 160 mm2 for area
                        Aerr(i0)= 1.8./SV_mean(i0)/NT; %set absolute error to 1.8 crystals when have no crystals
                        %then work out an error in the concentration from the mean sample volume    
                                                                                         
                        error_bar(1).dat=Aerr;
                        
                        
                        ylab = [num2str(NT) ' sec ' ylab];

                        
                        



                 case 'Total number size range (Jonny)'
                        %                icas_count=1;
                        ice_tot = 'ice';
%                        ice_tot = 'tot';
                        

                       
                        [i1,i2]=findheight_nearest(CIP_time_Jonny/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                                      
                        size_range_Jonny = [87.5 1600];
                        size_range_Jonny = [362.5 1600];     
%                        size_range_Jonny = [162.5 1600];                             
%                        size_range_Jonny = [0 1600];                        
                        i1 = find(CIP_size_bins_Jonny>=size_range_Jonny(1));
                        i2 = find(CIP_size_bins_Jonny<=size_range_Jonny(2))-1;                        
                        
                        inds_size_Jonny=[i1(1):i2(end)]; %for the PSD array
                        actual_sizes=num2str(CIP_size_bins_Jonny([i1(1) i2(end)]));
                        switch ice_tot
                            case 'ice'
                                ylab=['Ice number (L^{-1})'];
                                ydat(idat).y = 1000*sum(ice_PSD(inds_size_Jonny,inds2),1);
                            case 'tot'
                                ylab=['Total number (L^{-1})'];
                                ydat(idat).y = 1000*sum(tot_PSD(inds_size_Jonny,inds2),1);
                        end
                        
                        ice_accepts = sum(ice_no_particles_PSD(inds_size_Jonny,inds2),1);

%                        disp('**** WARNING - APPLYING 0.4 factor ****');

                        xdat(idat).x = CIP_time_Jonny(inds2)/3600;

                        labs(idat).l=[ylab actual_sizes(1,:) ' to ' actual_sizes(2,:) ' \mum'];

                        ismooth(idat)=0;  
                        
                        NT=5; %no. secs
                        time_bins = [xdat(1).x(1):NT/3600:xdat(1).x(end)];
                        [meanvals_Nice,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(xdat(1).x,ydat(1).y,time_bins);
                        [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(xdat(1).x,ydat(1).y,time_bins);
                        num_dat(1).dat = sum_vals;
                        
                        
                        iadd_nums_above=1;
                        num_dat(1).dat = ice_accepts;

                         
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;

            case 'Binned ice number size range (Jonny)'
                        %                icas_count=1;
                        ice_tot = 'ice';
%                        ice_tot = 'tot';
%                        ice_tot = 'round';       
                        

                    if exist('override_binned_options') &  override_binned_options==1
                    else
                        nsmooth_steps=60;
%                        nsmooth_steps=5;
                        %smoothed_plot=0; %smooth in windows and plot the error above and below

                        ismooth(idat)=1;
                        NT=0;
                        %                            NT=60; %no. secs
                        %                            NT=30; %no. secs
                    end
        
                         
                       
                        [i1,i2]=findheight_nearest(CIP_time_Jonny/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                                      
                        size_range_Jonny = [87.5 1600];  
                        size_range_Jonny = [112.5 1600];                          
%                        size_range_Jonny = [0 1600];   
%                        size_range_Jonny = [362.5 1600]; %min size draft
%                        size_range_Jonny = [387.5 1600]; %min size draft                        
%                        size_range_Jonny = [212.5 1600];
%                        size_range_Jonny = [162.5 1600];
%                        size_range_Jonny = [512.5 1600]; %min size draft

                        i1 = find(CIP_size_bins_Jonny>=size_range_Jonny(1));
                        i2 = find(CIP_size_bins_Jonny<=size_range_Jonny(2))-1;                        
                        
                        inds_size_Jonny=[i1(1):i2(end)]; %for the PSD array
                        actual_sizes=num2str(CIP_size_bins_Jonny([i1(1) i2(end)]));
                        switch ice_tot
                            case 'ice'
                                ylab=['Binned ice number (L^{-1})'];
                                ice_dat = 1000*sum(ice_PSD(inds_size_Jonny,inds2),1);
                                ice_accepts = sum(ice_no_particles_PSD(inds_size_Jonny,inds2),1);                                
                            case 'round'
                                ylab=['Number of round particles (L^{-1})'];
                                ice_dat = 1000*sum(round_PSD(inds_size_Jonny,inds2),1);
                                ice_accepts = sum(round_no_particles_PSD(inds_size_Jonny,inds2),1);                                
                            case 'tot'
                                ylab=['Total number (L^{-1})'];
                                ice_dat = 1000*sum(tot_PSD(inds_size_Jonny,inds2),1);
                                ice_accepts = sum(tot_no_particles_PSD(inds_size_Jonny,inds2),1);
                        end
                        


%                        disp('**** WARNING - APPLYING 0.4 factor ****');

                        CIP_time_plot = CIP_time_Jonny(inds2)/3600;

                        labs(idat).l=[ylab actual_sizes(1,:) ' to ' actual_sizes(2,:) ' \mum'];

                        tas_dat = TAS_Jonny(inds2);    
                        

                        
                 
                        

                        
                        if NT==0
                            num_dat(1).dat = ice_accepts;
                            mid_points = CIP_time_plot;
                            meanvals_Nice = ice_dat;
                            meanvals_tas = tas_dat;
                        else
                            time_bins = [CIP_time_Jonny(inds2(1))/3600:NT/3600:CIP_time_Jonny(inds2(end))/3600];
                            [meanvals_Nice,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny(inds2)/3600,ice_dat,time_bins);
                            [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny(inds2)/3600,ice_accepts,time_bins);
                            num_dat(1).dat = sum_vals;
                            [meanvals_tas,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny(inds2)/3600,tas_dat,time_bins);
                        end

                        xdat(idat).x = mid_points;
                        ydat(idat).y = meanvals_Nice;

                        if smoothed_plot==1
                            %doing a sum over each window instead of a mean
                            nfilter=nsmooth_steps;
                            [xdat_smooth,Ncry]=window_average(mid_points,num_dat(1).dat,nfilter,'sum');

                            %                                 bfilter=ones([1 nfilter])*1;
                            %
                            %                                 Ncry=filter(bfilter,1,num_dat(1).dat);
                            %                                 %ydat(idat).y(end-nfilter+1:end)=[]; xdat(idat).x(end-nfilter+1:end)=[];
                            %                                 Ncry(1:nfilter-1)=[];
                            %                                 if mod(nfilter,2)==0  %even number
                            %                                     xdat_smooth =0.5*(mid_points(nfilter/2:end-nfilter/2) + mid_points(nfilter/2+1:end-nfilter/2+1));
                            %                                 else
                            %                                     xdat_smooth = mid_points(nfilter/2+0.5:end-nfilter/2+0.5);
                            %                                     %                        xdat(idat).x(1:nfilter)=[];
                            %                                 end

                            %now do the actual mean of the
                            %concentrations
                            %                                 bfilter2=ones([1 nfilter])*1/nfilter;
                            %                                 ydat(2).y = filter(bfilter2,1,meanvals_Nice);
                            %                                 ydat(2).y(1:nfilter-1)=[];
                            
                            [xdat_smooth,ydat(2).y]=window_average(mid_points,meanvals_Nice,nfilter,'mean');
                            [xdat_smooth,tas_smoothed]=window_average(mid_points,meanvals_tas,nfilter,'mean');

                        
                            %override some values for line plotting
                            lwidth=[1.5 3 1.5];
                            override_patterns=1;
                            patterns(1).p=':';
                            patterns(2).p='-';
                            patterns(3).p=':';


                            %%%
                            %ydat(2).y = Ncry;
                            %%%
                            error = ydat(2).y .* sqrt(Ncry)./Ncry; %absolute error



                            i0=find(Ncry==0);
                            SV_mean = tas_smoothed*160/1000; %mean sample volume in litres assuming 160 mm2 for area
                            error(i0)= 1.8./SV_mean(i0)/nfilter; %set absolute error to 1.8 crystals when have no crystals
                            %then work out an error in the concentration from
                            %the mean sample volume (1.8 comes from
                            %Poisson stats)

                            ydat(1).y = ydat(2).y - error;
                            ydat(1).y(i0) = 0; %no lower bound when have zero counts
                            ydat(3).y = ydat(2).y + error;


                            xdat(1).x = xdat_smooth;
                            xdat(2).x = xdat_smooth;
                            xdat(3).x = xdat_smooth;

                            labs(2).l = labs(1).l;
                            labs(1).l = 'Lower error bound';
                            labs(3).l = 'Upper error bound';
                            

                        end


              
                        
                        
                       
                        if smoothed_plot==0
                             ylab = [num2str(NT) ' sec ' ylab];
                             iadd_nums_above=1;
                             iplot_error_bar=1;
                        else
                            ylab = [num2str(nsmooth_steps) ' sec ' ylab];
                        end
                        
                        if length(strfind(file_name_h5,'CENTRE'))>0
                            ylab=[ylab ' CENTRE-IN'];
                        else
                            ylab=[ylab ' ALL-IN'];
                        end
                        
                        if adjust_STP==1
                            ylab=[ylab ' (STP adjusted)'];
                        end
                        
                        i0=find(num_dat(1).dat==0);
%                        Aerr = meanvals_Nice.*sqrt(1./num_dat(1).dat +0.04); %combined absolute error due to Poisson and 20% in TAS
                        Aerr = meanvals_Nice.*sqrt(1./num_dat(1).dat); %combined absolute error due to Poisson and 20% in TAS                        
                        SV_mean = meanvals_tas*160/1000; %mean sample volume in litres assuming 160 mm2 for area
                        Aerr(i0)= 1.8./SV_mean(i0)/NT; %set absolute error to 1.8 crystals when have no crystals
                        %then work out an error in the concentration from the mean sample volume    
                                                                                         
                        error_bar(1).dat=Aerr;

                                                                                                
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;

              case 'N crystals per bin ice number size range (Jonny)'
                                                ice_tot = 'ice';
%                        ice_tot = 'tot';
                        

                       
                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                                      
                        size_range_Jonny = [87.5 1600];
                        i1 = find(CIP_size_bins_Jonny>=size_range_Jonny(1));
                        i2 = find(CIP_size_bins_Jonny<=size_range_Jonny(2))-1;                        
                        
                        inds_size_Jonny=[i1(1):i2(end)]; %for the PSD array
                        actual_sizes=num2str(CIP_size_bins_Jonny([i1(1) i2(end)]));
                        switch ice_tot
                            case 'ice'
                                ylab=['N crystals per bin'];
                                ice_dat = 1000*sum(ice_PSD(inds_size_Jonny,inds2),1);
                            case 'tot'
                                ylab=['Total number (L^{-1})'];
                                ice_dat = 1000*sum(tot_PSD(inds_size_Jonny,inds2),1);
                        end
                        
                        ice_accepts = sum(ice_no_particles_PSD(inds_size_Jonny,inds2),1);

%                        disp('**** WARNING - APPLYING 0.4 factor ****');

                        CIP_time_plot = CIP_time_Jonny2(inds2)/3600;

                        labs(idat).l=[ylab actual_sizes(1,:) ' to ' actual_sizes(2,:) ' \mum'];

                        ismooth(idat)=0;  
                        
                        NT=10; %no. secs
                        time_bins = [CIP_time_Jonny2(inds2(1))/3600:NT/3600:CIP_time_Jonny2(inds2(end))/3600];
                        [meanvals_Nice,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny2(inds2)/3600,ice_dat,time_bins);
                        [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny2(inds2)/3600,ice_accepts,time_bins);
                        
                        xdat(idat).x = mid_points;
                        ydat(idat).y = sum_vals;
                        
                        ylab = [num2str(NT) ' sec bins - ' ylab];
                        
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;

              case 'Number of ice accepts (Jonny)'
                          %                icas_count=1;
                        ylab=['Number'];
                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                                      
                        size_range_Jonny = [137.5 1600];
                        i1 = find(CIP_size_bins_Jonny>=size_range_Jonny(1));
                        i2 = find(CIP_size_bins_Jonny<=size_range_Jonny(2))-1;                        
                        
                        inds_size_Jonny=[i1(1):i2(end)]; %for the PSD array
                        actual_sizes=num2str(CIP_size_bins_Jonny([i1(1) i2(end)]));
                        
                        ydat(idat).y = sum(ice_no_particles_PSD(inds_size_Jonny,inds2),1);
%                        ydat(idat).y = 1000*sum(ice_no_particles_PSD(inds_size_Jonny,inds2),1);                        

%                        disp('**** WARNING - APPLYING 0.4 factor ****');

                        xdat(idat).x = CIP_time_Jonny2(inds2)/3600;

                        labs(idat).l=['Number particles used for ' actual_sizes(1,:) ' to ' actual_sizes(2,:) ' \mum'];

                        ismooth(idat)=0;   
                        
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;

               case 'Number of ice rejects (Jonny)'
                          %                icas_count=1;
                        ylab=['Number'];
                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt(inds)),max(time_flt(inds)));
                        inds2=[i1:i2];
                                      
                        size_range_Jonny = [137.5 1600];
                        i1 = find(CIP_size_bins_Jonny>=size_range_Jonny(1));
                        i2 = find(CIP_size_bins_Jonny<=size_range_Jonny(2))-1;                        
                        
                        inds_size_Jonny=[i1(1):i2(end)]; %for the PSD array
                        actual_sizes=num2str(CIP_size_bins_Jonny([i1(1) i2(end)]));
                        
                        ydat(idat).y = CIP_NRej_Jonny(inds2);


                        xdat(idat).x = CIP_time_Jonny2(inds2)/3600;

%                        labs(idat).l=['Number particles used for ' actual_sizes(1,:) ' to ' actual_sizes(2,:) ' \mum'];
                         labs(idat).l='Number of rejects';
                         
                        ismooth(idat)=0;   
                        
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                        
                        izlim=0;
                        zmin=0;
                        zmax=0.5;
%                        zmax=12;


               case 'Ice mass'
                        %                icas_count=1;
                        ylab=['Ice mass (mg m^{-3})'];
                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt),max(time_flt));
                        inds=[i1:i2];
                        
                        ydat(idat).y = 1000*ice_mass_Jonny(inds);

                        xdat(idat).x = CIP_time_Jonny2(inds)/3600;

                        labs(idat).l=['Total ice mass'];

                        ismooth(idat)=0;   
                        
                        
                        lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

         
                        izlim=0;
                        zmin=0;
                        zmax=1;
                        
                        
                     case 'LWC mode'
                        %                icas_count=1;
                        LWC_min=0.2;
                                                
                        ylab=['Mode for LWC.GT.' num2str(LWC_min) ' g m^{-3}'];
                        [i1,i2]=findheight_nearest(CAS_time_all/3600,min(time_flt),max(time_flt));
                        inds=[i1:i2];

                        
                        [LWC_mode_mid LWC_mode_lower LWC_mode_upper,LWC_percent_contribution]=LWC_mode(LWC_size_dist,CAS_bins,LWC_min);
                        
                        ydat(idat).y = LWC_mode_mid(indsCAS);


                        xdat(idat).x = CAS_time_all(indsCAS)/3600;

                        labs(idat).l=[instrument{idat}];

                        ismooth(idat)=0;   
                        
                        nmark=-1;
                        
                        
                        lor=3;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

         
                        izlim=0;
                        zmin=0;
                        zmax=1;
    



                end  %switch time_graph{idat}
                
               
                

                if ismooth(idat)==1
                    
                    if isum_smooth==0
                        nfilter=nsmooth_steps; bfilter=ones([1 nfilter])*1/nfilter;
                        [xdat(idat).x,ydat(idat).y] = window_average(xdat(idat).x,ydat(idat).y,nfilter,'mean');
                    else
                        nfilter=nsmooth_steps; bfilter=ones([1 nfilter])*1;
                        [xdat(idat).x,ydat(idat).y] = window_average(xdat(idat).x,ydat(idat).y,nfilter,'sum');
                    end
                     %note - the way that this was being done before was wrong for a smoothing window
                %the x values weren't being properly assigned as the
                %mid-points of the window. E.g. the last entry in the
                %smoothed data is actually the average of the last nfilter
                %datapoints and should not be ignored as before   
                
% %                 bfilter=ones([1 nfilter])*1/nfilter;
% %                 ydat(idat).y=filter(bfilter,1,ydat(idat).y);
% %                 ydat(idat).y(end-nfilter+1:end)=[]; xdat(idat).x(end-nfilter+1:end)=[];
% %                 ydat(idat).y(1:nfilter)=[]; xdat(idat).x(1:nfilter)=[];
% 
%                     ydat(idat).y=filter(bfilter,1,ydat(idat).y);
%                     %ydat(idat).y(end-nfilter+1:end)=[]; xdat(idat).x(end-nfilter+1:end)=[];
%                     ydat(idat).y(1:nfilter-1)=[];
%                     if mod(nfilter,2)==0  %even number
%                         xdat(idat).x =0.5*(xdat(idat).x(nfilter/2:end-nfilter/2) + xdat(idat).x(nfilter/2+1:end-nfilter/2+1)); 
%                     else
%                         xdat(idat).x = xdat(idat).x(nfilter/2+0.5:end-nfilter/2+0.5);
% %                        xdat(idat).x(1:nfilter)=[];
%                     end
                end
                
                if ix_distance==1
                    xdat(idat).x=interp1(dat_flt(:,1)/1e3/3600,dist_flt,xdat(idat).x);
                    xlab=['Distance along flight track (km)'];
                    idatetick=0;
                else
                    xlab=['UTC Time'];
                end
                
                if idatetick==1 %make sure that all xdat data has been changed in case they weren't created
                    %in the idat loop
                    for idat2=idat:length(xdat)                        
                        xdat(idat2).x=xdat(idat2).x/24; %convert to days
                    end
                end
                
                
                
                

        end   %for idat=1:length(time_graph)
        
        
%    eval(['X_flt = X_flt' flight_no ';']);
%    eval(['Y_flt = Y_flt' flight_no ';']);
                


%        titlenam=[ylab ' for flight ' flight_no ' on ' date_str ' at ' CAS_start_time_str];

if exist('flight_no')
        titlenam=[ylab ' for flight ' flight_no ' on ' date_str];
else
        titlenam=[ylab];
end

if ismooth(1)==1
    titlenam=[titlenam ' SMOOTHING OVER ' num2str(nsmooth_steps) ' data points'];
end

        figname=[titlenam];
        savename=titlenam;
        
        
        


        if ~exist('ioverride_xlims') | ioverride_xlims==0
            xlims=1;       
            xlimits=[min(xdat(1).x) max(xdat(1).x)]; %whole range of the flight
        else
            clear ioverride_xlims
        end
        
        

        maxx=-1e99;
        for idat=1:length(ydat)
            if size(ydat(idat).y,1)>1
                maxx=max([ydat(idat).y' maxx]);
            else
                maxx=max([ydat(idat).y maxx]);
            end
        end

        
   
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'AWS_antjan06'
        %run read_AWS_data.m first (or the Wisconsin script)
        idatetick=1; %flag to say the want the xaxis in proper time format rather than decimal time
        datetick_type=0; %specify the type with datetick_type (see help datetick) 15= HH:MM 13=HH:MM:SS
        
        i_set_dateticks=1;
        date_ticks=[5:0.25:8];

        iaxis_square=0;
            
        nmark=0;
        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        wis_aws=idir; %
        wis_aws=0;  %best to choose wis_aws=0 for Larsen as this is the original (not wisconsin AWS) data
        %wisconsin data seems to be smoothed?? Or lower time res? Be careful about teh ilat ilon locations though
        %choose the Wisconsin data for the smoothed wind speed data
        %BUT don't want to smooth the wind direction - then choose wis_aws=0
        
%        wis_aws=17;

        xlims=1;
        xlimits=[20.1 20.8407]; %
        xlimits=[19.2 22.75]; %
        xlimits=[5 8]; %
%        xlimits=[0 31];
%        xlimits=[32 AWStime_wis(7).dat(end)]; %
        

        ih_met=0;


        %timeseries comparing to AWS data
        %read_AWS_data.m reads in the data
        %fields are
        %1)Year
        %2)Month
        %3)Day
        %4)Hour
        %5)TEMPERATURE
        %6)PRESSURE
        %7)WIND_SPEED
        %8)WIND_DIRECTION
        %HUMIDITY values were all null.....



        % (1) Kirkwood Island - Lat : 68.34S  Long :  69.01W    Elev :   30 M
        % (2) Dismal Island - Lat : 68.09S  Long :  68.82W      Elev :   10 M
        % (3) Bonaparte Point -  Lat : 64.78S  Long :  64.07W   Elev :    8 M
        % (4) Sky Blu   -   Lat : 74.79S  Long :  71.49W        Elev : 1510 M
        % (5) Limbert   -   Lat : 75.91S  Long :  59.26W        Elev :   40 M
        % (6) Butler Island  -  Lat : 72.21S  Long :  60.17W    Elev :   91 M
        % (7) Larsen BAS     -  Lat : 67.01S  Long :  61.55W    Elev :   17 M


             flt_graph = 'Temp';
             flt_graph = 'Temp 2m';             
        %    flt_graph = 'Height';
 %           flt_graph = 'Vapour';
        %    flt_graph = 'Lat';
        %    flt_graph = 'Lon';
%            flt_graph = 'Pressure';
%%%            flt_graph = 'Wind';
            flt_graph = 'Wind 10m';
%            flt_graph = 'Wind dir';
            flt_graph = 'Wind dir 10m';



        incep=0;
        if incep==1

            pressure_corr(1)=0;  %pressure corrections (mb) to take into account the analysis coarse terrain
            pressure_corr(2)=0;  %these are done by sight rather than calculation to asses the pressure trends
            pressure_corr(3)=0;
            pressure_corr(4)=-95;
            pressure_corr(5)=20;  %according to pressure calculation this should be 12.5234 hPa
            pressure_corr(6)=44;
            pressure_corr(7)=0;
        else
            pressure_corr(1)=0;  %pressure corrections (mb) to take into account the analysis coarse terrain
            pressure_corr(2)=0;
            pressure_corr(3)=0;
            pressure_corr(4)=0;
            pressure_corr(5)=-4.897;  %this particular one was calculated from the actual profile data to take inot account the height diff
            %between the height of the coarse met_em terrain and the actual AWS height (40 m). Done by solving hydrostatic
            %equation from soilhgt upwards using PRES and TT array - note new ecmwf with changed Vtable doesn't give SOILHGT.
            %produces a good match - ideally should be done for each time in the timeseries
            pressure_corr(6)=2;
            pressure_corr(7)=0;
        end


        if is_met_em  %set the directory for the AWS positions here - actually makes little difference
            cd(['Y:/WRF/' rundir '/AWS_positions_d01/']);      %whether using d02 vs d03 for analysis
            loadname=['pressure_alt_corr_timser_level_' num2str(ih_met) '.mat'];
            %             loadname=['pressure_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'pressure');
            loadname=['wind_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'wind');
            loadname=['temperature_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'temperature');
            loadname=['winddir_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'winddir');

            if ih_met~=0
                loadname=['vapour_timser_level_' num2str(ih_met) '.mat'];
                load(loadname,'vapour');

            end

            ilat = pressure(wis_aws).ilat;
            ilon = pressure(wis_aws).ilon;


        else
            %Lat : 67.01S  Long :  61.55W  (also says 61.33 on same webiste!)
            %Elev :   17 M  Taken from website - but not
            %sure how accurate - 17m abmsl?

            %LAT=[-67.01]; %location of Larsen AWS
            %LON=[-61.55];

            %LAT=[-67.01 -68.4]; %location of Larsen AWS and point on the west side below Rothera
            %LON=[-61.55 -68.4];


            if wis_aws==0
%                LAT=[-68.34]; %location of Kirkwood AWS
%                LON=[-69.01];
%Larsen AWS
                LAT = -67.0100;
                LON = -61.5500;
  
            else
                LAT = LAT_AWS(wis_aws).dat;
                LON = LON_AWS(wis_aws).dat;
            end


            [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);
        end




        start_time = 5*24; %in hrs
        end_time = 8*24;
        wrf_Tinds = [1:24];
        wrf_Tinds = [1:12];

        wrf_Tinds = [1:size(Times,1)];

        %[ind,ind2] = findheight(AWShrs,start_time,end_time);
        %inds = ind:ind2;

        if i_set_dateticks==1
            xlab = 'Time (UTC)';
        else
            xlab = 'Time (days)';        
        end

        izlim=0;




        switch flt_graph
            case {'Temp','Temp 2m'}

                ylab=['Temperature (^{o}C)'];

                if wis_aws==0
%                    ydat(1).y = AWSdat2(5,inds)';
                    ydat(1).y = AWSdat2(5,:)';
%                    xdat(1).x = AWShrs(inds)'/24;
                    xdat(1).x = AWShrs2(:)'/24; %revised time points
                    labs(1).l='Larsen AWS';
                else
                    ydat(1).y = ( AWSdat_wis(wis_aws).dat(3,~isnan(AWSdat_wis(wis_aws).dat(3,:))) )';
                    xdat(1).x = ( AWStime_wis(wis_aws).dat(~isnan(AWSdat_wis(wis_aws).dat(3,:))) )';
%                    xdat(1).x=datenum(2009,1,1)+xdat(1).x;                    
%                    xlimits=datenum(2009,1,1)+xlimits;
                    labs(1).l=[aws_name(wis_aws).name ' AWS'];
                    izlim=0;
                    zmin=980;
                    zmax=1010;
                end


                switch wis_aws
                    case 1
                        izlim=1;
                        zmin=-5;
                        zmax=7;
                        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                    case 2
                        izlim=1;
                        zmin=-5;
                        zmax=7;
                        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                    case 7
                        izlim=0;
                        zmin=-5;
                        zmax=7;
                        lor=3; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane    

                        
                      
                end


            case 'Height'

                ylab=['Height (km)'];
                ydat(1).y = dat(:,11);
                xdat(1).x = time_flt19;
                labs(1).l='Height';

            case 'Vapour'

                ylab=['Water vapour mixing ratio (g kg^{-1})'];
                ydat(1).y = qv2*1000;
                xdat(1).x = time_flt19;
                labs(1).l='Vapour';

            case 'Lat'

                ylab=['Latitude (degrees)'];
                ydat(1).y = dat(:,2);
                xdat(1).x = time_flt19;
                labs(1).l='Latitude';

            case 'Lon'

                ylab=['Longitude (degrees)'];
                ydat(1).y = dat(:,3);
                xdat(1).x = time_flt19;
                labs(1).l='Longitude';

            case 'Pressure'

                ylab=['Pressure (mb)'];

                if wis_aws==0
                    ydat(1).y = AWSdat2(6,inds)';
                    xdat(1).x = AWShrs(inds)'/24;
                    labs(1).l='Larsen AWS';
                else
                    ydat(1).y = ( AWSdat_wis(wis_aws).dat(4,~isnan(AWSdat_wis(wis_aws).dat(4,:))) )';
                    xdat(1).x = ( AWStime_wis(wis_aws).dat(~isnan(AWSdat_wis(wis_aws).dat(4,:))) )';
                    labs(1).l=[aws_name(wis_aws).name ' AWS'];
                    izlim=1;
                    zmin=980;
                    zmax=1010;

                    lor=3;

                    switch wis_aws
                        case 2
                            zmin=980;
                            zmax=998;
                        case 3
                            zmin=955;
                            zmax=985;
                        case 4
                            zmin=800;
                            zmax=822;
                        case 5
                            zmin=975;
                            zmax=990;
                            lor=1;
                        case 6
                            zmin=950;
                            zmax=990;
                        case 10
                            zmin=950;
                            zmax=990;  
                        case 20
                            zmin=970;
                            zmax=995;      
                    end


                end



            case {'Wind','Wind 10m'}

                ylab=['Wind speed (m s^{-1})'];

                if wis_aws==0
                    ydat(1).y = AWSdat2(7,:)';
%                    xdat(1).x = AWShrs(:)'/24;
                    xdat(1).x = AWShrs2(:)'/24;
                    
                    labs(1).l='Larsen AWS';
                else
                    ydat(1).y = ( AWSdat_wis(wis_aws).dat(5,~isnan(AWSdat_wis(wis_aws).dat(5,:))) )';
                    xdat(1).x = ( AWStime_wis(wis_aws).dat(~isnan(AWSdat_wis(wis_aws).dat(5,:))) )';
                    labs(1).l=[aws_name(wis_aws).name ' AWS'];
                    izlim=1;
                    zmin=0;
                    zmax=15;
                end

                lor = 1;

                switch wis_aws
                    case 5
                        lor=4;
                end   
                
                iexecute_script=1;
                script_name='L_shaped_leg_std_devs';



            case {'Wind dir','Wind dir 10m'}
                
                if wis_aws==0
                    ydat(1).y = AWSdat2(8,:)';
%                    xdat(1).x = AWShrs(:)'/24;
                    xdat(1).x = AWShrs2(:)'/24;
                    labs(1).l='Larsen AWS';
                else                    
                    ydat(1).y = AWSdat_wis(wis_aws).dat(6,~isnan(AWSdat_wis(wis_aws).dat(5,:)))';
                    xdat(1).x = ( AWStime_wis(wis_aws).dat(~isnan(AWSdat_wis(wis_aws).dat(5,:))) )';
                    labs(1).l= [aws_name(wis_aws).name ' AWS'];
                end
                
                %move all the times with direction less than 20 to avoid too much "flipping"
                move_val=50;
                imove=find(ydat(1).y'<move_val & xdat(1).x<7);
                ydat(1).y(imove)=ydat(1).y(imove)+360;
                
                move_val=350;  %make the directions look a bit nicer by selective adding or removing of 360
                imove=find(ydat(1).y'>move_val & xdat(1).x>7);
                ydat(1).y(imove)=ydat(1).y(imove)-360;
                               
                ylab=['Wind direction (degrees)'];
                
                
                iexecute_script=1;
                script_name='L_shaped_leg_std_devs';
                



        end
        
        titlenam=[ylab ' for AWS data'];
        figname=[titlenam ' ' labs(1).l];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wrf/analysis plots %%%%%%%%%%%%%%%%%%%%%%
        for idat=2:length(ilon)+1;
            iloc=idat-1;
            if isnan(ilat(iloc))==1 | isnan(ilon(iloc))==1
                xdat(idat).x=NaN;
                ydat(idat).y=NaN;
                labs(idat).l='No WRF data';
                continue  %continue skips forward to the next iteration in the loop
            end

            switch flt_graph
                case 'Pressure'

                    labs(idat).l = filestr;





                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = pressure(wis_aws).y(i_order)'/100;
                    else

                        i_interp=1;
                        if i_interp==1

                            terr = nc{'HGT'}(1,:);
                            for ipres=1:length(wrf_Tinds)
                                PSFC=nc{'PSFC'}(wrf_Tinds(ipres),:,:)/100;
                                pres=( nc{'P'}(wrf_Tinds(ipres),:,ilat,ilon) + nc{'PB'}(wrf_Tinds(ipres),:,ilat,ilon) )/100;
                                H=WRFUserARW(nc,'Z',wrf_Tinds(ipres),ilat(iloc),ilon(iloc));
                                if ALT_AWS(wis_aws).dat>terr(ilat,ilon)
                                    ydat(idat).y(ipres) = interp1([terr(ilat,ilon) H],[PSFC(ilat,ilon) pres],ALT_AWS(wis_aws).dat);
                                else
                                    ydat(idat).y(ipres) = interp1([terr(ilat,ilon) H],[PSFC(ilat,ilon) pres],ALT_AWS(wis_aws).dat,'linear','extrap');
                                    labs(idat).l = [filestr '_extrapolation from ' num2str(terr(ilat,ilon)) ' m'];
                                end

                            end

                            ydat(idat).y=ydat(idat).y';

                        else

                            %ydat(idat).y = (nc{'P'}(wrf_Tinds,1,ilat(iloc),ilon(iloc)) + nc{'PB'}(wrf_Tinds,1,ilat(iloc),ilon(iloc)) )/100;
                            ydat(idat).y = nc{'PSFC'}(wrf_Tinds,ilat(iloc),ilon(iloc))/100;

                        end

                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
                    end

                    xdat(idat).x = xdat(idat).x';


                case 'Temp'
                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = temperature(wis_aws).y(i_order)'-273.15;
                    else
                        potemp = nc{'T'}(wrf_Tinds,1,ilat(iloc),ilon(iloc)) + 300;
                        P = nc{'P'}(wrf_Tinds,1,ilat(iloc),ilon(iloc)) + nc{'PB'}(wrf_Tinds,1,ilat(iloc),ilon(iloc));
                        ydat(idat).y = potemp ./ ( (1e5./P).^0.286 ) - 273.15;
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
%                        xdat(idat).x=datenum(2009,1,1)+xdat(idat).x;
                        %                    xdat(idat).x = xdat(idat).x';
                    end

                    labs(idat).l = filestr;
                    xdat(idat).x = xdat(idat).x';
                    
                case 'Temp 2m'
                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = temperature(wis_aws).y(i_order)'-273.15;
                    else                       
                        ydat(idat).y = nc{'T2'}(wrf_Tinds,ilat(iloc),ilon(iloc)) - 273.15;
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
%                        xdat(idat).x=datenum(2009,1,1)+xdat(idat).x;
                        %                    xdat(idat).x = xdat(idat).x';
                    end

                    labs(idat).l = [filestr ' at 2m'];
                    xdat(idat).x = xdat(idat).x';    


                case 'vapour'
                    if is_met_em
                        rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                        T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                        P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                        qsat = satvappress(T,'goff','liq',P,1)/f;
                        xdat(i).x = 1000 * rh/100 .* qsat;
                    else
                        xdat(i).x = 1000*nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                    end

                case 'Wind'
                    ih_wrf=1;
                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = wind(wis_aws).y(i_order)';                        
                    else
                        u=0.5* ( nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc)+1,ilon(iloc)) );
                        v=0.5* ( nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)+1) );
                        ydat(idat).y= sqrt( u.^2 + v.^2 );

                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end

                    end

                    xdat(idat).x = xdat(idat).x';
                    labs(idat).l = filestr;    

                    
                 case 'Wind 10m'
                    ih_wrf=1;
                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = wind(wis_aws).y(i_order)';                        
                    else
                        u=nc{'U10'}(wrf_Tinds,ilat(iloc),ilon(iloc));
                        v=nc{'V10'}(wrf_Tinds,ilat(iloc),ilon(iloc));
                        ydat(idat).y= sqrt( u.^2 + v.^2 );

                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end

                    end

                    xdat(idat).x = xdat(idat).x';
                    labs(idat).l = [filestr ' at 10m'];
                    
                    izlim=1;
                    zmin=0;
                    zmax=15;


                case 'Wind dir'
                    ih_wrf=1;
                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = winddir(wis_aws).y(i_order)';
                    else
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
                        xdat(idat).x = xdat(idat).x';



                        u=0.5* ( nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc)+1,ilon(iloc)) );
                        v=0.5* ( nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)+1) );


                        jnorth = ilat(iloc) + 10;
                        lons_north = lon2d.var(jnorth,:);
                        [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

                        %angle of the local north line relative to the grid
                        thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );




                        for iuv=1:length(u)

                            theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                            if u(iuv)==0 & v(iuv)==0
                                ydat(idat).y(iuv) = 0;
                            elseif u(iuv)>=0 & v(iuv)>=0
                                ydat(idat).y(iuv) = theta2;
                            elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                                ydat(idat).y(iuv) = 180 + theta2;
                            elseif u(iuv)<=0 & v(iuv)<=0
                                ydat(idat).y(iuv) = 180 + theta2;
                            elseif u(iuv)<0 & v(iuv)>0
                                ydat(idat).y(iuv) = 360 + theta2; %theta2 is negative
                            end

                        end

                        ydat(idat).y = ydat(idat).y' + 180 - thetaN*180/pi; %add 180 to make it the direction wind is coming from
                        % take away thetaN to give direction relative to north
                        i360 = find(ydat(idat).y>=360);
                        ydat(idat).y(i360) = ydat(idat).y(i360) - 360;
                        
                    end
                        
                        izlim=1;
                        zmin=0;
                        zmax=450;

                        labs(idat).l = filestr;
                        xdat(idat).x = xdat(idat).x';
                        
                case 'Wind dir 10m'

                    if is_met_em
                        [xdat(idat).x,i_order] = sort(pressure(1).times);
                        ydat(idat).y = winddir(wis_aws).y(i_order)';
                    else
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
                        xdat(idat).x = xdat(idat).x';

                       u=nc{'U10'}(wrf_Tinds,ilat(iloc),ilon(iloc));
                       v=nc{'V10'}(wrf_Tinds,ilat(iloc),ilon(iloc));
                       
                       ydat(idat).y=wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat(iloc),ilon(iloc),DX,DY);
                       
                       move_val=2;
                       imove=find(ydat(idat).y<move_val);
                       ydat(idat).y(imove)=ydat(idat).y(imove)+360;
                       
                        
                       labs(idat).l = [filestr ' 10m'];
                       xdat(idat).x = xdat(idat).x'; 

                    end
                   
                case 'cloud'
                    cloud=nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                    xdat(i).x=1000*cloud;


                    figname=['Cloud mixing ratio profile at ' tstr ' for ' filestr];
                    xlab='Total condensed water mixing ratio (g kg^{-1})';


            end

            if is_met_em==1
                labs(idat).l=[labs(idat).l ' analysis'];
            end


        end

        if idatetick==1
            for i=1:length(xdat)
                xdat(i).x=xdat(i).x;
            end
            xlimits=xlimits;
        end


        if is_met_em==1
            savename=[figname ' ' labs(1).l ' ' filestr ' analysis'];
        else
            savename=[figname ' ' labs(1).l ' ' filestr];
        end

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y' maxx]);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'wrf_timser'

        flt_graph = 'Temp';
        %    flt_graph = 'Height';
        %    flt_graph = 'Vapour';
        %    flt_graph = 'Lat';
        %    flt_graph = 'Lon';
            flt_graph = 'Pressure';
        %    flt_graph = 'Wind';
        %    flt_graph = 'Wind dir';
%        flt_graph = 'Surface fluxes';
        
        mean_at_lat=0; %flag to say we want a timseries of mean values across ice shelf at various latitudes

        %rundir='ant_jan06_sfRUC_v3';

        ih_wrf=0;
        ih_met=0;

        nmark=0;
        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        xlims=1;
        xlimits=[20.1 20.8407]; %some constant level and the descent
        xlimits=[19.2 22.75]; %whole range of flight
        xlimits=[5 8];


        wrf_Tinds = [1:size(Times,1)];
        wrf_Tinds = 2:25;

        %%%%%%%%% setting up of the location of the points for timeseries %%%%%%%%%

        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';

        %    LAT=[-67.5702 -67.1420];
        %    LON=[-68.1297 -61.6650];



        %----   points to plot on the flight track based on the time along the flight track
        times_flight_loc = [19.84565 20.056 20.12 20.2755]; %20.12 is the time where the max was seen for aircraft data - others are just points
        %along the aircraft track
        times_flight_loc = [20.3755]; %20.3755 is the time where the max was seen for aircraft data
        clear it_flt;
        for iflt=1:length(times_flight_loc)
            it_flt(iflt) = findheight(time_flt19,times_flight_loc(iflt));
        end
        LAT_extra = dat_flt19(it_flt,2);
        LON_extra = dat_flt19(it_flt,3);

        if mean_at_lat==1

            np=4; %number of points along the line of longitude to make averages for
            latA=-69;
            latB=-64.8;
            LAT_extra=latA:(latB-latA)/np:latB;
            LAT_extra=-69.5:1:-65.5;
            LON_extra = -60 * ones([1 length(LAT_extra)]);

        end

        

        %LAT_extra =[];
        %LON_extra =[];
        
        LAT_extra =[lat2d.var(125,240)];
        LON_extra =[lon2d.var(125,240)];
        
        LAT_extra2 =[];
        LON_extra2 =[];
        
        

        extra_x=[];
        extra_y=[];

        %----   extra points to plot - give as x and y km
        % extra_x = [575]; %for 12UTC, 6th Jan, ncep polar
        % extra_y = [351];

        % extra_x = [556]; %for 03UTC, 7th Jan, ncep polar
        % extra_y = [355];

        % extra_x = [570]; %for 12UTC, 6th Jan, ecmwf
        % extra_y = [355];

        %  extra_x = [600]; %random
        %  extra_y = [400];

        %  extra_x = [275]; %lefthand side of the equiv cross sections
        %  extra_y = [380];

        %  extra_x = [200 300]; %extra points to the westside of the peninsula
        %  extra_y = [600 500];

%        extra_x = [587 654.48 625 525]; %for 12UTC, 6th Jan, ecmwf_ml_0.5_nudging
%        extra_y = [325 296.87 120 225];
        
%        extra_x = [587 654.48 625 525 675 592]; %for 12UTC, 6th Jan, ecmwf_ml_0.5_nudging
%        extra_y = [325 296.87 120 225 325 222];
        
%        extra_x = [560]; %for 15UTC, 6th Jan, ncep lowest wind level
%        extra_y = [315];
        
%        extra_x = [620]; %for 15UTC, 7th Jan, ncep wind@10m
%        extra_y = [267];

%        extra_x = [654.48]; %Larsen AWS location
%        extra_y = [296.87];

        extra_x = [540 600]; %for 15UTC, 7th Jan, ncep wind@10m
        extra_y = [350 200];



        nlat = size(lat2d.var,1);
        nlon = size(lat2d.var,2);
        dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
        dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

        i_grid = dx_grid * [1:nlon];
        j_grid = dy_grid * [1:nlat];

        for iflt=1:length(extra_x)
            i_extra = findheight(i_grid,extra_x(iflt));
            j_extra = findheight(j_grid,extra_y(iflt));
            LAT_extra2(iflt) = [lat2d.var(j_extra,i_extra)];
            LON_extra2(iflt) = [lon2d.var(j_extra,i_extra)];
        end


        % ascent_str='034'; %'4' is the L-shaped segments - variation of L-shaped segs is minimal
        ascent_str='03'; %0 and 3 are the first and last ascent
        % ascent_str='1'; %upwind side

        if strfind(ascent_str,'0')  %descent after going over the peninsula
            %     LAT=[LAT_extra LAT_extra2];% lat2d(1).var(140,240)];  %places during descent plus other locations
            %     LON=[LON_extra  LON_extra2t];% lon2d(1).var(140,240)];
            LAT=[LAT_extra LAT_extra2];% lat2d(1).var(140,240)];  %places during descent plus other locations
            LON=[LON_extra  LON_extra2];% lon2d(1).var(140,240)];
            as_ds_str='descent';
        elseif strfind(ascent_str,'1')
            LAT=[-67.55 -67.62 -67.55 -66.8 -67.2 LAT_extra2]; %first ascent from Rothera
            LON=[-68.1 -67.8 -67.5 -73.9 -70.9 LON_extra2];
            as_ds_str='ascent';
        elseif strfind(ascent_str,'3')  %final ascent on the way back
            LAT=[LAT_extra LAT_extra2];% lat2d(1).var(140,240)];  %places during descent plus other locations
            LON=[LON_extra LON_extra2];% lon2d(1).var(140,240)];
            as_ds_str='ascent2';
        elseif strfind(ascent_str,'4')  %cumulative ascent during L-shaped flight legs
            LAT=[LAT_extra  ];% lat2d(1).var(140,240)];  %places during descent plus other locations
            LON=[LON_extra  ];% lon2d(1).var(140,240)];
            as_ds_str='ascent2';
        end




        if is_met_em
            cd(['Y:/WRF/' rundir '/AWS_positions/']);
            loadname=['pressure_alt_corr_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'pressure');
            loadname=['temperature_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'temperature');
            loadname=['wind_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'wind');
            if ih_met~=0
                loadname=['vapour_timser_level_' num2str(ih_met) '.mat'];
                load(loadname,'vapour');
            end
            loadname=['winddir_timser_level_' num2str(ih_met) '.mat'];
            load(loadname,'winddir');
            for itemp=1:length(pressure)
                ilat(itemp) = pressure(itemp).ilat;
                ilon(itemp) = pressure(itemp).ilon;
            end
        else
            if mean_at_lat==0
                [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);
            else
                ilat=[];
                ilon=[];
            end
        end

        %%%%% end of setting of location %%%%%


        xlab = 'Time (days)';

        izlim=0;

        multi_var=0;  %set to one if want to plot several variables for one location
        if multi_var==1
            l_plot=1;
        else
            l_plot=length(LAT);
        end

        lor=1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

        %%%%%%%%%%%%%%%%%%%%%% wrf plots %%%%%%%%%%%%%%%%%%%%%%
        for idat=1:l_plot

            if multi_var==1
                iloc=4; %%%%%%%%%%%  set to what location we require %%%%%%%%%%
            else
                iloc=idat;
            end

            switch flt_graph
                case 'Surface fluxes'
                    if is_met_em==1
                        %                    ydat(idat).y = temperature(iloc).y';
                        %                   xdat(idat).x = pressure(1).times;
                    else
                        

                        if mean_at_lat==1

                            flux_type='GRDFLX';
                            flux_type='SH';
                            
                            landmask=nc{'LANDMASK'}(1,:);
                            seaice=nc{'SEAICE'}(1,:); %
                            hgt=nc{'HGT'}(1,:); %


                            icons_inds = get_inds_constant_lat(LAT(iloc),lat2d,lon2d); %get indices for a constant latitude slice

                            ieast=find(lon2d.var(icons_inds)>-67.5); %to remove the mountains to the west of peninsula from being included
                            [peak_height ipeak]=max(hgt(icons_inds(ieast))); %find the position of the peninsula mountain and keep to the east of it
                            lon_peak=lon2d.var(icons_inds(ieast(ipeak)));

                            i2=find( abs(hgt(icons_inds)-0)<1e-5 & abs(landmask(icons_inds)-1)<1e-5 & abs(seaice(icons_inds)-0)<1e-5...
                                & lon2d(1).var(icons_inds)>lon_peak ); % 
                            
                            switch flux_type
                                case 'GRDFLX'
                                    GRDFLX=nc{'GRDFLX'}(wrf_Tinds,:);
                                    GRDFLX=mean(GRDFLX(:,icons_inds(i2)),2);
                                case 'SH'
                                    SH = - nc{'HFX'}(wrf_Tinds,:);
                                    SH = mean(SH(:,icons_inds(i2)),2);
                            end
                            
                            
                            
                        else

                            LW=get_wrf_point_surface(nc,'GLW',wrf_Tinds,ilat(iloc),ilon(iloc));
                            SW=get_wrf_point_surface(nc,'SWDOWN',wrf_Tinds,ilat(iloc),ilon(iloc)); %downwelling SW
                            SH= -( get_wrf_point_surface(nc,'HFX',wrf_Tinds,ilat(iloc),ilon(iloc)) );
                            %SH approximately varies with the difference between TSk and T2 - as would expect - i.e. if air (T2) warmer than surface (TSK) then get flux
                            % of heat from air to surface (positive SH) and vice versa
                            LH= -(get_wrf_point_surface(nc,'LH',wrf_Tinds,ilat(iloc),ilon(iloc)) ); %negative as WRF convention is that these are fluxes into air
                            ALBEDO=get_wrf_point_surface(nc,'ALBEDO',wrf_Tinds,ilat(iloc),ilon(iloc));

%                            ALBEDO=0.78;

                            EMISS=get_wrf_point_surface(nc,'EMISS',wrf_Tinds,ilat(iloc),ilon(iloc));
                            TSK=get_wrf_point_surface(nc,'TSK',wrf_Tinds,ilat(iloc),ilon(iloc));
                            T2=get_wrf_point_surface(nc,'T2',wrf_Tinds,ilat(iloc),ilon(iloc));
                            TSLB=nc{'TSLB'}(wrf_Tinds,:,ilat(iloc),ilon(iloc));
                            GRDFLX=get_wrf_point_surface(nc,'GRDFLX',wrf_Tinds,ilat(iloc),ilon(iloc)); %looks like positive ground flux indicates flux of heat to
                            % the surface from the ground below. E.g. when the ground starts to cool at night due to LW radiative losses it takes heat from the (probably) warmer air and also
                            % from the ground below (which was heated to zero during the daytime) - then the GRDFLX value is positive - indicating that positive means flux
                            % to the surface from below.

                            U10=get_wrf_point_surface(nc,'U10',wrf_Tinds,ilat(iloc),ilon(iloc));
                            V10=get_wrf_point_surface(nc,'V10',wrf_Tinds,ilat(iloc),ilon(iloc));
                            sp10=sqrt(U10.^2+V10.^2);

                            TSLB1=TSLB(:,1);
                            
                                                    %                   ALBEDO=0.78;
                        %                    ALBEDO=0.08;
                        %sign convention is postive means energy going into ground
                        SW_UP=ALBEDO.*SW; %think albedo means this for WRF
                        SW_NET=SW-SW_UP;
                        LW_UP=5.67e-8.*TSK.^4; %Boltzmann W/m2 using skin temperature (need emissivity? - from King paper seems doesn't include)
                        %                   LW_UP=EMISS*5.67e-8.*TSK.^4;
                        LW_NET=LW-LW_UP;


                        MNET=LW_NET + SW_NET + LH+SH;


                        end



                        %



                        nmark=-1; %all markers=-1 --plot markers at all data points


                        if multi_var==1
                            

                            %                        ydat(idat).y = SW_NET + LW_NET;
                            %                        ydat(idat+1).y = -(LH+SH);
                            %                        labs(idat).l = ['Net SW+LW'];
                            %                        labs(idat+1).l = ['Net SH+LH (-ve)'];

                            %                        labs(idat).l = ['SH'];
                            %                        labs(idat+1).l = ['GRDFLX'];

%                            labs(idat).l = ['SH'];
%                            labs(idat+1).l = ['LH'];
                            
                           % labs(idat).l = ['MNET+GRDFLX'];
%                            labs(idat+1).l = ['SH+LH'];

                            %                        labs(idat).l = ['TSLB1'];
                            %                        labs(idat+1).l = ['TSK'];
%                                                   labs(idat+2).l = ['T2'];

                            extras='';
                            %                        extras='-273.15'; %extra statements on the end (e.g. an offset)


                            for imulti=1:length(labs)
                                eval(['ydat(imulti).y = ' labs(imulti).l extras ';']);
                                    
                                for it=1:length(wrf_Tinds)
                                    iit=wrf_Tinds(it);
                                    xdat(imulti).x(it) = str2num(Times(iit,9:10)) + str2num(Times(iit,12:13))/24;
                                end
                            end
                            
                            


                            ylab='Temperature (^{o}C)';
%                            ylab='Flux (W m^{-2})' ;

                            %                        xdat(idat+1).x = xdat(idat).x';
                        else

                            for it=1:length(wrf_Tinds)
                                iit=wrf_Tinds(it);
                                xdat(idat).x(it) = str2num(Times(iit,9:10)) + str2num(Times(iit,12:13))/24;
                            end

                            %                        ydat(idat).y = MNET;
                            %                        ydat(idat).y = T2-273.15;
                                 ydat(idat).y = GRDFLX;
                               ydat(idat).y = SW_NET; 
%                                 ydat(idat).y = SH;
%                                  ydat(idat).y = LW_NET+SW_NET;
%                               ydat(idat).y = T2-273.15;
%                                ydat(idat).y = TSK-273.15;
                                ydat(idat).y = GRDFLX+MNET;  %overall melt flux rate
                            %     ydat(idat).y = LH+SH;
                            
%                                  ydat(idat).y = 
                            %                        ydat(idat).y = sp10;


                            labs(idat).l = filestr;

                            %ylab='LH+SH Flux (W m^{-2})';
                            ylab='Ground Flux (W m^{-2})';
                            ylab='Sensible Heat Flux (W m^{-2})';  
                            ylab='Melt heat flux (W m^{-2})';  
%                                            ylab='Net Flux (W m^{-2})';
%                                            ylab='Temperature at 2m (^{o}C)';
%                                        ylab='Skin temperature (^{o}C)';
                            %                ylab='GRDFLX + SH (W m^{-2})';
                            %                ylab='10m wind speed (m s^{-1})';
                           % ylab='Net radiation (W m^{-2})';
%                            ylab='LW (W m^{-2})';
%                            ylab='Net SW flux (W m^{-2})';
                           %  ylab='LW up radiation (W m^{-2})';


                        end


                    end

                    xdat(idat).x = xdat(idat).x';





                    if wrf_Tinds(1)==1;
                        if multi_var==1
                            for imulti=1:length(ydat)
                                xdat(imulti).x = xdat(imulti).x(2:end);  %as the fluxes for the first time are all zero
                                ydat(imulti).y = ydat(imulti).y(2:end);
                            end
                        else
                            xdat(idat).x = xdat(idat).x(2:end);  %as the fluxes for the first time are all zero
                            ydat(idat).y = ydat(idat).y(2:end);
                        end
                    end



                case 'Pressure'
                    if is_met_em
                        ydat(idat).y = pressure(iloc).y'/100;
                        xdat(idat).x = pressure(1).times;
                    else
                        ydat(idat).y = (nc{'P'}(wrf_Tinds,1,ilat(iloc),ilon(iloc)) + nc{'PB'}(wrf_Tinds,1,ilat(iloc),ilon(iloc)) )/100;
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
                    end

                    xdat(idat).x = xdat(idat).x';

                case 'Temp'
                    if is_met_em==1
                        ydat(idat).y = temperature(iloc).y';
                        xdat(idat).x = pressure(1).times;
                    else
                        potemp = nc{'T'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + 300;
                        P = nc{'P'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'PB'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc));
                        ydat(idat).y = potemp ./ ( (1e5./P).^0.286 ) - 273.15;
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
                        labs(idat).l = filestr;
                    end

                    xdat(idat).x = xdat(idat).x';

                    ylab='Temperature (^{o}C)';


                case 'vapour'
                    if is_met_em
                        rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                        T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                        P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                        qsat = satvappress(T,'goff','liq',P,1)/f;
                        xdat(idat).x = 1000 * rh/100 .* qsat;
                    else
                        xdat(idat).x = 1000*nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                    end

                    ylab


                case 'Wind'
                    if is_met_em
                        ydat(idat).y = wind(iloc).y';
                        xdat(idat).x = pressure(1).times;
                    else
                        u=0.5* ( nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc)+1,ilon(iloc)) );
                        v=0.5* ( nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)+1) );
                        ydat(idat).y= sqrt( u.^2 + v.^2 );
                        for it=wrf_Tinds
                            xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                        end
                    end

                    xdat(idat).x = xdat(idat).x';

                    ylab='Wind speed (m s^{-1})';

                    lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane



                case 'Wind dir'
                    for it=wrf_Tinds
                        xdat(idat).x(it) = str2num(Times(it,9:10)) + str2num(Times(it,12:13))/24;
                    end
                    xdat(idat).x = xdat(idat).x';
                    labs(idat).l = filestr;

                    if is_met_em
                        u=nc{'UU'}(ih_wrf,:,ilat(iloc),ilon(iloc));
                        v=nc{'VV'}(ih_wrf,:,ilat(iloc),ilon(iloc));
                    else
                        u=0.5* ( nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'U'}(wrf_Tinds,ih_wrf,ilat(iloc)+1,ilon(iloc)) );
                        v=0.5* ( nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)) + nc{'V'}(wrf_Tinds,ih_wrf,ilat(iloc),ilon(iloc)+1) );
                    end

                    jnorth = ilat(iloc) + 10;
                    lons_north = lon2d.var(jnorth,:);
                    [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

                    %angle of the local north line relative to the grid
                    thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );




                    for iuv=1:length(u)

                        theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                        if u(iuv)==0 & v(iuv)==0
                            ydat(idat).y(iuv) = 0;
                        elseif u(iuv)>=0 & v(iuv)>=0
                            ydat(idat).y(iuv) = theta2;
                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                            ydat(idat).y(iuv) = 180 + theta2;
                        elseif u(iuv)<=0 & v(iuv)<=0
                            ydat(idat).y(iuv) = 180 + theta2;
                        elseif u(iuv)<0 & v(iuv)>0
                            ydat(idat).y(iuv) = 360 + theta2; %theta2 is negative
                        end

                    end

                    ydat(idat).y = ydat(idat).y' + 180 - thetaN*180/pi; %add 180 to make it the direction wind is coming from
                    % take away thetaN to give direction relative to north
                    i360 = find(ydat(idat).y>=360);
                    ydat(idat).y(i360) = ydat(idat).y(i360) - 360;

                    izlim=1;
                    zmin=0;
                    zmax=450;




                case 'cloud'
                    cloud=nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                    cloud=cloud+nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                    xdat(i).x=1000*cloud;


                    figname=['Cloud mixing ratio profile at ' tstr ' for ' filestr];
                    xlab='Total condensed water mixing ratio (g kg^{-1})';


            end
            
            if multi_var==0
                iloc=idat; %otherwise is set earlier to that wanted
            end

            abc=['ABCDEFGHIJKLM'];

            if (length(strfind(ascent_str,'1'))>0)
                loc_lab=['W' abc(iloc)];
            elseif (length(strfind(ascent_str,'0'))>0 | length(strfind(ascent_str,'3'))>0 | length(strfind(ascent_str,'4') )>0)
                loc_lab=['E' abc(iloc)];
            end

            if multi_var==0
                if mean_at_lat==0
                    labs(idat).l=[num2str(lat2d.var(ilat(iloc),ilon(iloc)),3) ' , ' num2str(lon2d.var(ilat(iloc),ilon(iloc)),3) ];
                else
                    labs(idat).l=[loc_lab ' ' num2str(LAT(iloc),3)];
                end
            end


            %             if length(strfind(ascent_str,'1'))>0 & multi_var==0
            %                 labs(idat).l=['W' abc(idat) ' ' num2str(lat2d.var(ilat(iloc),ilon(iloc)),3) ' , ' num2str(lon2d.var(ilat(iloc),ilon(iloc)),3) ];
            %             elseif (length(strfind(ascent_str,'0'))>0 | length(strfind(ascent_str,'3'))>0 | length(strfind(ascent_str,'4') )>0) & multi_var==0
            %                 labs(idat).l=['E' abc(idat) ' ' num2str(lat2d.var(ilat(iloc),ilon(iloc)),3) ' , ' num2str(lon2d.var(ilat(iloc),ilon(iloc)),3) ];
            %             end

        end



        if is_met_em
            titlenam=[ylab ' for ' num2str(pressure(1).y(1)/100) ' hPa for ' strrep(rundir,'_',' ') ' analysis'];if ~exist('subplotting'); subplotting=0; end
        else
            if multi_var==0
                titlenam=[ylab ' for level ' num2str(ih_wrf) ' for ' filestr];if ~exist('subplotting'); subplotting=0; end
            else
                titlenam=[ylab ' for level ' num2str(ih_wrf) ' for ' filestr ' for ' loc_lab];if ~exist('subplotting'); subplotting=0; end
            end
        end
        figname=titlenam;
        savename=[titlenam];


        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y' maxx]);
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'antjan06_flt'

        %profiles from aircraft data
        %descent starts at ~20.8 UTC and finishes at ~21 UTC

        %[ibeg iend]=findheight(dat(:,1)/3.6e6,20.8,21.0);
        %inds_prof=ibeg:iend;
        %plot(qv2(inds_prof),dat(inds_prof,11));


        % 1) time - milliseconds (I think - divide by 3.6e6 to get decimal hours)
        % 2) lat
        % 3) lon
        % 4) static pressure  (mb)
        % 5) surface temp (from infra red thermometer)
        % 6) air temperature
        % 7) dew point (from frost point hygrometer)  (degC) *** wasn’t working properly according to the paper and they used the humicap data ***
        % 8) dew point (from humicap)		           (degC)
        % 9) wind speed
        % 10) wind direction (need to add 180 to get conventional wind direction - sorry)
        % 11) altitude (gps - not very good)
        % 12) short wave - upwelling
        % 13) short wave - downwelling
        % 14)long wave - upwelling
        % 15) long wave - down welling
        

        %flight_no should be automatically picked out in the reading of the CAS file
% flight_no='19';


iaxis_square=0; %switch to make axis square
fsize=18;
iplot_error_bar=0;


eval(['dat_flt = dat_flt' flight_no ';']);  %put the data for the required flight here


         if ~exist('man_choose_flt_graph')                  
            flt_graph = 'Temp';
%            flt_graph = 'Height';
%            flt_graph = 'Vapour'
%            flt_graph = 'Lat';
%            flt_graph = 'Lon';
%            flt_graph = 'Pressure';
%            flt_graph = 'Wind';
            flt_graph = 'Wind dir';
%            flt_graph = 'SW up';   
%            flt_graph = 'SW down';   
%            flt_graph = 'LW up';   
%            flt_graph = 'LW down';   
%            flt_graph = 'SW net down';               
%            flt_graph = 'LW net down';               
%            flt_graph = 'SW+LW net down';                           
%            flt_graph = 'Altitude';
            flt_graph = 'Radar Altitude';
%            flt_graph = 'Altitude3D';    
%            flt_graph = 'Temperature3D';                
%            flt_graph = 'Potemp';
%             flt_graph = 'Frost point humicap'
%              flt_graph = 'Frost point hygrometer'
%              flt_graph = 'RHi hygrometer'
%            flt_graph = 'Airspeed CAS';
%            flt_graph = 'Airspeed';
%            flt_graph = 'Ice number';
%            flt_graph = 'Ice mass';
%            flt_graph = 'Mean ice diameter';
%             flt_graph = 'Estimated displacement';
             
              iplot_3D=0; 
              iaxis_square=1;
              
              ismooth=[1 1 1 1 1 1];
              ismooth=[0 0 0 0 0 0];
              nfilter=40;
              nfilter=150; %smoothing filter no. points
              
              ihighlight_cloud=1; %flag to say whether want to highlight cloud in flight plots
                
         else
             clear man_choose_flt_graph
         end
         

%%%         %set the required ix_distance value below if are not using set_ix_distance  %%%
         if ~exist('set_ix_distance')
             switch flt_graph
                 case {'Altitude3D','Temperature3D'}
                     ix_distance=1;  %%set whether want distance from a point or time as the x-axis %%%%
                 otherwise
                     ix_distance=0;   %%default value
                     %ix_distance=1;   %%set the default value
             end
         else
             clear set_ix_distance
         end



switch flight_no
    case '99'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
%        CameraPosition = [350-500 250-500 31241.7];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngle = [10.1035];        
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;
        
    case '100'
        X_pos = 350;
        Y_pos = 350; %position to calculate distance from in km

        X_pos = 285;
        Y_pos = 450; %NW part of flight transect for flight 100
        
        
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
%        CameraPosition = [350-500 250-500 31241.7];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngle = [10.1035];        
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;
        
        

%     case '117'
%         X_pos = 365;
%         Y_pos = 382; %Rothera
%         
%         X_pos = 340;
%         Y_pos = 345; %Position at the southern extremity of the flight path
        
    case '104'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;
        
     case '105'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;
        
       case '108'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;  
        
    case '113'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;      
        
 case '117'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;  
        
    case '120'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;      
        
    case '122'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;   
        
case '123'
        X_pos = 365;
        Y_pos = 382; %Rothera
        
        X_pos=0;
        Y_pos=0;
        
        %settings for the best 3D view
        View = [25.7034 87.4156];
        CameraPosition = [622.66+365 -1307.89+382 31241.7];
        CameraPositionMode = 'manual';
        CameraTarget = [42.7639+365 -103.139+382 1619.62];
        CameraTargetMode = 'manual';
        CameraUpVector = [0 0 1];
        CameraUpVectorMode = 'manual';
        CameraViewAngle = [6.71035];
        CameraViewAngleMode = 'manual';
        
        iset_3D_view_properties=1;           
        
        
    case '102'
        X_pos = 365;
        Y_pos = 382; %Rothera 
        
        X_pos=0;
        Y_pos=0;
        
        iset_3D_view_properties=1;
        
        %settings for the best 3D view
        view_3d='cloud top (altitude)';    
%        view_3d='cloud top 20:20-21:00 (altitude)'; 
%         view_3d='cloud top (altitude) WRF grid';
%        view_3d='cloud top (temperature)';
        view_3d='top down';
%        view_3d='side on (altitude)';
        
        switch view_3d
            case 'cloud top (altitude)'
            %close up view of the cloud top
            View = [40.3771 87.3751];
        	CameraPosition = [393.455 -605.214 18284.1];
            CameraPositionMode = 'manual';
            CameraTarget = [-83.7947 -43.9931 2214.55];
            CameraTargetMode = 'manual';
            CameraUpVector = [0 0 1];
            CameraUpVectorMode = 'manual';
            CameraViewAngle = [3.40713];         
            CameraViewAngleMode = 'manual';
            
        X_pos = 365;
        Y_pos = 382; %Rothera 
        
        case 'cloud top (altitude) WRF grid'
            %close up view of the cloud top
            View = [40.3771 87.3751];
        	CameraPosition = [393.455+365 -605.214+382 18284.1];
            CameraPositionMode = 'manual';
            CameraTarget = [-83.7947+365 -43.9931+382 2214.55];
            CameraTargetMode = 'manual';
            CameraUpVector = [0 0 1];
            CameraUpVectorMode = 'manual';
            CameraViewAngle = [9.40713];         
            CameraViewAngleMode = 'manual';
            
            X_pos = 365;
            Y_pos = 382; %Rothera 
            
            X_pos = 0;
            Y_pos = 0; %

        
       
            
        
        
            
            case 'cloud top 20:20-21:00 (altitude)'
                CameraPosition = [146.274 -346.849 20392.9]
                CameraPositionMode = 'manual'
                CameraTarget = [-54.088 -36.9182 2845.07]
                CameraTargetMode = 'manual'
                CameraUpVector = [0 0 1]
                CameraUpVectorMode = 'manual'
                CameraViewAngle = [3.40713]
                CameraViewAngleMode = 'manual'
                
                
            X_pos = 365;
            Y_pos = 382; %Rothera             
        
            case 'side on (altitude)'
                CameraPosition = [400 200 2845.07];
%                CameraPosition = [350 300 2845.07];                
                CameraPositionMode = 'manual'
                CameraTarget = [325 325  2845.07]
                CameraTargetMode = 'manual'
                CameraUpVector = [0 0 1]
                CameraUpVectorMode = 'manual'
                CameraViewAngle = [89.8794];
                CameraViewAngleMode = 'manual'
                
                
%            X_pos = 365;
%            Y_pos = 382; %Rothera             
        
            
            
            case 'cloud top (temperature)'
            %close up view of the cloud top
        iset_3D_view_properties=1;            
            View = [26.9381 -2.25919];
%            View = [22.9381 -2.25919];
            
        	CameraPosition = [378.882 -906.53 -43.7944]
            CameraPositionMode = 'manual';
            CameraTarget = [-58.7939 -45.2448 -5.68062]
            CameraTargetMode = 'manual';
            CameraUpVector = [0 0 -1]
            CameraUpVectorMode = 'manual';
            CameraViewAngle = [7.11849]
            CameraViewAngleMode = 'manual';
            
            case 'top down'
            %top down view
            View = [-0.382523 90];
            CameraPosition = [-42.215 -32.3087 30639.7];
            CameraPosition = [300 300 30639.7];            
            CameraPositionMode = 'auto';
            CameraTarget = [-42.215 -32.3087 1694.46];
            CameraTarget = [300 300 1694.46];            
            CameraTargetMode = 'auto';
            CameraUpVector = [0 1 0];
            CameraUpVectorMode = 'auto';
            CameraViewAngle = [7.24041];
            CameraViewAngleMode = 'manual';
            
        end
end





if ix_distance==1
    eval(['X_flt = X_flt' flight_no ';']);
    eval(['Y_flt = Y_flt' flight_no ';']);
    time_flt = sqrt( (X_flt-X_pos).^2 + (Y_flt-Y_pos).^2 ); %distance from position (X_pos,Y_pos)
    xlab=['Distance (km)'];
    idatetick=0; %flag to say the want the xaxis in proper time format rather than decimal time    

elseif ix_distance==2
    eval(['X_flt = X_flt' flight_no ';']);
    eval(['Y_flt = Y_flt' flight_no ';']);
    time_flt = [0; cumsum(sqrt( (diff(X_flt)).^2 + (diff(Y_flt)).^2 ))]; %distance from position (X_pos,Y_pos)
    xlab=['Distance along flight track (km)'];
    idatetick=0; %flag to say the want the xaxis in proper time format rather than decimal time    

else
    eval(['time_flt = time_flt' flight_no ';']);
    xlab=['Time (UTC)'];
    idatetick=1; %flag to say the want the xaxis in proper time format rather than decimal time
    datetick_type=15; %specify the type with datetick_type (see help datetick) 15= HH:MM
end

eval(['time_flt2 = time_flt' flight_no ';']); %the actual time again





% if idatetick==1;
%     time_flt = time_flt/24; %convert to days
% end




        nmark=0;
        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane


if ~exist('iset_xlimits') | iset_xlimits==0
    
        xlims=1;
%        xlimits=[20.9 21.5]; %spike in wind speed after 2nd L-shaped leg
%        xlimits=[20.1 20.8407]; %some constant level and the descent
%        xlimits=[19.2 22.75]; %whole range of flight
        %         xlimits=[19.2 20.5]; %ascent and descent
        
%        xlimits=[20.3972 20.7435]; %1st constant level L-shaped leg - approx. Range (uncorrected from aircraft data = 67-89m). 
                                %Quoted as 15m abv ice shelf surface in King paper        
%        xlimits=[20.76 21.176]; %2nd constant level L-shaped leg - approx. Range (uncorrected from aircraft data = 197-272m). 
                                %Quoted as 152m abv ice shelf surface in King paper
%        xlimits=[21.1912 21.55]; %3rd constant level L-shaped leg - approx. Range (uncorrected from aircraft data = 346-439m)
                                %Quoted as 305m abv ice shelf surface in King paper
%        xlimits=[21.577 21.959]; %4th constant level L-shaped leg - approx. Range (uncorrected from aircraft data = 652-728m)
        %Quoted as 610m abv ice shelf surface in King paper
        
%        xlimits=[19.4817 20.2664] %approx the constant level flight at 3000 m across the mountain


%using max and min here as the x axis is not necessarily always increasing
%e.g. as when x is distance

        xlimits=[min(time_flt2) max(time_flt2)]; %whole range of the flight
        dxlims=xlimits(2)-xlimits(1);
        percent_either_side=5;
        xlimits(1)=xlimits(1)-dxlims*percent_either_side/100;
        xlimits(2)=xlimits(2)+dxlims*percent_either_side/100;
        
else
%xlimits are set in plot_flight_segments script
    xlims=1;
    clear iset_xlimits
end




[inds_1 inds_2] = findheight_nearest(time_flt2,xlimits(1),xlimits(2));
inds=inds_1:inds_2;




ilegs=0; %flag to say whether to plot all the legs at the different heights on 
% one plot. But, is split into the N-S transect and the E-W one
% ilegs=1 for N-S and ilegs=2 for E-W. ilegs=0 for normal plot.

            




%%% these are set in plotTimeHeightVap3.m - might want to override here
            % or to make sure they are the same if requried
            if ~exist('cut_off_size')
                cut_off_size=1; %size (microns) below which to ignore counts for particle concentration
            end
            if ~exist('air_speed_type')
                %air_speed_type = 'aircraft';
                %air_speed_type = 'constant 60m/s';
                air_speed_type = 'CIP probe';
            end
            
            air_speed_type = 'constant';
            air_speed_type = 'CIP probe';
    %%%    ------------------------------------   %%%
    
            CAS_LWC_cut_off_sizes = [10 30]; %size outside of which to ignore LWC from CAS
            
            
            switch flight_no
                case {'19','09'}
                    %do nothing
                otherwise
                    if exist('CAS_time_all')
                        time_timeseries=CAS_time_all;
                        data_particle=CAS_counts_all;
                        
           [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number...
                    ,CAS_total_number_cutoff ...                    
                    ,CIP_total_number,LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                    =cas_sample_volume_and_stats2...
                    (dat_flt,time_timeseries,...
                   CAS_bins,data_particle,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all...
                    ,CAS_LWC_cut_off_sizes,0);
                
%            %get the sample volume and total concentrations, plus air speed if required.           
%             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
%                 ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip...
%                 ,CAS_mode_diameter,CAS_mean_diameter,LWC_dist_cas_cutoff]...
%                 =cas_sample_volume_and_stats(dat_flt,CAS_time_all,...
%                 CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all...
%                 ,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes);   
                    end
            
            end
            
            


if size(dat_flt,2)==15
    %for flt_19
    col_temp=6;
    col_alt=11;
    col_lat=2;
    col_lon=3;
    col_press=4;
    col_wind=9;
    col_winddir=10;
elseif size(dat_flt,2)==49
    mpace_column_numbers
else
    %for Feb2010 flights
    col_temp=5;
    col_alt=12;
    col_radalt=11;
    col_lat=2;
    col_lon=3;
    col_press=6;
    col_wind=9;
    col_winddir=10;    
    col_frostpoint_hygro=7;    
    col_frostpoint_humi=8;
    col_airspeed=4;
end




        switch flt_graph
            case 'Number aerosol 610-1030 nm'
                %calculte the number of particles within a size limit        
        
                Nlim_sizes=[0.6 1.04]; %calculate the number in between these limits  
%                Nlim_sizes=[0.6 10.0]; %calculate the number in between these limits                  
                Nlim_sizes=[0.6 3.04]; %calculate the number in between these limits  
                Nlim_sizes=[0.6 50]; %calculate the number in between these limits  
                
                if Nlim_sizes(1)==0
                    ilims_lower=1;
                else
                    ilims01 = find(CAS_bins>=Nlim_sizes(1));
                    ilims_lower = ilims01(1)+1; %add one as the first number bin is for D<0.61
                end
                ilims02 = find(CAS_bins<=Nlim_sizes(2));                
                ilims_upper = ilims02(end); %bin N+1 is for sizeN to sizeN+1. Want N for 
                                            %sizes<sizeN so is bin N that we require
                N_scaled = CAS_counts_all./sample_volume_CAS;    
                

                                            
                ydat(1).y = interp1(CAS_time_all,sum(N_scaled(:,ilims_lower:ilims_upper),2),dat_flt(:,1)/1e3);
                ylab=['N_{' num2str(Nlim_sizes(1)) '-' num2str(Nlim_sizes(2)) ' \mum} (cm^{-3})'];       
                      
                if adjust_STP==1
                    disp('**** WARNING - concentrations are adjusted to STP values *****');
                    rho_factor_aircraft = rho_aircraft / rho_stp;  %=Vol_stp / Vol_altitude 
                                                                   % for constant mass
                    ydat(1).y = ydat(1).y./rho_factor_aircraft;     
                    
                    ylab = [ylab ' at STP'];
                end
                
                      
                xdat(1).x = time_flt;
                                
%                iydir=-1;
                
                labs(1).l='N';
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                izlim=0;
                zmin=-13;
                zmax=-8;
                
                ismooth(1)=1;
                nfilter=8;


            case 'Estimated displacement'
               %is an LWC point at 4 g/m3 that must be wrong                                
                LWC = interp1(CIP_time_all,LWC_CAS_all',dat_flt(:,1)/1e3);
                i_ignore = find(LWC<-0.05 & LWC>2.5); 
                LWC(i_ignore)=NaN;
                iLWC=find(LWC>=0.08);
                
                T=273.15+dat_flt(:,col_temp);
                P=100*dat_flt(:,col_press);
% for qv will use satuartion qv when have LWC and vapour measurement when do no
% since the qv measurement is likely to be unreliable in cloud
                clear qv
                qv=qv_flt_humi;
                qv(iLWC) = SatVapPress(T(iLWC),'goff','liq',P(iLWC),1)/f;                
                equiv = equivalent_potemp(T,P,qv);
                heights=dat_flt(:,col_alt);
                               
                w_obs = w2_turb;
                w_obs(i_ignore)=NaN;
                
                %idealised soundings from wave_temperature... Matlab model
                eq_use = eq_h_fp;
                eq_use = eq_h_humi; %the above two are very similar
                h_use = h_wave_model;
            %actual sounding
%                eq_use = equiv_sound_humi;
%                h_use = Zsound;
                
                min_LWC=-0.05;
                
                for ieq=1:length(equiv)
                    if equiv(ieq)>=min(eq_use) & equiv(ieq)<=max(eq_use) & LWC(ieq)>min_LWC
%                        isource = findheight(eq_h_fp,equiv(i));
%                        xdat(1).x(i) = heights(i) - interp1(eq_use,h_use,equiv(i));
                        ydat(1).y(ieq) = heights(ieq) - interp1(eq_use,h_use,equiv(ieq));                        
                    elseif equiv(ieq)<min(eq_use)
                        ydat(1).y(ieq)=800;
                    elseif equiv(ieq)>max(eq_use)
                        ydat(1).y(ieq)=-800;  %came from higher up than sounding
                    else
                        ydat(1).y(ieq)=NaN;
                    end
                end
                
                
                ylab=['Displacement (m)'];             
                xdat(1).x = time_flt;
                                
%                iydir=-1;
                
                labs(1).l='Displacement';
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                izlim=0;
                zmin=-13;
                zmax=-8;
                
                
            case 'Temp'

                ylab=['Temperature (^{o}C)'];
                ydat(1).y = dat_flt(:,col_temp);                
                xdat(1).x = time_flt;
                                
                iydir=-1;
                
                labs(1).l='Temperature';
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                izlim=0;
                zmin=-13;
                zmax=-8;
                
            case 'Vertical wind speed'

                ylab=['Vertical wind speed (m s^{-1})'];
                ydat(1).y = w2_turb;                
                xdat(1).x = time_flt;
                                
%                iydir=-1;
                
                labs(1).l='W';
                
                lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                izlim=1;
                zmin=-6;
                zmax=10;
                
                

            case 'Ice number'

                ylab=['Ice number (L^{-1})'];
%                ydat(1).y = interp1(CIP_time_Jonny2/3600,1000*ice_no_Jonny,time_flt);                
%                xdat(1).x = time_flt;

                xdat(1).x=CIP_time_Jonny2/3600;
                ydat(1).y=1000*ice_no_Jonny;

                labs(1).l='Ice number';
                
                lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                

                
                ihighlight_cloud=0;
                
             case 'Ice mass'

                ylab=['Ice mass (mg m^{-3})'];

                xdat(1).x=CIP_time_Jonny2/3600;
                ydat(1).y=1000*ice_mass_Jonny;

                labs(1).l='Ice mass';
                
                lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                

                
                ihighlight_cloud=0; 
                
                izlim=0;
                zmin=0;
                zmax=0.5;

                
            case 'Mean ice diameter'

                ylab=['Diameter (\mum)'];
%                ydat(1).y = interp1(CIP_time_Jonny2/3600,1000*ice_no_Jonny,time_flt);                
%                xdat(1).x = time_flt;

                xdat(1).x=CIP_time_Jonny/3600;
                ydat(1).y=mean_ice_size;

                labs(1).l='Mean ice diameter';
                
                lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                

                
                ihighlight_cloud=0;    
                
                
            case 'Height'

                ylab=['Height (km)'];
                ydat(1).y = dat_flt(:,col_alt);
                xdat(1).x = time_flt;
                labs(1).l='Height';
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

            case 'Vapour'

                ylab=['Water vapour mixing ratio (g kg^{-1})'];
                ydat(1).y = qv2*1000;
                xdat(1).x = time_flt;
                labs(1).l='Vapour';

            case 'Lat'

                ylab=['Latitude (degrees)'];
                ydat(1).y = dat_flt(:,col_lat);
                xdat(1).x = time_flt;
                labs(1).l='Latitude';

            case 'Lon'

                ylab=['Longitude (degrees)'];
                ydat(1).y = dat_flt(:,col_lon);
                xdat(1).x = time_flt;
                labs(1).l='Longitude';

            case 'Pressure'

                ylab=['Pressure (mb)'];
                ydat(1).y = dat_flt(:,col_press);
                xdat(1).x = time_flt;
                labs(1).l='Pressure';
                
                lor=-1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                

            case 'Vertical speed of aircraft'

                ylab=['Aircraft vertical speed (m s^{-1})'];
                alt=dat_flt(:,col_alt);
                %rate of aircraft ascent
                ydat(1).y = [0; diff(dat_flt(:,col_alt))./ diff(time_flt2*3600) ];
                xdat(1).x = time_flt;
                labs(1).l='Ascent speed';
                
                izlim=0;
                zmin=0;
                zmax=20;
                
            case 'Wind'

                ylab=['Wind speed (m s^{-1})'];
                ydat(1).y = dat_flt(:,col_wind);
                xdat(1).x = time_flt;
                labs(1).l='Wind speed';
                
                izlim=1;
                zmin=0;
                zmax=30;
                
                iNaN = find(ydat(1).y>25);
                ydat(1).y(iNaN)=NaN;

                
            case 'Wind dir'

                ylab=['Wind direction (degrees)'];
                ydat(1).y = dat_flt(:,col_winddir)+180; %need to add 180 according to TLC
                ydat(1).y(ydat(1).y>360)=ydat(1).y(ydat(1).y>360)-360;
                xdat(1).x = time_flt;
                labs(1).l='Wind direction';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=1;
                zmin=0;
                zmax=360;
                
            case 'SW up'

                ylab=['SW up (W m^{-2})'];
                ydat(1).y = dat_flt(:,12);               
                xdat(1).x = time_flt;
                labs(1).l='SW up';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=1;
                zmin=260;
                zmax=460;
                
            case 'SW down'

                ylab=['SW down (W m^{-2})'];
                ydat(1).y = dat_flt(:,13);               
                xdat(1).x = time_flt;
                labs(1).l='SW down';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=1;
                zmin=300;
                zmax=800;

            case 'LW up'

                ylab=['LW up (W m^{-2})'];
                ydat(1).y = dat_flt(:,14);               
                xdat(1).x = time_flt;
                labs(1).l='LW up';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=1;
                zmin=290;
                zmax=350;

                
            case 'LW down'

                ylab=['LW down (W m^{-2})'];
                ydat(1).y = dat_flt(:,15);               
                xdat(1).x = time_flt;
                labs(1).l='LW down';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=1;
                zmin=200;
                zmax=260;
                
              case 'SW net down'

                ylab=['SW net down (W m^{-2})'];
                ydat(1).y = dat_flt(:,13)-dat_flt(:,12);                  
                xdat(1).x = time_flt;
                labs(1).l='SW net down';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=0;
                zmin=200;
                zmax=260;   
                
             case 'LW net down'

                ylab=['LW net down (W m^{-2})'];
                ydat(1).y = dat_flt(:,15)-dat_flt(:,14);                
                xdat(1).x = time_flt;
                labs(1).l='LW net down';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=0;
                zmin=200;
                zmax=260;
                
             case 'SW+LW net down'

                ylab=['SW+LW net down(W m^{-2})'];
                ydat(1).y = dat_flt(:,13)-dat_flt(:,12)+dat_flt(:,15)-dat_flt(:,14);               
                xdat(1).x = time_flt;
                labs(1).l='SW+LW net down';
                
                lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
                izlim=1;
                zmin=0;
                zmax=300;
                
                
             case 'Altitude'

                ylab=['Altitude (m)'];
                ydat(1).y = dat_flt(:,col_alt);
                xdat(1).x = time_flt;
                labs(1).l='Altitude';
               
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane    
                
            case 'Radar Altitude'

                ylab=['Radar Altitude (m)'];
                ydat(1).y = dat_flt(:,col_radalt);
                xdat(1).x = time_flt;
                labs(1).l='Altitude';
               
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane    
   
                
             case 'Aircraft T pert'

                ylab=['Estimated T pert (K)'];
                alt = dat_flt(:,col_alt);
                press = dat_flt(:,col_press); %mb
                T = dat_flt(:,col_temp)+273.15; %K
                pot = T.*(1000./press).^0.286;
                pot_mean = mean(pot(inds));
                press_mean = mean(press(inds));
                %the exact values used for these baselines shouldn't matter too much
                %as the relative changes in potemp for a 1 or two degree change are quite small
                %e.g. 3/293 = 0.01, i.e. 1% error. similarly for pressure
                
                ydat(1).y = pot_mean * ( (press/1000).^0.286 - (press_mean/1000).^0.286 );                
                
                xdat(1).x = time_flt;
                

                
                
                labs(1).l='Aircraft T pert';
                
                iydir=-1;
                               
                lor=3;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane       
                
             case 'Altitude3D'

                ylab=['Y (km)']; 
                yylab='Altitude (m)';
                xlab=['X (km)'];

                labs(1).l='Altitude';
                
                lor=2;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                xdat(1).x = X_flt(:)-X_pos;
                ydat(1).y = dat_flt(:,col_alt);
                yydat(1).y = Y_flt(:)-Y_pos;

                
                iplot_3D=1;    
                
                xlims=0;
                izlim=0;
                
                
            case 'Temperature3D'

                ylab=['Y (km)']; 
                yylab='Temperature (^{o}C)';
                xlab=['X (km)'];

                labs(1).l='Temperature';
                
                lor=3;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                xdat(1).x = X_flt-X_pos;
                ydat(1).y = dat_flt(:,col_temp);
                yydat(1).y = Y_flt-Y_pos;

                
                iplot_3D=1;    
                
                xlims=0;
                izlim=0;  
                
                iydir=-1;
                
            case 'Potemp'

                ylab=['Potential temperature (K)'];
                T=dat_flt(:,col_temp)+273.15;
                P=dat_flt(:,col_press)*100;
                
                ydat(1).y = T.*(1000e2./P).^0.286; 
                xdat(1).x = time_flt;
                labs(1).l='Potential temperature';  
                
                lor=3;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'Equivalent Potemp'

                ylab=['Equivalent Potential temperature (K)'];
                T=dat_flt(:,col_temp)+273.15;
                P=dat_flt(:,col_press)*100;
                qv=qv_flt_humi;
%                qv(iLWC) = SatVapPress(T(iLWC),'goff','liq',P(iLWC),1)/f;                

                
                ydat(1).y = equivalent_potemp(T,P,qv); 
                xdat(1).x = time_flt;
                labs(1).l='Equivalent Potential temperature';  
                
                lor=3;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                

                
                
             case 'Frost point humicap'

                ylab=['Frost point humicap (^{o}C)'];
                
                ydat(1).y = dat_flt(:,col_frostpoint_humi); 
                xdat(1).x = time_flt;
                labs(1).l='Frost point humicap'; 
                
            case 'Frost point hygrometer'

                ylab=['Frost point hygrometer (^{o}C)'];
                
                ydat(1).y = dat_flt(:,col_frostpoint_hygro);
                xdat(1).x = time_flt;
                labs(1).l='Frost point hygrometer'; 
                
            case 'RHi hygrometer'

                ylab=['RHi hygrometer (%)'];
                eval(['qv_dat = qv_flt' flight_no '_fp;'])
                qv_sat = satvappress(dat_flt(:,col_temp)+273.15,'goff','ice',dat_flt(:,col_press)*100,1)/f;
                
                ydat(1).y = 100*qv_dat(:)./qv_sat(:);
                xdat(1).x = time_flt;
                labs(1).l='RHi hygrometer';   
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'RHi humi'

                ylab=['RHi humi (%)'];
                eval(['qv_dat = qv_flt' flight_no '_humi;'])
                qv_sat = satvappress(dat_flt(:,col_temp)+273.15,'goff','ice',dat_flt(:,col_press)*100,1)/f;
                
                ydat(1).y = 100*qv_dat(:)./qv_sat(:);
                xdat(1).x = time_flt;
                labs(1).l='RHi humicap';    
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane                
                
                
            case 'RH hygrometer'

                ylab=['RH hygrometer (%)'];
                eval(['qv_dat = qv_flt' flight_no '_fp;'])
                qv_sat = satvappress(dat_flt(:,col_temp)+273.15,'goff','liq',dat_flt(:,col_press)*100,1)/f;
                
                ydat(1).y = 100*qv_dat(:)./qv_sat(:);
                xdat(1).x = time_flt;
                labs(1).l='RH hygrometer';   
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
            case 'RH humi'

                ylab=['RH humi (%)'];
                eval(['qv_dat = qv_flt' flight_no '_humi;'])
                qv_sat = satvappress(dat_flt(:,col_temp)+273.15,'goff','liq',dat_flt(:,col_press)*100,1)/f;
                
                ydat(1).y = 100*qv_dat(:)./qv_sat(:);
                xdat(1).x = time_flt;
                labs(1).l='RH humi'; 
                
                lor=4;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
                
                
           case 'Airspeed'

                ylab=['Airspeed (m s^{-1})'];
                
                ydat(1).y = dat_flt(:,col_airspeed);
                xdat(1).x = time_flt;
                labs(1).l='Airspeed'; 
                
                ydat(2).y = interp1(CAS_time_all/3600,air_speed_1D'/100,time_flt2);
                xdat(2).x = time_flt;
                labs(2).l='Airspeed CAS'; 
                
                izlim=1;
                zmin=0;
                zmax=90;
                
           case 'Airspeed CAS'

                ylab=['Airspeed (m s^{-1})'];
                
                ydat(1).y = air_speed_1D'/100;
                xdat(1).x = CAS_time_all/3600;
                labs(1).l='Airspeed CAS';    
                
                

        end
        
        idat=length(xdat);                
        
        if iplot_3D==0
            yydat(1).y = ydat(1).y(inds);  %dummy data
        else
            yydat(1).y=yydat(1).y(inds);
        end
        
        for irestrict=1:length(xdat)
                ydat(irestrict).y=ydat(irestrict).y(inds);
                xdat(irestrict).x=xdat(irestrict).x(inds);
        end
        
        xlimits = [time_flt(inds(1)) time_flt(inds(end))];
        
            
        
        if ihighlight_cloud==1
           highlight_cloud_method = {'frost point'};
           highlight_cloud_method = {'CAS'};
           highlight_cloud_method = {'CAS LWC'};
%           highlight_cloud_method = {'CIP'};
           highlight_cloud_method = {'CAS LWC and hotwire','CIP','CAS and hotwire LWC and CIP'};
%           highlight_cloud_method = {'CAS LWC and hotwire'};
            highlight_cloud_method = {'Cloud highlight new'};
            highlight_cloud_method = {'Cloud highlight new ice concs'};            
           
           nmark=zeros([1 length(xdat)+length(highlight_cloud_method)]);
                   
           xdat(3:idat+1)=xdat(2:idat);
           ydat(3:idat+1)=ydat(2:idat); 
%           if ix_distance==1 | ix_distance==2
            if iplot_3D==1
               yydat(3:idat+1)=yydat(2:idat);            
           end
           labs(3:idat+1)=labs(2:idat);
           idat=1;
           
           threshold_HOT=0.02; %g/m3
           lwc_threshold=0.05; %CAS threshold - note these are for unsmoothed data whereas the plots are smoothed
           threshold=0.1*1e-3; %CIP threshold
           threshold=0.01*1e-3; %CIP threshold
           
           for ihighlight=1:length(highlight_cloud_method)
           
            switch highlight_cloud_method{ihighlight}
                case 'frost point'

                    %if either instrument suggests cloud then include it for now
                        icloud = find(dat_flt(:,col_frostpoint_humi)>=dat_flt(:,col_temp));
                        switch flight_no
                            case {'120','122'}  %humicap says that were always in cloud for this flight - faulty or perhaps is true?
                                        %Is very different to the hygrometer reading though.
                                icloud=[];
                        end

                        icloud2 = find(dat_flt(:,col_frostpoint_hygro)>=dat_flt(:,col_temp));
                        
                        ydat(idat+1).y = NaN*ones(size(xdat(1).x)); %default of NaNs (so have blanks for non-cloud)
                        ydat(idat+1).y(icloud) = ydat(1).y(icloud);
                        ydat(idat+1).y(icloud2) = ydat(1).y(icloud2);
                        xdat(idat+1).x = time_flt;
                        
                        labs(idat+1).l='Cloud (frostpoint)';
                        ismooth(idat+1)=0;
                        nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                        
                        idat=idat+1;
                        
                case 'CAS'
                    ncc_threshold=10;
                    icloud = find( CAS_total_number{icas_count} > ncc_threshold );
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    xdat(idat+1).x = CAS_time_all/3600;
                    labs(idat+1).l=['Cloud (CAS>' num2str(ncc_threshold) ' cm^{-3})'];
                    ismooth(idat+1)=0;
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;
                    
                case 'CAS LWC'

                    icloud = find( LWC_dist_cas > lwc_threshold );
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2,ydat(1).y,CAS_time_all(icloud)/3600);                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    
                    test=1;
                    if test==1
                        xdat(idat+1).x = interp1(time_flt2,xdat(1).x,CAS_time_all/3600);
                        yydat(idat+1).y = interp1(time_flt2,yydat(1).y,CAS_time_all/3600);                        
                    else
                        xdat(idat+1).x = interp1(time_flt2,time_flt,CAS_time_all/3600);
                    end
                    
                    labs(idat+1).l=['Cloud (CAS LWC>' num2str(lwc_threshold) ' g m^{-3})'];  
                    ismooth(idat+1)=0;
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;
                    
                 case 'CAS LWC and hotwire'

                                                                                                              
                    icloud_CAS = find( LWC_dist_cas > lwc_threshold );                                           
                    
                    %for hotwire data will probably want to smooth the data as it is quite noisy
                    nfilter=40; bfilter=ones([1 nfilter])*1/nfilter;
                    HOT_smooth = filter(bfilter,1,LWC_CAS_all);
                    HOT_smooth(HOT_smooth<-0.5)=NaN;
                    
                    timeHOT_smooth = CIP_time_all;
                    HOT_smooth(end-nfilter+1:end)=[];
                    timeHOT_smooth(end-nfilter+1:end)=[];
                    HOT_smooth(1:nfilter)=[]; 
                    timeHOT_smooth(1:nfilter)=[];
                    
                    
                     
                    icloud_HOT = find( HOT_smooth > threshold_HOT );                    
                    time_HOT = timeHOT_smooth(icloud_HOT);
                    icloud_HOT = findheight_nearest( CAS_time_all, time_HOT ); %finds indices for all time_CIP values in CAS_time_all array
                                        
                    icloud = intersect(icloud_CAS,icloud_HOT); %finds the common value between the two
                    
                                                                                                                            
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2,ydat(1).y,CAS_time_all(icloud)/3600);                                      
                    yydat(idat+1).y(icloud) = interp1(time_flt2,yydat(1).y,CAS_time_all(icloud)/3600);                     
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    test=1;
                    if test==1
                        xdat(idat+1).x = interp1(time_flt2,xdat(1).x,CAS_time_all/3600);
                    else
                        xdat(idat+1).x = interp1(time_flt2,time_flt,CAS_time_all/3600);
                    end
                    
                    labs(idat+1).l=['Cloud (CAS LWC>' num2str(lwc_threshold) ' g m^{-3})'];  
                    ismooth(idat+1)=0;
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;    
                    
                case 'CIP'
                    %threshold=0.0005;
                    %icloud_CIP = find( CIP_total_number{icas_count} > threshold );
                    

                    CIP_large_sum = sum(CIP_counts_all(:,6:end)./sample_volume_CIP(:,6:end),2);
                    icloud_CIP = find( CIP_large_sum > threshold );
                    
                    time_CIP = CIP_time_all(icloud_CIP);
                    icloud = findheight_nearest( CAS_time_all, time_CIP ); %finds indices for all time_CIP values in CAS_time_all array
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));                    
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2,ydat(1).y,CAS_time_all(icloud)/3600);                    
                    yydat(idat+1).y(icloud) = interp1(time_flt2,yydat(1).y,CAS_time_all(icloud)/3600); 
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2,time_flt,CAS_time_all/3600);  
                    xdat(idat+1).x = interp1(time_flt2,xdat(1).x,CAS_time_all/3600);  
                    
                    labs(idat+1).l=['Cloud (CIP N_{D>137.5 \mum} > ' num2str(1000*threshold) ' L^{-1})'];  
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;        
                    
                case 'CAS LWC and CIP'
                    icloud_CAS = find( LWC_dist_cas > lwc_threshold );
                    
                    
%                    threshold=1e-4; %set earlier
                    CIP_large_sum = sum(CIP_counts_all(:,6:end)./sample_volume_CIP(:,6:end),2);
                    icloud_CIP = find( CIP_large_sum > threshold );
                    
                    time_CIP = CIP_time_all(icloud_CIP);
                    icloud_CIP = findheight_nearest( CAS_time_all, time_CIP ); %finds indices for all time_CIP values in CAS_time_all array
                    
                    
                    icloud = intersect(icloud_CAS,icloud_CIP); %finds the common value between the two
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2,ydat(1).y,CAS_time_all(icloud)/3600);
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2,time_flt,CAS_time_all/3600);


                    labs(idat+1).l=['Both'];  
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;       
                    
                    
               case 'CAS and hotwire LWC and CIP'

                    icloud_CAS = find( LWC_dist_cas > lwc_threshold );
                    
                                        %for hotwire data will probably want to smooth the data as it is quite noisy
                    nfilter=40; bfilter=ones([1 nfilter])*1/nfilter;
                    HOT_smooth = filter(bfilter,1,LWC_CAS_all);
                    HOT_smooth(HOT_smooth<-0.5)=NaN;
                    
                    timeHOT_smooth = CIP_time_all;
                    HOT_smooth(end-nfilter+1:end)=[];
                    timeHOT_smooth(end-nfilter+1:end)=[];
                    HOT_smooth(1:nfilter)=[]; 
                    timeHOT_smooth(1:nfilter)=[];
                                                             
                    icloud_HOT = find( HOT_smooth > threshold_HOT );                    
                    time_HOT = timeHOT_smooth(icloud_HOT);
                    icloud_HOT = findheight_nearest( CAS_time_all, time_HOT ); %finds indices for all time_CIP values in CAS_time_all array
                    
                    
%                    threshold=1e-4; %set earlier
                    CIP_large_sum = sum(CIP_counts_all(:,6:end)./sample_volume_CIP(:,6:end),2);
                    icloud_CIP = find( CIP_large_sum > threshold );                    
                    time_CIP = CIP_time_all(icloud_CIP);
                    icloud_CIP = findheight_nearest( CAS_time_all, time_CIP ); %finds indices for all time_CIP values in CAS_time_all array
                                                                                
                    icloud = intersect(icloud_CAS,icloud_CIP); %finds the common value between the two
                    icloud = intersect(icloud,icloud_HOT); %and when also have hotwire LWC
                    
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2,ydat(1).y,CAS_time_all(icloud)/3600);
                    yydat(idat+1).y(icloud) = interp1(time_flt2,yydat(1).y,CAS_time_all(icloud)/3600);
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2,xdat(1).x,CAS_time_all/3600);


                    labs(idat+1).l=['Both'];  
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;            
                    
             case 'Cloud highlight new'
                 
                    nmark=zeros([1 5]);
                    
                    

                    icloud_CAS = find( LWC_dist_cas > lwc_threshold );
                    
                    nfilter=40; bfilter=ones([1 nfilter])*1/nfilter;
                    HOT_smooth = filter(bfilter,1,LWC_CAS_all);
                    HOT_smooth(HOT_smooth<-0.5)=NaN;
                    
                    timeHOT_smooth = CIP_time_all;
                    HOT_smooth(end-nfilter+1:end)=[];
                    timeHOT_smooth(end-nfilter+1:end)=[];
                    HOT_smooth(1:nfilter)=[]; 
                    timeHOT_smooth(1:nfilter)=[];
                                                             
                    icloud_HOT = find( HOT_smooth > threshold_HOT );                    
                    time_HOT = timeHOT_smooth(icloud_HOT);
                    icloud_HOT = findheight_nearest( CAS_time_all, time_HOT ); %finds indices for all time_CIP values in CAS_time_all array
                    
                    %smoothed hotwire dataz - is this right?
                    
%                    threshold=1e-4; %set earlier
                    CIP_large_sum = sum(CIP_counts_all(:,6:end)./sample_volume_CIP(:,6:end),2);
                    icloud_CIP = find( CIP_large_sum > threshold );                    
                    time_CIP = CIP_time_all(icloud_CIP);
                    icloud_CIP = findheight_nearest( CAS_time_all, time_CIP ); %finds indices for all time_CIP values in CAS_time_all array

% all 3 (i.e. CAS, hotwire and CIP cloud)                    
                    icloud = intersect(icloud_CAS,icloud_CIP); %finds the common values between the two
                    icloud = intersect(icloud,icloud_HOT); %and when also have hotwire LWC                    
                                       
                    
% just CAS and hotwire but not when have all 3                    
                    icloud_CASHOT = intersect(icloud_CAS,icloud_HOT); %finds the common values between the two     
                    [C,IA,IB] = intersect(icloud,icloud_CASHOT);
                    icloud2=icloud_CASHOT;
                    icloud2(IB)='';
                    
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all)); %could just make it the size of xdat(1).x?
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
                    
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);

                     ydat(idat+1).y(icloud2) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(icloud2)/3600);
                     yydat(idat+1).y(icloud2) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(icloud2)/3600);
                    

                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                    labs(idat+1).l=['Cloud (CAS LWC>' num2str(lwc_threshold) ' g m^{-3})'];   
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;    
                    

% just CIP and not when have all 3                    
                    %icloud_CIP = intersect(icloud_CAS,icloud_HOT); %finds the common values between the two     
                    [C,IA,IB] = intersect(icloud,icloud_CIP);
                    icloud2=icloud_CIP;
                    icloud2(IB)='';
                    
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud2) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(icloud2)/3600);
                    yydat(idat+1).y(icloud2) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(icloud2)/3600);
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                    labs(idat+1).l=['Cloud (CIP N_{D>137.5 \mum} > ' num2str(1000*threshold) ' L^{-1})'];   
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;     
                    
                    
                    
                    
                    % plot the all 3 points last
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(icloud)/3600);
                    yydat(idat+1).y(icloud) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(icloud)/3600);
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                    labs(idat+1).l=['Both'];  
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;         
                    
                    
                    if length(time_highlight_path)>0 ...%i.e. i_highlight_path==1
                            switch flt_graph
                                case {'Altitude3D','Temperature3D'}
                                    % highlight a section of the path
                                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));

                                    %                    ihighlight = find( CAS_time_ALL > lwc_threshold );
                                    [ihighlight(1), ihighlight(2)] = findheight(CAS_time_all,time_highlight_path(1)*3600 ...
                                        ,time_highlight_path(2)*3600);

                                    ihighlight=[ihighlight(1):ihighlight(2)];

                                    ydat(idat+1).y(ihighlight) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(ihighlight)/3600);
                                    yydat(idat+1).y(ihighlight) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(ihighlight)/3600);


                                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                                    labs(idat+1).l=['Highlight ' datestr(time_highlight_path(1)/24,13) ' to ' datestr(time_highlight_path(2)/24,13)];
                                    ismooth(idat+1)=0;
                                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)

                                    idat=idat+1;
                                    ismooth(idat)=0;

                            end

                    end
                    
                case 'Cloud highlight new ice concs' %actually decided to base it on when have more than 
                    %one ice crystal - so don't need concentrations at all then
                 
                    nmark=zeros([1 5]);
                    
                    

                    icloud_CAS = find( LWC_dist_cas > lwc_threshold );
                    
                     nfilter=40;
                     [time_HOT_smooth,HOT_smooth]=window_average(CIP_time_all,LWC_CAS_all,nfilter,'mean');
%                    HOT_smooth = filter(bfilter,1,LWC_CAS_all);
                    HOT_smooth(HOT_smooth<-0.5)=NaN;
                    
%                     timeHOT_smooth = CIP_time_all;
%                     HOT_smooth(end-nfilter+1:end)=[];
%                     timeHOT_smooth(end-nfilter+1:end)=[];
%                     HOT_smooth(1:nfilter)=[]; 
%                     timeHOT_smooth(1:nfilter)=[];
                                                             
                    icloud_HOT = find( HOT_smooth > threshold_HOT );                    
                    time_HOT = time_HOT_smooth(icloud_HOT);
                    icloud_HOT = findheight_nearest( CAS_time_all, time_HOT ); %finds indices for all time_CIP values in CAS_time_all array
                    
                    %smoothed hotwire data - is this right? Yes, otherwise
                    %is noisy.
                    
             %now for ice concs
%                    threshold=1e-4; %set earlier
%                     CIP_large_sum = sum(CIP_counts_all(:,6:end)./sample_volume_CIP(:,6:end),2);
%                     icloud_CIP = find( CIP_large_sum > threshold );                    
%                     time_CIP = CIP_time_all(icloud_CIP);
%                     icloud_CIP = findheight_nearest( CAS_time_all, time_CIP ); %finds indices for all time_CIP values in CAS_time_all array

                nsmooth_steps = 60;
                smoothed_plot = 0;
                calc_ice_concs_for_markers  %script

                    CIP_large_sum = ice_dat; %from script above
%                    icloud_CIP = find( CIP_large_sum > threshold*1e3 );       
                    icloud_CIP = find( Ncry_smooth > 0 );
                    time_CIP = time_ice_smooth(icloud_CIP)*3600;
                    icloud_CIP = findheight_nearest( CAS_time_all, time_CIP ); %finds indices for all time_CIP values in CAS_time_all array


                    
                    

% all 3 (i.e. CAS, hotwire and CIP cloud)                    
                    icloud = intersect(icloud_CAS,icloud_CIP); %finds the common values between the two
                    icloud = intersect(icloud,icloud_HOT); %and when also have hotwire LWC                    
                                       
                    
% just CAS and hotwire but not when have all 3                    
                    icloud_CASHOT = intersect(icloud_CAS,icloud_HOT); %finds the common values between the two     
                    [C,IA,IB] = intersect(icloud,icloud_CASHOT);
                    icloud2=icloud_CASHOT;
                    icloud2(IB)='';
                    
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all)); %could just make it the size of xdat(1).x?
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
                    
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);

                     ydat(idat+1).y(icloud2) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(icloud2)/3600);
                     yydat(idat+1).y(icloud2) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(icloud2)/3600);
                    

                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                    labs(idat+1).l=['Liquid (CAS LWC>' num2str(lwc_threshold) ' g m^{-3})'];   
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;    
                    

% just CIP and not when have all 3                    
                    %icloud_CIP = intersect(icloud_CAS,icloud_HOT); %finds the common values between the two     
                    [C,IA,IB] = intersect(icloud,icloud_CIP);
                    icloud2=icloud_CIP;
                    icloud2(IB)='';
                    
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud2) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(icloud2)/3600);
                    yydat(idat+1).y(icloud2) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(icloud2)/3600);
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


%                    labs(idat+1).l=['Cloud (CIP N_{D>112.5 \mum} > ' num2str(1000*threshold) ' L^{-1})'];   
                    labs(idat+1).l=['Ice (CIP N_{D>112.5 \mum} >= 1 ice particle)'];                       
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;     
                    
                    
                    
                    
                    % plot the all 3 points last
                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));
                    
%                    ydat(idat+1).y(icloud) = interp1(xdat(1).x,ydat(1).y,CAS_time_all(icloud)/3600);
                    ydat(idat+1).y(icloud) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(icloud)/3600);
                    yydat(idat+1).y(icloud) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(icloud)/3600);
                    
%                    xdat(idat+1).x = CAS_time_all/3600;
                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                    labs(idat+1).l=['Both'];  
                    ismooth(idat+1)=0;   
                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)
                    
                    idat=idat+1;         
                    
                    
                    if length(time_highlight_path)>0 ...%i.e. i_highlight_path==1
                            switch flt_graph
                                case {'Altitude3D','Temperature3D'}
                                    % highlight a section of the path
                                    ydat(idat+1).y = NaN*ones(size(CAS_time_all));
                                    yydat(idat+1).y = NaN*ones(size(CAS_time_all));

                                    %                    ihighlight = find( CAS_time_ALL > lwc_threshold );
                                    [ihighlight(1), ihighlight(2)] = findheight(CAS_time_all,time_highlight_path(1)*3600 ...
                                        ,time_highlight_path(2)*3600);

                                    ihighlight=[ihighlight(1):ihighlight(2)];

                                    ydat(idat+1).y(ihighlight) = interp1(time_flt2(inds),ydat(1).y,CAS_time_all(ihighlight)/3600);
                                    yydat(idat+1).y(ihighlight) = interp1(time_flt2(inds),yydat(1).y,CAS_time_all(ihighlight)/3600);


                                    xdat(idat+1).x = interp1(time_flt2(inds),xdat(1).x,CAS_time_all/3600);


                                    labs(idat+1).l=['Highlight ' datestr(time_highlight_path(1)/24,13) ' to ' datestr(time_highlight_path(2)/24,13)];
                                    ismooth(idat+1)=0;
                                    nmark(idat+1)=-1; %put markers on the cloud data to make it stand out (thicker line)

                                    idat=idat+1;
                                    ismooth(idat)=0;

                            end

                    end
                    
                    
                    
            end
            
            
           end
            

                
                
                


        end
        

        
        
        leg_choice={'15m','152m','305m','610m'};
        leg_choice={'15m'};        

        if ilegs==1
            idatetick=0;
            iydir=0;
            ixdir=0;
            
            dat=ydat(1).y;
            idat=0;
            
            if sum(strcmp(leg_choice,'15m'))>0
            
            idat=idat+1;
            [it it2]=findheight(time_flt19,20.4,20.57);
            xdat(idat).x = dat_flt(it:it2,col_lat);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '15 m';
            iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
            aws_times = time_flt19(it:it2);
            aws_time(idat) = aws_times(iaws_time(idat));
            
            end
            
            
            switch flt_graph
                case {'Wind','Wind dir'} %cut some of the leg off as there were problems with the GPS lock
                    [it it2]=findheight(time_flt19,20.97,21.07);            
                otherwise                    
                    [it it2]=findheight(time_flt19,20.97,21.17);
            end
            
    if sum(strcmp(leg_choice,'152m'))>0
            idat=idat+1;            
            xdat(idat).x = dat_flt(it:it2,col_lat);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '152 m';
            iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
            aws_times = time_flt19(it:it2);
            aws_time(idat) = aws_times(iaws_time(idat));
    end

    if sum(strcmp(leg_choice,'305m'))>0
            idat=idat+1;
            [it it2]=findheight(time_flt19,21.21,21.38);
            xdat(idat).x = dat_flt(it:it2,col_lat);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '305 m';
            iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
            aws_times = time_flt19(it:it2);
            aws_time(idat) = aws_times(iaws_time(idat));   
    end
            
    if sum(strcmp(leg_choice,'610m'))>0
            idat=idat+1;
            [it it2]=findheight(time_flt19,21.77,21.96);
            xdat(idat).x = dat_flt(it:it2,col_lat);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '610 m';
            iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
            aws_times = time_flt19(it:it2);
            aws_time(idat) = aws_times(iaws_time(idat));     
    end
            
             xlimits=[min(xdat(1).x) max(xdat(1).x)]; %whole range of the flight
%             izlim=1;
%             zmin=0;
%             zmax=20;
             xlab=['Latitude'];
             
             ylab = [ylab ' leg N-S'];
        end
        
        if ilegs==2
            idatetick=0;
            iydir=0;
            ixdir=-1;
            
            dat=ydat(1).y;
            idat=0;
            
        if sum(strcmp(leg_choice,'15m'))>0            
            idat=idat+1;
            [it it2]=findheight(time_flt19,20.57,20.74);
            xdat(idat).x = dat_flt(it:it2,col_lon);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '15 m';
        end
            
        if sum(strcmp(leg_choice,'152m'))>0
            idat=idat+1;
            [it it2]=findheight(time_flt19,20.77,20.97);
            xdat(idat).x = dat_flt(it:it2,col_lon);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '152 m';
        end
            
        if sum(strcmp(leg_choice,'305m'))>0
            idat=idat+1;
            [it it2]=findheight(time_flt19,21.38,21.56);
            xdat(idat).x = dat_flt(it:it2,col_lon);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '305 m';
        end
        
        
        if sum(strcmp(leg_choice,'610m'))>0    
            idat=idat+1;
            [it it2]=findheight(time_flt19,21.59,21.77);
            xdat(idat).x = dat_flt(it:it2,col_lon);
            ydat(idat).y = dat(it:it2);
            labs(idat).l = '610 m';
        end
            
            xlimits=[min(xdat(1).x) max(xdat(1).x)]; %whole range of the flight
%              izlim=1;
%              zmin=0;
%              zmax=20;
             xlab=['Longitude'];
             
              ylab = [ylab ' leg E-W'];
        end
        
%        titlenam=[ylab ' for 6thJan06 Antarctica aircraft data'];
        titlenam=[ylab ' for flight ' flight_no ' (Antarctica aircraft data)'];
        figname=[titlenam];
        savename=titlenam;
        

        for idat=1:length(xdat)

            if ismooth(idat)==1
                [xdat(idat).x,ydat(idat).y] = window_average(xdat(idat).x,ydat(idat).y,nfilter,'mean');
                
                %note - the way that this was being done below was wrong for a smoothing window
                %the x values weren't being properly assigned as the
                %mid-points of the window. E.g. the last entry in the
                %smoothed data is actually the average of the last nfilter
                %datapoints and should not be ignored as below                
%                 bfilter=ones([1 nfilter])*1/nfilter;
%                 ydat(idat).y=filter(bfilter,1,ydat(idat).y);
%                 ydat(idat).y(end-nfilter+1:end)=[]; xdat(idat).x(end-nfilter+1:end)=[];
%                 ydat(idat).y(1:nfilter)=[]; xdat(idat).x(1:nfilter)=[];
            end

        end
        





        maxx=-1e99;
        
        if idatetick==1
                xlimits=xlimits/24;
        end
            
        
        for idat=1:length(ydat)
            if idatetick==1
                xdat(idat).x = xdat(idat).x/24;
            end
            if size(ydat(idat).y,1)==1
                maxx=max([ydat(idat).y maxx]);
            else
                maxx=max([ydat(idat).y' maxx]);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_mac3'

        ihm=2; %q-field code for HM to be plotted

        ylab='Max HM content (g m^{-3})';
        ylab='Max HM content (m^{-3})';

        xlab='Time (mins)';

        dumprange=1:12;
        idirsave=idir;



        dirs=[1];
        for idat=1:length(dirs)
            idir=dirs(idat);

            % hmaxs=SerDan(idir).SER(:,37+nqp+ihm); %heights of the max HM contents
            % smooth_h=hmaxs;

            %         nfilter=3; bfilter=ones([1 nfilter])*1/nfilter;
            %         smooth_h=filter(bfilter,1,hmaxs);

            %         for ihs=1:length(hmaxs)
            %
            %             ih=findheight(GridDan(idir).Z , smooth_h(ihs) ); %find height index of the height of the max HM content
            %             ydat(idat).y(ihs) = SerDan(idir).SER(ihs,37+ihm) *1000 * GridDan(idir).RHO(ih); % convert to g/m3 from kg/kg
            %
            %
            %         end

            for it=dumprange
                ydat(idat).y(it) = maxALL(mac3(it).ni);
            end

            labs(idat).l=[runName(idir).nam];
            xdat(idat).x = dumprange * 5;


        end

        %    titlenam=['Max HM=' num2str(ihm)];
        titlenam=['Max ice no.'];

        figname=[titlenam];

        idir=idirsave;



        yys=[1 2];

        xlims=0;
        izlim=0;


        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=2;
        fsize=16;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 'max_hm_gm3'

        ihm=7; %q-field code for HM to be plotted

        ylab='Max HM content (gm^{-3})';
        ylab='Max HM content (cm^{-3})';
        %	ylab='Max HM content (kg^{-1})';

        xlab='Time (mins)';

        dumprange=1:44;
        idirsave=idir;



        dirs=[1:7];
        for idat=1:length(dirs)
            idir=dirs(idat);

            hmaxs=SerDan(idir).SER(:,37+nqp+ihm); %heights of the max HM contents
            smooth_h=hmaxs;

            %         nfilter=3; bfilter=ones([1 nfilter])*1/nfilter;
            %         smooth_h=filter(bfilter,1,hmaxs);

            for ihs=1:length(hmaxs)

                ih=findheight(GridDan(idir).Z , smooth_h(ihs) ); %find height index of the height of the max HM content
                %            ydat(idat).y(ihs) = SerDan(idir).SER(ihs,37+ihm) *1000 * GridDan(idir).RHO(ih); % convert to g/m3 from kg/kg   - mixing ratio
                ydat(idat).y(ihs) = SerDan(idir).SER(ihs,37+ihm) /1e6 * GridDan(idir).RHO(ih); % convert to #/cm3 from #/kg     - number conc
                %   ydat(idat).y(ihs) = SerDan(idir).SER(ihs,37+ihm);                                                              %  number conc per kg


            end

            labs(idat).l=[runName(idir).nam];
            xdat(idat).x = SerDan(idir).SER(:,1)/60;


        end

        titlenam=['Max HM=' num2str(ihm)];
        figname=[titlenam];

        idir=idirsave;



        yys=[1 2];

        xlims=0;
        izlim=0;


        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=-1;
        fsize=16;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'av_up'

        ylab='Average Updraught (m s^{-1})';

        dumprange=1:44;
        idirsave=idir;


        aind=280; %for ALu
        H=10; %height to plot (km)

        dirs=[1 2];
        for idat=1:length(dirs)
            idir=dirs(idat);

            for ihs=1:1

                ih=findheight((GridDan(idir).Z+620)/1000 , H)+ihs-1;

                area=icediagsALL(idir).i(ih,:,aind);
                area(area==0)=1e99;

                ydat((ihs-1)*2+idat).y = icediagsALL(idir).i(ih,:,[137])./area;   %av updraught ALu_W
                labs((ihs-1)*2+idat).l=[runName(idir).nam ' at ' num2str(GridDan(idir).Z(ih)/1000+0.62,4) ' km'];
                xdat((ihs-1)*2+idat).x = GridDan(idir).t(1:size(icediagsALL(idir).i,2))-16.75;
            end


        end



        titlenam=['Updraught'];
        figname=['Updraught ' direc(idir).dir];

        idir=idirsave;



        yys=[1 2];

        xlims=0;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];


        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'latent_heating'

        ylab='Latent heating rate (K s^{-1})';


        dumprange=1:44;
        idirsave=idir;

        f=1e6*28.97/18;
        Cp=1005;
        rlvap_on_cp = 2.501e6 / Cp;
        rlsub_on_cp = 2.834e6 / Cp ; %diving by Cp so when multiply by mass get temp change (from E=mL=MCpdT, M=mass of air=1 kg as q in kg/kg)

        H=9; %height for timeseries in km

        dirs=[1 2];
        for idat=1:length(dirs)
            idir=dirs(idat);

            for ihs=1:1

                ih=findheight((GridDan(idir).Z+620)/1000 , H)+ihs-1;
                dq_liq=sum(icediagsALL(idir).i(ih,:,[29:30]),3); %sum of liquid and rain sources, summed over time (using ALL_DQ02 etc.)
                dq_ice=sum(icediagsALL(idir).i(ih,:,[31:33]),3); %sum of ice, snow graupel sources, summed over time


                ydat((ihs-1)*2+idat).y = rlvap_on_cp*(dq_liq) + rlsub_on_cp*(dq_ice) ;   %temp change per sec due to latent heat (K/s)

                labs((ihs-1)*2+idat).l=[runName(idir).nam ' at ' num2str(GridDan(idir).Z(ih)/1000+0.62,4) ' km'];
                xdat((ihs-1)*2+idat).x = GridDan(idir).t(1:size(icediagsALL(idir).i,2))-16.75;

            end


        end



        titlenam=['Latent heat timeseries'];
        figname=['Latent heat timeseries' direc(idir).dir];

        idir=idirsave;



        yys=[1 2];

        xlims=0;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];


        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 'mean_vapour'

        ylab='Vapour mixing ratio (ppmv)';





        dumprange=1:44;
        idirsave=idir;

        f=1e6*28.97/18;

        dirs=[1 2];
        for idat=1:length(dirs)
            idir=dirs(idat);

            for ihs=1:5

                ih=findheight((GridDan(idir).Z+620)/1000 , 16)+ihs-1;
                ydat((ihs-1)*2+idat).y = f*squeeze( icediagsALL(idir).i(ih,:,37) );

                labs((ihs-1)*2+idat).l=[runName(idir).nam ' at ' num2str(GridDan(idir).Z(ih)/1000+0.62,4) ' km'];
                xdat((ihs-1)*2+idat).x = GridDan(idir).t(1:size(icediagsALL(idir).i,2));

            end


        end

        titlenam=['Vapour mixing ratio'];
        figname=['Vapour mixing ratio ' direc(idir).dir];

        idir=idirsave;



        yys=[1 2];

        xlims=0;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];


        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 'max_MR_totprc'



        dumprange=1:16;
        idirsave=idir;

        f=1e6*28.97/18;

        dirs=[1 2 3];

        massnum='tot';
        %massnum='vap';

        H=16.5;
        H=16;
        H=17;

        switch massnum
            case 'tot'
                for idat=1:length(dirs)
                    idir=dirs(idat);

                    ih=findheight((GridDan(idir).Z+620)/1000 , H);
                    ydat(idat).y = squeeze( tot_prctiles(idir).t(ih,dumprange,end) );
                    labs(idat).l=runName(idir).nam;
                    xdat(idat).x = GridDan(idir).t(dumprange)+3;
                end

                units=' (kg kg^{-1})';
                ylab='Mixing ratio';
            case 'vap'
                for idat=1:length(idirs)
                    idir=dirs(idat);

                    ih=findheight((GridDan(idir).Z+620)/1000 , H);
                    ydat(idat).y = squeeze( vap_prctiles(1).t(ih,dumprange,end) );
                    labs(idat).l='Vap';
                    xdat(idat).x = GridDan(idir).t(dumprange)+3;
                end

                units=' (kg kg^{-1})';
                ylab='Mixing ratio';
        end


        ylab=[ylab units];
        titlenam=['Max at ' num2str(GridDan(1).Z(ih)/1000+0.62,2) ' km'];
        figname=['Max ' lower(ylab) direc(idir).dir];



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        % xlimits=[GridDan(idir).t(dumprange(23))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    case 'max_MR_hslice'




        dumprange=1:44;
        idirsave=idir;

        f=1e6*28.97/18;

        dirs=[1];
        for idat=1:length(dirs)
            idir=dirs(idat);
            % idir=1;

            massnum='num';
            % massnum='mass';
            %  massnum='maxtemp';

            switch massnum
                case 'mass'
                    ydat(idat).y = squeeze( max(max(iceMR_hslice(idir).dat(:,:,dumprange))) )';
                    labs(idat).l=runName(idir).nam;
                    xdat(idat).x = GridDan(idir).t(dumprange)+3;
                    units=' (kg kg^{-1})';
                    ylab='Mixing ratio';
                case 'num'
                    ydat(idat).y = squeeze( max(max(iceNC_hslice(1).dat(:,:,dumprange))) )';
                    labs(idat).l='Ice';
                    xdat(idat).x = GridDan(idir).t(dumprange)+3-19.75;
                    idat=idat+1;
                    units=' (kg^{-1})';
                    ylab='Number concentration';
                case 'maxtemp'
                    td=potemp_hslice(1).dat(2:end-1,:,dumprange);
                    [sx sy st]=size(td);

                    ih=findheight((GridDan(idir).Z+620)/1000 , 17.8);
                    tref=repmat(GridDan(idir).THREF(ih),[sx sy st]);
                    td=td+tref;

                    P=pressure_hslice(1).dat(2:end-1,:,dumprange);
                    RHOref=repmat(GridDan(idir).RHON(ih),[sx sy st]);
                    Pref=repmat(GridDan(idir).PREFN(ih),[sx sy st]);
                    P=P.*RHOref + Pref;

                    T=td./(1e5./P).^0.286;

                    %             ydat(idat).y = squeeze( max(max(potemp_hslice(1).dat(:,:,dumprange))) )';
                    %             labs(idat).l='Max';
                    %             xdat(idat).x = GridDan(idir).t(dumprange)+3;
                    %             idat=idat+1;
                    %
                    %             ydat(idat).y = squeeze( min(min(potemp_hslice(1).dat(:,:,dumprange))) )';
                    %             labs(idat).l='Min';
                    %             xdat(idat).x = GridDan(idir).t(dumprange)+3;
                    %             idat=idat+1;

                    ydat(idat).y = squeeze( max(max(T)) )';
                    labs(idat).l='Max';
                    xdat(idat).x = GridDan(idir).t(dumprange)+3;
                    idat=idat+1;

                    ydat(idat).y = squeeze( min(min(T)) )';
                    labs(idat).l='Min';
                    xdat(idat).x = GridDan(idir).t(dumprange)+3;
                    idat=idat+1;

                    ydat(idat).y = squeeze( mean(mean(T)) )';
                    labs(idat).l='Mean';
                    xdat(idat).x = GridDan(idir).t(dumprange)+3;
                    idat=idat+1;

                    units=' (K)';
                    ylab='Temperature';
            end

        end

        ylab=[ylab units];
        titlenam=['Max ' lower(ylab)];
        figname=['Max ' lower(ylab) direc(idir).dir];


        %   ydat(idat).y = squeeze( max(max(snowNC_hslice(1).dat)) )';
        %   labs(idat).l='Snow';
        %   xdat(idat).x = GridDan(idir).t(dumprange)+3;
        %   idat=idat+1;


        %    end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3]-19.75;
        % xlimits=[GridDan(idir).t(dumprange(23))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 'max_supersat'

        ylab='Supersaturation (%)';

        titlenam=['Max supersaturation'];
        figname=['Max supersat ' direc(idir).dir];



        dumprange=1:44;
        idirsave=idir;

        f=1e6*28.97/18;

        dirs=[1];
        for idat=1:length(dirs)
            idir=dirs(idat);

            td=potemp_hslice(1).dat(2:end-1,:,dumprange);
            [sx sy st]=size(td);

            ih=findheight((GridDan(idir).Z+620)/1000 , 17.8);
            tref=repmat(GridDan(idir).THREF(ih),[sx sy st]);
            td=td+tref;

            P=pressure_hslice(1).dat(2:end-1,:,dumprange);
            RHOref=repmat(GridDan(idir).RHON(ih),[sx sy st]);
            Pref=repmat(GridDan(idir).PREFN(ih),[sx sy st]);
            P=P.*RHOref + Pref;

            T=td./(1e5./P).^0.286;
            qsi=satvapPress(T,'lem','ice',P,1)/f; %satvappress gives in ppmv if 5th argument=1
            si=100*(vap_hslice(1).dat(2:end-1,:,dumprange)-qsi)./qsi;

            ydat(idat).y = squeeze( max(max(si)) )';

            labs(idat).l=[runName(idir).nam];
            xdat(idat).x = GridDan(idir).t(dumprange)-16.75;


        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3]-19.75;
        xlimits=[GridDan(idir).t(dumprange(25))+3 GridDan(idir).t(dumprange(end))+3]-19.75;

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'echo_top'

        rads=[30 20 15 10 40 35];
        irad=4;


        ylab=['Height of ' num2str(rads(irad)) ' dBZ echotop (km)'];
        titlenam=['Height of ' num2str(rads(irad)) ' dBZ echotop'];

        figname=[titlenam direc(idir).dir];

        dumprange=1:16;
        %    dumprange=7;


        dirs=[1];
        for idat=1:length(dirs)
            idir=dirs(idat);

            dy=diff(GridDan(1).Y1(1:2))/1000;
            for it=dumprange
                itop=find(n10dbz(idir).n(:,irad,it)*dy>=1);  %points where 10 dbz exceeds 1 km distance over domain
                [imax emax]=max(itop);
                if length(imax)>0
                    ydat(idat).y(it) = (GridDan(idir).Z(imax)+620)/1000;
                else
                    ydat(idat).y(it) = 0;
                end
            end

            labs(idat).l=[runName(idir).nam];
            xdat(idat).x = GridDan(idir).t(dumprange)+3;
            xdat(idat).x = GridDan(idir).t(dumprange)-16.75; %change so that is just from model start


        end




        yys=[1 2];

        xlims=0;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(1)+3 GridDan(idir).t(dumprange(end))+3];

        izlim=1;
        zmin=0;
        zmax=22;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

    case 'height_dqvapmax'



        ylab='Height of Max of Vapour Deficit (km)';

        titlenam=['Height of Max DqVap'];
        figname=['Height of Max DqVap ' direc(idir).dir];

        dirs=[1 2];
        ndirs=length(dirs);

        dumprange=1:58;

        idirsave=idir;

        dirs=[3];
        for idat=1:length(dirs)
            idir=dirs(idat);

            [dqmax ih]=max(dq_vaps(idir).d(izmin:izmax,dumprange,2),[],1);
            ydat(idat).y = GridDan(idir).Z(ih)'+620;

            labs(idat).l=[runName(idir).nam];
            xdat(idat).x = GridDan(idir).t(dumprange)+3;


        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        % xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'eddy_flux'

        H=16.25; %height above which to average
        H=16.9;
        H=16;
        %    H=16.5;

        H2=17;

        datind=[355 356]; %v'w'
        datind=[353 354]; %w'th' (ad + sg)

        idir=1; %for purposes of finding height index

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        H2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;

        ylab='Potential Temperature Flux K m s^{-1}';

        titlenam=['Potential Temperature Flux from ' num2str(H,4) ' - ' num2str(H2,4) ' km'];
        titlenam=['Potential Temperature Flux at ' num2str(H,4) ' km'];

        figname=['Timeseries of Potential Temperature Flux ' direc(idir).dir];

        dirs=[1 2];
        ndirs=length(dirs);

        dumprange=1:60;
        idirsave=idir;

        for idat=1:ndirs
            idir=dirs(idat);

            ydat(idat).y = sum(mean(icediagsALL(idir).i(ih,dumprange,datind),1),3) / npess2(idir);

            labs(idat).l=[runName(idir).nam];
            xdat(idat).x = GridDan(idir).t(dumprange)+3;


        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        % xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;


    case 'cumatH_tracer'

        H=16.25; %height above which to average
        H=16.9;
        H=17.2;

        datind=151; %ALL_Q10
        datind=157; %ALu_WQ10
        datind=139; %ALL_WQ10
        %datind=[145]; %ALL_WQSG10
        datind=[145 139]; %ALL_WQSG10

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Tracer at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Tracer Mixing Ratio Sources at H' direc(idir).dir];

        dirs=[1 2];
        ndirs=length(dirs);

        dumprange=1:60;

        for idat=1:1*ndirs
            xdat(idat).x = GridDan(idir).t(dumprange)+3;
        end

        idirsave=idir;


        for idat=1:ndirs
            idir=dirs(idat);

            %ydat(idat).y = cumsum(icediagsALL(idir).i(ih,dumprange,145),2) / npess2(idir);

            dat=-f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,datind),3),GridDan(idir).t,ih-1,ih)/npess2(idir);

            ydat(idat).y=cumsum(dat);

            labs(idat).l=['Tracer Gain ' runName(idir).nam];

            % 		ydat(idir+ndirs).y = -cumsum(fallrate(ih,:),2)*300;
            %         labs(idir+ndirs).l=['Fall Speed Loss ' runName(idir).nam];
            %
            % %		ydat(idir+2*ndirs).y = cumsum(microrate(ih,:),2)*300;
            % %       labs(idir+2*ndirs).l=['Microphysical Gain ' runName(idir).nam];

        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        % xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;


    case 'ice_proc_rates'

        iser=1; %flag to say that want to use SER values instead of icediag values

        if iser==0
            time=GridDan(idir).t(t1:t2)+3;
        else
            time=SerDan(idir).SER(:,1)/3600 + GridDan(1).t(1)+3;
        end


        idir=1;
        H=16.25; %height above which to average

        height=GridDan(idir).Z + ground_heights(idir);


        H0=8.5;

        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        H=GridDan(idir).Z(ih)/1000+0.62;

        ih0=findheight(GridDan(idir).Z/1000+0.62,H0);
        H0=GridDan(idir).Z(ih0)/1000+0.62;

        ih=1:length(GridDan(idir).Z);

        ylab='Rate (ppmv s^{-1})';

        titlenam=['Average updraught in cloudy air at ' num2str(H,3) ' km'];



        vapsour=[24:27];
        vapsink=[999 1 9 30 31];
        %vapsink=999;
        %vapsink=[];

        imrsour=[];
        imrsink=[];

        imrsour=[29:34]; %sources of ice mixing ratio
        imrsink=[8 11 12 27 7 21 36];
        %       imrsour=[29 30 33]; %34 = PIFRW
        %   imrsink=[11 12 36 21];
        %   imrsour=[];
        %  imrsink=[45 46]; %RSAUT RIACI
        %             imrsink=imrncsink;

        smrsour=[8 9 18 11 35 36]; %sources of snow mixing ratio
        smrsink=[25 14 17 19 6];

        gmrsour=[1 12:17 10 19:21]; %sources of graupel mixing ratio
        gmrsink=[2 24 4];

        allsour=[1 13 15 20 16 10 imrsour 9 18 35]; %sources of all ice species - i.e. not including ones that convert from one ice to another
        allsink=[2 4 6 7 24 25 27]; %as above but sinks

        liqsink=[13 3 18 5 32 29 33 34];

        icencsour=[59:62];
        icencsink=[52 53 63 64 45 46];


        snowncsour=[39 45 48];
        snowncsink=[58 50 56 43 40 47];

        snowncsour2=[39 45 48];
        snowncsink2=[43 40 47];

        gncsour=[56 40 41 42];
        gncsink=[57 49];

        % vapsour=[1 9 29 31];
        % vapsour=[24 25 26 27];

        hmssour2=allsour;   %[imrsour smrsour gmrsour];
        hmssink2=allsink;   %[imrsink smrsink gmrsink];

        rainsour_str={'PGMLT','PRAUT','PGSHD','PRACW','PSMLT','PIMLT'};
        rainsink_str={'PGACR','PREVP','PGFR','PSACR','PIACR-G','PIACR-S'};

        rainsour=get_prnum(rainsour_str);
        rainsink=get_prnum(rainsink_str);



        % icencsour=[29 30 33];
        %icencsour=[];
        %icencsink=[11 12 36 21];
        %  icencsink=[];


        liqsour=999;

        clear indexarray
        indexarray(1).i=imrsour;
        indexarray(2).i=imrsink;
        indexarray(3).i=smrsour;
        indexarray(4).i=smrsink;
        indexarray(5).i=gmrsour;
        indexarray(6).i=gmrsink;
        indexarray(7).i=hmssour2;
        indexarray(8).i=hmssink2;
        indexarray(9).i=liqsink;
        indexarray(10).i=vapsour;
        indexarray(11).i=vapsink;
        indexarray(12).i=rainsour;
        indexarray(13).i=rainsink;

        imass_nc=14; %boundary between mass and number processes

        indexarray(imass_nc+1).i=icencsour;
        indexarray(imass_nc+2).i=icencsink;
        indexarray(imass_nc+3).i=snowncsour;
        indexarray(imass_nc+4).i=snowncsink;
        indexarray(imass_nc+5).i=gncsour;
        indexarray(imass_nc+6).i=gncsink;
        indexarray(imass_nc+7).i=snowncsour2;
        indexarray(imass_nc+8).i=snowncsink2;


        hmssour=[];
        hmssink=[];
        indexsour=[]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
        indexsink=[];

        hm='Snow';
        hm='Snownc';
        hm='Ice';
        %       hm='Icenc';
        %       hm='Liquid';
        %    hm='Total ice';
        %       hm='Vapour';
        hm='Graupel';
        %   hm='Graupel NC';

        %   hm='Rain';


        lab=hm;
        switch hm
            case 'Ice'
                indexsour=[1]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[2]; %same for sinks
            case 'Icenc'
                indexsour=[12]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[13]; %same for sinks
            case 'Snownc'
                lab='Snow';
                indexsour=[18]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[19]; %same for sinks
                lor=2;
            case 'Snow'
                lab='Snow';
                indexsour=[3]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[4]; %same for sinks
            case 'Liquid'
                indexsink=[9]; %same for sinks
                hmssour=[999];
            case 'Total ice'
                indexsour=[7]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[8]; %same for sinks
            case 'Vapour'
                indexsour=[10]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[11]; %same for sinks
            case 'Graupel'
                indexsour=[5]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[6]; %same for sinks
            case 'Graupel NC'
                indexsour=[imass_nc+5]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                indexsink=[imass_nc+6]; %same for sinks
                lor=2;
            case 'Rain'
                indexsour=12;
                indexsink=[13];

        end


        if iser==1
            indexarray(indexsour).i=get_pname_col2( { dgs{indexarray(indexsour).i} } ); %takes indices for icediag array and converts to indices
            indexarray(indexsink).i=get_pname_col2( { dgs{indexarray(indexsink).i} } ); %for in SER
        end

        a=indexarray(indexsour).i;
        if length(a)>0; indexarray(indexsour).i(a==0)=[]; end
        a=indexarray(indexsink).i;
        if length(a)>0; indexarray(indexsink).i(a==0)=[]; end



        for iind=1:length(indexsour)
            hmssour=[hmssour indexarray(indexsour(iind)).i];
        end
        for iind=1:length(indexsink)
            hmssink=[hmssink indexarray(indexsink(iind)).i];
        end


        ii=[hmssour hmssink];
        if length(ih)>1
            titlenam=['Microphysical Process Rates for ' num2str(height(ih(1))/1000) ' to ' num2str(height(ih(end))/1000) ' km'];
        else
            titlenam=['Microphysical Process Rates for ' num2str(height(ih(1))/1000) ' km'];
        end

        if max([indexsour indexsink])<imass_nc
            %ylab=[lab ' Microphysical Source Rate (ppmv)'];
            ylab=['Contribution to ' lab ' Mixing Ratio (g kg^{-1} km)'];

            fact=1000;
        else
            %                   xlab=[lab ' Microphysical Number Source Rate (kg^{-1})'];
            ylab=[lab ' No. Source Rate (kg^{-1} km)'];

            fact=1;
        end




        figname=[xlab ' ' titlenam]
        savename=figname;

        set_dgs_numrates; %defines dgs values for numrates option - may need to change, see icediags_5thSept_2005_32
        dgs{999}='PCOND';

        domfact=length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/fact; %multiply by domain size (in km) and covert to grams so is g/kg km
        domfactcont=length(GridDan(idircont).Y1).*diff(GridDan(idircont).Y1(1:2))/fact; %multiply by domain size (in km) and covert to grams so is g/kg km

        %  domfact=1;
        %  domfactcont=1;

        maxval=0;
        for i=1:length(ii)      %contribution to ice from the microphysical rates rather than rates themselves
            if (i==1 & hmssour==999) | ii(i)==999
                if iser==0
                    ydat(i).y = domfact*fact*( TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(ih,t1:t2,liqsink),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) )... %microphys - sums rate over the times given and multiplies by 300s
                        + ( domfact*fact*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(ih,t1:t2,29),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) );
                else
                    ydat(i).y = domfact*fact*( TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(ih,t1:t2,liqsink),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) )... %microphys - sums rate over the times given and multiplies by 300s
                        + ( domfact*fact*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(ih,t1:t2,29),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) );
                end
            else
                if iser==0
                    ydat(i).y = domfact*fact*sum(icediag(idir).i(ih,t1:t2,ii(i)),1)/npess2(idir); %microphys - sums rate over the times given and multiplies by 300s
                else
                    ydat(i).y = domfact*fact*SerDan(idir).SER(:,ii(i));
                end
            end

            xdat(i).x = time;

            if iser==0
                labs(i).l=upper(dgs{ii(i)});
            else
                pnames
                labs(i).l=upper( pname( ii(i) ).p );
            end

        end






        figname=[titlenam direc(idir).dir];


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(t1))+3 GridDan(idir).t(dumprange(t2))+3];
        xlimits=[18 18.6];

        izlim=0;
        zmin=0
        zmax=6e-5;;

        nmark=0;


    case 'emm_ncw'


        ylab='Number concentration (m^{-3})';

        titlenam=['Total droplet number conc'];
        figname=titlenam;

        idirs=[1:length(emmdat)];
        % idirs=[16 17];
        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = sum(emmdat(idirs(idat)).ncw(:,:,1),1);
            ydat(idat).y = max(emmdat(idirs(idat)).ncw(:,:,1),[],1);
            labs(idat).l=run_name_emm{ idirs(idat) };
            if length(ydat(idat).y)>length(xdat(idat).x)
                ydat(idat).y(end)=[];
            end
        end


        xlims=1;
        xlimits=[20.05 20.3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;


    case 'emm_nrwc'


        ylab='Number concentration (m^{-3})';

        titlenam=['Total rain number conc'];
        figname=titlenam;

        idirs=[1:length(emmdat)];
        idirs=[18 19];
        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = sum(emmdat(idirs(idat)).nr(:,:,1),1);
            labs(idat).l=run_name_emm{ idirs(idat) };
        end


        xlims=1;
        xlimits=[20.05 20.3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;


    case 'emm_isg'

        ylab='Mixing ratio (g m^{-3})';
        xlab='UTC time';

        titlenam=['QISG timeseries'];
        figname=titlenam;

        idirs=1:length(emmdat)   %[12 13];
        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = sum(emmdat(idirs(idat)).qisg(:,:,1),1);

            labs(idat).l=run_name_emm{ idirs(idat) };
        end


        xlims=1;
        xlimits=[20.05 20.3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;

    case 'emm_iwc'

        ylab='Mixing ratio (g m^{-3})';
        xlab='UTC time';

        titlenam=['IWC timeseries'];
        figname=titlenam;

        idirs=1:length(emmdat)   %[12 13];
        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = sum(emmdat(idirs(idat)).iwczt(:,:,1),1);

            labs(idat).l=run_name_emm{ idirs(idat) };
        end


        xlims=1;
        xlimits=[20.05 20.3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;


    case 'emm_lwc'

        ih1=findheight(GridDan(idir).Z/1000+0.62,11);
        ih2=findheight(GridDan(idir).Z/1000+0.62,19);

        ylab='Mixing ratio (g m^{-3})';

        titlenam=['Total liquid water content'];
        figname=titlenam;

        idirs=[12 13];
        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = sum(emmdat(idirs(idat)).lwc(:,:,1),1);
            labs(idat).l=run_name_emm{ idirs(idat) };
        end


        xlims=1;
        xlimits=[20.05 20.3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;

    case 'emm_rwc'


        ylab='Mixing ratio (g m^{-3})';

        titlenam=['Total rain content'];
        figname=titlenam;

        idirs=[1:length(emmdat)];
        %idirs=[1 2];

        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = sum(emmdat(idirs(idat)).rwc(:,:,1),1);
            labs(idat).l = run_name_emm{ idirs(idat) };
        end


        xlims=1;
        xlimits=[20.05 20.7];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;


    case 'emm_dcw'

        ylab='Diameter (microns)';
        xlab='UTC time';

        titlenam=['Max droplet diameter'];
        figname=titlenam;

        idirs=1:length(emmdat)   %[12 13];
        idirs=[16 17];
        for idat=1:length(idirs)
            xdat(idat).x = vec(idirs(idat)).time;
            ydat(idat).y = max(emmdat(idirs(idat)).dcw(:,:,1),[],1);
            %		ydat(idat).y = mean(emmdat(idirs(idat)).dcw(:,:,1),1);

            labs(idat).l=run_name_emm{ idirs(idat) };
        end


        xlims=1;
        xlimits=[20.05 20.3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;



    case 'min_icesatMR'

        ih1=findheight(GridDan(idir).Z/1000+0.62,11);
        ih2=findheight(GridDan(idir).Z/1000+0.62,19);


        ylab='Ice saturation mixing ratio (ppmv)';

        titlenam=['Min ice saturation mixing raito'];
        figname=titlenam;

        ndirs=1;
        for idat=1:ndirs
            xdat(idat).x = time(1:42);
        end

        for idir=1:ndirs
            tref=repmat(GridDan(idir).THREF,[1 42]); %ref potemp
            tref=tref./(1e5./pref).^0.286; %ref temp
            pref=repmat(GridDan(idir).PREFN,[1 42]); %ref p
            TT=tref+tpertTimH(1).t;
            ei=SatVapPress(TT,'goff','ice'); %Pa
            sat=0.622*ei./(pref-ei);

            ydat(idir).y=min(f*sat);
            labs(idir).l=[runName(idir).nam];
        end

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        nmark=0;
        lor=1;

    case 'av_rad'

        ih1=findheight(GridDan(idir).Z/1000+0.62,11);
        ih2=findheight(GridDan(idir).Z/1000+0.62,19);


        ylab='Forcing (K day^{-1})';

        titlenam=['Domain average radiative heating rate'];
        figname=['Domain average radiative heating rate'];

        ndirs=1;
        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            rad=icediagsRAD(idir).i(ih1:ih2,dumprange,[1]);
            ydat(idir).y=mean(rad);
            labs(idir).l=[runName(idir).nam];
        end

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;

    case 'LNB_max_dqtot_path'

        H=17.4; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H);


        titlenam=['Domain average cloud ice deposition rate at ' num2str(H,3) ' km'];
        figname=['Domain average cloud ice deposition rate at ' num2str(H,3) direc(idir).dir];



        ndirs=1;




        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            dqtotTimH=length(GridDan(idir).Y1)*( dq_tot(idir).d(1:ih,dumprange,2) ) *dy;
            [maxdq imax]=max(dqtotTimH,[],1); %get height indices for the max dqtot at each time (below 17.4 km)
            for t=dumprange
                %   ydat(idir).y(t)=meanlnb_bel(idir).m(imax(t),t);
                ydat(idir).y(t)=meanlnb_bel_tot(idir).m(imax(t),t);
                %ydat(idir).y(t)=minlnb_vap(idir).m(imax(t),t);
            end
            labs(idir).l=[runName(idir).nam];
        end

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;

    case 'ice_dep_rate'

        idir=1;
        H=16.25; %height above which to average
        H=10.5; %height above which to average
        H=17.9;

        H0=8.5;
        %    H=11.4;
        %  H=17;
        %    H=12.1;
        %    H=6.5;
        %    H=6.9;
        %    H=7.2;

        %    ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        %    H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        H=GridDan(idir).Z(ih)/1000+0.62;

        ih0=findheight(GridDan(idir).Z/1000+0.62,H0);
        H0=GridDan(idir).Z(ih0)/1000+0.62;

        ylab='Rate (ppmv s^{-1})';
        ylab='Rate (g kg^{-1} s^{-1} km)';
        %	ylab='Rate (kg^{-1} s^{-1} km)';
        ylab='Mean Vertical Velocity (m s^{-1})';


        titlenam=['Domain average cloud ice deposition rate at ' num2str(H,3) ' km'];
        %   titlenam=['Domain average Hallet Mossop process rate at ' num2str(H,3) ' km'];
        %   titlenam=['Domain average ice number source rate at ' num2str(H,3) ' km'];
        %   titlenam=['Domain average total ice source rate at ' num2str(H,3) ' km'];

        % titlenam=['Domain average total ice source rate averaged from ' num2str(H0,3) ' to ' num2str(H,3) ' km'];

        %   titlenam=['Average updraught in cloudy air at ' num2str(H,3) ' km'];
        %     titlenam=['Max updraught at ' num2str(H,3) ' km'];

        %    titlenam=['Average ice production over height ' num2str(H,3) ' km'];


        ndirs=1;





        ihm=31; %ice dep
        %  ihm=29; %Hallet Mossop
        ihm=34; %PIFRW
        ihm=3;

        ihm2=[40:42]; %total ice  %q07 =43, 34=ALL_DQ07
        ihm2=34;
        ihm2=[31:33]; %DQ for all ice
        ihm2=[33]; %DQ for all ice

        % ihm2=[76:78]; %= ALu_dq
        % ihm2=302; %ALL_ALu=137 302=Acu_W  303= W>1_W


        aind=285; %285=ACu_A 280=ALu_A 283=W>1_A
        aind=[];



        figname=[titlenam direc(idir).dir];

        idirs=[1];
        for idat=1:length(idirs)
            idir=idirs(idat);
            xdat(idat).x = GridDan(idir).t(dumprange)+3;

            %151 = low tracer
            ydat(idat).y = f*icediag(idir).i(ih,t1:t2,ihm)/npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors

            if length(aind)==1
                area=icediagsALL(idir).i(ih,t1:t2,aind)/npess2(idir);
                ilow=find(area<1e-5);
                %area(area==0)=1;
                area(area==0)=1e99;
            else
                area=1;
            end

            %area=ones(size(area));
            %area(ilow)=1e99;




            ydat(idat).y = f*sum(icediagsALL(idir).i(ih,dumprange,ihm2),3)./area /npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors

            %         ydat(idat).y = f*mean(sum(icediagsALL(idir).i(ih0:ih,t1:t2,ihm2),3),1)./area /npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors
            %
            %         ydat(idat).y = sum(icediagsALL(idir).i(ih,t1:t2,ihm2),3)./area /npess2(idir);
            %
            %       %  ydat(idat).y=MaxW(idir).w(ih,t1:t2);
            %         ydat(idat).y = f*mean(sum(icediag(idir).i(:,t1:t2,ihm),3),1)./area /npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors


            labs(idat).l=[runName(idir).nam];

        end


        %     labs(1).l=['HM ' runName(1).nam];
        %     labs(2).l=['HM ' runName(3).nam];
        %     labs(3).l=['Primary ' runName(1).nam];
        %     labs(4).l=['Primary ' runName(3).nam];


        %     ihm=30; %PIPRM
        %     idir=1;
        %      ydat(3).y = f*icediag(idir).i(ih,t1:t2,ihm)/npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors
        %
        %      idir=3;
        %      ydat(4).y = f*icediag(idir).i(ih,t1:t2,ihm)/npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];


        %  xlimits=[20.05 20.3];

        izlim=0;
        zmin=0
        zmax=6e-5;;

        nmark=0;
        lor=1;

    case 'ice_mic_rate'
        idir=1;
        H=16.25; %height above which to average
        H=10; %height above which to average
        H=11.05;
        %    H=6.5;
        %    H=6.9;
        %    H=7.2;

        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        H=GridDan(idir).Z(ih)/1000+0.62;

        ylab='Rate (ppmv s^{-1})';
        ylab='Rate (g kg^{-1} s^{-1} km)';
        %	ylab='Rate (kg^{-1} s^{-1} km)';

        titlenam=['Domain average cloud ice deposition rate at ' num2str(H,3) ' km'];
        titlenam=['Domain average Hallet Mossop process rate at ' num2str(H,3) ' km'];
        titlenam=['Domain average microphysical mxing ratio source rate at ' num2str(H,3) ' km'];

        figname=[titlenam direc(idir).dir];

        idirs=[1 2];



        ihm=31; %ice dep
        %  ihm=29; %Hallet Mossop
        %  ihm=34; %PIFRW
        ihm=[34];

        for idat=1:length(idirs)*length(ihm)
            xdat(idat).x = time;
        end

        for idat=1:length(idirs)
            idir=idirs(idat)
            %151 = low tracer

            for ihmdat=1:length(ihm)
                iuse=length(ihm)*(idat-1)+ihmdat;
                ydat(iuse).y = 1000*f*icediag(idir).i(ih,t1:t2,ihm(ihmdat))/npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors
                labs(iuse).l=[runName(idir).nam ' ' dgs{ihm(ihmdat)}];
            end
            %  ydat(idat).y = f*icediagsALL(idir).i(ih,t1:t2,34)/npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/f;   %sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors
            %q07 =43, 34=ALL_DQ07



        end


        %     labs(1).l=['HM ' runName(1).nam];
        %     labs(2).l=['HM ' runName(3).nam];
        %     labs(3).l=['Primary ' runName(1).nam];
        %     labs(4).l=['Primary ' runName(3).nam];


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;



    case 'ice_mass'

        H=16.25; %height above which to average
        H=10; %height above which to average
        H=11.5;
        % H=6.5;
        % H=6.9;

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Mixing Ratio (ppmv)';

        titlenam=['Domain average cloud ice mixing ratio at ' num2str(H,3) ' km'];
        figname=['Domain average cloud ice mixing ratio at ' num2str(H,3) direc(idir).dir];

        ndirs=4;

        for idat=1:ndirs
            xdat(idat).x = time;
        end

        ihm=42; %ice mixing ratio
        ihm=[31:33];

        for idir=1:ndirs
            ydat(idir).y = f*sum(icediagsALL(idir).i(ih,dumprange,[ihm]),3)/npess2(idir)* length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000/f; %dividing by no. processors; %dividing by no. processors
            ydat(idir).y=MaxW(idir).w(ih,dumprange);
            labs(idir).l=[runName(idir).nam];
        end

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=2;

    case 'ice_num'

        H=16.25; %height at which to do timeseries
        H=10; %height above which to average
        %  H=11.5;
        H=6.5;
        H=17.9;

        dumprange=1:44;


        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Number Concentration (# kg^{-1})';

        titlenam=['Domain average cloud ice no. conc. at ' num2str(H,3) ' km'];
        figname=['Domain average cloud ice no. conc. at ' num2str(H,3) direc(idir).dir];


        dirs=1;
        for idat=1:length(dirs)
            idir=idirs(idat);
            xdat(idat).x = GridDan(idir).t(dumprange)+3;
            ydat(idat).y = sum(icediagsALL(idir).i(ih,dumprange,[43]),3)/npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; %dividing by no. processors
            %ydat(idat).y = sum(icediagsALL(idir).i(ih,dumprange,[43]),3)/npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; %dividing by no. processors
            labs(idat).l=[runName(idir).nam];
        end

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;



    case 'low_tracer'

        H=16.25; %height above which to average
        H=15.9; %height above which to average

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Mean low tracer (microns)';

        titlenam=['Mean low tracer at ' num2str(H,4) ' km'];
        figname=['Mean low tracer at ' num2str(H,4) direc(idir).dir];

        ndirs=2;

        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            %151 = low tracer
            ydat(idir).y = sum(icediagsALL(idir).i(ih,dumprange,[42]),3)/npess2(idir); %dividing by no. processors
            labs(idir).l=[runName(idir).nam];
        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;

    case 'fall_comp'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Fall speed flux contribution (microns)';

        titlenam=['Fall speed flux contribution at ' num2str(H,4) ' km'];
        figname=['Fall speed flux contribution at ' num2str(H,4) direc(idir).dir];

        ndirs=2;

        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            FallSpeedTimHcalc       %make sure that height is correct in gamdistTimH

            ydat(idir).y = ice_flux(1).i(ih,dumprange); %dividing by no. processors
            labs(idir).l=[runName(idir).nam];
        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;

    case 'mode_diam'
        %    H=16.5; %height above which to average
        H=16.25; %height above which to average
        %  H=16.8; %height above which to average
        %  H=15.6; %height above which to average

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Mass Mode Diameter (microns)';

        titlenam=['Mass mode diameter at ' num2str(H,4) ' km'];
        figname=['Mass mode diameter at ' num2str(H,4) direc(idir).dir];


        idirs=[1 2];
        ndirs=length(idirs);



        for idat=1:ndirs
            idir=idirs(idat);

            xdat(idat).x = GridDan(idir).t(dumprange) - 16.75;


            gamdistTimH       %make sure that height is correct in gamdistTimH
            iend=2800;
            iend=2500;
            %iend=3500;
            %iend=3400;
            d=[D(1):D(iend)/500:D(iend)]*1e6;
            sum_dm=0;
            for it=1:3
                sum_dm=sum_dm+distIce(it).dm(:,dumprange);
            end
            sum_dm=1e-6*f*interp1(D*1e6,sum_dm,d);
            [ac,bc]=max(sum_dm,[],1);

            ydat(idat).y = d(bc);
            labs(idat).l=[runName(idir).nam];
        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3]-19.75;
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;

    case 'av_w'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Average Updraught (ms^{-1})';

        titlenam=['Average updraught at ' num2str(H,4) ' km'];
        figname=['Timeseries of Average updraught at ' num2str(H,4) direc(idir).dir];

        ndirs=2;

        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            ydat(idir).y = icediagsALL(idir).i(ih,dumprange,[137])/npess2(idir); %dividing by no. processors
            labs(idir).l=[runName(idir).nam];
        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=4;


    case 'max_w'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Max Updraught (ms^{-1})';

        titlenam=['Maximum updraught at ' num2str(H,4) ' km'];
        figname=['Timeseries of Maximum updraught at ' num2str(H,4) direc(idir).dir];

        ndirs=2;

        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            ydat(idir).y = MaxW(idir).w(ih,dumprange);  %microicerate is gain of ice or loss of vapour
            labs(idir).l=[runName(idir).nam];
        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        nmark=0;
        lor=1;


    case 'tot_dist'

        ixtime=0;

        ylab='Relative frequency (ppmv^{-1})';
        xlab='Total water mixing ratio (ppmv)';

        titlenam=['Total water frequency distribution'];
        figname=['Total water frequency distribution' direc(idir).dir];

        idirs=[1 2 3];

        for idat=1:length(idirs)
            ydat(idat).y=mean(totdist(idirs(idat)).v(1:200,:),2);  %/max(mean(pxx(idat).p,2));
            xdat(idat).x=[0.1:0.1:20];
            labs(idat).l=[runName(idat).nam];
        end

        logflag=2;

        xlims=1;
        xlimits=([0 10]);

        izlim=0;
        zmin=(1e-3);
        zmax=(1.2);


        nmark=0;
        lor=1;

    case 'vap_dist'

        ixtime=0;

        ylab='Relative frequency (ppmv^{-1})';
        xlab='Vapour mixing ratio (ppmv)';

        titlenam=['Vapour frequency distribution'];
        figname=['Vapour frequency distribution' direc(idir).dir];

        idirs=[1 2 3];

        for idat=1:length(idirs)
            ydat(idat).y=mean(vapdist(idirs(idat)).v,2);  %/max(mean(pxx(idat).p,2));
            xdat(idat).x=[0.1:0.1:20];
            labs(idat).l=[runName(idat).nam];
        end

        logflag=2;

        xlims=0;
        xlimits=([0 10]);

        izlim=0;
        zmin=(1e-3);
        zmax=(1.2);


        nmark=0;
        lor=3;

    case 'gwave_w-spectra'

        ixtime=0;

        H=16.25; %height for spectra
        % H=17.1;
        % H=18.1;

        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        H=GridDan(idir).Z(ih)/1000+0.62;

        ylab='Normalised PSD';
        ylab='PSD (m^2 s^{-1})';
        xlab='w (rad s^{-1})';

        titlenam=['Vertical velocity power spectrum at ' num2str(H,4) ' km'];
        figname=['Vertical velocity power spectrum at H' direc(idir).dir];



        itend=19;

        idirs=[1 2 3];

        for idat=1:length(idirs)
            idat
            dy=300; %data sampled every 300 secs

            if idat==2
                ihminus=100;
            else
                ihminus=0;
            end

            pxsum=0;
            L=length(GridDan(idirs(idat)).Y1);
            for i=1:L
                [px,fx]=periodogram(TwoD_alltim(idirs(idat)).W(ih-ihminus,i,1:itend),[],[],1/dy);
                pxsum=pxsum+px;
            end
            %psd is power per unit frequency so here is (m/s)^2 / (m^-1) since frequency is 1/m (grid data in metres)

            ydat(idat).y=pxsum/L;  %/max(mean(pxx(idat).p,2));
            xdat(idat).x = 2*pi*fx;
            labs(idat).l=[runName(idat).nam];
        end

        logflag=12;

        xlims=1;
        xlimits=([7e-5 2e-2]);

        izlim=0;
        zmin=(1e-3);
        zmax=(1.2);


        nmark=0;
        lor=3;

    case 'gwave_k-spectra'

        ixtime=0;

        H=16.25; %height for spectra
        % H=17.1;
        % H=18.1;

        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        H=GridDan(idir).Z(ih)/1000+0.62;

        ylab='Normalised PSD';
        ylab='PSD (m^3 s^{-2})';
        xlab='k (rad m^{-1})';
        xlab='Horizontal wavelength (km)';

        titlenam=['Vertical velocity power spectrum at ' num2str(H,4) ' km'];
        figname=['Vertical velocity power spectrum at H' direc(idir).dir];



        itend=19;

        idirs=[1 2 3];

        for idat=1:length(idirs)
            dy=diff(GridDan(idirs(idat)).Y1(1:2));

            if idat==2
                ihminus=100;
            else
                ihminus=0;
            end

            for i=1:itend
                [pxx(idat).p(:,i),fxx(idat).f(:,i)]=periodogram(TwoD_alltim(idirs(idat)).W(ih-ihminus,:,i),[],[],1/dy);
            end
            %psd is power per unit frequency so here is (m/s)^2 / (m^-1) since frequency is 1/m (grid data in metres)

            ydat(idat).y=mean(pxx(idat).p,2)';  %/max(mean(pxx(idat).p,2));
            xdat(idat).x = 2*pi*fxx(idat).f(:,1);
            %xdat(idat).x = 1/fxx(idat).f(:,1);
            labs(idat).l=[runName(idat).nam];
        end

        logflag=12;

        xlims=1;
        xlimits=([5e-5 10e-3]);

        izlim=0;
        zmin=(1e-3);
        zmax=(1.2);

        ixdir=0;


        nmark=0;
        lor=3;

    case 'cumatHchangeiceNC_resdiff'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Change (kg^{-1})';

        titlenam=['Change in Number Concentration at ' num2str(H,4) ' km'];
        figname=['Timeseries of Number Concentration Changes at H' direc(idir).dir];

        ndirs=2;

        for idat=1:ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            ad_calcs4timeseries;

            ydat(idir).y = changenc(ih,dumprange);  %microicerate is gain of ice or loss of vapour
            labs(idir).l=[runName(idir).nam];

        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;


    case 'cumatHchangeice_resdiff'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Change (ppmv)';

        titlenam=['Change in Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Mixing Ratio Changes at H' direc(idir).dir];

        ndirs=2;

        for idat=1:2*ndirs
            xdat(idat).x = time;
        end


        runName2(1).nam='Control';
        runName2(2).nam='960 cm^{-3}';
        for idir=1:ndirs
            ad_calcs4timeseries;

            ydat(idir).y = changevap(ih,dumprange);  %microicerate is gain of ice or loss of vapour
            %        labs(idir).l=['Vapour, ' runName(idir).nam];
            labs(idir).l=['Vapour, ' runName2(idir).nam];

            ydat(idir+ndirs).y = changeice(ih,dumprange);  %microicerate is gain of ice or loss of vapour
            %        labs(idir+ndirs).l=['Ice, ' runName(idir).nam];
            labs(idir+ndirs).l=['Ice, ' runName2(idir).nam];

        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %  xlimits=[GridDan(idir).t(dumprange(36))+3 GridDan(idir).t(dumprange(end))+3];
        %  xlimits=[GridDan(idir).t(dumprange(20))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=1;
        zmin=-0.5
        zmax=0.75;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;

    case 'cumatHdepsub_resdiff'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        H=GridDan(idir).Z(ih)/1000+0.62;


        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Vapour Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Vapour Mixing Ratio Sources at H' direc(idir).dir];

        ndirs=2;

        for idat=1:2*ndirs
            xdat(idat).x = time;
        end

        for idir=1:ndirs
            ad_calcs4timeseries;

            ydat(idir).y = f*cumsum( sum(icediag(idir).i(ih,dumprange,[1 9 31]) ,3) ,2)*300/npes;  %microicerate is gain of ice or loss of vapour
            labs(idir).l=['Ice Dep. ' runName(idir).nam];

            ydat(idir+ndirs).y = f*cumsum( sum(icediag(idir).i(ih,dumprange,[24 25 27]) ,3) ,2)*300/npes; %vapadcum is the cumlative advective loss of vapour
            labs(idir+ndirs).l=['Ice Sub. ' runName(idir).nam];

        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=2;

    case 'cumatHvap_resdiff'

        H=16.25; %height above which to average
        H=16.2;
        %    H=17;

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        H2=17.5;
        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        H2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;


        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Vapour Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Vapour Mixing Ratio Sources at H' direc(idir).dir];

        ndirs=2;


        for idat=1:2*ndirs
            xdat(idat).x = GridDan(1).t(dumprange)+3;
        end

        for idir=1:ndirs
            npes=npess2(idir);
            ad_calcs4timeseries;

            ydat(idir).y = cumsum(microicerate(ih,:),2)*300; %microicerate is gain of ice or loss of vapour
            ydat(idir).y = cumsum(sum(microicerate(ih:ih2,:)),2)*300; %microicerate is gain of ice or loss of vapour

            labs(idir).l=['Microphysical Loss' runName(idir).nam];

            ydat(idir+ndirs).y = -vapadcum(ih,:); %vapadcum is the cumlative advective loss of vapour
            ydat(idir+ndirs).y = -sum(vapadcum(ih:ih2,:)); %vapadcum is the cumlative advective loss of vapour
            %       ydat(idir+ndirs).y = f*cumsum(icediagsALL(idir).i(ih,dumprange,1),2) / npess2(idir) ;

            labs(idir+ndirs).l=['Advective Gain' runName(idir).nam];

        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=-1;

    case 'cumatHtotiM_resdiff'

        H=16.25; %height above which to average
        H=17.5;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        H2=17.5;
        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        H2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Total Ice Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Ice Mixing Ratio Sources at H' direc(idir).dir];

        time=GridDan(1).t(dumprange)+3;

        ndirs=1;

        for idat=1:3*ndirs
            xdat(idat).x = time;
        end

        idirsave=idir;

        for idir=1:ndirs
            ad_calcs4timeseries;
            ydat(idir).y = iceadcum(ih,:);
            %      ydat(idir).y = sum(iceadcum(ih:ih2,:));

            labs(idir).l=['Adv. Gain ' runName(idir).nam];

            ydat(idir+ndirs).y = -cumsum(fallrate(ih,:),2)*300;
            %		ydat(idir+ndirs).y = -cumsum(sum(fallrate(ih:ih2,:)),2)*300;

            labs(idir+ndirs).l=['Fall Speed Loss ' runName(idir).nam];

            ydat(idir+2*ndirs).y = cumsum(microicerate(ih,:),2)*300;
            labs(idir+2*ndirs).l=['Microphysical Gain ' runName(idir).nam];

        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;


    case 'cumatHtotw_resdiff'

        H=16.25; %height above which to average
        %  H=17.2;
        H=16.5;
        H=16.9;

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        H2=17.5;
        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        H2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;


        % ih=ih-1+izmin;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['ALd Contributions to Total Water Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Water Mixing Ratio Sources at H' direc(idir).dir];

        time=GridDan(1).t(dumprange)+3;

        ndirs=2;

        for idat=1:2*ndirs
            xdat(idat).x = time;
        end

        idirsave=idir;

        for idir=1:ndirs
            ad_calcs4timeseries;
            %        ydat(idir).y = iceadcum(ih,:);
            ydat(idir).y = sum( ad(ih:ih2,:) );     %[10:14]
            %    ydat(idir).y = [0 diff( ad(ih,:) )  / 300 ];     %[10:14]

            %    ydat(idir).y = cumsum(  f * (  sum( icediagsALL(idir).i(ih,dumprange,[1:6]),3)   ) / npess2(idir)   ,2);

            %    ydat(idir).y = cumsum(  f * (  sum( icediagsALL(idir).i(ih,dumprange,[46:51]),3)   ) / npess2(idir)   ,2);

            %    ydat(idir).y = cumsum(  f * (  sum( icediagsALL(idir).i(ih,dumprange,[91:96]),3)   ) / npess2(idir)   ,2);

            labs(idir).l=['Adv. Gain ' runName(idir).nam];

            ydat(idir+ndirs).y = -cumsum(  sum ( fallrate(ih:ih2,:) ) ,2)*300 ;
            %ydat(idir+ndirs).y = -cumsum(fallrate(ih,:),2)*300;
            %  ydat(idir+ndirs).y = cumsum( f * (sum( icediagsALL(idir).i(ih,dumprange,[19:24]),3)    ) / npess2(idir)  ,2);

            labs(idir+ndirs).l=['Fall Speed Loss ' runName(idir).nam];

            %		ydat(idir+2*ndirs).y = cumsum(microrate(ih,:),2)*300;
            %       labs(idir+2*ndirs).l=['Microphysical Gain ' runName(idir).nam];

        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;

    case 'cumatHtotwAD_resdiff'

        H=16.25; %height above which to average
        %   H=17.2;
        %   H=16.5;
        H=16.9;

        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        H2=17.5;
        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        H2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;


        ih=ih-1+izmin;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Total Water Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Water Mixing Ratio Sources at H' direc(idir).dir];

        time=GridDan(1).t(dumprange)+3;

        ndirs=2;


        ndiags=4;

        for idat=1:ndiags*ndirs
            xdat(idat).x = time;
        end

        idirsave=idir;

        for idir=1:ndirs
            ad_calcs4timeseries;
            %        ydat(idir).y = iceadcum(ih,:);
            ydat(idir).y = sum(ad(ih:ih2,:));                                   %[10:14]
            ydat(idir).y = cumsum(  f * (  sum( icediagsALL(idir).i(ih,dumprange,[1:6]),3)   ) / npess2(idir)   ,2);

            ydat(idir).y = cumsum(  f * (  sum( icediagsALL(idir).i(ih,dumprange,[46:51]),3)   ) / npess2(idir)   ,2); %Alu

            ydat(idir).y =   f * (  sum( icediagsALL(idir).i(ih,dumprange,[91:96 100:105]),3)   ) / npess2(idir)  ;   %ALd
            %        ydat(idir).y =cumsum(   f * (  sum( icediagsALL(idir).i(ih,dumprange,[91:96 100:105]),3)   ) / npess2(idir)  ,2);   %ALd


            labs(idir).l=['ALd ' runName(idir).nam];

            ydat(idir+ndirs).y = -cumsum(fallrate(ih,:),2)*300;
            %ydat(idir+ndirs).y = -cumsum(fallrate(ih,:),2)*300;
            ydat(idir+ndirs).y = cumsum( f * (sum( icediagsALL(idir).i(ih,dumprange,[19:24]),3)    ) / npess2(idir)  ,2);

            ydat(idir+ndirs).y =f * (  sum( icediagsALL(idir).i(ih,dumprange,[46:51 55:60]),3)   ) / npess2(idir)  ;   %ALd
            %   ydat(idir+ndirs).y =cumsum(   f * (  sum( icediagsALL(idir).i(ih,dumprange,[46:51 55:60]),3)   ) / npess2(idir)   ,2);   %ALd

            labs(idir+ndirs).l=['ALu ' runName(idir).nam];

            ydat(idir+2*ndirs).y = f * (  sum( icediagsALL(idir).i(ih,dumprange,[19:24]),3)   ) / npess2(idir)  ;   %ALd
            labs(idir+2*ndirs).l=['FALL ' runName(idir).nam];


            ydat(idir+3*ndirs).y = f * (  sum( icediagsALL(idir).i(ih,dumprange,[1:6]),3)   ) / npess2(idir)  ;   %ALd
            labs(idir+3*ndirs).l=['ALL ' runName(idir).nam];


            %		ydat(idir+2*ndirs).y = cumsum(microrate(ih,:),2)*300;
            %       labs(idir+2*ndirs).l=['Microphysical Gain ' runName(idir).nam];

        end

        idir=idirsave;



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(30))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;




    case 'cumatHtotiN_resdiff'

        H=16.25; %height above which to average
        H=17.5;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (kg^{-1})';

        titlenam=['Contributions to Total Ice Number at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Number Sources at H' direc(idir).dir];

        ndirs=1;

        for idat=1:3*ndirs
            xdat(idat).x = GridDan(1).t(dumprange)+3;
        end

        for idir=1:ndirs
            ad_calcs4timeseries;
            ydat(idir).y = ncadcum(ih,:);
            labs(idir).l=['Advective Gain ' runName(idir).nam];

            ydat(idir+ndirs).y = cumsum(fallnc(ih,:),2)*300;
            labs(idir+ndirs).l=['Fall Speed Gain ' runName(idir).nam];

            ydat(idir+2*ndirs).y = cumsum(micronc(ih,:),2)*300;
            labs(idir+2*ndirs).l=['Microphysical Gain ' runName(idir).nam];

        end



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215;
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=-1;


    case 'cumatHicetypes'

        H=16.95; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Ice Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Mixing Ratios of Different Ice Types Sources at H' direc(idir).dir];

        for idat=1:6
            xdat(idat).x = time;
        end

        %     idir=1;
        ad_calcs4timeseries;
        %     ncadcum_cont = ncadcum;
        %     fallnc_cont = fallnc;
        %     micronc_cont = micronc;

        ydat(1).y = HMiceadcum(ih,t1:t2);
        ydat(2).y = HMsnowadcum(ih,t1:t2);
        ydat(3).y = HMgraadcum(ih,t1:t2);

        ydat(4).y = -cumsum(HMicefallrate(ih,t1:t2),2)*300;
        ydat(5).y = -cumsum(HMsnowfallrate(ih,t1:t2),2)*300;
        ydat(6).y = -cumsum(HMgrafallrate(ih,t1:t2),2)*300;


        labs(1).l='Advective Gain Ice';
        labs(2).l='Advective Gain Snow';
        labs(3).l='Advective Gain Graupel';

        labs(4).l='Fall Speed Loss Ice';
        labs(5).l='Fall Speed Loss Snow';
        labs(6).l='Fall Speed Loss Graupel';

        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=4;

    case 'cumatHicetypesMicro'

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Ice Mixing Ratio at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Mixing Ratios of Different Ice Types Sources at H' direc(idir).dir];

        for idat=1:3
            xdat(idat).x = time;
        end



        ydat(1).y = cumsum(HMicerate(ih,t1:t2),2)*300;
        ydat(2).y = cumsum(HMsnowrate(ih,t1:t2),2)*300;
        ydat(3).y = cumsum(HMgrarate(ih,t1:t2),2)*300;



        labs(1).l='Microphysical Gain Ice';
        labs(2).l='Microphysical Gain Snow';
        labs(3).l='Microphysical Gain Graupel';



        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;


    case 'minchangetot'
        ad_calcs4timeseries;

        H=16.0; %height above which to averageN
        H2=17;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;
        H2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;

        ylab='Change in Total Water (ppmv)';

        titlenam=['Min Change in Tot Water from ' num2str(H,4) ' to ' num2str(H2,4) ' km'];
        figname=titlenam;

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(1).y = min(change(ih:ih2,t1:t2));

        labs(1).l='Min Change in Total Water';

        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(dumprange(t1))+3 GridDan(idir).t(t2)+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;

    case 'changenc'
        ad_calcs4timeseries;

        H=16.25; %height above which to averageN
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Ice Number (kg^{-1})';

        titlenam=['Contributions to Total Ice Number at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Number Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(1).y = changenc(ih,:);

        labs(1).l='Ice Number';

        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;

    case 'cumatHtotiN'
        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (kg^{-1})';

        titlenam=['Contributions to Total Ice Number at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Number Sources at H' direc(idir).dir];

        for idat=1:3
            xdat(idat).x = time;
        end

        ydat(1).y = ncadcum(ih,:);
        ydat(2).y = cumsum(fallnc(ih,:),2)*300;
        ydat(3).y = cumsum(micronc(ih,:),2)*300;

        labs(1).l='Advective Gain';
        labs(2).l='Fall Speed Gain';
        labs(3).l='Microphysical Gain';

        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;

    case 'cumatHtotiN_diff'

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (kg^{-1})';

        titlenam=['Contributions to Total Ice Number at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Number Sources at H' direc(idir).dir];

        for idat=1:6
            xdat(idat).x = time;
        end

        idir=1;
        ad_calcs4timeseries;
        ncadcum_cont = ncadcum;
        fallnc_cont = fallnc;
        micronc_cont = micronc;

        idir=2;
        ad_calcs4timeseries;

        ydat(1).y = ncadcum(ih,:);
        ydat(2).y = ncadcum_cont(ih,:);

        ydat(3).y = cumsum(fallnc(ih,:),2)*300;
        ydat(4).y = cumsum(fallnc_cont(ih,:),2)*300;

        ydat(5).y = cumsum(micronc(ih,:),2)*300;
        ydat(6).y = cumsum(micronc_cont(ih,:),2)*300;


        labs(1).l='Advective Gain CCN 960';
        labs(2).l='Advective Gain';

        labs(3).l='Fall Speed Gain CCN 960';
        labs(4).l='Fall Speed Gain';

        labs(5).l='Microphysical Gain CCN 960';
        labs(6).l='Microphysical Gain';


        yys=[1 2];

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];

        izlim=0;
        zmin=215
        zmax=220;

        maxx=-1e99;
        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;
        lor=1;


    case 'microvapsoursinks'

        %      logflag=2;
        %      lor=4;

        % ad_calcs4timeseries;

        H=16.15; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        ylab='Microphysical Rate (ppmv s^{-1})';

        titlenam=['Vapour Microphysical Source and Sink'];
        figname=['Vapour Micropohysics Source Sink ' direc(idir).dir];

        for idat=1:2
            xdat(idat).x=GridDan(idir).t(dumprange)+3;
        end

        pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[24 25 27]),3);
        ydat(1).y = pdat(1).p(ih,dumprange);

        pdat(1).p=f*sum(icediag(idir).i(izmin:izmax,dumprange,[1 9 31]),3);
        ydat(2).y = pdat(1).p(ih,dumprange);

        labs(1).l='Source of Vapour';
        labs(2).l='Sink of Vapour';

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

    case 'dqtotsum'
        H=1; %height above which to average
        %    ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        ih=findheight(GridDan(1).Z/1000+0.62,H);
        HH=GridDan(1).Z(ih+izmin-1)/1000+0.62;

        H2=25; %height above which to average
        %    ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        ih2=findheight(GridDan(1).Z/1000+0.62,H2);
        HH2=GridDan(1).Z(ih2+izmin-1)/1000+0.62;

        titlenam=['Sum of Tot Water Points from ' num2str(HH,4) ' to ' num2str(HH2,4) ' km'];
        titlenam=['Sum of Tot Water Points Less Than 5ppmv up to ' num2str(HH2,4) ' km'];
        figname=['Sum of Total Water Points ' direc(idir).dir];
        titlenam='';

        idirecs=[1 2];
        for idat=1:length(idirecs)
            idir=idirecs(idat);
            dy=diff(GridDan(idir).Y1(1:2));
            xdat(idat).x=GridDan(idir).t(dumprange)+3;
            air=repmat(GridDan(idir).RHON(2:end).*diff(GridDan(1).Z),[1 length(dumprange)]);
            pdat(idat).p=air.*length(GridDan(idir).Y1).*( dq_tot(idir).d(2:end,dumprange,2) ) .*dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more

            ydat(idat).y = sum(pdat(idat).p(ih:ih2,dumprange),1);


            air=repmat(GridDan(idir).RHON(ih:ih2).*diff(GridDan(idir).Z(ih-1:ih2)),[1 length(dumprange)]);

            %2d case
            pdat(1).p= air .* length(GridDan(idir).Y1).*( dq_vaps(idir).d(ih:ih2,dumprange,2) )/f *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more

            %3d case
            pdat(1).p= air .* length(GridDan(idir).Y1).*( dq_tot(idir).d(ih:ih2,dumprange,2) )/f *dy*dy; %

            ydat(idat).y = sum(pdat(1).p(:,dumprange),1);

            labs(idat).l=runName(idir).nam;
        end


        %labs(2).l='Fall Speed Loss';

        labs(2).l='2km 2d (for 1km in 3rd dim)';



        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];

        ylab='Sum (ppmv km)';
        ylab='Total water mass removed from below 5 ppmv (kg)';

        lor=2;

    case 'dqvapsum'
        idir=2;

        H=1; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H);
        HH=GridDan(idir).Z(ih)/1000+0.62;

        H2=17.7; %height above which to average
        H2=25; %height above which to average

        ih2=findheight(GridDan(idir).Z/1000+0.62,H2);
        HH2=GridDan(idir).Z(ih2)/1000+0.62;

        titlenam=['Sum of Vapour Points from ' num2str(HH,4) ' to ' num2str(HH2,4) ' km'];
        titlenam=['Sum of Vapour Points up to ' num2str(HH2,4) ' km'];
        figname=['Timerseries of Sum of Deficit Vapour Points ' direc(idir).dir];
        titlenam='';



        if ismember(loadselect(idir),jshifts)
            it1=3;
        else
            it1=1;
        end

        idirecs=[1 2];
        for idat=1:length(idirecs)
            idir=idirecs(idat);
            xdat(idat).x=GridDan(idir).t(dumprange)+3;
            dy=diff(GridDan(idir).Y1(1:2));
            air=repmat(GridDan(idir).RHON(ih:ih2).*diff(GridDan(idir).Z(ih-1:ih2)),[1 length(dumprange)]);

            %2d case
            pdat(1).p= air .* length(GridDan(idir).Y1).*( dq_vaps(idir).d(ih:ih2,dumprange,2) )/f *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more

            %3d case
            pdat(1).p= air .* length(GridDan(idir).Y1).*( dq_vaps(idir).d(ih:ih2,dumprange,2) )/f *dy*dy; %

            ydat(idat).y = sum(pdat(1).p(:,dumprange),1);
            %ydat(2).y = -cumsum(fallrate(ih,:),2)*300;

            labs(idat).l=runName(idir).nam;
            %labs(2).l='Fall Speed Loss';
        end

        labs(2).l='2km 2d (for 1km in 3rd dim)';




        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];

        ylab='Total vapour mass removed from below 5 ppmv (kg m^{-1})';
        ylab='Total vapour mass removed from below 5 ppmv (kg)';

        lor=1;

    case 'rhopert_vap'
        idir=2;

        H=16; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;
        HH=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        H2=17; %height above which to average
        ih2=findheight(GridDan(idir).Z/1000+0.62,H2)-izmin+1;
        HH2=GridDan(idir).Z(ih2+izmin-1)/1000+0.62;

        titlenam=['Sum of Density Perturbations ' num2str(HH,4) ' to ' num2str(HH2,4) ' km'];
        titlenam=['Sum of Density Perturbations up to ' num2str(HH2,4) ' km'];
        figname=['Sum of Vapour Points ' direc(idir).dir];

        if ismember(loadselect(idir),jshifts)
            it1=3;
        else
            it1=1;
        end

        clear diff

        idirecs=[2 4];
        for idat=1:length(idirecs)
            idir=idirecs(idat);
            xdat(idat).x=GridDan(idir).t(dumprange)+3;
            dy=diff(GridDan(idir).Y1(1:2));
            air=repmat(GridDan(idir).RHON(ih:ih2).*diff(GridDan(idir).Z(ih:ih2+1)),[1 length(dumprange)]);

            air=1; %should the mean density perturbation (not total amount of cold air) be weighted according to air mass variations
            %with height?
            %not multiplying by dy as is the mean density pert in < 5 ppmv air
            pdat(1).p= air .* ( rho5ppmv_vap(idir).r(ih:ih2,dumprange) );
            pdat(1).p(isnan(pdat(1).p))=0;
            ydat(idat).y = sum(pdat(1).p,1);
            %ydat(2).y = -cumsum(fallrate(ih,:),2)*300;

            labs(idat).l=runName(idir).nam;
            %labs(2).l='Fall Speed Loss';
        end




        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];

        ylab='Sum (ppmv km)';

        lor=0;


    case 'rhochange'

        %      logflag=2;
        %      lor=4;

        % ad_calcs4timeseries;

        H=16.75; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        ylab='Density Change (kgm^{-3})';

        titlenam=['Density Change at ' num2str(H) ' km'];
        figname=['Timeseries of denisty change' direc(idir).dir];

        for idat=1:1
            xdat(1).x=GridDan(idir).t(dumprange)+3;
        end

        ydat(1).y = pdat(1).p(ih,dumprange);
        %ydat(2).y = -cumsum(fallrate(ih,:),2)*300;

        labs(1).l='Advective Gain';
        %labs(2).l='Fall Speed Loss';


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];


    case 'cumabvHtoti'

        ad_calcs4timeseries;

        H=17; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Average Above ' num2str(H) ' km for Total Ice'];
        figname=['Timeseries of Cumlative Sources abv H' direc(idir).dir];

        for idat=1:2
            xdat(idat).x = time;
        end

        ydat(1).y = mean(iceadcum(ih:end,:),1);
        ydat(2).y = -mean(cumsum(fallrate(ih:end,:),2),1)*300;

        labs(1).l='Advective Gain';
        labs(2).l='Fall Speed Loss';


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %

        nmark=0;

    case 'cumabvHtotw'

        ad_calcs4timeseries;

        H=16; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H); %index in original Grid.Z
        ih2=ih-izmin+1;   %index for H in change, etc. that run from izmin:izmax
        z=GridDan(idir).Z(ih:end);

        ylab='Cumulative Contribution to Total Water (ppmv)';

        titlenam=['Average Above ' num2str(H) ' km for Total Water'];
        figname=['Timeseries of Cumlative Sources abv H' direc(idir).dir];

        for idat=1:2
            xdat(idat).x = time;
        end

        dat=ad(ih2+1:end,:); %ad is a cumulative sum

        rho=repmat(GridDan(idir).RHON(ih+1:end),[1 size(dat,2)]);
        dzz=repmat(diff(z),[1 size(dat,2)]);
        mtot=sum(rho.*dzz);

        ydat(1).y = mean(dat.*rho.*dzz,1)./mtot;

        dat2=cumsum(fallrate(ih2+1:end,:),2)*300;
        %ydat(1).y = sum(change(ih:end,:),1);
        ydat(2).y = -mean(dat2.*rho.*dzz,1)./mtot;

        labs(1).l='Advective Gain';
        labs(2).l='Fall Speed Loss';


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(dumprange(35))+3 GridDan(idir).t(dumprange(end))+3];
        %

        nmark=0;
    case 'cumatHtotw_diff'

        %      logflag=2;
        lor=4;


        H=16.25; %height above which to average
        %    H=16;
        % H=16.25;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution to Total Water (ppmv)';

        titlenam=['Contributions to Total Water at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        idir=1 ;
        ad_calcs4timeseries;
        ad_cont=ad;
        fallrate_cont=fallrate;

        idir=2;
        ad_calcs4timeseries;


        for idat=1:4
            xdat(idat).x = time;
        end

        %ydat(1).y = ad(ih,:)-ad_cont(ih,:);
        %ydat(2).y = -cumsum(fallrate(ih,:),2)*300 + cumsum(fallrate_cont(ih,:),2)*300;

        ydat(1).y = ad(ih,:);
        ydat(2).y=  ad_cont(ih,:)
        ydat(3).y = -cumsum(fallrate(ih,:),2)*300;
        ydat(4).y = -cumsum(fallrate_cont(ih,:),2)*300;


        labs(1).l='Advective Gain CCN 960';
        labs(2).l='Advective Gain Control';
        labs(3).l='Fall Speed Loss CCN 960';
        labs(4).l='Fall Speed Loss Control';


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(34)+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];


    case 'cumatHtotw'

        %      logflag=2;
        %      lor=4;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        H=16.7;
        %H=16.25;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution to Total Water (ppmv)';

        titlenam=['Contributions to Total Water at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:2
            xdat(idat).x = time;
        end

        ydat(1).y = ad(ih,:);
        ydat(2).y = -cumsum(fallrate(ih,:),2)*300;

        labs(1).l='Advective Gain';
        labs(2).l='Fall Speed Loss';


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(34)+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

    case 'cumatHvapflux'

        %      logflag=2;
        %      lor=4;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        % H=16.7;
        %H=16.25;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution to Vapour (ppmv)';

        titlenam=['Contributions to Vapour at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Advective Source of vapour at H' direc(idir).dir];

        idirs=[1:3]
        for idat=1:length(idirs)
            idir=idirs(idat);

            xdat(idat).x = GridDan(idir).t(dumprange)+3;
            ydat(idat).y = cumsum( sum(icediagsALL(idir).i(ih,dumprange,[1 10]),3) * 300 ) / npess2(idir) ;
            labs(idat).l = runName(idat).nam;


        end

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %         xlimits=[GridDan(idir).t(34)+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];




    case 'cumatHfallcomp'

        %      logflag=2;
        lor=4;

        ad_calcs4timeseries;

        H=15.5;
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Rate of Total Water Change (ppmv s^{-1})';

        titlenam=['Sensitivity to Graupel Fall Speed Scheme at ' num2str(H,4) ' km'];
        figname=['Timeseries of Fall Speed Source at H' direc(idir).dir];

        Lfall=length(fall_from_mean);
        for idat=1:Lfall
            xdat(idat).x = time;
            %ydat(idat).y = -cumsum(fall_from_mean(idat).i(ih,t1:t2),2)*300;
            ydat(idat).y = fall_from_mean(idat).i(ih,:);
            labs(idat).l=['Scheme ' num2str(idat)];
        end




        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %xlimits=[GridDan(idir).t(34)+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];



    case 'cumatHtoti'

        %      logflag=2;
        %      lor=4;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Total Ice at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:2
            xdat(idat).x = time;
        end

        ydat(1).y = iceadcum(ih,:);
        ydat(2).y = -cumsum(fallrate(ih,:),2)*300;

        labs(1).l='Advective Gain';
        labs(2).l='Fall Speed Loss';


        yys=[1 2];


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(dumprange(33))+3 GridDan(idir).t(dumprange(end))+3];
        xlimits=[GridDan(idir).t(dumprange(28))+3 GridDan(idir).t(dumprange(end))+3];


        izlim=0;
        zmin=215
        zmax=220;
        %
        maxx=-1e99;
        lor=0;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        %         ylims=[2 maxx];



    case 'cumatHmicro'

        %      logflag=2;
        %      lor=4;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Total Ice at ' num2str(H) ' km'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(1).y = cumsum(microicerate(ih,:),2)*300;
        %     ydat(2).y = -cumsum(fallrate(ih,:),2)*300;

        labs(1).l='Microphysics';
        %   labs(2).l='Fall Speed Loss';

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        %         ylims=[2 maxx];

        lor=2;


    case 'cumatHvap'

        %      logflag=2;
        %      lor=4;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        H=GridDan(idir).Z(ih+izmin-1)/1000+0.62;

        ylab='Cumulative Contribution (ppmv)';

        titlenam=['Contributions to Vapour at ' num2str(H,4) ' km'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:2
            xdat(idat).x = time;
        end

        ydat(1).y = cumsum(microicerate(ih,:),2)*300; %microicerate is gain of ice or loss of vapour
        ydat(2).y = -vapadcum(ih,:); %vapadcum is the cumlative advective loss of vapour

        labs(1).l='Microphysical Loss';
        labs(2).l='Advective Gain';

        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([max(ydat(idat).y) maxx]);
        end

        nmark=0;

        %         ylims=[2 maxx];

        lor=2;

    case 'maxW'

        xlab=['Time (hrs)'];
        ylab='Maximum Updraught (m s^{-1})';

        titlenam=['Timeseries'];
        figname=['Timeseries of Max W' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = SerDan(idir).SER(:,1)/3600;
        end

        ydat(1).y = SerDan(idir).SER(:,4);

        labs(1).l='Miles City';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([max(ydat(idat).y) maxx]);
        end

        nmark=0;
        lor=1;
        lwidth=2;

    case 'citop'

        xlab=['Time (hrs)'];
        ylab='Max Cloud Top Height (km)';

        titlenam=['Timeseries'];
        figname=['Timeseries of Max CITOP' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = SerDan(idir).SER(:,1)/3600;
        end

        ydat(1).y = SerDan(idir).SER(:,17)/1000;

        labs(1).l='Miles City';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([max(ydat(idat).y) maxx]);
        end

        nmark=0;
        lor=1;
        lwidth=2;

    case 'clbase'

        xlab=['Time (hrs)'];
        ylab='Min Liquid Cloud Base Height (km)';

        titlenam=['Timeseries'];
        figname=['Timeseries of Max CLBAS' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = SerDan(idir).SER(:,1)/3600;
        end

        ydat(1).y = SerDan(idir).SER(:,14)/1000;

        labs(1).l='Miles City';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([max(ydat(idat).y) maxx]);
        end

        nmark=0;
        lor=1;
        lwidth=2;

    case 'cltop'

        xlab=['Time (hrs)'];
        ylab='Max Liquid Cloud Top Height (km)';

        titlenam=['Timeseries'];
        figname=['Timeseries of Max CLTOP' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = SerDan(idir).SER(:,1)/3600;
        end

        ydat(1).y = SerDan(idir).SER(:,13)/1000;

        labs(1).l='Miles City';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([max(ydat(idat).y) maxx]);
        end

        nmark=0;
        lor=1;
        lwidth=2;

    case 'cibas'

        xlab=['Time (hrs)'];
        ylab='Min Ice Cloud Base Height (km)';

        titlenam=['Timeseries'];
        figname=['Timeseries of Max CIBAS' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = SerDan(idir).SER(:,1)/3600;
        end

        ydat(1).y = SerDan(idir).SER(:,18)/1000;

        labs(1).l='Miles City';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([max(ydat(idat).y) maxx]);
        end

        nmark=0;
        lor=1;
        lwidth=2;

    case 'siatH'

        %      logflag=2;
        %      lor=4;

        %ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H);

        ylab='Max Supersaturation wrt Ice (%)';

        titlenam=['Contributions for ' num2str(H) ' km for Total Water'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(1).y = simaxTimH(idir).s(ih,:);

        labs(1).l='Max Supersaturation wrt Ice';


        xlims=1;
        xlimits=[GridDan(idir).t(dumprange(1))+3 GridDan(idir).t(dumprange(end))+3];
        %
        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=0;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

        % nmark=-1;

    case 'maxwatH'

        %ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H);

        ylab='Max Updraught (m^s{-1})';

        titlenam=['Contributions for ' num2str(H) ' km for Total Water'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(1).y = maxW(idir).w(ih,:);
        labs(1).l='Max Updraught Speed (m/s)';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=-1;
        lwidth=1;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

        % nmark=-1;

    case 'upfluxatH'

        %ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H);

        ylab='Mean Upwards Flux (kg m^{-2} s^{-1})';

        titlenam=['Contributions for ' num2str(H) ' km for Total Water'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(idat).y=icediagsALL(idir).i(ih,dumprange,137).*GridDan(idir).RHON(ih);
        labs(1).l='Mean Upwards Flux (kg m^{-2} s^{-1})';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=-1;
        lwidth=1;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

        % nmark=-1;

    case 'meaniceatH'

        logflag=2;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        ylab='Mean Ice Mixing Ratio (kg kg^{-1}))';

        titlenam=['Contributions for ' num2str(H) ' km for Total Water'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(idat).y=changeice(ih,dumprange);
        labs(1).l='Mean Ice Mixing Ratio (kg kg^{-1}))';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=-1;
        lwidth=1;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

        % nmark=-1;

    case 'changevap'

        logflag=0;

        ad_calcs4timeseries;

        H=16.25; %height above which to average
        ih=findheight(GridDan(idir).Z/1000+0.62,H)-izmin+1;

        ylab='Change in Vapour Mixing Ratio (kg kg^{-1}))';

        titlenam=['Contributions for ' num2str(H) ' km for Total Water'];
        figname=['Timeseries of Cumlative Sources at H' direc(idir).dir];

        for idat=1:1
            xdat(idat).x = time;
        end

        ydat(idat).y=changevap(ih,dumprange);
        labs(1).l='Change in Vapour Mixing Ratio (kg kg^{-1}))';

        maxx=-1e99;

        for idat=1:length(ydat)
            maxx=max([ydat(idat).y maxx]);
        end

        nmark=-1;
        lwidth=1;

        izlim=0;
        zmin=215;
        zmax=220;

        savename=[titlenam];

        % nmark=-1;


end

if exist('override_binned_options') &  override_binned_options==1
    clear override_binned_options
end

clear ioverride_lwc_cutoffs



% savename=[titlenam ' ' ylab];