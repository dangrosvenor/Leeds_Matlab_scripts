itimser='cumatHtotiN_resdiff';
itimser='cumatHtotiM_resdiff';
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
itimser = 'echo_top'; %echotop timeseries
%itimser = 'max_supersat';
%itimser = 'max_MR_hslice';
%itimser = 'max_MR_totprc';
%itimser = 'mean_vapour';
%itimser = 'latent_heating';
%itimser = 'av_up'; %check out cloud top too from profiles
itimser = 'max_hm_gm3';
%itimser = 'max_mac3';


multi={'ice_dep_rate','ice_mass','ice_num'};

if subplotting==1
    itimser=multi{iplot};
end


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
xlab=['Time (hrs)'];

    
ititle=1;
    
switch itimser
    %%%%   N.B. - don't forget that for 2D runs you need to divide icediagsALL by npess2(idir) - don't need to worry for 3D
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 
 savename=[ylab ' ' titlenam];