%24

updown=35; %microphysical source rates
%updown=24; %advection, fall speed, microphysics and change for various HMs
% %updown=25; %advection, microphysics and change for vapour
% %updown=18; %change in ice for different runs
% updown=38; %mean HM contents over time
% updown=382; %mean HM contents over time
% %updown=388; %max totwater contents over time
% %updown=42; %comparison between case of individual process rates 
%updown=43; %comparison between cases of advective vapour flux
%updown=44; %profile of mode mass at end of sim
%updown=45; %latent heat release
%updown=46; %mean updraught


idirlabel=0; %flag to label with directory name

    logflag=0;
    ititle=1;
    
runs='ccn';
switch runs
  case 'ccn'
    dire(1).dir='Control';
    dire(2).dir='CCN = 960 cm^{-3}';
    dire(3).dir='IN x 10';
    dire(4).dir='500m res';
case 'res'
    dire(1).dir='250m res';
    dire(2).dir='500m res';
    dire(3).dir='1km res';
    dire(4).dir='250mres 13.4TH';
end
    
hrange=12;
switch hrange
    case 1
        zmin=0.01;
        zmax=18;
    case 2
        zmin=0.01;
        zmax=20;
    case 3
        zmin=6;
        zmax=12;
    case 4
        zmin=15;
        zmax=19;
    case 5
        zmin=14.8;
        zmax=19;
    case 6
        zmin=13;
        zmax=19;    
    case 7
        zmin=0.1;
        zmin=1.2;
        zmax=9;     
     case 8
        zmin=1.2;
        zmax=20.5;    
	case 9
        zmin=16.5;
        zmax=20;  
	case 10
        zmin=0.01;
        zmax=15;        
end

 times=12;
    switch times
    case 1
        time1=19.75; %time (UTC)
        time2=22.4869;
    case 2
		time1=19.75; %time (UTC)
        time2=20.5;
    case 3
        time1=19.75;
        time2=23.75;
    case 4
        time1=19.75;
        time2=23.46;
    case 5
        time1=19.75;
        time2=21.75;
    case 6
        time1=19.75;
        time2=24.6667;    
    case 7
        time1=19.75;
        time2=25.6667;        
    case 8
        time1=19.75;
        time2=20.25;       
    case 9
        time1=19.75;
        time2=21.25;  
    case 10
        time1=19.75;
        time2=23.4167;   
    case 11
        time1=19.75;
        time2=25.1667;  
    case 12
        time1=0;
        time2=2; 
    case 13
        time1=1.2;
        time2=3.5;
    case 14
        time1=0 + 19.75;
        time2=2.67 + 19.75 ;
    case 15
        time1=0;
        time2=2.67 ;    %time used for simulation end in paper (should be 3:35 ideally - need to get extra diags)
    end
    
        mins=(time1-floor(time1))*60;
        minstr=num2str(mins,'%2.0f');
    
        hrs=mod(floor(time1),24);
        hrstr=num2str(hrs,'%2.0f');
        if mins==0; minstr='00';end
        if hrs==0; hrstr='00';end
        t1str=[hrstr ':' minstr];
        
        mins=(time2-floor(time2))*60;
        minstr=num2str(mins,'%2.0f');
    
        hrs=mod(floor(time2),24);
        hrstr=num2str(hrs,'%2.0f');
        if mins==0; minstr='00';end
        if hrs==0; hrstr='00';end
        t2str=[hrstr ':' minstr];

    
%    t1str=num2str(mod(time1,24),'%2.2f');
%    t2str=num2str(mod(time2,24),'%2.2f');

    
    
%     
%     xlims=0;
%     %xlimits=[-0.67 0.63];
%     xlimits=[-0.5 1.2];
%     xlimits=[-500 50];
%     
%     zmin=13;  %(km)
%     %zmin=1;
%     %zmax=22;   %z(end)/1000;
%     
%     %zmin=0.01;
%     %zmax=17.5;
%     zmax=22;
%     
%     zmin=16;
%     zmax=19;
    
    dz=GridDan(idir).Z(2:end)-GridDan(idir).Z(1:end-1); 
    rho=GridDan(idir).RHO;
    
   
    
    t1=findheight(GridDan(1).t,time1);
    t2=findheight(GridDan(1).t,time2);
    
    if t2>size(icediagsALL(1).i,2)
        disp('******* Warning: t2 is larger than time array in icediags(1).i *********');
    end
    
    iz=2;
    iz2=length(GridDan(idir).Z);
        
    z=GridDan(idir).Z(iz+1:iz2)+620;    
    secyA=z/1000;
    secyB=GridDan(idir).PREFN(iz+1:iz2)/100;
        
        switch updown                   
            
        case 46    %profiles of advective flux of vap_neg, etc.

			
            idircont=1;
            
            lor=1; %3=bottom left, 4= bottom right, 1=top right, 2=top left 0=auto, -1= on right
            
            xlimits=[-50 50];
            xlims=0;

%            zmin=15;
%            zmax=19.5;
            izlim=1;                                                                               
            
            
          %  domfact=1;
	            
            maxval=0;  
            
            ind=137; %ALu_W
            aind=280; %ALu_A
            
            ind=302; %ACu_W
            aind=285; %ACu_A
            
             ind=303; %W>1_W
             aind=283; %W>1_A
            


            
            cf=0;
            idirs=[1 2];
            for idat=1:length(idirs)
                idir=idirs(idat);   
                
                area=icediagsALL(idir).i(iz+1:iz2,t1:t2,aind);
                area(area==0)=1e99;                           
                
                if cf==1
                    areacont=icediagsALL(idircont).i(iz+1:iz2,t1:t2,aind);
                    areacont(areacont==0)=1e99;
                    cont=mean(icediagsALL(idircont).i(iz+1:iz2,t1:t2,[ind])./areacont,2);
                    xdat(idat).x= 100*( mean(icediagsALL(idir).i(iz+1:iz2,t1:t2,[ind])./area,2) - cont ) ./cont ;
                    xlab=['Increase in mean updraught (%)'];
                else
                    xdat(idat).x= mean(icediagsALL(idir).i(iz+1:iz2,t1:t2,[ind])./area,2);
                    
                    xlab=['Mean updraught (m s^{-1})'];
                end                    
                     
                dy=diff(GridDan(idir).Y1(1:2)); %multiply by domain size (in km) and covert to grams so is g/kg km
                labs(idat).l=runName(idir).nam; 
                ydat(idat).y=(GridDan(idir).Z(iz+1:iz2)+620) / 1000;

                titlenam=['Average profile ' num2str(mod(time1-19.75,24),4) ' to ' num2str(mod(time2-19.75,24),4) ' hours'];                                            
                
				figname=[xlab ' ' titlenam]  %check that multiplied by rho
                savename=figname;
                
                

            end
            
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        case 45    

			
            idircont=1;
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT - make sure to switch off smoothing if not required!! %%%%%           
            ismooth=1;
			windowSize = 5; %size of the smoothing window
			smooth_nums=[1];
      %%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT - make sure to switch off smoothing if not required!! %%%%%


            
            lor=1; %3=bottom left, 4= bottom right, 1=top right, 2=top left 0=auto, -1= on right
            
            xlimits=[-50 50];
            xlims=0;

%            zmin=15;
%            zmax=19.5;
            izlim=0;                                                                               
            
            
          %  domfact=1;
	            
            maxval=0;
            
            Cp=1005;
            rlvap_on_cp = 2.501e6 / Cp;
            rlsub_on_cp = 2.834e6 / Cp ; %diving by Cp so when multiply by mass get temp change (from E=mL=MCpdT, M=mass of air=1 kg as q in kg/kg)
            
            cf=0;
            idirs=[1];
            for idat=1:length(idirs)
                idir=idirs(idat);    
                dq_liq=mean(sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[29:30]),3),2); %sum of liquid and rain sources, summed over time (using ALL_DQ02 etc.)
                dq_ice=mean(sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[31:33]),3),2); %sum of ice, snow graupel sources, summed over time
                
                if cf==1
                    dq_liq_cont=mean(sum(icediagsALL(idircont).i(iz+1:iz2,t1:t2,[29:30]),3),2); %sum of liquid and rain sources, summed over time (using ALL_DQ02 etc.)
                    dq_ice_cont=mean(sum(icediagsALL(idircont).i(iz+1:iz2,t1:t2,[31:33]),3),2); %sum of ice, snow graupel sources, summed over time
                    cont=rlvap_on_cp*(dq_liq_cont) + rlsub_on_cp*(dq_ice_cont);
                    xdat(idat).x= (rlvap_on_cp*(dq_liq) + rlsub_on_cp*(dq_ice) - cont ); %*100./cont ;
                    
                    xlab=['Increase in latent heating rate (%)'];
                    xlab=['Increase in latent heating rate (K s^{-1})'];
                    
                else
                    xdat(idat).x= rlvap_on_cp*(dq_liq) + rlsub_on_cp*(dq_ice);                    
                    xlab=['Latent heating rate (K s^{-1})'];
                end
                   %temp change per sec due to latent heat (K/s)                     
                     
                dy=diff(GridDan(idir).Y1(1:2)); %multiply by domain size (in km) and covert to grams so is g/kg km
                labs(idat).l=runName(idir).nam; 
                ydat(idat).y=(GridDan(idir).Z(iz+1:iz2)+620) / 1000;

                titlenam=['Mean profile for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' hours'];                                            
                
				figname=[xlab ' ' titlenam]  %check that multiplied by rho
                savename=figname;
                
                
                if ismooth==1
                    for ism=1:length(smooth_nums)
                        xdat(smooth_nums(ism)).x=filter(ones(1,windowSize)/windowSize,1,xdat(smooth_nums(ism)).x); %smooth the data
                    end
                end
                
                

            end
            
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            
        case 44	
            
           % secyA=GridDan(1).Z;
%		    secyB=GridDan(1).PREFN;
            
            ix=45; %time index of 45
            
            lor=1; %3=bottom left, 4= bottom right, 1=top right, 2=top left 0=auto, -1= on right
            
            xlimits=[-500 2500];
            xlims=0;

            zmin=15;
            zmax=18;
            
            
            
            dual=1;
            
            
            titlenam=['Microphysical Process Rates for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' UTC'];
            titlenam='Mass mode diameter profile';
            
            
            xlab=['Mode diameter of total ice mass distribution (microns)']; 

            
            figname=[xlab ' ' titlenam]
            savename=figname;
	            
            maxval=0;
            
            idirs=[1 2];
            for idat=1:length(idirs)
                idir=idirs(idat);   
                
%                gam_dist_vs_height
                
                xdat(idat).x=massMODE(idir).m';
                xdat(idat).x(139:end)=0;
                ydat(idat).y=(GridDan(idir).Z+620)/1000;
                labs(idat).l=runName(idir).nam; 
            end

            
            
        case 43     %profiles of advective flux of vap_neg, etc.

			
            
            lor=1; %3=bottom left, 4= bottom right, 1=top right, 2=top left 0=auto, -1= on right
            
            xlimits=[-500 2500];
            xlims=0;

            zmin=15;
            zmax=19.5;
            
            idirs=[1 2]; %indices of cases to plot (as in direcDan)
            
            adcase='vapneg_tot';
            adcase='vappos_tot';            
           % adcase='micro';
            
            switch adcase
            case 'vapneg_tot'
                titlenam=['Neg Advective Flux of vapour for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' UTC'];
                
	        case 'vappos_tot'
                titlenam=['Pos Advective Flux of vapour for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' UTC'];
	                                                   
            case 'micro'
                titlenam=['Microphyscical source of vapour for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' UTC'];
                
            end
            
            
				figname=[xlab ' ' titlenam]  %check that multiplied by rho
                savename=figname;
                xlab=['Flux (kg m^{-2} s^{-1})']; 
            
            
            
          %  domfact=1;
	            
            maxval=0;
            
            for idat=1:length(idirs)
                idir=idirs(idat);    
                
                dy=diff(GridDan(idir).Y1(1:2)); %multiply by domain size (in km) and covert to grams so is g/kg km
                labs(idat).l=runName(idir).nam; 
                ydat(idat).y=(GridDan(idir).Z+620) / 1000;

                
                switch adcase
                case 'vapneg_tot'
                
                    xdat(idat).x= dy * mean( vapneg_up(idir).dat(:,t1:t2) + vapneg_d(idir).dat(:,t1:t2) ...
                    , 2 ) ; %microphys - sums rate over the times given and multiplies by 300s
                
                case 'vappos_tot'
                
                    xdat(idat).x= dy * mean( vappos_up(idir).dat(:,t1:t2) + vappos_d(idir).dat(:,t1:t2) ...
                    , 2 ) ; %microphys - sums rate over the times given and multiplies by 300s
                
              %      xdat(idat).x= dy * mean( vappos_up(idir).dat(:,t1:t2) ...
              %      , 2 ) ; %microphys - sums rate over the times given and multiplies by 300s

                
                case 'micro'
% 					xdat(idat).x= f*mean( cumsum(  flipud(sum(icediag(idir).i(:,t1:t2,[31 9 1]),3)) , 1) ,2) / npess2(idir) ;   
%                     xdat(idat).x = flipud(xdat(idat).x);
% 
%                     labs(idat).l=['PIDEP' runName(idir).nam];
%                     
%                     xdat(idat+2).x= f*mean( cumsum(  flipud(sum(icediag(idir).i(:,t1:t2,[24 25 27]),3)) , 1) ,2) / npess2(idir) ;   
%                     xdat(idat+2).x = flipud(xdat(idat+2).x);
%                     
%                     labs(idat+2).l=['PISUB' runName(idir).nam];                    
%                     ydat(idat+2).y=(GridDan(idir).Z+620) / 1000;

                    
                    xdat(idat).x= f*mean( cumsum(  flipud(sum(icediagsALL(idir).i(:,t1:t2,[28]),3)) , 1) ,2) / npess2(idir) ;   
                    xdat(idat).x = flipud(xdat(idat).x);

                    labs(idat).l=['ALL_DQ01 ' runName(idir).nam];
                    
                end

            end
            
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        case 42                     %%%%%   individual microphysical process rates for case comparisons %%%%%%

			
            
            lor=1; %3=bottom left, 4= bottom right, 1=top right, 2=top left 0=auto, -1= on right
            
            xlimits=[-500 2500];
            xlims=0;

            zmin=1.1;
            zmax=18.5;
            
            idirs=[1 2]; %indices of cases to plot (as in direcDan)
            
            set_dgs_numrates; %defines dgs values for numrates option - may need to change, see icediags_5thSept_2005_32
            
            pnames={'PIHAL','PIPRM'};
            pnames={'PSACW','PGACW'};
            pnames={'PIFRW'};

            pnam='';
            for i=1:length(pnames)
                pnam=[pnam ' ' pnames{i}];
            end
            iproc=get_prnum(pnames);
  
            iproc_nc=39; %index of first number process rate
            
            titlenam=['Microphysical Process Rates for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' UTC'];
            if min(iproc)<iproc_nc
                xlab=['Mean Microphysical Mixing Ratio Source Rate (kg kg^{-1} s^{-1})']; 
                xlab=['Mixing Ratio Source for ' pnam ' (kg kg^{-1} s^{-1})']; 
                fact=f;
            else
				xlab=['Number Source for ' pnames ' (kg^{-1})']; 
                fact=1;
            end
            
            figname=[xlab ' ' titlenam]
            savename=figname;
            
            
            domfact=length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/fact; %multiply by domain size (in km) and covert to grams so is g/kg km
            
          %  domfact=1;
	            
            maxval=0;
            
            for idat=1:length(idirs)
                idir=idirs(idat);    
                
                xdat(idat).x=domfact*fact*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,iproc),3),GridDan(idir).t,iz+1,iz2)/npess2(idir); %microphys - sums rate over the times given and multiplies by 300s
                ydat(idat).y=z/1000;
                labs(idat).l=runName(idat).nam; 
            end
             

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               

        case 1
        figname=['Vap Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        
        xdat(1).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
        xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
        %xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2)*300;
        xdat(3).x=xdat(1).x(1:end)-xdat(2).x(1:end);
        
        
        xlab='Source of Vapour (ppmv)';
        titlenam=figname;
        
             ydat(1).y=z/1000;
            labs(1).l='Change';
            
            ydat(2).y=z/1000;
            labs(2).l='Microphysics';
            
            ydat(3).y=z/1000;
            labs(3).l='Advection';
            
            xlims=0;
            xlimits=[-2 1.2];
            
            zmin=16;
            zmax=19;
            
            
        
        case 2
        figname=['Total Water Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

                    
         xdat(1).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37:42]),3) ); %difference in end and start vapour diags
       % xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
        xdat(2).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2);
        xdat(3).x=xdat(1).x(1:end)-xdat(2).x(1:end);
        
        
        xlab='Source of Total Water (ppmv)';
        titlenam=figname;
        
             ydat(1).y=z/1000;
            labs(1).l='Change';
            
            ydat(2).y=z/1000;
            labs(2).l='Fall Speed Flux';
            
            ydat(3).y=z/1000;
            labs(3).l='Advection';
            
            xlims=1;
            xlimits=[-2 1.2];
            
            zmax=22;
            
            lor=3;
            
        
        case 3
        titlenam=['Ice Sources for Times ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        
            
         xdat(1).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags
         xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         %xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         
         xdat(3).x=xdat(1).x(1:end)-xdat(2).x(1:end); %fall speed + advection
         
         %xdat(4).x=xdat(1).x(1:end)-xdat(2).x(1:end)-xdat(3).x(1:end); %advection
         
         %xdat(5).x=-sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[4:6 13:15]),3),GridDan(idir).t,iz,iz2),2);
        
        xlab='Source of Ice (ppmv)';        
        savename=[titlenam ' ' xlab];
        figname=savename;
        
            ydat(1).y=z/1000;
            labs(1).l='Change';
            
            ydat(2).y=z/1000;
            labs(2).l='Microphysics';
            
            ydat(3).y=z/1000;
            labs(3).l='Advection + Fall Speed Flux';
            
            %ydat(3).y=z/1000;
            %labs(3).l='Fall Speed Flux';
            
            %ydat(4).y=z/1000;
            %labs(4).l='Advection';
            
            
            
            %xlims=1;
            xlimits=[-2 1.2];
            
            zmax=22;
            zmin=0.01;
            
        case 4
%              case 1
%             idiff=[85:87];
%             imicro=[76:78];
%             ifall=[
         
        figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC in up and down draughts'];  
        titlenam=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Microphysical Source of Ice (ppmv)';

         xdat(1).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,76:78),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,121:123),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         xdat(3).x=xdat(1).x+xdat(2).x;           
            
            ydat(1).y=z/1000;
            labs(1).l='Updraughts';
           
            ydat(2).y=z/1000;
            labs(2).l='Downdraughts';
            
            ydat(3).y=z/1000;
            labs(3).l='Total';
            
            
            case 5
%              case 1
%             idiff=[85:87];
%             imicro=[76:78];
%             ifall=[
         
        figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC in up and down draughts'];  
        titlenam=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Microphysical Source of Ice (ppmv)';

         xdat(1).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,76:78),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,121:123),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         xdat(3).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[85:87]),3) );
         xdat(4).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[130:132]),3) );
         xdat(5).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[67:69]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(6).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[112:114]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(7).x=xdat(3).x(1:end)-xdat(1).x(1:end)-xdat(5).x(1:end);
         xdat(8).x=xdat(4).x(1:end)-xdat(2).x(1:end)-xdat(6).x(1:end);
            
            ydat(1).y=z/1000;
            labs(1).l='Microphysics in Updraughts';
                      
            ydat(2).y=z/1000;
            labs(2).l='Microphysics in Downdraughts';
            
            ydat(3).y=z/1000;
            labs(3).l='Change in Updraughts';
                        
            ydat(4).y=z/1000;
            labs(4).l='Change in Downdraughts';
           
            
            ydat(5).y=z/1000;
            labs(5).l='Fall Speed Flux in Updraughts';
            
            ydat(6).y=z/1000;
            labs(6).l='Fall Speed Flux in Downdraughts';
            
            ydat(7).y=z/1000;
            labs(7).l= 'Advection in Updraughts';
            
            ydat(8).y=z/1000;
            labs(8).l='Advection in Downdraughts';
            
         case 6
%              case 1
%             idiff=[85:87];
%             imicro=[76:78];
%             ifall=[
         
        figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC in up and down draughts'];  
        titlenam=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Source of Ice from Advection (ppmv)';

         upmic=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,76:78),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         downmic=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,121:123),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         upchange=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[85:87]),3) );
         downchange=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[130:132]),3) );
         upfall=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[67:69]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         downfall=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[112:114]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(1).x=upchange-upmic-upfall;
         xdat(2).x=downchange-downmic-downfall;
         xdat(3).x=xdat(1).x+xdat(2).x;
         
            ydat(1).y=z/1000;
            labs(1).l='Updraughts';
            
            ydat(2).y=z/1000;
            labs(2).l='Downdraughts';
            
            ydat(3).y=z/1000;
            labs(3).l='Total';
            
        
		ylab='Height (km)';    
        
        case 7
%              case 1
%             idiff=[85:87];
%             imicro=[76:78];
%             ifall=[
         
        figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC in up and down draughts'];  
        titlenam=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Source of Ice from Fall Speed Flux (ppmv)';

        xdat(1).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[67:69]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
        xdat(2).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[112:114]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
        xdat(3).x=xdat(1).x+xdat(2).x;
        
            ydat(1).y=z/1000;
            labs(1).l='Updraughts';
            
            ydat(2).y=z/1000;
            labs(2).l='Downdraughts';
            
            ydat(3).y=z/1000;
            labs(3).l='Total';
  
       case 8
%              case 1
%             idiff=[85:87];
%             imicro=[76:78];
%             ifall=[
         
        figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC in up and down draughts'];  
        titlenam=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Source of Ice from Microphysics (ppmv)';

        xdat(1).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,[10 16]),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
        xdat(2).x=-f*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,[15 19]),3),GridDan(idir).t,iz+1,iz2); 
        xdat(3).x=xdat(1).x+xdat(2).x;
        
            ydat(1).y=z/1000;
            labs(1).l='Deposition';
            
            ydat(2).y=z/1000;
            labs(2).l='Sublimation';
            
            ydat(3).y=z/1000;
            labs(3).l='Total';
            
            
            
         case 9
              figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];    
                
             xdat(1).x=( sum(icediagsALL(idir).i(iz+1:iz2,t2,[43:45]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[43:45]),3) ); %difference in end and start vapour diags
             xdat(2).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,34:36),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             xdat(3).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[25:27]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
             xdat(4).x=xdat(1).x(1:end)-xdat(2).x(1:end)-xdat(3).x(1:end);
             
             %xdat(5).x=-sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[4:6 13:15]),3),GridDan(idir).t,iz,iz2),2);
            
              xlab='Source of Ice (kg^{-1})';
              titlenam=figname;
            
                ydat(1).y=z/1000;
                labs(1).l='Change';
                
                ydat(2).y=z/1000;
                labs(2).l='Microphysics';
                
                ydat(3).y=z/1000;
                labs(3).l='Fall Speed Flux';
                
                ydat(4).y=z/1000;
                labs(4).l='Advection';
                
                %ydat(5).y=z/1000;
                %labs(5).l='Flux Sources';
            
             case 10
              figname=['Ice Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
              
              MI0=1e-15;
              
              n=sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[43]),3);
		      q=sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[40]),3);
              
                %q(q<1e-10)=1e60;
                n(n==1)=0;
				Rmass_I=n./q;
                Rmass_I(Rmass_I>1/MI0)=1/MI0;
				
             xdat(1).x=300*sum(icediag(idir).i(iz+1:iz2,t1:t2,3).*Rmass_I,2);
             xdat(2).x=300*sum(icediag(idir).i(iz+1:iz2,t1:t2,4).*Rmass_I,2); 
             xdat(3).x=300*sum(icediag(idir).i(iz+1:iz2,t1:t2,5).*Rmass_I,2); 
             xdat(4).x=300*sum(icediag(idir).i(iz+1:iz2,t1:t2,6).*Rmass_I,2);

             
             xdat(5).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,7),3),GridDan(idir).t,iz+1,iz2)./MI0; 
             xdat(6).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,8),3),GridDan(idir).t,iz+1,iz2)./MI0; %microphys - sums rate over the times given and multiplies by 300s
              xdat(7).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,9),3),GridDan(idir).t,iz+1,iz2)./MI0; 
             xdat(8).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,12),3),GridDan(idir).t,iz+1,iz2)./MI0; 
             
             xdat(9).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,13),3),GridDan(idir).t,iz+1,iz2); %RSAUT
             xdat(10).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,14),3),GridDan(idir).t,iz+1,iz2); %RIACI 
             
             xdat(11).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,34),3),GridDan(idir).t,iz+1,iz2); 

             
             %xdat(5).x=-sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[4:6 13:15]),3),GridDan(idir).t,iz,iz2),2);
            
              xlab='Source of Ice (kg^{-1})';
              titlenam=figname;
              
              dgs{1}='PIMLT';
		dgs{2}='PSAUT';
		dgs{3}='PSACI';
		dgs{4}='PRACI_S';
		dgs{5}='PGACI';
		dgs{6}='PRACI_G';
		dgs{7}='PIHAL';
		dgs{8}='PIPRM';
		dgs{9}='PICNT';
		dgs{10}='PIDEP';
		dgs{11}='PIACW';
		dgs{12}='PIFRW';
        dgs{13}='RSAUT';
        dgs{14}='RIACI';
        dgs{15}='PISUB';
        
            
                ydat(1).y=z/1000;
                labs(1).l=dgs{3};
                
                ydat(2).y=z/1000;
                labs(2).l=dgs{4};
                
                ydat(3).y=z/1000;
                labs(3).l=dgs{5};
                
                ydat(4).y=z/1000;
                labs(4).l=dgs{6};
                
                ydat(5).y=z/1000;
                labs(5).l=dgs{7};
                
                ydat(6).y=z/1000;
                labs(6).l=dgs{8};
                
                ydat(7).y=z/1000;
                labs(7).l=dgs{9};
                
                ydat(8).y=z/1000;
                labs(8).l=dgs{12};
                
                ydat(9).y=z/1000;
                labs(9).l=dgs{13};
                
                ydat(10).y=z/1000;
                labs(10).l=dgs{14};
                
                ydat(11).y=z/1000;
                labs(11).l='Total';

                
                
                %ydat(5).y=z/1000;
                %labs(5).l='Flux Sources';
                
        case 11
         
        figname=['Fall speed isg number breakdown ' num2str(time1) ' to ' num2str(time2) ' UTC '];  
        titlenam=['Fall Speed Number Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Fall Speed Number Source (kg^{-1})';

         xdat(1).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[25]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(2).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[27]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(3).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[26]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         
            ydat(1).y=z/1000;
            labs(1).l='Ice';
           
            ydat(2).y=z/1000;
            labs(2).l='Snow';
            
            ydat(3).y=z/1000;
            labs(3).l='Graupel';
                
         
        case 12
            
            ismrsour=[8 9 18 11 36 35];
            ismrsink=[25 14 17 19 6];
    
            figname=['Fall speed isg breakdown ' num2str(time1) ' to ' num2str(time2) ' UTC '];  
            titlenam=['Fall Speed Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
            xlab='Snow Mixing Ratio Process Rate (kg/kg/s)';
              
            for ip=1:length(ismrsour)
                 xdat(ip).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediag(idir).i(:,t1:t2,[ismrsour(ip)]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
                 ydat(ip).y=z/1000;
                 labs(ip).l=dgs{ismrsour(ip)};
            end
            
            for ip2=ip+1:ip+length(ismrsink)
                 xdat(ip2).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediag(idir).i(:,t1:t2,[ismrsink(ip2-ip)]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
                 ydat(ip2).y=z/1000;
                 labs(ip2).l=dgs{ismrsink(ip2-ip)};
            end
      
        case 13
         
        figname=['Fall speed isg mixing ratio breakdown ' num2str(time1) ' to ' num2str(time2) ' UTC '];  
        titlenam=['Fall Speed Mixing Ratio Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Fall Speed Source (ppmv)';

         xdat(1).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(2).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[23]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         
         
            ydat(1).y=z/1000;
            labs(1).l='Ice';
           
            ydat(2).y=z/1000;
            labs(2).l='Snow';
            
            ydat(3).y=z/1000;
            labs(3).l='Graupel';
            
            
        case 14
%              case 1
%             idiff=[85:87];
%             imicro=[76:78];
%             ifall=[
         
        figname=['Ice Number Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC in up and down draughts'];  
        titlenam=['Ice Number Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];  
        
        xlab='Source of Ice Number from Advection (kg^{-1})';

         upmic=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,79:81),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         downmic=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,124:126),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         upchange=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[88:90]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[88:90]),3) );
         downchange=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[133:135]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[133:135]),3) );
         upfall=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[70:72]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         downfall=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[115:117]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         xdat(1).x=upchange-upmic-upfall;
         xdat(2).x=downchange-downmic-downfall;
         xdat(3).x=xdat(1).x+xdat(2).x;
         
            ydat(1).y=z/1000;
            labs(1).l='Updraughts';
            
            ydat(2).y=z/1000;
            labs(2).l='Downdraughts';
            
            ydat(3).y=z/1000;
            labs(3).l='Total';
            
        
		ylab='Height (km)';   
        
        
        case 15
        figname=['Fall Speed Intercomparison for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
            xdat(idir).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2);
            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Fall Speed Source of Total Water (ppmv)';
        titlenam=figname;
        
        lor=2; %put legend in top left
        
        xlimits=[-500 50];
        
        case 16
        figname=['Advection Intercomparison for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
            a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37:42]),3) ); %difference in end and start vapour diags
       % xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
            b=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2);
            xdat(idir).x=a-b;
            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Advection Source of Total Water (ppmv)';
        titlenam=figname;
        
        lor=1; %put legend in top right
        
        xlimits=[-50 500];
        
        case 17
        figname=['Change in Total Water Intercomparison for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
            xdat(idir).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37:42]),3) ); %difference in end and start vapour diags
       
            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Change of Total Water (ppmv)';
        titlenam=figname;
        
        lor=1; %put legend in top right
        
        xlimits=[-0.5 2];
        xlims=1;
        
        
        case 18
                lor=4; %put legend in top right

                
        figname=['Intercomparison for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        
        %snow, graupel, ice - order of MR in q fields
        hm='allice';
        hm='icenum';
        hm='ice';
     %  hm='snow';
     %  hm='graupel';
        
     aind=[];

        switch hm
        case 'allice'
            inds=[40:42];
            inds=[226:228]; %ACC_Q0 4:6
            aind=[284];
            lab='Total Ice Mixing Ratio (ppmv)';
            fact=f;
        case 'ice'
            inds=[42];
            inds=[228];
            aind=[284];
            lab='Ice Mixing Ratio (ppmv)';
            fact=f;   
        case 'snow'
            inds=[40];
            lab='Snow Mixing Ratio (ppmv)';
            fact=f;       
        case 'graupel'
            inds=[41];
            lab='Graupel Mixing Ratio (ppmv)';
            fact=f;    
            lor=1;
        case 'icenum'
            inds=229;
            aind=[284]; %ACC_A            
            lab='Ice Number Concentration (kg^{-1})';
            fact=1;
        end
        
        
        idirs=[1:2];
        for idat=1:length(idirs)          %direcDan)  
            idir=idirs(idat);
            
        if length(aind)==1
            area2=icediagsALL(idir).i(iz+1:iz2,t2,aind)/npess2(idir);
            area2(area2==0)=1;
            
            area=icediagsALL(idir).i(iz+1:iz2,t1,aind)/npess2(idir);
            area(area==0)=1;
        else
            area=1;
            area2=1;
        end
            
            %area=1;
            %area2=1;
            %inds=42;
            
             xdat(idat).x=fact*( sum(icediagsALL(idir).i(iz+1:iz2,t2,inds),3)./area2 - sum(icediagsALL(idir).i(iz+1:iz2,t1,inds),3)./area )...
                 /npess2(idir); %difference in end and start vapour diags
       
            ydat(idat).y=z/1000;
            labs(idat).l=runName(idir).nam;    
        end

        xlab=['Change in ' lab];
        titlenam=figname;
        
        
        xlimits=[-0.1 2];
        xlims=0;
        
        savename=[figname lab];

        
        
        case 19
        figname=['Microphysical Source of Ice Intercomparison for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
%            xdat(idir).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
            xdat(idir).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
   
            ydat(idir).y=z/1000;
            labs(idir).l=runName(idir).nam;    
        end

        xlab='Microphysical Source of Ice Mixing Ratio (ppmv)';
        titlenam=figname;
        
        lor=1; %put legend in top right
        
        xlimits=[-0.5 10];
        xlims=0;
        
        
        case 20
        figname=['Advection Source of Ice Intercomparison for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags
             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
            
            xdat(idir).x=a-b-c;
            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Advection Source of Ice Mixing Ratio (ppmv)';
        titlenam=figname;
        
        lor=1; %put legend in top right
        
        xlimits=[-0.5 500];
        xlims=1;
        
        xlimits=[-0.5 1];
        zmin=16.5;
        
        case 21
        figname=['Change in Vapour Mixing Ratio Intercomparison ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
            xdat(idir).x=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Change in Vapour Mixing Ratio (ppmv)';
        titlenam=figname;
        
        lor=1; %put legend in top right
        
        xlimits=[-0.5 500];
        xlims=1;
        
        xlimits=[-0.5 1.3];
        %zmin=16.5;
        
        case 22
        figname=['Advection Source of Vapour MR Intercomparison ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
            a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
            b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
            %xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2)*300;
            xdat(idir).x=a-b;
            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Advection Source of Vapour Mixing Ratio (ppmv)';
        titlenam=figname;
        
        lor=1; %put legend in top right
        
        xlimits=[-0.5 500];
        xlims=0;
        
        xlimits=[-1 15];
        %zmin=16.5;
        
        case 23 %N.B. - is the same as the ice microphysical gain with opposite sign
        figname=['Microphysical Source of Vapour Mixing Ratio Intercomparison ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        for idir=1:length(direcDan)                    
            xdat(idir).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
            %xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2)*300;

            ydat(idir).y=z/1000;
            labs(idir).l=dire(idir).dir;    
        end

        xlab='Microphysical Source of Vapour Mixing Ratio (ppmv)';
        titlenam=figname;
        
        lor=2; %legend position 1=right 2=left
        
        xlimits=[-0.5 500];
        xlims=0;
        
        xlimits=[-5 0.3];
        %zmin=16.5;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        case 24
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idircont=1;
            idir=2;
            
            casestr=runName(idir).nam;
            
            lor=2; %1=right, 2=left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
            lor=2;
        
            
        xlims=0;
        xlimits=[-1 1]*1e8;
%        xlimits=[-200 150];

        
        cf=1;
        
        
%        xlimits=[-600 700];


        
            
        ihm=[1 2 3]; %set as 1 for ice, 2 for snow, and 3 for graupel or combinations - e.g. [1 2 3] for all 3. 4=liquid, 567=ice,snow, graupel no. conc
      %  ihm=[1];
       % ihm=[5];
   %    ihm=[5 6 7];
   %    ihm=4;
        tithm={'Ice MR','Snow MR','Graupel MR','Liquid MR','Ice NC','Snow NC','Grapel NC'};    
        
        if max(ihm)<=4
            units='(g kg^{-1} km)';
            units='(g kg^{-1})';            
            units='(ppmv)';            

            f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
           % f=1; %1/1000 cancelled by conversion to grams
            if length(ihm)==3
                xtit='total ice mixing ratio ';
            else
                xtit=tithm{ihm};
            end
        else
            units='(kg^{-1})';
%            units='(kg^{-1} km)';

            f=1/1000;%so that is per km
            f=1;
            
            if length(ihm)==3
                xtit='total ice no. conc. ';
            else
                xtit=tithm{ihm};
            end
            
        end
        

        
        
            switch cf
            case 1
                titlenam=[casestr ' cf. Control for ' t1str ' to ' t2str ' UTC'];
                titlenam=[casestr ' cf. Control'];
                
                xlab=['Difference of ' xtit ' Source from Control ' units];
                
                if length(tithm)>1
                    xlab=['Difference of total ice source from control ' units];
                end
                    
            case 0
                titlenam=[casestr ' case up to ' t2str ' UTC '];
                titlenam=['Domain average contributions up to ' t2str ' UTC '];
                titlenam=[runName(idir).nam];
                
                xlab=['Contribution to ' xtit units];
               % xlab=['Contribution to total ice number ' units];

            end
                
            savename=[titlenam xlab];
            figname=savename;
            
%         xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
%         titlenam=figname;
%         
%         figname=[casestr tithm{ihm} ' Sources cf. Control ' num2str(time1) ' to ' num2str(time2) ' UTC'];
%         savename=figname;

        
        ifall=[24 22 23 20 25 27 26]; %q456=snow, gr, ice q789=ice, graupel, snow no. concs
        imic=[33 31 32 29 34 36 35]; %for microphysics
        iq=[42 40 41 38 43 45 44]; %for q values
        
        domfact=f*length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2)); %multiply by domain size (in km) and covert to grams so is g/kg km
        domfactcont=f*length(GridDan(idircont).Y1).*diff(GridDan(idircont).Y1(1:2)); %multiply by domain size (in km) and covert to grams so is g/kg km
        
        domfact=f;
        domfactcont=f;
        
             a=domfact*( sum(icediagsALL(idir).i(iz+1:iz2,t2,iq(ihm)),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,iq(ihm)),3) )/npess2(idir); %difference in end and start values
             acont=domfactcont*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,iq(ihm)),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,iq(ihm)),3) )/npess2(idircont); %difference in end and start values

             b=domfact*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,imic(ihm)),3),GridDan(idir).t,iz+1,iz2)/npess2(idir); %microphys - sums rate over the times given and multiplies by 300s
             bcont=domfactcont*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,imic(ihm)),3),GridDan(idircont).t,iz+1,iz2)/npess2(idircont); %microphys - sums rate over the times given and multiplies by 300s
                          
             c=domfact*sum(TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,ifall(ihm)),3),GridDan(idir).t,iz,iz2),2)/npess2(idir); %fall speed flux source
             ccont=domfactcont*sum(TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,ifall(ihm)),3),GridDan(idircont).t,iz,iz2),2)/npess2(idircont); %fall speed flux source
             
             ihh=findheight(GridDan(idir).Z+620,16e3)-iz;             
             ihh2=findheight(GridDan(idir).Z+620,17e3)-iz;
             
             meanchange=mean(a(ihh:ihh2));
             
          
             if ihm==4 %for liquid case where no fall speed flux
                 c=0;
                 ccont=0;
             end
             
             d=a-b-c;
             dcont=acont-bcont-ccont;


    
    if cf==1
        xdat(1).x=a-acont;
        xdat(2).x=b-bcont;
        xdat(3).x=c-ccont;
        xdat(4).x=- (-d+dcont);
    else
        xdat(1).x=a;
        xdat(2).x=b;
        xdat(3).x=c;
        xdat(4).x=d;
    end

        
        ydat(1).y=z/1000;
        labs(1).l='Change';
        
        ydat(2).y=z/1000;
        labs(2).l='Microphysics';
        
        ydat(3).y=z/1000;
        labs(3).l='Fall Speed';
        
        ydat(4).y=z/1000;
        labs(4).l='Advection';
        
        if ihm==4
            xdat(3)=xdat(4);
            ydat(3)=ydat(4);
            labs(3)=labs(4);
            xdat(4)=[];
            ydat(4)=[];
            labs(4)=[];
        end
        
            ih2=findheight(GridDan(1).Z+620,17.1e3);
		%    ih1=findheight(GridDan(1).Z+620,15.335e3);
            ih1=findheight(GridDan(1).Z+620,16.1e3);
    
    
        
        for idat=1:length(xdat)
            %ydat(idat).y = GridDan(idir).Z/1000 + 0.62; 
            
            rho=GridDan(idir).RHON(ih1:ih2); %convert to kg/km3 as xdat in g/kg km 
            dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;
            
            air1=mean( dz.*rho );
   
           % i0=find(xdat(idat).x>6.5);            

            tot(idat).t=mean( xdat(idat).x(ih1:ih2) .*dz.*rho )   ./ air1;
            
             
		end

        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        case 25
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=4;
            casestr=dire(idir).dir;
            
        figname=[casestr ' Vapour Sources Compared With Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        savename=figname;

        idircont=3;
        
        
         a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
         acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
  
         b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,28),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s

         c=a-b;
         ccont=acont-bcont;



        xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
        titlenam=figname;
        
        xdat(1).x=a-acont;
        xdat(2).x=- (b-bcont);
        xdat(3).x=c-ccont;
        

        
        ydat(1).y=z/1000;
        labs(1).l='Change';
        
        ydat(2).y=z/1000;
        labs(2).l='Microphysics (-ve)';
        
        ydat(3).y=z/1000;
        labs(3).l='Advection';
        
        %ydat(4).y=z/1000;
        %labs(4).l='Ice MR';
        
        
        lor=1; %1=right, 2=left
        
        xlimits=[-0.5 500];
        xlims=0;
        
        xlimits=[-0.5 1];
        %zmin=6.2;
        
        
        case 26
        
        %casestr='500m res';
        %casestr='CCN=960 cm^{-3}';
        idir=2;
        casestr=dire(idir).dir;
            
        figname=[casestr ' Ice Sources cf. Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        savename=['updown ' figname];
        
            idircont=1;
        
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[85:87]),3) ); %difference in end and start
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[85:87]),3) ); %difference in end and start vapour diags

             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,76:78),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,76:78),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[67:69]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
             ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,67:69),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
              
             dup=a-b-c;
             dcontup=acont-bcont-ccont;
             
             xdat(1).x=a-acont;
             
             
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[130:132]),3) ); %difference in end and start
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[130:132]),3) ); %difference in end and start vapour diags

             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,121:123),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,121:123),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[112:114]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
             ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[112:114]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
              
             ddown=a-b-c;
             dcontdown=acont-bcont-ccont;
             
             
             xdat(2).x=a-acont;
             
             
%             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags
%             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%             xdat(3).x=a-b-c;
         
         
         


        xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
        titlenam=figname;
        
%         xdat(1).x=dup-dcontup;
%         xdat(2).x=ddown-dcontdown;

        
        
        ydat(1).y=z/1000;
        labs(1).l='Updraught';
        
        ydat(2).y=z/1000;
        labs(2).l='Downdraught';
        
%         ydat(3).y=z/1000;
%         labs(3).l='Total';
        

        
        lor=1; %1=right, 2=left
        
        xlimits=[-0.5 500];
        xlims=0;
        
        xlimits=[-0.5 1];
        zmin=0.01;
        
        
        case 27
        figname=['CCN=960 cm^{-3} Ice MR Compared With Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        savename=figname;

        idir=2;
        idircont=1;
        
        
  
         b=(1/300)*f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,40:42),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         bcont=(1/300)*f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,40:42),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s

         xdat(1).x=b-bcont;



        xlab='Difference of CCN=960 cm^{-3} Case from Control (ppmv)';
        titlenam=figname;
        
        ydat(1).y=z/1000;
        labs(1).l='Ice MR';
        
        
        lor=1; %1=right, 2=left
        
        xlimits=[-0.5 500];
        xlims=0;
        
        xlimits=[-0.5 1];
        zmin=0.1;
        
        case 28
        figname=['Advective Ice Sources in Control Case, ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        savename=['updown ' figname];

        idir=2;
        idircont=1;
        
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[85:87]),3) ); %difference in end and start vapour diags

             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,76:78),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,67:69),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
              
             xdat(6).x=acont-bcont-ccont;
             
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[130:132]),3) ); %difference in end and start vapour diags

             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,121:123),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[112:114]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
             
             
             
             
             xdat(1).x=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[130:132]),3) ); %difference in end and start vapour diags

             xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,121:123),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[112:114]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
             
              
             %xdat(1).x=acont-bcont-ccont;
             
           xdat(4).x=acont-bcont-ccont;
           
           
           a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags
            b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
            c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
            xdat(5).x=a-b-c;


        xlab='Advective Ice Source (ppmv)';
        titlenam=figname;
        
        ydat(1).y=z/1000;
        labs(1).l='Change';
        
        ydat(2).y=z/1000;
        labs(2).l='Microphysics';
        
        ydat(3).y=z/1000;
        labs(3).l='Fall Speed Flux';
        
        ydat(4).y=z/1000;
        labs(4).l='Advection (-ve)';
        
        ydat(5).y=z/1000;
        labs(5).l='Total';
        
        ydat(6).y=z/1000;
        labs(6).l='Up adv';
        
        
%         ydat(1).y=z/1000;
%         labs(1).l='Updraught';
%         
%         ydat(2).y=z/1000;
%         labs(2).l='Downdraught';
        
        xlims=0;
        zmin=0.1;
        
        
        case 29
            
        
            %casestr='500m res';
            %idir=2;
            casestr=dire(idir).dir;
            
        figname=[casestr ' Total Water Sources cf. Control ' num2str(time1,'%2.2') ' to ' num2str(time2,'%2.2') ' UTC'];
        savename=figname;

        idircont=1;
        
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37:42]),3) );%difference in end and start vapour diags
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37:42]),3) );%difference in end and start vapour diags

             b=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2);
             bcont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2);
                          
             d=a-b;
             dcont=acont-bcont  ;


        xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
        titlenam=figname;
        
        xdat(1).x=a-acont;
        xdat(2).x= - (b-bcont);
        xdat(3).x= d-dcont;
        
      
        
            ydat(1).y=z/1000;
            labs(1).l='Change';
            
            ydat(2).y=z/1000;
            labs(2).l='Fall Speed Flux (-ve)';
            
            ydat(3).y=z/1000;
            labs(3).l='Advection';
            
            xlims=0;
            %zmin=0.1;
            
            xlimits=[-2 2];
            
            lor=2; %1=right, 2=left
            
                        
            case 30
            
            casestr='500m res';
            idir=4;

            %casestr='CCN=960 cm^{-3}';
            %idir=2;

            
            ivapsource=[24:27]; %source and sink od vapour indices in icediag
            ivapsink=  [9 1 31 30]; 
            
        
        figname=[casestr ' Microphysical Vapour Sources cf. Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        savename=figname;

        idircont=1;
        

      
             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,28),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             xdat(1).x=b-bcont; %overall difference in microphysics term
             
             %note these use icediag not icediagsALL since are using the process rates for sources and sinks of vapour
             c=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,ivapsource),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             ccont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediag(idircont).i(:,t1:t2,ivapsource),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             xdat(2).x=c-ccont; %overall difference in microphysics term
             
             d=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,ivapsink),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             dcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediag(idircont).i(:,t1:t2,ivapsink),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             xdat(3).x=d-dcont; %overall difference in microphysics term
         
         
         


        xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
        titlenam=figname;
        
      
        
            ydat(1).y=z/1000;
            labs(1).l='Total';
            
            ydat(2).y=z/1000;
            labs(2).l='Sources';
            
            ydat(3).y=z/1000;
            labs(3).l='Sinks';
            
            xlims=0;
            %zmin=0.1;
            
            xlimits=[-5 90];
            
            lor=1; %1=right, 2=left
            
        case 31
            
            
            idir=4;
            dire(idir).dir;
            
            %casestr='CCN=960 cm^{-3}';
            %idir=2;

            
            ivapsource=[24:27]; %source and sink of vapour indices in icediag
            ivapsink=  [9 1 31 30]; 
            
        
        %figname=[casestr ' Mean Updraught cf. Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        figname=['Mean Updraught Intercomparison ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];

        savename=figname;

        idircont=1;
  
        %c=sum(TotMassBudgetALL(GridDan(idir),sum(fluxes(idir).f(:,t1:t2,[2]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
        %ccont=sum(TotMassBudgetALL(GridDan(idircont),sum(fluxes(idircont).f(:,t1:t2,[2]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
        
        c=mean(fluxes(idir).f(iz+1:iz2,t1:t2,2),2);
        ccont=mean(fluxes(idircont).f(iz+1:iz2,t1:t2,2),2);
        
        
%        xdat(1).x=(c-ccont)./ccont * 100;     
%        xdat(1).x=(c-ccont);     

        %xdat(1).x=ccont;
        %xdat(2).x=c;   
        

            

        %xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
        xlab=['Mean Updraught Speed (m/s)'];
        titlenam=figname;

      
%         
%             ydat(1).y=z/1000;
%             labs(1).l='Control';
%             
%             ydat(2).y=z/1000;
%             labs(2).l=casestr;
            
            for idir=1:length(direcDan)                    
                xdat(idir).x=mean(fluxes(idir).f(iz+1:iz2,t1:t2,2),2);
           
                ydat(idir).y=z/1000;
                labs(idir).l=dire(idir).dir;    
            end
        
        
            
            
            xlims=1;
            zmin=0.1;
            zmax=z(end)/1000;
            
            xlimits=[0 0.8];
            
            lor=1; %1=right, 2=left

            
            case 32
            %casestr='CCN=960 cm^{-3}';
           % idir=4;
            casestr=dire(idir).dir;
                
            figname=[casestr ' Water Changes cf. Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
            savename=figname;
	
            idircont=1; %number in direcDan for control run
            
             %change in vapour
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
                
             %change in total water
             b=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37:42]),3) ); %difference 
             bcont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37:42]),3) ); %difference 
             
      
	
	
	
            xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
            titlenam=figname;
            
            xdat(1).x=a-acont;
            xdat(2).x= (b-bcont);
            
        
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags

             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
             ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
              
             d=a-b-c;
             dcont=acont-bcont-ccont;
             
             
             xdat(3).x= - (b-bcont); %difference in microphysics
             xdat(4).x=(d-dcont) + (c-ccont) ; 
             
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37]),3) ); %difference in end and start vapour diags
      
             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,28),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,28),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
	
             c=a-b;
             ccont=acont-bcont;

             xdat(5).x= (c-ccont) 
             
   
            
	
            
            ydat(1).y=z/1000;
            labs(1).l='Vapour';
            
            ydat(2).y=z/1000;
            labs(2).l='Total Water';
            
            ydat(3).y=z/1000;
            labs(3).l='Mphys (vap)';
            
            ydat(4).y=z/1000;
            labs(4).l='Adv + Fall (ice)';
            
            ydat(5).y=z/1000;
            labs(5).l='Adv (vapour)';
            
            
            

            
            %ydat(4).y=z/1000;
            %labs(4).l='Ice MR';
            
            
            lor=1; %1=right, 2=left
            
            %xlimits=[-2 2];
            xlims=1;
            
            xlimits=[-1 1];
            zmin=15;
            
             case 322
            %casestr='CCN=960 cm^{-3}';
           % idir=4;
            casestr=dire(idir).dir;
                
            figname=[casestr ' Water Changes cf. Control ' t1str ' to ' t2str ' UTC'];
            savename=figname;
	
            idircont=1; %number in direcDan for control run
            
             %change in vapour
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37]),3) )/npes; %difference in end and start vapour diags
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37]),3) )/npes; %difference in end and start vapour diags
                
             %change in total water
             b=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[37:42]),3) )/npes; %difference 
             bcont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[37:42]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[37:42]),3) )/npes; %difference 
             
             %change in ice
             c=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[39:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[39:42]),3) )/npes; %difference 
             ccont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[39:42]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[39:42]),3) )/npes; %difference 
      
	
	
	
            xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
            titlenam=figname;
            
            xdat(1).x=a-acont;
            xdat(2).x= (b-bcont);
            xdat(3).x= (c-ccont);
            
            
            ydat(1).y=z/1000;
            labs(1).l='Vapour Change';
            
            ydat(2).y=z/1000;
            labs(2).l='Total Water Change';
            
            ydat(3).y=z/1000;
            labs(3).l='Ice Change';
            
            
            
            

            
            %ydat(4).y=z/1000;
            %labs(4).l='Ice MR';
            
            
            lor=1; %1=right, 2=left
            
            %xlimits=[-2 2];
            xlims=1;
            
            xlimits=[-0.05 0.18];
            zmin=15;
            
            
            case 33
            figname=['Homogeneous freezing for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        
		xdat(1).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediags(1).i(:,t1:t2,34),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s        
        
        xlab='Source of Ice (ppmv)';
        titlenam=figname;
        
            ydat(1).y=z/1000;
            labs(1).l=dire(2).dir;
            
       case 34
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=2;
            casestr=dire(idir).dir;
            
        figname=['Control : ' num2str(time1) ' to ' num2str(time2) ' UTC'];
        savename=figname;

        idircont=1;
        
             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start values
             acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start values

             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
             bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                          
             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
             ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
              
             d=a-b-c;
             dcont=acont-bcont-ccont;


    xlab=['Ice Sources in Control Case (ppmv)'];
        titlenam=figname;
        
%         xdat(1).x=(a-acont)./acont*100;
%         xdat(2).x=(b-bcont)./bcont*100;
%         xdat(3).x=(c-ccont)./ccont*100;
%         xdat(4).x=(-d+dcont)./dcont*100;

        xdat(1).x=acont;
        xdat(2).x=bcont;
        xdat(3).x=ccont;
        xdat(4).x=dcont;

        
        ydat(1).y=z/1000;
        labs(1).l='Change';
        
         ydat(2).y=z/1000;
         labs(2).l='Microphysics';
        
        ydat(3).y=z/1000;
        labs(3).l='Fall Speed Flux';
        
        ydat(4).y=z/1000;
        labs(4).l='Advection (-ve)';
        
        lor=1; %1=right, 2=left 0=auto
        
        xlimits=[-0.5 0.75]*1e4;
        xlims=1;
        
        %xlimits=[-600 700];
        zmin=0.01;
        
        case 35                     %%%%%   microphysical process rates for ice, snow or graupel %%%%%%

            cf=0;  %set to say whether a difference run (=1) or not (=0)
            
            idircont=1; %index in direcDan for control run *******
            idir=1; %%% NOTE - when cf=1 make sure idir and idircont are different!!!!
            
            overall=0; %whether to plot overall sum of process rates

            
 %%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT - make sure to switch off smoothing if not required!! %%%%%           
 
            ismooth=1; %0=switch off all smoothing
			windowSize = 20; %size of the smoothing window (5)
			smooth_nums=[1]; %vector of the lines to be smoothed
            smooth_heights=[10]; %heights (km) below which to apply smoothing - for domain top enter a height above domain top
            
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT - make sure to switch off smoothing if not required!! %%%%%
            
            
            %casestr='500m res';
            casestr='CCN=960 cm^{-3}';
%            casestr=dire(idir).dir;

			casestr=runName(idir).nam;
            
            lor=2; %3=bottom left, 4= bottom right, 1=top right, 2=top left 0=auto, -1= on right
            
            xlimits=[-500 2500];
            xlims=0;
            
            %xlimits=[-600 700];
            %zmin=9;
            zmin=1.1;
            zmax=18.5;
       %     zmax=20.5;
%            zmax=9;
        
            
 
            
            
           % casestr=direcDan(idir).dir;
            
           vapsour=[24:27];
           vapsink=[999 1 9 30 31];
           %vapsink=999;
           %vapsink=[];
           
           imrsour=[];
           imrsink=[];
           
            imrsour=[29:34]; %sources of ice mixing ratio 
            imrsink=[8 11 12 27 7 21 36];
            
        %    imrsour=[29 30]; %PIHAL, PIPRM
      %       imrsour=[29 30 33]; %34 = PIFRW
          %   imrsink=[11 12 36 21];
          %   imrsour=[];
           %  imrsink=[45 46]; %RSAUT RIACI
%             imrsink=imrncsink;

       %     imrsour=[20]; %PIHAL, PIPRM
            
            smrsour=[8 9 18 11 35 36]; %sources of snow mixing ratio 
            smrsink=[25 14 17 19 6];
            
            gmrsour=[1 12:17 10 19:21]; %sources of graupel mixing ratio 
            gmrsink=[2 24 4];
            
            allsour=[1 13 15 20 16 10 imrsour 9 18 35]; %sources of all ice species - i.e. not including ones that convert from one ice to another
            allsink=[2 4 6 7 24 25 27]; %as above but sinks
            
            liqsink=[13 3 18 5 32 29 33 34]; 
            
            icencsour=[59:62]; %62=PIFRW
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
            
            
%%%%%%%%%%%%%%%%%%%%%%%%         Choose the required hydrometeor here         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            hm='Snow';
            hm='Snownc';
            hm='Ice';
%            hm='Icenc';           
            hm='Liquid';
      %      hm='Total ice';
     %       hm='Vapour';
      %      hm='Graupel';
%            hm='Graupelnc';
           
%            hm='Rain';
            
            
            lab=hm;
            switch hm
            case 'Ice'
                  indexsour=[1]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                  indexsink=[2]; %same for sinks
             case 'Icenc'
                  indexsour=[15]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                  indexsink=[16]; %same for sinks      
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
            case 'Graupelnc'
                  indexsour=[imass_nc+5]; %choose index groups (from indexarray) for required sources, e.g. [1 3 5] for ice, snow & graupel MR sources
                  indexsink=[imass_nc+6]; %same for sinks        
              case 'Rain'
                  indexsour=12;
                  indexsink=[13];
                  
            end
        
            
                
            for iind=1:length(indexsour)
                hmssour=[hmssour indexarray(indexsour(iind)).i];
            end
            for iind=1:length(indexsink)
                hmssink=[hmssink indexarray(indexsink(iind)).i];
            end
            

            
            ii=[hmssour hmssink];
            
            
            
            
            
          
            
            if cf==1
                titlenam=[casestr ' cf. Control for ' num2str(mod(time1,24),'%2.2f') ' to ' num2str(mod(time2,24),'%2.2f') ' UTC'];
                titlenam=[casestr ' cf. Control'];
%                xlab=['Difference of Ice Microphysical Source Rate from Control (ppmv)']; 
                if max([indexsour indexsink])<imass_nc
                    xlab=['Diff in Mphys Contibutions to ' lab ' from Control (ppmv)']; 
                    xlab=['Diff in Mphys Contibutions to ' lab ' from Control (g kg^{-1} km)']; 

                    fact=f;
                else
                    xlab=['Diff in Mphys Contibutions to ' lab ' No. from Control (kg^{-1} km)']; 
                    fact=1;
                end

            else %not a difference from control graph
                titlenam=['Microphysical Process Rates for ' num2str(mod(time1,24),4) ' to ' num2str(mod(time2,24),4) ' UTC'];
                if max([indexsour indexsink])<imass_nc
                    xlab=[lab ' Microphysical Source Rate (ppmv)']; 
                    xlab=['Contribution to ' lab ' Mixing Ratio (g kg^{-1} km)']; 
                    
                    fact=f;
                else
                   xlab=[lab ' Microphysical Number Source Rate (kg^{-1})']; 
                   xlab=[lab ' No. Source Rate (kg^{-1} km)']; 
                   
                    fact=1;
                end

                
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
                 if cf==1
                    
                  xdat(i).x= fact*( domfact*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,liqsink),3),GridDan(idir).t,iz+1,iz2)/npess2(idir)... %microphys - sums rate over the times given and multiplies by 300s        
                            - domfactcont*TotMassBudgetProfRate2(GridDan(idircont),sum(icediag(idircont).i(:,t1:t2,liqsink),3),GridDan(idircont).t,iz+1,iz2)/npess2(idircont) )...
                         + ( domfact*fact*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,29),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) - ...
                               domfactcont*fact*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,29),3),GridDan(idircont).t,iz+1,iz2)/npess2(idircont) );
                           %have done all the sinks minus the change due to the microphysics term
                           
                           if hmssink(1)==999
                             %  xdat(i).x = - xdat(i).x;
                           end
                           
                   else
                       xdat(i).x= domfact*fact*( TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,liqsink),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) )... %microphys - sums rate over the times given and multiplies by 300s        
                     + ( domfact*fact*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,29),3),GridDan(idir).t,iz+1,iz2)/npess2(idir) );
                   end
             else
                 
                if cf==1 %if want a comparison (difference from control run)
                    xdat(i).x=fact*(domfact*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,ii(i)),3),GridDan(idir).t,iz+1,iz2)/npess2(idir)... %microphys - sums rate over the times given and multiplies by 300s        
                            - domfactcont*TotMassBudgetProfRate2(GridDan(idircont),sum(icediag(idircont).i(:,t1:t2,ii(i)),3),GridDan(idircont).t,iz+1,iz2)/npess2(idircont) );
                else %if want just an absolute plot of a run
                    xdat(i).x=domfact*fact*TotMassBudgetProfRate2(GridDan(idir),sum(icediag(idir).i(:,t1:t2,ii(i)),3),GridDan(idir).t,iz+1,iz2)/npess2(idir); %microphys - sums rate over the times given and multiplies by 300s
                end
                        %     xdat(i).x=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediag(idircont).i(:,t1:t2,ii(i)),3),GridDan(idircont).t,iz+1,iz2) ;
                        
             end
             
        %      xdat(i).x=xdat(i).x * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/fact; %multiply by domain size (in km) and covert to grams so is g/kg km

                
                ydat(i).y=z/1000;
                
                labs(i).l=upper(dgs{ii(i)});  
                
                
                iund=findstr(labs(i).l,'_'); %replace underscores with a dash as matlab makes text subscript next to them
                if length(iund)>0; labs(i).l(iund)='-'; end;
                
                if i>length(hmssour); 
                    xdat(i).x=-xdat(i).x; 
                    labs(i).l=[labs(i).l '(-)'];
                else %if it is a sink then invert the sign as positive change will reduce the ice hm ir  
                    labs(i).l=[labs(i).l '(+)'];
                end
                %so is showing the changes caused to the ice and not the difference in the sources/sinks
                
                maxval=max([maxval abs(max(xdat(i).x))]);    
            end
	        
            summ=zeros(size(xdat(1).x));
            ic=0;
            for i=1:length(ii)   %remove the rates with very low contributions for clarity
                summ=summ+xdat(i).x;
	            if max(abs(xdat(i).x))>maxval/20;
                    ic=ic+1;
                    xdat(ic).x=xdat(i).x;
                    ydat(ic).y=ydat(i).y;
                    labs(ic).l=labs(i).l;
                end  
            end
                
%                 ic=2;                             
                xdat(ic+1:length(ii))=[];
                ydat(ic+1:length(ii))=[];
                labs(ic+1:length(ii))=[];
%                 
%                 
%                 xdat(1).x=summ;
%                 
                 b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
                 bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                 
%                 xdat(2).x=b-bcont;
	
            
%            xdat(ic+1).x=b-bcont; %add these to put the actual change in the ice hm mr on as well
%            ydat(ic+1).y=z/1000;
%            labs(ic+1).l='Change';

if ismooth==1
    
    if length( ydat ) == 0 
        disp('*************  Error! - ydat has zero length! ***********');
        return
    end
    
    for ism=1:length(smooth_nums)
        ihsmooth=findheight(ydat(smooth_nums(ism)).y,smooth_heights(ism));
        xdat(smooth_nums(ism)).x(1:ihsmooth)=filter(ones(1,windowSize)/windowSize,1,xdat(smooth_nums(ism)).x(1:ihsmooth)); %smooth the data
    end
end

        
switch overall
case 1
%            xdat(ic+1).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,34:36),3),GridDan(idir).t,iz+1,iz2);
        xdat(ic+1).x=0;
        for i=1:ic
            xdat(ic+1).x=xdat(ic+1).x+xdat(i).x;
        end
        labs(ic+1).l='Overall';
        ydat(ic+1).y=z/1000;
end
            
            
            
            
            
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        case 36
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=2;
            casestr=dire(idir).dir;
            
            figname=['Ice Number Change for ' num2str(time1) ' to ' num2str(time2) ' UTC'];
            savename=figname;
	
            idircont=1;
        
            for idir=1:length(dire)
             xdat(idir).x=( sum(icediagsALL(idir).i(iz+1:iz2,t2,[43:45]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[43:45]),3) ); %difference in end and start values
             ydat(idir).y=z/1000;
             labs(idir).l=dire(idir).dir;
            end

%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              d=a-b-c;
%              dcont=acont-bcont-ccont;


            xlab=['Change in Ice Number Concentration (#/kg)'];
            titlenam=figname;
            
	%         xdat(1).x=(a-acont)./acont*100;
	%         xdat(2).x=(b-bcont)./bcont*100;
	%         xdat(3).x=(c-ccont)./ccont*100;
	%         xdat(4).x=(-d+dcont)./dcont*100;
	

%             xdat(3).x=ccont;
%             xdat(4).x=dcont;
	

%             
%             ydat(3).y=z/1000;
%             labs(3).l='Fall Speed Flux';
%             
%             ydat(4).y=z/1000;
%             labs(4).l='Advection (-ve)';
            
            lor=1; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            xlims=1;
            
            %xlimits=[-600 700];
            zmin=14.5;
        
        case 37
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=2;
            casestr=dire(idir).dir;
            
            figname=['Average Ice Particle Mass for ' num2str(time1) ' to ' num2str(time2) ' UTC'];
            savename=figname;
	
            idircont=1;
        
            for idir=1:length(direcDan)
                qdivn=sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[40:42]),3)./sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[43:45]),3);
                qdivn=sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[42]),3)./sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[43]),3);
                 xdat(idir).x=f*mean(qdivn,2);
                 ydat(idir).y=z/1000;
                 labs(idir).l=dire(idir).dir;
            end

%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              d=a-b-c;
%              dcont=acont-bcont-ccont;


            xlab=['Average Ice Particle Mass (kg)'];
            titlenam=figname;
            
	%         xdat(1).x=(a-acont)./acont*100;
	%         xdat(2).x=(b-bcont)./bcont*100;
	%         xdat(3).x=(c-ccont)./ccont*100;
	%         xdat(4).x=(-d+dcont)./dcont*100;
	

%             xdat(3).x=ccont;
%             xdat(4).x=dcont;
	

%             
%             ydat(3).y=z/1000;
%             labs(3).l='Fall Speed Flux';
%             
%             ydat(4).y=z/1000;
%             labs(4).l='Advection (-ve)';
            
            lor=1; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            %xlims=1;
            
            %xlimits=[-600 700];
            zmin=0.01;
            zmax=17.5;
            
            logflag=1;
            
            case 377
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            
            figname=['Average Surface Area for ' num2str(time1) ' to ' num2str(time2) ' UTC'];
            savename=figname;
	        
            rhoice=1000;
            
            for idir=1:length(direcDan)
                Vair=diff(GridDan(idir).Z(iz:iz2)).*GridDan(idir).RHON(iz+1:iz2)*diff(GridDan(1).Y1([1 end])); %total volume of air in each layer
                Vair=repmat(Vair,[1 length(dumprange(t1:t2))]);
                m=sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[42]),3).*Vair;
                surfmean=( (3*m/rhoice).^(2/3)*(4*pi)^(1/3) ) ./sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,[43]),3)./Vair;
                 xdat(idir).x=f*mean(surfmean,2);
                 ydat(idir).y=z/1000;
                 labs(idir).l=runName(idir).nam;
            end

%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              d=a-b-c;
%              dcont=acont-bcont-ccont;


            xlab=['Average Surface Area (m^2)'];
            titlenam=figname;
            
	%         xdat(1).x=(a-acont)./acont*100;
	%         xdat(2).x=(b-bcont)./bcont*100;
	%         xdat(3).x=(c-ccont)./ccont*100;
	%         xdat(4).x=(-d+dcont)./dcont*100;
	

%             xdat(3).x=ccont;
%             xdat(4).x=dcont;
	

%             
%             ydat(3).y=z/1000;
%             labs(3).l='Fall Speed Flux';
%             
%             ydat(4).y=z/1000;
%             labs(4).l='Advection (-ve)';
            
            lor=1; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            %xlims=1;
            
            %xlimits=[-600 700];
            zmin=0.01;
            zmax=17.5;
            
            logflag=1;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            case 38 %%%%%%%%%%%%%%%%%%%%           mean contents   %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=1;
            casestr=dire(idir).dir;
            
            ihm=[2 3 4];
            
            hmstr={'Liquid Water Content','Ice MR','Snow MR','Graupel MR','Ice NC','Snow NC','Graupel NC','Vapour MR','Updraught','Rain'};
            iq=[38 42 40 41 43 45 44 37 303 39]; %302  137 303=W>1_W
            
            
            aind=283; %285=ACu_A 280=ALu_A 283=W>1_A
            area=icediagsALL(idir).i(iz+1:iz2,t1:t2,aind)/npess2(idir);
            ilow=find(area<1e-7);
            area(area==0)=1e99;
            
          %  area=ones(size(area));
%            area(ilow)=1e99;
            
            iarea=0;
            if iarea==0
                area=1;
            end
            
            
            xlab=['Mean Mixing Ratio (ppmv)'];
            xlab=['Time Mean Mixing Ratio x Domain Size (g kg^{-1} km)'];
          %  xlab=['Time Mean Number Concentration x Domain Size (kg^{-1} km)'];
            xlab=['Time Mean Mixing Ratio x Domain Size (g kg^{-1})'];

            titlenam=['Time Average ' hmstr{ihm} ' for ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
            figname=[xlab titlenam];
            savename=figname;
	
            idircont=1;
            
            logflag=1;
            
            clear diff tot
        
            for idir=1:length(direcDan)
             xdat(idir).x=f*( mean( sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,iq(ihm)),3) ,2) ); %mean value of liq water
             
             xdat(idir).x=1000*( mean( sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,iq(ihm)),3)./area ,2) )  ./npess2(idir);  % * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; %mean value of liq water
             
         %    xdat(idir).x=( mean( sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,iq(ihm)),3)./area ,2) )  ./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000; %mean value of liq water
             
             %     pdat(1).p=1000*sum(icediagsALL(idir).i(izmin:izmax,dumprange,[38]),3)./npess2(idir) * length(GridDan(idir).Y1).*diff(GridDan(idir).Y1(1:2))/1000;

             rho=GridDan(idir).RHON(iz+1:iz2)*1e9; %convert to kg/km3 as xdat in g/kg km 
             dz=diff(GridDan(idir).Z(iz:iz2))/1000;
             
             tot(idir)=sum(xdat(idir).x.*dz.*rho);
                          
             ydat(idir).y=z/1000;
             labs(idir).l=runName(idir).nam;
            end

%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              d=a-b-c;
%              dcont=acont-bcont-ccont;


               
            
	%         xdat(1).x=(a-acont)./acont*100;
	%         xdat(2).x=(b-bcont)./bcont*100;
	%         xdat(3).x=(c-ccont)./ccont*100;
	%         xdat(4).x=(-d+dcont)./dcont*100;
	

%             xdat(3).x=ccont;
%             xdat(4).x=dcont;
	

%             
%             ydat(3).y=z/1000;
%             labs(3).l='Fall Speed Flux';
%             
%             ydat(4).y=z/1000;
%             labs(4).l='Advection (-ve)';
            
            lor=4; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            xlims=0;
            
            %xlimits=[-600 700];
        %    zmin=6;
        %    zmax=18;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            case 382 %%%%%%%%%%%%%%%%%%%%           mean contents from TwoD  %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=1;
            casestr=dire(idir).dir;
            
            ihm=[2 3 4]; %all ice MR
            %ihm=[5 6 7]; %all ice NC
            
            hmstr={'Liquid Water Content','Ice MR','Snow MR','Graupel MR','Ice NC','Snow NC','Graupel NC','Vapour MR','Updraught','Rain'};
            iq=[2 6 4 5 7 9 8 1]; 
            
            if min(ihm)<=4                    
                xlab=['Mean mixing ratio (ppmv)'];
                fact=f;
            else
          %  xlab=['Time Mean Mixing Ratio x Domain Size (g kg^{-1} km)'];
                xlab=['Mean number concentration (kg^{-1})'];
          %  xlab=['Time Mean Mixing Ratio x Domain Size (g kg^{-1})'];
              fact=1;
            end

%            titlenam=['Time Average ' hmstr{ihm} ' for ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
            titlenam=xlab;
            figname=[xlab titlenam];
            savename=figname;
	
            idircont=1;            
            logflag=0;
            
            clear diff tot
        
            for idir=1:length(direcDan)
             xdat(idir).x=fact*( mean( sum(TwoDDan(idir).Q(iz+1:iz2,:,iq(ihm)),3) ,2) ); %mean value
             
             rho=GridDan(idir).RHON(iz+1:iz2);   %*1e9; %convert to kg/km3 as xdat in g/kg km 
             dz=diff(GridDan(idir).Z(iz:iz2))/1000;
             
             tot(idir)=sum(xdat(idir).x.*dz.*rho);
                          
             ydat(idir).y=z/1000;
             labs(idir).l=runName(idir).nam;
            end

            
            lor=4; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            xlims=0;
            
            %xlimits=[-600 700];
        %    zmin=6;
        %    zmax=18;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            case 388 %%%%%%%%%%%%%%%%%%%%           max contents   %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=1;
            casestr=dire(idir).dir;
            
            ihm=[2];            
            hmstr={'vapour MR','tot MR'};            
            
            switch hmstr{ihm}
            case 'vapour MR'
                
                xlab=['Max mixing ratio (g kg^{-1})'];
                
                titlenam=['Max ' hmstr{ihm} ' for ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];                                
                
                for idir=1:length(direcDan)
                    xdat(idir).x=f* max( vap_prctiles(idir).t(iz+1:iz2,t1:t2,end),[],2 ); %mean value of liq water
                    ydat(idir).y=(GridDan(idir).Z(iz+1:iz2)+620)/1000;    
                    labs(idir).l=runName(idir).nam;
                end
                
            case 'tot MR'
                
                xlab=['Max mixing ratio (g kg^{-1})'];
                
                titlenam=['Max ' hmstr{ihm} ' for ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];

                for idir=1:length(direcDan)
                    xdat(idir).x=f* max( tot_prctiles(idir).t(iz+1:iz2,t1:t2,end),[],2 ); %mean value of liq water
                    ydat(idir).y=(GridDan(idir).Z(iz+1:iz2)+620)/1000;
                    labs(idir).l=runName(idir).nam;
                end
                
            end
            


%
            
            lor=1; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            xlims=0;
            
            figname=[xlab titlenam];
            savename=figname;
            
            %xlimits=[-600 700];
        %    zmin=6;
        %    zmax=18;
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            
            case 39
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=3;
            casestr=dire(idir).dir;
            
            titlenam=[casestr ' cf. Control for ' num2str(time1) ' to ' num2str(time2) ' UTC'];
            
           
            
            	
            idircont=1;
            
            iii=[43 45 44];
        
            for ii=1:3
                
             xdat(ii).x=( mean( sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,iii(ii)),3) ,2) - mean( sum(icediagsALL(idircont).i(iz+1:iz2,t1:t2,iii(ii)),3) ,2) ); %mean value of liq water
             ydat(ii).y=z/1000;
             
            end
            
            labs(1).l='ice';
            labs(2).l='snow';
            labs(3).l='graupel';
%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              d=a-b-c;
%              dcont=acont-bcont-ccont;


            xlab=['Difference of Time Mean Number Concs. from Control (kg^{-1})'];
            savename=[titlenam xlab];
            figname=savename;
            
	%         xdat(1).x=(a-acont)./acont*100;
	%         xdat(2).x=(b-bcont)./bcont*100;
	%         xdat(3).x=(c-ccont)./ccont*100;
	%         xdat(4).x=(-d+dcont)./dcont*100;
	

%             xdat(3).x=ccont;
%             xdat(4).x=dcont;
	

%             
%             ydat(3).y=z/1000;
%             labs(3).l='Fall Speed Flux';
%             
%             ydat(4).y=z/1000;
%             labs(4).l='Advection (-ve)';
            
            lor=0; %1=right, 2=left 0=auto
            
            xlimits=[-0.1 2]*1e7;
            xlims=0;
            
            %xlimits=[-600 700];
            zmin=0.01;
            zmax=10;
            
            logflag=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            case 40
            %casestr='500m res';
            %casestr='CCN=960 cm^{-3}';
            idir=2;
            casestr=dire(idir).dir;
            
            figname=['Diff from Control for ' num2str(time1) ' to ' num2str(time2) ' UTC'];
            savename=figname;
            	
            idircont=1;
            
            iii=[42 40 41];
        
            for ii=1:3
                
             xdat(ii).x=f*( mean( sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,iii(ii)),3) ,2) - mean( sum(icediagsALL(idircont).i(iz+1:iz2,t1:t2,iii(ii)),3) ,2) ); %mean value of liq water
             ydat(ii).y=z/1000;
             
            end
            
            labs(1).l='ice';
            labs(2).l='snow';
            labs(3).l='graupel';
%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,31:33),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[22 23 24]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              d=a-b-c;
%              dcont=acont-bcont-ccont;


            xlab=['Mixing Ratio (ppmv)'];
            titlenam=figname;
            
	%         xdat(1).x=(a-acont)./acont*100;
	%         xdat(2).x=(b-bcont)./bcont*100;
	%         xdat(3).x=(c-ccont)./ccont*100;
	%         xdat(4).x=(-d+dcont)./dcont*100;
	

%             xdat(3).x=ccont;
%             xdat(4).x=dcont;
	

%             
%             ydat(3).y=z/1000;
%             labs(3).l='Fall Speed Flux';
%             
%             ydat(4).y=z/1000;
%             labs(4).l='Advection (-ve)';
            
            lor=4; %1=right, 2=left 0=auto 4=bottom right
            
            xlimits=[-0.1 2]*1e7;
            xlims=0;
            
            %xlimits=[-600 700];
            zmin=0.01;
            
            logflag=0;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            case 41
        
        %casestr='500m res';
        %casestr='CCN=960 cm^{-3}';
        idir=2;
        casestr=dire(idir).dir;
            
        figname=[casestr ' Mean Ice cf. Control ' num2str(time1,'%2.2f') ' to ' num2str(time2,'%2.2f') ' UTC'];
        savename=['updown ' figname];
        
            idircont=1;
        
%              a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[85:87]),3) ); %difference in end and start
%              acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[85:87]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[85:87]),3) ); %difference in end and start vapour diags
% 
%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,76:78),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,76:78),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[67:69]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,67:69),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              dup=a-b-c;
%              dcontup=acont-bcont-ccont;
%              
%              
%              
%              
%              a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[130:132]),3) ); %difference in end and start
%              acont=f*( sum(icediagsALL(idircont).i(iz+1:iz2,t2,[130:132]),3) - sum(icediagsALL(idircont).i(iz+1:iz2,t1,[130:132]),3) ); %difference in end and start vapour diags
% 
%              b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,121:123),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%              bcont=f*TotMassBudgetProfRate2(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,121:123),3),GridDan(idircont).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%                           
%              c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[112:114]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%              ccont=sum(f*TotMassBudgetALL(GridDan(idircont),sum(icediagsALL(idircont).i(:,t1:t2,[112:114]),3),GridDan(idircont).t,iz,iz2),2); %fall speed flux source
%               
%              ddown=a-b-c;
%              dcontdown=acont-bcont-ccont;
             
             

             iup=[49:51]; %upwards flux of ice 
             iup=[67:69]; %downwards fall speed flux
             a=f*( mean(sum(icediagsALL(idir).i(iz+1:iz2,t1:t2,iup),3) ,2) );  %average amounts of ice in updraughts
             acont=f*mean(sum(icediagsALL(idircont).i(iz+1:iz2,t1:t2,iup),3) ,2);  %downdraughts
             
             xdat(1).x=a-acont;
             
             
             idown=[94:96];
             idown=[112:114];
             a=f*mean(sum(icediagsALL(idir).i(iz+1:iz2,t2,idown),3) ,2); 
             acont=f*mean(sum(icediagsALL(idircont).i(iz+1:iz2,t2,idown),3) ,2); 
             
             xdat(2).x=a-acont;
             
             
%             a=f*( sum(icediagsALL(idir).i(iz+1:iz2,t2,[40:42]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[40:42]),3) ); %difference in end and start vapour diags
%             b=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,31:33),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
%             c=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
%             xdat(3).x=a-b-c;
         
         
         


        xlab=['Difference of ' casestr ' Case from Control (ppmv)'];
        titlenam=figname;
        
%         xdat(1).x=dup-dcontup;
%         xdat(2).x=ddown-dcontdown;

        
        
        ydat(1).y=z/1000;
        labs(1).l='Updraught';
        
        ydat(2).y=z/1000;
        labs(2).l='Downdraught';
        
%         ydat(3).y=z/1000;
%         labs(3).l='Total';
        

        
        lor=1; %1=right, 2=left
        
        xlimits=[-0.5 500];
        xlims=0;
        
        xlimits=[-0.5 1];
        zmin=0.01;
        
        case 42
        titlenam=['Ice Number Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];    
            
         xdat(1).x=( sum(icediagsALL(idir).i(iz+1:iz2,t2,[43:45]),3) - sum(icediagsALL(idir).i(iz+1:iz2,t1,[43:45]),3) ); %difference in end and start vapour diags
         xdat(2).x=TotMassBudgetProfRate2(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,34:36),3),GridDan(idir).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
         xdat(3).x=sum(TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[25:27]),3),GridDan(idir).t,iz,iz2),2); %fall speed flux source
         
         %xdat(3).x=xdat(1).x(1:end)-xdat(2).x(1:end);
         
         xdat(4).x=xdat(1).x(1:end)-xdat(2).x(1:end)-xdat(3).x(1:end);
         
         %xdat(5).x=-sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,t1:t2,[4:6 13:15]),3),GridDan(idir).t,iz,iz2),2);
        
        xlab='Source of Ice Number (kg^{-1})';
        fignam=[titlenam '_' xlab];
        savename=fignam;
        
            ydat(1).y=z/1000;
            labs(1).l='Change';
            
            ydat(2).y=z/1000;
            labs(2).l='Microphysics';
            
            %ydat(3).y=z/1000;
            %labs(3).l='Advection + Fall Speed Flux';
            
            ydat(3).y=z/1000;
            labs(3).l='Fall Speed Flux';
            
            ydat(4).y=z/1000;
            labs(4).l='Advection';
            
            %ydat(5).y=z/1000;
            %labs(5).l='Flux Sources';
            
            xlims=1;
            xlimits=[-1 1]*1e9;
            
            zmin=0.01;
            zmax=22;
            
            lor=3;
            idirlabel=1;
            
        
        case 43 %plots the process rates for several cases on one graph
            iproc=30;
            
        figname=[dgs{iproc} ' Ice Number Sources for Times ' num2str(time1) ' to ' num2str(time2) ' UTC'];    
        
        icencsour=[59:62]; 
            icencsink=[52 53 63 64 45 46];
            
            snowncsour=[39 45 48];
            snowncsink=[58 50 56 43 40 47];
            
            gncsour=[56 40 41 42];
            gncsink=[57 49];
         
         files=[1:5];
            
         for idir=1:length(files)
            id=files(idir); 
            xdat(idir).x=TotMassBudgetProfRate2(GridDan(id),sum(icediag(id).i(:,t1:t2,iproc),3),GridDan(id).t,iz+1,iz2); %microphys - sums rate over the times given and multiplies by 300s
            ydat(idir).y=z/1000;
            labs(idir).l=direcDan(id).dir;    
         end     
         

        
        xlab='Source of Ice Number (kg^{-1})';
        titlenam=figname;
        

            %xlims=1;
            xlimits=[-2 1.2];
            
            zmin=0.1;
            zmax=22;
            
            lor=3; %1=right, 2=left
            
            
            
        end
        
        
        
  
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=GridDan(idir).PREFN(iz+1:iz2)/100;
    lab2='Pressure (hPa)';  
    dual=0;