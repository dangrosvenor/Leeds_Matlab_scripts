if izmin==1|izmin==2; izmin=3; end

fallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[22 23 24]),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 
fallnc=1/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[25:27]),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in ice number due to fall speed flux 

clear fall_from_mean
if exist('tot_fallflux');
    for ifall=1:length(tot_fallflux)
		fall_from_mean(ifall).i=f/300*TotMassBudgetALL(GridDan(idir),tot_fallflux(ifall).i(:,:),GridDan(idir).t,izmin-1,izmax); %rate of increase in tot water due to fall speed flux 
	end
end

%fallrate=fall_from_mean;
% 
%     fluxrate= -f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[1:6]),3),GridDan(idir).t,izmin-1,izmax); %rate of tot increase due to flux of total water (wq)
    
 %   pdat(1).p=f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[22 23 24]),3),GridDan(idir).t,izmin-1,izmax)...
  %      - f*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,[1:6]),3),GridDan(idir).t,izmin-1,izmax);
    
    init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[37:42]),3),[1 length(dumprange)]);
    change=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[37:42]),3) - init )/npes; %difference in end and start total water diags
    
    init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[40:42]),3),[1 length(dumprange)]);
    changeice=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40:42]),3) - init )/npes;
    
    init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[37]),3),[1 length(dumprange)]);
    changevap=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[37]),3) - init )/npes;
    
    clear diff;
    dz=diff(GridDan(idir).Z)/1000;
    dzz=repmat(dz,[1 length(dumprange)]);    
    rho=repmat(GridDan(idir).RHON(2:end),[1 length(dumprange)]);
    init=repmat(sum(icediagsALL(idir).i(2:end,1,[37:42]),3),[1 length(dumprange)]);
    change2=f*( sum(icediagsALL(idir).i(2:end,dumprange,[37:42]),3) - init )/npes; %difference in end and start total water diags
    
    change2=f*( sum(icediagsALL(idir).i(2:end,dumprange,[22:24]),3))/npes; 
    topdown=1000*flipud(cumsum(flipud(change2.*rho.*dzz),1));
    %topdown=topdown(izmin:izmax,:);
    
    %microrate=sum(icediagsALL(idir).i(izmin:izmax,dumprange,28:33),3);
    %microchange=cumsum(microrate,2)*300; %sum rate cumulativley and multiply by dt for change due to microphysics
    
   %pdat(1).p=change-cumsum(fallrate,2); %cumulative sum of fallrate to give the contribution to change of fall flux at each time
  % pdat(1).p=cumsum(fallrate,2); %cumulative sum of fallrate to give the contribution to change of fall flux at each time
   %pdat(1).p=change;
   
   ad=change-cumsum(fallrate,2)*300; %change due to advection sum(dy)=sum(dy/dt *dt) = sum(dy/dt) * dt
      %cumulative total water from ad          %note sum(dy/dt*dt) NE sum(dy/dt)*sum(dt)
   
   fluxrate2=(ad(:,2:end)-ad(:,1:end-1))/300;
   fluxrate2(:,dumprange(2):dumprange(end))=fluxrate2(:,1:end);
   fluxrate2(:,1)=0;
%    
%    changerate=(change(:,2:end)-change(:,1:end-1))/300;
%    changerate(:,2:dumprange(end))=changerate(:,1:end);
%    changerate(:,1)=0;
   
   microicerate=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,31:33),3)/npes; %ice mixing ratio source rate
   microicerate=-f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,28),3)/npes; %ice mixing ratio (of all ice) source rate
   
   init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[42]),3),[1 length(dumprange)]);
   HMicechange=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[42]),3) - init )/npes;
   
   init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[40]),3),[1 length(dumprange)]);
   HMsnowchange=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[40]),3) - init )/npes;
   
   init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[41]),3),[1 length(dumprange)]);
   HMgrachange=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[41]),3) - init )/npes;
   
	HMicefallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,24),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 
	HMsnowfallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,22),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 
	HMgrafallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,23),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 

    HMicefallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,24),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 
	HMsnowfallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,22),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 
	HMgrafallrate=f/300*TotMassBudgetALL(GridDan(idir),sum(icediagsALL(idir).i(:,dumprange,23),3),GridDan(idir).t,izmin-1,izmax)/npes; %rate of increase in tot water due to fall speed flux 
    
    
    
   HMicerate=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,33),3)/npes; %ice mixing ratio source rate
   HMsnowrate=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,31),3)/npes; %ice mixing ratio source rate
   HMgrarate=f*sum(icediagsALL(idir).i(izmin:izmax,dumprange,32),3)/npes; %ice mixing ratio source rate
   
   micronc=sum(icediagsALL(idir).i(izmin:izmax,dumprange,34:36),3)/npes; %ice number source rate
   
   initvap=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[37]),3),[1 length(dumprange)]);
   changevap=f*( sum(icediagsALL(idir).i(izmin:izmax,dumprange,[37]),3) - initvap )/npes; %difference in end and start vapour diags
   
   vapadcum=-( changevap + cumsum(microicerate,2)*300 ); %is plus here since the microphysical rate used is that for ice
      %made negative so is the advective loss of vapour
   vapad=(vapadcum(:,2:end)-vapadcum(:,1:end-1))/300; %calculate rate of change of advective source
   vapad(:,dumprange(2):dumprange(end))=vapad(:,1:end);
   vapad(:,1)=0;
   
   iceadcum=changeice - cumsum(microicerate,2)*300 - cumsum(fallrate,2)*300;
   icead=adratef(iceadcum,dumprange); %this function does the same as the 3 vapad= lines above
   
   HMiceadcum=HMicechange - cumsum(HMicerate,2)*300 - cumsum(HMicefallrate,2)*300;
   HMsnowadcum=HMsnowchange - cumsum(HMsnowrate,2)*300 - cumsum(HMsnowfallrate,2)*300;
   HMgraadcum=HMgrachange - cumsum(HMgrarate,2)*300 - cumsum(HMgrafallrate,2)*300;
   
   
   
   init=repmat(sum(icediagsALL(idir).i(izmin:izmax,1,[43:45]),3),[1 length(dumprange)]);
   changenc=(sum(icediagsALL(idir).i(izmin:izmax,dumprange,[43:45]),3) - init )/npes;
%    changenc=(changenc(:,2:end)-changenc(:,1:end-1))/300;
%    changenc(:,2:dumprange(end))=changenc(:,1:end);
%    changenc(:,1)=0;
   
   ncadcum=changenc - cumsum(micronc,2)*300 - cumsum(fallnc,2)*300; %is plus here since the microphysical rate used is that for ice
   
   adnc=(ncadcum(:,2:end)-ncadcum(:,1:end-1))/300; %calculate rate of change of advective source
   adnc(:,dumprange(2):dumprange(end))=adnc(:,1:end);
   adnc(:,1)=0;
 
  init=repmat(sum(icediagsALL(idir).i(:,1,[151]),3),[1 length(dumprange)]);
  changelow=(sum(icediagsALL(idir).i(:,dumprange,[151]),3) - init )/npes;
  adlow=adratef(changelow,dumprange); 
   
   liqsink=[13 3 18 5 32 29 33 34]; %dq02 = 29 for ALL