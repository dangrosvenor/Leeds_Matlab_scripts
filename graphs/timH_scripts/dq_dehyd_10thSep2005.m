fact=1e6*28.97/18;
%dumprange=[1:55];

	logflag=0;
    
 	tit(1).tit='Loss Due to Total Water Areas Below 3.67 ppmv (ppbv)';
%	tit(1).tit='Loss Due to Total Water Areas Below 5 ppmv (ppbv)';
    tit(1).tit='Sum of Total Water Below 5 ppmv (ppmv)';
    tit(1).tit='Sum of Water Vapour Points Below 5 ppmv * Grid Length (ppmv km)';

    nplots2d=1;
    
    clines=0; %makes black contour lines appear
    clab=0;
    manclab=0;
    
    minZ=15.3e3;
    minZ=13e3;
    maxZ=19e3;
%    maxZ=22e3;
    ncont=20;
    
    dlogflag=0;
    dlogmin=100;
    
    imaxovr=0;
    iminovr=0;
    
    mincovOvr = dlog(1.000000,dlogmin);
maxcovOvr = dlog(380.000000,dlogmin);

    mincovOvr = dlog(0.000000,dlogmin);
    maxcovOvr = dlog(500.000000,dlogmin);

    
%    mincovOvr=log10(1);
    %mincovOvr=log10(0.03); %mass
    
 %   mincovOvr = 0.000000;
%maxcovOvr = 60.000000;
    
    z=GridDan(idir).Z;
    time=GridDan(idir).t+3;
    
    savename=['dq_dehyd_1-' num2str(dumprange(end))];
    
    %three=repmat(3.67,size(dq(idir).d(izmin:izmax,dumprange,1));
    
%     pdat(i).p=1000*(dq(idir).d(izmin:izmax,dumprange,1)...
%         - nn(idir).n(izmin:izmax,dumprange)/502 .* (3.67 - f*sum(icediag4(idir).i(izmin:izmax,dumprange,[35:36]),3) ) );

    % pdat(i).p=1000*( dq_tot(idir).d(izmin:izmax,dumprange,1) ) ;	%reduction of layer tot water if points < 3.67 were mixed
    
    dy=(GridDan(idir).Y1(2)-GridDan(idir).Y1(1))/1000;
	pdat(i).p=length(GridDan(idir).Y1)*( dq_tot(idir).d(izmin:izmax,dumprange,2) ) *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more
	pdat(i).p=length(GridDan(idir).Y1)*( dq_vaps(idir).d(izmin:izmax,dumprange,2) ) *dy; %multiply by dy so is in ppmv*km since otherwise high res will mean there are more

    %*** WHY multiplying by length(GridDan(idir).Y1???                                 %points counted
     %*** - is because was divided by this in Allimp...m      
    
   % pdat(i).p=1000*(dq_tot(idir).d(izmin:izmax,dumprange,1)...
   %      - nn(idir).n(izmin:izmax,dumprange)/502 .* (3.67 - 3.8 ) );

%    initialq=repmat(f*sum(icediag4(idir).i(izmin:izmax,dumprange(1),[35:36]),3),[1 length(dumprange)]);
%    
%    pdat(i).p=1000*(dq_tot(idir).d(izmin:izmax,dumprange,1)...
%          - nn(idir).n(izmin:izmax,dumprange)/502 .* (3.67 - initialq ) );
     
     %pdat(i).p=3.67 - 502./nn(idir).n(izmin:izmax,dumprange) .* dq_tot(idir).d(izmin:izmax,dumprange,1);
     %pdat(i).p(isnan(pdat(i).p))=0;
    
%     dz=repmat(z(2:end)-z(1:end-1),[1 length(dumprange)]); 
%     rho=repmat(GridDan(idir).RHON(izmin:izmax),[1 length(dumprange)]);
% 	pdat(i).p=( dq_tot(idir).d(izmin:izmax,dumprange,1) .*rho .*dz(izmin:izmax,:) ) ;

    
    %icontovr=1;
    %conts=[1 2 4 8 