%cumulative top-down time height plot

clear dz

logflag=0;
    fact=1e6*28.97/18;
    fact=1;
    
    z=GridDan(idir).Z(2:end); %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
    %len=length(icediag4(idir).i(:,1,35))-1;
    %len=findheight(z,18.25e3)+1;
    %len=length(z)-1;
    
	ncont=24;
    
    iminovr=1;
%    mincovOvr=(1e-7);
    mincovOvr=(-0.1e-3)*1000;

    
    imaxovr=1;
%    maxcovOvr=(0.05);
    maxcovOvr=(0.1E-3)*1000;

        
	sig=2;

    
    minZ=10e3;
    maxZ=z(end);
    
    iylim=1;
    iylims=[14.7 30];
    iylims=[15.1965   22.6869];
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    %i2d=2; %tells it to label x axis in km   
   
  if (~exist('phase')); phase=1; end
      
    if (phase==2) %if are in part of plot program that wants the data in pdat
        izmax=izmax;
        len=izmax+1; %set izmin so that is one below max
        hstr=num2str( round2(GridDan(idir).Z(len-1)/1000,1) );
        
        tit(1).tit=['Cumulative Net Gain Down in Vapour from Domain Top (g/m^{2})'];
        %tit(1).tit=['Cumulative Net Gain in Total Water Down from Domain Top (g/m^{2})'];
        
        
        %tit(1).tit=['Averaged Net Gain Down from ' hstr 'km (g/m^{2})'];
        %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	    %tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
	    figlab=['Cumulative net gain down from ' hstr 'km'];
    
        
        dz(2:len,1:size(icediag4(idir).i,2))=repmat(GridDan(idir).Z(2:end)-GridDan(idir).Z(1:end-1),[1 size(icediag4(idir).i,2)]); 
        rho=repmat(GridDan(idir).RHO(:),[1 size(icediag4(idir).i,2)]);
        
        start=repmat(sum(icediag4(idir).i(:,dumprange(1),[35:36]),3),[1 size(icediag4(idir).i,2)]);
        
        if  izmin==1; izmin=2; end
        
        for icu=izmin+1:len
            %me=mean(icediag3(1).i(len-i+1:len,end,11));
            %xdat(1).x(len-i+1) = sum(rho(len-i+1:len).*dz(len-i+1:len).*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ); %difference in end and start vapour diags                                 
            %airmass=sum( rho(len-i+1:len) .*dz(len-i+1:len) );
            airmass=1;
            
            
            
            %converted into mixing ratio
                %pdat(i).p(len-icu+1,1:length(dumprange)) = sum( rho(len-icu+1:len,dumprange) .*dz(len-icu+1:len,dumprange) .* ( sum(icediag4(idir).i(len-icu+1:len,dumprange,[35:36]),3) - start(len-icu+1:len,dumprange) )  ,1) ./ airmass;

            pdat(i).p(len-icu+1,1:length(dumprange)) = 1000 * sum( rho(len-icu+izmin:len-1,dumprange) .*dz(len-icu+izmin+1:len,dumprange) .* ( sum(icediag4(idir).i(len-icu+izmin+1:len,dumprange,[35 36]),3) - start(len-icu+izmin+1:len,dumprange) )  ,1) ./ airmass;
            %pdat(i).p(len-icu+1,1:length(dumprange)) = 1000 * sum( rho(len-icu+izmin:len-1,dumprange) .*dz(len-icu+izmin+1:len,dumprange) .* ( sum(icediag4(idir).i(len-icu+izmin+1:len,dumprange,[35:42]),3) - start(len-icu+izmin+1:len,dumprange) )  ,1) ./ airmass;

            %xdat(2).x(len-i+1) = f* sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediag4(idir).i(len-i+1:len,end,[35:36]),3) - sum(icediag4(idir).i(len-i+1:len,3,[35:36]),3) )  ) ./ airmass;
	
            
            %mass gain
            %xdat(1).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .* ( sum(icediag3(1).i(len-i+1:len,end,[2 5 8 11]),3) - sum(icediag3(1).i(len-i+1:len,1,[2 5 8 11]),3) )  ) ;
            %xdat(2).x(len-i+1) =  sum( rho(len-i+1:len) .*dz(len-i+1:len) .*( icediag3(1).i(len-i+1:len,end,11) - icediag3(1).i(len-i+1:len,1,11) ) ) ; %difference in end and start vapour diags    
        end
    
end
    
  