function [lwc]=adLWC(pdat,potemp,zdat,qdat)
% calculates the adiabatic LWC profile (LWC in g/kg)


        tdat=potemp./(1e5./pdat).^0.286; %temperature in K        
		
        f=1e6*28.97/18;
		qsat=satvappress(tdat,'goff','liq',pdat,1)/f;
		
        [CAPE,CIN,HLCL,TLCL,PLCL]=calc_cape(pdat,tdat,qdat,qsat,zdat);                      
        ipos=findheight(zdat,HLCL);  
        
        th_start=potemp(ipos);
        
		[tad,th_grid,p_grid]=moist_adiabat2(th_start,pdat(ipos),pdat(end)); %temperaure, potemp and pressure during moist saturated rise
        zdat2=interp1(pdat,zdat,p_grid); %find corresponding altitudes

            qsat2=satvappress(tad,'goff','liq',p_grid,1)/f; %the saturation mixing ratio during the adiabatic ascent - i.e. the vapour value
                                                            %of the parcel as are assuming it's at saturation all the way up
            
            lwc=1000*(qsat(ipos)-qsat2); %the adiabatic LWC is the initial vapour value of the parcel at LCL minus  
                                                  %value maintained during ascent (=saturation mixing ratio)
                                                  
            lwc = interp1(p_grid,lwc,pdat); %interpolate back onto the original pdat grid