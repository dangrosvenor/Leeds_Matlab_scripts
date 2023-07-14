function [CAPE,CIN,HLCL,TLCL]=CALC_CAPE_v1(P,T,R,ALT)
% This routine calculates the CAPE from a sounding
%P in Pa, T in K, R=mixing ratio (kg/kg)

N=max(size(P));

TLP=zeros(N,N);
TLVP=zeros(N,N);
TVPDIF=zeros(N,N);
CAPEP=zeros(N);
NAP=zeros(N);
PAP=zeros(N);


CPD=1005.7;	% specific heat dry air
CPV=1870.0;	% ?
CL=4190.0;	% specific heat water at 0deg
CPVMCL=2320.0;	
RV=461.5; % gas constant water vapor = R*/M_H2O
RD=287.04;	% gas constant dry air = R*/M_dry
EPS=RD/RV;
ALV0=2.501E6;	% "Verdampfungswaerme"
T0=273.15;

% water vapour pressure ev and sat vap press
TC=T-T0; % celsius
EV=R.*P./(EPS+R); % vapour press

%ES=GGEW(T); %sat vapour press (Pa)

RS=esLES(P,T);
ES=P.*RS/(RS+EPS);



%ES=RSAT.*P./(EPS+RSAT); % vapour press

% begin outer loop
for I=1:N./3
    %RS=EPS.*ES(I)./(P(I)-ES(I));
    ALV=ALV0-CPVMCL.*TC(I);
    EM=(EV(I) > 1E-6).*EV(I) +(EV(I) <1E-6).*1E-6;
    SP=CPD.*real(log(T(I)))-RD.*real(log(P(I)-EV(I))) + ...
        ALV.*R(I)./T(I)-R(I).*RV.*real(log(EM./ES(I)));
    
    % Find lifted condensation pressure PLCL 
    RH=R(I)./RS(I);     % relative humidity
    RH=(RH < 1.0).*RH+(RH > 1.0);
    CHI=T(I)./(1669.0-122.0.*RH-T(I));
    PLCL=(RH > 0.0).*P(I).*RH.^CHI+(RH < 0.);
    
    %   ***  Begin updraft loop   ***

    SUM=0.0;
    RG0=R(I);
    TG0=T(I);
    for J=I:N 

        %***  Calculate estimates of the rates of change of the entropies  ***
        %***           with temperature at constant pressure               ***
  
        %RS=EPS.*ES(J)./(P(J)-ES(J)); % saturation mixing ratio
        ALV=ALV0-CPVMCL.*TC(J);
        SLP=(CPD+RS(J).*CL+ALV.*ALV.*RS(J)./(RV.*T(J)*T(J)))./T(J);
   
        %***  Calculate lifted parcel temperature below its LCL   ***

        if P(J) >= PLCL 
           TLP(I,J)=T(I).*(P(J)./P(I)).^(RD./CPD);
           TLVP(I,J)=TLP(I,J).*(1.+R(I)./EPS)./(1.+R(I));
           TVPDIF(I,J)=TLVP(I,J)-T(J).*(1.+R(J)./EPS)./(1.+R(J));
        else

       %***  Iteratively calculate lifted parcel temperature and mixing   ***
       %***    ratios for pseudo-adiabatic ascent     ***

           TG=T(J);
           RG=RS(J);
           for K=1:7
	       CPW=(J > 0).*( SUM+CL.*0.5.*(RG0+RG).*(real(log(TG))-real(log(TG0))));
               EM=RG.*P(J)./(EPS+RG);
               ALV=ALV0-CPVMCL.*(TG-273.15);
               SPG=CPD.*real(log(TG))-RD.*real(log(P(J)-EM))+CPW+ALV.*RG./TG;
               TG=TG+(SP-SPG)./SLP;  
	       %ENEW=GGEW(TG);
               %RG=EPS.*ENEW./(P(J)-ENEW);
               RG=esLES(P(J),TG);
           end % K
           TLP(I,J)=TG;
           TLVP(I,J)=TG.*(1.+RG./EPS)./(1.+RG);
           TVPDIF(I,J)=TLVP(I,J)-T(J).*(1.+R(J)./EPS)./(1.+R(J));
           RG0=RG;
           TG0=TG;
           SUM=CPW;
       end
   end % J
   
   
   %***  Find positive and negative areas  PA and NA and
   %     CAPE (=PA-NA) from pseudo-adiabatic ascent ***



   %***  Find lifted condensation level and maximum level   ***
   %***               of positive buoyancy                  ***

   ICB=N;	% index of lifted cond level
   INBP=0;	% index of maximum level of positive buoyancy
   for J=(N):-1:I
        if P(J) < PLCL 
            ICB=min([ICB,J]);
        end
        if TVPDIF(I,J) > 0.0
           INBP=(J > INBP).*J+(J < INBP).*INBP;
	       break;
       end
   end % J
    
   IMAX=(INBP > I).*INBP+(INBP <= I).*I;
   TVPDIF(I,IMAX:N)=0.;	% set to zero for levels above IMAX

   %***  Do updraft loops        ***

   if INBP > I
       for J=(I+1):INBP
           TVM=0.5.*(TVPDIF(I,J)+TVPDIF(I,J-1));
           PM=0.5.*(P(J)+P(J-1));
           if TVM < 0.0 
              NAP(I)=NAP(I)-RD.*TVM.*(P(J-1)-P(J))./PM;
           else
              PAP(I)=PAP(I)+RD.*TVM.*(P(J-1)-P(J))./PM;
           end 
       end %J
       CAPEP(I)=PAP(I)-NAP(I);
   end % else cape=0 if no free convection is possible
 
end %I	; loop over air parcel origins

CAPE=max(max(CAPEP));
CIN=max(max(NAP));
max(max(PAP));





    RS=EPS.*ES(1)./(P(1)-ES(1));
    ALV=ALV0-CPVMCL.*TC(1);
    EM=(EV(1) > 1E-6).*EV(1) +(EV(1) <1E-6).*1E-6;
    SP=CPD.*real(log(T(1)))-RD.*real(log(P(1)-EV(1))) + ...
        ALV.*R(1)./T(1)-R(1).*RV.*real(log(EM./ES(1)));
    
    % Find lifted condensation pressure PLCL 
    RH=R(1)./RS;     % relative humidity
    RH=(RH < 1.0).*RH+(RH > 1.0);
    CHI=T(1)./(1669.0-122.0.*RH-T(1));
    PLCL=(RH > 0.0).*P(1).*RH.^CHI+(RH < 0.);
    
    % find the height
    %HLCL=interp1(P,ALT,PLCL); %taken out because P data has some points the same as previous ones
    %TLCL=interp1(P,T-273,PLCL);
    
