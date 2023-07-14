function [LWC,DENS,VISC,VSCW,CND,CNDW,RE,PRF,PRW,TB,P,DRYP,FACT]=lwc_calc_from_DMT_PACS_comparison(T,PMB,TAS,V,TB)
%function LWC=lwc_calc_from_DMT(T,PMB,TAS,V)
%DMT liquid water content code from the document 
%LWC_Manual_DOC-0067_Rev-A.pdf in folder
%\My Documents\logbook\Antarctica\Flights and instruments_Feb2010\
%T = air temp in celsius
%PMB = pressure in mb
%TAS = air speed in m/s
%V = voltage output from LWC probe (0-10V)


%TW = wire temperature in sensor (Celsius) (125 degC)
%L = length of sensor (2 cm)
%D = diameter of sensor (0.18 cm)


%TW=125;
%TW=157;
%TW=115.2;
%TW=110.5;   %the approx value when using CIP airspeed
%TW=105.5;
%TW=100.5; %when are using the aircraft airspeed (TAS_aircraft in read_CAS_data_textscan2)
%TW=150;
TW=125;

L=2;
D=0.18;

%TW=100; %Tom LC's script gives TW and L values of 100 degC and 1.5 cm
%L=1.5;


TK=T+273.16;
TWK=TW+273.16;
TFLM=(TWK+TK)/2;

%convert volts to watts
P=10*V; %must be I=10A? No, since P not in Watts  - see below FACT calculation

%calculate the thermal conductivity
CND = 5.8e-5*(398./(125+TFLM)).*(TFLM/273).^1.5;  %is 125 here the wire temperature??
CNDW = 5.8e-5*(398./(125+TWK)).*(TWK/273).^1.5;

%calculate the viscosity
VISC = 1.718E-4*(393./(120+TFLM)).*(TFLM/273).^1.5;
VSCW = 1.718E-4*(393./(120+TWK)).*(TWK/273).^1.5;

%calculate the density
DENS = PMB./(2870.5*TFLM);
FCT = 3.14159*L*CND.*(TWK-TK);

%The Reynolds number
RE = 100*DENS.*TAS*D./VISC;

%The Prandtl numbers
PRF = 0.24*VISC./CND;
PRW = 0.24*VSCW./CNDW;

%TB=373.16; %boiling point of water - in DMT spreadsheet this varies (presumably according to pressure?) 
%TB=372.734741; 
%TB=364.48996;

%Calculate the dry air loss
DRYP = 0.26*RE.^0.6.*PRF.^0.37.*(PRF./PRW).^0.25.*FCT/0.239;
FACT = 1.238E6*0.239./(L*D*TAS*100.*(597.3+TB-TK)); %373.16 = boiling point of water
%1e6 is here to convert from cm^3 to m^3 for LDV. 100 converts V from m/s to cm/s
%the 597.3 doesn't quite look right for ratio of Lv/SHC, which it works out as
%P is probably not in Watts as Lv and SHC not in SI

%FACT=FACT/1.137; %factor that seems to be indicated by spreasheet differences
%but the formulae here are the same as those from Darrel Baumgardner's script

%the FACT below gets a similar answer - chice of Lv and SHC water prob the difference
%FACT = 1./(1e-4*L*D.*TAS*1e-3.*(2260e3 + 4260*(373.16-TK)));

%Calculate the LWC
LWC = (P-DRYP).*FACT;

%N.B. - DMT use a non-constant boiling point of water - varies with pressure I guess?
iprint=0;
if iprint==1
    fprintf(1,'\n%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n'...
        ,LWC,0,0,0,T,PMB,TAS,DENS,VISC,VSCW,CND,CNDW,RE,PRF,PRW,TB,P,DRYP,0,0,FACT);    
end
%   Columns from the DMT PADS hotwire LWC output (Excel spreadsheet - 02Hotwire_LWC20090218105850_matlab_comparison.csv)
%    in folder C:\Documents and Settings\dan\My Documents\logbook\Antarctica\Flights and instruments_Feb2010
% LWC hotwire(V)	
% LWC Slave(V)	
% LWC raw (g/m^3)           - Think this is FACT*P
% DAT Calculated (g/m^3)	- Dry Air Term - think this is FACT*DRYP
% DAT Observed (g/m^3)	- seems to be fairly constant all through the flight
%                       - actually changes in steps
% LWC-DAT Calc (g/m^3)	-
% LWC-DAT Observed(g/m^3)	
% Ambient Temp(C)	
% Pressure(mB)	
% Airspeed(m/s)	
% Air Density(kg/m^3)	
% Viscosity Dry(g/sec-cm)	
% Viscosity Wet(g/sec-cm)	
% Thermal Cond Dry(cal/sec-cm-K)	
% Thermal Cond Wet(cal/sec-cm-K)	
% Reynold's Number	
% Pradtl Numb Dry	
% Pradtl Numb Wet	
% Boiling Point Water (K)	
% P Total (W)	
% P Dry Calculated (W)

