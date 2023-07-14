function sVol=SampleVolumeCIP(tasVal,bin,sample_time)  
%function sVolGetSampleVolumeCIP(tasVal,bin,sample_time)  
%sVol is the sample volume per second of sample time in cm^3
%tasVal is the true airspeed (m/s)
%bin is the width of the particles in number of bins minus 1
%sample_time is the length of time that the sample was taken for in seconds
%
%This calculates the sample volume for the CIP instrument based on the statistical likelihood of
%having sampled a particle as a function of height. So the sample volume is scaled based on the size
%with larger particles being more likely to be rejected by the instrument (if touches sides) and thus 
%the sample volume is scaled down more than for small particles.


CIPRES = 25;    %CIPRES is the resolution in microns (25 for BAS)
max_dof = 100;   %max depth of field (mm) needs to be below the arm width (=10 cm)
arrayWidth = 64; %width of diode array (pixels)
Taur = 0.1;     %constant regarding electronics recovery time
CIPlambda = 0.0000658;	% wavelength of CIP laser (units? - think are centrimetres - see dof line)
    %this would make it 658 nm, which sounds right

radius=(bin+1) * CIPRES * 1e-4 * 0.5; 	%radius is in cm                 
f=100*min(0.5./(1-exp(-1000*2*10*radius./(tasVal*Taur))),1);
z=30.46-(0.628*f)+(0.003246*f.*f); %z is a correction factor for instrument dead time - must be dimensionless as
    %has f and f^2 terms.
dof = z*10.*radius.*radius/CIPlambda; %depth of field is converted to mm by the factor of 10 here
    %can assume that z has no units (think is just a scaling factor) then this means that CIPlambda is in centimetres
    %to allow dof to be in mm (as says in the CIP manaul in Y:/dataAnalysisUsersGuideChapterII.doc)
    %Nmetres = Nm = 1e-2*Ncm = 1e-3*Nmm  so dof = 10 * 1e2 Nm *1e2 Nm / CIPlambda = Nmm = Nm*1e3
    %this leads to CIPlambda = 1e2 * Nm = Ncm, so CIPlambda is in centimetres
dof = min(dof, max_dof);
eaw=CIPRES*1.0e-3*((arrayWidth-2)-bin); %effective sample width in mm (or rather the free sample width 
            %left available to the side of the particle image - the less this is the less likely it is that 
            %the particle would have been sampled and so it gets scaled down (and thus concentrations increased)
sVol=eaw.*dof.*tasVal*sample_time;  %cm^3 = mm * mm * m/s * s
	