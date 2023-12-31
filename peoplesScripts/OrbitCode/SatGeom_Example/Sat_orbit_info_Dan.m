%Demonstrate Calculation of Beacon Satellite Propagation Geometry
%
%    Calculates the propagation geometry and magnetic field
%    parameters from spg10sys orbit and IRGF11 B-field predictions.  For 
%    satellites other than object 32765, tsince parameters and sampling
%    intervals should be adjusted
%
%    The outputs are stored in a .mat file for subsequent display compuation
%
%Required Libraries: GPS_CoordinateXforms, SPG4, IGRF, utilities 
%

%Written by
%Charles L. Rino
%Rino Consulting
%crino@mindspring.com
%Version: September 28, 2010

%modified by Dan G, 18th April, 2012

use_own_files=1; %flag to say that want to bypass the demo files

dtr=pi/180;
min_per_day=60*24;

%***************Get TLE Orbital Elements**********************
%path2data='..\SGP4\DemoData';

path2data='C:\Users\Dan\Documents\MATLAB\work\peoplesScripts\OrbitCode\SGP4\DemoData\';
TLEfiles=dir([path2data,'\*.mat']);

%set filename for satellite TLE here
if use_own_files==1 %overwrite the pathname for the satellite TLE
    %Two-line Element file for TERRA (although TLE files may change with time...)
    TLEfiles(2).name = 'TERRA-TLE.mat';
end

if isempty(TLEfiles)
    error('No TLE file')
else
    %DemoData currently has only a single set of elements
    %ADD library and selection options here for more general use
    load(fullfile(path2data,TLEfiles(2).name));
    longstr1=line1;
    longstr2=line2;
    fprintf('USING 2-Line Elements: \n')
    fprintf('%s \n',longstr1)
    fprintf('%s \n',longstr2)
%    load(fullfile(path2data,TLEfiles(1).name));  %
end

%****************Get Station Location**************************
stationfiles=dir([path2data,'\*.mat']);

if isempty(stationfiles) & use_own_files==0
   error('No station files file')
else
    % *** set station/location lat, lon and altitude here ***
    if use_own_files==1
%        lct.station.rx_latitude=79.5848;
%        lct.station.rx_longitude=106.2550;
        lct.station.rx_latitude=50.5571;
        lct.station.rx_longitude=-37.1856;
        lct.station.rx_altitude=0;
        lct.station.rx_name=[num2str(lct.station.rx_latitude) '-lat,'  num2str(lct.station.rx_longitude) '-lon,' num2str(lct.station.rx_altitude) '-alt'];
    else
        load([path2data,'\',stationfiles(1).name]);
    end

   origin_llh=[lct.station.rx_latitude*dtr;...
               lct.station.rx_longitude*dtr;...
               lct.station.rx_altitude];
   rx_name=lct.station.rx_name;
end

%********************Initialize spg4****************************
%set in sgp4init
%global tumin mu radiusearthkm xke j2 j3 j4   
%global opsmode
satrec = twoline2rvMOD(longstr1,longstr2);
rx_name=sscanf(rx_name,'%s');
fprintf('\n')
fprintf('Satellite ID %5i \n',satrec.satnum)
fprintf('Station %s: Lon=%6.4f Lat %6.4f Alt=%6.2f m \n',...
            rx_name,origin_llh(2)/dtr,origin_llh(1)/dtr,origin_llh(3))

%********************Run SPG4 for Multiple Orbits****************            
%NOTE:  dt, npts, and tsince offset tailored to object 32765 parameters 
%specify the time after the epoch (satrec.epochdays) that we want orbit
%data for
dt=10/60;  %10 sec (dt is therfore in minutes)
nhrs = 24;
npts=ceil(nhrs/(dt/60)); %no. of steps of dt that we want orbit info for
tsince_offset=0; %tsince_offset is the start time for the orbit info we want (in minutes)
%Expand dt, npts, and tsince_offset for broader selection capabilities

tsince=tsince_offset+[0:npts-1]*dt; %minutes
if (satrec.epochyr < 57)
    Eyear= satrec.epochyr + 2000;
else
    Eyear= satrec.epochyr + 1900;
end
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,0);
UTsec=Ehr*3600+Emin*60+Esec;
fprintf('EPOCH: YEAR MO  DAY UTSEC \n')
fprintf('      %5i %2i %4i %5.2f \n\n',Eyear,Emon,Eday,UTsec);



xsat_ecf=zeros(3,npts);
vsat_ecf=zeros(3,npts);
for n=1:npts
   [satrec, xsat_ecf(:,n), vsat_ecf(:,n)]=spg4_ecf(satrec,tsince(n));
end
%Scale state vectors to mks units
xsat_ecf=xsat_ecf*1000;  %m
vsat_ecf=vsat_ecf*1000;  %mps

%ecf is Earth Centred Earth Fixed
%The Earth-centered earth-fixed (ECEF or ECF) or conventional terrestrial 
%coordinate system rotates with the Earth and has its origin at the centre of the Earth. 
%The X axis passes through the equator at the prime meridian.
%The Z axis passes through the north pole but it does not exactly coincide 
%with the instantaneous Earth rotational axis.
%The Y axis can be determined by the right-hand rule to be passing through
%the equator at 90� longitude.


%llh means  [LAT,LON,HEIGHT]. ecf2llhT produces a [3xtime] array of this for
%the satellite
sat_llh=ecf2llhT(xsat_ecf);            %ECF to geodetic (llh) 
% -------------------------------------------------------------------------
% *** if just want [LAT,LON,HEIGHT] of the satellite then can stop here ***
% *** LLH is independent of the ground location (origin) ***
% -------------------------------------------------------------------------
%function llh=ecf2llhT(ecf)
%
% Description:
%   This function returns the [lat lon hgt] LLH coordinates
%   of the given earth center fixed (ECF) coordinates.
%
% Usage:
%   llh=ecf2llhT(ecf)
%
% Inputs
%   ecf   - 3xN Matrix of Vectors in ECF coordinates (meters)
%           where
%              x(1,:) is Greenwich meridan; (0 lon, 0 lat)
%              y(2,:) is Easterly
%              z(3,:) is North Pole
%
% Outputs:
%   llh   - 3xN Matrix of Vectors
%           where
%              llh(1,:) is latitude positive (radians)
%              llh(2,:) is longitude positive East (radians)
%              llh(3,:) is height (meters)
%


% So, if wanted to calculate elevation angles (and therefore visibility of
% MODIS then would just need to calculate the llh of the orbit above and do
% the following for each lat lon
sat_tcs=llh2tcsT(sat_llh,origin_llh);  %llh to tcs at origin_llh
%function [tcs]=llh2tcsT(llh,origin);
%
% Description:
% 	This function returns the x, y, and z topocentric (TCS)
%	coordinates of the point specified by llh [lat lon hgt],
%   relative to the input origin [lat lon alt].


%The satellite elevation as viewed from the origin (ground location)
%I.e., presumably the angle between the surface and the satellite
sat_elev=atan2(sat_tcs(3,:),sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2));
% -------------------------------------------------------------------------
% *** if also just want elevation then can stop here ***
% -------------------------------------------------------------------------


%Identify visible segments: 
%notVIS=find(sat_tcs(3,:)<0); %i.e. when z is < 0 it is not visible - also the elevation
%will < 0 for these points. I guess that if we are above sea level then we
%might be able to see the satellite with a negative elevation angle -
%although I think that the presumption is that this will not be the case

%all the points that didn't satisfy the above
%VIS=setdiff([1:npts],notVIS);
%SETDIFF Set difference.
%    SETDIFF(A,B) when A and B are vectors returns the values
%    in A that are not in B.  The result will be sorted.  A and B
%    can be cell arrays of strings.

VIS = find(sat_tcs(3,:)>=0);

%sets all the lat and lons and tcs of non-visible points to NaN - may not
%want to do this
%sat_llh(:,notVIS)=NaN;
%sat_tcs(:,notVIS)=NaN;

%Extract start and end times for visible segment:
%When the index jumps by more than one we are at the start of a new
%consequetive segment
iend = find(diff(VIS)>1);
t_end= tsince(VIS(iend));  iends=VIS(iend);
t_start = tsince(VIS(iend+1)); istarts=VIS(iend+1);
t_start = [tsince(VIS(1)) t_start]; istarts=[VIS(1) istarts];
t_end =[t_end tsince(VIS(end))]; iends =[iends VIS(end)];


%t_start=tsince(notVIS(diff(notVIS)>1));  istarts=notVIS(diff(notVIS)>1);
%t_end  =tsince(VIS(diff(VIS)>1));        iends = VIS(diff(VIS)>1);
% if length(t_end)<length(t_start)
%    %t_end  =[t_end,tsince(VIS(end))];     istarts  =[istarts,VIS(end)];
% elseif length(t_start)<length(t_end)
%   % t_start=[tsince(notVIS(1)),t_end];    iends=[notVIS(1),iends];
% end
%addition by Dan - if the sat is visible at the start then make this a
%start point


% t_end = tsince(VIS(istart-1));
% if VIS(1)==1
%     t_start=[tsince(1) t_start]; istarts=[1 istarts];
% end
% %Also add the last end point
% t_end  =[t_end,tsince(VIS(end))];  iends=[iends,VIS(end)];

%t_mid  =(t_start+t_end)/2;

figure
plot(tsince(:)/min_per_day,sat_elev(:)/dtr,'r')
hold on
plot(t_start/min_per_day,zeros(length(t_start)),'b^')
hold on
plot(t_end/min_per_day,zeros(length(t_end)),'b^')






clear max_elev lat_max_elev d aob sens_angle
fprintf('PASS#  N START/END: MONTH DAY HOUR:MIN\n');   
for npass=1:length(t_start)
    inds = istarts(npass):iends(npass);
    [max_elev(npass),imax]=max(sat_elev(inds));
    max_elev(npass)=max_elev(npass)/dtr; %convert to degrees
    lat_max_elev(npass) = sat_llh(1,inds(imax))/dtr;
    lon_max_elev(npass) = sat_llh(2,inds(imax))/dtr;
    [d(npass),aob(npass)]=distlatlon(lct.station.rx_latitude,lct.station.rx_longitude,lat_max_elev(npass),lon_max_elev(npass));
    sens_angle(npass) = 90 - aob(npass) - max_elev(npass);
    
   
    
    
    hold on
    npass_str=num2str(npass);
    text(t_end(npass)/min_per_day,10,npass_str)
    [Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays+t_start(npass)/min_per_day);    
    fprintf('PASS#%3i START: %5i %2i %4i:%2.2d ',npass,Emon,Eday,Ehr,Emin);    
    [Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays+t_end(npass)/min_per_day);
    fprintf('END: %5i %2i %4i:%2.2d ',Emon,Eday,Ehr,Emin);
    fprintf('MAX_ELEV, SENS ANGLE: %2.2f %2.2f\n',max_elev(npass),sens_angle(npass));
    if t_start(npass)==tsince(1)
        fprintf(' -- N.B. above is START of requested period (not elev=0) --\n');
    end
    
   
end



grid on
xlabel('UT--days')
ylabel('elevation--deg')
title(['Satellite ID: ',num2str(satrec.satnum),'--Station: ',rx_name])







%  ****** nPass and 300km value set here ******
%Select Pass for summary
nPass=2;
h_intercept=300000;   %300 km

SatID.satnum=satrec.satnum;
SatID.year=Eyear;
SatID.mon =Emon;
SatID.day =Eday;
plotTitle=['Object ',num2str(satrec.satnum),' Pass ',num2str(nPass)];

t_start=t_start(nPass);
t_end=  t_end(nPass);
npts=ceil((t_end-t_start)*60);
tsince_min=linspace(t_start,t_end,npts);

%get the satellite geometry for the a single pass

satGEOM_struct=satGEOM(satrec,Eyear,tsince_min,origin_llh,h_intercept,plotTitle);
%function  satGEOM_struct=satGEOM(satrec,Eyear,tsince_min,origin_llh,h_intercept,PassID)
    %Generate Satellite Geometry
    %USAGE:
    %
    %INPUTS:
    %   satrec      = structure generated by twoline2rvMOD
    %   Eyear       = TLE Epoch fractional years 
    %   tsince_M    = time since epoch in minutes
    %   origin_llh  = llh column vector [lat_deg, lon_deg, altitude_met]
    %   PassID      = ID string 
    %
    %OUTPUTS:
    %*************ecf Coordinates**********************************************************
    %xsat_ecf, vsat_ecf = satellite state vector in ecf coordinates from spg4         (3XN)
    %************GPS Coordinates***********************************************************
    %sat_llh            = satellite geodetic coordinates                              (3XN)
    %************Station TCS coordinates***************************************************
    %sat_tcs,  vsat_tcs = satellite state vector in receiver tcs system at origin_llh (3XN)
    %sat_rng,  sat_rdot = satellite rang & range rate (=> sat_tcs)                    (3XN)           (NX1)
    %sat_elev, sat_phi  = satellite elevation & true bearing                          (NX1)
    %*************Propagation Reference Coordinates at penetration point*******************
    %xyzp               = propagation coordinates with origin at h_intercept          (3xN)
    %thetap, phip       = polar angles wrt x  (phip  cw from y-axis                   (NX1)
    %                                          theta cw from x-axis
    %rngp               = range from receiver to intercept point                      (Nx1)
    %uk_xyzp            = unit vector pointing along propagation direction            (3XN)
    %s                  = unit magnetic field vector xp,yp,zp system                  (3xN)
    %thetaB,psiB        = polar angles wrt xp                                         (Nx1)
    %vp                 = penetration point velocity <= satellite motion
    %vk                 = apparent velocity in measurement plane
    
    

% -------------------------------------------------------------------------
% ***          STOPPING here        ***
% -------------------------------------------------------------------------
return


%SummarizeSatelliteGeometry(satGEOM_struct,SatID,origin_llh,...
%                                      tsince_min,plotTitle);
