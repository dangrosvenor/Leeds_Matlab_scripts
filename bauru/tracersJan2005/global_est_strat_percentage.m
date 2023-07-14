%global estimate of moistening based on the percentage of stratospheric flux of air

%area of globe from -L to +L latitude
L=15; %latitude for tropics (degrees)
R=6400e3; %radius of the Earth (m)

A=4*pi*R^2*sin(pi*L/180); %( integral from -L to +L of dA=2*pi*R'*R*dthi
                   % where R' is the radius of circle based on rotation axis
                   % at lat thi, R'=Rcos(thi)
                   % dA=circumference of this cicle * arc length R*dthi )
                   
%volume of air entering the strat per year at uplift rate of W m per month
W=1e3; %mean stratospheric uplift (m per month)
V=A*12*W; %volume in one year (12 months) (m^3)

%vapour entering the stratosphere in one year
rho=0.2; %mean density of air at tropopause (kg/m^3)
qvap=3.8; %mean mixing ratio of air at trop (ppmv)
%% should actually use 3.8 ppmv acc. to Dessler 1998

f=1e6*28.97/18; %conversion factor for kg/kg to ppmv
Mvap = V * rho * qvap/f %mass of vapour into strat in one year (kg)


%NOTE - Holton (1995) quotes an annual average mass flux across 100 hPa of 85e8 kg/s (for L=-15 to +15)
%This is based on Rosenlof and Holton (1993), based on UKMO measurements (radiosondes I think)
%The estimate here for this lat range at w=1e3 gives 100e8 kg/s - so pretty close
M=85e8 * 24*365*3600; %total mass of air crossing 100 hPa in a year based on Rosenlof and Holton (1995) 

Mvap2=M*qvap/f %total water vapour mass in a year (kg) based on above and qvap ppmv mean entry mixing ratio

Mvap2=1e9*365; %Thomas Peter's quoted value (1e9 kg per day) - is lower than the one above (1.73e9 kg per day)


% Mvap2 = 6.6e11 kg/year based on the above 85e8 kg/s and 4 ppmv stratospheric entry mixing ratio. 
% This is compared with 7.9e11 kg/year estimated from the mean uplift of 1 km/month, rho=0.2 L=15, qvap=4 ppmv.
% 
% Holton (1995) exercises a word of caution about using the 100 hPa surface though but should be suitable for this purpose.
% 
% Now just remains to compare this with LEM flux of vapour/ice into stratosphere.


%so are going to use the total amount of mass (vapour or total water) transported into the stratosphere
%in the various model runs. total will just then be this multiplied by the number of events per year
%estimated from Bauru radar and number of regions likely to be similar. Plus can use satellite stats
%such as the TRMM 20 dBZ data

%dndt is the number per m^2 per month
dndt=(0.29+0.22+0.23)/100*10198 / (pi*(240e3)^2) / (51/30); %10 dbz > 18 km - 3D + 3D high ccn
%dndt = ((1.94+0.48+0.23+0.22+0.29))/100 * 10198  / (pi*(240e3)^2) / (51/30);  %no. 10 dbz echotop clouds reaching >16 km per area per month (as in 3d-weak)
%dndt = ((0.23+0.22+0.29))/100 * 10198  / (pi*(240e3)^2) / (51/30); %no. 10 dbz echotop clouds reaching >17 km per area per month (as in 3D-med)

A=1000e3*1000e3; %total area of convection in the tropics km^2 (*** need a better estimate of this!)

%dndt_tot=A*dndt*12; %to give the total number per year to compare with the slow uplift estimate

Nmonths=59;

%dndt_tot=9641*12/Nmonths; %from Liu and Zipser paper - they observed 9641 events that penetrated the tropopause over ~46 months (-20 to +20)
%dndt_tot=71989*12/Nmonths; %as above but for > 14 km
                    % perhaps the 20 dbZ threshold is too high? Perhaps not useful to look at 20 dBZ from the model
                    % given the problems with the higher dbz reflectivities. But radar data may give an idea of ratio of 20 dBZ storms above 
                    %tropopause to the number of 10 dBZ ones. But don't have any stats for 20 dBZ.
%dndt_tot=5512*12/Nmonths; %same paper - above 380K
%actually looks like Liu and Zipser stats are over 58 or 59 months (says 5 years - Jan 98 to Nov 2000 = 34 or 35 months 
%depending on it's up to the start or end of November (prob end I would guess) and then from December 2001 to Dec 2003 
%(=24 months) giving 59 months in total. Not sure where 46 months came from?

% NEW ESTIMATE - the value used above would be a large underestimate of the TOTAL number of overshooting events over the 5 years
% since the number recorded is not a sample from all locations at all times just a snapshot of what the satellite happened to see on its way
% around - it only sampled each region for about 80 seconds every 11 hours and so we can only build up statistics based on these
% samples. Overshoots are likely to not move much during 80 seconds and so it is unlikely to record any new events during that time
% But since we only have a snapshot we need to assume a lifetime for the events to extrapolate the total number recorded.  
% We know the total number from all samples (e.g. N=9641 globally for > trop). Can work out the total number of samples globally
% S=59/12*365*24/11 (one every 11 hours for almost 5 years). So then have an average of N/S overshoots per sample. If we assume this
% many new events occur every T hours globally (i.e. each has lifetime of T hours) then we have
% N/S/T events per hour. This gives...

N=9641;
S=59/12*365*24/11;
T=1;

Nhour=N/S/T; % = 2.46 per hour over the globe when T=1 - same ballpark as the Nhour=12 figure quoted by Thomas Peter. Seems likely that the lifetime will be
             % shorter for 20 dbZ above trop. 12/2.46 = 4.9 - so would need a lifetime of 12 mins to get to his number - seems reasonable I guess
             % In addition Thomas Peter used a different value for the slow uplift value. He quotes a value of 1e9 kg 
             % per day from large scale ascent - this approximately agrees with value that I used in the estimates (85e8*3600*24*3.8/f = 1.7e9 kg)
             % but is a little lower
             
Nhour=12.33/T; %from Thomas Peter (Thierry) figure of 12.33 events at any one time - divide this by the average lifetime to get number per hour             
Nhour=Nhour * 9641/5512;

             
% Tonnes (1e3 kg) of water            
%Vapour
% (1) 1115.753475   Normal case 
% (2) 86.435065     Weak 
% (3) 194.267679    Medium
% (4) 1303.933967   High CCN

% Total water
% (5) 1247.471101  Normal case
% (6) 87.146444    Weak
% (7) 197.309389   Medium
% (8) 1808.750781  High CCN           

%Vapour
mass(1).dat = 1115.753475e3  % Normal case 
mass(2).dat = 86.435065e3    % Weak 
mass(3).dat = 194.267679e3   % Medium
mass(4).dat = 1303.933967e3  % High CCN

% Total water
mass(5).dat = 1247.471101e3  % Normal case
mass(6).dat = 87.146444e3    % Weak
mass(7).dat = 197.309389e3   % Medium
mass(8).dat = 1808.750781e3  % High CCN    

% to get masses procedure is to:
% connect to y:/ drive
% loadvapdata with the following files for just the names (justname=1)
% loadselect=[118 112 113 133]; %strong, medium, weak and high CCN cases
% mimp dump 44 for full 3D fields of vapour, ice, snow and graupel
% then loadvapdata again for the icediagsALL for the first file (normal case)
% then again for grids for all 4
% then run watervapourMay2005 with graph=4876

             
dndt_tot=Nhour*24*365;
             
for i=1:length(mass)
    m_inject(i).dat=mass(i).dat*dndt_tot;
    frac(i).dat=100*m_inject(i).dat/Mvap2; %expressed as a percentage of the total stratospheric slow transport
    fprintf(1,'\n%d) %f',i,frac(i).dat);
end

%Using these values we get percentages of:

% T=1/5 hours, Mvap2 is as estimated from Holton paper (using 85e8 kg/s mass flux of air and 3.8 ppmv entry)
% 1) 19.012933
% 2) 1.472892
% 3) 3.310407
% 4) 22.219612
% 5) 21.257460
% 6) 1.485014
% 7) 3.362239
% 8) 30.821914

% T=1/5 hours, Mvap2 = Thomas Peter's value of 1e9 kg/day
% 1) 32.967756
% 2) 2.553942
% 3) 5.740130
% 4) 38.528024
% 5) 36.859686
% 6) 2.574962
% 7) 5.830005
% 8) 53.444113


% T=1 hour, Mvap2 is as estimated from Holton paper (using 85e8 kg/s mass flux of air and 3.8 ppmv entry)
% 1) 3.802587
% 2) 0.294578
% 3) 0.662081
% 4) 4.443922
% 5) 4.251492
% 6) 0.297003
% 7) 0.672448
% 8) 6.164383
%                   Values are much closer to those in thesis (e.g. 8.6% for high CCN case total water)
                    % factor of 1.4 too high - 40 % error

% T=1 hour, Mvap2 = Thomas Peter's value of 1e9 kg/day
% 1) 6.593551
% 2) 0.510788
% 3) 1.148026
% 4) 7.705605
% 5) 7.371937
% 6) 0.514992
% 7) 1.166001
% 8) 10.688823
                  % Again close to thesis values - factor of 1.2 too low - 20% error