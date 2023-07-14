

function solang=zenith(doy,TIME,LATin,longitude)
%  function solang=zenith(doy,TIME,LATin,longitude)
%  This routine takes the date, time and latitude and
%  returns the relevant solar zenith angle.
%

%  TIME      = Universal time (decimal hours)
%  LATin       = Latitude (degrees north)
%  longitude = Longitude (degrees east)
%  solang    = (Solar Zenith Angle)
%
%                      A. Collard.  Univ. of Wisconsin
%
PI = 3.141592654;
TWOPI = 2*PI;
%
%  The following is an array of the day of year of the first day
%  of each month (365 day year) minus 1.
%
doy_mm=[0,31,59,90,120,151,181,212,243,273,304,334];
%
% Calculate the day of year from YYMMDD
%

%
%  Calculate the solar declination
%

longitude=-longitude;

GAM=2.*PI*(doy-1)/365.;
dec=0.006918-0.399912*cos(GAM)+0.070257*sin(GAM)-...
0.006758*cos(2.*GAM)+0.000907*sin(2.*GAM)-...
0.002697*cos(3.*GAM)+0.00148*sin(3.*GAM);
%
%  Calculate the equation of time
%  True Solar Time = Mean Solar Time + eqntim
%  From "Spherical Astronomy" by W.M. Smart
%
SOL_LONG = (doy-80.)*0.017202;
eqntim = -103.9*sin(SOL_LONG) - 429.6*cos(SOL_LONG) + ...
596.3*sin(2.*SOL_LONG) - 2.0*cos(2.*SOL_LONG) + ...
4.3*sin(3.*SOL_LONG) + 19.3*cos(3*SOL_LONG) - 12.7*cos(4.*SOL_LONG);
% Convert to hours
eqntim = eqntim/3600.;
%
%  Find local time from UT and longitude
%
TIM = TIME - longitude/15.;
%
%  Calculate solar hour angle
%
hang = (TIM - 12. + eqntim)*TWOPI/24.;
%
%  Convert latitude to radians
%
LAT = LATin*0.017453293;
%
%  Calculate the cosine of the solar zenith angle.
%
solang=sin(dec).*sin(LAT)+cos(dec).*cos(LAT).*cos(hang);
solang=(180./pi).*acos(solang);


