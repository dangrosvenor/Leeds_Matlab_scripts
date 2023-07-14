function lst_out = LST(lon,date)
%function lst_out = LST(lon,date)

[Y,M,D] = datevec(date);

%See Wikipedia entry for the method (I think). Or in OneNote
lstm = 15*lon/360*24; %the LST standard meridion
d =-( datenum(date)-datenum(Y,1,1) ); %number of days since start of year 
B=360/365*(d-81);
eot=9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B); % "Equation of time"
tc = (4*(lstm-lon)+eot)/60; 
lst_out = lon/360*24 + tc;