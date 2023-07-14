function cost=elec_charges(days,kwh,dayA,dayB,night,kwh_night)
%function cost=elec_charges(days,kwh,dayA,dayB,night,kwh_night)

years=days/365;

ann=kwh/years;

if nargin==6
    
    cost = ( 1000*dayA + (ann-1000)*dayB )*years + kwh_night*night;
    
else
    
    cost = ( 900*dayA + (ann-900)*dayB )*years;
    
end