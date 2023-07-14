%ds/lon, where s=min sensor ZA
dx3 = diff(ydat(3).y)./diff(xdat(3).x);
sig=sign(dx3);
sigd=diff(sig);
%where the sign of the gradient changes from 1 (pos) to -1 neg - i.e. the
%triangle peaks in s
is=find(sigd==-2);


%ds/lon, where s=min sensor ZA
dx3 = diff(ydat(3).y(3:end))./diff(xdat(3).x(1:end-2));
lon3 = MLON(2:end-1);
sig=sign(dx3);
sigd=diff(sig);
%where the sign of the gradient changes from 1 (pos) to -1 neg - i.e. the
%triangle peaks in s
is=find(sigd==-2);



%maximum spacing between positive to negative gradient changes in ds/dlon
[maxval imax]=max(diff(lon3(is)));

peaks = lon3(is(imax)):maxval:lon3(is(imax))+360;
%peaks(peaks>180)=peaks(peaks>180)-360;
peaks(peaks>180)=[];


%max value of minSZA
[maxval2 imax2]=max(ydat(1).y);
overlap_boundaryR = MLON(imax2);

iloc=find(peaks<MLON(imax2));
peak_match = peaks(iloc(end));

overlap_boundaryL = peak_match - (overlap_boundaryR - peak_match);


