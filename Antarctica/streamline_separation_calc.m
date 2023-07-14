po=291;
%po=289;
po=292.4;  %one that works for H=3.2 km
po=292.0;  %one that works for H=3.2 km
%po=290.9;
%po=288;
po=291; %perhaps the best one
po=292;
%po=294; %for equiv. potemp

po=292;
po2=293;

%po=293;
%po2=294;

Y=zz(1).z;
X=timesTH(1).t;

np=100;
dx=( LATS(end)-X(1) )/np; %spacing in x
X3=[X(1):dx:LATS(end)];

X3=LATS;




ipo=findheight_nearest(potemp,po);
po=potemp(ipo);

[cA.c cB.c]=contour(timesTH(1).t,zz(1).z,pdat(1).p,[po po],'b'); %x,y data for the 291 K potemp contour (run plot.. first
[cA2.c cB2.c]=contour(timesTH(1).t,zz(1).z,pdat(1).p,[po2 po2],'b'); %x,y data for the 291 K potemp contour (run plot.. first
    %for vertical cross section
    
 
%cA contains the useful data cA(1).c(1,1) is the contour value (291 K here). N=cA(1).c(2,1) is then the number of x,y points that follow
%cA(1).c(1,2:N) will then be all the x values and (2,2:N) all the y ones

X_po=cA(1).c(1,2:end);
Y_po=cA(1).c(2,2:end);

X_po2=cA2(1).c(1,2:end);
Y_po2=cA2(1).c(2,2:end);



%discard some contour areas
%istream=find(Y_po<3.3e3 & X_po<800e3);  %limit the region of contour x,y values to reject the ones we don't want
istream=[1:length(X_po)];
X2_po=X_po(istream);
Y2_po=Y_po(istream);
%sort
[X2_po isort]=sort(X2_po);
Y2_po=Y2_po(isort);

%istream2=find(Y_po2<3.3e3 & X_po2<800e3);  %limit the region of contour x,y values to reject the ones we don't want
istream2=[1:length(X_po2)];
X2_po2=X_po2(istream2);
Y2_po2=Y_po2(istream2);

[X2_po2 isort]=sort(X2_po2);
Y2_po2=Y2_po2(isort);

%make unique
[X2_po,I,J]=unique(X2_po);
Y2_po=Y2_po(I);

[X2_po2,I,J]=unique(X2_po2);
Y2_po2=Y2_po2(I);


%X3_po = [X2_po(1):dx/1000:X2_po(end)];
Y3_po = interp1(X2_po,Y2_po,X3);  %interpolate to the finer grid of dx=1 km (like X3)
inan = isnan(Y3_po);
Y3_po(inan==1)=[];
X3_po=X3;
X3_po(inan==1)=[];

Y3_po2 = interp1(X2_po2,Y2_po2,X3);  %interpolate to the finer grid of dx=1 km (like X3)
inan = isnan(Y3_po2);
Y3_po2(inan==1)=[];
X3_po2=X3;
X3_po2(inan==1)=[];

diff_stream = Y3_po2 - Y3_po;


fprintf(1,'Finished calc streamline diff\n');