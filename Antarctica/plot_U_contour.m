Uplot=9.57; %the best one for a dividing line for U contours
Uplot=9.74;
Uplot=8.75;
Uplot=8;
Uplot=9.5;

data_type='cross_section';
data_type='XY_pot_cross_data';

switch data_type
    case 'XY_pot_cross_data'
        U_cross = U_cross_6th_Jan_6UTC;
        Y=XY_pot_cross_data.Y_cross;
%        X=XY_pot_cross_data.X_cross(1:pend)*1000;
        X=XY_pot_cross_data.X_cross;  %
        Udat=U_cross(:,1);
        
    case 'cross_section'        
        Y=zz(1).z'*1000;
        X=timesTH(1).t(1:pend);
        Udat=pdat(1).p(1,1:pend);
end

ipo=findheight_nearest(Udat,Uplot);
%Uplot=Udat(ipo);

[cA.c cB.c]=contour(X,Y,U_cross,[Uplot Uplot],'b','linewidth',3); %x,y data for the 291 K potemp contour (run plot.. first
    %for vertical cross section
       

%cA contains the useful data cA(1).c(1,1) is the contour value (291 K here). N=cA(1).c(2,1) is then the number of x,y points that follow
%cA(1).c(1,2:N) will then be all the x values and (2,2:N) all the y ones

% X_pot=cA(1).c(1,2:end);
% Y_pot=cA(1).c(2,2:end);
% 
% y_lim = 3.3; %km
% x_lim = 800;
% x_lim=-63; %longitude
% 
% istream=find(Y_pot<y_lim & X_pot<x_lim);  %limit the region of contour x,y values to reject the ones we don't want
% X2_pot=X_pot(istream);
% Y2_pot=Y_pot(istream);
% [X2_pot isort]=sort(X2_pot);
% Y2_pot=Y2_pot(isort);
% 
% [X2_pot,I_pot,J]=unique(X2_pot);
% Y2_pot=Y2_pot(I_pot);
