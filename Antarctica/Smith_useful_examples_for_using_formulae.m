
%to find the critical mountain height for a given Ha & Hb as well as delA and delB at the peak
[hpeak,dA,dB]=Smith_dual_layer_find_hpeak(Ha*LS,Hb*LS,LS);

%to find the delA and delB values at mountain height lower than hpeak
[delA,delB]=Smith_solve_dual_layer_solve_del(800*LS,1900*LS,500*LS);

%continuously stratified
z_surf=1100;
[hmax,del] = mountain_smith_solve_hhat_for_given_critical_H0_hat((H0-zsurf)*LS,LS)

del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
H0_hat_strat=hhat - del_hat + acos(hhat./del_hat);

%solve for del up and down the mountain, as well as for U etc.
[h,hx,hx2,d,d2,ddel_dz,Ahat,Bhat,U]=mountain_smith_solve(LS,hhat,3.1*pi/6/LS);


%this script loads in the cross section data
load_cross_section_data_from_streamlines.m

%average windspeed
ipos=1;
zz(1).z = XY_pot_cross_data.Y_cross;
iz_find=findheight_nearest(zz(1).z,3.7);
%istart=24; %1.8 km
istart=22; %1.6 km
%iend=29; %2.3 km
%iend=35; %2.9 km
iend=iz_find;

%weighted average - BUT zz(1).z is has a regula z grid
%U_mean = sum(squeeze(zz(1).z(istart+1:iend,ipos))'.*diff(zz(1).z(istart:iend)))./sum(diff(zz(1).z(istart:iend)));
%so can just do... as includes all the data points
U_mean = mean(U_cross_6th_Jan_6UTC(istart:iend,ipos));


%% calculating N
ipos=1;
%istart=29;
%iend=35;
N=calc_N(pot_cross_6th_Jan_6UTC(istart,ipos),pot_cross_6th_Jan_6UTC(iend,ipos),XY_pot_cross_data(1).Y_cross(istart)*1000,XY_pot_cross_data(1).Y_cross(iend)*1000)
%N=calc_N(pot_cross(istart,ipos),pot_cross(iend,ipos),zz(1).z(istart)*1000,zz(1).z(iend)*1000)



%setting effective surfaces and z values, etc.
%3710 gives the best position if use the wind contours (U=9.57 m/s)
z=1000;za=1800;zb=3710;H0=zb-z
z=1000;za=1600;zb=3740;H0=zb-z
z=1000;za=1600;zb=3850;H0=zb-z;hm=640
z=1000;za=1600;zb=3149;H0=zb-z;hm=370
z=1257;za=1600;zb=3550;H0=zb-z;hm=1640-z  


ipos=findheight_nearest(XY_pot_cross_data(1).X_cross,-68.0461);
z=1260;za=1600;zb=3550;H0=zb-z;hm=1640-z  %when using ipos=1
z=1260;za=1600;zb=3500;H0=zb-z;hm=1640-z  %when using ipos for -68.0461

istart=findheight_nearest(zz(1).z,za/1000);
iend=findheight_nearest(zz(1).z,zb/1000);
N=calc_N(pot_cross_6th_Jan_6UTC(istart,ipos),pot_cross_6th_Jan_6UTC(iend,ipos),XY_pot_cross_data(1).Y_cross(istart)*1000,XY_pot_cross_data(1).Y_cross(iend)*1000)
U_mean = mean(U_cross_6th_Jan_6UTC(istart:iend,ipos));
LS=N/U_mean;

[hpeak,dA,dB]=Smith_dual_layer_find_hpeak((za-z)*LS,H0*LS,LS)
%hm=640;
percent_error=100*(hpeak-hm)/hm   %this should be consistent with hpeak
%can work out a percentage error:
%100*(1445-z)/hpeak-100 %positive means heff>hpeak  (heff=Hmountain-z)


%hpeak=630;
%calculate delA and delB for the whole approach up the mountain and back down
%mountain needs to be the critical mountain height for Hb as calculated above
[delA,delB,hx,delA2,delB2,hx2,U,U2]=Smith_solve_dual_layer_del_for_h_up_and_down2((za-z)*LS,H0*LS,hpeak*LS,U_mean,1)


%terr_cross2 is calculated from the load_cross_section_data_from_streamlines.m
%using the locations of the NaNs - perhaps not the best measure of terrain
%or use streamline as contour - use potemp_contour_height.m

terrain_type='real';
terrain_type='theta';

switch terrain_type
    case 'real'
        terr_del = terr_cross2;
        X_del=XY_pot_cross_data.X_cross;
    case 'theta'
        terr_del = Y2_pot; %after running potemp_contour_height.m
        X_del=X2_pot;
end


%potemp_contour_height also plots the effective surface contour
%in bold
%to plot the height of the streamline use the iadd_streamline_z0 option
%in plotTimeHeightVap3. Can also plot the effective surface there.

[maxp ipeak]=max(terr_del);
terr_eff=terr_del*1000-z;
terr_eff(1:ipeak) = max(0,terr_eff(1:ipeak));

clear delA_x delB_x
for i=1:ipeak
    delA_x(i) = interp1(hx/LS,delA,terr_eff(i)); %MAKE SURE z is SET CORRECTLY!!
    delB_x(i) = interp1(hx/LS,delB,terr_eff(i));
end

for i=ipeak+1:length(terr_eff)
    delA_x(i) = interp1(hx2/LS,delA2,terr_eff(i)); %MAKE SURE z is SET CORRECTLY!!
    delB_x(i) = interp1(hx2/LS,delB2,terr_eff(i));    
end


%remove the NaN values
inan=find(isnan(delB_x)==1);
delB_x2=delB_x;
delB_x2(inan)=[];
delA_x2=delA_x;
delA_x2(inan)=[];
X2_pot2=X_del;
X2_pot2(inan)=[];
terr_eff2=terr_eff;
terr_eff2(inan)=[];
%plotting the del data as a streamline
%streamline as terrain
%plot(X2_pot2,(delA_x2/LS+za)/1000,'wx-');
%plot(X2_pot2,(delB_x2/LS+zb)/1000,'wx-');

plot(X2_pot2,(delA_x2/LS+za)/1000,'wo-','markersize',10,'linewidth',1.5);
plot(X2_pot2,(delB_x2/LS+zb)/1000,'wo-','markersize',10,'linewidth',1.5);
ii=find(terr_eff2==0); %want to keep the potential temp contour 
%so only plot up to where the new effective surface meets it - otherwise it
%makes a solid rather than dotted line
plot(X2_pot2(ii),(terr_eff2(ii)+z)/1000,'w--','linewidth',2);


%actual terrain
%plot(XY_pot_cross_data.X_cross,(delA_x/LS+za)/1000,'wo');
%plot(XY_pot_cross_data.X_cross,(delB_x/LS+zb)/1000,'kx-');

%grid_xy=meshgrid(XY_pot_cross_data.X_cross,XY_pot_cross_data.Y_cross);
UA_wrf = interp2(XY_pot_cross_data.X_cross,XY_pot_cross_data.Y_cross,U_cross_6th_Jan_6UTC,X2_pot2, (delA_x2/LS+za)/1000);

Ha_hat = (za-z)*LS;
Hb_hat = (zb-z)*LS;
%calculating the speed of streamline delA2
%note, doesn't matter that are using values*LS for AA2 
%since LS values will all cancel. Still gives A such that del_1=Az+B
A2=(delA2-hx2)./((za-z)*LS + delA2 - hx2);
AA2=(delA_x2-terr_eff2*LS)./((za-z)*LS + delA_x2 - terr_eff2*LS);
U0=7; %what to choose?
%U0=U_mean;
Udown = (1-A2)*U0;
Udown2 = (1-AA2)*U0;

%find C_hat and D_hat from eqn set (22) in Smith and Sun - these are independent of z
%first solve using the 5th and 6th equations
D_hat = AA2 ./ ( cos(Ha_hat+delA_x2) - cot(Hb_hat+delB_x2).*sin(Ha_hat+delA_x2) );
C_hat = D_hat .* cot(Hb_hat+delB_x2);
%BB2 = terr_eff2*LS .* (1-AA2);
BB2 = terr_eff2 .* (1-AA2); %there is a mistake in Smith and Sun eqn 22 for the 
%1st 2 equations : should be B_hat rather than B - so remove LS from the previous
%version to get the "non-hatted" B

%test by solving using the 3rd and 4th equations - get the same answer
XA=Ha_hat+delA_x2;
D_hat2 = ( AA2 + delA_x2.*tan(XA) ) ./ ( tan(XA).*sin(XA) + cos(XA) );
C_hat2 = ( delA_x2 - D_hat2.*sin(XA) ) ./ cos(XA);


%trying again...
XA=Ha_hat+delA_x2;
XB=Hb_hat+delB_x2;
C_hat = (delB_x2.*cos(XA) - AA2.*sin(XB) ) ./ (cos(XA).*cos(XB) + sin(XA).*sin(XB)) ;
D_hat = C_hat .* tan(XB);



%so now can calculate del for any position by knowing h(x) and z 
%(and delA(x) and delB(x), A(x), B(x) and C(x) and D(x))
%from del_1(x,z) = LS*(Az + B) for the lower layer
%and  del_2(x,z) = C_hat*cos(Lz) + D_hat*sin(LZ) (dels are hatted here)
%for each x-position can work out a column of dels from these equations
%will have to look out for where del_1 meets del_2 ; where del_1/dz=del_2/dz so know which to apply
dz=10;
Hb=zb-z;
z_d = [-Hb:dz:Hb]; %include negative values for negative mountain (downslope)
z_d_arr = repmat(z_d,[length(AA2) 1]); %creates 2_d array dimensions [X,Z]

%calculate dels for the lower layer (those that start below Ha at x=0)
A_arr = repmat(AA2,[length(z_d) 1])'; %creates 2_d array dimensions [X,Z]
B_arr = repmat(BB2,[length(z_d) 1])'; %creates 2_d array dimensions [X,Z]
del_1 = LS*(A_arr.*(z_d_arr) + B_arr); %calculate del_hat
del_1_dz = LS*A_arr; %calculate del_hat/dz

%calculate dels for the upper layer (those that start between Ha and Hb at x=0)
C_arr = repmat(C_hat,[length(z_d) 1])'; %creates 2_d array dimensions [X,Z]
D_arr = repmat(D_hat,[length(z_d) 1])'; %creates 2_d array dimensions [X,Z]
del_2 = C_arr.*cos(LS*z_d_arr) + D_arr.*sin(LS*z_d_arr); %calculate del_hat
del_2_dz = -LS*C_arr.*sin(LS*z_d_arr) + LS*D_arr.*cos(LS*z_d_arr); %calculate del_hat/dz

clear temp iz del_all
for i=1:size(del_1,1)
    z_a=(za-z)+delA_x2(i)/LS; %height that the lower A streamline should be at
    [temp(i),iz(i)]=min(abs(z_a-z_d_arr(i,:)));
    lower = del_1(i,1:iz(i));
    upper = del_2(i,iz(i)+1:end);
    del_all(i,1:length(lower)) = lower;
    del_all(i,length(lower)+1:length(lower)+length(upper)) = upper;
end

z0=2650*LS; %the z to draw a streamline for
[z_stream]=find_streamline_from_dels_Smith_and_Sun(z0,del_all,z_d_arr*LS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculating wind speedup factor such that U=U0*Ufac  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ufac2=(1-del_2_dz/LS);
%2D color plot of lon-height for 1-ddel/dz
figure
pcolor(X2_pot2,z_d,Ufac2');
shading interp
colorbar
  %plot the delA and delB streamlines
plot(X2_pot2,za+delA_x2/LS-z,'k')
plot(X2_pot2,zb+delB_x2/LS-z,'k')
  %plot the maximum windspeed line for the flow below the dividing streamline
ilim=findheight(z_d,1000);
[temp,imax]=max(Ufac2(:,1:330),[],2);
plot(X2_pot2,z_d(imax),'k--');
  %plot the terrain
plot(X2_pot2,terr_eff2,'k','linewidth',2)
set(gca,'tickdir','out')
set(gca,'clim',[0 5.5]);
contour(X2_pot2,z_d,Ufac2',[1 1],'w'); %plot on the Ufac=1 line
xlabel('Longitude');
ylabel('Height from z=1260 m');



%U profiles at LON = -66
ii=findheight_nearest(X2_pot2,-66);
ii2=findheight_nearest(XY_pot_cross_data.X_cross,-66);
figure
U0=5;
plot(Ufac2(ii,:)*U0,(z_d+z)/1000,'r'); %plot the Smith windspeed for given U0 at -66
hold on; grid
plot(U_cross_6th_Jan_6UTC(:,ii2),XY_pot_cross_data.Y_cross); %plot model cross section profile
set(gca,'fontsize',14);
xlabel('Wind speed (m/s)','fontsize',14);
ylabel('Height (km)','fontsize',14);
title('Predicted and modelled wind speed profiles at LON=-66','fontsize',14);

%now check what happens if use the actual upstream windspeed for each streamline
z0s=z_d+z-del_2(ii,:)/LS; %the upstream intial heights for the streamlines ( since z=H0(i)+delB(i)  )
zA=(delA_x2(ii)/LS+za)/1000; %height of the lower streamline
zB=(delB_x2(ii)/LS+zb)/1000; %height of delB dividing streamline
ir=find(z_d+z>=zA*1000 & z_d+z<=zB*1000); %find all points between delA and delB streamlines
z0s2=z0s(ir); %limit upstream heights to these
U0s=interp1(XY_pot_cross_data.Y_cross,U_cross_6th_Jan_6UTC(:,1),z0s2/1000); %interpolate
  %the wind speed from x=0 profile (upstream profile)
Us=Ufac2(ii,ir).*U0s; %multiply by the speedup factor
plot(Us,(z_d(ir)+z)/1000,'g','linewidth',2);
line([-15 30],[zA zA],'color','k'); %plot horizontal lines to show the positions
line([-15 30],[zB zB],'color','k'); %of the delA and delB streamlines at LON=-66

%calculate ratio between predicted and modelled
%need to interpolate one onto the grid of the other to do this
Umod=interp1(XY_pot_cross_data.Y_cross,U_cross_6th_Jan_6UTC(:,ii),(z+z_d(ir))/1000);
figure
plot(z_d(ir)+z,Us./Umod);
grid

%Richardson number
%Ri = N^2/ (dU/dz)^2
%know that U=U0(1-ddel/dz) so dU/dz = -U0*d2del/dz2
%from (1) of Smith and Sun know that d2del/dz2 = -L^2*del for layer i=2
%so dU/dz = U0*L^2*del and
%Ri= N^2/(U0^2*[L^2*del]^2) = 1/(L^2*del^2) = 1/del_hat^2
%work out for delB2 - lowest Ri will be on z=zb
Ri=1./(delB_x2).^2;
[temp,ii]=min(abs(Ri-0.25)); %find where Ri is closest to 0.25 as below
%this value will start to get instabilities
plot(X2_pot2(ii),(delB_x2(ii)/LS+zb)/1000,'ko-','markerfacecolor','k','markersize',20);



%cont strat


for i=ipeak+1:length(terr_eff)
    delA_x(i) = interp1(hx2/LS,d,terr_eff(i)); %MAKE SURE z is SET CORRECTLY!!
end

inan=find(isnan(delA_x)==1);
delA_x2=delA_x;
delA_x2(inan)=[];
X2_pot2=X_del;
X2_pot2(inan)=[];
terr_eff2=terr_eff;
terr_eff2(inan)=[];
%plotting the del data as a streamline
%streamline as terrain
%plot(X2_pot2,(delA_x2/LS+za)/1000,'wx-');
%plot(X2_pot2,(delB_x2/LS+zb)/1000,'wx-');

plot(X2_pot2,(delA_x2/LS+z+3.1*pi/6/LS)/1000,'ko-','markersize',10,'linewidth',1.5);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  interpolating wind speed along a potential temperature contour %%%%%%%
%run potemp_contour_height.m first to get the potemp contour required

U_potemp = interp2(XY_pot_cross_data.X_cross,XY_pot_cross_data.Y_cross,U_cross_6th_Jan_6UTC,X2_pot,Y2_pot);
figure;

plot(XY_pot_cross_data.X_cross,terr_cross2,'r'); %plot a scale plot of the terrain to get an idea of where the speed occurs
%relative to it 
labs(1).l=['Terrain height'];

hold on
plot(X2_pot,U_potemp,'b');
labs(4).l=['Wind speed for ' num2str(round(po)) ' K contour'];

grid
ylabel('Wind speed (m s^{-1}) / Mountain height (km)');
xlabel('Longitude');

legend(labs(1).l,labs(2).l,labs(3).l,labs(4).l,'location','northwest');
%saved some of these graphs in
cd('C:\Documents and Settings\dan\My Documents\Antarctica_mydocs\Mountain_wave_work_Smith_etc')
