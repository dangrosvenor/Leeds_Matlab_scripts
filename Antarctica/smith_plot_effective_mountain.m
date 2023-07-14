po=291;
%po=289;
po=292.4;  %one that works for H=3.2 km
po=292.0;  %one that works for H=3.2 km
%po=290.9;
%po=288;
po=291; %perhaps the best one
po=292;
%po=294; %for equiv. potemp

po=283;

U=8;
gd=0.63;
gd=0.55;
gd=0.4;

data_type='cross_section';
data_type='XY_pot_cross_data';

switch data_type
    case 'XY_pot_cross_data'
        pot_cross = pot_cross_6th_Jan_6UTC;
        Y=XY_pot_cross_data.Y_cross*1000;
%        X=XY_pot_cross_data.X_cross(1:pend)*1000;
        X=XY_pot_cross_data.X_cross*1000;  %still need to multiply by 1000 if is longitude
        potemp=pot_cross(:,1);
        
    case 'cross_section'        
        Y=zz(1).z'*1000;
        X=d_line*1000;
        potemp=pdat(1).p(:,1);
end

ipo=findheight_nearest(potemp,po);
H=Y(ipo)/1000;
po=potemp(ipo);

[cA.c cB.c]=contour(X/1000,Y/1000,pot_cross,[po po],'b'); %x,y data for the 291 K potemp contour (run plot.. first
    %for vertical cross section
    
    
    
%HM=1.5; %mountain height in km
dx=1e3; %spacing in x
X3=[X(1):dx:X(end)];
terr3=interp1(X,terr_slice,X3);
[HM,ihm]=max(terr3/1000);


%cA contains the useful data cA(1).c(1,1) is the contour value (291 K here). N=cA(1).c(2,1) is then the number of x,y points that follow
%cA(1).c(1,2:N) will then be all the x values and (2,2:N) all the y ones

X_291=cA(1).c(1,2:end)*1000;
Y_291=cA(1).c(2,2:end)*1000;

y_lim = 3.3e3;
x_lim = 800e3;
x_lim=-64;

istream=find(Y_291<x_lim & X_291<y_lim);  %limit the region of contour x,y values to reject the ones we don't want
X2_291=X_291(istream);
Y2_291=Y_291(istream);
[X2_291 isort]=sort(X2_291);
Y2_291=Y2_291(isort);

[X2_291,I,J]=unique(X2_291);
Y2_291=Y2_291(I);


%X3_291 = [X2_291(1):dx/1000:X2_291(end)];
Y3_291 = interp1(X2_291,Y2_291,X3);  %interpolate to the finer grid of dx=1 km (like X3)
inan = isnan(Y3_291);
Y3_291(inan==1)=[];
X3_291=X3;
X3_291(inan==1)=[];

%now try to match those to a h value
fprintf(1,'\n******  NOTE - need to set L value *******\n');
L=0.0016; %set L
L=0.002; %set L
L=LS;


del_aim=(Y3_291-H*1000);  %displacement from the chosen potential temp contour

%find location of the mountain peak in X3_291
ipeak=findheight_nearest(X3_291,X3(ihm));
del_peak=del_aim(ipeak);
[hhat_peak,H0_hat_peak]=mountain_smith_solve_hhat_for_given_critical_del(del_peak*L);





%NL=6.77;
NL=H0_hat_peak*6/pi
H0=NL*pi/6/L; %set H0

clear h heff2
for i=1:length(del_aim);
    heff2(i)=mountain_smith_solve_fun_h_for_given_del(del_aim(i),L,H0);
end


heff2=interp1(X3_291,heff2,X3);
heff2(isnan(heff2))=0;
heff=HM-max(heff2)/1000+heff2/1000+0.2; %match to fit to the top of the mountain



%%%%%%%%%%%%%%%%%% Houghton
%K4=uA*hA;
%thi_crit=K4^(2/3)*gd^(-1/3); %critical value at mountain crest

thi_crit = H*1000 - HM*1000 + del_peak;
u_crit = sqrt(gd*thi_crit); 
K4 = u_crit*thi_crit;
H0_Houghton = K4/U;
    %if know del at the mountain top then we know thi_crit and so u_crit from u_crit=sqrt(g*thi_crit)
    %then can calculate K4=u_crit*thi_crit, which gives hA, if we know uA, from uA*hA=K4
    %Want to try and get heff from the displacement of the streamline and also consider
    %the displacement at the mountain crest to calculated hA?

F0=U/sqrt(gd*H0_Houghton);
hpeak_Houghton = H0_Houghton .* ( 1 + 0.5*F0.^2 - 1.5*F0.^(2/3) );
off_Houghton = HM - hpeak_Houghton/1000;

heff_Houghton=zeros([1 length(del_aim)]);
%i0=find(del_aim(1:ipeak))>=0;
%heff_Houghton(i0)=0;
inds=find(del_aim<-2);

iwind=find(inds<ipeak-10);
for i=inds(iwind)    
  heff_Houghton(i)=Houghton_find_heff_from_streamline_del(del_aim(i),U,gd,H0_Houghton,hpeak_Houghton,'windward');  
end

heff_Houghton(ipeak-1:ipeak) = hpeak_Houghton;

ilee=find(inds>ipeak);
for i=inds(ilee)
  heff_Houghton(i)=Houghton_find_heff_from_streamline_del(del_aim(i),U,gd,H0_Houghton,hpeak_Houghton,'lee');
end


%%%%



'Finished effective mountain calculation'