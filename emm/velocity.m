z=GridDan(1).Z/1000 + 0.62;

z0=3;  %cloud base for w
w0=16; %m/s for cloud base
zb=10.8; %start of braking height
zt=16.8; %cloud top
wb=20.1; %m/s at start of braking height

%linear change of wpeak from LEM
w=interp1([z0 zb],[w0 wb],z,'linear','extrap');

%ft=1 up to braking height, zb, and then reduces to 0.1 linearly up to zt
bb=1;
bt=0.1;
ft=interp1([zb zt],[bb bt],z,'linear','extrap');
ft(ft>1)=1; %b=1 below braking height



%reference velocity as in EMM for ft=1 bsaed on linear LEM wpeak profile
f=0;
wbar=5; %f and wbar can change for each thermal = same as wdeb
wr=2*(w+3/2*f*wbar)/(3/2*f+1); 

%actual wpeak that will result from multiplication with z dependant ft
wp=ft.*wr/2*(3/2*f+1) - 3/2*f*wbar;
%N.B. f and wbar don't affect wp in the linear part as cancel - only affect in braking zone


top=findheight(z,19); %find index for 19km 

%max w for LEM between times 19.75 & 21.2 UTC
maxw=max(MaxW(1).w(:,1:18),[],2);

figure;
plot(wp(1:top),z(1:top));
hold on;
plot(maxw(1:top),z(1:top),'r');

[maxwt maxwti]=max(MaxW(1).w(:,1:18),[],1);
maxh=GridDan(1).Z(maxwti); %height of the max w (thermal centre)

%now need to input wr at zb and zt into EMM (linear part)
wr0=2*(w0+3/2*f*wbar)/(3/2*f+1);
wrti=findheight(z,zt);
wrt=wr(wrti); %linear extrapolation of wr up to zt