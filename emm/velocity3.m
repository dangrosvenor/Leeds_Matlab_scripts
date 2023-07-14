%z=GridDan(1).Z/1000 + 0.62;

z=-10:0.01:20;



dt=5;
time=[0:dt:3600]; %timesteps as in EMM dynamics

%remember these are thermal velocities wth - wr=wth*2 is entered into the EMM
z0=3.5;  %cloud base for w
w0=8.5; %m/s for cloud base
zb=10.2; %start of braking height (of thermal centre so add tdepth/2 to this for braking height- set as height of cldtop_max above braking height 
zt=18.5; %cloud top
wb=10.5; %m/s at start of braking height 

f=0.8; %set in thermal_amplification subroutine in EMM
wbar=2; %f and wbar can change for each thermal = same as wdeb

tdepth=12; %depth of thermal (km) - set in generalinit (so is ceiling)

thinit=4+(tdepth/2); %start of thermal top

%linear change of wref from LEM
w=interp1([z0 zb],[w0 wb],z,'linear','extrap');

%ft=1 up to braking height, zb, and then reduces to 0.1 linearly up to zt
bb=1;
bt=0.00;
ft=interp1([zb zt],[bb bt],z,'linear','extrap');
ft(ft>1)=1; %b=1 below braking height



%reference velocity as in EMM for ft=1 bsaed on linear LEM wpeak profile

%wr=2*(w+3/2*f*wbar)/(3/2*f+1); 

wr=2*w; %w is the linear part of the profile

wr(wr<w0*2)=w0*2; %min of wr is w0*2

wth=ft.*wr/2;

%actual wpeak that will result from multiplication with z dependant ft
%wp=ft.*wr/2*(3/2*f+1) - 3/2*f*wbar;


%N.B. f and wbar don't affect wp in the linear part as cancel - only affect in braking zone

wp=wth + (3/2)*f*(wth-wbar);


clear top bot zmid wmid wmax

top(1)=thinit;
bot(1)=max([thinit-tdepth 0]);

for i=2:length(time)
    itop=findheight(z,top(i-1));
    vtop(i)=wr(itop)/2;
    
    if bot(i-1)<z0;
        bot(i-1)=top(i-1)-tdepth;
        ibot=findheight(z,bot(i-1)); %if thermal is below cloud base make = top-depth
    else
        ibot=findheight(z,bot(i-1));
    end
    vbot(i)=wr(ibot)/2;
 
    
    
    zbt=zb+tdepth/2;
    grad=(1-bt)/(zbt-zt);
    if top(i-1)>zbt
        vtop(i)=vtop(i)*(1-grad*(zbt-top(i-1)));
    end
    if bot(i-1)>zbt
        vbot(i)=vbot(i)*(1-grad*(zbt-bot(i-1)));
    end
        
    top(i)=top(i-1)+vtop(i)*dt/1000;
    bot(i)=bot(i-1)+vbot(i)*dt/1000;
    
    %bottemp=max([z0 bot(i)]);
    %i0=findheight(z,bottemp);
    %zmid(i)=z(round((itop+i0)/2));
    zmid(i)=(bot(i)+top(i))/2;
    wmid(i)=(vtop(i)+vbot(i))/2;
   
  vmidh=1;  
  if vmidh==1
    imid=findheight(z,zmid(i-1));
    wmid(i)=wr(imid)/2;
    grad=(1-bt)/(zb-zt);
    if zmid(i-1)>zb
        wmid(i)=wmid(i)*(1-grad*(zb-zmid(i-1)));
    end
  end
    
    %wcrit a function of height of thermal top above cloud base
    %fcrit=1.5+(top(i)-z0)*0.6;
    
    fcrit=1.5+(top(i)-zb)*0.9;
    wcrit=wmid(i) + (3/2)*fcrit;
    wmax(i)=wmid(i) + (3/2)*f*(wmid(i)-wbar);
    
    
    
    %wmax(i)=max([wmax(i) wcrit]); %replace with the wcrit term if is bigger
    
    if top(i)>=zt; break; end;
    
end


%top=findheight(z,19); %find index for 19km 

%max w for LEM between times 19.75 & 21.2 UTC
maxw=max(MaxW(1).w(:,1:18),[],2);

% figure;
% 
% plot(time(1:length(wmax))/60,wmax,'r');
% hold on;
% plot(Wmax(:,3),Wmax(:,1));
% 
% figure;
% 
% plot(time(1:length(bot))/60,top,'r');
% hold on;
% plot(thermalheight(:,1),thermalheight(:,2));
% plot(time(1:length(bot))/60,bot,'r');
% plot(thermalheight(:,1),thermalheight(:,3));

figure;
set(gcf,'name','wmax profile');
plot(wmax,zmid);
hold on;
plot(Wmax(:,1),Wmax(:,2),'r');

tim=time(1:length(bot));


% figure;
% plot(wp(1:top),z(1:top));
% hold on;
% plot(maxw(1:top),z(1:top),'r');
% 
 [maxwt maxwti]=max(MaxW(1).w(:,1:18),[],1);
 maxh=GridDan(1).Z(maxwti); %height of the max w (thermal centre)
 
%max w for LEM between times 19.75 & 21.2 UTC
maxw=max(MaxW(1).w(:,1:18),[],2);
top=findheight(GridDan(1).Z,19e3); %find index for 19km 
plot(maxw(1:top),GridDan(1).Z(1:top)/1000,'k--');


% 
% %now need to input wr at zb and zt into EMM (linear part)
% wr0=2*(w0+3/2*f*wbar)/(3/2*f+1);
% wrti=findheight(z,zt);
% wrt=wr(wrti); %linear extrapolation of wr up to zt
% 
% %wcrit a function of height of thermal top above cloud base
% zbraketop=zb+tdepth/2;
% wcrit=1.5+(z+tdepth/2-z0)*0.6;
% wpcrit=wth + (3/2)*wcrit;
% plot(wpcrit(1:top),z(1:top),'k--');

%plot(maxwt,maxh);

