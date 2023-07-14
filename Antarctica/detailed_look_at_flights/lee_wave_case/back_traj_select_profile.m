x_prof=1.46e5; %location of profile
ix_prof=findheight_nearest(X2,x_prof);

%z2_sound=z;
%T_h_sound=T_h; %temperature (K)
%qv2_sound=qv_prof_humi;
%P2_sound=P_wave;

T_h2 = temp3(ix_prof,:);
inan=isnan(T_h2);
inan=find(inan==0);


lapse_rate=6; %K/km - doesn't really matter?
inds_top=[inan(end)+1:length(T_h2)];
T_h2(inds_top)=T_h2(inan(end)) -  lapse_rate/1000*(Zsound_pressure(inds_top)-Zsound_pressure(inan(end)));

inds=[inan(1):length(T_h2)];
T_h=T_h2(inds);
z = Zsound_pressure(inds);
qv_prof_humi = tot3(ix_prof,inds) - L3(ix_prof,inds);
P_wave = interp1(Zsound_pressure,Psound_extended,z);

Zsound=Zsound_pressure(inan);  %just the valid points - the rest will be made
        %through assumptions such as constant (low) RH, etc.
        
iuse_mean_extrap=0;

vert_vel = W2(ix_prof:end);
time_wave = ( X2(ix_prof:end)-X2(ix_prof) ) / Uconstant;

