%%%% ****  get the environmental (ascent & descent) sounding for flight 102 ****
dat_flt=dat_flt102;
%get the profile performed at the end of the flight from the flight data
times=[21.6772 22.0281]/24; %flight 102 end profile - perhaps the better one to use because
   %the ascent profile occurred quite close to the Peninsula and so may have been affected by 
   %terrain effects quite strongly.
times2=[19.6861 20.3949]/24; %flight 102 start profile   
[time_0,time_1]=findheight_nearest(dat_flt(:,1)/1000/3600/24,times(1,1),times(1,2));
[time_0_2,time_1_2]=findheight_nearest(dat_flt(:,1)/1000/3600/24,times2(1,1),times2(1,2));
inds=time_0:time_1;
inds2=time_0_2:time_1_2;

[Zsound,isort]=sort(dat_flt(inds,12)); %sort in height order  - end descent
Tsound=dat_flt(inds,5); Tsound=Tsound(isort);
Psound=dat_flt(inds,6); Psound=Psound(isort);
qfp_sound=qv_flt102_fp(inds); qfp_sound=qfp_sound(isort);
qhumi_sound=qv_flt102_humi(inds); qhumi_sound=qhumi_sound(isort);
Usound=dat_flt(inds,col_wind); Usound=Usound(isort);

[Zsound2,isort2]=sort(dat_flt(inds2,12)); %sort in height order - start, ascent
Tsound2=dat_flt(inds2,5); Tsound2=Tsound2(isort2);
Psound2=dat_flt(inds2,6); Psound2=Psound2(isort2);
qfp_sound2=qv_flt102_fp(inds2); qfp_sound2=qfp_sound2(isort2);
qhumi_sound2=qv_flt102_humi(inds2); qhumi_sound2=qhumi_sound2(isort2);

equiv_sound_humi = equivalent_potemp(273.15+Tsound,100*Psound,qhumi_sound);
equiv_sound_fp = equivalent_potemp(273.15+Tsound,100*Psound,qfp_sound);

equiv_sound_humi2 = equivalent_potemp(273.15+Tsound2,100*Psound2,qhumi_sound2);
equiv_sound_fp2 = equivalent_potemp(273.15+Tsound2,100*Psound2,qfp_sound2);

Wsound=w2_turb(inds); Wsound=Wsound(isort);
Pot_sound = (273.15+Tsound) .* (1000./Psound).^0.286;
