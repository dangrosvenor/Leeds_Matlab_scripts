

%% UKESM ensemble
i=1; %first trend period
i=2; %second trend period

fprintf(1,'UKESM ens mean\n')
t = trend_dat_box;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))

fprintf(1,'UKESM ens min\n')
t = trend_dat_box_ens;
minval=1e99;
for ii=1:size(t,3)
   val = t{1,i,ii}.coeffs(2)*diff(t{1,i,ii}.x([1 end]));
   if val<minval
      ival = ii;
      minval = val;
   end
end
fprintf(1,'%f\n',100*minval)
fprintf(1,'%f\n',100*t{1,i,ival}.uncer_max*diff(t{1,i,ival}.x([1 end])))

fprintf(1,'UKESM ens max\n')
t = trend_dat_box_ens;
maxval=-1e99;
for ii=1:size(t,3)
   val = t{1,i,ii}.coeffs(2)*diff(t{1,i,ii}.x([1 end]));
   if val>maxval
      ival = ii;
      maxval = val;
   end
end
fprintf(1,'%f\n',100*maxval)
fprintf(1,'%f\n',100*t{1,i,ival}.uncer_max*diff(t{1,i,ival}.x([1 end])))


%%
fprintf(1,'Hist-GHG\n')
t = trend_dat_box_hist_GHG;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))

fprintf(1,'Hist-nat\n')
t = trend_dat_box_hist_nat;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))


fprintf(1,'Sum\n')
t = trend_dat_box_hist_linear;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))


%% HADEGM and DAMIP
i=1; %first trend period
i=2; %second trend period

fprintf(1,'HADGEM\n')
t = trend_dat_box;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))

fprintf(1,'Hist-aer\n')
t = trend_dat_box_hist_aer;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))

fprintf(1,'Hist-GHG\n')
t = trend_dat_box_hist_GHG;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))

fprintf(1,'Hist-nat\n')
t = trend_dat_box_hist_nat;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))


fprintf(1,'Sum\n')
t = trend_dat_box_hist_linear;
fprintf(1,'%f\n',100*t{i}.coeffs(2)*diff(t{i}.x([1 end])))
fprintf(1,'%f\n',100*t{i}.uncer_max*diff(t{i}.x([1 end])))