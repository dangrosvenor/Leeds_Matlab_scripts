function [time_100,var_100,anom_100,N_100] = RyanE_traj_analysis_100_day_anomalies(time,var,nthresh_days)
% Calculates the 100-day running mean of a variable.
% Require a minimum of nthesh non-NaN points within the 100-day running mean before being accepted.

Ndays=100;

c=permute(var,[3 1 2]); %put the time dimension first
%calculate the number of time points required
dtime = diff(time(1:2)); %in days
Nwindow = round(Ndays/dtime);
%Make sure it's an odd number since makes things easier with choosing the
%part of the array we want.
if mod(Nwindow,2)==0
    Nwindow=Nwindow+1;
end
nthresh = nthresh_days/Ndays * Nwindow;
[time_100,var_100,N_100]=window_average_NaN(time,c,Nwindow,'mean',nthresh);
istart=(Nwindow+1)/2;
nT=size(var_100,1);
var_100 = permute(var_100,[2 3 1]);
N_100 = permute(N_100,[2 3 1]);
dat_anom = var(:,:,istart:istart+nT-1);
anom_100 = dat_anom - var_100;
