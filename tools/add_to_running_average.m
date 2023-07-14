function [tot_out,Nvals_out] = add_to_running_average(dat,tot,Nvals)
%adds dat to a running total (tot), but only where have non-NaN data
%also returns the running total number of datapoints in the total.

%tot is a running total
tot_out = tot;
Nvals_out = Nvals;

inan=find(isnan(dat)==1);

%make all the NaN values 0 and the others 1
Nvals = ones(size(Nvals));
Nvals(inan) = 0;

%increment the number of datapoints in the non-NaN total
Nvals_out = Nvals_out + Nvals;

%increment the data total, but only for non-NaN values
dat(inan) = 0;
tot_out = tot_out + dat;

