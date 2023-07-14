%calculate the N,H values for the MODIS-L3 joint histogram bins
%so this gives an N, H and W calculated from each combination of tau and
%reff given by the histogram intervals (using mid-points of the bin edges)

tau_bins=Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Ybins;
reff_bins=Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Xbins*1e-6;

tau = mid_vals(tau_bins);
reff = mid_vals(reff_bins); %convert to metres

%experiment with log variation within bins - different mean values for the
%bins. Doesn't change the bin means by much at all
%tau = (tau_bins(2:end)-tau_bins(1:end-1)) ./ log(tau_bins(2:end)./tau_bins(1:end-1)) ;
%reff = (reff_bins(2:end)-reff_bins(1:end-1)) ./ log(reff_bins(2:end)./reff_bins(1:end-1)) ;



%WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
[tau2d,reff2d] = meshgrid(tau,reff);  %make a 2D grid of all the tau,reff comibations

itau_lim=find(tau2d>50);
tau_lim_zeros = ones(size(tau2d));
tau_lim_zeros(itau_lim)=0;




if ~exist('set_MODIS_NH_flags') | set_MODIS_NH_flags==0

    Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
    %Wflag='MODIS'; %use the MODIS LWP

else
    clear set_MODIS_NH_flags
end

%WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
WMOD=NaN;

%calculate N, H and W values for each combination of tau and reff in the
%histogram
[N2d,H2d,W2d,k,Q,cw]=MODIS_N_H_func(tau2d,reff2d,Wflag,WMOD);

size_data = size(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.data);

%create 4D arrays of replicas of N2d,H2d for each lat,lon point
N4d=repmat(N2d,[1 1 size_data(3),size_data(4)]);
H4d=repmat(H2d,[1 1 size_data(3),size_data(4)]);
W4d=repmat(W2d,[1 1 size_data(3),size_data(4)]);
Fzero4d=repmat(tau_lim_zeros,[1 1 size_data(3),size_data(4)]);

F = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.data;
F(isnan(F))=0;
%F=F.*Fzero4d; %remove frequency values when tau is too high
%fprintf(1,'\n**** WARNING - am limiting tau histogram frequencies to <=50 ***\n');

%total number of pixels for each lat,lon 
totN = squeeze( sum(sum(F,1),2) );

N_histo_mean = squeeze(  sum(sum(F.*N4d,1),2)  ) ./totN;
H_histo_mean = squeeze(  sum(sum(F.*H4d,1),2)  ) ./totN;
W_histo_mean = squeeze(  sum(sum(F.*W4d,1),2)  ) ./totN;

%now find the mode value of Nd, etc. - i.e. the value at the histgram
%location where the most pixels were contained
%first reshape to roll out all of the 2D histos into a linear 110 vector
%(from a 10 by 11 array)

F2=F(:); %roll out completely
F2=reshape(F2,[110 180*360]); %reshape so just the histogram parts and the lat lon parts
%are reshaped separately
%I guess this only works if the elements are the first two dimensions of
%the array

%find the maximum frequency (i.e. the mode) of the histo
[Fmax,iFmax]=max(F2,[],1);
Fmode = reshape(Fmax,[180 360]); %save the Fmax values as might be useful
%for testing robustness - e.g. if the mode value has only a low percentage
%of the number of points, perhaps it is less reliable - can screen later
%for this

%now find N, H and W at the mode position from the 2d histo arrays
Nmode = N2d(iFmax);
Nmode = reshape(Nmode,[180 360]);

Hmode = H2d(iFmax);
Hmode = reshape(Hmode,[180 360]);

Wmode = W2d(iFmax);
Wmode = reshape(Wmode,[180 360]);








