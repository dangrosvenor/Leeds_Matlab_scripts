function [histo_mean,histo_std,std_norm]=histo_mean_calc_func_MODIS(hist_vals2D,freq)

%calculates the mean of a 2D histogram of data
%Need to supply the values that correspond to each position in the histo
%(histo_vals2D). MODvar is the name of the MODIS variable and dataset_var is
%the name of the dataset within that variable given in string form
%e.g. 'data' for the full lat,lon data, 'timeseries' for the timeseries,
%etc. inds are the indices of the data that is required within that dataset

freq=squeeze(freq);
 
size_data = size(freq);
rep_inds=[1 1 size_data(3:end)];

%create 4D arrays of replicas of the histo_vals data for each lat,lon point
%vals4d=repmat(hist_vals2D,rep_inds);
vals4d = hist_vals2D;


%the NaN values within the histogram for each location
%are kinda meaningless - think they should just mean they are zeros
%(i.e. no pixels with those values). So set to zero here. Otherwise the
%whole histogram will give a sum of NaN and won't be used.
F = freq;
F(isnan(F))=0;

%total number of pixels for each lat,lon 
totN = squeeze( sum(sum(F,1),2) );

histo_mean = squeeze(  sum(sum(F.*vals4d,1),2)  ) ./totN;


%%%now calculate the std dev
rep_inds2=[ones(size(size_data(3:end))) size_data(1) size_data(2)];
%mean4d=repmat(histo_mean,[1 1 size_data(1) size_data(2)]);
mean4d=repmat(histo_mean,rep_inds2);

%F4d2=permute(F,[3 4 1 2]);  %re-arrange so that the next step works
F4d2=permute(F,[3:length(size_data) 1 2]);  %re-arrange so that the next step works

vals4d2=permute(vals4d,[3:length(size_data) 1 2]); %N4d are the N values as a funciton of the Reff-Tau histo values

strc=':';
for ic=1:length(size_data)-2
    strc=[strc ',:'];
end
%then replicated for each lat/lon point
%"expand" the last two dimensions to make them into a vector array
Npix3d=eval(['F4d2(' strc ')']); %Number of pixels with particular Reff-Tau values replicated over lat-lon
vals3d=eval(['vals4d2(' strc ')']);
mean3d = eval(['mean4d(' strc ')']); 

sumdim=length(size_data)-1;

histo_std=sqrt ( sum(Npix3d.*(vals3d - mean3d).^2 , sumdim) ./ (totN-1) );
histo_std(totN==0)=NaN;

std_norm=histo_std./histo_mean;

