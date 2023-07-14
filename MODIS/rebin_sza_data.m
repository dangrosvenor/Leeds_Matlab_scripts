    %%%%% rebin the data into larger sza bins.... %%%%%
    
clear Nd_rebin Np_rebin
Nwindow=2;

[a,b] = window_average(sza1D,Nd_sza(:,1,1),Nwindow,'mean');

Nd_rebin = NaN*ones([length(a) size(Np_sza,2) size(Np_sza,3)]);
Np_rebin = NaN*ones([length(a) size(Np_sza,2) size(Np_sza,3)]);



    for ilat=1:size(Nd_sza,2)
        ilat
        for iday_rebin=1:size(Nd_sza,3)
            [new_sza_bins,Nd_rebin(:,ilat,iday_rebin)] =  window_average(sza1D,Nd_sza(:,ilat,iday_rebin),Nwindow,'mean');
            [new_sza_bins,Np_rebin(:,ilat,iday_rebin)] =  window_average(sza1D,Np_sza(:,ilat,iday_rebin),Nwindow,'sum');
        end
    end