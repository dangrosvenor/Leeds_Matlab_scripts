function [histo_output]=histo_mean_calc_MODIS_run(histo_var,histo_type_str,histo_freqs,CTT)
%function histo_mean_calc_func_MODIS(histo_var,histo_type_str,dataset_str,inds)
%calculate the N,H values for the MODIS-L3 joint histogram bins
%For tau-reff histos this gives an N, H and W calculated from each combination of tau and
%reff given by the histogram intervals (using mid-points of the bin edges)
%histo_var is the name of the structure containing the MODIS histogram variable
%histo_type_str is a string saying what type of histogram it is (each one will
%allow different quantities to be calculated from it)
%histo_type_str='tau-reff' for tau reff histograms that allow N, H and W to
%be calculated
%dataset_str is the string for the name of the dataset within the structure
%that contains the desired frequencies (e.g. could be 'data' or
%'timeseries'
%inds_str are the indices within that dataset in string format (e.g. '(:,1:10)' - for the whole array set to []


% histo_freqs = eval(['histo_var.' dataset_str]); %all the data
% if length(inds_str)~=0
%     histo_freqs = eval(['histo_freqs' inds_str]);
% end


%NaN values are removed from histo_freqs in histo_mean_calc_func_MODIS


switch histo_type_str
    case 'tau-reff'

        tau_bins=histo_var.Ybins;
        reff_bins=histo_var.Xbins*1e-6; %convert to metres

        tau = mid_vals(tau_bins);
        reff = mid_vals(reff_bins); 

        %experiment with log variation within bins - different mean values for the
        %bins. Doesn't change the bin means by much at all
        %tau = (tau_bins(2:end)-tau_bins(1:end-1)) ./ log(tau_bins(2:end)./tau_bins(1:end-1)) ;
        %reff = (reff_bins(2:end)-reff_bins(1:end-1)) ./ log(reff_bins(2:end)./reff_bins(1:end-1)) ;

         %WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
        [tau2d,reff2d] = meshgrid(tau,reff);  %make a 2D grid of all the tau,reff comibations

%        itau_lim=find(tau2d>50);
%        tau_lim_zeros = ones(size(tau2d));
%        tau_lim_zeros(itau_lim)=0;


         size_data = size(squeeze(histo_freqs));
         rep_inds=[1 1 size_data(3:end)];
         rep_inds2=[1 1 size_data(1:2)];

        tau4d=repmat(tau2d,rep_inds);
        reff4d=repmat(reff2d,rep_inds);
        
        CTT_4d=repmat(CTT,rep_inds2);
        LL=length(size(CTT_4d));
        CTT_4d = permute(CTT_4d,[LL-1 LL 1 2]);
        
        [N4d,H4d,W4d,k,Q,cw4d]=MODIS_N_H_func(tau4d,reff4d,'calc',NaN,CTT_4d);
        
        LWC_CT_4d = cw4d .* H4d; %cloud top LWC
        


        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)

        %calculate Nd, H and W for the each square of the tau-reff histogram
%        [N2d,H2d,W2d,k,Q,cw]=MODIS_N_H_func(tau2d,reff2d,Wflag,0);
  
%note - had W and H the wrong way around when passing to these functions as
%of 6th Nov, 2011 - meant that W_timeseries was actually H_timeseries and
%vice versa
        [histo_output.N_histo_mean,histo_output.N_histo_std,histo_output.N_std_norm]=histo_mean_calc_func_MODIS(N4d,histo_freqs);
        [histo_output.W_histo_mean,histo_output.W_histo_std,histo_output.W_std_norm]=histo_mean_calc_func_MODIS(W4d,histo_freqs);
        [histo_output.H_histo_mean,histo_output.H_histo_std,histo_output.H_std_norm]=histo_mean_calc_func_MODIS(H4d,histo_freqs);
        [histo_output.LWC_histo_mean,histo_output.LWC_histo_std,histo_output.LWC_std_norm]=histo_mean_calc_func_MODIS(LWC_CT_4d,histo_freqs);        
       
end