function [coeffs_out,t_out,bint_out,stats_out,p_Durbin_out,d_Durbin_out,t_out2,uncer_out2,uncer_out,Neff_out] = trend_linear_fit_nd_data(x1,dat,dim)


ndim=ndims(dat);
permute_inds=[1:ndim];
%remove the requested dimension from the list
permute_inds(dim)=[];
permute_inds = [dim permute_inds];
dat = permute(dat,permute_inds);
siz=size(dat);
%Putting the required dim at the start allows us to compact the array into
%a 2d one, with all the other dims together
dat = dat(:,:);

coeffs_out = NaN*ones([2 size(dat,2)]);
t_out = NaN*ones([1 size(dat,2)]);

for i=1:size(dat,2)
    y = dat(:,i)';   
    x = [ones(size(y')) x1];
    %Catch cases with NaNs or no data and make the output NaN
    if length(find(isnan(y)==1)) > 0 | length(x)==0 | length(y)==0 
        coeffs=[NaN; NaN]; bint=[NaN NaN; NaN NaN]; residuals=NaN*ones([size(x,1) 1]); rint=NaN*ones(size(x)); t=NaN;
        stats=[NaN NaN NaN NaN];
        t2 = NaN; p_Durbin=NaN; d_Durbin=NaN; uncer=NaN; uncer2=NaN; Neff=NaN;
    else        
        [coeffs,bint,residuals,rint,stats] = regress(y',x);
        %Got the following from Chapter 1 of some lecture notes I found -
        %Chapter 1 - Data Analysis :- 
        %<file://C:\Users\eardgro\Documents\logbook_OLD\lecture notes etc\Statistical_significance_of_linear_trend.pdf>
        
        %Is based on Santer et al. - Statistical significance of trends and trend
        %differences in layer-average atmospheric temperature time series,
        %JGR, 2000        
        [t,serr,t2,serr2,Neff] = t_value(coeffs(2),residuals,x(:,2),y);
        %x here has a column of ones in first column, so just pass the
        %second column
        
        %serr is the standard error of the fit. Uncertainty is usually taken
        %as 2* this ; approx matches the bint values from Matlab (Matlab
        %came out as 2.05 * serr in one test)
        %serr2 is the standard error when the effective number of samples
        %is reduced according to Santer, JGR, 2000 formula
        %Santer uses 1.96*serr
        uncer = 2*serr;
        uncer2 = 2*serr2;
              
        %Durbin-Watson test for degree of auto-correlation
        [p_Durbin,d_Durbin] = dwtest(residuals,x);
            %For no auto-correlation we want a high p-value (range 0 to 1)
            %and a d-value near 2.
            %p is the probability that there is no auto-correlation
            
        
    end
    
    coeffs_out(:,i) = coeffs;
    bint_out(:,:,i) = bint; %The 95% confidence intervals of the coeff values (y-int and slope)
    stats_out(:,i) = stats; %Stats contains R^2, the F-statistic and it's p-value,
        %and an estimate of its error variance
    t_out(i) = t; 
    t_out2(i) = t2; 
    p_Durbin_out(i) = p_Durbin;
    d_Durbin_out(i) = d_Durbin;
    uncer_out(i) = uncer;
    uncer_out2(i) = uncer2;
    Neff_out(i) = Neff;
    
end

coeffs_out = reshape(coeffs_out,[2 siz(2:end)]);
bint_out = reshape(bint_out,[2 2 siz(2:end)]);
stats_out = reshape(stats_out,[4 siz(2:end)]);
t_out = squeeze(reshape(t_out,[1 siz(2:end)]));
t_out2 = squeeze(reshape(t_out2,[1 siz(2:end)]));
p_Durbin_out = squeeze(reshape(p_Durbin_out,[1 siz(2:end)]));
d_Durbin_out = squeeze(reshape(d_Durbin_out,[1 siz(2:end)]));
uncer_out = squeeze(reshape(uncer_out,[1 siz(2:end)]));
uncer_out2 = squeeze(reshape(uncer_out2,[1 siz(2:end)]));
Neff_out = squeeze(reshape(Neff_out,[1 siz(2:end)]));
