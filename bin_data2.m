function [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals,meanvals_ignore_zeros,std_dev,mean_error]=bin_data2(X,Y,bins,A)
%bin data in bins of X and calcualte the means, max and indices of the max
%[meanvals,maxvals,max_inds]=bin_data2(X,Y,bins,rerr)
%A is the absolute error for each point for a calculation of the error of the bin means

if bins(1)>=bins(end)
    disp('*** Bins must be increasing ***');
    return
end

        a=isnan(Y);
        b=find(a==0);  
        
        Y3=Y(b);
        X3=X(b);  %ignore all the NaN data in Y
        A3=A(b);
        
        a=isnan(X3);
        b=find(a==0);  
        
        Y2=Y3(b);
        X2=X3(b);  %ignore all the NaN data in X
        A2=A3(b);  %ignore all the NaN data in X
        
        b=find(Y2~=0);          
        Y0=Y2(b);
        X0=X2(b);  %ignore all the zero data for a special mean
        
for i=1:length(bins)-1
    ii=find(X2>=bins(i) & X2<bins(i+1));
    if length(ii)>0      
        meanvals(i)=mean(Y2(ii));        
        med_vals(i)=median(Y2(ii));
        [maxvals(i) max_ind_ii]=max(Y2(ii));
        max_inds(i)=ii(max_ind_ii);
        nvals(i)=length(ii);
        sum_vals(i)=sum(Y2(ii));
        std_dev(i)=std(Y2(ii));
        %calculate the combined error of the mean value given all the relative errors of the individual
        %values (rerr)
        mean_error(i) = meanvals(i) * sqrt(sum((A2(ii)).^2))/sum(Y2(ii));
        
    else
        meanvals(i)=NaN;
        med_vals(i)=NaN;
        maxvals(i)=NaN;
        max_inds(i)=NaN;
        nvals(i)=NaN;    
        sum_vals(i)=NaN;
    end
    
    
    ii=find(X0>=bins(i) & X0<bins(i+1));
    if length(ii)>0      
       meanvals_ignore_zeros(i)=mean(Y0(ii));
    else
       meanvals_ignore_zeros(i)=NaN;
    end

    
    
    
end

mid_points = 0.5 *( bins(1:end-1) + bins(2:end) );





