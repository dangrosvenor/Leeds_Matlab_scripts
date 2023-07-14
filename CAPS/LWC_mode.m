function [LWC_mode_mid,LWC_mode_lower,LWC_mode_upper,LWC_percent_contribution]=LWC_mode(LWC_size_dist,CAS_bins,LWC_min)

if nargin<3
    LWC_min=0.05;
end
    
for i=1:size(LWC_size_dist,1)
    LWC=sum(LWC_size_dist(i,:));    
    [a,b]=max(LWC_size_dist(i,:));
    
%    b=24;    
%    a=LWC_size_dist(i,b);
    
    if LWC>LWC_min
        LWC_mode_lower(i)=CAS_bins(b);
        LWC_mode_upper(i)=CAS_bins(b+1);
        LWC_mode_mid(i)=0.5*(LWC_mode_lower(i)+LWC_mode_upper(i));
        LWC_percent_contribution(i)=100*a/LWC;
    else
        LWC_mode_lower(i)=0;
        LWC_mode_upper(i)=0;
        LWC_mode_mid(i)=0;
        LWC_percent_contribution(i)=0;
    end
    
end