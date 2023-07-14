try

%tests using Latin hypercube sampling on tau and re to assess the effect on
%Nd. Can add correlation between the tau and re - doesn't make a lot of
%difference.
%Results are fairly similar to that based on just using the mean values
%Sensitive to the cut-off value that needs to be imposed for re to ensure
%non-negative values (should prob do for tau too).
%The sampler used assumes a normal distribution - can introduce other
%distributions into it. Would be best to introduce the actual distribution
%(with zero prob for values below a certain tau and re value).
%This would circumvent the cut-off problem.
%See make_sza_rel_increases_table_vals2 for script that runs this several times and 
%makes the table


if ~exist('override_attr') | override_attr==0
    use_pdfs = 1; %flag to say that we want to use specific PDFs rather than normal gaussians
    band_str='16'; %string to define which Re band to use
    band_str='21';
    band_str='37';
end

facNd = 1.2034e+04;

%order is tau then re
%vals={[20 11],[30 11],[20 10],[30 10]};
base = [30 15];
xsd=[3 3];

ftau=1; fre=1;
ftau_std=1; fre_std=1; %factors by which to change the std dev at low and high sza
%i.e. by changing just fre_std can assess the effect of keeping the mean
%constant, but having a change in the width of the re distribution - has a
%large effect

if use_pdfs==0
    vals={base,[base(1)*ftau  base(2)],[base(1)  base(2)*fre],[base(1)*ftau  base(2)*fre]};
    vals_std={xsd,[xsd(1)*ftau_std  xsd(2)],[xsd(1)  xsd(2)*fre_std],[xsd(1)*ftau_std  xsd(2)*fre_std]};
else
%     [me_Tau01,std_Tau01] = pdf_me_std(xTau(1).x,yTau(1).y);
%     [me_Tau02,std_Tau02] = pdf_me_std(xTau(2).x,yTau(2).y);
%     [me_Re01,std_Re01] = pdf_me_std(xRe37(1).x,yRe37(1).y);
%     [me_Re02,std_Re02] = pdf_me_std(xRe37(2).x,yRe37(2).y);
        
%    vals={[me_Tau01 me_Re01],[me_Tau02 me_Re01],[me_Tau01 me_Re02],[me_Tau02 me_Re02]};
%    vals_std={[std_Tau01 std_Re01],[std_Tau02 std_Re01],[std_Tau01 std_Re02],[std_Tau02 std_Re02]};    
    
    %need to also create arrays of the PDFs. Need to be in the for
    %[nsample,nvar] where nsample is the number of datapoints and nvar the
    %number of variables (e.g. Tau and Re).
clear xvals_orig fvals_orig
        xvals_orig{1} = xTau(1).x;
        eval(['xvals_orig{2} = xRe' band_str '(1).x;']);
        fvals_orig{1} = yTau(1).y;
        eval(['fvals_orig{2} = yRe' band_str '(1).y;']);
        
        xvals_tau{1} = xTau(2).x;
        eval(['xvals_tau{2} = xRe' band_str '(1).x;']);
        fvals_tau{1} = yTau(2).y;
        eval(['fvals_tau{2} = yRe' band_str '(1).y;']);        
        
        xvals_re{1} = xTau(1).x;
        eval(['xvals_re{2} = xRe' band_str '(2).x;']);        
        fvals_re{1} = yTau(1).y;
        eval(['fvals_re{2} = yRe' band_str '(2).y;']);
        
        xvals_both{1} = xTau(2).x;
        eval(['xvals_both{2} = xRe' band_str '(2).x;']);        
        fvals_both{1} = yTau(2).y;
        eval(['fvals_both{2} = yRe' band_str '(2).y;']);   
        
        corr_actual_01 = eval(['corr(Tau_vals(1).dat,Re' band_str '_vals(1).dat);']);
        corr_actual_02 = eval(['corr(Tau_vals(2).dat,Re' band_str '_vals(2).dat);']);        

end
%note xvals are the bins edges - required to get the cumulative PDF right

re_cutoff = 1;
nsample = 1e5;

%correlation determined by e.g. corr_choose(1,2) is correlation between
%variable 1 and 2 and vice versa for corr_choose(2,1)
nocorr = [1 0;0 1]; %no correlation
corr2 = [1 -0.99; -0.99 1]; %negative correlation
%corr2 = [1 0.99; 0.99 1];
corr2 = [1 -0.22; -0.22 1];
%corr2 = [1 0.37; 0.37 1];

corr_actual_01_mat = [1 corr_actual_01; corr_actual_01 1];
corr_actual_02_mat = [1 corr_actual_02; corr_actual_02 1];

corr_choose=nocorr;
%corr_choose=corr2;
use_corr_actual = 0; %but only for the base values - i.e. orig and both calculations
%for the delta_tau and delta_re operations would need to choose a
%correlation to use between the low and high re and tau dataset. Since the
%correlation changes with sza this is difficult.

%calling N1 the value at low SZA (original tau and re values)
%N2 is the value with both tau and re changed (i.e. the actual Nd that
%would be obtained assuming that the LHS worked as in reality)
%N_orig_mean is the original N value (low sza from teh mean values)
%N_orig_lhs is the low sza value obtained using the PDFs (LHS)
%Ntau_mean is the N value from just changing the tau value and using only
%the means in the calculation. Similarly for Nre_mean
%Ntau_lhs and Nre_lhs are the corresponding values using the LHS method
% si is the scaling factor = N_orig_lhs/N_orig_mean



%-----   using Latin Hypercube sampling -----------------------------------------


%original values

if use_pdfs==1
     if use_corr_actual==1
         z_orig=lhs_iman_owndist_Dan(corr_actual_01_mat,nsample,xvals_orig,fvals_orig); 
     else
         z_orig=lhs_iman_owndist_Dan(corr_choose,nsample,xvals_orig,fvals_orig); %
    %function z=lhs_iman_owndist(xmean,xsd,corr_choose,nsample,x_pdf,f_pdf,ntry)
     end
else
    xmean = vals{1}; %orig values
    xsd2 = vals_std{1};   
    z_orig=lhs_iman(xmean,xsd2,corr_choose,nsample); %
    icut = find(z_orig(:,2)<re_cutoff);
    z_orig(icut,:)=[];
end



if use_pdfs==1
    z_tau=lhs_iman_owndist_Dan(corr_choose,nsample,xvals_tau,fvals_tau); %
else
    xmean = vals{2}; %just tau increased
    xsd2 = vals_std{2};
    z_tau=lhs_iman_owndist_Dan(xmean,xsd2,corr_choose,nsample); %
    icut = find(z_tau(:,2)<re_cutoff);
    z_tau(icut,:)=[];
end


if use_pdfs==1
    z_re=lhs_iman_owndist_Dan(corr_choose,nsample,xvals_re,fvals_re); %
else
    xmean = vals{3}; %just re decreased
    xsd2 = vals_std{3};
    z_re=lhs_iman_owndist_Dan(xmean,xsd2,corr_choose,nsample); %
    icut = find(z_re(:,2)<re_cutoff);
    z_re(icut,:)=[];
end


if use_pdfs==1
    if use_corr_actual==1
        z_both=lhs_iman_owndist_Dan(corr_actual_02_mat,nsample,xvals_both,fvals_both); %
    else
        z_both=lhs_iman_owndist_Dan(corr_choose,nsample,xvals_both,fvals_both); %
    end
else
    xmean = vals{4};  %both tau and re changed
    xsd2 = vals_std{4};          
    icut = find(z_both(:,2)<re_cutoff);
    z_both(icut,:)=[];
end

% ----- end of the sampling ------

%     [me_Re02,std_Re02] = pdf_me_std(xRe37(2).x,yRe37(2).y);

me_Tau01=mean(z_orig(:,1)); std_Tau01=std(z_orig(:,1));
me_Tau02=mean(z_both(:,1)); std_Tau02=std(z_both(:,1));
me_Re01=mean(z_orig(:,2)); std_Re01=std(z_orig(:,2));
me_Re02=mean(z_both(:,2)); std_Re02=std(z_both(:,2));

me_Nactual_pdf_01 = eval(['sum(0.5*(xNd' band_str '(1).x(1:end-1)+xNd' band_str '(1).x(2:end)).*yNd' band_str '(1).y);']);
me_Nactual_pdf_02 = eval(['sum(0.5*(xNd' band_str '(2).x(1:end-1)+xNd' band_str '(2).x(2:end)).*yNd' band_str '(2).y);']);

%the actual means from all teh X values (loaded in)
me_Nactual_01 = eval(['meanN' band_str '(1);']);
me_Nactual_02 = eval(['meanN' band_str '(2);']);

vals={[me_Tau01 me_Re01],[me_Tau02 me_Re01],[me_Tau01 me_Re02],[me_Tau02 me_Re02]};
vals_std={[std_Tau01 std_Re01],[std_Tau02 std_Re01],[std_Tau01 std_Re02],[std_Tau02 std_Re02]};
    

% ----- Results -----------------------------
%-----   using the means -----------------------------------------

N_orig_mean = (facNd*vals{1}(1)^0.5 .* vals{1}(2)^(-2.5));
%varying both tau and re
r_both=( facNd*vals{4}(1)^0.5 .* vals{4}(2)^(-2.5) ) / N_orig_mean;
% r_both =
%     1.5543
%vary just tau
r_tau=( facNd*vals{2}(1)^0.5 .* vals{2}(2)^(-2.5) ) / N_orig_mean;
% r_tau =
%     1.2247
%vary juse re
r_re=( facNd*vals{3}(1)^0.5 .* vals{3}(2)^(-2.5) ) / N_orig_mean;
% r_re =
% 
%     1.2691
N_tau_mean = r_tau * N_orig_mean;
N_re_mean = r_re * N_orig_mean;
N_both_mean = r_both * N_orig_mean;

fprintf(1,'\nr_both=%f r_tau=%f r_re=%f (means)',r_both,r_tau,r_re)


% --------------  Results using the LHS ---------------------------

%z_both is an [N,i] sized array of values sampled from the input PDFS - i=1 is tau and i=2 is re


%N_orig2 = facNd*mean(z_orig(:,1)).^0.5 .* mean(z_orig(:,2)).^(-2.5); %this
%is just akin to N_orig since the mean vlaues of the PDFs should be those
%specified (or close to)
N_orig_lhs = facNd*mean(z_orig(:,1).^0.5 .* z_orig(:,2).^(-2.5));

r_both2 = ( facNd*mean(z_both(:,1).^0.5 .* z_both(:,2).^(-2.5)) ) / N_orig_lhs;
r_tau2 = ( facNd*mean(z_tau(:,1).^0.5 .* z_tau(:,2).^(-2.5)) ) / N_orig_lhs;
r_re2 = ( facNd*mean(z_re(:,1).^0.5 .* z_re(:,2).^(-2.5)) ) / N_orig_lhs;

N_both_lhs = r_both2 * N_orig_lhs;
N_tau_lhs = r_tau2 * N_orig_lhs;
N_re_lhs = r_re2 * N_orig_lhs;

fprintf(1,'\nr2_both=%f r2_tau=%f r2_re=%f (LHS)\n',r_both2,r_tau2,r_re2);


%looking at the f factor as used in the paper draft

%for mean values. orig means low sza, both means high sza.
siA = N_orig_lhs/N_orig_mean; %can calc the s factor using either the low or high sza values...
siB = N_both_lhs/N_both_mean;

si=siA;
si=siB;
%si = mean([siA siB]);

f_tau_mean =  (N_tau_mean*si - N_orig_lhs) / (N_both_lhs - N_orig_lhs);
f_re_mean =  (N_re_mean*si - N_orig_lhs) / (N_both_lhs - N_orig_lhs);
f_both_mean =  (N_both_mean*si - N_orig_lhs) / (N_both_lhs - N_orig_lhs);

fprintf(1,'\nf_both_mean=%f f_tau_mean=%f f_re_mean=%f f_rem_mean=%f\n',f_both_mean,f_tau_mean,f_re_mean,1-f_tau_mean-f_re_mean);

%for LHS
f_tau_lhs =  (N_tau_lhs - N_orig_lhs) / (N_both_lhs - N_orig_lhs);
f_re_lhs =  (N_re_lhs - N_orig_lhs) / (N_both_lhs - N_orig_lhs);
f_both_lhs =  (N_both_lhs - N_orig_lhs) / (N_both_lhs- N_orig_lhs);

fprintf(1,'f_both_lhs=%f f_tau_lhs=%f f_re_lhs=%f f_rem_lhs=%f\n',f_both_lhs,f_tau_lhs,f_re_lhs,1-f_tau_lhs-f_re_lhs);

clear override_attr
catch error_attr
   clear override_attr
   rethrow(error_attr)
end




