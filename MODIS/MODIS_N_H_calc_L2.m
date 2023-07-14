%calculate N and H from MOD06 optical thickness, effective radius and cloud fraction

use_QA_vals=0;

if use_QA_vals==1
    %using MOD06 only here
    
    fprintf(1,'\n**** USING QA means for tau and reff ***\n');
    title_info = 'QA Means';

   
    tau = Cloud_Optical_Thickness_Liquid_QA_Mean.data;
    %tau(tau>50)=NaN; %remove large tau values
    %fprintf(1,'\n**** WARNING - am removing data when tau >50 ***\n');

    reff = Cloud_Effective_Radius_Liquid_QA_Mean.data*1e-6; %convert to metres
    WMOD=Cloud_Water_Path_Liquid_QA_Mean.data/1000; %convert to kg/m2
else
    %using MOD06 only here

    tau = Cloud_Optical_Thickness_Liquid_Mean.data;
    %tau(tau>50)=NaN; %remove large tau values
    %fprintf(1,'\n**** WARNING - am removing data when tau >50 ***\n');

    reff = Cloud_Effective_Radius_Liquid_Mean.data*1e-6; %convert to metres
    WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
end

    sangle = Scattering_Angle_Mean.data;
    cf = Cloud_Fraction_Liquid.data;
    
    
    


if ~exist('set_MODIS_NH_flags') | set_MODIS_NH_flags==0

    Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
    %Wflag='MODIS'; %use the MODIS LWP

else
    clear set_MODIS_NH_flags
end

[N,H,W,k,Q,cw]=MODIS_N_H_func(tau,reff,Wflag,WMOD);



