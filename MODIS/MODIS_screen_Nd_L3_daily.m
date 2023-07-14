%save_file = ['/home/disk/eos1/d.grosvenor/modis_work/saved_data_L3/Nd_from_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_meanSZA_65/'...
%    modis_data_case{1}{1} '.mat'];


thresh_SZA=[45 55];
thresh_SZA=[45 55];
thresh_SZA=[0 90];
thresh_SZA=[0 65];
%thresh_SZA=[50 55];
%thresh_SZA=[75 82];
%thresh_SZA=[46 56];
%thresh_SZA=[75 90];

switch CF_str
    case '0.8'
        thresh_CF=[0.8 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
    case '0.0'
        thresh_CF=[-0.01 1.000001];
end
%thresh_CF=[0.99 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.01 0.8];
%thresh_CF=[0.0 1.000001];
thresh_NP=10;
thresh_NP=50;
%thresh_sensZA=45;
%thresh_sensZA=50;
%thresh_sensZA=[0 41.4];
thresh_sensZA=[0 90];
%thresh_sensZA=[55 65];
%thresh_sensZA=[30 70];
%thresh_sensZA=58;

%thresh_CTH = [-20 20];
thresh_CTH = [-20 3.2];

%thresh_AZ = [50 130]; %not used

thresh_relAZ = [0 180];

%thresh_CTT = [273-5 273+100];
%thresh_CTT = [273 273+100];
%thresh_CTT = [273 273+100];
thresh_CTT = [273-100 273+100];
%thresh_CTP = 800; %Cloud top pressure (hPa)

%    thresh_sensZA=80;
%thresh_maxSZA=200;

% thresh_Nd_per_error = 100;
% thresh_Reff_per_error = 50;
% thresh_Reff_abs_error = 4;
% thresh_Reff = 30;

thresh_stdW = [0 1e9];




%screen_type='NP + CF + MAX sensZA';
%                                    screen_type='NP + MAX sensZA';
screen_type='NP + CF + MEAN sensZA';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT + noice';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT';
screen_type ='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT';
%screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH';
screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';


%Run screening script
iplot_global=1;
igcm_screen=0;
modisL3_screening_timeseries3
%This calculates the indices to screen with



vars={'N_time3','N_time3_37','Cloud_Optical_Thickness_Liquid_Mean.timeseries3','Cloud_Effective_Radius_Liquid_Mean.timeseries3','Cloud_Effective_Radius_37_Liquid_Mean.timeseries3',...
    'Cloud_Fraction_Liquid.timeseries3'};
for ivar=1:length(vars)
    var = vars{ivar};
    eval_str = [var '_screened = ' var ';'];
    eval(eval_str);
    eval_str = [var '_screened(ihtot)=NaN;'];
    eval(eval_str);
    
    
    
%     var2=[var '_screened'];
    
%     if ivar==1
%         iappend=0;
%     else
%         iappend=1;
%     end
%     eval_str = ['save_vars_mat_func(save_file,var,' var2 ',iappend);']
%     eval(eval_str);
end


