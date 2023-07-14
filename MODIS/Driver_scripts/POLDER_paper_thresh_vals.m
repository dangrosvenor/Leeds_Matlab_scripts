thresh_SZA=[75 90]; %upper range for SZA paper
thresh_SZA=[50 55];
%thresh_SZA=[45 55];
thresh_SZA=[-1 90];
%thresh_SZA=[0 70];
%thresh_SZA=[70 90];
%thresh_SZA=[0 65];
%thresh_SZA=[50 55];  %lower range for SZA paper
%thresh_SZA=[48.78 52.82];
%thresh_SZA=[75 82];
%thresh_SZA=[46 56];
%thresh_SZA=[55 60];
%thresh_SZA=[50 55];
%thresh_SZA=[45 60];
%thresh_SZA=[45 85];

thresh_CF=[0.799 1.0000001];
%thresh_CF=[0.47 1.0000001];
%thresh_CF=[0.99 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.1 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.99 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
%thresh_CF=[0.01 0.8];
thresh_CF=[-0.1 1.000001];

thresh_NP=0;
%thresh_NP=10;
thresh_NP=25;
thresh_NP=50;

thresh_ndays=15;

thresh_zeroCF = 0.05;



thresh_sensZA=45;
%thresh_sensZA=50;
thresh_sensZA=[0 20];
thresh_sensZA=[0 41.4];
%thresh_sensZA=[41.4 90];
thresh_sensZA=[-1 90];
%thresh_sensZA=[0 65];
%thresh_sensZA=[55 65];
%thresh_sensZA=[30 70];
%thresh_sensZA=58;

thresh_CTH = [-20 1e9];
%thresh_CTH = [-20 2];
thresh_CTH = [-0.01 3.2];

thresh_AZ = [50 130];

thresh_relAZ = [-1 181];
%thresh_relAZ = [0 80];
%thresh_relAZ = [100 180];

thresh_CTT = [273 273+100];
thresh_CTT = [268 273+100];
%thresh_CTT = [173 273+100];
thresh_CTT = [273 273+100];
%thresh_CTT = [273 273+100];
%thresh_CTT = [273-100 273+100];
%thresh_CTT = [265 273+100];
%thresh_CTT = [250 273+100];
thresh_CTP = 800; %Cloud top pressure (hPa)

thresh_reff = [0 30];
%thresh_reff = [0 14];

thresh_sigCTT = [0 1e9];

%    thresh_sensZA=80;
thresh_maxSZA=200;

thresh_Nd_per_error = 100;
thresh_Reff_per_error = 50;
thresh_Reff_abs_error = 4;
thresh_Reff = 30;

thresh_tau = [-1 300];

%this is set up for homog = (LWP/std_LWP).^2  - so is HOMOGENEITY. Lower
%numbers are more heterogeneous
thresh_stdW = [25 1e9];
thresh_stdW = [0 1e9];


minfrac_CF = 0.9; 
%minfrac_CF = 0.99;  
%minfrac_CF = 0.79;  
%minfrac_CF = 0.47;  
%minfrac_CF = 0.1;  
%minfrac_CF = 0.0;  
%minimum fraction of the sampled points that had successful cloudy/clear/phase
      %determination (i.e. Npix/Nptot_mockL3 =
      %Cloud_Fraction_Liquid_Pixel_Counts./Cloud_Fraction_Liquid./Total_pixels
      % - restriction (2) as presented in the paper
  

minfrac_NpNd = 0.9;
%minfrac_NpNd = 0.99;
%minfrac_NpNd = 0.47;
%minfrac_NpNd = 0.1;
%minfrac_NpNd = 0.0;
%Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3./Cloud_Fraction_Liquid_Pixel_Counts.timeseries3
%Fraction of points that remain after all previous filtering for
%which we have an Nd retrieval. Restriction (4) in the SZA paper.

thresh_dreff = 0.6; %maximum allowed difference between 1.6, 2.1 and 3.7 um retrievals
thresh_dreff = 1e9; %maximum allowed difference between 1.6, 2.1 and 3.7 um retrievals - is set as a percentage at the moment
%thresh_dreff = 10; %maximum allowed difference between 1.6, 2.1 and 3.7 um retrievals - is set as a percentage at the moment

%max allowed height of the upper cloud layer
thresh_maxlayerH = [0 1e5];
%thresh_maxlayerH = [1800 1e5];
thresh_maxlayerH = [0 800];
thresh_nlayers = [1 1]; %max number of allowed layers
%thresh_nlayers = [1 4]; %max number of allowed layers

thresh_calipso_highCF = [-0.01 0.3]; %desired range of Calipso mid+high cloud