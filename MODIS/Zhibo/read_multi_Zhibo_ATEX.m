try
case_str_zhibo = 'ATEX';

if ~exist('ioverride_iview') | ioverride_iview==0
    iview2=1; %counter for looping through the views
    iview=7; %index for viewing geometry (7 means VZA=0, i.e. nadir)
% iview=12 and 13 are messed up (Zhibo mentioned this in an email)
% So can't use backscatter for 50 or 60 degrees.
    iplot_zhibo=1;
end

if iview2==1  %only to read and calc means once (does for all view angles anyway)

    time_strs = {'003734','004155','004576','006260','006681'}; %Think these are (variable) timesteps
    clear times
    for i=1:length(time_strs)
        times(i) = str2num(time_strs{i});
    end

    %ATEX_CCN40_dharma_003734_retrieval_results_sz0.5000.nc
    filedir='/home/disk/eos1/d.grosvenor/modis_work/Zhibo_Marshak_work/Zhibo_netCDF_files_July2013/LES_Data/';
    %See end for a file list


    %SZA=20
    %DYCOMS file :-
    %file_read = 'Ackerman_DYCOMS2_25bins_dharma_003043_retrieval_results_sz0.9396.nc';
    %ATEX :-
    file_read = ['ATEX_CCN40_dharma_TIME_retrieval_results_sz0.9396.nc'];
    read_Zhibo_DYCOMS
    %filetag = 'SZA_20';
    SZA_20 = LES_struct;
    SZA_20_times = LES_struct_times;        


    %SZA=60
    %file_read = 'Ackerman_DYCOMS2_25bins_dharma_003043_retrieval_results_sz0.1736.nc';
    file_read = ['ATEX_CCN40_dharma_TIME_retrieval_results_sz0.5000.nc'];
    %filetag = 'SZA_80';
    read_Zhibo_DYCOMS
    SZA_80 = LES_struct;
    SZA_80_times = LES_struct_times;

end

%Also want to check that the reference re (1D, low SZA) does not change
%with vza

prc_bias_overall_sza20_3D_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(end,iview)./ SZA_20.timemean_Re21_1D_vs_res(1,iview) - 1)
prc_bias_sub_sza20_1D_21 = 100* ( SZA_20.timemean_Re21_1D_vs_res(end,iview)./ SZA_20.timemean_Re21_1D_vs_res(1,iview) - 1)
prc_bias_sub_sza20_3D_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(end,iview)./ SZA_20.timemean_Re21_3D_vs_res(1,iview) - 1)
prc_bias_3Dhet_sza20_100m_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(1,iview)./ SZA_20.timemean_Re21_1D_vs_res(1,iview) - 1)
prc_bias_3Dhet_sza20_800m_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(end,iview)./ SZA_20.timemean_Re21_1D_vs_res(end,iview) - 1)

prc_bias_overall_sza20_3D_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(end,iview)./ SZA_20.timemean_Re37_1D_vs_res(1,iview) - 1)
prc_bias_sub_sza20_1D_37 = 100* ( SZA_20.timemean_Re37_1D_vs_res(end,iview)./ SZA_20.timemean_Re37_1D_vs_res(1,iview) - 1)
prc_bias_sub_sza20_3D_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(end,iview)./ SZA_20.timemean_Re37_3D_vs_res(1,iview) - 1)
prc_bias_3Dhet_sza20_100m_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(1,iview)./ SZA_20.timemean_Re37_1D_vs_res(1,iview) - 1)
prc_bias_3Dhet_sza20_800m_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(end,iview)./ SZA_20.timemean_Re37_1D_vs_res(end,iview) - 1)


%high SZA (labelled SZA=80, but was 60 for ATEX)
prc_bias_overall_sza80_3D_21 = 100* ( SZA_80.timemean_Re21_3D_vs_res(end,iview)./ SZA_80.timemean_Re21_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza80_1D_21 = 100* ( SZA_80.timemean_Re21_1D_vs_res(end,iview)./ SZA_80.timemean_Re21_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza80_3D_21 = 100* ( SZA_80.timemean_Re21_3D_vs_res(end,iview)./ SZA_80.timemean_Re21_3D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza80_100m_21 = 100* ( SZA_80.timemean_Re21_3D_vs_res(1,iview)./ SZA_80.timemean_Re21_1D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza80_800m_21 = 100* ( SZA_80.timemean_Re21_3D_vs_res(end,iview)./ SZA_80.timemean_Re21_1D_vs_res(end,iview) - 1);

prc_bias_overall_sza80_3D_37 = 100* ( SZA_80.timemean_Re37_3D_vs_res(end,iview)./ SZA_80.timemean_Re37_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza80_1D_37 = 100* ( SZA_80.timemean_Re37_1D_vs_res(end,iview)./ SZA_80.timemean_Re37_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza80_3D_37 = 100* ( SZA_80.timemean_Re37_3D_vs_res(end,iview)./ SZA_80.timemean_Re37_3D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza80_100m_37 = 100* ( SZA_80.timemean_Re37_3D_vs_res(1,iview)./ SZA_80.timemean_Re37_1D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza80_800m_37 = 100* ( SZA_80.timemean_Re37_3D_vs_res(end,iview)./ SZA_80.timemean_Re37_1D_vs_res(end,iview) - 1);


ref_re_21 = SZA_20.timemean_Re21_1D_vs_res(1,iview);
ref_re_37 = SZA_20.timemean_Re37_1D_vs_res(1,iview);

%x=[50 100 400 800];
x=[100 400 800];

if iplot_zhibo==1



figure
plot(x,SZA_20.timemean_Re21_1D_vs_res(:,iview),'rs-');
hold on
plot(x,SZA_20.timemean_Re21_3D_vs_res(:,iview),'bo-');
plot(x,SZA_80.timemean_Re21_1D_vs_res(:,iview),'rs--');
plot(x,SZA_80.timemean_Re21_3D_vs_res(:,iview),'bo--');

legend({'1D, SZA=20','3D, SZA=20','1D, SZA=60','3D, SZA=60'});
ylabel('Mean R_e (\mum)');
xlabel('Averaging scale (m)');
title('2.1 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);


figure
plot(x,SZA_20.timemean_Re37_1D_vs_res(:,iview),'rs-');
hold on
plot(x,SZA_20.timemean_Re37_3D_vs_res(:,iview),'bo-');
plot(x,SZA_80.timemean_Re37_1D_vs_res(:,iview),'rs--');
plot(x,SZA_80.timemean_Re37_3D_vs_res(:,iview),'bo--');

legend({'1D, SZA=20','3D, SZA=20','1D, SZA=60','3D, SZA=60'});
ylabel('Mean R_e (\mum)');
xlabel('Averaging scale (m)');
title('3.7 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);


% ---- Plotting 2.1 and 3.7 minus their correct ref value (1D SZA=20 at high res)
% N.B. this is different for 2.1 and 3.7, presumably due to vertical
% variation (3.7 is larger than 2.1, which makes sense for an adiabatic
% cloud)
figure
plot(x,SZA_20.timemean_Re21_1D_vs_res(:,iview) - SZA_20.timemean_Re21_1D_vs_res(1,iview),'rs-');
hold on
plot(x,SZA_20.timemean_Re21_3D_vs_res(:,iview) - SZA_20.timemean_Re21_1D_vs_res(1,iview),'bo-');
plot(x,SZA_80.timemean_Re21_1D_vs_res(:,iview) - SZA_20.timemean_Re21_1D_vs_res(1,iview),'rs--');
plot(x,SZA_80.timemean_Re21_3D_vs_res(:,iview) - SZA_20.timemean_Re21_1D_vs_res(1,iview),'bo--');

grid
legend({'1D, SZA=20','3D, SZA=20','1D, SZA=60','3D, SZA=60'});
ylabel('Mean R_e minus ref (\mum)');
xlabel('Averaging scale (m)');
title('2.1 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);

figure
plot(x,SZA_20.timemean_Re37_1D_vs_res(:,iview) - SZA_20.timemean_Re37_1D_vs_res(1,iview),'rs-');
hold on
plot(x,SZA_20.timemean_Re37_3D_vs_res(:,iview) - SZA_20.timemean_Re37_1D_vs_res(1,iview),'bo-');
plot(x,SZA_80.timemean_Re37_1D_vs_res(:,iview) - SZA_20.timemean_Re37_1D_vs_res(1,iview),'rs--');
plot(x,SZA_80.timemean_Re37_3D_vs_res(:,iview) - SZA_20.timemean_Re37_1D_vs_res(1,iview),'bo--');

grid
legend({'1D, SZA=20','3D, SZA=20','1D, SZA=60','3D, SZA=60'});
ylabel('Mean R_e minus ref (\mum)');
xlabel('Averaging scale (m)');
title('3.7 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);



% ---- re2.1 minus 3.7 -----

figure
plot(x,SZA_20.timemean_Re21_1D_vs_res(:,iview) - SZA_20.timemean_Re37_1D_vs_res(:,iview),'rs-');
hold on
plot(x,SZA_20.timemean_Re21_3D_vs_res(:,iview) - SZA_20.timemean_Re37_3D_vs_res(:,iview),'bo-');
plot(x,SZA_80.timemean_Re21_1D_vs_res(:,iview) - SZA_80.timemean_Re37_1D_vs_res(:,iview),'rs--');
plot(x,SZA_80.timemean_Re21_3D_vs_res(:,iview) - SZA_80.timemean_Re37_3D_vs_res(:,iview),'bo--');

legend({'1D, SZA=20','3D, SZA=20','1D, SZA=60','3D, SZA=60'});
ylabel('Difference in Mean R_e (\mum)');
xlabel('Averaging scale (m)');
title('2.1 - 3.7 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);

end


%File list :-
%filedir='/home/disk/eos1/d.grosvenor/modis_work/Zhibo_Marshak_work/Zhibo_netCDF_files_July2013/LES_Data/';

%-rw-r--r-- 1 d.grosvenor atgstaff 74725816 Jul 31  2013 ATEX_CCN40_dharma_003734_Physics_Optics.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_003734_retrieval_results_sz0.5000.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_003734_retrieval_results_sz0.9396.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 74725816 Jul 31  2013 ATEX_CCN40_dharma_004155_Physics_Optics.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_004155_retrieval_results_sz0.5000.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_004155_retrieval_results_sz0.9396.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 74725816 Jul 31  2013 ATEX_CCN40_dharma_004576_Physics_Optics.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_004576_retrieval_results_sz0.5000.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_004576_retrieval_results_sz0.9396.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 74725816 Jul 31  2013 ATEX_CCN40_dharma_006260_Physics_Optics.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_006260_retrieval_results_sz0.5000.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_006260_retrieval_results_sz0.9396.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 74725816 Jul 31  2013 ATEX_CCN40_dharma_006681_Physics_Optics.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_006681_retrieval_results_sz0.5000.nc
% -rw-r--r-- 1 d.grosvenor atgstaff 13955576 Jul 31  2013 ATEX_CCN40_dharma_006681_retrieval_results_sz0.9396.nc
% -rw-r--r-- 1 d.grosvenor atgstaff     1262 Mar 26  2014 ncdump_Physical_Optics.out
% -rw-r--r-- 1 d.grosvenor atgstaff    10499 Mar 26  2014 ncdump_retrieval_results.out

clear ioverride_iview
catch error_zhibo
    clear ioverride_iview
    rethrow(error_zhibo);
end