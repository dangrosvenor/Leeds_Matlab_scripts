try
case_str_zhibo = 'DYCOMS';

if ~exist('ioverride_iview') | ioverride_iview==0
    iview2=1; %counter for looping through the views
    iview=7; %index for viewing geometry (7 means VZA=0, i.e. nadir)
% iview=12 and 13 are messed up (Zhibo mentioned this in an email)
% So can't use backscatter for 50 or 60 degrees.
    iplot_zhibo=1;
end

if iview2==1  %only to read and calc means once (does for all view angles anyway)

    time_strs = {'003043','003784','004637','005416','006263','007138','008036','008884'};
%    time_strs = {'003043'};
    clear times
    for i=1:length(time_strs)
        times(i) = str2num(time_strs{i});
    end

    filedir='/home/disk/eos1/d.grosvenor/modis_work/Zhibo_Marshak_work/Ackermann_DYCOMS2_25bins_from_Zhibo/';

    %SZA=20
    file_read = 'Ackerman_DYCOMS2_25bins_dharma_TIME_retrieval_results_sz0.9396.nc';
    read_Zhibo_DYCOMS
    %filetag = 'SZA_20';
    SZA_20 = LES_struct;
    SZA_20_times = LES_struct_times;    


    %SZA=80
    file_read = 'Ackerman_DYCOMS2_25bins_dharma_TIME_retrieval_results_sz0.1736.nc';
    %filetag = 'SZA_80';
    read_Zhibo_DYCOMS
    SZA_80 = LES_struct;
    SZA_80_times = LES_struct_times; 

end

%low SZA (labelled SZA=20)
prc_bias_overall_sza20_3D_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(end,iview)./ SZA_20.timemean_Re21_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza20_1D_21 = 100* ( SZA_20.timemean_Re21_1D_vs_res(end,iview)./ SZA_20.timemean_Re21_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza20_3D_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(end,iview)./ SZA_20.timemean_Re21_3D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza20_100m_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(1,iview)./ SZA_20.timemean_Re21_1D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza20_800m_21 = 100* ( SZA_20.timemean_Re21_3D_vs_res(end,iview)./ SZA_20.timemean_Re21_1D_vs_res(end,iview) - 1);

prc_bias_overall_sza20_3D_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(end,iview)./ SZA_20.timemean_Re37_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza20_1D_37 = 100* ( SZA_20.timemean_Re37_1D_vs_res(end,iview)./ SZA_20.timemean_Re37_1D_vs_res(1,iview) - 1);
prc_bias_sub_sza20_3D_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(end,iview)./ SZA_20.timemean_Re37_3D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza20_100m_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(1,iview)./ SZA_20.timemean_Re37_1D_vs_res(1,iview) - 1);
prc_bias_3Dhet_sza20_800m_37 = 100* ( SZA_20.timemean_Re37_3D_vs_res(end,iview)./ SZA_20.timemean_Re37_1D_vs_res(end,iview) - 1);


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


x=[50 100 400 800];
%x=[100 400 800];


if iplot_zhibo==1



figure
plot(x,SZA_20.mean_Re21_1D_vs_res(:,iview),'rs-');
hold on
plot(x,SZA_20.mean_Re21_3D_vs_res(:,iview),'bo-');
plot(x,SZA_80.mean_Re21_1D_vs_res(:,iview),'rs--');
plot(x,SZA_80.mean_Re21_3D_vs_res(:,iview),'bo--');

legend({'1D, SZA=20','3D, SZA=20','1D, SZA=80','3D, SZA=80'});
ylabel('Mean R_e (\mum)');
xlabel('Averaging scale (m)');
title('2.1 \mum, DYCOMS');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);

figure
plot(x,SZA_20.mean_Re37_1D_vs_res(:,iview),'rs-');
hold on
plot(x,SZA_20.mean_Re37_3D_vs_res(:,iview),'bo-');
plot(x,SZA_80.mean_Re37_1D_vs_res(:,iview),'rs--');
plot(x,SZA_80.mean_Re37_3D_vs_res(:,iview),'bo--');

legend({'1D, SZA=20','3D, SZA=20','1D, SZA=80','3D, SZA=80'});
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
plot(x,SZA_20.mean_Re21_1D_vs_res(:,iview) - SZA_20.mean_Re21_1D_vs_res(1,iview),'rs-');
hold on
plot(x,SZA_20.mean_Re21_3D_vs_res(:,iview) - SZA_20.mean_Re21_1D_vs_res(1,iview),'bo-');
plot(x,SZA_80.mean_Re21_1D_vs_res(:,iview) - SZA_20.mean_Re21_1D_vs_res(1,iview),'rs--');
plot(x,SZA_80.mean_Re21_3D_vs_res(:,iview) - SZA_20.mean_Re21_1D_vs_res(1,iview),'bo--');

grid
legend({'1D, SZA=20','3D, SZA=20','1D, SZA=80','3D, SZA=80'});
ylabel('Mean R_e minus ref (\mum)');
xlabel('Averaging scale (m)');
title('2.1 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);

figure
plot(x,SZA_20.mean_Re37_1D_vs_res(:,iview) - SZA_20.mean_Re37_1D_vs_res(1,iview),'rs-');
hold on
plot(x,SZA_20.mean_Re37_3D_vs_res(:,iview) - SZA_20.mean_Re37_1D_vs_res(1,iview),'bo-');
plot(x,SZA_80.mean_Re37_1D_vs_res(:,iview) - SZA_20.mean_Re37_1D_vs_res(1,iview),'rs--');
plot(x,SZA_80.mean_Re37_3D_vs_res(:,iview) - SZA_20.mean_Re37_1D_vs_res(1,iview),'bo--');

grid
legend({'1D, SZA=20','3D, SZA=20','1D, SZA=80','3D, SZA=80'});
ylabel('Mean R_e (\mum)');
xlabel('Averaging scale (m)');
title('3.7 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);



% ---- re2.1 minus 3.7 -----

figure
plot(x,SZA_20.mean_Re21_1D_vs_res(:,iview) - SZA_20.mean_Re37_1D_vs_res(:,iview),'rs-');
hold on
plot(x,SZA_20.mean_Re21_3D_vs_res(:,iview) - SZA_20.mean_Re37_3D_vs_res(:,iview),'bo-');
plot(x,SZA_80.mean_Re21_1D_vs_res(:,iview) - SZA_80.mean_Re37_1D_vs_res(:,iview),'rs--');
plot(x,SZA_80.mean_Re21_3D_vs_res(:,iview) - SZA_80.mean_Re37_3D_vs_res(:,iview),'bo--');

legend({'1D, SZA=20','3D, SZA=20','1D, SZA=80','3D, SZA=80'});
ylabel('Difference in Mean R_e (\mum)');
xlabel('Averaging scale (m)');
title('2.1 - 3.7 \mum');
set(gca,'xlim',[0 1000]);
%set(gca,'xlim',[0 1000]);

end

clear ioverride_iview
catch error_zhibo
    clear ioverride_iview
    rethrow(error_zhibo);
end