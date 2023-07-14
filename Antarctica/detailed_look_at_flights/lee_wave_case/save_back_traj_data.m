filedir_bt = 'Y:\BAS_flights\flight102\Matlab_back_traj_data\';
savename = 'W plus 1.5';
savename = 'W plus 1.25';
savename = 'W plus 1';
savename = 'W plus 1.5';
savename = 'W plus 1.5 fixed LWC'; %fixed as in fixed the bug
%savename = 'W plus 1.5 fixed LWC liq dew point'; %using 'liq' flag to calculate qv from the dew point (humicap)
savename='tot_back_traj_Wplus_1.5_qvliq_cas_RH_thresh0.85';
%savename='tot_back_traj_Wplus_0_qvliq_cas_RH_thresh0.85';
%savename='tot_back_traj_Wplus_0_qvliq_cas_RH_thresh0.85_U_22';
%savename='tot_back_traj_Wplus_1.25_qvliq_cas_RH_thresh0.85';
savename='tot_back_traj_Wplus_0_upward_tilt_68.2_constantU_15';

save_load='load';
save_load='save';

switch save_load
    case 'save'

        save([filedir_bt savename '.mat'],'X2','Z2','Zsound_pressure','tot3','equiv3','temp3',...
            'zz','tot2','equiv2','temp2');

    case 'load'
        load([filedir_bt savename '.mat'],'X2','Z2','Zsound_pressure','tot3','equiv3','temp3',...
            'zz','tot2','equiv2','temp2');


end

disp('Done save/load');