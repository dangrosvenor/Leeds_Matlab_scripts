function [dat,dat_annual_ens_mean,dat_ens_mean,dat_annual_ens] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP)



mat_obj = matfile(load_file);
if isingle_ens==1
    dat = squeeze(fscale * mat_obj.dat_ens(iens,t_inds,:,:)); %monthly, individual ensemble members
    switch var_DAMIP
        case 'rsut'; %For rsut we also need this (annual mean for iens)
            dat_annual_ens_mean = squeeze( mat_obj.dat_annual_ens(iens,t_inds_annual,:,:) ); %annual, individual ensemble members
    end
else
    %1980 is all 12 months of all 165 years.
    dat = fscale * mat_obj.dat_ens_mean(t_inds,:,:); %convert from % to 0 to 1 form.
    switch var_DAMIP
        case 'rsut'; %For rsut we also need this (annual mean for iens)
            dat_annual_ens_mean = mat_obj.dat_annual(t_inds_annual,:,:);
    end
end

switch var_DAMIP
    case 'rsut';
        dat_ens_mean = mat_obj.dat_ens_mean(t_inds,:,:); %monthly ensemble mean
        dat_annual_ens = mat_obj.dat_annual_ens(:,t_inds_annual,:,:); %annual individual ensemble values
end


switch expt_str
    case 'AerChemMIP_hist-AerProxy'
        mat_obj2 = matfile(load_file2);
        if isingle_ens==1
            dat2 = squeeze(fscale * mat_obj2.dat_ens(iens,t_inds,:,:));
            switch var_DAMIP
                case 'rsut'; %For rsut we also need this (annual mean for iens)
                    dat_annual_ens_mean2 = squeeze( mat_obj2.dat_annual_ens(iens,t_inds_annual,:,:) ); %annual, individual ensemble members
            end
        else
            %1980 is all 12 months of all 165 years.
            dat2 = fscale * mat_obj2.dat_ens_mean(t_inds,:,:); %convert from % to 0 to 1 form.
            switch var_DAMIP
                case 'rsut'; %For rsut we also need this (annual mean for iens)
                    dat_annual_ens_mean2 = mat_obj2.dat_annual(t_inds_annual,:,:);
            end
        end
        
        switch var_DAMIP
            case 'rsut';
                dat_ens_mean2 = mat_obj2.dat_ens_mean(t_inds,:,:); %monthly ensemble mean
                dat_annual_ens2 = mat_obj2.dat_annual_ens(:,t_inds_annual,:,:); %annual individual ensemble values
        end
                       
        
        %Take average over first 20 years.
        tav_inds = [1:20*12];
        [dat] = sw_diff_AerChemMIP(dat,dat2,tav_inds);
        
        switch var_DAMIP
            case 'rsut';
                %Take average over first 20 years.
                tav_inds = [1:20];
                [dat_annual_ens_mean] = sw_diff_AerChemMIP(dat_annual_ens_mean,dat_annual_ens_mean2,tav_inds);
                
                 tav_inds = [1:20*12];
                [dat_ens_mean] = sw_diff_AerChemMIP(dat_ens_mean,dat_ens_mean2,tav_inds);
                
                tav_inds = [1:20];
                for iens=1:size(dat_annual_ens,1)
                    [dat_annual_ens(iens,:,:,:)] = sw_diff_AerChemMIP(squeeze(dat_annual_ens(iens,:,:,:)),squeeze(dat_annual_ens2(iens,:,:,:)),tav_inds);
                end
                
     
        end
        
        
%         dat_t0_01 = meanNoNan(dat(1:20*12,:,:),1);
%         %Do same for dat2
%         dat_t0_02 = meanNoNan(dat2(1:20*12,:,:),1);
%         %take the average
%         dat0 = repmat( 0.5*(dat_t0_01+dat_t0_02) , [1 1 size(dat,1)]);
%         dat0 = permute(dat0,[3 1 2]);
%         
%         %load_file is the control and load_file2 is piAer (GHG-only proxy)
%         %So do control minus GHG for the deltas in the aerosol proxy
%         %The deltas are relative to a common baseline (1850-1870 average) when it is
%         %assumed that the control and piAer runs are the same.
%         
%         dat = dat0 + dat - dat2;
%         %dat = SW_up_TOA_dat_ens_monthly  - SW_up_TOA_dat_ens_monthly2;
        
end


function [dat_out] = sw_diff_AerChemMIP(dat_01,dat_02,tav_inds)
%Take average over first 20 years.
dat_t0_01 = meanNoNan(dat_01(tav_inds,:,:),1);
%Do same for dat2
dat_t0_02 = meanNoNan(dat_02(tav_inds,:,:),1);
%take the average
dat0 = repmat( 0.5*(dat_t0_01+dat_t0_02) , [1 1 size(dat_01,1)]);
dat0 = permute(dat0,[3 1 2]);

%load_file is the control and load_file2 is piAer (GHG-only proxy)
%So do control minus GHG for the deltas in the aerosol proxy
%The deltas are relative to a common baseline (1850-1870 average) when it is
%assumed that the control and piAer runs are the same.

dat_out = dat0 + dat_01 - dat_02;
%dat = SW_up_TOA_dat_ens_monthly  - SW_up_TOA_dat_ens_monthly2;
        
        


