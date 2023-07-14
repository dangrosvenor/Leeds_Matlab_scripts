% This is run several times from :- ACSIS_Robson_paper_TABLE_stats_noobs3_FUNC.m

%deltaT=1;
%var_str_tab='totCF'; fscale=1e5;

%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
%trends = load(tr_file);

%will assume the file is loaded

switch iregion{1}
    case '0'
        table_str='table_vals_region0.';
    otherwise
        table_str='table_vals.';
end
    


var_str_tab2=[var_str_tab run_str]; %e.g. CFHAD run_str='HAD'
trend_str=[var_str_tab2 'trend']; %e.g. CFtrend

if iobs==0
    eval([table_str trend_str period_str ' = deltaT*fscale*trends.trend_dat_box' hist_str '{1,itr}.coeffs(2);']); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
else
    %Obs are only for one period here (usually itr=3).
    eval([table_str trend_str period_str ' = deltaT*fscale*trends.trend_dat_box_obs2{1}.coeffs(2);']); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.    
end
eval([table_str trend_str period_str 'exp = -log10(fscale);']); %the exponent for the data when quoting. E.g. Y x 10^{exp}
eval([table_str trend_str period_str 'exp_precision = ''%d'';']);
eval([table_str trend_str period_str 'un = deltaT*fscale*trends.trend_dat_box' hist_str '{1,itr}.uncer_max;']); %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
if ido_ens==1
    Nens = size(trends.trend_dat_box_ens,3);
    clear dat
    for iens=1:Nens
        dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values
    end
    [minval imin]=min(dat);
    %ensemble member with the smallest trend
    eval([table_str trend_str period_str 'min = deltaT*minval;']);
    %The uncertainyy from that ensemble member.
    eval([table_str trend_str period_str 'minUN = deltaT*fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;']);
    [maxval imax]=max(dat);
    %ensemble member with the smallest trend
    eval([table_str trend_str period_str 'max = deltaT*maxval;']);
    %The uncertainyy from that ensemble member.
    eval([table_str trend_str period_str 'maxUN = deltaT*fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;']);
end
%Not sure if the following is correct? Since i0 is just the single index of
%the start of the time period?
eval(['x5 = meanNoNan(trends.dat_annual_box_ukesm' hist_str '(i0:i0+4),2);']);
eval([table_str var_str_tab2 'avFirstFive' period_str ' = x5;']);
eval([table_str trend_str 'Pct' period_str ' = 100*' table_str trend_str period_str '/' table_str var_str_tab2 'avFirstFive' period_str '/deltaT;']);
eval([table_str trend_str 'Pct' period_str '_precision=''%.9f'';']);

%% Caclulation of fc and LWPic change expected from change in Nd based on G20 paper (i.e., the local aerosol effect only).

%dNd = 56;
dNd = 59; %(dNd from \NASWdNdregionalmeanprcbias in /home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/Percentage_change_in_Nd_(PD_minus_PI)_ocean_only,_no_sea-ice_STATS.tex)
%S_CF = 0.0232; %dfc/dNd from G20 paper in percentage form (%change in Fc per % change in Nd).
S_CF = 1.1/dNd; %dfc (dfc from \NASWdlowCFregionalmeanprcbias in /home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/Percentage_change_in_cloud_fraction_(PD_minus_PI)_ocean_only,_no_sea-ice_STATS.tex)
%S_LWPic = 0.0155;
S_LWPic = 0.82 / dNd;
%S_SW =   0.0446; % dSW = 2.5%, dNd = 56%
S_SW =   2.7 / dNd; % /home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/Percentage_change_in_SW_TOA_(PD_minus_PI)_ocean_only,_no_sea-ice_STATS.tex
% 2.2 % change from indirect effect - will assume that direct effect
% makes up the difference from the total SW - although presumably do have
% an ARI dSW value - however, perhaps harder to scale? Could use AOD
% perhaps. Or just scale by Nd too? I think the direct+indirect combined to
% something close to the total SW anyway.
S_SWind = 2.2 / dNd; %/home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/Percentage_change_in_SW_TOA_due_to_indirect_forcing_(PD_minus_PI)_ocean_only,_no_sea-ice_STATS.tex
%For direct effect :- % SW change was 0.57%, dSW was 0.53 W/m2; /home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/Percentage_change_in_SW_TOA_due_to_direct_forcing_(PD_minus_PI)_ocean_only,_no_sea-ice_STATS.tex
%for indirect dSW was 2.1 W/m2. For total SW dSW was 2.5 W/m2. So
%direct+indirect = 2.63 W/m2, which is just a little higher than the actual
%total. Calculating the direct effect from the difference gives 2.5-2.1 =
%0.4 W/m2, which is 100*(1-0.4/0.53) = 25% too low.
%


% Useful to know the contribution of Nd etc. to SWTOA indirect effect
% relative to the total (CF+Nd+LWPic).
%		§ Further down from the above UM_ACSIS*.m there are sections that calculate the contributions for TOA (previously did it for the surface).
% 			§ Then look in the table .tex files for the regional avearges. E.g. :-
% 				? /home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/Change_in_SWTOA_due_to_Nd_PI+PD_average_ocean_only,_no_sea-ice_STATS.tex
% 				? To get the \NASWSWTOAforcingNdPIPDmetrue etc. values.
% 				? Shows that for the NASW region when doing the ocean and sea-ice screening the Nd contribution was -1.1, CF was -0.46 and LWPic was -0.66 W/m2.
% 				? So Nd contribution percentage is -1.1/(1.1+0.46+0.66) = -0.4955
% 				? So, around 50%.

%Also do for the other time period from the ACSIS paper.

switch period_str
    case 'PC'
        exp_Nd=NaN;
        T_Nd=NaN;
        
    otherwise
        switch iregion{1}
            case '4'                                
                switch run_str
                    %case {'HistAer2','AerChemAerosol','AerChemGHG','UKESMAerChemMIP'}
                    case {'AerChemAerosol','AerChemGHG','UKESMAerChemMIP',}
                        run_str_temp = 'HistAer'; %until I get Nd for hist-aer2
                        exp_Nd = eval([table_str 'Nd' run_str_temp 'trend' period_str 'exp;']);
                        T_Nd = 10.^exp_Nd * eval([table_str 'Nd' run_str_temp 'trendPct' period_str]);
                        
                    otherwise
                        exp_Nd = eval([table_str 'Nd' run_str 'trend' period_str 'exp;']);
                        T_Nd = 10.^exp_Nd * eval([table_str 'Nd' run_str 'trendPct' period_str]); %Percentage trend in Nd for run in question
                        %I.e., here we are scaling the response in the ACSIS run to the change in
                        %Nd for the UKESM historical, DAMIP, etc.
                end
                
            otherwise
                exp_Nd=NaN;
                T_Nd=NaN;
                
        end
        
end

switch var_str_mat            
    case {'SW_TOA_up','F_SW_TOA_upwelling','totCF','LWPic'}
        %T_Xaero = S_SW * T_Nd * x5 /100;
        eval_str = ['T_Xaero = S_' var_str_tab ' * T_Nd * x5 /100;']; eval(eval_str);
        eval([table_str var_str_tab2 'TrendAeroLocal' period_str ' = deltaT*fscale * T_Xaero;']);                  
        
        %Estimate the SW feedback to temperature
        %eval([table_str trend_str period_str ' = deltaT*fscale*trends.trend_dat_box' hist_str '{1,itr}.coeffs(2);']); %N.B. - latex newcommand names can't have numbers in...
        
        var_str_tab_TS2=[var_str_tab run_str]; %e.g. CFHAD run_str='HAD'
        trend_str_TS = ['TS' run_str 'trend' period_str]; %e.g. 
        trend_str_Y = [var_str_tab run_str 'trend' period_str]; %e.g. SWHistGhgtrendPA - where Y refers to SW, CF, LWPic, etc.
        %trend_str_SW = ['SW' run_str 'trend' period_str]; %e.g. SWHistGhgtrendPA        
        dTS = eval([table_str trend_str_TS]);
        dY = eval([table_str trend_str_Y]);
        eval([table_str var_str_tab 'feedback' run_str period_str ' = dY./dTS;']);
        %eval([table_str 'SWfeedback' run_str period_str ' = dSW./dTS;']);
        
end

switch var_str_mat
    case {'SW_TOA_up','F_SW_TOA_upwelling'}
        T_Xaero = S_SWind * T_Nd * x5 /100;
        eval([table_str var_str_tab2 'IndirectTrendAeroLocal' period_str ' = deltaT*fscale * T_Xaero;']);
end
    

% 
% switch var_str_mat
%     case 'totCF'               
%         T_Xaero = S_CF * T_Nd * x5 /100;
%         eval([table_str var_str_tab2 'TrendAeroLocal' period_str ' = deltaT*fscale * T_Xaero;']);
%         
%     case 'LWPic'               
%         T_Xaero = S_LWPic * T_Nd * x5 /100;
%         eval([table_str var_str_tab2 'TrendAeroLocal' period_str ' = deltaT*fscale * T_Xaero;']);        
%         
%     case {'SW_TOA_up','F_SW_TOA_upwelling'}
%         T_Xaero = S_SW * T_Nd * x5 /100;
%         eval([table_str var_str_tab2 'TrendAeroLocal' period_str ' = deltaT*fscale * T_Xaero;']);  
%         
%         T_Xaero = S_SWind * T_Nd * x5 /100;
%         eval([table_str var_str_tab2 'IndirectTrendAeroLocal' period_str ' = deltaT*fscale * T_Xaero;']); 
%         
%         %Estimate the SW feedback to temperature
%         %eval([table_str trend_str period_str ' = deltaT*fscale*trends.trend_dat_box' hist_str '{1,itr}.coeffs(2);']); %N.B. - latex newcommand names can't have numbers in...
%         
%         var_str_tab_TS2=[var_str_tab run_str]; %e.g. CFHAD run_str='HAD'
%         trend_str_TS=['TS' run_str 'trend' period_str]; %e.g. 
%         trend_str_SW=['SW' run_str 'trend' period_str]; %e.g. SWHistGhgtrendPA
%         dTS = eval([table_str trend_str_TS]);
%         dSW = eval([table_str trend_str_SW]);
%         eval([table_str 'SWfeedback' run_str period_str ' = dSW./dTS;']);
%         
% end






