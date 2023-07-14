% models = {'BCC-ESM1',...
% 'CNRM-CM6-1',...
% 'CNRM-ESM2-1',...
% 'EC-Earth3-AerChem',...
% 'MPI-ESM-1-2-HAM',... 
% 'IPSL',...
% 'MIROC6',...
% 'MIROC-ES2L',...
% 'CESM2',...
% 'CESM2-FV2',...
% 'CESM2-WACCM',...
% 'CESM2-WACCM-FV2',...
% 'NorESM2-MM',...
% 'NIMS-KMAUKESM1-0-LL',...
% 'GFDL-CM4',...
% 'GFDL-ESM4'};
% %'MRI-ESM2-0',... 

istart_multi_model=1;

if iadd_trend_str==1
    trend_str_multi_model=[', trend=(' num2str(fscale*trend_dat_box_multi_model{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_multi_model{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
else
    trend_str_multi_model='';
end

baseline_type = 'zero';
%Sum of perturbations from baseline
switch baseline_type
    case 'common'
        baseline_hist_aer = baseline;
        baseline_hist_GHG = baseline;
        baseline_hist_nat = baseline;
    case 'individual'
        baseline_hist_aer = mean([dat_annual_box_ukesm_hist_aer(istart_multi_model:istart_multi_model+9)]); %use the hist-aer one for the hist-aer line, etc.
        baseline_hist_GHG = mean([dat_annual_box_ukesm_hist_GHG(istart_multi_model:istart_multi_model+9)]); %use the hist-aer one for the hist-aer line, etc.
        baseline_hist_nat = mean([dat_annual_box_ukesm_hist_nat(istart_multi_model:istart_multi_model+9)]); %use the hist-aer one for the hist-aer line, etc.
    case 'zero'
        for imodel=1:length(models)
            expt_str = models{imodel}; expt_str2=remove_character(expt_str,'-','_');
            eval_str =  ['baseline_'  expt_str2 '= 0'];
            eval(eval_str);
        end

end

%cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]}; %10 colours
%cols = {[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]}; %9 colours
cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 1 0],[0 1 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]}; %save black [0 0 0] for obs
cols = {[0.75 0.75 1.0],[0.5 0.5 0.75],[0 1 0],[0 1 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]}; %save black [0 0 0] for obs
line_patts = {'-','--','-.'};

icol=0;
ipatt=1;
for imodel=1:length(models)
    icol=icol+1;
    if icol>length(cols)
        ipatt=ipatt+1;
        icol=1;
    end
    col = cols{icol};
    patt = line_patts{ipatt};
    expt_str = models{imodel}; expt_str2=remove_character(expt_str,'-','_');
    eval_str = ['perts = dat_annual_box_ukesm_' expt_str2 '(istart_multi_model:end) - baseline_' expt_str2 ';']; eval(eval_str);
    eval_str = ['years_multi_model = dat_ukesm_' expt_str2 '.years_ukesm_1d(istart_multi_model:end)']; ; eval(eval_str);
    
%     switch expt_str
%                 case {'CESM2'};
%                     perts = perts*3;                    
%     end
    
    h_multi_model{imodel} = plot(years_multi_model,perts,'color',col,'linestyle',patt);
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=[remove_character(expt_str2,'_',' ') trend_str_multi_model]; ileg=ileg+1;
    set(h_multi_model{imodel},'linewidth',4);
    set(h_multi_model{imodel},'markerfacecolor',col);
    
end









