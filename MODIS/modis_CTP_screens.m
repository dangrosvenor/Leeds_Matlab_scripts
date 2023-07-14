thresh_CF=[0.8 1.0001];
thresh_NP=50;
%thresh_CTPs_1 = [880];
%thresh_CTPs_2 = [1000];

thresh_CTPs_1 = [NaN];
thresh_CTPs_2 = [NaN];


%thresh_CTPs_1 = [880  830  780  730  680  440 50];
%thresh_CTPs_2 = [1000 1000 1000 1000 1000 680 440];

screen_types{1}='NP + CF';
%screen_types{1}='NP + CF + mean pressure';
%screen_types{2}='NP + CF + (mean pressure - std_dev)';
%screen_types{3}='NP + CF + (mean pressure - 2*std_dev)';

%for plot_global
if exist('iCTP_screen')
    CTP_str = screen_types{iCTP_screen};
    CTP_thresh_str = [num2str(min(thresh_CTPs_1(iCTPs))) ' to ' num2str(max(thresh_CTPs_2(iCTPs))) ' hPa CF.GT.' num2str(thresh_CF(1))];
end

month_str='';
if exist('imonths')
    for imonth=1:length(imonths)
        month_str=[month_str ' ' datestr(datenum(2000,imonths(imonth),1),'mmm')];
    end
end