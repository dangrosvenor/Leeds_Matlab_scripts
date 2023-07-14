% Make a NetCDF dataset for Odran
% Run read_multi_view_Zhibo first

%Just do for ATEX low SZA case at for now - SZA_20_times_ATEX(itime)
%structure

case_strs = {'SZA_20_times_DYCOMS','SZA_20_times_ATEX'};
iview = 7; %Nadir


for icase=1:length(case_strs)

    case_str = case_strs{icase};

    save_name = ['/home/disk/eos1/d.grosvenor/modis_work/Zhibo_Marshak_work/' case_str '_SAVED_' datestr(now,30) '.mat'];

    vars_to_save = {'R86_1km_3D','R21_1km_3D','R37_1km_3D'...
        ,'R86_std_1km_3D','R21_std_1km_3D','R37_std_1km_3D'...
        ,'H86_1km_3D','H21_1km_3D','H37_1km_3D','Hcov_R86_R21_1km_3D','Hcov_R86_R37_1km_3D'...
        ,'R86_1km_1D','R21_1km_1D','R37_1km_1D'...
        ,'R86_std_1km_1D','R21_std_1km_1D','R37_std_1km_1D'...
        ,'H86_1km_1D','H21_1km_1D','H37_1km_1D','Hcov_R86_R21_1km_1D','Hcov_R86_R37_1km_1D'...
        ,'Re21_1km_1D','Re37_1km_1D','tau_1km_1D'...
        ,'Re21_800m_1D','Re37_800m_1D','tau_800m_1D'...
        ,'Re21_800m_3D','Re37_800m_3D','tau_800m_3D'...
        };

    clear time

    for it=1:eval(['length(' case_str ')'])
        


        time(it)=it; %store a time variable to make it clearer

        for ivar=1:length(vars_to_save)
            var = vars_to_save{ivar};
            if it==1
                eval(['clear ' var]);
            end
        
            eval([var '(:,:,it) = ' case_str '(' num2str(it) ').' var '(iview,:,:);']);
        end


    end


    %save variables to file

    save(save_name,'time','-V7.3');

    for ivar=1:length(vars_to_save)
        var = vars_to_save{ivar};
        eval(['save(save_name,''' var ''',''-V7.3'',''-APPEND'');']);
    end

    %convert to NetCDF
    mat2nc_Dan(save_name,[save_name '.nc'],1); %,1 tells it to give unique dimensions


end




