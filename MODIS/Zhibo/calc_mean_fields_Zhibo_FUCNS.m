

switch case_str_zhibo
    case 'ATEX'
        res={'100m','400m','800m'};
        %ATEX native LES res is 100m and the LES and 100m fields are the
        %same - so ignore the LES fields.
    otherwise
        %FOr Dycoms the LES resolution is 50m
        res={'LES','100m','400m','800m'};
end

rt={'1D','3D'};

clear R86_250m N_vals_R86_250m R86_std_250m R21_250m R21_std_250m cov_R86_R21_250m R37_250m R37_std_250m cov_R86_R37_250m

for ires=1:length(res)
    for irt=1:length(rt)      
        var_str = [res{ires} '_' rt{irt}]; CM_var_str = ['CM_' var_str];
        rt_str = rt{irt};
        %Use only the points where the cloud mask indicates cloud
%        icloud = find( LES_struct.(CM_var_str) >=1 ); %Changed to >= here since in some files 
        % (e.g.
        % Ackerman_DYCOMS2_25bins_dharma_TIME_retrieval_results_sz0.9396.nc
        % ) the cloud mask file is all 1e36 values - perhaps should not
        % trust? Still seem to get rerievals where there is no cloud mask -
        % how can this produce an re value?
        tau_str = ['Tau_' var_str];
        tau = LES_struct.(tau_str);
       
        
      %2.1 um
        re21_str = ['Re21_' var_str];
        re_read = LES_struct.(re21_str);
        R86_str = ['R86_' var_str];
        R86_read = LES_struct.(R86_str);
        R21_str = ['R21_' var_str];
        R21_read = LES_struct.(R21_str);
        R37_str = ['R37_' var_str];
        R37_read = LES_struct.(R37_str);
        %Remove non-sensible values (perhaps occuring due to cloud mask not
        %working very well).
        lower_re_thresh=4; %4
        upper_re_thresh=30; %60; %30
        lower_tau_thresh=0.3; %4
        upper_tau_thresh=150; %30        
        icloud = find( tau>=lower_tau_thresh & tau<=upper_tau_thresh & re_read>=lower_re_thresh & re_read<=upper_re_thresh); %Trying with tau a re limits

        
        re = NaN*ones( size(LES_struct.(re21_str)) );
        R86=re;
        R21=re;
        R37=re;
        re(icloud)=re_read(icloud);
        R86(icloud)=R86_read(icloud);
        R21(icloud)=R21_read(icloud);
        R37(icloud)=R37_read(icloud);
        
        if ires==1 %Just to std dev and cov for the highest resolution base data
            for iview=1:size(R86,1)
                R86_2=squeeze(R86(iview,:,:));
                R21_2=squeeze(R21(iview,:,:));
                R37_2=squeeze(R37(iview,:,:));
                N_av = 20;
                M_av = 20;
                [R86_250m(iview,:,:),N_vals_R86_250m(iview,:,:),R86_std_250m(iview,:,:),R21_250m(iview,:,:),R21_std_250m(iview,:,:),cov_R86_R21_250m(iview,:,:)]=reduce_matrix_subsample_mean(R86_2,N_av,M_av,'covariance',R21_2);
                [R86_250m(iview,:,:),N_vals_R86_250m(iview,:,:),R86_std_250m(iview,:,:),R37_250m(iview,:,:),R37_std_250m(iview,:,:),cov_R86_R37_250m(iview,:,:)]=reduce_matrix_subsample_mean(R86_2,N_av,M_av,'covariance',R37_2);
            end
        end
        

        
%         re(re<lower_re_thresh)=NaN;
%         re(re>upper_re_thresh)=NaN;
        LES_struct.(['mean_' re21_str]) = MeanNoNan(re(:,:),2);
      %3.7 um
        re37_str = ['Re37_' var_str];
        re = NaN*ones( size(LES_struct.(re37_str)) );
        re(icloud)=LES_struct.(re37_str)(icloud);        
        LES_struct.(['mean_' re37_str]) = MeanNoNan(re(:,:),2);  
        
        %Put all of the values from different resolutions in order in a [nres
        %nvza] sized array
        LES_struct.(['mean_Re21_' rt_str '_vs_res'])(ires,:) = LES_struct.(['mean_' re21_str]);
        LES_struct.(['mean_Re37_' rt_str '_vs_res'])(ires,:) = LES_struct.(['mean_' re37_str]);
    
                
    end  
    

      
end
        


% eval(['mean_Re21_3D_' filetag ' = NaN*ones([1 4]);']);
% eval(['mean_Re21_1D_' filetag ' = NaN*ones([1 4]);']);



% icloud = eval(['find(CM_100m_3D_' filetag '==1);']);
% Re21_100m_3D = NaN*ones(size(eval(['Re' wavelength '_100m_3D_' filetag])));
% Re21_100m_3D(icloud)=eval(['Re21_100m_3D_' filetag '(icloud);']);
% 
% icloud = eval(['find(CM_100m_1D_' filetag '==1);']);
% Re21_100m_1D = NaN*ones(size(eval(['Re21_100m_1D_' filetag])));
% Re21_100m_1D(icloud)=eval(['Re21_100m_1D_' filetag '(icloud);']);
% 
% icloud = eval(['find(CM_400m_3D_' filetag '==1);']);
% Re21_400m_3D = NaN*ones(size(eval(['Re21_400m_3D_' filetag])));
% Re21_400m_3D(icloud)=eval(['Re21_400m_3D_' filetag '(icloud);']);
% 
% icloud = eval(['find(CM_400m_1D_' filetag '==1);']);
% Re21_400m_1D = NaN*ones(size(eval(['Re21_400m_1D_' filetag])));
% Re21_400m_1D(icloud)=eval(['Re21_400m_1D_' filetag '(icloud);']);
% 
% icloud = eval(['find(CM_800m_3D_' filetag '==1);']);
% Re21_800m_3D = NaN*ones(size(eval(['Re21_800m_3D_' filetag])));
% Re21_800m_3D(icloud)=eval(['Re21_800m_3D_' filetag '(icloud);']);
% 
% icloud = eval(['find(CM_800m_1D_' filetag '==1);']);
% Re21_800m_1D = NaN*ones(size(eval(['Re21_800m_1D_' filetag])));
% Re21_800m_1D(icloud)=eval(['Re21_800m_1D_' filetag '(icloud);']);
% 
% eval(['mean_Re21_3D_' filetag ' = NaN*ones([1 4]);']);
% eval(['mean_Re21_1D_' filetag ' = NaN*ones([1 4]);']);