%quick version of gcm_process for creating Nd fields etc. with different CF
%screenings

cf_thresh=[0.4 1.01];
cf_thresh=[0.8 1.01];
%cf_thresh=[0.01 1.01];
%cf_thresh=[-0.01 1.01];
%cf_thresh=[0.0 0.8];
%cf_thresh=[0.89 0.91];

cf_type='COSP 2D';
%cf_type='3D CF';
%cf_type='3D CF2';
%cf_type='COSP 2D maxliq';
%cf_type='COSP 2D & 3D maxliq';
%cf_type='';

switch cf_type
    case ''
    otherwise
        Nd = eval(['gcm_drop_read_' gcm_str '*1e-6']);
        %also removing low Nd points because ACWNC sets Nd to zero when the cloud
        %threshold fails (when there is not enough QL or CF). So we want to remove
        %those points. Should really use the FREQL variable, but I don't think that
        %it was output.
        iNd = find(Nd<1);
        Nd(iNd) = NaN;
        cf=eval(['gcm_cf_' gcm_str]);
end



switch cf_type           
    case '3D CF'
       
        icf = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2))']);
        Nd(icf) = NaN;

        %find max over height - if using CAM ANCWC then this is already thresholded
%for reasonable LWC

        Nd = squeeze(max(Nd,[],2));
        
    case '3D CF2' %doing the max Nd first and using the max CF in an attempt to be similar to the COSP approach

        cf=max(cf,[],2);
        icf = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2))']);
        
        Nd=squeeze(max(Nd,[],2));
        Nd(icf) = NaN;

        %find max over height - if using CAM ANCWC then this is already thresholded
%for reasonable LWC

%        Nd = squeeze(max(Nd,[],2));
        
         case 'COSP 2D'
             Nd = squeeze(max(Nd,[],2));
             
             cf = eval(['liqCF_modis_' gcm_str '/100']);
             icf = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2))']);     
             Nd(icf)=NaN;                 
             
    case 'COSP 2D maxliq'     
        cf2 = eval(['gcm_cf_' gcm_str]);        
        icf = find(cf2<0.05);
        
         gcm_liq=eval(['gcm_liq_av_' gcm_str ' ./ gcm_cf_' gcm_str]);
         gcm_liq(icf)=NaN;         
         gcm_liq = permute(gcm_liq,[2 1 3 4]);
         
         
         Nd = permute(Nd,[2 1 3 4]);         
         s_gcm_low = size(gcm_liq);
                  
         [max_liq,imax_liq] = max(gcm_liq,[],1);
         imax_lin = sub2ind([s_gcm_low(1) s_gcm_low(2)*s_gcm_low(3)*s_gcm_low(4)],imax_liq(:),[1:s_gcm_low(2)*s_gcm_low(3)*s_gcm_low(4)]');

         Nd = reshape(Nd(imax_lin),s_gcm_low(2),s_gcm_low(3),s_gcm_low(4)); 


%         Nd = squeeze(max(Nd,[],2));
         
             cf = eval(['liqCF_modis_' gcm_str '/100;']);
             icf = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2));']);     
             Nd(icf)=NaN;   
             
           
             
             
    case 'COSP 2D & 3D maxliq'
        cf2 = eval(['gcm_cf_' gcm_str]);        
        icf = find(cf2<0.05);
        icf2 = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2))']);        
        
        
         gcm_liq=eval(['gcm_liq_av_' gcm_str ' ./ gcm_cf_' gcm_str]);
         gcm_liq(icf)=NaN;  
         gcm_liq(icf2)=NaN;
         gcm_liq = permute(gcm_liq,[2 1 3 4]);
         
         
         Nd = permute(Nd,[2 1 3 4]);         
         s_gcm_low = size(gcm_liq);
                  
         [max_liq,imax_liq] = max(gcm_liq,[],1);
         imax_lin = sub2ind([s_gcm_low(1) s_gcm_low(2)*s_gcm_low(3)*s_gcm_low(4)],imax_liq(:),[1:s_gcm_low(2)*s_gcm_low(3)*s_gcm_low(4)]');

         Nd = reshape(Nd(imax_lin),s_gcm_low(2),s_gcm_low(3),s_gcm_low(4)); 


%         Nd = squeeze(max(Nd,[],2));
         
             cf = eval(['liqCF_modis_' gcm_str '/100']);
             icf = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2))']);     
             Nd(icf)=NaN;    

end

cf = eval(['liqCF_modis_' gcm_str '/100;']);
icf = eval(['find(cf<=cf_thresh(1) | cf>cf_thresh(2));']);

  %using the MODIS COSP CF to scale for in-cloud (maybe this is
             %dubious?)
%             eval(['gcm_lwp_COSPCF_' gcm_str ' = gcm_lwp_' gcm_str './cf;']);
             eval(['gcm_lwp_COSPCF_' gcm_str ' = gcm_lwp_' gcm_str ';']);             
             %remove points outside of the requested cloud range
             eval(['gcm_lwp_COSPCF_' gcm_str '(icf)=NaN;']);
             
             %N.B. there are some points for which the COSP CF is very
             %high, but have zero gcm_lwp??? Perhaps COSP is using the rain
             %field. I.e. rain could still be falling, but cloud has gone?
             
             
              tau = eval(['liqTau_modis_' gcm_str './cf']);
              re = eval(['liqRe_modis_' gcm_str './cf']);
              
              CF_gcm_thresh=0.01;
              
        tau(cf<CF_gcm_thresh)=NaN;
%        tau=permute(tau,[2 3 1]);
        re(cf<CF_gcm_thresh)=NaN;
%        re=permute(re,[2 3 1]);
        
        
%        X = X(ilat,ilon,itime);

%multiply by CF to get the grid-box average as also done for MODIS
%        eval(['LWP_COSP_' gcm_str '= 5/9*1000.*re.*tau.*cf;']);
        %or not
        eval(['LWP_COSP_' gcm_str '= 5/9*1000.*re.*tau;']);        
        eval(['LWP_COSP_' gcm_str '(icf)=NaN;']);
             
             

disp('Done gcm_process_lite');
