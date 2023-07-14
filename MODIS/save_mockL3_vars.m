
Npix_str=num2str(Npix);

savedir_mockL3 = '/home/disk/eos8/d.grosvenor/';
now_str=datestr(now,30);
save_filename = [savedir_mockL3 remove_character(file_name_h5,'.','_') '_SAVED_mockL3_' box_type '_' Npix_str 'km_' now_str];

ivar=1; clear mockL3_var

%mean cloud fraction with the lat/lon cell
mockL3_var{ivar}='CF_mockL3'; ivar=ivar+1;

mockL3_var{ivar}='Nice_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='Nreject_mockL3'; ivar=ivar+1;


%1D histograms within the cell
mockL3_var{ivar}='NdPDF_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='WPDF_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='RePDF_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='TauPDF_mockL3'; ivar=ivar+1;

mockL3_var{ivar}='Tau_un_PDF_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='Re_un_PDF_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='Nd_un_PDF_mockL3'; ivar=ivar+1;

% *** calculate a 1D histogram of Nd vs CTP.
%first create a 2D histo of Nd and CTP
mockL3_var{ivar}='Nd_CTP_2D_PDF'; ivar=ivar+1; 
mockL3_var{ivar}='Nd_bins_rep2D'; ivar=ivar+1; 
mockL3_var{ivar}='Nd_vs_CTP_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='NP_vs_CTP_mockL3'; ivar=ivar+1;

%2D tau-reff histograms
mockL3_var{ivar}='Tau_Re_2D_PDF'; ivar=ivar+1;



mockL3_var{ivar}='meanNd_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanTau_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanRe_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanNd_un_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanTau_un_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanRe_un_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanCTT_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='meanSolarZA_mockL3'; ivar=ivar+1;

mockL3_var{ivar}='totNP_vs_CTP'; ivar=ivar+1; 
mockL3_var{ivar}='meanNd_vs_CTP_mockL3'; ivar=ivar+1; 
mockL3_var{ivar}='midCTP_2D_bins'; ivar=ivar+1;

mockL3_var{ivar}='meanNd_mockL3_16'; ivar=ivar+1;
mockL3_var{ivar}='meanNd_mockL3_37'; ivar=ivar+1;
mockL3_var{ivar}='meanRe_mockL3_16'; ivar=ivar+1;
mockL3_var{ivar}='meanRe_mockL3_37'; ivar=ivar+1;
mockL3_var{ivar}='meanW_mockL3'; ivar=ivar+1;
mockL3_var{ivar}='logW_mockL3'; ivar=ivar+1;

mockL3_var{ivar}='LAT_mockL3_edge';
mockL3_var{ivar}='LON_mockL3_edge';


for ivar=1:length(mockL3_var)
    fprintf(1,'Saving %d of %d\n',ivar,length(mockL3_var));
    
    if ivar==1
        eval_str = ['save(save_filename,''' mockL3_var{ivar} ''',''-V7.3'');'];
    else
        eval_str = ['save(save_filename,''' mockL3_var{ivar} ''',''-V7.3'',''-APPEND'');'];
    end
    
    eval(eval_str);
    
end









