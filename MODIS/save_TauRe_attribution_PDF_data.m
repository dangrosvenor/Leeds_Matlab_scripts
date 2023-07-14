%isave=0;

%save_pdfs_filename='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/s
%aved_Arctic_tau_re_Nd_PDFs.mat';
% - specified in watervap, so won't run here too

if isave==0

if length(strfind(xlab,'Mean Optical Depth'))>0
    clear xTau yTau
    
    for i=1:length(xdat)
        xTau(i).x = Xbins';
        yTau(i).y = ydat(i).y.*bin_widths';
    end
    meanTau = meanXvals;
    stdTau = stdXvals;
    Tau_vals = Xvals_save;
    fprintf(1,'\nSaved Tau data\n');    
end

if length(strfind(xlab,'R_{eff 1.6 \mum}'))>0
    clear xRe16 yRe16
    
    for i=1:length(xdat)
        xRe16(i).x = Xbins';
        yRe16(i).y = ydat(i).y.*bin_widths';
    end
    meanRe16 = meanXvals;
    stdRe16 = stdXvals;   
    Re16_vals = Xvals_save;    
    fprintf(1,'\nSaved 1.6 micron Re data\n');    
end

if length(strfind(xlab,'R_{eff 2.1 \mum}'))>0
    clear xRe21 yRe21
    
    for i=1:length(xdat)
        xRe21(i).x = Xbins';
        yRe21(i).y = ydat(i).y.*bin_widths';
    end
    meanRe21 = meanXvals;
    stdRe21 = stdXvals;   
    Re21_vals = Xvals_save;    
    fprintf(1,'\nSaved 2.1 micron Re data\n');    
end

if length(strfind(xlab,'R_{eff 3.7 \mum}'))>0
    clear xRe37 yRe37
    
    for i=1:length(xdat)
        xRe37(i).x = Xbins';
        yRe37(i).y = ydat(i).y.*bin_widths';
    end
    meanRe37 = meanXvals;
    stdRe37 = stdXvals;   
    Re37_vals = Xvals_save;    
    fprintf(1,'\nSaved 3.7 micron Re data\n');
end

if length(strfind(xlab,'N_{d Re 1.6}'))>0
    clear xN16 yNd16
    
    for i=1:length(xdat)
        xNd16(i).x = Xbins';
        yNd16(i).y = ydat(i).y.*bin_widths';
    end
    meanN16 = meanXvals;
    stdN16 = stdXvals;    
    fprintf(1,'\nSaved 1.6 micron Nd data\n');
end

if length(strfind(xlab,'N_d (cm^{-3})'))>0
    clear xN21 yNd21
    
    for i=1:length(xdat)
        xNd21(i).x = Xbins';
        yNd21(i).y = ydat(i).y.*bin_widths';
    end
    meanN21 = meanXvals;
    stdN21 = stdXvals;    
    fprintf(1,'\nSaved 2.1 micron Nd data\n');
end

if length(strfind(xlab,'N_{d Re 3.7}'))>0
    clear xNd37 yNd37
    
    for i=1:length(xdat)
        xNd37(i).x = Xbins';
        yNd37(i).y = ydat(i).y.*bin_widths';
    end
    meanN37 = meanXvals;
    stdN37 = stdXvals;    
    fprintf(1,'\nSaved 3.7 micron Nd data\n');
end


end


if isave==1
    save(save_pdfs_filename,'xTau','yTau','xRe16','yRe16','xRe21','yRe21','xRe37','yRe37','xNd16','yNd16','xNd21','yNd21','xNd37','yNd37','meanN16','stdN16','meanN21','stdN21','meanN37','stdN37','meanRe16','stdRe16','meanRe21','stdRe21','meanRe37','stdRe37','meanTau','stdTau','Tau_vals','Re16_vals','Re21_vals','Re37_vals');
    fprintf(1,'\nWrote data to file\n');
end


