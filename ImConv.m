scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

% pat='C:\matlabR12\work\bauru\tracersJan2005\EGUgraphs\'; 
% pat='g:\runs\sdlavap\results\processrates\iceNCprocesses\MPCratio_iceNC_corr'; 
% %force+3_3th3qv\TotalMR';
% patout='C:\matlabR12\work\bauru\tracersJan2005\EGUgraphs\jpgformat\';

pat='c:/documents and settings/login/my documents/leeds_mmocca/microphysics_results_stewart/images/'; 
patout='c:/documents and settings/login/my documents/leeds_mmocca/microphysics_results_stewart/images/emf'; 

lis=dir(pat);

for i=3:length(lis)
    
    im=imread(strcat(pat,lis(i).name),'bmp');
    
    %nam=strcat(patout,lis(i).name);
    nam=strcat(patout,lis(i).name);
    
    %figure('position',posit);
    
    %image(im);
    
    %set(gcf,'paperpositionmode','auto');
    %print(gcf,'-epsc',pat2);
    imwrite(im,nam,'jpg','quality',25);
    %fclose(gcf);    
end