scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

% pat='C:\matlabR12\work\bauru\tracersJan2005\EGUgraphs\'; 
% pat='g:\runs\sdlavap\results\processrates\iceNCprocesses\MPCratio_iceNC_corr'; 
% %force+3_3th3qv\TotalMR';
% patout='C:\matlabR12\work\bauru\tracersJan2005\EGUgraphs\jpgformat\';

pat='c:/documents and settings/login/my documents/leeds_mmocca/microphysics_results_stewart/images/'; 
patout='c:/documents and settings/login/my documents/leeds_mmocca/microphysics_results_stewart/images/emf'; 

pat='C:\Documents and Settings\Login\My Documents\logbook\vapour_paper\pics\radar_slices\Fig5b_cpz_24feb04_1753_final';
pat='C:\Documents and Settings\Login\My Documents\logbook\vapour_paper\pics\radar_slices\Radar reflectivity (dBZ) 2d plot radar_1_2';  
pat2=[pat '.png'];
%patout='C:\Documents and Settings\Login\My Documents\logbook\vapour_paper\pics\radar_slices\Fig5b_cpz_24feb04_1753_finalJPG.jpg';
patout=[pat 'JPG.jpg']


    
    [im,map]=imread(pat2);

    imwrite(im,patout,'jpg','quality',100);
  
