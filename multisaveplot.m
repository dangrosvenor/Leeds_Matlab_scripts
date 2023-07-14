% direc(1).dir='c:/cygwin/home/user/runs/dmidamp_2/';
% direc(2).dir='c:/cygwin/home/user/runs/damp_ccn960/';
% %direc(4).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv/';
% direc(3).dir='c:/cygwin/home/user/runs/damp_inx10/';
% direc(4).dir='g:runs/500mres/';

tag='totbelow5_14-22km_14-Feb06';
%tag='_15.3-18km_29-Feb06';
tag='_13-22km_06-Mar06';
tag='_15.3-17km_06-Mar06_advapscale_clines_clab';


isavemulti=1;


minVal=9e99;
maxVal=-9e99;

for idir=1:4
     plottimeheightvap3;    
     %waterVapourMay2005;
     picname=[direcDan(idir).dir 'results/' savename '_' tag '.emf'];
     set(gcf,'paperpositionmode','auto');
     
     if isavemulti==1
         print(gcf,picname,'-dmeta');
         picname=['c:/documents and settings/login/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         print(gcf,picname,'-dmeta');
     end
     
     if dlogflag==1
         minVal=min([minVal idlog(conts(1),dlogmin)]);
         maxVal=max([maxVal idlog(conts(end),dlogmin)]);
     else
            minVal=min([minVal conts(1)]);
            maxVal=max([maxVal conts(end)]);
     end    
 end
 
 if dlogflag==1
     fprintf(1,'mincovOvr = dlog(%f,dlogmin);\nmaxcovOvr = dlog(%f,dlogmin);',minVal,maxVal);
 else
     fprintf(1,'mincovOvr = %f;\nmaxcovOvr = %f;',minVal,maxVal);
 end
         