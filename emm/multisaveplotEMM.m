%DON'T USE THIS ONE! use plotEMMtimH.m

% direc(1).dir='c:/cygwin/home/user/runs/dmidamp_2/';
% direc(2).dir='c:/cygwin/home/user/runs/damp_ccn960/';
% %direc(4).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv/';
% direc(3).dir='c:/cygwin/home/user/runs/damp_inx10/';
% direc(4).dir='g:runs/500mres/';

% tag='totbelow5_14-22km_14-Feb06';
% %tag='_15.3-18km_29-Feb06';
% tag='_13-22km_06-Mar06';
tag='nonMinus50';


isavemulti=1;


minVal=9e99;
maxVal=-9e99;

files={'naero' 'ssat' 'lwc' 'ncw' 'rwc' 'nr' 'qi' 'qs' 'qg' 'qi2'};

cat{1}='Cloudwater';
cat{2}='Condensation Freezing';
cat{3}='Contact Nucleation';
cat{4}='Homogeneous Nucleation';
cat{5}='Splinters 1';
cat{6}='Frozen Rain 1';
cat{7}='Splinters 2';
cat{8}='RF Splinters 1';
cat{9}='Frozen Rain 2';
cat{10}='Splinters 3';
cat{11}='RF Splinters 2';

type='isg';
type='inc';

switch type
case 'inc';
ifiles=[1 2 4 6 7 9 10];
for j=1:length(ifiles)
    files{j}=[cat{ifiles(j)}];
end
cols=[1 4:13]; %index in iwczt

case 'other'
    ifiles=[1 2 3 4 5 6];
case 'isg'
    files={'ice' 'snow' 'graupel'};
    ifiles=[1:3];end




ndir=1;

icat2=6;

comp='lacieLap';
comp='readEMM2';

switch comp
case 'uni'
     savedir=['c:/documents and settings/login/my documents/emm_0012/'];
case 'lacieLap'
      savedir=['c:/documents and settings/g/my documents/emm_0012/'];
case 'readEMM2'
      savedir=[rootdir];    
end

for icat=1:length(ifiles)     %icat2:icat2

for idir=1:ndir
     plottimeheightvap3;    
     %waterVapourMay2005;
%     picname=[savedir direcDan(idir).dir '/' savename '_' tag '.emf'];
     picname=[savedir '/' savename '_' tag '.emf'];
     
     set(gcf,'paperpositionmode','auto');
     
     if isavemulti==1
         print(gcf,picname,'-dmeta');
         switch comp
         case 'uni'
             picname=['c:/documents and settings/login/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         case 'lacieLap'
              picname=['c:/documents and settings/g/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         end
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
 
 

end