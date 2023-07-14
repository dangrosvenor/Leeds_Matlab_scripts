pat='c:/documents and settings/tom/my documents/dan/johannes';

filen='sf4sdlab';
icempc=load([pat,'/outdata/',filen,'-ICE.DAT']);
icempc=reshape(icempc,[ 703     1   123    83]); %number density 1e-6/kg

icemrmpc=load([pat,'/outdata/',filen,'-IWC.DAT']);
icemrmpc=reshape(icemrmpc,[ 703     1   123    83]); %number density 1e-6/kg
 
%xdim=703
%ydim=1
%zdim=123
%tdim=83
%Sample times written in
realstarttime=15.67;
TimeMPC =load([pat,'/outdata/',filen,'-Time.DAT'])+realstarttime;
%  tdim=size(Time,1)

rhompc=load([pat,'/outdata/',filen,'-RHOAIR.DAT']);
rhompc=reshape(rhompc,[ 703     1   123    83]);

exname=[pat '\sdla_Grid_Johannes.mat'];
load(exname); %loads x,z,t


szmpc=123;
zmpc=(z(1):(z(end)-z(1))/(szmpc-1):z(end));
sxmpc=703;
xmpc=(x(1):(x(end)-x(1))/(sxmpc-1):x(end));



load([pat '/sdla_icemr']);
load([pat '/sdla_snowmr']);
load([pat '/sdla_graupelmr']);