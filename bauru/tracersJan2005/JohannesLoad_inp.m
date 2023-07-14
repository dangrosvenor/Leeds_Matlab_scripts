pat='c:/documents and settings/tom/my documents/dan/johannes';

loadcase=[0 0 0 1];

filen='sf4sdlab';

if loadcase(1)==1
	tmpmpc=load([pat '/inpdata/2/sf4sdlab.tmp']);
	tmpmpc=reshape(tmpmpc,[83 703 123]); 
end

if loadcase(2)==1
    wmvmpc=load([pat '/inpdata/2/sf4sdlab.mrwv']);
	wmvmpc=reshape(wmvmpc,[83 703 123]); 
end

if loadcase(3)==1
    prsmpc=load([pat '/inpdata/2/sf4sdlab.prs']);
	prsmpc=reshape(prsmpc,[83 703 123]); 
end

load([pat '/sdla_Grid_Johannes.mat']); %loads x,z,time

szmpc=123;
zmpc=(z(1):(z(end)-z(1))/(szmpc-1):z(end));
sxmpc=703;
xmpc=(x(1):(x(end)-x(1))/(sxmpc-1):x(end));

realstarttime=15.67;
TimeMPC =load([pat,'/', filen,'-Time.DAT'])+realstarttime;

if loadcase(4)==1
    load([pat '/sdla_potemp_johannes_1-83.mat']);
    load([pat '/sdla_pressure_johannes_1-83.mat']);
end


% icemrmpc=load([pat,'/outdata/',filen,'-IWC.DAT']);
% icemrmpc=reshape(icemrmpc,[ 703     1   123    83]); %number density 1e-6/kg
%  
% %xdim=703
% %ydim=1
% %zdim=123
% %tdim=83
% %Sample times written in



% %  tdim=size(Time,1)
% 
% rhompc=load([pat,'/outdata/',filen,'-RHOAIR.DAT']);
% rhompc=reshape(rhompc,[ 703     1   123    83]);
% 
% exname='C:\matlabR12\work\bauru\tracersJan2005\force+3_3th3qv\diags\sdla_Grid_Johannes.mat';
% load(exname); %loads x,z,t
% 
% 
% szmpc=123;
% zmpc=(z(1):(z(end)-z(1))/(szmpc-1):z(end));
% sxmpc=703;
% xmpc=(x(1):(x(end)-x(1))/(sxmpc-1):x(end));