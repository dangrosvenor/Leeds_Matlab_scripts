pat='c:/documents and settings/login/my documents/hibiscus/johannes';

filen=strcat(pat,'/sf4sdlab');
ext='';
ice=load(['outdata/',filen,ext,'-ICE.DAT']);
c=reshape(ice,[ 703     1   123    83]); %number density 1e-6/kg
 
%xdim=703
%ydim=1
%zdim=123
%tdim=83
%Sample times written in
Time = load(['outdata/',filen,ext,'-Time.DAT'])+realstarttime;
%  tdim=size(Time,1)
