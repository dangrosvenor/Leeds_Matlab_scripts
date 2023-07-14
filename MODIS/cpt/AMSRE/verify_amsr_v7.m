%   This program reads amsr-e files (version-7) 
%   and writes out the data as shown in the verify_amsre_v7.txt file
%
%	Unzip the data file, if necessary, use matlab function gunzip
%	Change the file_name below to point to your data file
%   Remove or comment out code if you do not have that file

%*****************daily*******************
clear all

file_name = 'your drive:\your directory\amsre_20020715v7';
[time,sst,wspdLF,wspdMF,vapor,cloud,rain]=read_amsr_day_v7(file_name);

time(170:175,274:278,1)'
sst(170:175,274:278,1)'
wspdLF(170:175,274:278,1)'
wspdMF(170:175,274:278,1)'
vapor(170:175,274:278,1)'
cloud(170:175,274:278,1)'
rain(170:175,274:278,1)'


%****************3-day*******************
clear all

file_name = 'your drive:\your directory\amsre_20020715v7_d3d';
[sst,wspdLF,wspdMF,vapor,cloud,rain]=read_amsr_averaged_v7(file_name);

sst(170:175,274:278)'
wspdLF(170:175,274:278)'
wspdMF(170:175,274:278)'
vapor(170:175,274:278)'
cloud(170:175,274:278)'
rain(170:175,274:278)'


%*************weekly********************
clear all

file_name = 'your drive:\your directory\amsre_20020720v7';
[sst,wspdLF,wspdMF,vapor,cloud,rain]=read_amsr_averaged_v7(file_name);

sst(170:175,274:278)'
wspdLF(170:175,274:278)'
wspdMF(170:175,274:278)'
vapor(170:175,274:278)'
cloud(170:175,274:278)'
rain(170:175,274:278)'



%*************monthly**********************
clear all

file_name = 'your drive:\your directory\amsre_200207v7';
[sst,wspdLF,wspdMF,vapor,cloud,rain]=read_amsr_averaged_v7(file_name);

sst(170:175,274:278)'
wspdLF(170:175,274:278)'
wspdMF(170:175,274:278)'
vapor(170:175,274:278)'
cloud(170:175,274:278)'
rain(170:175,274:278)'
