filename='/home/mbexddg5/manchester_flt19.txt';

dat=dlmread(filename);
dat=sortrows(dat); %the times in the file are not in order so this puts them in order