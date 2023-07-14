%Read in the zlevs (height) file for UM

function [n,z,dz] = read_zlevs_UM(file)


fid=fopen(file,'rt'); %is quicker to use fscanf that dlmread
tmp = fgetl(fid);
dat = textscan(fid,'%f %s %f %s %f %s');

n = dat{1};
z = dat{3};
dz = dat{5};

fclose(fid);

