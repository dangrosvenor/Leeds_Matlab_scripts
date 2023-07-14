%read Marshak data for 2006 paper (JGR - Impact of 3D radiative effects..)

filedir = '/home/disk/eos1/d.grosvenor/modis_work/Zhibo_Marshak_work/';

file_read = 'seg_A.txt';
M = 14;
filename = [filedir file_read];
fid = fopen(filename,'rt');
headersA = textscan(fid,'%s',M);
N = length(headersA{1});
segA = fscanf(fid,'%f',[N Inf]);
fclose(fid);

file_read = 'seg_B.txt';
M = 9;
filename = [filedir file_read];
fid = fopen(filename,'rt');
headersB = textscan(fid,'%s',M);
N = length(headersB{1});
segB = fscanf(fid,'%f',[N Inf]);
fclose(fid);








