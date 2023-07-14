function height=find_height_from_p_ant_d03(pres)

file = 'C:\Documents and Settings\dan\My Documents\WRF\pressure vs. height approx from wrfout_text3.txt';

fid=fopen(file,'rt');
fgetl(fid); %skip first line
fgetl(fid); %skip first line
dat = fscanf(fid,'%f',[4 inf]);

height = interp1(dat(3,:),dat(2,:),pres,'','extrap');

