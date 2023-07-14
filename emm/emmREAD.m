function a=emmread(emmdir,filestr)

fid=fopen([emmdir filestr],'rt');
cloudtop=dlmread(fn,' ');
fclose(fid);