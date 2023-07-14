function [dat,time_hours,X_flt,Y_flt]=read_mpace_met_files_func(filedir,filename,ncol,NaN_val)

fn=[filedir filename];

fid=fopen(fn,'rt');

for i=1:66
    fgetl(fid);
end

dat=fscanf(fid,'%f',[ncol inf]);

dat(dat>=NaN_val)=NaN;

dat=dat';


time_hours = dat(:,1)/3600;
