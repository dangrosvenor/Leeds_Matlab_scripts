function [time,phase,temp,height,cwc,lwc,iwc,rew,rei,Nd,Ni]=read_mpace_mphys_files_func(filedir,filename,NaN_val)

fid=fopen([filedir filename],'rt');

go=1;
idat=0;
while go==1
    idat=idat+1;

    [line_time{idat}]=textscan(fid,'%s %s %f',1); %time in secs from midnight
    [line_phase{idat}]=textscan(fid,'%s %s %f',1); %
    [line_temp{idat}]=textscan(fid,'%s %f',1'); %
    [line_height{idat}]=textscan(fid,'%s %f',1'); %
    [line_cwc{idat}]=textscan(fid,'%s %f',1'); %
    [line_lwc{idat}]=textscan(fid,'%s %f',1'); %
    [line_iwc{idat}]=textscan(fid,'%s %f',1'); %
    [line_rew{idat}]=textscan(fid,'%s %f',1'); %
    [line_rei{idat}]=textscan(fid,'%s %f',1'); %
    [line_Nd{idat}]=textscan(fid,'%s %f',1'); %Nd per litre
    [line_Ni{idat}]=textscan(fid,'%s %f',1'); %Ni per litre
    
    if size(line_time{idat}{1},1)==0 %when runs out of data to read it returns a zero length variable
        go=0; %exit the loop
    end
    
%    fprintf(1,'%f\n',line_time{idat}{3});
    
    for i=1:27
        fgetl(fid);
    end


end

for j=1:idat-1
    time(j)=line_time{j}{3}/3600; %convert to hours
    phase(j)=line_phase{j}{3};    
    temp(j)=line_temp{j}{2};
    height(j)=line_height{j}{2};
    cwc(j)=line_cwc{j}{2};
    lwc(j)=line_lwc{j}{2};
    iwc(j)=line_iwc{j}{2};
    rew(j)=line_rew{j}{2};
    rei(j)=line_rei{j}{2};
    Nd(j)=line_Nd{j}{2}/1000; %convert to per cc
    Ni(j)=line_Ni{j}{2};
    
end






