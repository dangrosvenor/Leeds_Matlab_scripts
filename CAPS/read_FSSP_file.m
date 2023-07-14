function [bins_FSSP,data_FSSP,time_FSSP]=read_FSSP_file(fileName)

fid=fopen(fileName,'rt');

%skip some lines
for i=1:2
    textline=fgetl(fid);
end

header_line=textscan(fid,'%s ',44);

clear bins_FSSP data_FSSP data_FSSP_temp time_FSSP
%make bins - these are the upper diameters of the bins
for i=1:20
    bins_FSSP(i) = str2num(char(header_line{1}(i+24)));
end




idat=0;
go=1;
while go==1
    idat=idat+1;
    %data = textscan(fid,'%f %f %f',1,'delimiter',',');
    [str] = textscan(fid,'%s ',3);
    
    data_FSSP_temp(idat)=textscan(fid,'%f ',42);

    L = length(data_FSSP_temp{idat});
    
    if L==0
        go=0;
        break
    end
    
    date_FSSP = char(str{1}(1));
    time_FSSP_str=char(str{1}(2));
    ampm = char(str{1}(3));

    
    data_FSSP(1:L,idat)=double(data_FSSP_temp{idat}(:));
    
    icolon=findstr(time_FSSP_str,':');
    
    hour_FSSP = str2num(time_FSSP_str(1:icolon(1)-1));
    
    time_FSSP(idat) = hour_FSSP + str2num(time_FSSP_str(icolon(1)+1:icolon(1)+2))/60 ...
        + str2num(time_FSSP_str(icolon(2)+1:icolon(2)+2))/3600 + data_FSSP(1,idat)/3600e3;

    if strcmp(ampm,'PM')==1 & hour_FSSP~=12
        time_FSSP(idat)=time_FSSP(idat)+12;
    end

    counts=data_FSSP(23:42,idat);
    if length(find(isnan(counts)==1))>0
        data_FSSP(:,idat)=[];
        time_FSSP(idat)=[];
        idat=idat-1;
    end
        
end

fclose(fid);

if size(data_FSSP,2)==0
    clear bins_FSSP data_FSSP data_FSSP_temp time_FSSP
end

disp('read_FSSP.m finished');

