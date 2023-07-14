function data=readwelas(filename)
dir1='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons/';
fid=fopen(filename,'r');
fgetl(fid);fgetl(fid);fgetl(fid);

str1=fgetl(fid);
str2=fgetl(fid);
ind1=findstr(str1,'/');
ind2=findstr(str2,':');

% start time of the data gathering
data.startTime=datenum(str2num(str1(ind1(end)+1:end)),...% date
    str2num(str1(6:ind1(1)-1)),...
    str2num(str1(ind1(1)+1:ind1(2)-1)),... 
    str2num(str2(6:ind2(1)-1)),... % time
    str2num(str2(ind2(1)+1:ind2(2)-1)),...
    str2num(str2(ind2(2)+1:ind2(2)+3)));

if(strcmp(str2(end-1:end),'PM')& ~strcmp(str2(6:ind2(1)-1),'12'))
    data.startTime=data.startTime+12./24;
end

fgetl(fid);fgetl(fid);fgetl(fid);
fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);

str=fgetl(fid);
data.distributionmeasuredelay=str2num(str(findstr(str,'=')+1:end));

fgetl(fid);fgetl(fid);
str=fgetl(fid);
data.name=str(findstr(str,'=')+1:end);
fgetl(fid);fgetl(fid);

str=fgetl(fid);
data.h=str2num(str(findstr(str,'=')+1:end-1));
str=fgetl(fid);
data.l=str2num(str(findstr(str,'=')+1:end-1));
str=fgetl(fid);
data.w=str2num(str(findstr(str,'=')+1:end-1));
str=fgetl(fid);
data.density=str2num(str(findstr(str,'=')+1:findstr(str,'m^')-1));
str=fgetl(fid);
data.sensorflow=str2num(str(findstr(str,'=')+1:findstr(str,'s^')-1));
str=fgetl(fid);
data.particlespeed=str2num(str(findstr(str,'=')+1:findstr(str,'s^')-1));
str=fgetl(fid);
data.dilution=str2num(str(findstr(str,'=')+1:end));
str=fgetl(fid);
data.Tmin=str2num(str(findstr(str,'=')+1:end));
str=fgetl(fid);
data.sensor=str(findstr(str,'=')+1:end);
str=fgetl(fid);
data.calibrationfunction=str(findstr(str,'=')+1:end);
str=fgetl(fid);
data.timescalefactor=str2num(str(findstr(str,'=')+1:findstr(str,'us')-1));

% read few lines we dont need
for i=32:66
    fgetl(fid);
end
% calibration information
for i=1:3
    str=fgetl(fid);
    data.calib{:,i}=str(2:end-1);
    str=fgetl(fid);
    data.name1{:,i}=str(findstr(str,'=')+1:end);
    str=fgetl(fid);
    data.file{:,i}=str(findstr(str,'=')+1:end);
    str=fgetl(fid);
    data.minsize(:,i)=str2num(str(findstr(str,'=')+1:findstr(str,'um')-1));
    str=fgetl(fid);
    data.maxsize(:,i)=str2num(str(findstr(str,'=')+1:findstr(str,'um')-1));
    fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
    
end

for i=94:115
    fgetl(fid);
end
i=1;
data.system1=NaN.*ones(4096,100);
data.system2=NaN.*ones(4096,100);

while 1
    try
        fgetl(fid);fgetl(fid);fgetl(fid);
        str=fgetl(fid);
        ind=findstr(str,':');
        data.time(i)=datenum(str2num(str1(ind1(end)+1:end)),...% date
            str2num(str1(6:ind1(1)-1)),...
            str2num(str1(ind1(1)+1:ind1(2)-1)),... 
            str2num(str(6:ind(1)-1)),... % time
            str2num(str(ind(1)+1:ind(2)-1)),...
            str2num(str(ind(2)+1:ind(2)+3)));

        if(strcmp(str(end-1:end),'PM')& ~strcmp(str(6:ind2(1)-1),'12'))
            data.time(i)=data.time(i)+12./24;
        end
        
        fgetl(fid);
        str=fgetl(fid);
        data.calibrationfunctionact{:,i}=str(findstr(str,'=')+1:end);
        fgetl(fid);
        fgetl(fid);
        str=fgetl(fid);
        data.duration(i)=str2num(str(findstr(str,'=')+1:end-1));
        fgetl(fid);
        str=fgetl(fid);
        data.dkp(i)=str2num(str(findstr(str,'=')+1:end));
        str=fgetl(fid);
        data.sumdN(i)=str2num(str(findstr(str,'=')+1:end));
        str=fgetl(fid);
        try
            data.system1(:,i)=eval(['[',(str(findstr(str,'=')+1:end)),']']);
        catch
            data.system1(:,i)=NaN.*ones(4096,1);
        end
        str=fgetl(fid);
        %data.system2(:,i)=eval(['[',(str(findstr(str,'=')+1:end)),']']);
        try
            data.system2(:,i)=eval(['[',(str(findstr(str,'=')+1:end)),']']);
        catch
            data.system2(:,i)=NaN.*ones(4096,1);
        end
        str=fgetl(fid);
        data.timedistribution(:,i)=eval(['[',(str(findstr(str,'=')+1:end)),']']);
        fgetl(fid);
        
        %data.calibrationfunctionact tells us which size calibration (which size bins)
        %were used - then the appropriate file is loaded
        ind=find(strcmp(data.calib,data.calibrationfunctionact{:,i}));
        
        try
        %%%  need a calibration file here containing the bin sizes or this won't work%%%
        overrule_calib_file='no';        
%        overrule_calib_file='yes';        
        switch overrule_calib_file
            case 'no'
                dat=readcalib([dir1,data.file{:,ind}]);
            case 'yes'
                dir_calib='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\from Stefan Benz\';
                file_calib = 'latex0u25_17.txt';
%                file_calib = 'WELAS_Wasser_kalib_Dan.TXT';                

                dat=readcalib([dir1,data.file{:,3}]);
                
%                dat=readcalib([dir_calib file_calib]);
                

                disp('****** WARNING ******** - Are overriding the proper calibration file');
        end
        data.size(:,i)=dat(2,:);
        end
        
        if(i==29)
            i
        end
        i=i+1;
        
    catch
        break;
    end
end
        
fclose(fid);
data.conc1=data.system1./data.h./data.w./data.particlespeed./data.duration(1);
data.conc2=data.system2./data.h./data.w./data.particlespeed./data.duration(1);

function dat=readcalib(filename)
fid=fopen(filename,'r');
for i=1:8
    fgetl(fid);
end
dat=fscanf(fid,'%f %f',[2 inf]);
fclose(fid);
