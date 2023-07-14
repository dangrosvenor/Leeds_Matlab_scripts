%read in the PACS output hotwire LWC data from APPRAISE
fileName='\\Data\Database3\APPRAISE_Clouds\BAe146\RawData\B433\CAPS\20090303\20090303105417\01CAS20090303105417.csv';
fileName='Y:\BAS_flights\29thJune2010_CAS_comparisons\Manchester_CAS_29Jun2010\20100629\20100629102411\01CAS_POL_PBP20100629102411.csv'

fid=fopen(fileName,'rt');

%skip 12 lines to the beam area info
for i=1:23
    textline=fgetl(fid);
end

bin_1=textscan(fid,'%s',1,'delimiter',',');
bin_rest=textscan(fid,'%f',29,'delimiter',',');
bin_1_val = str2num(bin_1{1}{1}(9:end));
bins_PACS = [bin_1_val; bin_rest{1}];

textline=fgetl(fid);
textline=fgetl(fid);
textline=fgetl(fid);

bin_back_1=textscan(fid,'%s',1,'delimiter',',');
bin_back_rest=textscan(fid,'%f',29,'delimiter',',');
bin_1_val = str2num(bin_back_1{1}{1}(9:end));
bins_back_PACS = [bin_1_val; bin_back_rest{1}];


for i=1:8
    textline=fgetl(fid);
end

clear data_CAS_PACS data_CAS

idat=0;
go=1;
while go==1
    idat=idat+1;
    %data = textscan(fid,'%f %f %f',1,'delimiter',',');
    data_CAS(idat)=textscan(fid,'%f',168,'delimiter',',');

    L = length(data_CAS{idat});
    
    if L==0
        go=0;
        break
    end
    
    data_CAS_PACS(1:L,idat)=data_CAS{idat}(:);

end

fclose(fid);

disp('read_PACS_CAS.m finished');

