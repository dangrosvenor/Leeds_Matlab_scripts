%read in the PACS output hotwire LWC data from APPRAISE
fileName='C:\Documents and Settings\dan\My Documents\logbook\Antarctica\Flights and instruments_Feb2010\02Hotwire_LWC20090218105850.csv';
fileName='\\Data\Database3\APPRAISE_Clouds\BAe146\RawData\B433\CAPS\20090303\20090303105417\02Hotwire_LWC20090303105417.csv';
fid=fopen(fileName,'rt');

%skip 12 lines to the beam area info
for i=1:17
    textline=fgetl(fid);
end

clear data_LWC_PACS

idat=0;
go=1;
while go==1
    idat=idat+1;
    %data = textscan(fid,'%f %f %f',1,'delimiter',',');
    data_LWC(idat)=textscan(fid,'%f',24,'delimiter',',');

    L = length(data_LWC{idat});
    
    if L==0
        go=0;
        break
    end
    
    data_LWC_PACS(1:L,idat)=data_LWC{idat}(:);

end

fclose(fid);

disp('read_PACS_LWC.m finished');

