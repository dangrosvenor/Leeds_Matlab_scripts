function [CPC_time,CPC_counts,CPC_conc]=read_CPC_data(fileName)

flow_rate = 1e3/60; %flow rate in cm^3/sec =1 litre per minute = 1000/60 cm^3 / sec

%fileName='Y:\BAS_flights\flight122\cpc201002261832.txt';
%fileName='Y:\BAS_flights\flight123\cpc201002271757.txt';


% istr=findstr(fileName,'\flight'); %finds the location of the string '\flight'
% flight_no=fileName(istr(2)+7:istr(2)+9); %this should be followed by the flight number
% if strcmp(flight_no(end),'_')==1
%     flight_no=flight_no(1:end-1);
% end

fid=fopen(fileName,'rt');

%CPC_counts=[];

idat=0;
go=1;
ready=0;
while go==1
    
    %%%%%%% Read in the time data - also tells us the probe number   %%%%%%
    cpc_data = textscan(fid,'%f %f %s %f %f %s %s %f',1,'delimiter',',');
    
    if length(cpc_data{1})==0
        break
    end


    if strcmp(cpc_data{7},'READY')
        idat=idat+1;
        CPC_counts(idat)=cpc_data{2};
        CPC_time(idat) = cpc_data{1}/1e3/3600;
%        ready=ready+1;
%        CPC_conc(idat)=NaN;
%        if ready>=2            
%            CPC_conc(idat) = CPC_counts(idat)/(CPC_time(idat)-CPC_time(idat-1));
%        end
%    else
%        ready=0;
    end
    
    
end

CPC_conc = CPC_counts/flow_rate; %counts are the number of counts per second
%CPC_conc = CPC_conc/flow_rate; %counts are the number of counts per second
%so divide by flow rate to get numbers per cc

