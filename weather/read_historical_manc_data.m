clear filename

%Ringway
%Location: 3814E 3844N, 69 metres amsl
%Estimated data is marked with a * after the value.
%Missing data (more than 2 days missing in month) is marked by  ---.
%Sunshine data taken from an automatic Kipp & Zonen sensor marked with a #, otherwise sunshine data taken from a Campbell Stokes recorder.
%   yyyy  mm   tmax    tmin      af    rain     sun
%              degC    degC    days      mm   hours



idat=1;
ALT_stat(idat).dat = 134; %height in metres (according to http://www.antarctica.ac.uk/met/stations/manned_stations.html)
aws_name(idat).name='Bradford';filename(idat).name='C:\Documents and Settings\dan\My Documents\weather\historical_bradford_data.txt';
num_cols_aws(idat)=7;

idat=idat+1;
ALT_stat(idat).dat = 67; %height in metres (according to http://www.antarctica.ac.uk/met/stations/manned_stations.html)
aws_name(idat).name='Ross-on_Wye';filename(idat).name='C:\Documents and Settings\dan\My Documents\weather\historical_ross_glasto_data.txt';
num_cols_aws(idat)=7;

%idat=idat+1;
%ALT_stat(idat).dat = 69; %height in metres (according to http://www.antarctica.ac.uk/met/stations/manned_stations.html)
%aws_name(idat).name='Ringway';filename(idat).name='C:\Documents and Settings\dan\My Documents\weather\historical_manchester_data.txt';
%num_cols_aws(idat)=7;

for i=1:length(filename)

    fid=fopen(filename(i).name,'rt');
    for ii=1:7
    station_text=fgetl(fid);
    end

    station_dat = [];

    go=1;
    while go==1

        [stationdat_temp,count]=fscanf(fid,'%f'); %reads in until hits some text (or if have reached the end of the file)
        if count>0
            station_dat = [station_dat; stationdat_temp]; %add the data to a 1D vector
        end
        
        [text,count]=fscanf(fid,'%s',[1 1]); %try to read in text in case have 'null'
        if count>0
            if strcmp(text,'*')~=1 %if have an asterix than ignore
                station_dat = [station_dat; NaN]; %null data flag
            end
        end

        if count==0  %if doesn't read any text then must have reached end as previous read reads in all numbers until reaching text
            go=0;
        end

    end




    num_cols=num_cols_aws(i);
    station_dat=reshape(fliplr(station_dat),[num_cols length(station_dat)/num_cols]);


end













'done read of station data'



