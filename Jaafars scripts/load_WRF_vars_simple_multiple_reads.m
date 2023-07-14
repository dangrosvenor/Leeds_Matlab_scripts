%This program for loading WRF output files and then you plot it using other
%program (profile or contour scripts)
%N.B. the mexdcf53 bit I wrote doesn't work for larger file names so are
%changing directory to keep it as short as possible

if exist('cd')==1
    fprintf(1,'variable cd exists - clear it!');
end

is_ecmwf=0;
is_met_em=0;
file_type = 'wrfout';

file_num=[1]; %set file_num for file required


%%%%%%%% set the path and filename here %%%%%%%
%fileWRF(1).file=['C:/Yanbu_modelling/wrfout_d03_Jan0110_for'];   %filename and path of WRF output file
fileWRF(1).file=['Y:\WRF\ecmwf_ml_0.5_nudging\d02'];   %filename and path of WRF output file
fileWRF(2).file=['Y:\WRF\ecmwf_ml_0.5_nudging\d03'];   %filename and path of WRF output file
fileWRF(3).file='Y:\WRF\ecmwf_ml_0.5\met_em_d03-2006-01-06_18-00-00';  %
%%%%%%%%


try
    close(nc);
catch
end




for file_num=1:length(fileWRF)

    nc=netcdf(fileWRF(file_num).file);     
    nca(file_num).nc=netcdf(fileWRF(file_num).file);
    
        try

            if length(strfind(fileWRF(file_num).file,'met_em'))>0 ...
                    | length(strfind(fileWRF(file_num).file,'geo_em'))>0
                file_type = 'met_em';
                file_type2{file_num} = 'met_em';
                is_met_em=1;
            else
                file_type = 'wrfout';
                file_type2{file_num} = 'wrfout';
            end


            if strcmp(fileWRF(file_num).file(1:5),'ecmwf')==1
                is_ecmwf==1;
                file_type = 'ecmwf';
            end

        catch
        end
        
end
 

add_ground_height=0; %value for some plotting scripts

% Times array
Times=nca(1).nc{1}(:);
% Times=datenum(str2num(Times(:,1:4)),str2num(Times(:,6:7)),...
%     str2num(Times(:,9:10)),str2num(Times(:,12:13)),str2num(Times(:,15:16)),...ncnclslfalsdfj
%     str2num(Times(:,18:19)))+8.5./24-tShift./24;
% The times array
time=1;

filestr=fileWRF(1).file;
iund=findstr('_',filestr);
filestr(iund)='-';

iund=findstr('/',filestr);
filestr(iund)='';






switch file_type
    case 'wrfout'
        lat2d=WRFUserARW(nc,'XLAT',time);
        lon2d=WRFUserARW(nc,'XLONG',time);
    case 'met_em'
        lat2d.var=nc{'XLAT_M'}(:);
        lon2d.var=nc{'XLONG_M'}(:);
    case 'ecwmf'
        lon2d.var=nc{'longitude'}(:);
        lat2d.var=nc{'latitude'}(:);
end

if size(Times,2)==1  %if there is only one time then it gets put the wrong way around so change it
    Times=Times';
end



disp('done load WRF');


